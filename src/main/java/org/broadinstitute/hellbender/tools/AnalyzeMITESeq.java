package org.broadinstitute.hellbender.tools;

import htsjdk.samtools.*;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.BetaFeature;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.barclay.help.DocumentedFeature;
import org.broadinstitute.hellbender.cmdline.programgroups.CoverageAnalysisProgramGroup;
import org.broadinstitute.hellbender.engine.GATKTool;
import org.broadinstitute.hellbender.engine.ReferenceDataSource;
import org.broadinstitute.hellbender.engine.filters.ReadFilterLibrary;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.tools.spark.utils.HopscotchMap;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.gcs.BucketUtils;
import org.broadinstitute.hellbender.utils.read.GATKRead;

import java.io.BufferedOutputStream;
import java.io.IOException;
import java.io.OutputStreamWriter;
import java.text.DecimalFormat;
import java.util.*;
import java.util.stream.Stream;

@DocumentedFeature
@CommandLineProgramProperties(
        summary = "(Experimental) Processes reads from a MITESeq experiment.",
        oneLineSummary = "(EXPERIMENTAL) Processes reads from a MITESeq experiment.",
        programGroup = CoverageAnalysisProgramGroup.class
)
@BetaFeature
public class AnalyzeMITESeq extends GATKTool {
    @Argument(doc = "minimum quality score for analyzed portion of read", fullName = "min-q")
    private static int minQ = 30;

    @Argument(doc = "minimum size of high-quality portion of read", fullName = "min-length")
    private static int minLength = 15;

    @Argument(doc = "minimum number of wt calls flanking variant", fullName = "min-flanking-length")
    private static int minFlankingLength = 18;

    @Argument(doc = "reference indices of the ORF (1-based, closed), for example, '134-180,214-238'", fullName = "orf")
    private static String orfCoords;

    @Argument(doc = "minimum number of observations of reported variants", fullName = "min-variant-obs")
    private static long minVariantObservations = 0;

    @Argument(doc = "codon translation (a string of 64 amino acid codes", fullName = "codon-translation")
    private static String codonTranslation = "KNKNTTTTRSRSIIMIQHQHPPPPRRRRLLLLEDEDAAAAGGGGVVVVZYZYSSSSZCWCLFLF";

    @Argument(doc = "output file prefix", fullName = "output-file-prefix", shortName = "O")
    private static String outputFilePrefix;

    @Argument(doc = "paired mode", fullName = "paired-mode")
    private static boolean pairedMode = true;

    private byte[] refSeq; // the amplicon -- all bytes are converted to upper-case 'A', 'C', 'G', or 'T', no nonsense
    private long[] refCoverage; // number of molecules aligning to each reference position -- same length as above
    private long[] refCoverageSizeHistogram; // number of molecules having a given reference coverage size

    private CodonTracker codonTracker;

    // a map of SNV sets onto number of observations of that set of variations
    private final HopscotchMap<SNVCollectionCount, Long, SNVCollectionCount> variationCounts =
            new HopscotchMap<>(10000000);

    private long nReadsTotal = 0; // total number of reads in input bam
    private long nReadsUnmapped = 0; // number of unmapped reads in bam
    private long nReadsLowQuality = 0; // number of reads where trimming made the read disappear

    private long nTotalBaseCalls = 0; // number of base calls over all reads in bam (mapped or not, call trimmed or not)

    // a bunch of mutually exclusive counts of molecules
    private long nWildTypeMolecules = 0; // number of molecules in which no variation from reference was detected
    private long nInconsistentPairs = 0; // number of molecules where mates with conflicting variants in overlap region
    private long nInsufficientFlankMolecules = 0; // number of molecules where variation was too close to end of region
    private long nLowQualityVariantMolecules = 0; // number of molecules where a variation was called with low quality
    private long nCalledVariantMolecules = 0; // number of molecules with a least one variant

    private IntervalCounter intervalCounter; // count of molecules having a particular [start, stop) on reference

    private GATKRead read1; // a place to stash the first read of a pair during pairwise processing of the read stream

    private static final int UPPERCASE_MASK = 0xDF; // e.g., 'a' & UPPERCASE_MASK == 'A'
    private final static byte NO_CALL = (byte)'-';

    private static final int N_REGULAR_CODONS = 64;
    private static final int FRAME_PRESERVING_INDEL_INDEX = 64;
    private static final int FRAME_SHIFTING_INDEL_INDEX = 65;
    private static final int CODON_COUNT_ROW_SIZE = 66;

    private static final String[] labelForCodonValue = {
        "AAA", "AAC", "AAG", "AAT", "ACA", "ACC", "ACG", "ACT", "AGA", "AGC", "AGG", "AGT", "ATA", "ATC", "ATG", "ATT",
        "CAA", "CAC", "CAG", "CAT", "CCA", "CCC", "CCG", "CCT", "CGA", "CGC", "CGG", "CGT", "CTA", "CTC", "CTG", "CTT",
        "GAA", "GAC", "GAG", "GAT", "GCA", "GCC", "GCG", "GCT", "GGA", "GGC", "GGG", "GGT", "GTA", "GTC", "GTG", "GTT",
        "TAA", "TAC", "TAG", "TAT", "TCA", "TCC", "TCG", "TCT", "TGA", "TGC", "TGG", "TGT", "TTA", "TTC", "TTG", "TTT"
    };

    @Override
    public boolean requiresReads() {
        return true;
    }

    @Override
    public boolean requiresReference() {
        return true;
    }

    @Override
    public void onTraversalStart() {
        super.onTraversalStart();
        initializeRefSeq();
        codonTracker = new CodonTracker(orfCoords, refSeq);
        if ( codonTranslation.length() != N_REGULAR_CODONS ) {
            throw new UserException("codon-translation string must contain exactly 64 characters");
        }
    }

    @Override
    public void traverse() {
        // ignore unaligned reads and non-primary alignments
        final Stream<GATKRead> reads = getTransformedReadStream(ReadFilterLibrary.PRIMARY_LINE);

        if ( !pairedMode ) {
            reads.forEach(read -> {
                try {
                    applyReport(getReadReport(read));
                } catch ( final Exception e ) {
                    final String readName = read.getName();
                    throw new GATKException("Caught unexpected exception on read " + nReadsTotal + ": " + readName, e);
                }
            });
        } else {
            reads.forEach(read -> {
                try {
                    if ( !read.isPaired() ) {
                        applyReport(getReadReport(read));
                        return;
                    }
                    if ( read1 == null ) {
                        read1 = read;
                        return;
                    }
                    if ( !read1.getName().equals(read.getName()) ) {
                        System.err.println("Read " + read1.getName() + " has no mate.");
                        applyReport(getReadReport(read));
                        read1 = read;
                        return;
                    }
                    final ReadReport report1 = getReadReport(read1);
                    final ReadReport report2 = getReadReport(read);
                    final ReadReport combinedReport = combineReports(report1, report2);
                    if ( combinedReport != null ) {
                        applyReport(combinedReport);
                    } else {
                        applyReport(report1);
                        applyReport(report2);
                    }
                    read1 = null;
                } catch ( final Exception e ) {
                    final String readName = read.getName();
                    throw new GATKException("Caught unexpected exception on read " + nReadsTotal + ": " + readName, e);
                }
            });
        }
        if ( read1 != null ) {
            if ( read1.isPaired() ) {
                System.err.println("Read " + read1.getName() + " has no mate.");
            }
            try {
                applyReport(getReadReport(read1));
            } catch ( final Exception e ) {
                final String readName = read1.getName();
                throw new GATKException("Caught unexpected exception on read " + nReadsTotal + ": " + readName, e);
            }
        }
    }

    @Override
    public Object onTraversalSuccess() {
        writeVariationCounts(getVariationEntries());
        writeRefCoverage();
        writeCodonCounts();
        writeCodonFractions();
        writeAACounts();
        writeAAFractions();
        writeReadCounts();
        writeCoverageSizeHistogram();
        return null;
    }

    private List<SNVCollectionCount> getVariationEntries() {
        final long outputSize =
                variationCounts.stream().filter(entry -> entry.getCount() >= minVariantObservations).count();
        final List<SNVCollectionCount> variationEntries = new ArrayList<>((int)outputSize);
        for ( final SNVCollectionCount entry : variationCounts ) {
            final long count = entry.getCount();
            if ( count >= minVariantObservations ) {
                variationEntries.add(entry);
            }
        }
        variationEntries.sort(Comparator.naturalOrder());
        return variationEntries;
    }

    private void writeVariationCounts( final List<SNVCollectionCount> variationEntries ) {
        final String variantsFile = outputFilePrefix + ".variantCounts";
        try ( final OutputStreamWriter writer =
                      new OutputStreamWriter(new BufferedOutputStream(BucketUtils.createFile(variantsFile))) ) {
            final DecimalFormat formatter = new DecimalFormat("0.0");
            for ( final SNVCollectionCount entry : variationEntries ) {
                writer.write(Long.toString(entry.getCount()));
                writer.write('\t');
                final List<SNV> snvs = entry.getSNVs();
                final int start = snvs.get(0).getRefIndex() - minFlankingLength;
                final int end = snvs.get(snvs.size() - 1).getRefIndex() + minFlankingLength;
                writer.write(Long.toString(intervalCounter.countSpanners(start, end)));
                writer.write('\t');
                writer.write(Integer.toString(snvs.stream().mapToInt(SNV::getQuality).sum()));
                writer.write('\t');
                writer.write(formatter.format(entry.getMeanRefCoverage()));
                writer.write('\t');
                writer.write(Integer.toString(snvs.size()));
                String sep = "\t";
                for ( final SNV snv : snvs ) {
                    writer.write(sep);
                    sep = ", ";
                    writer.write(snv.toString());
                }
                describeVariantsAsCodons(writer, snvs);
                writer.write('\n');
            }
        } catch ( final IOException ioe ) {
            throw new UserException("Can't write " + variantsFile, ioe);
        }
    }

    private void describeVariantsAsCodons( final OutputStreamWriter writer, final List<SNV> snvs ) throws IOException {
        final List<CodonVariation> variations = codonTracker.encodeSNVsAsCodons(snvs);
        final int[] refCodonValues = codonTracker.getRefCodonValues();

        final long nVariations = variations.stream().filter(cv -> !cv.isFrameshift()).count();
        writer.write('\t');
        writer.write(Long.toString(nVariations));

        String sep = "\t";
        for ( final CodonVariation variation : variations ) {
            if ( variation.isFrameshift() ) continue;
            writer.write(sep);
            sep = ", ";
            final int codonId = variation.getCodonId();
            writer.write(Integer.toString(codonId));
            writer.write(':');
            writer.write(variation.isInsertion() ? "---" : labelForCodonValue[refCodonValues[codonId]]);
            writer.write('>');
            writer.write(variation.isDeletion() ? "---" : labelForCodonValue[variation.getCodonValue()]);
        }

        sep = "\t";
        for ( final CodonVariation variation : variations ) {
            if ( variation.isFrameshift() ) continue;
            writer.write(sep);
            sep = ", ";
            final int codonId = variation.getCodonId();
            if ( variation.isInsertion() ) {
                writer.write("I:->");
                writer.write(codonTranslation.charAt(variation.getCodonValue()));
            } else if ( variation.isDeletion() ) {
                writer.write("D:");
                writer.write(codonTranslation.charAt(refCodonValues[codonId]));
                writer.write(":-");
            } else {
                final char fromAA = codonTranslation.charAt(refCodonValues[codonId]);
                final char toAA = codonTranslation.charAt(variation.getCodonValue());
                final char label = fromAA == toAA ? 'S' : CodonTracker.isStop(variation.getCodonValue()) ? 'N' : 'M';
                writer.write(label);
                writer.write(':');
                writer.write(fromAA);
                writer.write('>');
                writer.write(toAA);
            }
        }
    }

    private void writeRefCoverage() {
        final String refCoverageFile = outputFilePrefix + ".refCoverage";
        try ( final OutputStreamWriter writer =
                      new OutputStreamWriter(new BufferedOutputStream(BucketUtils.createFile(refCoverageFile))) ) {
            writer.write("RefPos\tCoverage\n");
            int refPos = 1;
            for ( final long coverage : refCoverage ) {
                writer.write(Integer.toString(refPos++));
                writer.write('\t');
                writer.write(Long.toString(coverage));
                writer.write('\n');
            }
        } catch ( final IOException ioe ) {
            throw new UserException("Can't write " + refCoverageFile, ioe);
        }
    }

    private void writeCodonCounts() {
        final String codonCountsFile = outputFilePrefix + ".codonCounts";
        try ( final OutputStreamWriter writer =
                      new OutputStreamWriter(new BufferedOutputStream(BucketUtils.createFile(codonCountsFile))) ) {
            final long[][] codonCounts = codonTracker.getCodonCounts();
            final int nCodons = codonCounts.length;
            writer.write("AAA\tAAC\tAAG\tAAT\tACA\tACC\tACG\tACT\tAGA\tAGC\tAGG\tAGT\tATA\tATC\tATG\tATT\t" +
                    "CAA\tCAC\tCAG\tCAT\tCCA\tCCC\tCCG\tCCT\tCGA\tCGC\tCGG\tCGT\tCTA\tCTC\tCTG\tCTT\t" +
                    "GAA\tGAC\tGAG\tGAT\tGCA\tGCC\tGCG\tGCT\tGGA\tGGC\tGGG\tGGT\tGTA\tGTC\tGTG\tGTT\t" +
                    "TAA\tTAC\tTAG\tTAT\tTCA\tTCC\tTCG\tTCT\tTGA\tTGC\tTGG\tTGT\tTTA\tTTC\tTTG\tTTT\t" +
                    "NFS\tFS\tTotal\n");
            for ( int codonId = 0; codonId != nCodons; ++codonId ) {
                final long[] rowCounts = codonCounts[codonId];
                long total = 0;
                for ( final long count : rowCounts ) {
                    writer.write(Long.toString(count));
                    writer.write('\t');
                    total += count;
                }
                writer.write(Long.toString(total));
                writer.write('\n');
            }
        } catch ( final IOException ioe ) {
            throw new UserException("Can't write " + codonCountsFile, ioe);
        }
    }

    private void writeCodonFractions() {
        final String codonFractFile = outputFilePrefix + ".codonFractions";
        try ( final OutputStreamWriter writer =
                      new OutputStreamWriter(new BufferedOutputStream(BucketUtils.createFile(codonFractFile))) ) {
            final long[][] codonCounts = codonTracker.getCodonCounts();
            final int nCodons = codonCounts.length;
            writer.write("Codon" + "" +
                    "   AAA   AAC   AAG   AAT   ACA   ACC   ACG   ACT" +
                    "   AGA   AGC   AGG   AGT   ATA   ATC   ATG   ATT" +
                    "   CAA   CAC   CAG   CAT   CCA   CCC   CCG   CCT" +
                    "   CGA   CGC   CGG   CGT   CTA   CTC   CTG   CTT" +
                    "   GAA   GAC   GAG   GAT   GCA   GCC   GCG   GCT" +
                    "   GGA   GGC   GGG   GGT   GTA   GTC   GTG   GTT" +
                    "   TAA   TAC   TAG   TAT   TCA   TCC   TCG   TCT" +
                    "   TGA   TGC   TGG   TGT   TTA   TTC   TTG   TTT" +
                    "   NFS    FS    Total\n");
            final DecimalFormat formatter = new DecimalFormat("##0.00");
            for ( int codonId = 0; codonId != nCodons; ++codonId ) {
                writer.write(String.format("%5d", codonId + 1));
                final long[] rowCounts = codonCounts[codonId];
                final long total = Arrays.stream(rowCounts).sum();
                for ( final long count : rowCounts ) {
                    writer.write(formatter.format(100. * count / total));
                }
                writer.write(String.format("%9d\n", total));
            }
        } catch ( final IOException ioe ) {
            throw new UserException("Can't write " + codonFractFile, ioe);
        }
    }

    private void writeAACounts() {
        final String aaCountsFile = outputFilePrefix + ".aaCounts";
        try ( final OutputStreamWriter writer =
                      new OutputStreamWriter(new BufferedOutputStream(BucketUtils.createFile(aaCountsFile))) ) {
            final long[][] codonCounts = codonTracker.getCodonCounts();
            final int nCodons = codonCounts.length;
            for ( int codonId = 0; codonId != nCodons; ++codonId ) {
                final long[] rowCounts = codonCounts[codonId];
                final SortedMap<Character, Long> aaCounts = new TreeMap<>();
                for ( int codonValue = 0; codonValue != N_REGULAR_CODONS; ++codonValue ) {
                    aaCounts.merge(codonTranslation.charAt(codonValue), rowCounts[codonValue], Long::sum);
                }
                if ( codonId == 0 ) {
                    String prefix = "";
                    for ( final char chr : aaCounts.keySet() ) {
                        writer.write(prefix);
                        prefix = "\t";
                        writer.write(chr);
                    }
                    writer.write('\n');
                }
                String prefix = "";
                for ( final long count : aaCounts.values() ) {
                    writer.write(prefix);
                    prefix = "\t";
                    writer.write(Long.toString(count));
                }
                writer.write('\n');
            }
        } catch ( final IOException ioe ) {
            throw new UserException("Can't write " + aaCountsFile, ioe);
        }
    }

    private void writeAAFractions() {
        final String aaFractFile = outputFilePrefix + ".aaFractions";
        try ( final OutputStreamWriter writer =
                      new OutputStreamWriter(new BufferedOutputStream(BucketUtils.createFile(aaFractFile))) ) {
            final DecimalFormat formatter = new DecimalFormat("##0.00");
            final long[][] codonCounts = codonTracker.getCodonCounts();
            final int nCodons = codonCounts.length;
            for ( int codonId = 0; codonId != nCodons; ++codonId ) {
                final long[] rowCounts = codonCounts[codonId];
                final SortedMap<Character, Long> aaCounts = new TreeMap<>();
                for ( int codonValue = 0; codonValue != N_REGULAR_CODONS; ++codonValue ) {
                    aaCounts.merge(codonTranslation.charAt(codonValue), rowCounts[codonValue], Long::sum);
                }
                if ( codonId == 0 ) {
                    writer.write("Codon");
                    for ( final char chr : aaCounts.keySet() ) {
                        writer.write("     ");
                        writer.write(chr);
                    }
                    writer.write("    Total\n");
                }
                writer.write(String.format("%5d", codonId + 1));
                final long total = Arrays.stream(rowCounts).sum();
                for ( final long count : aaCounts.values() ) {
                    writer.write(formatter.format(100. * count / total));
                }
                writer.write(String.format("%9d\n", total));
            }
        } catch ( final IOException ioe ) {
            throw new UserException("Can't write " + aaFractFile, ioe);
        }
    }

    private void writeReadCounts() {
        final String readCountsFile = outputFilePrefix + ".readCounts";
        try ( final OutputStreamWriter writer =
                      new OutputStreamWriter(new BufferedOutputStream(BucketUtils.createFile(readCountsFile))) ) {
            final DecimalFormat df = new DecimalFormat("0.000");
            writer.write("Total Reads:\t" + nReadsTotal + "\n");
            writer.write("Unmapped Reads:\t" + nReadsUnmapped + "\t" +
                    df.format(100. * nReadsUnmapped / nReadsTotal) + "%\n");
            writer.write("LowQ Reads:\t" + nReadsLowQuality + "\t" +
                    df.format(100. * nReadsLowQuality / nReadsTotal) + "%\n");
            final long totalMolecules = nInconsistentPairs + nWildTypeMolecules + nInsufficientFlankMolecules +
                    nLowQualityVariantMolecules + nCalledVariantMolecules;
            writer.write("Number of inconsistent pair molecules:\t" + nInconsistentPairs + "\t" +
                    df.format(100. * nInconsistentPairs / totalMolecules) + "%\n");
            writer.write("Number of wild type molecules:\t" + nWildTypeMolecules + "\t" +
                    df.format(100. * nWildTypeMolecules / totalMolecules) + "%\n");
            writer.write("Number of insufficient flank molecules:\t" + nInsufficientFlankMolecules + "\t" +
                    df.format(100. * nInsufficientFlankMolecules / totalMolecules) + "%\n");
            writer.write("Number of low quality variation molecules:\t" + nLowQualityVariantMolecules + "\t" +
                    df.format(100. * nLowQualityVariantMolecules / totalMolecules) + "%\n");
            writer.write("Number of called variant molecules:\t" + nCalledVariantMolecules + "\t" +
                    df.format(100. * nCalledVariantMolecules / totalMolecules) + "%\n");
            writer.write("Base calls evaluated for variants:\t" +
                    df.format(100. * Arrays.stream(refCoverage).sum() / nTotalBaseCalls) + "%\n");
        } catch ( final IOException ioe ) {
            throw new UserException("Can't write " + readCountsFile, ioe);
        }
    }

    private void writeCoverageSizeHistogram() {
        if ( refCoverageSizeHistogram == null ) return;

        final String trimCountsFile = outputFilePrefix + ".coverageLengthCounts";
        try ( final OutputStreamWriter writer =
                      new OutputStreamWriter(new BufferedOutputStream(BucketUtils.createFile(trimCountsFile))) ) {
            final int len = refCoverageSizeHistogram.length;
            for ( int idx = minLength; idx != len; ++idx ) {
                writer.write(Integer.toString(idx));
                writer.write('\t');
                writer.write(Long.toString(refCoverageSizeHistogram[idx]));
                writer.write('\n');
            }
        } catch ( final IOException ioe ) {
            throw new UserException("Can't write " + trimCountsFile, ioe);
        }
    }

    // describes an interval on some sequence as a pair of offsets (0-based, half-open).
    private final static class Interval {
        private final int start;
        private final int end;

        public Interval( final int start, final int end ) {
            this.start = start;
            this.end = end;
        }

        public int getStart() { return start; }
        public int getEnd() { return end; }
        public int size() { return end - start; }

        public static Interval nullInterval = new Interval(0, 0);
    }

    // a description of a single-base deviation from reference
    private static final class SNV implements Comparable<SNV> {
        private final int refIndex;
        private final byte refCall;
        private final byte variantCall;
        private final byte qual; // not a part of equality or hash

        public SNV( final int refIndex, final byte refCall, final byte variantCall, final byte qual ) {
            this.refIndex = refIndex;
            this.refCall = refCall;
            this.variantCall = variantCall;
            this.qual = qual;
        }

        public int getRefIndex() { return refIndex; }
        public byte getRefCall() { return refCall; }
        public byte getVariantCall() { return variantCall; }
        public byte getQuality() { return qual; }

        @Override
        public int hashCode() {
            return 47 * (47 * (47 * refIndex + refCall) + variantCall);
        }

        @Override
        public boolean equals( final Object obj ) {
            return obj instanceof SNV && equals((SNV)obj);
        }

        public boolean equals( final SNV that ) {
            return this.refIndex == that.refIndex &&
                    this.refCall == that.refCall &&
                    this.variantCall == that.variantCall;
        }

        @Override
        public int compareTo( final SNV that ) {
            int result = Integer.compare(this.refIndex, that.refIndex);
            if ( result == 0 ) result = Byte.compare(this.refCall, that.refCall);
            if ( result == 0 ) result = Byte.compare(this.variantCall, that.variantCall);
            return result;
        }

        @Override
        public String toString() {
            return (refIndex + 1) + ":" + (char)refCall + ">" + (char)variantCall;
        }
    }

    // a count of molecules that start and stop at particular places on the reference
    // this allows us to calculate the number of molecules that span any given interval
    private static final class IntervalCounter {
        // triangular matrix indexed first by starting position on reference, and second by the size of the interval
        final long[][] counts;

        public IntervalCounter( final int refLen ) {
            counts = new long[refLen][];
            for ( int rowIndex = 0; rowIndex != refLen; ++rowIndex ) {
                counts[rowIndex] = new long[refLen - rowIndex + 1];
            }
        }

        public void addCount( final int refStart, final int refEnd ) {
            counts[refStart][refEnd - refStart] += 1;
        }

        public long countSpanners( final int refStart, final int refEnd ) {
            long total = 0;
            for ( int rowIndex = 0; rowIndex <= refStart; ++rowIndex ) {
                final long[] row = counts[rowIndex];
                for ( int spanIndex = refEnd - rowIndex; spanIndex < row.length; ++spanIndex ) {
                    total += row[spanIndex];
                }
            }
            return total;
        }
    }

    // an array of SNVs that serves as a key (comparison, equality, and hashCode depend only on this part)
    // and a count of the number of observations of those SNVs, plus the total reference coverage over all observations
    // implements a Map.Entry as a single object to conserve memory
    private static final class SNVCollectionCount
            implements Map.Entry<SNVCollectionCount, Long>, Comparable<SNVCollectionCount> {
        private static final SNV[] emptyArray = new SNV[0];
        private final SNV[] snvs;
        private long count; // number of observations of this set of SNVs
        private int totalRefCoverage; // the sum of the reference coverage over all observations
        private final int hash;

        public SNVCollectionCount( final List<SNV> snvs, final int refCoverage ) {
            this.snvs = snvs.toArray(emptyArray);
            this.count = 1;
            this.totalRefCoverage = refCoverage;
            int hashVal = 0;
            for ( final SNV snv : snvs ) {
                hashVal = 47 * hashVal + snv.hashCode();
            }
            hash = 47 * hashVal;
        }

        @Override
        public SNVCollectionCount getKey() { return this; }

        @Override
        public Long getValue() { return count; }

        @Override
        public Long setValue( final Long value ) {
            final Long result = count;
            count = value;
            return result;
        }

        public List<SNV> getSNVs() { return Arrays.asList(snvs); }
        public long getCount() { return count; }

        public void bumpCount( final int refCoverage ) {
            count += 1;
            totalRefCoverage += refCoverage;
        }

        public float getMeanRefCoverage() {
            return 1.f * totalRefCoverage / count;
        }

        @Override
        public boolean equals( final Object obj ) {
            return obj instanceof SNVCollectionCount && equals((SNVCollectionCount)obj);
        }

        public boolean equals( final SNVCollectionCount that ) {
            return this.hash == that.hash && Arrays.equals(this.snvs, that.snvs);
        }

        @Override
        public int hashCode() {
            return hash;
        }

        @Override
        public int compareTo( final SNVCollectionCount that ) {
            final int minSize = Math.min(this.snvs.length, that.snvs.length);
            int result = 0;
            for ( int idx = 0; idx != minSize; ++idx ) {
                result = this.snvs[idx].compareTo(that.snvs[idx]);
                if ( result != 0 ) break;
            }
            if ( result == 0 ) result = Integer.compare(this.snvs.length, that.snvs.length);
            return result;
        }
    }

    private enum CodonVariationType {
        FRAMESHIFT,
        INSERTION,
        DELETION,
        MODIFICATION
    }

    private static final class CodonVariation {
        private final int codonId;
        private final int codonValue;
        private final CodonVariationType variationType;

        public CodonVariation( final int codonId, final int codonValue, final CodonVariationType variationType ) {
            this.codonId = codonId;
            this.codonValue = codonValue;
            this.variationType = variationType;
        }

        public int getCodonId() { return codonId; }
        public int getCodonValue() { return codonValue; }
        public CodonVariationType getVariationType() { return variationType; }
        public boolean isFrameshift() { return variationType == CodonVariationType.FRAMESHIFT; }
        public boolean isInsertion() { return variationType == CodonVariationType.INSERTION; }
        public boolean isDeletion() { return variationType == CodonVariationType.DELETION; }
        public boolean isModification() { return variationType == CodonVariationType.MODIFICATION; }
    }

    private static final class CodonTracker {
        private final byte[] refSeq;
        private final List<Interval> exonList;
        private final long[][] codonCounts;
        private final int[] refCodonValues;

        public CodonTracker( final String orfCoords, final byte[] refSeq ) {
            this.refSeq = refSeq;
            exonList = getExons(orfCoords);

            codonCounts = new long[exonList.stream().mapToInt(Interval::size).sum() / 3][];
            for ( int codonId = 0; codonId != codonCounts.length; ++codonId ) {
                codonCounts[codonId] = new long[CODON_COUNT_ROW_SIZE];
            }

            refCodonValues = parseReferenceIntoCodons(refSeq, exonList);
        }

        public int[] getRefCodonValues() { return refCodonValues; }
        public long[][] getCodonCounts() { return codonCounts; }

        public List<CodonVariation> encodeSNVsAsCodons( final List<SNV> snvs ) {
            final List<CodonVariation> codonVariations = new ArrayList<>();
            final Iterator<SNV> snvIterator = snvs.iterator();
            SNV snv = null;
            while ( snvIterator.hasNext() ) {
                final SNV testSNV = snvIterator.next();
                if ( isExonic(testSNV.getRefIndex()) ) {
                    snv = testSNV;
                    break;
                }
            }

            while ( snv != null ) {
                int refIndex = snv.getRefIndex();
                final Iterator<Interval> exonIterator = exonList.iterator();
                Interval curExon;
                do {
                    curExon = exonIterator.next();
                } while ( curExon.getEnd() < refIndex );
                if ( curExon.getStart() > refIndex ) {
                    throw new GATKException("can't find current exon, even though refIndex should be exonic.");
                }

                int codonId = exonicBaseCount(refIndex);
                int codonPhase = codonId % 3;
                codonId /= 3;

                int codonValue = refCodonValues[codonId];
                if ( codonPhase == 0 ) codonValue = 0;
                else if ( codonPhase == 1 ) codonValue >>= 4;
                else codonValue >>= 2;

                int leadLag = 0;
                do {
                    boolean codonValueAltered = false;
                    boolean bumpRefIndex = false;
                    if ( snv == null || snv.getRefIndex() != refIndex ) {
                        codonValue = (codonValue << 2) | "ACGT".indexOf(refSeq[refIndex]);
                        codonValueAltered = true;
                        bumpRefIndex = true;
                    } else {
                        if ( snv.getVariantCall() == NO_CALL ) {
                            if ( --leadLag == -3 ) {
                                codonVariations.add(
                                        new CodonVariation(codonId, -1, CodonVariationType.DELETION));
                                if ( ++codonId == refCodonValues.length ) {
                                    return codonVariations;
                                }
                                leadLag = 0;
                            }
                            bumpRefIndex = true;
                        } else if ( snv.getRefCall() == NO_CALL ) {
                            leadLag += 1;
                            codonValue = (codonValue << 2) | "ACGT".indexOf(snv.getVariantCall());
                            codonValueAltered = true;
                        } else {
                            codonValue = (codonValue << 2) | "ACGT".indexOf(snv.getVariantCall());
                            codonValueAltered = true;
                            bumpRefIndex = true;
                        }
                        snv = null;
                        while ( snvIterator.hasNext() ) {
                            final SNV testSNV = snvIterator.next();
                            if ( isExonic(testSNV.getRefIndex()) ) {
                                snv = testSNV;
                                break;
                            }
                        }
                    }
                    if ( bumpRefIndex ) {
                        if ( ++refIndex == curExon.getEnd() ) {
                            curExon = exonIterator.next();
                            if ( curExon.getStart() != Integer.MAX_VALUE ) {
                                refIndex = curExon.getStart();
                            }
                        }
                        if ( refIndex == refSeq.length ) {
                            return codonVariations;
                        }
                    }
                    if ( codonValueAltered ) {
                        if ( ++codonPhase == 3 ) {
                            if ( leadLag == 3 ) {
                                codonVariations.add(
                                        new CodonVariation(codonId, codonValue, CodonVariationType.INSERTION));
                                leadLag = 0;
                                codonId -= 1;
                            } else if ( codonValue != refCodonValues[codonId] ) {
                                codonVariations.add(
                                        new CodonVariation(codonId, codonValue, CodonVariationType.MODIFICATION));
                            }
                            if ( isStop(codonValue) ) {
                                return codonVariations;
                            }
                            if ( ++codonId == refCodonValues.length ) {
                                return codonVariations;
                            }
                            codonPhase = 0;
                            codonValue = 0;
                        }
                    }
                } while ( leadLag != 0 || codonPhase != 0 );
            }

            return codonVariations;
        }

        public void reportVariantCodonCounts( final Interval refCoverage,
                                              final List<CodonVariation> variantCodons ) {
            final int startingCodonId = (exonicBaseCount(refCoverage.getStart()) + 2) / 3;
            final int endingCodonId = exonicBaseCount(refCoverage.getEnd()) / 3;
            final Iterator<CodonVariation> variantCodonIterator = variantCodons.iterator();
            CodonVariation codonVariation = variantCodonIterator.hasNext() ? variantCodonIterator.next() : null;
            for ( int codonId = startingCodonId; codonId < endingCodonId; ++codonId ) {
                while ( codonVariation != null && codonVariation.getCodonId() < codonId ) {
                    codonVariation = variantCodonIterator.hasNext() ? variantCodonIterator.next() : null;
                }
                if ( codonVariation == null || codonVariation.getCodonId() != codonId ) {
                    codonCounts[codonId][refCodonValues[codonId]] += 1;
                } else {
                    final int NO_INDEL = -1;
                    int indelColumn = NO_INDEL;
                    do {
                        switch ( codonVariation.getVariationType() ) {
                            case FRAMESHIFT:
                                indelColumn = FRAME_SHIFTING_INDEL_INDEX;
                                break;
                            case DELETION:
                            case INSERTION:
                                if ( indelColumn == NO_INDEL ) indelColumn = FRAME_PRESERVING_INDEL_INDEX;
                                break;
                            case MODIFICATION:
                                codonCounts[codonId][codonVariation.getCodonValue()] += 1;
                                break;
                        }
                        codonVariation = variantCodonIterator.hasNext() ? variantCodonIterator.next() : null;
                    } while ( codonVariation != null && codonVariation.getCodonId() == codonId );
                    if ( indelColumn != NO_INDEL ) {
                        codonCounts[codonId][indelColumn] += 1;
                    }
                }
            }
        }

        public void reportWildCodonCounts( final Interval refCoverage ) {
            final int startingCodonId = (exonicBaseCount(refCoverage.getStart()) + 2) / 3;
            final int endingCodonId = exonicBaseCount(refCoverage.getEnd()) / 3;
            for ( int codonId = startingCodonId; codonId < endingCodonId; ++codonId ) {
                codonCounts[codonId][refCodonValues[codonId]] += 1;
            }
        }

        private static List<Interval> getExons( final String orfCoords ) {
            final List<Interval> exonList = new ArrayList<>();
            for ( final String coordPair : orfCoords.split(",") ) {
                final String[] coords = coordPair.split("-");
                if ( coords.length != 2 ) {
                    throw new UserException("Can't interpret ORF as list of pairs of coords: " + orfCoords);
                }
                try {
                    final int start = Integer.valueOf(coords[0]);
                    if ( start < 1 ) {
                        throw new UserException("Coordinates of ORF are 1-based.");
                    }
                    final int end = Integer.valueOf(coords[1]);
                    if ( end < start ) {
                        throw new UserException("Found ORF end coordinate less than start: " + orfCoords);
                    }
                    // convert 1-based, inclusive intervals to 0-based, half-open
                    final Interval exon = new Interval(start - 1, end);
                    exonList.add(exon);
                } catch ( final NumberFormatException nfe ) {
                    throw new UserException("Can't interpret ORF coords as integers: " + orfCoords);
                }
                for ( int idx = 1; idx < exonList.size(); ++idx ) {
                    if ( exonList.get(idx - 1).getEnd() >= exonList.get(idx).getStart() ) {
                        throw new UserException("ORF coordinates are not sorted: " + orfCoords);
                    }
                }
            }

            // it's helpful to have this 0-length sentinel at the end of the list
            exonList.add(new Interval(Integer.MAX_VALUE, Integer.MAX_VALUE));

            final int orfLen = exonList.stream().mapToInt(Interval::size).sum();
            if ( (orfLen % 3) != 0 ) {
                throw new UserException("ORF length must be divisible by 3.");
            }

            return exonList;
        }

        private static int[] parseReferenceIntoCodons( final byte[] refSeq, final List<Interval> exonList ) {
            final int nCodons = exonList.stream().mapToInt(Interval::size).sum() / 3;
            final int[] refCodonValues = new int[nCodons];
            int codonId = 0;
            int codonPhase = 0;
            int codonValue = 0;
            for ( final Interval exon : exonList ) {
                final int exonEnd = exon.getEnd();
                for ( int refIndex = exon.getStart(); refIndex != exonEnd; ++refIndex ) {
                    codonValue = (codonValue << 2) | "ACGT".indexOf(refSeq[refIndex]);
                    if ( ++codonPhase == 3 ) {
                        if ( isStop(codonValue) && codonId != nCodons - 1 ) {
                            final int idx = refIndex + 1;
                            throw new UserException("There is an upstream stop codon at reference index " + idx + ".");
                        }
                        refCodonValues[codonId] = codonValue;
                        codonValue = 0;
                        codonPhase = 0;
                        codonId += 1;
                    }
                }
            }

            final int START_CODON = 0x0E;
            if ( refCodonValues[0] != START_CODON ) {
                System.err.println("WARNING:  Your ORF does not start with the expected ATG codon.");
            }
            final int lastCodon = refCodonValues[nCodons - 1];
            if ( !isStop(lastCodon) ) {
                System.err.println("WARNING:  Your ORF does not end with the expected stop codon.");
            }

            return refCodonValues;
        }

        public static boolean isStop( final int codonValue ) {
            final int STOP_OCH = 0x30;
            final int STOP_AMB = 0x32;
            final int STOP_OPA = 0x38;
            return codonValue == STOP_OCH || codonValue == STOP_AMB || codonValue == STOP_OPA;
        }

        private boolean isExonic( final int refIndex ) {
            for ( final Interval exon : exonList ) {
                if ( exon.getStart() > refIndex ) return false;
                if ( exon.getEnd() > refIndex ) return true;
            }
            throw new GATKException("can't reach here");
        }

        private int exonicBaseCount( final int refIndex ) {
            int baseCount = 0;
            for ( final Interval exon : exonList ) {
                if ( refIndex >= exon.getEnd() ) {
                    baseCount += exon.size();
                } else {
                    if ( refIndex > exon.getStart() ) {
                        baseCount += refIndex - exon.getStart();
                    }
                    break;
                }
            }
            return baseCount;
        }
    }

    private final static class ReadReport {
        final List<Interval> refCoverage;
        final List<SNV> snvList;

        public ReadReport( final List<Interval> refCoverage,
                           final List<SNV> snvList ) {
            this.refCoverage = refCoverage;
            this.snvList = snvList;
        }

        public List<Interval> getRefCoverage() {
            return refCoverage;
        }
        public List<SNV> getVariations() {
            return snvList;
        }

        public static ReadReport nullReport = new ReadReport(new ArrayList<>(), new ArrayList<>());
    }

    private void initializeRefSeq() {
        final ReferenceDataSource reference = ReferenceDataSource.of(referenceArguments.getReferencePath());
        final SAMSequenceDictionary seqDict = reference.getSequenceDictionary();
        if ( seqDict.size() != 1 ) {
            throw new UserException("Expecting a reference with a single contig. " +
                    "The supplied reference has " + seqDict.size() + " contigs.");
        }
        final SAMSequenceRecord tig0 = seqDict.getSequence(0);
        final int refSeqLen = tig0.getSequenceLength();
        final SimpleInterval wholeTig = new SimpleInterval(tig0.getSequenceName(), 1, refSeqLen);
        refSeq = Arrays.copyOf(reference.queryAndPrefetch(wholeTig).getBases(), refSeqLen);
        for ( int idx = 0; idx < refSeqLen; ++idx ) {
            switch ( refSeq[idx] &= UPPERCASE_MASK ) { // make into upper case
                case 'A':
                case 'C':
                case 'G':
                case 'T':
                    break;
                default:
                    throw new UserException("Reference sequence contains something other than A, C, G, and T.");
            }
        }
        refCoverage = new long[refSeq.length];
        refCoverageSizeHistogram = new long[refSeq.length + 1];
    }

    private ReadReport getReadReport( final GATKRead read ) {
        nReadsTotal += 1;
        nTotalBaseCalls += read.getLength();

        if ( read.isUnmapped() ) {
            nReadsUnmapped += 1;
            return ReadReport.nullReport;
        }

        final Interval trim = calculateTrim(read.getBaseQualitiesNoCopy());
        if ( trim.size() == 0 ) {
            return ReadReport.nullReport;
        }
        return analyzeAlignment(read, trim.getStart(), trim.getEnd());
    }

    private Interval calculateTrim( final byte[] quals ) {
        // find initial end-trim
        int readStart = 0;
        int hiQCount = 0;
        while ( readStart < quals.length ) {
            if ( quals[readStart] < minQ ) {
                hiQCount = 0;
            } else if ( ++hiQCount == minLength ) {
                break;
            }
            readStart += 1;
        }
        if ( readStart == quals.length ) {
            nReadsLowQuality += 1;
            return Interval.nullInterval;
        }
        readStart -= minLength - 1;

        // find final end-trim
        int readEnd = quals.length - 1;
        hiQCount = 0;
        while ( readEnd >= 0 ) {
            if ( quals[readEnd] < minQ ) {
                hiQCount = 0;
            } else if ( ++hiQCount == minLength ) {
                break;
            }
            readEnd -= 1;
        }
        readEnd += minLength;

        return new Interval(readStart, readEnd);
    }

    private ReadReport analyzeAlignment( final GATKRead read, final int start, final int end ) {
        final Cigar cigar = read.getCigar();
        final Iterator<CigarElement> cigarIterator = cigar.getCigarElements().iterator();
        CigarElement cigarElement = cigarIterator.next();
        CigarOperator cigarOperator = cigarElement.getOperator();
        int cigarElementCount = cigarElement.getLength();

        final byte[] readSeq = read.getBasesNoCopy();
        final byte[] readQuals = read.getBaseQualitiesNoCopy();

        final List<SNV> variations = new ArrayList<>();

        int refIndex = read.getStart() - 1; // 0-based numbering
        int readIndex = 0;

        // pretend that soft-clips are matches
        if ( cigarOperator == CigarOperator.S ) {
            refIndex -= cigarElementCount;
        }

        final List<Interval> refCoverageList = new ArrayList<>();
        int refCoverageBegin = -1;
        int refCoverageEnd = -1;

        while ( true ) {
            if ( readIndex >= start && refIndex >= 0 ) {
                if ( refCoverageBegin == -1 ) {
                    refCoverageBegin = refIndex;
                    refCoverageEnd = refIndex;
                }
                if ( cigarOperator == CigarOperator.D ) {
                    variations.add(new SNV(refIndex, refSeq[refIndex], NO_CALL, readQuals[readIndex]));
                } else if ( cigarOperator == CigarOperator.I ) {
                    final byte call = (byte)(readSeq[readIndex] & UPPERCASE_MASK);
                    variations.add(new SNV(refIndex, NO_CALL, call, readQuals[readIndex]));
                } else if ( cigarOperator == CigarOperator.M || cigarOperator == CigarOperator.S ) {
                    final byte call = (byte)(readSeq[readIndex] & UPPERCASE_MASK);
                    if ( call != refSeq[refIndex] ) {
                        variations.add(new SNV(refIndex, refSeq[refIndex], call, readQuals[readIndex]));
                    }
                    if ( refIndex == refCoverageEnd ) {
                        refCoverageEnd += 1;
                    } else {
                        refCoverageList.add(new Interval(refCoverageBegin, refCoverageEnd));
                        refCoverageBegin = refIndex;
                        refCoverageEnd = refIndex + 1;
                    }
                } else {
                    throw new GATKException("unanticipated cigar operator: " + cigarOperator.toString());
                }
            }

            if ( cigarOperator != CigarOperator.D ) {
                if ( ++readIndex == end ) {
                    break;
                }
            }

            if ( cigarOperator != CigarOperator.I ) {
                if ( ++refIndex == refSeq.length )
                    break;
            }

            if ( --cigarElementCount == 0 ) {
                if ( !cigarIterator.hasNext() ) {
                    throw new GATKException("unexpectedly exhausted cigar iterator");
                }
                cigarElement = cigarIterator.next();
                cigarOperator = cigarElement.getOperator();
                cigarElementCount = cigarElement.getLength();
            }
        }

        if ( refCoverageBegin < refCoverageEnd ) {
            refCoverageList.add(new Interval(refCoverageBegin, refCoverageEnd));
        }

        return new ReadReport(refCoverageList, variations);
    }

    private static ReadReport combineReports( final ReadReport report1, final ReadReport report2 ) {
        final List<Interval> refCoverage1 = report1.getRefCoverage();
        final List<Interval> refCoverage2 = report2.getRefCoverage();
        if ( refCoverage1.isEmpty() ) {
            return report2;
        } else if ( refCoverage2.isEmpty() ) {
            return report1;
        }

        final int overlapStart = Math.max(refCoverage1.get(0).getStart(), refCoverage2.get(0).getStart());
        final int overlapEnd = Math.min(refCoverage1.get(refCoverage1.size() - 1).getEnd(),
                                        refCoverage2.get(refCoverage2.size() - 1).getEnd());
        if ( overlapEnd < overlapStart ) {
            return null; // disjoint reports can't be combined -- caller will apply reports separately
        }

        final List<Interval> combinedCoverage = combineIntervals(refCoverage1, refCoverage2);
        final List<SNV> combinedSNVs =
                combineSNVs(report1.getVariations(), report2.getVariations(), overlapStart, overlapEnd);
        return new ReadReport(combinedCoverage, combinedSNVs);
    }

    private static List<Interval> combineIntervals( final List<Interval> refCoverage1,
                                                    final List<Interval> refCoverage2 ) {
        final List<Interval> combinedCoverage = new ArrayList<>();
        final Iterator<Interval> refCoverageItr1 = refCoverage1.iterator();
        final Iterator<Interval> refCoverageItr2 = refCoverage2.iterator();
        Interval refInterval1 = refCoverageItr1.next();
        Interval refInterval2 = refCoverageItr2.next();
        Interval curInterval;
        if ( refInterval1.getStart() < refInterval2.getStart() ) {
            curInterval = refInterval1;
            refInterval1 = refCoverageItr1.hasNext() ? refCoverageItr1.next() : null;
        } else {
            curInterval = refInterval2;
            refInterval2 = refCoverageItr2.hasNext() ? refCoverageItr2.next() : null;
        }
        while ( refInterval1 != null || refInterval2 != null ) {
            final Interval testInterval;
            if ( refInterval1 == null ) {
                testInterval = refInterval2;
                refInterval2 = refCoverageItr2.hasNext() ? refCoverageItr2.next() : null;
            } else if ( refInterval2 == null ) {
                testInterval = refInterval1;
                refInterval1 = refCoverageItr1.hasNext() ? refCoverageItr1.next() : null;
            } else if ( refInterval1.getStart() < refInterval2.getStart() ) {
                testInterval = refInterval1;
                refInterval1 = refCoverageItr1.hasNext() ? refCoverageItr1.next() : null;
            } else {
                testInterval = refInterval2;
                refInterval2 = refCoverageItr2.hasNext() ? refCoverageItr2.next() : null;
            }
            if ( curInterval.getEnd() < testInterval.getStart() ) {
                combinedCoverage.add(curInterval);
                curInterval = testInterval;
            } else {
                curInterval =
                        new Interval(curInterval.getStart(), Math.max(curInterval.getEnd(), testInterval.getEnd()));
            }
        }
        combinedCoverage.add(curInterval);

        return combinedCoverage;
    }

    private static List<SNV> combineSNVs( final List<SNV> snvs1, final List<SNV> snvs2,
                                   final int overlapStart, final int overlapEnd ) {
        final List<SNV> combinedSNVs = new ArrayList<>();
        final Iterator<SNV> iterator1 = snvs1.iterator();
        final Iterator<SNV> iterator2 = snvs2.iterator();
        SNV snv1 = iterator1.hasNext() ? iterator1.next() : null;
        SNV snv2 = iterator2.hasNext() ? iterator2.next() : null;
        while ( snv1 != null || snv2 != null ) {
            final SNV next;
            if ( snv1 == null ) {
                next = snv2;
                snv2 = iterator2.hasNext() ? iterator2.next() : null;
                final int refIndex = next.getRefIndex();
                if ( refIndex >= overlapStart && refIndex < overlapEnd ) {
                    return null;
                }
            } else if ( snv2 == null ) {
                next = snv1;
                snv1 = iterator1.hasNext() ? iterator1.next() : null;
                final int refIndex = next.getRefIndex();
                if ( refIndex >= overlapStart && refIndex < overlapEnd ) {
                    return null;
                }
            } else {
                final int refIndex1 = snv1.getRefIndex();
                final int refIndex2 = snv2.getRefIndex();
                if ( refIndex1 < refIndex2 ) {
                    next = snv1;
                    snv1 = iterator1.hasNext() ? iterator1.next() : null;
                    if ( refIndex1 >= overlapStart && refIndex1 < overlapEnd ) {
                        return null;
                    }
                } else if ( refIndex2 < refIndex1 ) {
                    next = snv2;
                    snv2 = iterator2.hasNext() ? iterator2.next() : null;
                    if ( refIndex2 >= overlapStart && refIndex2 < overlapEnd ) {
                        return null;
                    }
                } else if ( !snv1.equals(snv2) ) {
                    return null;
                } else {
                    next = snv1.getQuality() > snv2.getQuality() ? snv1 : snv2;
                    snv1 = iterator1.hasNext() ? iterator1.next() : null;
                    snv2 = iterator2.hasNext() ? iterator2.next() : null;
                }
            }
            combinedSNVs.add(next);
        }
        return combinedSNVs;
    }

    private void applyReport( final ReadReport readReport ) {
        final List<Interval> refCoverageList = readReport.getRefCoverage();
        if ( refCoverageList.isEmpty() ) {
            return;
        }

        int coverage = 0;
        for ( final Interval refInterval : refCoverageList ) {
            final int refIntervalStart = refInterval.getStart();
            final int refIntervalEnd = refInterval.getEnd();
            coverage += refIntervalEnd - refIntervalStart;
            for ( int idx = refInterval.getStart(); idx != refIntervalEnd; ++idx ) {
                refCoverage[idx] += 1;
            }
        }

        refCoverageSizeHistogram[coverage] += 1;

        final int refStart = refCoverageList.get(0).getStart();
        final int refEnd = refCoverageList.get(refCoverageList.size() - 1).getEnd();
        if ( intervalCounter == null ) {
            intervalCounter = new IntervalCounter(refSeq.length);
        }
        intervalCounter.addCount(refStart, refEnd);

        final List<SNV> variations = readReport.getVariations();
        if ( variations == null ) {
            nInconsistentPairs += 1;
        } else if ( variations.isEmpty() ) {
            nWildTypeMolecules += 1;
            if ( readReport.getRefCoverage().size() != 1 ) {
                throw new GATKException("expecting a single coverage interval for a wild-type molecule");
            }
            codonTracker.reportWildCodonCounts(readReport.getRefCoverage().get(0));
        } else if ( variations.stream().anyMatch(
                snv -> snv.getQuality() < minQ || "-ACGT".indexOf(snv.getVariantCall()) == -1) ) {
            nLowQualityVariantMolecules += 1;
        } else if ( refEnd - variations.get(variations.size() - 1).getRefIndex() < minFlankingLength ||
                variations.get(0).getRefIndex() - refCoverageList.get(0).getStart() < minFlankingLength ) {
            nInsufficientFlankMolecules += 1;
        } else {
            nCalledVariantMolecules += 1;
            final List<Interval> refCoverage = readReport.getRefCoverage();
            final int refCoverageSize = refCoverage.size();
            final Interval totalCoverage = refCoverageSize == 1 ? refCoverage.get(0) :
                        new Interval(refCoverage.get(0).getStart(), refCoverage.get(refCoverageSize - 1).getEnd());
            codonTracker.reportVariantCodonCounts(totalCoverage, codonTracker.encodeSNVsAsCodons(variations));

            final SNVCollectionCount newVal = new SNVCollectionCount(variations, coverage);
            final SNVCollectionCount oldVal = variationCounts.find(newVal);
            if ( oldVal != null ) {
                oldVal.bumpCount(coverage);
            } else {
                variationCounts.add(newVal);
            }
        }
    }
}
