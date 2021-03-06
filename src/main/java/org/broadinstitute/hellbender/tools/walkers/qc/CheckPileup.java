package org.broadinstitute.hellbender.tools.walkers.qc;

import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.barclay.help.DocumentedFeature;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import picard.cmdline.programgroups.DiagnosticsAndQCProgramGroup;
import org.broadinstitute.hellbender.engine.*;
import org.broadinstitute.hellbender.engine.filters.ReadFilter;
import org.broadinstitute.hellbender.engine.filters.ReadFilterLibrary;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.utils.IntervalUtils;
import org.broadinstitute.hellbender.utils.codecs.sampileup.SAMPileupFeature;
import org.broadinstitute.hellbender.utils.pileup.ReadPileup;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.PrintStream;
import java.util.List;

/**
 * Compare GATK's internal pileup to a reference Samtools pileup
 *
 * <p>
 *     This tool compares the mpileup data (reference base, aligned base from each overlapping read, and quality score)
 *     generated internally by GATK to a reference pileup data generated by Samtools, for each position in the requested
 *     interval. Note that this only considers the single-sample mpileup format. See the <a href="http://samtools.sourceforge.net/pileup.shtml">Pileup format documentation</a>
 *     for more details on the format.
 * </p>
 *
 * <h3>Inputs</h3>
 * <ul>
 *     <li>A BAM or CRAM file containing aligned read data.</li>
 *     <li>An mpileup file generated by Samtools covering the interval of interest.
 *         It should be accompanied by an index, e.g. that generated with IndexFeatureFile.</li>
 * </ul>
 *
 * <h3>Output</h3>
 * <p>
 *     A text file listing mismatches between the input mpileup and the GATK's internal pileup. If there are no mismatches,
 *     the output file is empty.
 * </p>
 *
 * <h3>Usage example</h3>
 * <pre>
 *   gatk CheckPileup \
 *     -R reference.fasta \
 *     -I input.bam \
 *     -O output.txt \
 *     --pileup samtools.pileup
 * </pre>
 *
 */
@DocumentedFeature
@CommandLineProgramProperties(
    summary = "This tool compares the mpileup data (reference base, aligned base from each overlapping read, and quality score) " +
            "generated internally by GATK to a reference pileup data generated by Samtools, for each position in the requested " +
            "interval.",
    oneLineSummary = "Compare GATK's internal pileup to a reference Samtools mpileup",
    programGroup = DiagnosticsAndQCProgramGroup.class)
public final class CheckPileup extends LocusWalker {

    /**
     * This is the existing mpileup against which we'll compare GATK's internal pileup at each genome position in the desired interval.
     */
    @Argument(
            fullName = "pileup",
            doc="Pileup generated by Samtools"
    )
    public FeatureInput<SAMPileupFeature> mpileup;

    @Argument(
            fullName = StandardArgumentDefinitions.OUTPUT_LONG_NAME,
            shortName = StandardArgumentDefinitions.OUTPUT_SHORT_NAME,
            doc = "Output file (if not provided, defaults to STDOUT)",
            optional = true
    )
    public File outFile = null;

    // This is to check for previous versions of samtools without having overlap detection
    @Argument(
            fullName = "ignore-overlaps",
            doc = "Disable read-pair overlap detection",
            optional = true
    )
    public boolean ignoreOverlaps = false;

    /**
     * By default the program will quit if it encounters an error (such as missing truth data for a given position).
     * Use this flag to override the default behavior; the program will then simply print an error message and move on
     * to the next position.
     */
    @Argument(
            fullName="continue-after-error",
            doc="Continue after encountering an error",
            optional=true
    )
    public boolean continueAfterAnError = false;

    // These are the default read filters from samtools
    @Override
    public List<ReadFilter> getDefaultReadFilters() {
        final List<ReadFilter> defaultFilters = super.getDefaultReadFilters();
        defaultFilters.add(ReadFilterLibrary.NOT_DUPLICATE);
        defaultFilters.add(ReadFilterLibrary.PASSES_VENDOR_QUALITY_CHECK);
        defaultFilters.add(ReadFilterLibrary.NOT_SECONDARY_ALIGNMENT);
        return defaultFilters;
    }

    @Override
    public boolean requiresReference() {
        return true;
    }

    private long nLoci = 0;
    private long nBases = 0;
    private PrintStream out;

    @Override
    public void onTraversalStart() {
        try {
            out = (outFile == null) ? System.out : new PrintStream(outFile);
        } catch (FileNotFoundException e) {
            throw new UserException.CouldNotCreateOutputFile(outFile, e.getMessage());
        }
    }

    @Override
    public void apply(final AlignmentContext context, final ReferenceContext ref, final FeatureContext featureContext) {
        final ReadPileup pileup = context.getBasePileup();
        final SAMPileupFeature truePileup = getTruePileup(featureContext);

        if (!ignoreOverlaps) {
            pileup.fixOverlaps();
        }

        if ( truePileup == null ) {
            out.printf("No truth pileup data available at %s%n", pileup.getPileupString((char) ref.getBase()));
            if ( !continueAfterAnError) {
                throw new UserException.BadInput(
                        String.format("No pileup data available at %s given GATK's output of %s -- this walker requires samtools mpileup data over all bases",
                        context.getLocation(), new String(pileup.getBases())));
            }
        } else {
            final String pileupDiff = pileupDiff(pileup, truePileup);
            if ( pileupDiff != null ) {
                out.printf("%s vs. %s%n", pileup.getPileupString((char) ref.getBase()), truePileup.getPileupString());
                if ( !continueAfterAnError) {
                    throw new UserException.BadInput(String.format("The input pileup doesn't match the GATK's internal pileup: %s", pileupDiff));
                }
            }
        }
        nLoci++;
        nBases += pileup.size();
    }

    public String pileupDiff(final ReadPileup a, final SAMPileupFeature b) {
        // compare sizes
        if ( a.size() != b.size() ) {
            return String.format("Sizes not equal: %s vs. %s", a.size(), b.size());
        }
        // compare locations
        if (IntervalUtils.compareLocatables(a.getLocation(), b, getReferenceDictionary()) != 0 ) {
            return String.format("Locations not equal: %s vs. %s", a.getLocation(), b);
        }
        // compare bases
        final String aBases = new String(a.getBases());
        final String bBases = b.getBasesString();
        if ( ! aBases.toUpperCase().equals(bBases.toUpperCase()) ) {
            return String.format("Bases not equal: %s vs. %s", aBases, bBases);
        }

        // compare the qualities
        final String aQuals = new String(a.getBaseQuals());
        final String bQuals = new String(b.getBaseQuals());
        if ( ! aQuals.equals(bQuals) ) {
            return String.format("Quals not equal: %s vs. %s", aQuals, bQuals);
        }
        return null;
    }

    /**
     * Extracts the true pileup data from the given mpileup.
     * @param featureContext Features spanning the current locus.
     * @return true pileup data; {@code null} if not covered
     */
    private SAMPileupFeature getTruePileup( final FeatureContext featureContext ) {
        final List<SAMPileupFeature> features = featureContext.getValues(mpileup);
        return (features.isEmpty()) ? null : features.get(0);
    }

    @Override
    public Object onTraversalSuccess() {
        return String.format("Validated %d sites covered by %d bases%n", nLoci, nBases);
    }

    @Override
    public void closeTool() {
        out.close();
    }
}
