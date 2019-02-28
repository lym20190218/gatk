package org.broadinstitute.hellbender.utils.codecs.copynumber;

import htsjdk.tribble.AsciiFeatureCodec;
import htsjdk.tribble.readers.LineIterator;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.tools.copynumber.formats.CopyNumberFormatsUtils;
import org.broadinstitute.hellbender.tools.copynumber.formats.records.SimpleCount;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.tsv.TableUtils;

import java.util.Arrays;
import java.util.List;


public final class SimpleCountCodec extends AsciiFeatureCodec<SimpleCount> {

    private static final String COLUMN_HEADER_STRING = String.format("CONTIG%1$sSTART%1$sEND%1$sCOUNT", TableUtils.COLUMN_SEPARATOR_STRING);
    private static final List<String> COLUMN_HEADER = Arrays.asList(COLUMN_HEADER_STRING.split(TableUtils.COLUMN_SEPARATOR_STRING));

    private boolean havePassedHeader = false;

    public SimpleCountCodec() {
        super(SimpleCount.class);
    }

    @Override
    public SimpleCount decode(final String line) {
        if (line.startsWith(CopyNumberFormatsUtils.COMMENT_PREFIX) || line.startsWith(COLUMN_HEADER_STRING)) {
            if (!havePassedHeader) {
                return null;
            } else {
                throw new UserException.MalformedFile("Header lines must be at the beginning of the file.");
            }
        } else {
            havePassedHeader = true;
            final String[] split = line.split(TableUtils.COLUMN_SEPARATOR_STRING);
            try {
                return new SimpleCount(new SimpleInterval(split[0], Integer.parseInt(split[1]), Integer.parseInt(split[2])), Integer.parseInt(split[3]));
            } catch (final NumberFormatException e) {
                throw new UserException.MalformedFile("Line = " + line + " is not formatted correctly.");
            }
        }
    }

    @Override
    public List<String> readActualHeader(final LineIterator reader) {
        //we check that the SAM header lines and the column header line are present in the correct order, then return the mandatory column header
        boolean isSAMHeaderPresent = false;
        while (reader.hasNext()) {
            final String line = reader.peek();
            if (line.startsWith(CopyNumberFormatsUtils.COMMENT_PREFIX)) {
                isSAMHeaderPresent = true;
                reader.next();
            } else {
                if (!isSAMHeaderPresent) {
                    throw new UserException.MalformedFile("SAM header lines must be at the beginning of the file.");
                } else if (!line.startsWith(COLUMN_HEADER_STRING)) {
                    throw new UserException.MalformedFile("File does not have a column header.");
                } else {
                    break;
                }
            }
        }
        return COLUMN_HEADER;
    }

    @Override
    public boolean canDecode(final String path) {
        return path.toLowerCase().endsWith(".counts.tsv");
    }
}

