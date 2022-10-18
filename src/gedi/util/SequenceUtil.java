package gedi.util;

/*
 * The MIT License
 *
 * Copyright (c) 2009 The Broad Institute
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in
 * all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
 * THE SOFTWARE.
 */

import htsjdk.samtools.Cigar;
import htsjdk.samtools.CigarElement;
import htsjdk.samtools.CigarOperator;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMTag;
import java.util.Arrays;
import java.util.List;


/**
 * Class SequenceUtil from https://github.com/samtools/htsjdk/tree/master/src/main/java/htsjdk/samtools
 * and slightly modified
 */
public class SequenceUtil {
    /** Byte typed variables for all normal bases. */
    public static final byte a = 'a', c = 'c', g = 'g', t = 't', n = 'n', A = 'A', C = 'C', G = 'G', T = 'T', N = 'N';


    /**
     * A set of bases supported by BAM in reads, see http://samtools.github.io/hts-specs/SAMv1.pdf chapter 4.2 on 'seq' field.
     * Effectively these are upper cased IUPAC codes with equals sign ('=') and without dot ('.').
     */
    private static final byte[] BAM_READ_BASE_SET = "=ABCDGHKMNRSTVWY".getBytes();

    private static final int BASES_ARRAY_LENGTH = 127;
    private static final int SHIFT_TO_LOWER_CASE = a - A;
    /**
     * A lookup table to find a corresponding BAM read base.
     */
    private static final byte[] bamReadBaseLookup = new byte[BASES_ARRAY_LENGTH];
    static {
        Arrays.fill(bamReadBaseLookup, N);
        for (final byte base: BAM_READ_BASE_SET) {
            bamReadBaseLookup[base] = base;
            bamReadBaseLookup[base + SHIFT_TO_LOWER_CASE] = base;
        }
    }

    private static final byte A_MASK = 1;
    private static final byte C_MASK = 2;
    private static final byte G_MASK = 4;
    private static final byte T_MASK = 8;

    private static final byte[] bases = new byte[BASES_ARRAY_LENGTH];
    private static final byte NON_IUPAC_CODE = 0;
    /*
     * Definition of IUPAC codes:
     * http://www.bioinformatics.org/sms2/iupac.html
     */
    static {
        Arrays.fill(bases, NON_IUPAC_CODE);
        bases[A] = A_MASK;
        bases[C] = C_MASK;
        bases[G] = G_MASK;
        bases[T] = T_MASK;
        bases['M'] = A_MASK | C_MASK;
        bases['R'] = A_MASK | G_MASK;
        bases['W'] = A_MASK | T_MASK;
        bases['S'] = C_MASK | G_MASK;
        bases['Y'] = C_MASK | T_MASK;
        bases['K'] = G_MASK | T_MASK;
        bases['V'] = A_MASK | C_MASK | G_MASK;
        bases['H'] = A_MASK | C_MASK | T_MASK;
        bases['D'] = A_MASK | G_MASK | T_MASK;
        bases['B'] = C_MASK | G_MASK | T_MASK;
        bases['N'] = A_MASK | C_MASK | G_MASK | T_MASK;
        // Also store the bases in lower case
        for (int i = 'A'; i <= 'Z'; i++) {
            bases[(byte) i + SHIFT_TO_LOWER_CASE] = bases[(byte) i];
        }
        bases['.'] = A_MASK | C_MASK | G_MASK | T_MASK;
    };


    /**
     * Calculate the reverse complement of the specified sequence
     * (Stolen from Reseq)
     *
     * @param sequenceData
     * @return reverse complement
     */
    public static String reverseComplement(final String sequenceData) {
        final byte[] bases = htsjdk.samtools.util.StringUtil.stringToBytes(sequenceData);
        reverseComplement(bases);
        return htsjdk.samtools.util.StringUtil.bytesToString(bases);
    }


    /** Returns the complement of a single byte. */
    public static byte complement(final byte b) {
        switch (b) {
            case a:
                return t;
            case c:
                return g;
            case g:
                return c;
            case t:
                return a;
            case A:
                return T;
            case C:
                return G;
            case G:
                return C;
            case T:
                return A;
            default:
                return b;
        }
    }

    /** Reverses and complements the bases in place. */
    public static void reverseComplement(final byte[] bases) {
        reverseComplement(bases, 0, bases.length);
    }



    public static void reverseComplement(final byte[] bases, final int offset, final int len) {
        final int lastIndex = len - 1;

        int i, j;
        for (i = offset, j = offset + lastIndex; i < j; ++i, --j) {
            final byte tmp = complement(bases[i]);
            bases[i] = complement(bases[j]);
            bases[j] = tmp;
        }
        if (len % 2 == 1) {
            bases[i] = complement(bases[i]);
        }
    }


    /**
     * Calculate MD and NM similarly to Samtools, except that N->N is a match.
     *
     * @param record Input record for which to calculate NM and MD.
     *               The appropriate tags will be added/updated in the record
     * @param ref    The reference bases for the sequence to which the record is mapped
     * @param calcMD A flag indicating whether to update the MD tag in the record
     * @param calcNM A flag indicating whether to update the NM tag in the record
     */
    public static void calculateMdAndNmTags(final SAMRecord record, final byte[] ref,
                                            final boolean calcMD, final boolean calcNM) {
        if (!calcMD && !calcNM)
            return;

        final Cigar cigar = record.getCigar();
        final List<CigarElement> cigarElements = cigar.getCigarElements();
        final byte[] seq = record.getReadBases();
        final int alignmentStart = record.getAlignmentStart() - 1;
        int cigarIndex, blockRefPos, blockReadStart, matchCount = 0;
        int nmCount = 0;
        final StringBuilder mdString = new StringBuilder();

        final int nElements = cigarElements.size();

        for (cigarIndex = blockReadStart = 0, blockRefPos = alignmentStart; cigarIndex < nElements; ++cigarIndex) {
            final CigarElement ce = cigarElements.get(cigarIndex);
            int inBlockOffset;
            final int blockLength = ce.getLength();
            final CigarOperator op = ce.getOperator();
            if (op == CigarOperator.MATCH_OR_MISMATCH || op == CigarOperator.EQ
                    || op == CigarOperator.X) {
                for (inBlockOffset = 0; inBlockOffset < blockLength; ++inBlockOffset) {
                    System.out.println("XXXXXXXXXXXXXXXXXX");

                    final int readOffset = blockReadStart + inBlockOffset;

                    if (ref.length <= blockRefPos + inBlockOffset) break; // out of boundary

                    final byte readBase = seq[readOffset];
                    final byte refBase = ref[blockRefPos + inBlockOffset];

                    if ((bases[readBase] == bases[refBase]) || readBase == 0) {
                        // a match
                        ++matchCount;
                    } else {
                        mdString.append(matchCount);
                        mdString.appendCodePoint(refBase);
                        matchCount = 0;
                        ++nmCount;
                    }
                }
                if (inBlockOffset < blockLength) break;
                blockRefPos += blockLength;
                blockReadStart += blockLength;
            } else if (op == CigarOperator.DELETION) {
                mdString.append(matchCount);
                mdString.append('^');
                for (inBlockOffset = 0; inBlockOffset < blockLength; ++inBlockOffset) {
                    if (ref[blockRefPos + inBlockOffset] == 0) break;
                    mdString.appendCodePoint(ref[blockRefPos + inBlockOffset]);
                }
                matchCount = 0;
                if (inBlockOffset < blockLength) break;
                blockRefPos += blockLength;
                nmCount += blockLength;
            } else if (op == CigarOperator.INSERTION
                    || op == CigarOperator.SOFT_CLIP) {
                blockReadStart += blockLength;
                if (op == CigarOperator.INSERTION) nmCount += blockLength;
            } else if (op == CigarOperator.SKIPPED_REGION) {
                blockRefPos += blockLength;
            }
        }
        mdString.append(matchCount);

        if (calcMD) record.setAttribute(SAMTag.MD.name(), mdString.toString());
        if (calcNM) record.setAttribute(SAMTag.NM.name(), nmCount);
    }

    public static void calculateMdAndNmTags(final SAMRecord record, final byte[] ref, final int refOffset, final boolean calcMD, final boolean calcNM) {
        if (!calcMD && !calcNM)
            return;

        final Cigar cigar = record.getCigar();
        final List<CigarElement> cigarElements = cigar.getCigarElements();
        final byte[] seq = record.getReadBases();
        final int start = record.getAlignmentStart() - 1;
        int i, x, y, u = 0;
        int nm = 0;
        final StringBuilder str = new StringBuilder();

        final int size = cigarElements.size();
        for (i = y = 0, x = start; i < size; ++i) {
            final CigarElement ce = cigarElements.get(i);
            int j;
            final int length = ce.getLength();
            final CigarOperator op = ce.getOperator();
            if (op == CigarOperator.MATCH_OR_MISMATCH || op == CigarOperator.EQ || op == CigarOperator.X) {
                for (j = 0; j < length; ++j) {
                    final int z = y + j;

                    if (refOffset + ref.length <= x + j)
                        break; // out of boundary

                    int c1 = 0;
                    int c2 = 0;
                    // try {
                    c1 = seq[z];
                    c2 = ref[x + j - refOffset];

                    if ((c1 == c2) || c1 == 0) {
                        // a match
                        ++u;
                    } else {
                        str.append(u);
                        str.appendCodePoint(ref[x + j - refOffset]);
                        u = 0;
                        ++nm;
                    }
                }
                if (j < length)
                    break;
                x += length;
                y += length;
            } else if (op == CigarOperator.DELETION) {
                str.append(u);
                str.append('^');
                for (j = 0; j < length; ++j) {
                    if (ref[x + j - refOffset] == 0)
                        break;
                    str.appendCodePoint(ref[x + j - refOffset]);
                }
                u = 0;
                if (j < length)
                    break;
                x += length;
                nm += length;
            } else if (op == CigarOperator.INSERTION || op == CigarOperator.SOFT_CLIP) {
                y += length;
                if (op == CigarOperator.INSERTION)
                    nm += length;
            } else if (op == CigarOperator.SKIPPED_REGION) {
                x += length;
            }
        }
        str.append(u);

        if (calcMD)
            record.setAttribute(SAMTag.MD.name(), str.toString());
        if (calcNM)
            record.setAttribute("nM", nm);
    }
}
