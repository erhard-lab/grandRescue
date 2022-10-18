package gedi.util;

import gedi.core.reference.Strandness;
import gedi.util.functions.EI;
import gedi.util.functions.ExtendedIterator;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMRecordIterator;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;

import java.io.*;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.Map;
import java.util.zip.GZIPOutputStream;

public class ReadExtraction {

    /**
     * @param file .bam-File of all reads from 4sU-labelled experiment
     *             Writes 3 fastq-files: All mapped, all unmapped reads and unmapped reads where all T's are converted to C's.
     *             Saves the original read sequence (with T->C MM's from 4su labeling) in ReadID, then changes all T's to C's in Read Sequence.
     *             These reads can then be mapped to the pseudo-genome where all T's to C's are converted.
     */
    public static void extractUnmappedReadsToFastq(File file, Strandness strandness, String from, String to) {
        System.out.println("Extracting unmappable reads...");
        SamReader reader = SamReaderFactory.makeDefault().open(file);
        SAMRecordIterator pairedIT = reader.iterator();
        boolean pairedEnd = pairedIT.next().getReadPairedFlag();
        pairedIT.close();
        ExtendedIterator<SAMRecord> it = EI.wrap(reader.iterator());
        HashMap<String, SAMRecord[]> pairedReadMap = new HashMap<>();

        try {
            BufferedWriter idWriter = new BufferedWriter(new FileWriter(file.getPath().replace(".bam", ".idMap")));
            String unmappedT2CPath = file.getPath().replace(".bam", "_unmapped_" + from + "2" + to + ".fastq");
            BufferedWriter unmappedT2CWriter = null;

            if (!pairedEnd) {
                unmappedT2CWriter = new BufferedWriter(new FileWriter(unmappedT2CPath));

                for (SAMRecord rec : it.loop()) {
                    if (!rec.getReadUnmappedFlag()) {
                        continue;
                    }

                    String seq = rec.getReadString();
                    String quality = rec.getBaseQualityString();

                    //Reverse complement the Sequence & reverse quality string if ReadNegative
                    if (rec.getReadNegativeStrandFlag()) {
                        seq = SequenceUtils.getDnaReverseComplement(seq);
                        quality = StringUtils.reverse(quality).toString();
                    }

                    String tagsToAdd = "";
                    for (SAMRecord.SAMTagAndValue s : rec.getAttributes()) {
                        tagsToAdd = tagsToAdd + "*" + s.tag + "~" + s.value + ";";
                    }

                    idWriter.append(rec.getReadName() + "_" + seq + "_" + tagsToAdd + "\n");
                    unmappedT2CWriter.append("@" + rec.getReadName() + "\n");
                    unmappedT2CWriter.append(convertNucleotides(seq, strandness, pairedEnd, from, to) + "\n");
                    unmappedT2CWriter.append("+\n");
                    unmappedT2CWriter.append(quality + "\n");
                }
            } else {
                unmappedT2CWriter = new BufferedWriter(new FileWriter(unmappedT2CPath.replace(".fastq", "_1.fastq")));
                BufferedWriter unmappedT2CWriter2 = new BufferedWriter(new FileWriter(unmappedT2CPath.replace(".fastq", "_2.fastq")));

                for (SAMRecord rec : it.loop()) {
                    //If both read & mate are mapped, skip
                    if (!rec.getReadUnmappedFlag() && !rec.getMateUnmappedFlag()) {
                        continue;
                    }

                    if (!rec.getReadUnmappedFlag() && rec.getNotPrimaryAlignmentFlag()) {
                        continue;
                    }

                    if (!pairedReadMap.containsKey(rec.getReadName())) {
                        pairedReadMap.put(rec.getReadName(), new SAMRecord[2]);
                    }

                    SAMRecord[] current = pairedReadMap.get(rec.getReadName());

                    if (rec.getFirstOfPairFlag()) {
                        current[0] = rec;
                    } else if (rec.getSecondOfPairFlag()) {
                        current[1] = rec;
                    } else {
                        throw new IllegalArgumentException("Read neither first of pair nor second of pair \n" + rec.getSAMString());
                    }

                    if (current[0] != null && current[1] != null) {
                        writePairedReadsToFastq(current[0], current[1], idWriter, unmappedT2CWriter, unmappedT2CWriter2, pairedEnd, strandness, from, to);
                        pairedReadMap.remove(rec.getReadName());
                    }
                }

                unmappedT2CWriter2.close();
            }

            unmappedT2CWriter.close();
            idWriter.close();
        } catch (IOException e) {
            System.out.println(e.getMessage());
        }
    }

    public static void writePairedReadsToFastq(SAMRecord rec1, SAMRecord rec2, BufferedWriter idWriter, BufferedWriter unmappedT2CWriter, BufferedWriter unmappedT2CWriter2, boolean pairedEnd, Strandness strandness, String from, String to) throws IOException {
        String seq1 = rec1.getReadString();
        String seq2 = rec2.getReadString();
        String quality1 = rec1.getBaseQualityString();
        String quality2 = rec2.getBaseQualityString();

        if (rec1.getReadNegativeStrandFlag()) {
            seq1 = SequenceUtils.getDnaReverseComplement(seq1);
            quality1 = StringUtils.reverse(quality1).toString();
        }
        if (rec2.getReadNegativeStrandFlag()) {
            seq2 = SequenceUtils.getDnaReverseComplement(seq2);
            quality2 = StringUtils.reverse(quality2).toString();
        }

        String tagsToAdd1 = "";
        String tagsToAdd2 = "";

        for (SAMRecord.SAMTagAndValue s : rec1.getAttributes()) {
            tagsToAdd1 = tagsToAdd1 + "*" + s.tag + "~" + s.value + ";";
        }
        for (SAMRecord.SAMTagAndValue s : rec2.getAttributes()) {
            tagsToAdd2 = tagsToAdd2 + "*" + s.tag + "~" + s.value + ";";
        }

        //idWriter.append(rec1.getReadName()+"/"+ rec1.getReadName() + "_#" + seq1 + "_#" + seq2 + "_"  + tagsToAdd1 + "\\" + rec2.getReadName() + "_#" + seq1 + "_#" + seq2 + "_" + tagsToAdd2 + "\n");
        idWriter.append(rec1.getReadName() + "_" + seq1 + "_" + tagsToAdd1 + "\\" + rec2.getReadName() + "_" + seq2 + "_" + tagsToAdd2 + "\n");

        unmappedT2CWriter.append("@" + rec1.getReadName() + "\n");
        unmappedT2CWriter.append(convertNucleotides(seq1, strandness, pairedEnd, from, to) + "\n");
        unmappedT2CWriter.append("+\n");
        unmappedT2CWriter.append(quality1 + "\n");

        unmappedT2CWriter2.append("@" + rec2.getReadName() + "\n");
        unmappedT2CWriter2.append(convertNucleotides(seq2, strandness, pairedEnd, from, to) + "\n");
        unmappedT2CWriter2.append("+\n");
        unmappedT2CWriter2.append(quality2 + "\n");
    }


    /**
     * @param sequence   Read Sequence
     * @param strandness Sequencing protocol sense or antisense
     * @param pairedEnd  Paired End sequencing
     * @return
     */
    private static String convertNucleotides(String sequence, Strandness strandness, boolean pairedEnd, String from, String to) {
        String out = sequence;

        if (strandness.equals(Strandness.Sense)) {
            out = out.replace(from, to);
        } else if (strandness.equals(Strandness.Antisense)) {
            out = out.replace(SequenceUtils.getDnaComplement(from), SequenceUtils.getDnaComplement(to));
        } else {
            throw new IllegalArgumentException("Strandness must be Sense or Antisense");
        }
        return out;
    }

}
