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
    public static void extractUnmappedReadsToFastq(File file, ArrayList<String> tags, Strandness strandness) {
        System.out.println("Extracting unmappable reads...");
        SamReader reader = SamReaderFactory.makeDefault().open(file);
        SAMRecordIterator pairedIT = reader.iterator();
        boolean pairedEnd = pairedIT.next().getReadPairedFlag();
        pairedIT.close();
        ExtendedIterator<SAMRecord> it = EI.wrap(reader.iterator());
        HashMap<String, SAMRecord[]> pairedReadMap = new HashMap<>();

        try {
            BufferedWriter idWriter = new BufferedWriter(new FileWriter(file.getPath().replace(".bam", ".idMap")));
            String unmappedT2CPath = file.getPath().replace(".bam", "_unmapped_T2C.fastq");
            BufferedWriter unmappedT2CWriter = null;

            if (!pairedEnd) {
                unmappedT2CWriter = new BufferedWriter(new FileWriter(unmappedT2CPath));

                for (SAMRecord rec : it.loop()) {

                    String seq = rec.getReadString();
                    String quality = rec.getBaseQualityString();

                    //Reverse complement the Sequence & reverse quality string if ReadNegative
                    if (rec.getReadNegativeStrandFlag()) {
                        seq = SequenceUtils.getDnaReverseComplement(seq);
                        quality = StringUtils.reverse(quality).toString();
                    }

                    String tagsToAdd = "";
                    for (String s : tags) {
                        Object val = rec.getAttribute(s);
                        if (val == null) {
                            continue;
                        }
                        tagsToAdd = tagsToAdd + "*" + s + "~" + val + ";";
                    }

                    if (rec.getReadUnmappedFlag()) {
                        idWriter.append(rec.getReadName() + "/" + rec.getReadName() + "_#" + seq + "_" + tagsToAdd + "\n");
                        unmappedT2CWriter.append("@" + rec.getReadName() + "\n");
                        unmappedT2CWriter.append(convertNucleotides(seq, strandness, pairedEnd) + "\n");
                        unmappedT2CWriter.append("+\n");
                        unmappedT2CWriter.append(quality + "\n");
                    }
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
                }

                for (Map.Entry<String, SAMRecord[]> entry : pairedReadMap.entrySet()) {
                    SAMRecord[] current = entry.getValue();
                    SAMRecord rec1 = current[0];
                    SAMRecord rec2 = current[1];
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

                    for (String s : tags) {
                        Object val = rec1.getAttribute(s);
                        if (val == null) {
                            continue;
                        }
                        tagsToAdd1 = tagsToAdd1 + "*" + s + "~" + val + ";";

                        val = rec2.getAttribute(s);
                        if (val == null) {
                            continue;
                        }
                        tagsToAdd2 = tagsToAdd2 + "*" + s + "~" + val + ";";
                    }

                    idWriter.append(rec1.getReadName() + "/" + rec1.getReadName() + "_#" + seq1 + "_#" + seq2 + "_" + tagsToAdd1 + "\\" + rec2.getReadName() + "_#" + seq1 + "_#" + seq2 + "_" + tagsToAdd2 + "\n");
                    //unmappedT2CWriter.append("@" + rec1.getReadName() + "_#" + seq1 + "_#" + seq2 + tagsToAdd1 + "\n");
                    unmappedT2CWriter.append("@" + rec1.getReadName() + "\n");
                    unmappedT2CWriter.append(convertNucleotides(seq1, strandness, pairedEnd) + "\n");
                    unmappedT2CWriter.append("+\n");
                    unmappedT2CWriter.append(quality1 + "\n");
                    //unmappedT2CWriter2.append("@" + rec2.getReadName() + "_#" + seq1 + "_#" + seq2 + tagsToAdd2 + "\n");
                    unmappedT2CWriter2.append("@" + rec2.getReadName() + "\n");
                    unmappedT2CWriter2.append(convertNucleotides(seq2, strandness, pairedEnd) + "\n");
                    unmappedT2CWriter2.append("+\n");
                    unmappedT2CWriter2.append(quality2 + "\n");
                }

                unmappedT2CWriter2.close();
            }

            unmappedT2CWriter.close();
            idWriter.close();
        } catch (IOException e) {
            System.out.println(e.getMessage());
        }
    }

    public static void extractUnmappedReadsToFastq(String file, ArrayList<String> tags, Strandness strandness) {
        extractUnmappedReadsToFastq(new File(file), tags, strandness);
    }

    /**
     * @param file .sam-File of all reads from 4sU-labelled experiment
     *             Writes 3 bam-files: All mapped, all unmapped reads and unmapped reads where all T's are converted to C's.
     *             Saves the original read sequence (with T->C MM's from 4su labeling) in ReadID, then changes all T's to C's in Read Sequence.
     *             These reads can then be mapped to the pseudo-genome where all T's to C's are converted.
     */
    public static void extractUnmappedReadsToSAM(File file) {
        SamReader reader = SamReaderFactory.makeDefault().open(file);
        ExtendedIterator<SAMRecord> it = EI.wrap(reader.iterator()).progress();

        try {
            BufferedWriter mappedWriter = new BufferedWriter(new FileWriter(file.getPath().replace(".bam", "_mapped.sam")));
            BufferedWriter unmappedWriter = new BufferedWriter(new FileWriter(file.getPath().replace(".bam", "_unmapped.sam")));
            BufferedWriter unmappedT2CWriter = new BufferedWriter(new FileWriter(file.getPath().replace(".bam", "_unmapped_T2C.sam")));

            mappedWriter.append(reader.getFileHeader().getTextHeader());
            unmappedWriter.append(reader.getFileHeader().getTextHeader());
            unmappedT2CWriter.append(reader.getFileHeader().getTextHeader());

            for (SAMRecord rec : it.loop()) {
                if (rec.getReadUnmappedFlag()) {
                    unmappedWriter.append(rec.getSAMString());

                    rec.setReadName(rec.getReadName() + "_#" + rec.getReadString());
                    rec.setReadString(rec.getReadString().replace("T", "C"));
                    unmappedT2CWriter.append(rec.getSAMString());
                } else {
                    mappedWriter.append(rec.getSAMString());
                }
            }
            mappedWriter.close();
            unmappedWriter.close();


        } catch (IOException e) {
            System.out.println(e.getMessage());
        }
    }

    public static void extractUnmappedReadsToSAM(String file) {
        extractUnmappedReadsToSAM(new File(file));
    }

    public static void T2CFastq(String fastq) {
        try {
            BufferedWriter writer = new BufferedWriter(new FileWriter(fastq.replace(".fastq", "_T2C.fastq")));
            FunctorUtils.FixedBlockIterator<String, ArrayList<String>> it = EI.lines(fastq).block(4);
            for (ArrayList<String> s : it.loop()) {
                writer.append(s.get(0) + "\n");
                writer.append(s.get(1).replace("T", "C") + "\n");
                writer.append(s.get(2) + "\n");
                writer.append(s.get(3) + "\n");
            }
            writer.close();
        } catch (IOException e) {
            e.printStackTrace();
        }
    }

    /**
     * @param sequence   Read Sequence
     * @param strandness Sequencing protocol sense or antisense
     * @param pairedEnd  Paired End sequencing
     * @return
     */
    private static String convertNucleotides(String sequence, Strandness strandness, boolean pairedEnd) {
        String out = sequence;

        if (!pairedEnd) {
            if (strandness.equals(Strandness.Sense)) {
                out = out.replace("T", "C");
            } else if (strandness.equals(Strandness.Antisense)) {
                out = out.replace("A", "G");
            } else {
                throw new IllegalArgumentException("Strandness must be Sense or Antisense");
            }
        } else {
            if (strandness.equals(Strandness.Sense)) {
                out = out.replace("T", "C");
            } else if (strandness.equals(Strandness.Antisense)) {
                out = out.replace("A", "G");
            } else {
                throw new IllegalArgumentException("Strandness must be Sense or Antisense");
            }
        }

        return out;
    }

    private static void gzip(String file) {
        System.out.println("Gzipping " + file);
        try {
            FileInputStream in = new FileInputStream(file);
            FileOutputStream out = new FileOutputStream(file + ".gz");
            GZIPOutputStream gzOut = new GZIPOutputStream(out);
            byte[] buffer = new byte[1024];
            int len;
            while ((len = in.read(buffer)) != -1) {
                gzOut.write(buffer, 0, len);
            }

            gzOut.close();
            out.close();
            in.close();

            new File(file).delete();
        } catch (IOException e) {
            e.printStackTrace();
        }
    }


}
