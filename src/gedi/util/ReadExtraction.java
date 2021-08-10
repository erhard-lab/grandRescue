package gedi.util;

import gedi.util.functions.EI;
import gedi.util.functions.ExtendedIterator;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMRecordIterator;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;

import java.io.*;
import java.nio.file.Files;
import java.nio.file.Paths;
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
    public static void extractUnmappedReadsToFastq(File file, boolean writeAll, boolean compress, ArrayList<String> tags) {
        System.out.println("Extracting unmappable reads...");
        SamReader reader = SamReaderFactory.makeDefault().open(file);
        SAMRecordIterator pairedIT = reader.iterator();
        boolean pairedEnd = pairedIT.next().getReadPairedFlag();
        pairedIT.close();
        ExtendedIterator<SAMRecord> it = EI.wrap(reader.iterator()).progress();
        HashMap<String, SAMRecord[]> pairedReadMap = new HashMap<>();

        try {
            String mappedPath = file.getPath().replace(".bam", "_mapped.fastq");
            String unmappedPath = file.getPath().replace(".bam", "_unmapped.fastq");
            String unmappedT2CPath = file.getPath().replace(".bam", "_unmapped_T2C.fastq");


            BufferedWriter mappedWriter = null;
            BufferedWriter unmappedWriter = null;
            BufferedWriter unmappedT2CWriter = null;

            if (!pairedEnd) {
                mappedWriter = new BufferedWriter(new FileWriter(mappedPath));
                unmappedWriter = new BufferedWriter(new FileWriter(unmappedPath));
                unmappedT2CWriter = new BufferedWriter(new FileWriter(unmappedT2CPath));

                for (SAMRecord rec : it.loop()) {

                    String seq = rec.getReadString();
                    String quality = rec.getBaseQualityString();

                    //Reverse complement the Sequence & reverse quality string if ReadNegative
                    if (rec.getReadNegativeStrandFlag()) {
                        seq = SequenceUtils.getDnaReverseComplement(seq);
                        quality = StringUtils.reverse(quality).toString();
                    }

                    String tagsToAdd = "_";
                    for(String s : tags){
                        Object val = rec.getAttribute(s);
                        if(val == null){
                            continue;
                        }
                        tagsToAdd = tagsToAdd + "*" + s + "~" + val + ";";
                    }

                    if (rec.getReadUnmappedFlag()) {
                        unmappedWriter.append("@" + rec.getReadName() + "\n");
                        unmappedWriter.append(seq + "\n");
                        unmappedWriter.append("+\n");
                        unmappedWriter.append(quality + "\n");
                        unmappedT2CWriter.append("@" + rec.getReadName() + "_#" + seq + tagsToAdd + "\n");
                        unmappedT2CWriter.append(seq.replace("T", "C") + "\n");
                        unmappedT2CWriter.append("+\n");
                        unmappedT2CWriter.append(quality + "\n");
                    } else {
                        mappedWriter.append("@" + rec.getReadName() + "\n");
                        mappedWriter.append(seq + "\n");
                        mappedWriter.append("+\n");
                        mappedWriter.append(quality + "\n");
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

                    String tagsToAdd1 = "_";
                    String tagsToAdd2 = "_";

                    for(String s : tags){
                        Object val = rec1.getAttribute(s);
                        if(val == null){
                            continue;
                        }
                        tagsToAdd1 = tagsToAdd1 + "*" + s + "~" + val + ";";

                        val = rec2.getAttribute(s);
                        if(val == null){
                            continue;
                        }
                        tagsToAdd2 = tagsToAdd2 + "*" + s + "~" + val + ";";
                    }
                    unmappedT2CWriter.append("@" + rec1.getReadName() + "_#" + seq1 + "_#" + seq2 + tagsToAdd1 + "\n");
                    unmappedT2CWriter.append(seq1.replace("T", "C") + "\n");
                    unmappedT2CWriter.append("+\n");
                    unmappedT2CWriter.append(quality1 + "\n");
                    unmappedT2CWriter2.append("@" + rec2.getReadName() + "_#" + seq1 + "_#" + seq2 + tagsToAdd2 + "\n");
                    unmappedT2CWriter2.append(seq2.replace("T", "C") + "\n");
                    unmappedT2CWriter2.append("+\n");
                    unmappedT2CWriter2.append(quality2 + "\n");
                }

                unmappedT2CWriter2.close();
            }

            if (!pairedEnd) {
                mappedWriter.close();
                unmappedWriter.close();
            }
            unmappedT2CWriter.close();

            if (!writeAll) {
                Files.delete(Paths.get(mappedPath));
                Files.delete(Paths.get(unmappedPath));
            } else if (compress) {
                gzip(mappedPath);
                gzip(unmappedPath);
            }
        } catch (IOException e) {
            System.out.println(e.getMessage());
        }
    }

    public static void extractUnmappedReadsToFastq(String file, boolean writeAll, boolean compress, ArrayList<String> tags) {
        extractUnmappedReadsToFastq(new File(file), writeAll, compress, tags);
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
