package gedi.util;

import gedi.bam.tools.BamUtils;
import gedi.centeredDiskIntervalTree.CenteredDiskIntervalTreeStorage;
import gedi.core.data.annotation.Transcript;
import gedi.core.data.reads.DefaultAlignedReadsData;
import gedi.core.genomic.Genomic;
import gedi.core.reference.Chromosome;
import gedi.core.region.ArrayGenomicRegion;
import gedi.core.region.GenomicRegion;
import gedi.core.region.ImmutableReferenceGenomicRegion;
import gedi.core.region.ReferenceGenomicRegion;
import gedi.region.bam.FactoryGenomicRegion;
import gedi.util.dynamic.DynamicObject;
import gedi.util.functions.EI;
import gedi.util.functions.ExtendedIterator;
import gedi.util.functions.IterateIntoSink;
import htsjdk.samtools.*;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.nio.charset.StandardCharsets;
import java.util.*;
import java.util.function.Function;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

import static gedi.util.SequenceUtil.calculateMdAndNmTags;

public class BAMUtils {

    static int plusCount = 0;
    static int minusCount = 0;
    static int noInduced = 0;


    public static void checkReadMapping(String bam, String fastq, String origGenome, String mappedToGenome, boolean positionEquality) {
        try {
            //Two Mappings are equal either if mapped to same Position (+/- X bases) or to same gene
            boolean geneNameEquality = !positionEquality;

            SamReader reader = SamReaderFactory.makeDefault().open(new File(bam));
            ExtendedIterator<SAMRecord> it = EI.wrap(reader.iterator());
            BufferedWriter writer = new BufferedWriter(new FileWriter(bam.replace(".bam", ".mapStat")));
            Genomic originalGenome = Genomic.get(origGenome);
            Genomic mappedGenome = Genomic.get(mappedToGenome);

            //Mapping statistics
            int totalReads = 0;
            int unmappedReads = 0;
            int mappedReads = 0;    // = totalReads-unmappedReads
            int uniqueReads = 0;
            int notPrimaryReads = 0;
            int correctMap = 0;
            int correctUnique = 0;
            int correctPrimary = 0;
            int falseMap = 0;
            int totalTCMM = countGlobalTCMM(fastq);
            int retainedTCMM = 0;

            for (SAMRecord rec : it.loop()) {
                totalReads++;

                if (rec.getReadUnmappedFlag()) {
                    unmappedReads++;
                    continue;
                }

                if (rec.getIntegerAttribute("NH") == 1) {
                    uniqueReads++;
                }

                if (rec.getNotPrimaryAlignmentFlag()) {
                    notPrimaryReads++;
                    continue;
                }


                Pattern p = Pattern.compile("([\\w&&[^_]]+)_([\\w&&[^_]]+)_(\\d+)_(\\d+)(_\\d+)*");
                Matcher m = p.matcher(rec.getReadName());
                m.find();

                if (m.group(5) != null && rec.getIntegerAttribute("NH") == 1) {
                    retainedTCMM = retainedTCMM + Integer.valueOf(m.group(5).substring(1));
                }

                //Check for same chromosome
                if (!rec.getReferenceName().equals(m.group(1))) {
                    falseMap++;
                    continue;
                }

                //Check for same strand
                String strand = "+";
                if (rec.getReadNegativeStrandFlag()) {
                    strand = "-";
                }
                if (!rec.getReadNegativeStrandFlag() && !m.group(2).equals("plus")) {
                    falseMap++;
                    continue;
                } else if (rec.getReadNegativeStrandFlag() && !m.group(2).equals("minus")) {
                    falseMap++;
                    continue;
                }

                if (positionEquality) {
                    if (rec.getStart() != Integer.valueOf(m.group(3)) + 1) {
                        falseMap++;
                        continue;
                    }
                    if (rec.getEnd() != Integer.valueOf(m.group(4))) {
                        falseMap++;
                        continue;
                    }

                } else {

                    Chromosome chromosome = Chromosome.obtain("" + rec.getReferenceName() + strand);
                    GenomicRegion mappedRegion = new ArrayGenomicRegion(rec.getAlignmentStart() - 1, rec.getAlignmentEnd());
                    ExtendedIterator<ImmutableReferenceGenomicRegion<Transcript>> mappedEI = mappedGenome.getTranscripts().ei(chromosome, mappedRegion);
                    GenomicRegion origRegion = new ArrayGenomicRegion(Integer.valueOf(m.group(3)), Integer.valueOf(m.group(4)));
                    ExtendedIterator<ImmutableReferenceGenomicRegion<Transcript>> origEI = originalGenome.getTranscripts().ei(chromosome, origRegion);

                    if (mappedEI.hasNext() && origEI.hasNext()) {
                        String mappedGene = mappedEI.first().getData().getGeneId();
                        String origGene = origEI.first().getData().getGeneId();
                        if (!mappedGene.equals(origGene)) {
                            falseMap++;
                            continue;
                        }
                    } else if (origEI.hasNext() && !mappedEI.hasNext()) {
                        falseMap++;
                        continue;

                    } else if (!origEI.hasNext() && mappedEI.hasNext()) {
                        falseMap++;
                        continue;
                    }
                }

                if (!rec.getNotPrimaryAlignmentFlag()) {
                    correctPrimary++;
                }
                if (rec.getIntegerAttribute("NH") == 1) {
                    correctUnique++;
                }
                correctMap++;
            }

            int originalT2CMMs = countGlobalTCMM(fastq);
            mappedReads = totalReads - unmappedReads;

            writer.append("Read mapping statistics for " + bam + "\n\n");
            writer.append("Parameters:\n");
            writer.append("Original Genome: " + origGenome + " - mapped Genome: " + mappedToGenome + "\n");
            if (positionEquality) {
                writer.append("Correct mapping defined by read position");
            } else {
                writer.append("Correct mapping defined by mapping to same gene");
            }
            writer.append("\n\n");


            writer.append("Total reads: " + totalReads + "\n\n");

            writer.append("Mapped reads: " + mappedReads + "\n");
            writer.append("Unmapped reads: " + unmappedReads + "\n");
            writer.append("% unmapped reads: " + ((double) unmappedReads / totalReads) + "\n\n");

            writer.append("Unique reads: " + uniqueReads + "\n");
            writer.append("% Unique reads: " + ((double) uniqueReads / mappedReads) + "\n");
            writer.append("Not primary reads: " + notPrimaryReads + "\n");
            writer.append("Average # of mappings per non-unique read: " + ((double) (mappedReads - uniqueReads) / (notPrimaryReads)) + "\n\n");

            writer.append("% correctly mapped reads: " + ((double) correctMap / mappedReads) + "\n");
            writer.append("% correctly mapped primary reads: " + ((double) correctPrimary / (mappedReads - notPrimaryReads)) + "\n");
            writer.append("% correctly mapped unique reads: " + ((double) correctUnique / uniqueReads) + "\n\n");

            writer.append("Original T2C-MM count: " + originalT2CMMs + "\n");
            writer.append("still existing T->C MM's in unique reads: " + retainedTCMM + "\n");
            writer.append("% T->C MM's lost: " + (1.0 - (double) retainedTCMM / originalT2CMMs) + "\n");
            writer.append("TCMM's lost per unmapped read: " + ((double) (originalT2CMMs - retainedTCMM) / unmappedReads) + "\n");
            writer.append("TCMM's lost per nonprimary read: " + ((double) (originalT2CMMs - retainedTCMM) / notPrimaryReads));


            writer.close();

        } catch (IOException e) {
            System.out.println(e.getMessage());
        }
    }

    public static void pseudoMapping(String bamFile, String originalGenome, String mappedGenome, boolean mapToPlusStrand) {
        SamReader reader = SamReaderFactory.makeDefault().open(new File(bamFile));
        ExtendedIterator<SAMRecord> it = EI.wrap(reader.iterator());
        Map<String, Set<SAMRecord>> nonUniqueMap = new HashMap<>();
        List<SAMRecord> correctMaps = new ArrayList<>();
        Pattern p = Pattern.compile("([\\w&&[^_]]+)_([\\w&&[^_]]+)_(\\d+)_(\\d+)(_\\d+)*");


        Genomic origGenome = Genomic.get(originalGenome);
        Genomic mapGenome = Genomic.get(mappedGenome);

        int primaryReads = 0;
        int correctMap = 0;
        int notSameGene = 0;
        int notSamePos = 0;


        for (SAMRecord rec : it.loop()) {

            if (rec.getReadNegativeStrandFlag()) {
                continue;
            }
            if (rec.getStringAttribute("MD").contains("N")) {
                continue;
            }
            if (rec.getCigarString().contains("N")) {
                continue;
            }


            Matcher m = p.matcher(rec.getReadName());
            m.find();

            if (Integer.valueOf(m.group(4)) - Integer.valueOf(m.group(3)) > 70) {
                continue;
            }
            //Put non-unique reads in map to resolve later
            if (rec.getIntegerAttribute("NH") != 1) {
                if (nonUniqueMap.get(rec.getReadName()) == null) {
                    nonUniqueMap.put(rec.getReadName(), new HashSet<>());
                }

                nonUniqueMap.get(rec.getReadName()).add(rec);
                continue;
            }


            primaryReads++;


            Boolean geneMap = sameGeneMap(rec, origGenome, mapGenome, mapToPlusStrand);
            if (geneMap == null) {
                //orig mapping und new mapping both didn't map to gene -> ignore
                primaryReads--;
                continue;
            } else if (geneMap == false) {
                /*notSameGene++;
                continue;*/
            } else if (geneMap == true) {
                //just continue with position mapping
                ;
            }

            Boolean posMap = samePosMap(rec, origGenome, mapGenome, p);
            if (posMap) {
                correctMaps.add(rec);
                correctMap++;
                continue;
            } else {
                notSamePos++;
                continue;
            }


        }

        //Deal with non-unique Reads
        for (String s : nonUniqueMap.keySet()) {
            Set<SAMRecord> set = nonUniqueMap.get(s);
            primaryReads++;

            Set<SAMRecord> correctNonPrimary = new HashSet<>();

            SAMRecord bestFit = null;
            int shortestMD = 0;

            for (SAMRecord rec : set) {
                if (shortestMD == 0) {
                    bestFit = rec;
                    shortestMD = rec.getStringAttribute("MD").length();
                    continue;
                }
                if (rec.getStringAttribute("MD").length() < shortestMD) {
                    bestFit = rec;
                    shortestMD = rec.getStringAttribute("MD").length();
                }
            }
            Boolean mapPos = samePosMap(bestFit, origGenome, mapGenome, p);
            if (mapPos) {
                correctNonPrimary.add(bestFit);
                correctMap++;
            } else {
                continue;
            }
            correctMaps.add(bestFit);
        }

        System.out.println();
        System.out.println(primaryReads);
        System.out.println(correctMap);
        System.out.println((double) correctMap / primaryReads);
        System.out.println(notSamePos);
        System.out.println(plusCount);
        System.out.println(minusCount);
    }

    public static void checkPseudoMapping(String bam, String fastq, String indexFile, String origGenome, String mappedToGenome, boolean mapToPlusStrand) {
        try {
            SamReader reader = SamReaderFactory.makeDefault().open(new File(bam));
            ExtendedIterator<SAMRecord> it = EI.wrap(reader.iterator());
            BufferedWriter writer = new BufferedWriter(new FileWriter(new File(bam.replace(".bam", ".mapStat"))));
            Genomic originalGenome = Genomic.get(origGenome);
            Genomic mappedGenome = Genomic.get(mappedToGenome);
            Map<String, IndexEntry> map = getIndexMap(indexFile);


            //Mapping statistics
            int totalReads = 0;
            int unmappedReads = 0;
            int mappedReads = 0;    // = totalReads-unmappedReads
            int uniqueReads = 0;
            int notPrimaryReads = 0;
            int correctMap = 0;
            int correctUnique = 0;
            int correctPrimary = 0;
            int falseMap = 0;
            int falsePos = 0;
            int falseGene = 0;
            int mapToNeg = 0;
            int totalTCMM = countGlobalTCMM(fastq);
            int retainedTCMM = 0;

            for (SAMRecord rec : it.loop()) {
                totalReads++;

                if (rec.getReadUnmappedFlag()) {
                    unmappedReads++;
                    continue;
                }

                //Some reads still mapped to - Strand even though genome only contains +??
                if (mapToPlusStrand && rec.getReadNegativeStrandFlag()) {
                    mapToNeg++;
                    continue;
                }

                if (rec.getIntegerAttribute("NH") == 1) {
                    uniqueReads++;
                }

                if (rec.getNotPrimaryAlignmentFlag()) {
                    notPrimaryReads++;
                    continue;
                }


                Pattern p = Pattern.compile("([\\w&&[^_]]+)_([\\w&&[^_]]+)_(\\d+)_(\\d+)(_\\d+)*");
                Matcher m = p.matcher(rec.getReadName());
                m.find();

                if (m.group(5) != null && rec.getIntegerAttribute("NH") == 1) {
                    retainedTCMM = retainedTCMM + Integer.valueOf(m.group(5).substring(1));
                }

                //Check for same chromosome
                if (!rec.getReferenceName().equals(m.group(1))) {
                    falseMap++;
                    continue;
                }

                //Check for same strand
                String strand = "+";

                if (rec.getReadNegativeStrandFlag()) {
                    strand = "-";
                }
                if (!rec.getReadNegativeStrandFlag() && !m.group(2).equals("plus")) {
                    falseMap++;
                    continue;
                } else if (rec.getReadNegativeStrandFlag() && !m.group(2).equals("minus")) {
                    falseMap++;
                    continue;
                }

                //Test for same Gene Name
                Chromosome chromosome = Chromosome.obtain("" + rec.getReferenceName() + strand);
                GenomicRegion mappedRegion = new ArrayGenomicRegion(rec.getAlignmentStart() - 1, rec.getAlignmentEnd());
                ExtendedIterator<ImmutableReferenceGenomicRegion<String>> mappedEI = mappedGenome.getGenes().ei(chromosome, mappedRegion);
                GenomicRegion origRegion = new ArrayGenomicRegion(Integer.valueOf(m.group(3)), Integer.valueOf(m.group(4)));
                ExtendedIterator<ImmutableReferenceGenomicRegion<String>> origEI = originalGenome.getGenes().ei(chromosome, origRegion);

                String mappedGene = "";
                if (mappedEI.hasNext() && origEI.hasNext()) {
                    mappedGene = mappedEI.first().getData();
                    String origGene = origEI.first().getData();
                    if (!mappedGene.equals(origGene)) {
                        falseGene++;
                        continue;
                    }
                } else if (origEI.hasNext() && !mappedEI.hasNext()) {
                    falseGene++;
                    continue;

                } else if (!origEI.hasNext() && mappedEI.hasNext()) {
                    falseGene++;
                    continue;
                }


                //Test for same Position
                IndexEntry indexEntry = map.get(mappedGene);
                int inducedStart = indexEntry.origPosStart + (rec.getAlignmentStart() - 1 - indexEntry.pseudoPosStart);
                //int inducedEnd = inducedStart + (rec.getAlignmentEnd() - rec.getAlignmentStart() + 1);
                if (rec.getReadName().contains("minus")) {
                    inducedStart = indexEntry.origPosStart + (indexEntry.pseudoPosEnd - 1 - rec.getAlignmentEnd());
                    //inducedEnd = inducedStart + (rec.getAlignmentEnd() - rec.getAlignmentStart() + 1);
                }


                if (inducedStart != Integer.valueOf(m.group(3))) {
                    falsePos++;
                    continue;
                }


                if (!rec.getNotPrimaryAlignmentFlag()) {
                    correctPrimary++;
                }
                if (rec.getIntegerAttribute("NH") == 1) {
                    correctUnique++;
                }
                correctMap++;
            }

            int originalT2CMMs = countGlobalTCMM(fastq);
            mappedReads = totalReads - unmappedReads;

            writer.append("Read mapping statistics for " + bam + "\n\n");
            writer.append("Parameters:\n");
            writer.append("Original Genome: " + origGenome + " - mapped Genome: " + mappedToGenome + "\n");
            writer.append("Mapping to pseudogenome based on gene name & read position");

            writer.append("\n\n");


            writer.append("Total reads: " + totalReads + "\n\n");

            writer.append("Mapped reads: " + mappedReads + "\n");
            writer.append("Unmapped reads: " + unmappedReads + "\n");
            writer.append("% unmapped reads: " + ((double) unmappedReads / totalReads) + "\n\n");

            writer.append("Unique reads: " + uniqueReads + "\n");
            writer.append("% Unique reads: " + ((double) uniqueReads / mappedReads) + "\n");
            writer.append("Not primary reads: " + notPrimaryReads + "\n");
            writer.append("Average # of mappings per non-unique read: " + ((double) (mappedReads - uniqueReads) / (notPrimaryReads)) + "\n\n");

            writer.append("Total falsely mapped reads: " + (falseMap + falseGene + falsePos));
            writer.append("% Falsely mapped reads: " + ((double) (falseMap + falseGene + falsePos) / mappedReads));
            writer.append("Mapped to wrong gene: " + falseGene);
            writer.append("Mapped to right gene, wrong Pos: " + falsePos);

            writer.append("% correctly mapped reads: " + ((double) correctMap / mappedReads) + "\n");
            writer.append("% correctly mapped primary reads: " + ((double) correctPrimary / (mappedReads - notPrimaryReads)) + "\n");
            writer.append("% correctly mapped unique reads: " + ((double) correctUnique / uniqueReads) + "\n\n");

            writer.append("Original T2C-MM count: " + originalT2CMMs + "\n");
            writer.append("still existing T->C MM's in unique reads: " + retainedTCMM + "\n");
            writer.append("% T->C MM's lost: " + (1.0 - (double) retainedTCMM / originalT2CMMs) + "\n");
            writer.append("TCMM's lost per unmapped read: " + ((double) (originalT2CMMs - retainedTCMM) / unmappedReads) + "\n");
            writer.append("TCMM's lost per nonprimary read: " + ((double) (originalT2CMMs - retainedTCMM) / notPrimaryReads));


            writer.close();

        } catch (IOException e) {
            System.out.println(e.getMessage());
        }
    }

    public static Map<String, IndexEntry> getIndexMap(String file) {
        Map<String, IndexEntry> map = new HashMap<>();

        try {
            Scanner sc = new Scanner(new File(file));

            sc.nextLine();
            while (sc.hasNextLine()) {
                String line = sc.nextLine();
                String[] lineArr = line.split("\t");

                Pattern p = Pattern.compile("(\\d+)/(\\d+)");
                Matcher m1 = p.matcher(lineArr[1]);
                m1.find();
                Matcher m2 = p.matcher(lineArr[2]);
                m2.find();
                map.put(lineArr[0], new IndexEntry(Integer.valueOf(m1.group(1)), Integer.valueOf(m1.group(2)), Integer.valueOf(m2.group(1)), Integer.valueOf(m2.group(2))));
            }
        } catch (IOException e) {
            System.out.println(e.getMessage());
        }


        return map;
    }

    public static int countGlobalTCMM(String fastQFile) {
        int MMCount = 0;

        try {
            Scanner scanner = new Scanner(new File(fastQFile));

            while (scanner.hasNext()) {
                String nextLine = scanner.nextLine();
                if (!nextLine.startsWith("@")) {
                    continue;
                }

                Pattern p = Pattern.compile("([\\w&&[^_]]+)_([\\w&&[^_]]+)_(\\d+)_(\\d+)(_\\d+)*");
                Matcher m = p.matcher(nextLine);
                m.find();

                if (m.group(5) != null) {
                    MMCount = MMCount + Integer.valueOf(m.group(5).substring(1));
                }
            }
            scanner.close();
        } catch (IOException e) {
            System.out.println(e.getMessage());
        }
        return MMCount;


    }

    public static int countTCMMinSeq(List<SAMRecord> records, Genomic origGenome, Genomic mapGenome) {
        int overallT = 0;
        int TCMM = 0;

        for (SAMRecord rec : records) {

            String mappedStrand = "+";
            if (rec.getReadNegativeStrandFlag()) {
                mappedStrand = "-";
            }
            Chromosome chrom = Chromosome.obtain(rec.getReferenceName() + mappedStrand);

            //Saved the original ReadSequence before T2C-Conversion in ReadName as "@ATGC"
            CharSequence mappedSequence = rec.getReadName().substring(1);
            String[] inducedPositions = inducedPositions(rec, origGenome, mapGenome, true);
            CharSequence origSequence = origGenome.getSequence(chrom, new ArrayGenomicRegion(Integer.valueOf(inducedPositions[0]), Integer.valueOf(inducedPositions[1])));

            if (mappedSequence.length() != origSequence.length()) {
                throw new IllegalArgumentException("Mapped & orig sequence don't have same length");
            }

            for (int i = 0; i < origSequence.length(); i++) {
                if (origSequence.charAt(i) == 'T') {
                    overallT++;
                    if (mappedSequence.charAt(i) == 'C') {
                        TCMM++;
                    }
                }
            }

        }

        System.out.println((double) TCMM / overallT);

        return TCMM;
    }

    public static Boolean sameGeneMap(SAMRecord rec, Genomic origGenome, Genomic mapGenome, boolean mapToPlusStrand) {
        Pattern p = Pattern.compile("([\\w&&[^_]]+)_([\\w&&[^_]]+)_(\\d+)_(\\d+)(_\\d+)*");
        Matcher m = p.matcher(rec.getReadName());
        m.find();
        if (!m.group(1).equals(rec.getReferenceName())) {
            return false;
        }
        if (rec.getReadNegativeStrandFlag() && mapToPlusStrand) {
            return false;
        }

        String origStrand = "";
        if (m.group(2).equals("plus")) {
            origStrand = "+";
        } else if (m.group(2).equals("minus")) {
            origStrand = "-";
        } else {
            throw new IllegalArgumentException("Strand in Read Name is neither plus nor minus");
        }

        String mapStrand = "+";
        if (!mapToPlusStrand && rec.getReadNegativeStrandFlag()) {
            mapStrand = "-";
        }

        if (!mapToPlusStrand && !origStrand.equals(mapStrand)) {
            return false;
        }

        Chromosome origChromosome = Chromosome.obtain("" + rec.getReferenceName() + origStrand);
        Chromosome mapChromosome = Chromosome.obtain("" + rec.getReferenceName() + mapStrand);
        GenomicRegion mappedRegion = new ArrayGenomicRegion(rec.getAlignmentStart() - 1, rec.getAlignmentEnd());
        ExtendedIterator<ImmutableReferenceGenomicRegion<String>> mappedEI = mapGenome.getGenes().ei(mapChromosome, mappedRegion);
        GenomicRegion origRegion = new ArrayGenomicRegion(Integer.valueOf(m.group(3)), Integer.valueOf(m.group(4)));
        ExtendedIterator<ImmutableReferenceGenomicRegion<String>> origEI = origGenome.getGenes().ei(origChromosome, origRegion);

        String mapGene = "";
        String origGene = "";
        ImmutableReferenceGenomicRegion<String> mapData = null;
        ImmutableReferenceGenomicRegion<String> origData = null;
        try {
            mapData = mappedEI.first();
            origData = origEI.first();

            mapGene = mapData.getData();
            origGene = origData.getData();


        } catch (NullPointerException e) {
            if (origData == null && mapData == null) {
                return null;
            } else {
                return false;
            }
        }


        if (!mapGene.equals(origGene)) {
/*            System.out.println(mapGene);
            System.out.println(origGene);
            System.out.println(rec.getSAMString());*/
            return false;
        } else {
            return true;
        }
    }

    public static Boolean samePosMap(SAMRecord rec, Genomic origGenome, Genomic mapGenome, Pattern p) {
        Matcher m = p.matcher(rec.getReadName());
        m.find();

        String mapStrand = "+";
        if (rec.getReadNegativeStrandFlag()) {
            mapStrand = "-";
        }


        Chromosome mapChromosome = Chromosome.obtain("" + rec.getReferenceName() + mapStrand);
        GenomicRegion mappedRegion = new ArrayGenomicRegion(rec.getAlignmentStart() - 1, rec.getAlignmentEnd());
        ExtendedIterator<ImmutableReferenceGenomicRegion<String>> mappedEI = mapGenome.getGenes().ei(mapChromosome, mappedRegion);

        String mapGene = "";
        ImmutableReferenceGenomicRegion<String> mapData = null;
        mapData = mappedEI.first();
        mapGene = mapData.getData();

        Function<String, ReferenceGenomicRegion<String>> origGeneMapFunc = origGenome.getGeneMapping();
        Function<String, ReferenceGenomicRegion<String>> mapGeneMapFunc = mapGenome.getGeneMapping();

        ReferenceGenomicRegion origRGR = origGeneMapFunc.apply(mapGene);
        ReferenceGenomicRegion mapRGR = mapGeneMapFunc.apply(mapGene);


        int inducedStart = origRGR.map(mapRGR.induce(rec.getAlignmentStart() - 1));

        if (origRGR.toString().contains("-:")) {
            inducedStart = origRGR.map(mapRGR.induce(rec.getAlignmentEnd() - 1));
        }

        if (inducedStart == Integer.valueOf(m.group(3))) {
            return true;
        } else {
            if (rec.getReadName().contains("plus")) {
                plusCount++;
            } else if (rec.getReadName().contains("minus")) {
                minusCount++;
            }

            return false;
        }
    }

    public static String[] inducedPositions(SAMRecord rec, Genomic origGenome, Genomic mapGenome, boolean toPlusStrand) {
        String[] inducedPositions = new String[3];
        String mapStrand = "+";
        if (!toPlusStrand && rec.getReadNegativeStrandFlag()) {
            mapStrand = "-";
        }

        Chromosome mapChromosome = Chromosome.obtain("" + rec.getReferenceName() + mapStrand);
        GenomicRegion mappedRegion = new ArrayGenomicRegion(rec.getAlignmentStart() - 1, rec.getAlignmentEnd());
        ExtendedIterator<ImmutableReferenceGenomicRegion<String>> mappedEI = mapGenome.getGenes().ei(mapChromosome, mappedRegion);

        String mapGene = "";
        ImmutableReferenceGenomicRegion<String> mapData = null;
        mapData = mappedEI.first();
        mapGene = mapData.getData();

        Function<String, ReferenceGenomicRegion<String>> origGeneMapFunc = origGenome.getGeneMapping();
        Function<String, ReferenceGenomicRegion<String>> mapGeneMapFunc = mapGenome.getGeneMapping();

        ReferenceGenomicRegion origRGR = origGeneMapFunc.apply(mapGene);
        ReferenceGenomicRegion mapRGR = mapGeneMapFunc.apply(mapGene);

        int inducedStart = origRGR.map(mapRGR.induce(rec.getAlignmentStart() - 1));
        if (origRGR.toString().contains("-:")) {
            inducedStart = origRGR.map(mapRGR.induce(rec.getAlignmentEnd() - 1));
        }

        inducedPositions[0] = origRGR.getReference().getName();

        inducedPositions[1] = "+";
        if (origRGR.getReference().isMinus()) {
            inducedPositions[1] = "-";
        }

        inducedPositions[2] = "" + inducedStart;
        return inducedPositions;
    }

    public static ReferenceGenomicRegion inducedRGR(SAMRecord rec, Genomic origGenome, Genomic mapGenome, boolean toPlusStrand) {
        String mapStrand = "+";
        if (!toPlusStrand && rec.getReadNegativeStrandFlag()) {
            mapStrand = "-";
        }

        Chromosome mapChromosome = Chromosome.obtain("" + rec.getReferenceName() + mapStrand);
        GenomicRegion mappedRegion = new ArrayGenomicRegion(rec.getAlignmentStart() - 1, rec.getAlignmentEnd());
        ExtendedIterator<ImmutableReferenceGenomicRegion<String>> mappedEI = mapGenome.getGenes().ei(mapChromosome, mappedRegion);

        String mapGene = "";
        ImmutableReferenceGenomicRegion<String> mapData = null;
        mapData = mappedEI.first();
        mapGene = mapData.getData();

        Function<String, ReferenceGenomicRegion<String>> origGeneMapFunc = origGenome.getGeneMapping();
        Function<String, ReferenceGenomicRegion<String>> mapGeneMapFunc = mapGenome.getGeneMapping();

        ReferenceGenomicRegion origRGR = origGeneMapFunc.apply(mapGene);
        ReferenceGenomicRegion mapRGR = mapGeneMapFunc.apply(mapGene);


        int inducedStart = origRGR.map(mapRGR.induce(rec.getAlignmentStart() - 1));

        if (origRGR.toString().contains("-:")) {
            inducedStart = origRGR.map(mapRGR.induce(rec.getAlignmentEnd() - 1));
        }

        int inducedEnd = inducedStart + BamUtils.getArrayGenomicRegion(rec).length();


        return origRGR;
    }

    public static void samOutputFromPseudoMapping(String bamFile, String referenceBam, Genomic origGenome, Genomic mapGenome, boolean toPlusStrand) {
        SamReader reader = SamReaderFactory.makeDefault().open(new File(bamFile));
        boolean pairedEnd = reader.iterator().next().getReadPairedFlag();
        ExtendedIterator<SAMRecord> it = EI.wrap(reader.iterator());
        SAMFileHeader header = SamReaderFactory.makeDefault().getFileHeader(new File(referenceBam));
        SAMFileWriter samWriter = new SAMFileWriterFactory().makeBAMWriter(header, false, new File(bamFile.replace(".bam", "_reverted.bam")));
        HashMap<String, SAMRecord[]> pairedReadMap = new HashMap<>();
        Pattern p = Pattern.compile("\\S+#(\\w+)");
        Pattern tagPattern = Pattern.compile("\\*([\\S&&[^;\\*]]+)\\-([\\S&&[^;\\*]]+);");


        if (pairedEnd) {
            p = Pattern.compile("\\S+_#(\\w+)_#(\\w+)");
        }

        for (SAMRecord rec : it.loop()) {
            if (!pairedEnd && toPlusStrand && rec.getReadNegativeStrandFlag()) {
                continue;
            }
            if (rec.getStringAttribute("MD").contains("N")) {
                continue;
            }
            if (rec.getCigarString().contains("N")) {
                continue;
            }
            if (rec.getNotPrimaryAlignmentFlag()) {
                continue;
            }

            if (!pairedEnd) {
                rec = createRecFromUnmappable(rec, p, tagPattern, origGenome, mapGenome, toPlusStrand, pairedEnd);
                samWriter.addAlignment(rec);
            } else {
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
        }

        for (Map.Entry<String, SAMRecord[]> entry : pairedReadMap.entrySet()) {
            SAMRecord[] recs = entry.getValue();
            if (recs == null) {
                continue;
            }
            SAMRecord rec1 = entry.getValue()[0];
            SAMRecord rec2 = entry.getValue()[1];
            if (rec1 != null) {
                rec1 = createRecFromUnmappable(rec1, p, tagPattern, origGenome, mapGenome, toPlusStrand, pairedEnd);
            }
            if (rec2 != null) {
                rec2 = createRecFromUnmappable(rec2, p, tagPattern, origGenome, mapGenome, toPlusStrand, pairedEnd);
            }
            completePairedReads(rec1, rec2);

            if (rec1 != null) {
                samWriter.addAlignment(rec1);
            }
            if (rec2 != null) {
                samWriter.addAlignment(rec2);
            }
        }

        samWriter.close();

        System.out.println("-noInduced:" + noInduced);
    }

    public static SAMRecord createRecFromUnmappable(SAMRecord record, Pattern p, Pattern tagPattern, Genomic origGenome, Genomic mapGenome, boolean toPlusStrand, boolean pairedEnd) {
        SAMRecord rec = record.deepCopy();
        try {
            Matcher m = p.matcher(rec.getReadName());
            m.find();
            String readSequence = "";
            if (!pairedEnd) {
                readSequence = m.group(1);
            } else if (rec.getFirstOfPairFlag()) {
                readSequence = m.group(1);
            } else if (rec.getSecondOfPairFlag()) {
                readSequence = m.group(2);
            }
            rec.setReadString(readSequence);

            Matcher tagMatcher = tagPattern.matcher(rec.getReadName());
            while (tagMatcher.find()) {
                rec.setAttribute(tagMatcher.group(1), tagMatcher.group(2));
            }


            String[] inducedPos = inducedPositions(rec, origGenome, mapGenome, toPlusStrand);
            rec.setReadNegativeStrandFlag(inducedPos[1].equals("-"));
            rec.setAlignmentStart(Integer.valueOf(inducedPos[2]) + 1);
            rec.setReferenceName(inducedPos[0]);
            String origSequence = origGenome.getSequence(Chromosome.obtain(inducedPos[0] + inducedPos[1]), BamUtils.getArrayGenomicRegion(rec)).toString();

            if (rec.getReadNegativeStrandFlag()) {
                origSequence = SequenceUtils.getDnaReverseComplement(origSequence);
                rec.setReadString(SequenceUtils.getDnaReverseComplement(readSequence));
                Cigar newCig = new Cigar();
                Cigar oldCig = rec.getCigar();
                for (int i = oldCig.numCigarElements() - 1; i >= 0; i--) {
                    newCig.add(oldCig.getCigarElement(i));
                }
                rec.setCigar(newCig);
            }

            calculateMdAndNmTags(rec, origSequence.getBytes(StandardCharsets.UTF_8), rec.getAlignmentStart() - 1, true, true);
            //Remove reads with more than a certain amount of MM regarding readlength
            if (rec.getIntegerAttribute("NM") > rec.getReadString().length() * 0.75) {
                rec.setReadUnmappedFlag(true);
            }
            //Put negative strand read to unmapped, when its mate is already on negative strand
            if (pairedEnd && rec.getReadNegativeStrandFlag() && rec.getMateNegativeStrandFlag()) {
                rec.setReadUnmappedFlag(true);
            }

        } catch (NullPointerException | IndexOutOfBoundsException | IllegalArgumentException e) {
            rec.setReadUnmappedFlag(true);
            noInduced++;
        }
        return rec;
    }

    public static void completePairedReads(SAMRecord r1, SAMRecord r2) {
        if (r1 == null && r2 == null) {
            return;
        } else if (r1 != null && r2 == null) {
            r1.setMateAlignmentStart(0);
            r1.setMateReferenceName("*");
            r1.setMateUnmappedFlag(true);
        } else if (r1 == null && r2 != null) {
            r2.setMateAlignmentStart(0);
            r2.setMateReferenceName("*");
            r2.setMateUnmappedFlag(true);
        } else if (r1 != null && r2 != null) {
            if (r1.getReadUnmappedFlag()) {
                r2.setMateAlignmentStart(0);
                r2.setMateReferenceName("*");
                r2.setMateUnmappedFlag(true);
            } else {
                r2.setMateAlignmentStart(r1.getAlignmentStart());
                r2.setMateReferenceName(r1.getReferenceName());
                r2.setMateUnmappedFlag(false);
            }
            if (r2.getReadUnmappedFlag()) {
                r1.setMateAlignmentStart(0);
                r1.setMateReferenceName("*");
                r1.setMateUnmappedFlag(true);
            } else {
                r1.setMateAlignmentStart(r2.getAlignmentStart());
                r1.setMateReferenceName(r2.getReferenceName());
                r1.setMateUnmappedFlag(false);
            }
        }
    }

    public static void createMetadata(String prefix) {
        try {
            CenteredDiskIntervalTreeStorage<DefaultAlignedReadsData> out = new CenteredDiskIntervalTreeStorage<>(prefix + "_rescued.cit", DefaultAlignedReadsData.class);
            out.setMetaData(DynamicObject.from("conditions", DynamicObject.arrayOfObjects("name", prefix)));
        } catch (IOException e) {
            e.printStackTrace();
        }
    }

    public static void mergeBAMFilesInCIT(String bamFile1, String bamFile2, String prefix) {
        SamReader reader1 = SamReaderFactory.makeDefault().open(new File(bamFile1));
        SamReader reader2 = SamReaderFactory.makeDefault().open(new File(bamFile2));

        ExtendedIterator<SAMRecord> it1 = EI.wrap(reader1.iterator());
        ExtendedIterator<SAMRecord> it2 = EI.wrap(reader2.iterator());

        try {
            CenteredDiskIntervalTreeStorage<DefaultAlignedReadsData> out = new CenteredDiskIntervalTreeStorage<>(bamFile1.replace(".bam", "_merged.cit"), DefaultAlignedReadsData.class);
            out.setMetaData(DynamicObject.from("conditions", DynamicObject.arrayOfObjects("name", prefix)));
            IterateIntoSink<ImmutableReferenceGenomicRegion<DefaultAlignedReadsData>> sink = new IterateIntoSink<>(r -> out.fill(r));

            for (SAMRecord rec : it1.loop()) {
                try {
                    if (rec.getReadUnmappedFlag()) {
                        continue;
                    }
                    int[] cumNum = new int[1];
                    cumNum[0] = 1;
                    FactoryGenomicRegion fac = BamUtils.getFactoryGenomicRegion(rec, cumNum, false, false, null);
                    fac.add(rec, 0);
                    String strand = "+";
                    if (rec.getReadNegativeStrandFlag()) {
                        strand = "-";
                    }
                    ImmutableReferenceGenomicRegion reg = new ImmutableReferenceGenomicRegion(Chromosome.obtain(rec.getReferenceName() + strand), fac, fac.create());
                    sink.put(reg);
                } catch (SAMFormatException ex) {
                    continue;
                }
            }

            for (SAMRecord rec : it2.loop()) {
                try {
                    if (rec.getReadUnmappedFlag()) {
                        continue;
                    }
                    int[] cumNum = new int[1];
                    cumNum[0] = 1;
                    FactoryGenomicRegion fac = BamUtils.getFactoryGenomicRegion(rec, cumNum, false, false, null);
                    fac.add(rec, 0);
                    String strand = "+";
                    if (rec.getReadNegativeStrandFlag()) {
                        strand = "-";
                    }
                    ImmutableReferenceGenomicRegion reg = new ImmutableReferenceGenomicRegion(Chromosome.obtain(rec.getReferenceName() + strand), fac, fac.create());
                    sink.put(reg);
                } catch (SAMFormatException ex) {
                    continue;
                }
            }
            sink.finish();
        } catch (IOException | InterruptedException e) {
            System.out.println(e.getMessage());
        }
    }


}

class IndexEntry {

    int pseudoPosStart;
    int pseudoPosEnd;
    int origPosStart;
    int origPosEnd;

    public IndexEntry(int pseudoPosStart, int pseudoPosEnd, int origPosStart, int origPosEnd) {
        this.pseudoPosStart = pseudoPosStart;
        this.pseudoPosEnd = pseudoPosEnd;
        this.origPosStart = origPosStart;
        this.origPosEnd = origPosEnd;
    }

    public int getOrigPosStart() {
        return origPosStart;
    }

    public int getOrigPosEnd() {
        return origPosEnd;
    }

    public int getPseudoPosEnd() {
        return pseudoPosEnd;
    }

    public int getPseudoPosStart() {
        return pseudoPosStart;
    }

    @Override
    public String toString() {
        return "IndexEntry{" +
                "pseudoPosStart=" + pseudoPosStart +
                ", pseudoPosEnd=" + pseudoPosEnd +
                ", origPosStart=" + origPosStart +
                ", origPosEnd=" + origPosEnd +
                '}';
    }
}

