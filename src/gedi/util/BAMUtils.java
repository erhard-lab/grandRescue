package gedi.util;

import gedi.bam.tools.BamUtils;
import gedi.centeredDiskIntervalTree.CenteredDiskIntervalTreeStorage;
import gedi.core.data.annotation.Transcript;
import gedi.core.data.reads.DefaultAlignedReadsData;
import gedi.core.genomic.Genomic;
import gedi.core.reference.Chromosome;
import gedi.core.reference.Strandness;
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
import static gedi.util.SequenceUtil.reverseComplement;

public class BAMUtils {

    static int noInduced = 0;

    public static String[] inducedPositions(SAMRecord rec, Genomic origGenome, Genomic mapGenome) {
        String[] inducedPositions = new String[3];

        Chromosome mapChromosome = Chromosome.obtain("" + rec.getReferenceName() + "+");
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

        int inducedStart = -1;

        if (origRGR.getReference().isMinus()) {
            inducedStart = origRGR.map(mapRGR.induce(rec.getAlignmentEnd() - 1));
            Cigar newCig = new Cigar();
            Cigar oldCig = rec.getCigar();
            for (int i = oldCig.numCigarElements() - 1; i >= 0; i--) {
                newCig.add(oldCig.getCigarElement(i));
            }
            rec.setCigar(newCig);

            rec.setReadNegativeStrandFlag(!rec.getReadNegativeStrandFlag());
            if(rec.getReadPairedFlag()) {
                rec.setMateNegativeStrandFlag(!rec.getMateNegativeStrandFlag());
            }

        } else {
            inducedStart = origRGR.map(mapRGR.induce(rec.getAlignmentStart() - 1));
        }

        inducedPositions[0] = origRGR.getReference().getName();

        inducedPositions[1] = "+";
        if (rec.getReadNegativeStrandFlag()) {
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

    public static void samOutputFromPseudoMapping(String bamFile, String referenceBam, Genomic origGenome, Genomic mapGenome, boolean keepID, Strandness strandness, int maxMM) {
        SamReader reader = SamReaderFactory.makeDefault().open(new File(bamFile));
        SAMRecordIterator pairedIT = reader.iterator();
        SAMFileHeader header = SamReaderFactory.makeDefault().getFileHeader(new File(referenceBam));
        SAMFileWriter samWriter = new SAMFileWriterFactory().makeBAMWriter(header, false, new File(bamFile.replace(".bam", "_reverted.bam")));
        if (!pairedIT.hasNext()) {
            System.out.println("File " + bamFile + " has 0 reads");
            samWriter.close();
            return;
        }
        boolean pairedEnd = pairedIT.next().getReadPairedFlag();
        pairedIT.close();
        ExtendedIterator<SAMRecord> it = EI.wrap(reader.iterator());
        HashMap<String, SAMRecord[]> pairedReadMap = new HashMap<>();
        Pattern p = Pattern.compile("\\S+#(\\w+)_");
        Pattern tagPattern = Pattern.compile("\\*([\\S&&[^;\\*]]+)\\~([\\S&&[^;\\*]]+);");


        if (pairedEnd) {
            p = Pattern.compile("\\S+_#(\\w+)_#(\\w+)_");
        }

        for (SAMRecord rec : it.loop()) {
            if (!pairedEnd && !strandness.equals(Strandness.Antisense) && rec.getReadNegativeStrandFlag()) {
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
                rec = createRecFromUnmappable(rec, p, tagPattern, origGenome, mapGenome, pairedEnd, keepID, maxMM);
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
                rec1 = createRecFromUnmappable(rec1, p, tagPattern, origGenome, mapGenome, pairedEnd, keepID, maxMM);
            }
            if (rec2 != null) {
                rec2 = createRecFromUnmappable(rec2, p, tagPattern, origGenome, mapGenome, pairedEnd, keepID, maxMM);
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
    }

    public static SAMRecord createRecFromUnmappable(SAMRecord record, Pattern p, Pattern tagPattern, Genomic origGenome, Genomic mapGenome, boolean pairedEnd, boolean keepID, int maxMM) {
        SAMRecord rec = record.deepCopy();

        if(maxMM == 100){
            maxMM = (int) (rec.getReadString().length() * 0.75);
        }
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


            String[] inducedPos = inducedPositions(rec, origGenome, mapGenome);
            rec.setReadNegativeStrandFlag(inducedPos[1].equals("-"));
            rec.setAlignmentStart(Integer.valueOf(inducedPos[2]) + 1);
            rec.setReferenceName(inducedPos[0]);
            String origSequence = origGenome.getSequence(Chromosome.obtain(inducedPos[0] + "+"), BamUtils.getArrayGenomicRegion(rec)).toString();

            if (rec.getReadNegativeStrandFlag()) {
                rec.setReadString(reverseComplement(rec.getReadString()));
            }

            rec.setAttribute("nM", null);
            rec.setAttribute("NM", null);

            calculateMdAndNmTags(rec, origSequence.getBytes(StandardCharsets.UTF_8), rec.getAlignmentStart() - 1, true, true);


            //Remove reads with more than a certain amount of MM regarding readlength
            if (rec.getIntegerAttribute("NM") > maxMM) {
                rec.setReadUnmappedFlag(true);
            }

        } catch (NullPointerException | IndexOutOfBoundsException | IllegalArgumentException e) {
            e.printStackTrace();
            rec.setReadUnmappedFlag(true);
            noInduced++;
        }
        if (!keepID) {
            rec.setReadName(rec.getReadName().substring(0, rec.getReadName().indexOf("_")));
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

    public static String getPrefix(String path) {
        return path.substring(path.lastIndexOf("/") + 1, path.lastIndexOf(".bam"));
    }

}
