package executables;


import gedi.bam.tools.BamUtils;
import gedi.centeredDiskIntervalTree.CenteredDiskIntervalTreeStorage;
import gedi.core.data.annotation.Transcript;
import gedi.core.data.reads.DefaultAlignedReadsData;
import gedi.core.data.reads.ReadCountMode;
import gedi.core.genomic.Genomic;
import gedi.core.reference.Chromosome;
import gedi.core.region.*;
import gedi.region.bam.FactoryGenomicRegion;
import gedi.statistics.Rtbeta;
import gedi.util.ParseUtils;
import gedi.util.SequenceUtils;
import gedi.util.StringUtils;
import gedi.util.dynamic.DynamicObject;
import gedi.util.functions.EI;
import gedi.util.functions.ExtendedIterator;
import gedi.util.functions.IterateIntoSink;
import gedi.util.math.stat.RandomNumbers;
import htsjdk.samtools.*;
import java.io.*;
import java.nio.file.Paths;
import java.util.*;
import java.util.regex.Matcher;
import java.util.regex.Pattern;


import static gedi.javapipeline.createRescueBash.createRescueBash;
import static gedi.javapipeline.createRescuePipe.getPrefix;
import static gedi.util.BAMUtils.*;
import static gedi.util.MMSimulationUtil.*;

public class MismatchSimulation implements Runnable {


    String genomeString;
    String SAMpath;
    String NTRpath;
    double ex;
    double s;
    boolean rtBetaDist;
    boolean uniqueReads = true;
    boolean noMito = true;

    public static void main(String[] args) {

        String genomeString = "m.ens102";
        String SAMpath = "/home/kevin/Desktop/PhD/gedi/Mismatch_Simulation/res/MCMV_SiLAM/MCMV.4sU.bam";
        String NTRpath = "/home/kevin/Desktop/PhD/gedi/Mismatch_Simulation/res/MCMV_SiLAM/bulk_NTRMap.csv";
/*        String SAMpath = "/home/kevin/Desktop/PhD/gedi/Mismatch_Simulation/res/Spt6/Spt6+.no4sU.bam";
        String NTRpath = "/home/kevin/Desktop/PhD/gedi/Mismatch_Simulation/res/Spt6/geneNTRMap.csv";*/
        boolean rtBetaDist = false;
        boolean onlyUniqueReads = true;
        boolean noMito = true;

        double[] exP = new double[]{0,0.25};
        double[] shapes = new double[]{0.25};

        //pseudoMapping("/home/kevin/Desktop/PhD/gedi/Mismatch_Simulation/res/Spt6/checkPseudoMapping/run8_25p_wholeGenome/25p_unmappable.bam", "h.ens90", "h.ens90_oneLine", true);
        //System.out.print(countGlobalTCMM("/home/kevin/Desktop/PhD/gedi/Mismatch_Simulation/res/Spt6/checkPseudoMapping/run7_25p_chrom22/25p_unmapped.fastq"));
        //samOutputFromPseudoMapping("/home/kevin/Desktop/PhD/gedi/Mismatch_Simulation/res/Spt6/checkPseudoMapping/run8_25p_wholeGenome/25p_unmappable.bam","/home/kevin/Desktop/PhD/gedi/Mismatch_Simulation/res/Spt6/checkPseudoMapping/run8_25p_wholeGenome/25p.bam", Genomic.get("h.ens90"), Genomic.get("h.ens90_pseudoGenomeT2C"), true);
        //mergeBAMFilesInCIT("/home/kevin/Desktop/PhD/gedi/Mismatch_Simulation/res/Spt6/checkPseudoMapping/run8_25p_wholeGenome/25p.bam","/home/kevin/Desktop/PhD/gedi/Mismatch_Simulation/res/Spt6/checkPseudoMapping/run8_25p_wholeGenome/25p_unmappable_reverted.bam");
        System.out.println(new Date(System.currentTimeMillis()));
        Genomic.get(genomeString).getTranscripts();
        //simulateExperiment_checkMapping("h.ens90", "/home/kevin/Desktop/PhD/gedi/Mismatch_Simulation/res/Spt6/Spt6+.no4sU.bam", NTRpath.replace("MCMV_SiLAM", "Spt6"), 0.25, true, true);

        for (double ex : exP) {
            for (double s : shapes) {
                new Thread(new MismatchSimulation(genomeString, SAMpath, NTRpath, ex, s, rtBetaDist, onlyUniqueReads, noMito)).start();
            }
        }






        //extractUnmappedReadsToFastq("/home/kevin/Desktop/PhD/gedi/Mismatch_Simulation/res/MCMV_SiLAM/pseudoMapping/MCMV.0p.bam", false, true);
        //samOutputFromPseudoMapping("/home/kevin/Desktop/PhD/gedi/Mismatch_Simulation/res/MCMV_SiLAM/pseudoMapping/MCMV.25p_pseudo.bam","/home/kevin/Desktop/PhD/gedi/Mismatch_Simulation/res/MCMV_SiLAM/pseudoMapping/MCMV.no4sU.bam", Genomic.get("m.ens102"), Genomic.get("m.ens102_pseudoGenomeT2C"), true);
        //mergeBAMFilesInCIT("/home/kevin/Desktop/PhD/gedi/Mismatch_Simulation/res/MCMV_SiLAM/pseudoMapping/MCMV.25p.bam","/home/kevin/Desktop/PhD/gedi/Mismatch_Simulation/res/MCMV_SiLAM/pseudoMapping/MCMV.25p_pseudo_reverted.bam", "MCMV.no4sU");


        //String bamFile = "/home/kevin/Desktop/PhD/gedi/Mismatch_Simulation/res/Spt6_fullrun_no4sU/checkMappingT_toGenes/T2CMapping/Spt6+.25p.unmapped.T2C.bam";
        //String fastqFile = "/home/kevin/Desktop/PhD/gedi/Mismatch_Simulation/res/Spt6_fullrun_no4sU/checkMappingT_toGenes/origMapping/25p.fastq";
        //samOutputFromPseudoMapping("/home/kevin/Desktop/PhD/gedi/Mismatch_Simulation/res/Spt6/checkPseudoMapping/run10_complete+No4sU/mappedToPseudo/no4sU_unmappable.bam","/home/kevin/Desktop/PhD/gedi/Mismatch_Simulation/res/Spt6/checkPseudoMapping/run10_complete+No4sU/mappedToHuman/Spt6+.no4sU.bam", Genomic.get("h.ens90"), Genomic.get("h.ens90_pseudoGenomeT2C"), true);
        //mergeBAMFilesInCIT("/home/kevin/Desktop/PhD/gedi/Mismatch_Simulation/res/Spt6/checkPseudoMapping/run10_complete+No4sU/mappedToPseudo/no4sU_unmappable_reverted.bam","/home/kevin/Desktop/PhD/gedi/Mismatch_Simulation/res/Spt6/checkPseudoMapping/run10_complete+No4sU/mappedToHuman/Spt6+.no4sU.bam");

        String[] files = new String[]{"no4sU.A","no4sU.B","15min.A","15min.B","30min.A","30min.B","60min.A","60min.B","90min.A","90min.B","120min.A","120min.B"};
        files = new String[]{"no4sU.A"};

        for(String name : files){
            //extractUnmappedReadsToFastq("/home/kevin/Desktop/PhD/gedi/Mismatch_Simulation/res/nih_kinetics/mapToOriginal/"+name+".bam", false);
            //samOutputFromPseudoMapping("/home/kevin/Desktop/PhD/gedi/Mismatch_Simulation/res/nih_kinetics/mapToPseudo/"+name+"_unmappable.bam","/home/kevin/Desktop/PhD/gedi/Mismatch_Simulation/res/nih_kinetics/mapToOriginal/"+name+".bam", Genomic.get("m.ens102"), Genomic.get("m.ens102_pseudoGenomeT2C"), true);
            //mergeBAMFilesInCIT("/home/kevin/Desktop/PhD/gedi/Mismatch_Simulation/res/nih_kinetics/mapToPseudo/"+name+"_unmappable_reverted.bam","/home/kevin/Desktop/PhD/gedi/Mismatch_Simulation/res/nih_kinetics/mapToOriginal/"+name+".bam");
        }
        //createRescueBash(false, true, "h.ens90", "h.ens90_pseudo", "Spt6+.10p.bam", "/mnt/hilbert/projects/berg/h.ens90_pseudoT2C/", "temp", "Spt6+.10p");
    }

    public MismatchSimulation(String genomeString, String SAMpath, String NTRpath, double ex, double s, boolean rtBetaDist, boolean onlyUniqueReads, boolean noMito) {
        this.genomeString = genomeString;
        this.SAMpath = SAMpath;
        this.NTRpath = NTRpath;
        this.rtBetaDist = rtBetaDist;
        this.uniqueReads = onlyUniqueReads;
        this.noMito = noMito;
        this.ex = ex;
        this.s = s;
    }

    public static void simulateExperiment(String genomeString, String SAMpath, String NTRpath, int exchangePercentage, boolean onlyUniqueReads, boolean noMito) {
        Genomic genomic = Genomic.get(genomeString);
        SamReader samReader = SamReaderFactory.makeDefault().open(new File(SAMpath));
        Iterator<SAMRecord> it = samReader.iterator();
        Map<String, Double> map = getNTRMap(NTRpath);
        long seed = System.currentTimeMillis();
        RandomNumbers rnd = new RandomNumbers(seed);

        try {
            BufferedWriter writer = new BufferedWriter(new FileWriter(new File(SAMpath.replace("Spt6+.no4sU_2.bam", exchangePercentage + "p.fastq"))));
            BufferedWriter logWriter = new BufferedWriter(new FileWriter(new File(SAMpath.replace(".bam", ".log")), true));
            logWriter.append("Seed " + exchangePercentage + "p File: " + seed + "\n");


            //CIT-File
            CenteredDiskIntervalTreeStorage<DefaultAlignedReadsData> out = new CenteredDiskIntervalTreeStorage<>("/home/kevin/Desktop/PhD/gedi/Mismatch_Simulation/res/Spt6_fullrun_no4sU/Spt6+." + exchangePercentage + "p.cit", DefaultAlignedReadsData.class);
            IterateIntoSink<ImmutableReferenceGenomicRegion<DefaultAlignedReadsData>> sink = new IterateIntoSink<>(r -> out.fill(r));
            while (it.hasNext()) {
                SAMRecord record = it.next();
                String strand = "+";
                String gene = "";
                int start = record.getAlignmentStart() - 1;
                int end = record.getAlignmentEnd();
                String chrom = record.getReferenceName();
                String sequence = record.getReadString();


                int[] cumNum = new int[1];
                cumNum[0] = 1;
                FactoryGenomicRegion fac = BamUtils.getFactoryGenomicRegion(record, cumNum, false, false, null);
                fac.add(record, 0);
                if (onlyUniqueReads && record.getIntegerAttribute("NH") != 1) {
                    continue;
                }

                if (chrom.contains("ERCC")) {
                    continue;
                }

                if (noMito && chrom.contains("MT")) {
                    continue;
                }


                if (record.getReadNegativeStrandFlag()) {
                    strand = "-";
                    sequence = SequenceUtils.getDnaReverseComplement(sequence);
                }

                Chromosome chromosome = Chromosome.obtain("" + chrom + strand);
                GenomicRegion region = new ArrayGenomicRegion(start, end);
                ExtendedIterator<ImmutableReferenceGenomicRegion<Transcript>> ei = genomic.getTranscripts().ei(chromosome, region);
                ImmutableReferenceGenomicRegion<Transcript> current = null;
                if (ei.hasNext()) {
                    current = ei.first();
                    gene = current.getData().getGeneId();
                } else {
                    writer.append("@" + record.getReadName() + "\n");
                    writer.append(sequence + "\n");
                    writer.append("+\n");
                    writer.append(record.getBaseQualityString() + "\n");

                    ImmutableReferenceGenomicRegion reg = new ImmutableReferenceGenomicRegion(Chromosome.obtain(chrom + strand), fac, fac.create());
                    sink.put(reg);

                    continue;
                }

                try {
                    double ntrRatio = map.get(gene);
                    if (rnd.getBool(ntrRatio)) {
                        char[] nucleotides = sequence.toCharArray();
                        String cigarString = cigarToString(record.getCigarString());
                        if (record.getReadNegativeStrandFlag()) {
                            cigarString = StringUtils.reverse(cigarString).toString();
                        }

                        int seqCounter = 0;
                        int cigarCounter = 0;
                        for (int i = 0; i < fac.getTotalLength(); i++) {
                            if (cigarString.charAt(cigarCounter) == 'D' || cigarString.charAt(cigarCounter) == 'N') {
                                cigarCounter++;
                                continue;
                            } else if (cigarString.charAt(cigarCounter) == 'I') {
                                i--;
                                seqCounter++;
                                cigarCounter++;
                                continue;
                            } else if (cigarString.charAt(cigarCounter) == 'S') {
                                i--;
                                seqCounter++;
                                cigarCounter++;
                                continue;
                            }
                            if (nucleotides[seqCounter] != 'T' || !rnd.getBool((double) exchangePercentage / 100)) {
                                seqCounter++;
                                cigarCounter++;
                                continue;
                            }

                            if (exchangePercentage == 0) {
                                throw new IllegalStateException("0 exchange!");
                            }

                            ImmutableReferenceGenomicRegion reg = new ImmutableReferenceGenomicRegion(Chromosome.obtain(chrom + strand), fac.asRegion());
                            String facSeq = genomic.getSequence(reg).toString();
                            int mappedPos = reg.mapMaybeOutSide(i);
                            if (mappedPos < 0) {
                                seqCounter++;
                                cigarCounter++;
                                continue;
                            }
                            int inducedPos = -1;
                            try {
                                inducedPos = current.induceMaybeOutside(mappedPos);
                            } catch (IllegalArgumentException e) {
                                seqCounter++;
                                cigarCounter++;
                                continue;
                            }
                            String transcrSeq = genomic.getSequence(current).toString();

                            if (inducedPos < 0 || inducedPos >= transcrSeq.length()) {
                                seqCounter++;
                                cigarCounter++;
                                continue;
                            }

                            if (facSeq.charAt(i) == 'T' && transcrSeq.charAt(inducedPos) == 'T' && !fac.create().toString().contains("" + i)) {
                                nucleotides[seqCounter] = 'C';
                                fac.getFactory().addMismatch(i, 'T', 'C', false);
                            }

                            seqCounter++;
                            cigarCounter++;
                        }
                        String replaced = String.valueOf(nucleotides);

                        writer.append("@" + record.getReadName() + "\n");
                        writer.append(replaced + "\n");
                        writer.append("+\n");
                        writer.append(record.getBaseQualityString() + "\n");
                        ImmutableReferenceGenomicRegion reg = new ImmutableReferenceGenomicRegion(Chromosome.obtain(chrom + strand), fac, fac.create());
                        sink.put(reg);
                    } else {
                        writer.append("@" + record.getReadName() + "\n");
                        writer.append(sequence + "\n");
                        writer.append("+\n");
                        writer.append(record.getBaseQualityString() + "\n");

                        ImmutableReferenceGenomicRegion reg = new ImmutableReferenceGenomicRegion(Chromosome.obtain(chrom + strand), fac, fac.create());
                        sink.put(reg);
                    }
                } catch (NullPointerException e) {
                    writer.append("@" + record.getReadName() + "\n");
                    writer.append(sequence + "\n");
                    writer.append("+\n");
                    writer.append(record.getBaseQualityString() + "\n");

                    ImmutableReferenceGenomicRegion reg = new ImmutableReferenceGenomicRegion(Chromosome.obtain(chrom + strand), fac, fac.create());
                    sink.put(reg);
                    continue;
                }

            }
            sink.finish();
            logWriter.close();
            writer.close();
        } catch (IOException | InterruptedException e) {
            System.out.println(e.getMessage());
        }
    }

    public static void simulateExperiment_checkMapping(String genomeString, String SAMpath, String NTRpath, double exchangePercentage, boolean onlyUniqueReads, boolean noMito) {
        Genomic genomic = Genomic.get(genomeString);
        SamReader samReader = SamReaderFactory.makeDefault().open(new File(SAMpath));
        Iterator<SAMRecord> it = EI.wrap(samReader.iterator()).progress();
        Map<String, Double> map = getNTRMap(NTRpath);
        long seed = System.currentTimeMillis();
        RandomNumbers rnd = new RandomNumbers(seed);
        int mismatchCount = 0;
        int tcount = 0;
        int tcmm = 0;

        try {
            BufferedWriter writer = new BufferedWriter(new FileWriter(new File(SAMpath.replace("Spt6+.no4sU.bam", exchangePercentage + "p.fastq"))));
            BufferedWriter logWriter = new BufferedWriter(new FileWriter(new File(SAMpath.replace(".bam", ".log")), true));
            logWriter.append("Seed " + exchangePercentage + "p File: " + seed + "\n");


            //CIT-File
            CenteredDiskIntervalTreeStorage<DefaultAlignedReadsData> out = new CenteredDiskIntervalTreeStorage<>("/home/kevin/Desktop/PhD/gedi/Mismatch_Simulation/res/Spt6/Spt6+." + exchangePercentage + "p.cit", DefaultAlignedReadsData.class);
            IterateIntoSink<ImmutableReferenceGenomicRegion<DefaultAlignedReadsData>> sink = new IterateIntoSink<>(r -> out.fill(r));
            while (it.hasNext()) {
                SAMRecord record = it.next();
                String strand = "+";
                String gene = "";
                int start = record.getAlignmentStart() - 1;
                int end = record.getAlignmentEnd();
                String chrom = record.getReferenceName();
                String sequence = record.getReadString();


                int[] cumNum = new int[1];
                cumNum[0] = 1;
                FactoryGenomicRegion fac = BamUtils.getFactoryGenomicRegion(record, cumNum, false, false, null);
                fac.add(record, 0);
                if (onlyUniqueReads && record.getIntegerAttribute("NH") != 1) {
                    continue;
                }

                if (chrom.contains("ERCC")) {
                    continue;
                }

                if (noMito && chrom.contains("MT")) {
                    continue;
                }

                if(record.getIntegerAttribute("nM") >= 5){
                    continue;
                }


                if (record.getReadNegativeStrandFlag()) {
                    strand = "-";
                    sequence = SequenceUtils.getDnaReverseComplement(sequence);
                }

                Chromosome chromosome = Chromosome.obtain("" + chrom + strand);
                GenomicRegion region = new ArrayGenomicRegion(start, end);
                ExtendedIterator<ImmutableReferenceGenomicRegion<Transcript>> ei = genomic.getTranscripts().ei(chromosome, region);
                ImmutableReferenceGenomicRegion<Transcript> current = null;
                String recordName = chrom;
                if (strand == "+") {
                    recordName = recordName + "_plus_" + start + "_" + end;
                } else {
                    recordName = recordName + "_minus_" + start + "_" + end;
                }
                if (ei.hasNext()) {
                    current = ei.first();
                    gene = current.getData().getGeneId();
                } else {
                    writer.append("@" + recordName + "\n");
                    writer.append(sequence + "\n");
                    writer.append("+\n");
                    writer.append(record.getBaseQualityString() + "\n");

                    ImmutableReferenceGenomicRegion reg = new ImmutableReferenceGenomicRegion(Chromosome.obtain(chrom + strand), fac, fac.create());
                    sink.put(reg);

                    continue;
                }

                try {
                    double ntrRatio = map.get(gene);
                    if (rnd.getBool(ntrRatio)) {
                        char[] nucleotides = sequence.toCharArray();


                        int seqCounter = 0;
                        int cigarCounter = 0;
                        int currentTmismatches = 0;
                        for (int i = 0; i < fac.getTotalLength(); i++) {

                            String cigarString = record.getCigarString();
                            if(record.getReadNegativeStrandFlag()){
                                cigarString = reverseCigar(cigarString);
                            }
                            String getCigarNumber = "";
                            int cigarNumber = 0;
                            char type = cigarString.charAt(cigarCounter);

                            while (Character.isDigit(type)) {
                                getCigarNumber = getCigarNumber + type;
                                cigarCounter++;
                                type = cigarString.charAt(cigarCounter);
                            }
                            cigarCounter++;

                            cigarNumber = Integer.valueOf(getCigarNumber);


                            if (type == 'D') {
                                i = i + cigarNumber - 1;
                                continue;
                            } else if (type == 'I' || type == 'S') {
                                i--;
                                seqCounter = seqCounter + cigarNumber;
                                continue;
                            } else if (type == 'N'){
                                i--;
                                continue;
                            }

                            for (int c = 0; c < cigarNumber; c++) {
                                if(nucleotides[seqCounter] == 'T'){
                                    tcount++;
                                }
                                if (nucleotides[seqCounter] != 'T' || !rnd.getBool(exchangePercentage)) {
                                    seqCounter++;
                                    i++;
                                    continue;
                                }


                                if (exchangePercentage == 0) {
                                    throw new IllegalStateException("0 exchange!");
                                }

                                ImmutableReferenceGenomicRegion reg = new ImmutableReferenceGenomicRegion(Chromosome.obtain(chrom + strand), fac.asRegion());
                                String facSeq = genomic.getSequence(reg).toString();
                                int mappedPos = reg.mapMaybeOutSide(i);
                                if (mappedPos < 0) {
                                    seqCounter++;
                                    i++;
                                    continue;
                                }

                                int inducedPos = -1;
                                try {
                                    inducedPos = current.induceMaybeOutside(mappedPos);
                                } catch (IllegalArgumentException e) {
                                    seqCounter++;
                                    i++;
                                    continue;
                                }

                                String transcrSeq = genomic.getSequence(current).toString();

                                if (inducedPos < 0 || inducedPos >= transcrSeq.length()) {
                                    seqCounter++;
                                    i++;
                                    continue;
                                }


                                if (facSeq.charAt(i) == 'T' && transcrSeq.charAt(inducedPos) == 'T' && !fac.create().toString().contains("" + i)) {
                                    nucleotides[seqCounter] = 'C';
                                    fac.getFactory().addMismatch(i, 'T', 'C', false);
                                    currentTmismatches++;
                                    tcmm++;
                                }

                                seqCounter++;
                                i++;

                            }
                            i--;
                        }
                        String replaced = String.valueOf(nucleotides);

                        mismatchCount = mismatchCount + currentTmismatches;
                        writer.append("@" + recordName + "_" + currentTmismatches + "\n");
                        writer.append(replaced + "\n");
                        writer.append("+\n");
                        writer.append(record.getBaseQualityString() + "\n");
                        ImmutableReferenceGenomicRegion reg = new ImmutableReferenceGenomicRegion(Chromosome.obtain(chrom + strand), fac, fac.create());
                        sink.put(reg);
                    } else {
                        writer.append("@" + recordName + "\n");
                        writer.append(sequence + "\n");
                        writer.append("+\n");
                        writer.append(record.getBaseQualityString() + "\n");

                        ImmutableReferenceGenomicRegion reg = new ImmutableReferenceGenomicRegion(Chromosome.obtain(chrom + strand), fac, fac.create());
                        sink.put(reg);
                    }
                } catch (NullPointerException e) {
                    writer.append("@" + recordName + "\n");
                    writer.append(sequence + "\n");
                    writer.append("+\n");
                    writer.append(record.getBaseQualityString() + "\n");

                    ImmutableReferenceGenomicRegion reg = new ImmutableReferenceGenomicRegion(Chromosome.obtain(chrom + strand), fac, fac.create());
                    sink.put(reg);
                    continue;
                }

            }
            System.out.println("Sum t->c mismatches: " + mismatchCount);
            sink.finish();
            logWriter.close();
            writer.close();
        } catch (IOException | InterruptedException e) {
            System.out.println(e.getMessage());
        }
    }

    public static void simulateExperiment_rtBeta(String genomeString, String SAMpath, String NTRpath, boolean rtBetaDist, boolean onlyUniqueReads, boolean noMito, double lower, double upper, double shape) {
        if (lower > 1 || upper > 1) {
            throw new IllegalArgumentException("lower or upper cannot be > 1");
        }

        int tcmm = 0;
        int tcount = 0;
        double ntrs = 0;
        int counter = 0;

        double s1 = Math.exp(shape);
        double s2 = Math.exp(-shape);

        Genomic genomic = Genomic.get(genomeString);
        SamReader samReader = SamReaderFactory.makeDefault().open(new File(SAMpath));
        Iterator<SAMRecord> it = EI.wrap(samReader.iterator()).progress();
        Map<String, Double> map = getNTRMap(NTRpath);
        long seed = System.currentTimeMillis();
        RandomNumbers rnd = new RandomNumbers(seed);
        String shapeTag = "";
        if(rtBetaDist){
            shapeTag = "_"+(int)shape;
        }


        try {
            BufferedWriter writer = new BufferedWriter(new FileWriter(new File(SAMpath.replace(".4sU.bam", "_"+(int)(upper*100) + "p" + shapeTag + ".fastq"))));
            BufferedWriter logWriter = new BufferedWriter(new FileWriter(new File(SAMpath.replace(".bam", ".log")), true));
            logWriter.append("Seed " + upper + "p" + shapeTag + "_" + "File: " + seed + "\n");
            BufferedWriter exchangePwriter = new BufferedWriter(new FileWriter(new File(SAMpath.replace(".4sU.bam", "_"+(int)(upper*100) + "p" + shapeTag + ".exRates.tsv"))));


            //CIT-File
            CenteredDiskIntervalTreeStorage<DefaultAlignedReadsData> out = new CenteredDiskIntervalTreeStorage<>(SAMpath.replace(".4sU.bam", "_" + (int)(upper*100) + "p" + shapeTag + ".cit"), DefaultAlignedReadsData.class);
            IterateIntoSink<ImmutableReferenceGenomicRegion<DefaultAlignedReadsData>> sink = new IterateIntoSink<>(r -> out.fill(r));

            while (it.hasNext()) {
                SAMRecord record = it.next();
                String strand = "+";
                String gene = "";
                int start = record.getAlignmentStart() - 1;
                int end = record.getAlignmentEnd();
                String chrom = record.getReferenceName();
                String sequence = record.getReadString();


                int[] cumNum = new int[1];
                cumNum[0] = 1;
                FactoryGenomicRegion fac = BamUtils.getFactoryGenomicRegion(record, cumNum, false, false, null);
                fac.add(record, 0);
                if (onlyUniqueReads && record.getIntegerAttribute("NH") != 1) {
                    continue;
                }

                if (chrom.contains("ERCC")) {
                    continue;
                }

                if (noMito && chrom.contains("MT")) {
                    continue;
                }


                if (record.getReadNegativeStrandFlag()) {
                    strand = "-";
                    sequence = SequenceUtils.getDnaReverseComplement(sequence);
                }

                Chromosome chromosome = Chromosome.obtain("" + chrom + strand);
                GenomicRegion region = new ArrayGenomicRegion(start, end);
                ExtendedIterator<ImmutableReferenceGenomicRegion<Transcript>> ei = genomic.getTranscripts().ei(chromosome, region);
                ImmutableReferenceGenomicRegion<Transcript> current = null;
                if (ei.hasNext()) {
                    current = ei.first();
                    gene = current.getData().getGeneId();
                } else {
                    writer.append("@" + record.getReadName() + "\n");
                    writer.append(sequence + "\n");
                    writer.append("+\n");
                    writer.append(record.getBaseQualityString() + "\n");

                    ImmutableReferenceGenomicRegion reg = new ImmutableReferenceGenomicRegion(Chromosome.obtain(chrom + strand), fac, fac.create());
                    sink.put(reg);

                    continue;
                }

                try {
                    double ntrRatio = map.get(gene);
                    counter++;
                    if (!Double.isNaN(ntrRatio)) {
                        ntrs = ntrs + ntrRatio;
                    }
                    if (rnd.getBool(ntrRatio)) {
                        double mmCount = 0.0;
                        double exchangePercentage = upper;
                        if(rtBetaDist) {
                            exchangePercentage = Rtbeta.rtbeta(1, lower, upper, s1, s2)[0];
                        }
                        exchangePwriter.append(gene + "\t" + exchangePercentage + "\t" + ntrRatio + "\n");

                        char[] nucleotides = sequence.toCharArray();

                        int seqCounter = 0;
                        int cigarCounter = 0;

                        for (int i = 0; i < fac.getTotalLength(); i++) {

                            String cigarString = record.getCigarString();
                            if(record.getReadNegativeStrandFlag()){
                                cigarString = reverseCigar(cigarString);
                            }
                            String getCigarNumber = "";
                            int cigarNumber = 0;
                            char type = cigarString.charAt(cigarCounter);

                            while (Character.isDigit(type)) {
                                getCigarNumber = getCigarNumber + type;
                                cigarCounter++;
                                type = cigarString.charAt(cigarCounter);
                            }
                            cigarCounter++;

                            cigarNumber = Integer.valueOf(getCigarNumber);


                            if (type == 'D') {
                                i = i + cigarNumber - 1;
                                continue;
                            } else if (type == 'I' || type == 'S') {
                                i--;
                                seqCounter = seqCounter + cigarNumber;
                                continue;
                            } else if (type == 'N'){
                                i--;
                                continue;
                            }

                            for (int c = 0; c < cigarNumber; c++) {
                                if(nucleotides[seqCounter] == 'T'){
                                    tcount++;
                                }
                                if (nucleotides[seqCounter] != 'T' || !rnd.getBool(exchangePercentage)) {
                                    seqCounter++;
                                    i++;
                                    continue;
                                }

                                ImmutableReferenceGenomicRegion reg = new ImmutableReferenceGenomicRegion(Chromosome.obtain(chrom + strand), fac.asRegion());
                                String facSeq = genomic.getSequence(reg).toString();
                                int mappedPos = reg.mapMaybeOutSide(i);
                                if (mappedPos < 0) {
                                    seqCounter++;
                                    i++;
                                    continue;
                                }

                                int inducedPos = -1;
                                try {
                                    inducedPos = current.induceMaybeOutside(mappedPos);
                                } catch (IllegalArgumentException e) {
                                    seqCounter++;
                                    i++;
                                    continue;
                                }

                                String transcrSeq = genomic.getSequence(current).toString();

                                if (inducedPos < 0 || inducedPos >= transcrSeq.length()) {
                                    seqCounter++;
                                    i++;
                                    continue;
                                }


                                if (facSeq.charAt(i) == 'T' && transcrSeq.charAt(inducedPos) == 'T' && !fac.create().toString().contains("" + i)) {
                                    nucleotides[seqCounter] = 'C';
                                    fac.getFactory().addMismatch(i, 'T', 'C', false);
                                    tcmm++;
                                    mmCount++;
                                }

                                seqCounter++;
                                i++;

                            }
                            i--;
                        }
                        String replaced = String.valueOf(nucleotides);

                        writer.append("@" + record.getReadName() + "\n");
                        writer.append(replaced + "\n");
                        writer.append("+\n");
                        writer.append(record.getBaseQualityString() + "\n");
                        ImmutableReferenceGenomicRegion reg = new ImmutableReferenceGenomicRegion(Chromosome.obtain(chrom + strand), fac, fac.create());
                        //Simulating read dropout in CIT-File due to unmappable reads. The higher the mmCount in a read, the higher the chance for a read to not make it into CIT-File
                        //if(rnd.getBool(1.0-(mmCount/100))){
                            sink.put(reg);
                        //}
                    } else {
                        writer.append("@" + record.getReadName() + "\n");
                        writer.append(sequence + "\n");
                        writer.append("+\n");
                        writer.append(record.getBaseQualityString() + "\n");

                        ImmutableReferenceGenomicRegion reg = new ImmutableReferenceGenomicRegion(Chromosome.obtain(chrom + strand), fac, fac.create());
                        sink.put(reg);
                    }
                } catch (NullPointerException e) {
                    writer.append("@" + record.getReadName() + "\n");
                    writer.append(sequence + "\n");
                    writer.append("+\n");
                    writer.append(record.getBaseQualityString() + "\n");

                    ImmutableReferenceGenomicRegion reg = new ImmutableReferenceGenomicRegion(Chromosome.obtain(chrom + strand), fac, fac.create());
                    sink.put(reg);
                    continue;
                }

            }
            logWriter.append(upper + "p" + shapeTag + " NTR-Ratio: " + (ntrs / counter) + "\n");
            logWriter.append("total Ts: " + tcount + "\n");
            logWriter.append("TCMMs: " + tcmm + "\n");


            sink.finish();
            logWriter.close();
            writer.close();
            exchangePwriter.close();
            System.out.println(new Date(System.currentTimeMillis()));
        } catch (IOException | InterruptedException e) {
            System.out.println(e.getMessage());
        }
    }

    @Override
    public void run() {
        simulateExperiment_rtBeta(genomeString, SAMpath, NTRpath, rtBetaDist, uniqueReads, noMito, 3 * Math.pow(10, -4), ex, s);
    }
}
