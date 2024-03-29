package gedi.javapipeline;

import gedi.core.genomic.Genomic;
import gedi.core.reference.Strandness;
import gedi.util.program.GediParameter;
import gedi.util.program.GediParameterSet;
import gedi.util.program.parametertypes.*;

import java.io.File;

public class RescueParameterSet extends GediParameterSet {
    public GediParameter<Genomic> genome = new GediParameter<Genomic>(this,"genome", "The indexed GEDI genome.", true, new GenomicParameterType());
    public GediParameter<Genomic> pseudogenome = new GediParameter<Genomic>(this,"pseudogenome", "The indexed GEDI pseudogenome.", false, new GenomicParameterType());
    public GediParameter<String> origmaps = new GediParameter<String>(this,"origmaps", "The original STAR-mapped bam-file", false, new StringParameterType());
    public GediParameter<String> pseudomaps = new GediParameter<String>(this,"pseudomaps", "The bam-file mapped to pseudogenome", false, new StringParameterType());
    public GediParameter<String> idMap = new GediParameter<String>(this,"idMap", "Path to .idMap-file generated by ExtractUnmappedReads", false, new StringParameterType());
    public GediParameter<Integer> maxMM = new GediParameter<Integer>(this,"maxMM", "The maximum number of mismatches allowed after rescue", false, new IntParameterType(), 100, true);
    public GediParameter<Strandness> strandness = new GediParameter<Strandness>(this,"strandness", "Whether sequencing protocol was stranded (Sense) or opposite strand (Antisense).", false, new EnumParameterType<>(Strandness.class), Strandness.Sense,true);
    public GediParameter<Boolean> keepID = new GediParameter<Boolean>(this,"keepID", "Keep modified read ID from rescue procedure", false, new BooleanParameterType(), false, true);
    public GediParameter<Boolean> noMito = new GediParameter<Boolean>(this,"noMito", "Don't rescue mitochondrial reads", false, new BooleanParameterType(), false, true);
    //public GediParameter<String> prefix = new GediParameter<String>(this,"pseudomaps", "The bam-file mapped to pseudogenome", false, new StringParameterType(),  RescueReads.getPrefix(origmaps.toString()), true);
    public GediParameter<File> outFile = new GediParameter<File>(this,"rescued.bam", "BAM-file containing all originally mapped reads and rescued reads", false, new FileParameterType());
    //public GediParameter<File> paramFile = new GediParameter<File>(this,"${prefix}.param", "File containing the parameters used to call PRICE", false, new FileParameterType());
    public GediParameter<String> chrPrefix = new GediParameter<String>(this,"chrPrefix", "A Prefix/Pattern for chromosome names, * replaces chromosome number. E.g. chr*", false, new StringParameterType(), "");
}
