package gedi.javapipeline;

import gedi.core.genomic.Genomic;
import gedi.core.reference.Strandness;
import gedi.util.program.GediParameter;
import gedi.util.program.GediParameterSet;
import gedi.util.program.parametertypes.*;

import java.io.File;

public class RescuePipeParameterSet extends GediParameterSet {

    public GediParameter<File> paramFile = new GediParameter<File>(this,"${prefix}.param", "File containing the parameters used to call this program", false, new FileParameterType());

    public GediParameter<String> genome = new GediParameter<String>(this,"genome", "The indexed GEDI genome.", true, new StringParameterType());
    public GediParameter<String> pseudogenome = new GediParameter<String>(this,"pseudogenome", "The indexed GEDI pseudo-genome.", false, new StringParameterType());
    public GediParameter<String> prefix = new GediParameter<String>(this,"prefix", "The prefix used for all output files", false, new StringParameterType());
    public GediParameter<String> pseudoSTAR = new GediParameter<String>(this,"pseudoSTAR", "The path to the STAR index of the pseudo genome", false, new StringParameterType());
    public GediParameter<String> tmp = new GediParameter<String>(this,"tmp", "The path for temporary files to be placed", false, new StringParameterType());
    public GediParameter<String> files = new GediParameter<String>(this,"files", "The paths to the mapped 4sU-experiment bam files", true, new StringParameterType());
    public GediParameter<Strandness> strandness = new GediParameter<Strandness>(this,"strandness", "Whether sequencing protocol was stranded (Sensse), strand unspecific (Unspecific), or opposite strand (Antisense).", false, new EnumParameterType<>(Strandness.class), Strandness.Sense,true);

    public GediParameter<File> bashFile = new GediParameter<File>(this,"${prefix}.sh", "Final bash file", false, new FileParameterType());



}
