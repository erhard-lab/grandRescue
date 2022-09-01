package gedi.javapipeline;

import gedi.core.reference.Strandness;
import gedi.util.program.GediParameter;
import gedi.util.program.GediParameterSet;
import gedi.util.program.parametertypes.EnumParameterType;
import gedi.util.program.parametertypes.FileParameterType;
import gedi.util.program.parametertypes.StringParameterType;

import java.io.File;

public class ExtractionParameterSet extends GediParameterSet {
    public GediParameter<Strandness> strandness = new GediParameter<Strandness>(this,"strandness", "Whether sequencing protocol was stranded (Sense) or opposite strand (Antisense).", false, new EnumParameterType<>(Strandness.class), Strandness.Sense,true);
    public GediParameter<String> tags = new GediParameter<String>(this,"tags", "BAM-File tags to keep after rescue.", true, new StringParameterType(), "", true);
    public GediParameter<String> from = new GediParameter<String>(this,"from", "Nucleotide to be converted", false, new StringParameterType(), "T", true);
    public GediParameter<String> to = new GediParameter<String>(this,"to", "Nucleotide to be converted to", false, new StringParameterType(), "C", true);
    public GediParameter<String> f = new GediParameter<String>(this,"f", "The original STAR-mapped bam-file", false, new StringParameterType(), false);

    public GediParameter<File> outFile = new GediParameter<File>(this, "ExtractReads.param", "BAM-file containing all originally mapped reads and rescued reads", false, new FileParameterType());

}
