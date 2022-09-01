package gedi;

import gedi.core.reference.Strandness;
import gedi.javapipeline.ExtractionParameterSet;
import gedi.util.program.GediProgram;
import gedi.util.program.GediProgramContext;

import java.io.File;
import java.util.ArrayList;

import static gedi.util.ReadExtraction.extractUnmappedReadsToFastq;

public class ExtractUnmappedReads extends GediProgram {
    public ExtractUnmappedReads(ExtractionParameterSet params){
        addInput(params.strandness);
        addInput(params.tags);
        addInput(params.from);
        addInput(params.to);
        addInput(params.f);

        addOutput(params.outFile);
    }

    @Override
    public String execute(GediProgramContext context) {
        Strandness strandness = getParameter(0);
        ArrayList<String> tags = getParameters(1);
        String from = getParameter(2);
        String to = getParameter(3);
        String file = getParameter(4);

        extractUnmappedReadsToFastq(new File(file), tags, strandness, from, to);

        return null;
    }
}
