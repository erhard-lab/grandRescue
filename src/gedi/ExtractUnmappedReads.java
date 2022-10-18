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
        addInput(params.from);
        addInput(params.to);
        addInput(params.f);

        addOutput(params.outFile);
    }

    @Override
    public String execute(GediProgramContext context) {
        Strandness strandness = getParameter(0);
        String from = getParameter(1);
        String to = getParameter(2);
        String file = getParameter(3);

        extractUnmappedReadsToFastq(new File(file), strandness, from, to);

        return null;
    }
}
