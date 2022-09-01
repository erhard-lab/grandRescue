package executables;

import gedi.ExtractUnmappedReads;
import gedi.app.Gedi;
import gedi.javapipeline.ExtractionParameterSet;
import gedi.util.LogUtils;
import gedi.util.program.CommandLineHandler;
import gedi.util.program.GediProgram;

public class ExtractReads {
    private static String getChangelog() {
        return "1.0:\n"
                + "Extraction of unmapped reads\n\n";
    }

    public static void main(String[] args) {
        Gedi.startup(false, LogUtils.LogMode.Normal,"ExtractReads");

        ExtractionParameterSet params = new ExtractionParameterSet();
        GediProgram pipeline = GediProgram.create("ExtractReads", new ExtractUnmappedReads(params));
        pipeline.setChangelog(getChangelog());
        GediProgram.run(pipeline, new CommandLineHandler("ExtractReads", "Extraction of unmapped reads", args));
    }
}
