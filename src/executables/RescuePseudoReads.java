package executables;

import gedi.RescueReads;
import gedi.app.Gedi;
import gedi.javapipeline.RescueParameterSet;
import gedi.util.LogUtils;
import gedi.util.program.CommandLineHandler;
import gedi.util.program.GediProgram;


public class RescuePseudoReads {

    private static String getChangelog() {
        return "0.9.2:\n"
                + "Added handling of multiple genomes as input\n\n";
    }

    public static void main(String[] args) {
        Gedi.startup(false, LogUtils.LogMode.Normal,"readRescue");

        RescueParameterSet params = new RescueParameterSet();
        GediProgram pipeline = GediProgram.create("readRescue", new RescueReads(params));
        pipeline.setChangelog(getChangelog());
        GediProgram.run(pipeline, new CommandLineHandler("readRescue", "Rescue unmappable 4sU reads", args));
    }

}
