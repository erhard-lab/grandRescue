package executables;

import gedi.app.Gedi;
import gedi.javapipeline.RescuePipeParameterSet;
import gedi.javapipeline.createRescuePipe;
import gedi.util.LogUtils;
import gedi.util.program.CommandLineHandler;
import gedi.util.program.GediProgram;

public class ReadRescuePipe {

    public static void main(String[] args) {
        Gedi.startup(false, LogUtils.LogMode.Normal,"readRescue");

        RescuePipeParameterSet params = new RescuePipeParameterSet();
        GediProgram pipeline = GediProgram.create("ReadRescue-Pipe", new createRescuePipe(params));
        GediProgram.run(pipeline, params.paramFile, new CommandLineHandler("ReadRescuePipe", "Rescue unmappable 4sU reads", args));
    }

    private static void usage() {
        System.out.println("\nA pipeline to rescue unmappable reads in 4sU rna-seq experiments via mapping of unmappable reads to a pseudo genome and subsequent backtracking of new mappings to real genome.\n");
        System.out.println("\nCreates a bash-file {prefix}_rescue.sh that runs the whole pipeline and returns the {prefix}.cit file containing the whole dataset.\n");
        System.out.println("\nReadRescuePipe [-threads] [-genome] [-pseudogenome] [-pseudoSTAR] [-prefix] [-tmp] [-sense] [-files]  \n\n -threads Number of threads\n -genome The gedi indexed genome (e.g. h.ens90, m.ens90...) \n -pseudogenome The gedi indexed pseudo genome of the original genome\n -pseudoSTAR The directory of the STAR index of the pseudo genome \n -prefix Use Prefix of the -origmaps bam-file\n\n");
    }

}
