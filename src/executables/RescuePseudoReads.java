package executables;

import gedi.RescueReads;
import gedi.app.Gedi;
import gedi.javapipeline.RescueParameterSet;
import gedi.util.LogUtils;
import gedi.util.program.CommandLineHandler;
import gedi.util.program.GediProgram;


public class RescuePseudoReads {

    private static String getChangelog() {
        return "0.9.1:\n"
                + "Added -strandness param & handling of Sense/Antisense\n"
                + "added -maxMM param\n"
                + "changed RescuePseudoReads to GediProgram\n\n";
    }

    public static void main(String[] args) {
        Gedi.startup(false, LogUtils.LogMode.Normal,"readRescue");

        RescueParameterSet params = new RescueParameterSet();
        GediProgram pipeline = GediProgram.create("readRescue", new RescueReads(params));
        pipeline.setChangelog(getChangelog());
        GediProgram.run(pipeline, new CommandLineHandler("readRescue", "Rescue unmappable 4sU reads", args));

        /*String origGenome = "";
        String pseudoGenome = "";
        String pseudoMapped = "";
        String origMapped = "";
        int maxMM = 100;
        boolean keepID = false;
        Strandness strandness = Strandness.Sense;

        int i;
        for (i=0; i<args.length; i++) {
            if (args[i].equals("-genome")) {
                ArrayList<String> gnames = new ArrayList<>();
                i = checkMultiParam(args, ++i, gnames);
                origGenome = gnames.get(0);
            }
            else if(args[i].equals("-pseudogenome")) {
                ArrayList<String> gnames = new ArrayList<>();
                i = checkMultiParam(args, ++i, gnames);
                pseudoGenome = gnames.get(0);
            }
            else if (args[i].equals("-pseudomaps")) {
                ArrayList<String> gnames = new ArrayList<>();
                i = checkMultiParam(args, ++i, gnames);
                pseudoMapped = gnames.get(0);
            }
            else if(args[i].equals("-origmaps")){
                ArrayList<String> gnames = new ArrayList<>();
                i = checkMultiParam(args, ++i, gnames);
                origMapped = gnames.get(0);
            }
            else if(args[i].equals("-maxMM")){
                ArrayList<String> mm = new ArrayList<>();
                i = checkMultiParam(args, ++i, mm);
                maxMM = Integer.valueOf(mm.get(0));
            }
            else if(args[i].equals("-strandness")){
                ArrayList<String> strand = new ArrayList<>();
                i = checkMultiParam(args, ++i, strand);
                if(strand.get(0).equals("Sense")){
                    strandness = Strandness.Sense;
                } else if(strand.get(0).equals("Antisense")){
                    strandness = Strandness.Antisense;
                } else {
                    usage();
                    break;
                }
            }
            else if (args[i].equals("-keepID")) {
                keepID = true;
            }
            else if(args[i].equals("-h")){
                usage();
                return;
            }
            else
                break;
        }

        System.out.println("Rescue pseudo-reads...");
        samOutputFromPseudoMapping(pseudoMapped, origMapped, Genomic.get(origGenome), Genomic.get(pseudoGenome), keepID, strandness, maxMM);
        createMetadata(getPrefix(origMapped));*/
    }

/*    private static void usage() {
        System.out.println("\nInduces the positions of reads mapped to pseudo genome in the original genome and merges these rescued reads with the original bam-files. The output is a bam-file with originally mapped reads and rescued reads.\n");
        System.out.println("\nrescuePseudoReads [-genome] [-pseudogenome] [-pseudomaps] [-origmaps] [-strandness] [-maxMM]\n\n -genome The gedi indexed genome (e.g. h.ens90, m.ens90...) \n -pseudogenome The gedi indexed pseudo genome of the original genome\n -origmaps the original bam-file with all mapped reads \n -pseudomaps bam-file of reads mapped to pseudo genome\n -strandness Sense or Antisense\n -maxMM Maximum number of Mismatches allowed \n\n");
    }

    private static int checkMultiParam(String[] args, int index, ArrayList<String> re) {
        while (index<args.length && !args[index].startsWith("-"))
            re.add(args[index++]);
        return index-1;
    }*/

}
