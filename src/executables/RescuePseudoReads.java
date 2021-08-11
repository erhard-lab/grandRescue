package executables;

import gedi.app.Gedi;
import gedi.core.genomic.Genomic;

import java.util.ArrayList;

import static gedi.util.BAMUtils.*;

public class RescuePseudoReads {


    public static void main(String[] args) {
        Gedi.startup(false);

        String origGenome = "";
        String pseudoGenome = "";
        String pseudoMapped = "";
        String origMapped = "";
        String prefix = "";
        boolean mapToPlusStrand = true;

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
            else if (args[i].equals("-prefix")) {
                ArrayList<String> gnames = new ArrayList<>();
                i = checkMultiParam(args, ++i, gnames);
                prefix = gnames.get(0);
            }
            else if(args[i].equals("-h")){
                usage();
                return;
            }
            else
                break;
        }

        System.out.println("Rescue pseudo-reads...");
        samOutputFromPseudoMapping(pseudoMapped, origMapped, Genomic.get(origGenome), Genomic.get(pseudoGenome), mapToPlusStrand);
        //mergeBAMFilesInCIT(pseudoMapped.replace(".bam", "_reverted.bam"), origMapped, prefix);
        createMetadata(prefix);
    }

    private static void usage() {
        System.out.println("\nInduces the positions of reads mapped to pseudo genome in the original genome and merges these rescued reads with the original bam-files. The output is a bam-file with originally mapped reads and rescued reads.\n");
        System.out.println("\nrescuePseudoReads [-genome] [-pseudogenome] [-pseudomaps] [-origmaps]\n\n -genome The gedi indexed genome (e.g. h.ens90, m.ens90...) \n -pseudogenome The gedi indexed pseudo genome of the original genome\n -origmaps the original bam-file with all mapped reads \n -pseudomaps bam-file of reads mapped to pseudo genome\n\n");
    }

    private static int checkMultiParam(String[] args, int index, ArrayList<String> re) {
        while (index<args.length && !args[index].startsWith("-"))
            re.add(args[index++]);
        return index-1;
    }

}
