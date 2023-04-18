package executables;

import gedi.app.Gedi;
import gedi.core.genomic.Genomic;
import gedi.core.reference.Strandness;
import gedi.util.LogUtils;

import java.nio.file.Paths;
import java.util.ArrayList;

import static gedi.javapipeline.createRescueBash.createRescueBash;

public class ReadRescue {

    public static void main(String[] args) {
        Gedi.startup(false, LogUtils.LogMode.Normal, "readRescue");

        String origGenome = null;
        String pseudoGenome = null;

        String pseudoStarIndex = "";
        String tmpDir = Paths.get("").toString()+"tmp";
        String from = "";
        String to = "";
        String file = "";
        Strandness strandness = Strandness.Sense;
        boolean pe = false;
        boolean grandslam = false;
        int maxMM = 999;
        String chrPrefix = "";


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
            else if(args[i].equals("-pseudoSTAR")){
                ArrayList<String> gnames = new ArrayList<>();
                i = checkMultiParam(args, ++i, gnames);
                pseudoStarIndex = gnames.get(0);
            }
            else if(args[i].equals("-f")) {
                ArrayList<String> gnames = new ArrayList<>();
                i = checkMultiParam(args, ++i, gnames);
                file = gnames.get(0);
            }
            else if(args[i].equals("-tmp")){
                ArrayList<String> gnames = new ArrayList<>();
                i = checkMultiParam(args, ++i, gnames);
                tmpDir = gnames.get(0);
            }
            else if(args[i].equals("-prefix")){
                ArrayList<String> gnames = new ArrayList<>();
                i = checkMultiParam(args, ++i, gnames);
                file = gnames.get(0);
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
            else if(args[i].equals("-pe")){
                pe = true;
            }
            else if(args[i].equals("-grandslam")){
                grandslam = true;
            }
            else if(args[i].equals("-from")){
                ArrayList<String> gnames = new ArrayList<>();
                i = checkMultiParam(args, ++i, gnames);
                from = gnames.get(0);
            }
            else if(args[i].equals("-to")){
                ArrayList<String> gnames = new ArrayList<>();
                i = checkMultiParam(args, ++i, gnames);
                to = gnames.get(0);
            }
            else if(args[i].equals("-maxMM")){
                ArrayList<String> gnames = new ArrayList<>();
                i = checkMultiParam(args, ++i, gnames);
                maxMM = Integer.valueOf(gnames.get(0));
            }
            else if(args[i].equals("-chrPrefix")){
                ArrayList<String> gnames = new ArrayList<>();
                i = checkMultiParam(args, ++i, gnames);
                chrPrefix = gnames.get(0);
            }
            else if(args[i].equals("-h")){
                usage();
                return;
            }
            else
                break;
        }

        if(file == "" || origGenome == "" || pseudoGenome == ""){
            System.out.println("Necessary infomation: -f, -genome, -pseudogenome, -pseudoSTAR");

            usage();
            System.exit(1);
        }
        if(pseudoStarIndex == ""){
            System.out.println("IMPORTANT: You have not specified a STAR-index, make sure to edit the bash-file to add a proper call to a mapping tool.");
        }

        createRescueBash(origGenome, pseudoGenome, pseudoStarIndex, tmpDir, file.replace(".bam", ""), strandness, pe, grandslam, from, to, maxMM, chrPrefix);

    }


    private static int checkMultiParam(String[] args, int index, ArrayList<String> re) {
        while (index<args.length && !args[index].startsWith("-"))
            re.add(args[index++]);
        return index-1;
    }

    private static void usage() {
        System.out.println("\nRescue unmappable reads of a single file in 4sU rna-seq experiments via mapping of unmappable reads to a pseudo genome and subsequent backtracking of new mappings to real genome.\n");
        System.out.println("\nReadRescue [-genome] [-pseudogenome] [-pseudoSTAR] [-f] [-strandness] [-pe] [-from] [-to] [-maxMM] [-chrPrefix]\n\n -genome The gedi indexed genome (e.g. h.ens90, m.ens90...) \n -pseudogenome The gedi indexed pseudo genome of the original genome\n -pseudoSTAR The directory of the STAR index of the pseudo genome\n -origmaps Absolute path to the original bam-file with all mapped reads \n -prefix Use Prefix of the -origmaps bam-file\n\n");
    }
}
