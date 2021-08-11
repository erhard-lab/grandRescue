package executables;

import gedi.app.Gedi;
import gedi.core.genomic.Genomic;

import java.nio.file.Paths;
import java.util.ArrayList;

import static gedi.javapipeline.createRescueBash.createRescueBash;

public class ReadRescue {

    public static void main(String[] args) {
        Gedi.startup(false);

        boolean writeAll = false;

        String origGenome = null;
        String pseudoGenome = null;

        String origMapped = "";
        String pseudoStarIndex = "";
        String tmpDir = Paths.get("").toString()+"tmp";
        String prefix = "";
        ArrayList<String> tags = new ArrayList<>();


        int i;
        for (i=0; i<args.length; i++) {
            if (args[i].equals("-all"))
                writeAll = true;
            else if (args[i].equals("-genome")) {
                ArrayList<String> gnames = new ArrayList<>();
                i = checkMultiParam(args, ++i, gnames);
                origGenome = gnames.get(0);
            }
            else if(args[i].equals("-pseudogenome")) {
                ArrayList<String> gnames = new ArrayList<>();
                i = checkMultiParam(args, ++i, gnames);
                pseudoGenome = gnames.get(0);
            }
            else if(args[i].equals("-origmaps")){
                ArrayList<String> gnames = new ArrayList<>();
                i = checkMultiParam(args, ++i, gnames);
                origMapped = gnames.get(0);
            }
            else if(args[i].equals("-pseudoSTAR")){
                ArrayList<String> gnames = new ArrayList<>();
                i = checkMultiParam(args, ++i, gnames);
                pseudoStarIndex = gnames.get(0);
            }
            else if(args[i].equals("-tmp")){
                ArrayList<String> gnames = new ArrayList<>();
                i = checkMultiParam(args, ++i, gnames);
                tmpDir = gnames.get(0);
            }
            else if(args[i].equals("-prefix")){
                ArrayList<String> gnames = new ArrayList<>();
                i = checkMultiParam(args, ++i, gnames);
                prefix = gnames.get(0);
            }
            else if(args[i].equals("-tags")){
                ArrayList<String> tagnames = new ArrayList<>();
                i = checkMultiParam(args, ++i, tagnames);
                tags.addAll(tagnames);
            }else if(args[i].equals("-h")){
                usage();
                return;
            }
            else
                break;
        }

        if(prefix == "" || origGenome == "" || pseudoGenome == "" || pseudoStarIndex == "" || origMapped == ""){
            System.out.println("Necessary infomation: -prefix, -genome, -pseudogenome, -origmaps, -pseudoSTAR");

            usage();
            System.exit(1);
        }

        createRescueBash(writeAll, origGenome, pseudoGenome, origMapped, pseudoStarIndex, tmpDir, prefix, tags);

    }

    private static String checkParam(String[] args, int index) {
        if (index>=args.length || args[index].startsWith("-")) throw new RuntimeException("Missing argument for "+args[index-1]);
        return args[index];
    }
    private static int checkMultiParam(String[] args, int index, ArrayList<String> re) {
        while (index<args.length && !args[index].startsWith("-"))
            re.add(args[index++]);
        return index-1;
    }

    private static void usage() {
        System.out.println("\nRescue unmappable reads of a single file in 4sU rna-seq experiments via mapping of unmappable reads to a pseudo genome and subsequent backtracking of new mappings to real genome.\n");
        System.out.println("\nReadRescue [-genome] [-pseudogenome] [-pseudoSTAR] [-origmaps] [-prefix]\n\n -genome The gedi indexed genome (e.g. h.ens90, m.ens90...) \n -pseudogenome The gedi indexed pseudo genome of the original genome\n -pseudoSTAR The directory of the STAR index of the pseudo genome\n -origmaps Absolute path to the original bam-file with all mapped reads \n -prefix Use Prefix of the -origmaps bam-file\n\n");
    }
}
