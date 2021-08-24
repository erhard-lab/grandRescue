package executables;

import gedi.app.Gedi;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;

import static gedi.util.createOnlyExonGenome.*;

public class CreatePseudo {

    public static void main(String[] args) {
        Gedi.startup(false);

        String gtf = "";
        String fasta = "";
        String outPath = "";

        int i;
        for (i=0; i<args.length; i++) {
            if (args[i].equals("-gtf")) {
                gtf = checkParam(args, ++i);
            }
            else if(args[i].equals("-fasta")) {
                fasta = checkParam(args, ++i);
            }
            else if(args[i].equals("-o")){
                outPath = checkParam(args, ++i);
            }else if(args[i].equals("-h")){
                usage();
                return;
            }
            else
                break;
        }

        try {
            fastaSeparator(new File(fasta), outPath);
            createGenomeFiles(gtf, fasta, outPath, true);
        } catch(IOException e){
            e.printStackTrace();
        }
    }

    private static void usage() {
        System.out.println("\nCreate fasta- & gtf-Files of PseudoGenome for ReadRescue-Procedure, replacing all T's with C's and replacing all intronic regions with 100 N spacer.\nSTAR's genomeGenerate on these files might take several days, depending on the size of the genome!\n");
        System.out.println("\nCreatePseudo [-gtf] [-fasta] [-o]\n\n -gtf The gtf-file of the genome \n -fasta The fasta-File of the genome\n -o Output-Folder for the PseudoGenome\n\n");
    }

    private static String checkParam(String[] args, int index) throws IllegalArgumentException {
        if (index>=args.length) throw new IllegalArgumentException("Missing argument for "+args[index-1]);
        return args[index];
    }

}
