package executables;

import gedi.app.Gedi;

import java.util.ArrayList;

import static gedi.util.MMSimulationUtil.extractUnmappedReadsToFastq;

public class ExtractUnmappedReads {

    public static void main(String[] args) {
        Gedi.startup(false);

        boolean writeAll = false;
        boolean compress = false;
        boolean pairedEnd = false;

        int i;
        for (i=0; i<args.length; i++) {
            if (args[i].equals("-all"))
                writeAll = true;
            else if(args[i].equals("-c")) {
                compress = true;
            }
            else if(args[i].equals("-pairedEnd")) {
                pairedEnd = true;
            }
            else if (args[i].equals("-h")) {
                usage();
                return;
            }
            else
                break;
        }

        if (i+1!=args.length) {
            usage();
            System.exit(1);
        }

        extractUnmappedReadsToFastq(args[i], writeAll, compress, pairedEnd);
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
        System.out.println("\nA method to extract unmapped reads from a bam-file and convert all T's to C's for 4sU-read rescue via PseudoMapping. (Bam-files only contain unmapped reads if the STAR parameter outSAMmapped is set to Within)\n");
        System.out.println("\nextractUnmappedReads [-all] [-c] [-pe] <input>\n\n -all Output fastq-files for mapped reads and non-T-to-C conversed reads too\n -c compress intermediate files from 'all' param to fastq.gz \n -pe Data is paired end");
    }


}
