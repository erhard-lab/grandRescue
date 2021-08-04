package gedi.util;

import java.io.*;
import java.sql.SQLOutput;
import java.util.*;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

public class createOnlyExonGenome {

    public static void main(String[] args) {
        try {
            String gtf = "/home/kevin/Desktop/PhD/gedi/Mismatch_Simulation/res/m.ens102/mus_musculus.102.gtf";
            String fasta = "/home/kevin/Desktop/PhD/gedi/Mismatch_Simulation/res/m.ens102/mus_musculus.102.fasta";

            fastaSeparator(new File(fasta));
            createGenomeFiles(gtf, fasta, true);
            System.out.println(getSeq(new File(fasta),"20", 145415, 145751));
            //System.out.println(getSeq("20",13092, 13428, new File("/home/kevin/Desktop/PhD/gedi/Mismatch_Simulation/res/h.ens90/homo_sapiens.90_genes_corr_T2C.fasta")));


        }catch (IOException e){
            System.out.println(e.getMessage());
        }
    }


    public static HashMap<String, ArrayList<String>> gtfSeparator(File gtf) throws IOException{
        HashMap<String, ArrayList<String>> map = new HashMap<>();

        Scanner sc = new Scanner(gtf);
        while(sc.hasNextLine()){
            String current = sc.nextLine();
            if(current.startsWith("#")){
                continue;
            }

            Pattern p = Pattern.compile("(\\S+)\t");
            Matcher m = p.matcher(current);
            m.find();
            String key = m.group(1);

            /*System.out.println("--");
            System.out.println(current);
            System.out.println(key);*/

            if(!map.containsKey(key)){
                map.put(key, new ArrayList<String>());
            }

            map.get(key).add(current);
        }

        return map;
    }

    public static HashMap<String, String> fastaSeparator(File fasta) throws IOException{
        HashMap<String, String> map = new HashMap<>();
        BufferedReader reader = new BufferedReader(new FileReader(fasta));

        String current = reader.readLine();
        String chrom = "";
        BufferedWriter writer = null;
        while(current != null){
            if(current.startsWith(">")){
                if(writer!=null){
                    writer.close();
                }
                Pattern p = Pattern.compile("(\\S+)\\s");
                Matcher m = p.matcher(current);
                m.find();
                chrom = m.group(1).replace(">","");
                System.out.println(chrom);

                writer = new BufferedWriter(new FileWriter(fasta.getPath().replace(".fasta", "_chrom_"+chrom+".fasta")));
                current = reader.readLine();
                continue;
            }
            //writer.append(current+"\n");
            writer.append(current);
            current = reader.readLine();
        }
        writer.close();

        return map;
    }


    public static String getSeq(File fasta, String chrom, int startPos, int endPos) throws IOException{
        BufferedReader reader = new BufferedReader(new FileReader(fasta.getPath().replace(".fasta", "_chrom_"+chrom+".fasta")));
        reader.skip(startPos-1);
        String output = "";
        for(int i = 0; i < endPos-startPos+1; i++){
            char currChar = (char)reader.read();
            output = output+currChar;
        }
        return output;
    }


    public static String getGTFHeader(File gtf) throws IOException{
        String header = "";

        Scanner sc = new Scanner(gtf);
        while(sc.hasNextLine()){
            String current = sc.nextLine();
            if(current.startsWith("#")){
                if(current.contains("!genome-build ")){
                    current = current + " modified, genes";
                }
                header = header+current+"\n";
            } else {
                break;
            }
        }
        return header;
    }

    /**
     * Add \n after certain amount of chars
     * @param fasta
     * @param separateAfter
     * @throws IOException
     */
    public static void correctFASTA(File fasta, int separateAfter) throws IOException{
        System.out.println("Correcting FASTA-File to " + separateAfter + " nucleotides per line...");
        BufferedReader reader = new BufferedReader(new FileReader(fasta));
        BufferedWriter writer = new BufferedWriter(new FileWriter(fasta.getPath().replace(".fasta", "_corr.fasta")));

        int counter = 0;
        boolean notFirstLine = false;

        String currLine = reader.readLine();
        while(currLine != null){
            if(currLine.startsWith(">")){
                if(notFirstLine){
                    counter = 0;
                    writer.append("\n");
                } else {
                    notFirstLine = true;
                }
                writer.append(currLine+"\n");
                currLine = reader.readLine();
                continue;
            }
            for(int i = 0; i < currLine.length(); i++){
                writer.append(currLine.charAt(i));
                counter++;
                if(counter == separateAfter){
                    counter = 0;
                    writer.append("\n");
                }
            }
            currLine = reader.readLine();

        }

        reader.close();
        writer.close();
    }

    public static void createGenomeFiles(String gtfString, String fastaString, boolean toPlusStrand){

        System.out.println("Starting to create Genome Files...");

        File gtf = new File(gtfString);
        File fasta = new File(fastaString);

        System.out.println("Files located...");

        try {
            HashMap<String, ArrayList<String>> map = gtfSeparator(gtf);
            String fastaPath = fasta.getPath().replace(".fasta", "_pseudo.fasta");
            BufferedWriter fastaWriter = new BufferedWriter(new FileWriter(fastaPath));
            BufferedWriter gtfWriter = new BufferedWriter(new FileWriter(fastaPath.replace(".fasta", ".gtf")));
            BufferedWriter indexWriter = new BufferedWriter(new FileWriter(fastaPath.replace(".fasta", ".index")));

            gtfWriter.append(getGTFHeader(gtf));
            indexWriter.append("Gene\tpseudoPos\torigPos\n");

            System.out.println("GTF separated into chromosomal records...");
            System.out.println("Start writing FASTA & GTF-Files...");
            for(String key : map.keySet()){
                fastaWriter.append(">"+key+"\n");
                System.out.println("- Chromosome " + key);

                long fastaPointer = 1;
                long lastGeneStart = -1;
                long lastGeneStartFasta = -1;

                ArrayList<String> records = map.get(key);
                for(int i = 0; i < records.size(); i++){
                    String rec = records.get(i);

                    //Groups:
                    //1 Chromosome                          4 startPos  7 strand (+,-,.)
                    //2 Source                              5 endPos    8 frame (0,1,2..)
                    //3 Feature(gene, exon, transcript..)   6 score     9 attributes

                    Pattern p = Pattern.compile("(\\S+)\\t(\\S+)\\t(\\S+)\\t(\\d+)\\t(\\d+)\\t(\\S+)\\t(\\S)\\t(\\S)\\t(.+)");
                    Matcher m = p.matcher(rec);
                    m.find();

                    if(!m.group(1).equals(key)){
                        throw new IllegalStateException("Key != chrom in gtf: " + key + " <-> " + m.group(1));
                    }

                    if(!m.group(3).equals("gene") && !m.group(3).equals("exon") && !m.group(3).equals("CDS")){
                        continue;
                    }


                    if(m.group(3).equals("gene")){
                        fastaWriter.append("NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN");
                        fastaPointer += 100;
                        lastGeneStart = Integer.valueOf(m.group(4));

                        lastGeneStartFasta = fastaPointer;
                        String geneSeq = getSeq(fasta, key, Integer.valueOf(m.group(4)), Integer.valueOf(m.group(5)));
                        if(m.group(7).equals("-") && toPlusStrand){
                            geneSeq = SequenceUtils.getDnaReverseComplement(geneSeq);
                        }

                        fastaWriter.append(geneSeq);
                        if(toPlusStrand) {
                            gtfWriter.append(m.group(1) + "\t" + m.group(2) + "\t" + m.group(3) + "\t" + fastaPointer + "\t" + (fastaPointer + geneSeq.length() - 1) + "\t" + m.group(6) + "\t" + "+" + "\t" + m.group(8) + "\t" + m.group(9) + "\n");
                        } else {
                            gtfWriter.append(m.group(1) + "\t" + m.group(2) + "\t" + m.group(3) + "\t" + fastaPointer + "\t" + (fastaPointer + geneSeq.length() - 1) + "\t" + m.group(6) + "\t" + m.group(7) + "\t" + m.group(8) + "\t" + m.group(9) + "\n");
                        }
                        Pattern g = Pattern.compile("gene_id \"(\\S+)\"");
                        Matcher gm = g.matcher(m.group(9));
                        gm.find();
                        indexWriter.append(gm.group(1)+"\t"+fastaPointer+"/"+(fastaPointer+geneSeq.length()-1)+"\t"+m.group(4)+"/"+m.group(5)+"\n");

                        fastaPointer = fastaPointer + geneSeq.length();

                    } else {
                        long exonStart = lastGeneStartFasta+(Integer.valueOf(m.group(4))-lastGeneStart);
                        long exonEnd = exonStart+(Integer.valueOf(m.group(5))-Integer.valueOf(m.group(4)));
                        if(toPlusStrand) {
                            gtfWriter.append(m.group(1) + "\t" + m.group(2) + "\t" + m.group(3) + "\t" + exonStart + "\t" + exonEnd + "\t" + m.group(6) + "\t" + "+" + "\t" + m.group(8) + "\t" + m.group(9) + "\n");
                        } else {
                            gtfWriter.append(m.group(1) + "\t" + m.group(2) + "\t" + m.group(3) + "\t" + exonStart + "\t" + exonEnd + "\t" + m.group(6) + "\t" + m.group(7) + "\t" + m.group(8) + "\t" + m.group(9) + "\n");

                        }
                    }



                }
                fastaWriter.append("\n");

            }
            fastaWriter.close();
            gtfWriter.close();
            indexWriter.close();

/*            System.out.println("FASTA & GTF-Files complete...");

            correctFASTA(new File(fastaPath), 60);*/

            System.out.println("T2C conversion of genome...");

            T2CConversion(new File(fastaPath));

            System.out.println("Finished createGenomeFiles...");
        }catch (IOException e){
            System.out.println(e.getMessage());
        }
    }

    public static void T2CConversion(File fasta){
        try{
            BufferedReader reader = new BufferedReader(new FileReader(fasta));
            BufferedWriter writer = new BufferedWriter(new FileWriter(fasta.getPath().replace(".fasta", "_T2C.fasta")));
            int lineCount = 1;

            String line = reader.readLine();
            while(line != null){
                if(line.startsWith(">")){
                    writer.append(line+"\n");
                } else {
                    writer.append(line.replace("T", "C")+"\n");
                }
                lineCount++;
                line = reader.readLine();
            }
            reader.close();
            writer.close();
        }catch(IOException e){
            System.out.println(e.getMessage());
        }
    }

}
