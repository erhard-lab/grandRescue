package gedi.util;

import java.io.*;
import java.nio.file.Files;
import java.sql.SQLOutput;
import java.util.*;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

public class createOnlyExonGenome {

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

            if(!map.containsKey(key)){
                map.put(key, new ArrayList<String>());
            }

            map.get(key).add(current);
        }

        return map;
    }

    public static void fastaSeparator(File fasta, String outPath) throws IOException{
        BufferedReader reader = new BufferedReader(new FileReader(fasta));

        String current = reader.readLine();
        String chrom = "";
        BufferedWriter writer = null;
        outPath = outPath + fasta.getAbsolutePath().substring(fasta.getAbsolutePath().lastIndexOf("/"));
        System.out.println(outPath);

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

                File file = new File(outPath.replace(".fasta", "_chrom_"+chrom+".fasta.tmp"));
                file.getParentFile().mkdirs();
                file.deleteOnExit();
                writer = new BufferedWriter(new FileWriter(file));
                current = reader.readLine();
                continue;
            }

            writer.append(current);
            current = reader.readLine();
        }
        writer.close();
    }


    public static String getSeq(String fastaPath, String chrom, int startPos, int endPos) throws IOException{
        BufferedReader reader = new BufferedReader(new FileReader(fastaPath.replace(".fasta", "_chrom_"+chrom+".fasta.tmp")));
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

    public static void createGenomeFiles(String gtfString, String fastaString, String outPath, String from, String to, boolean toPlusStrand){

        System.out.println("Starting to create Genome Files...");

        File gtf = new File(gtfString);
        File fasta = new File(fastaString);

        System.out.println("Files located...");

        try {
            HashMap<String, ArrayList<String>> map = gtfSeparator(gtf);
            outPath = outPath + fasta.getAbsolutePath().substring(fasta.getAbsolutePath().lastIndexOf("/"));
            String fastaPath = outPath.replace(".fasta", "_pseudo.fasta");
            BufferedWriter fastaWriter = new BufferedWriter(new FileWriter(fastaPath));
            BufferedWriter gtfWriter = new BufferedWriter(new FileWriter(fastaPath.replace(".fasta", ".gtf")));
            BufferedWriter indexWriter = new BufferedWriter(new FileWriter(fastaPath.replace(".fasta", ".index")));
            System.out.println(fastaPath);

            gtfWriter.append(getGTFHeader(gtf));
            indexWriter.append("Gene\tpseudoPos\torigPos\n");

            System.out.println("GTF separated into chromosomal records...");
            System.out.println("Start writing FASTA & GTF-Files...");
            int counter = 0;
            for(String key : map.keySet()){
                fastaWriter.append(">"+key+"\n");
                System.out.println("- Chromosome " + key);

                long fastaPointer = 1;
                long lastGeneStart = -1;
                long lastGeneStartFasta = -1;

                ArrayList<String> records = map.get(key);
                for(int i = 0; i < records.size(); i++){
                    counter++;
                    if(counter == 11){
                        counter = 0;
                        //break;
                    }
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
                        String geneSeq = getSeq(outPath, key, Integer.valueOf(m.group(4)), Integer.valueOf(m.group(5)));
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

            System.out.println("Nucleotide conversion of genome ("+from+" -> "+to+") ...");

            NucleotideConversion(new File(fastaPath), from, to);

            System.out.println("Finished createGenomeFiles...");
        }catch (IOException e){
            System.out.println(e.getMessage());
        }
    }

    public static void NucleotideConversion(File fasta, String from, String to){
        try{
            BufferedReader reader = new BufferedReader(new FileReader(fasta));
            BufferedWriter writer = new BufferedWriter(new FileWriter(fasta.getPath().replace(".fasta", "_"+from+"2"+to+".fasta")));
            Files.delete(fasta.toPath());

            String line = reader.readLine();
            while(line != null){
                if(line.startsWith(">")){
                    writer.append(line+"\n");
                } else {
                    writer.append(line.replace(from, to)+"\n");
                }
                line = reader.readLine();
            }
            reader.close();
            writer.close();
        }catch(IOException e){
            System.out.println(e.getMessage());
        }
    }

}
