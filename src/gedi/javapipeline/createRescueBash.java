package gedi.javapipeline;

import executables.Pipeline;
import gedi.core.reference.Strandness;
import gedi.util.FileUtils;
import gedi.util.io.text.LineIterator;
import java.io.File;
import java.io.IOException;
import java.net.URL;
import java.nio.charset.Charset;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.ArrayList;

public class createRescueBash {


    public static String createRescueBash(String origGenome, String pseudoGenome, String pseudoStarIndex, String tmpDir, String prefix, ArrayList<String> tags, Strandness strandness, boolean pe, String from, String to, int maxMM, String chrPrefix) {
        String pathString = Paths.get(prefix+".sh").toAbsolutePath().toString();
        try {
            Charset charset = StandardCharsets.UTF_8;
            String name = "readRescue.sh";


            String tagnames = "";
            if(!tags.isEmpty()){
                tagnames = "-tags";
                for(String s : tags){
                    tagnames = tagnames + " "+s;
                }
            }
            String file = "";
            if (!new File(name).exists()) {
                URL res = Pipeline.class.getResource("/resources/"+name);
                if (res!=null) {
                    file = new LineIterator(res.openStream()).concat("\n");
                }
            }
            else {
                file = FileUtils.readAllText(new File(name));
            }

            file = file.replaceAll("\\{tmp}", tmpDir);
            if(!from.equals("")) {
                file = file.replaceAll("\\{from}", "-from " + from);
            } else {
                file = file.replaceAll("\\{from}", "");
            }
            if(!to.equals("")) {
                file = file.replaceAll("\\{to}", "-to " + to);
            } else {
                file = file.replaceAll("\\{to}", "");
            }
            if(pe){
                file = file.replaceAll("\\{pe_unmapped}", "samtools view -b -f 8 -F 4 {prefix}_unmappedMates.bam\n"+
                        "samtools merge unmappedMatesReads.bam {prefix}_unmapped.bam {prefix}_unmappedMates.bam\n"+
                        "samtools sort -n unmappedMatesReads.bam\n"+
                        "mv unmappedMatesReads.bam {prefix}_unmapped.bam\n");
            } else {
                file = file.replaceAll("\\{pe_unmapped}", "");
            }

            if(maxMM!=999) {
                file = file.replaceAll("\\{maxMM}", "-maxMM " + maxMM);
            } else {
                file = file.replaceAll("\\{maxMM}", "");
            }
            if(!chrPrefix.equals("")) {
                file = file.replaceAll("\\{chrPrefix}", "-chrPrefix " + chrPrefix);
            } else {
                file = file.replaceAll("\\{chrPrefix}", "");
            }
            file = file.replaceAll("\\{pseudoStarIndex}", pseudoStarIndex);
            file = file.replaceAll("\\{genome}", origGenome);
            file = file.replaceAll("\\{pseudogenome}", pseudoGenome);
            file = file.replaceAll("\\{files}", prefix+"_unmapped_T2C.fastq");
            file = file.replaceAll("\\{prefix}", prefix);
            file = file.replaceAll("\\{tags}", tagnames);
            file = file.replaceAll("\\{strandness}", strandness.name());



            Files.write(Paths.get(pathString), file.getBytes(charset));
            File finalFile = new File(pathString);
            finalFile.setExecutable(true);
        }catch (IOException e){
            e.printStackTrace();
        }

        return pathString;
    }
}
