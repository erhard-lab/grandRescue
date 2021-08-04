package gedi.javapipeline;

import executables.Pipeline;
import gedi.util.FileUtils;
import gedi.util.io.text.LineIterator;
import java.io.File;
import java.io.IOException;
import java.net.URL;
import java.nio.charset.Charset;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Paths;

public class createRescueBash {


    public static String createRescueBash(boolean writeAll, String origGenome, String pseudoGenome, String origMapped, String pseudoStarIndex, String tmpDir, String prefix, boolean pairedEnd) {
        String pathString = Paths.get(prefix+".sh").toAbsolutePath().toString();
        try {
            Charset charset = StandardCharsets.UTF_8;

            String name = "readRescue.sh";
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

            if(!writeAll){
                file = file.replaceAll("-all", "");
            }
            if(!pairedEnd){
                file = file.replaceAll("-pairedEnd", "");
            }
            file = file.replaceAll("\\{tmp}", tmpDir);
            file = file.replaceAll("\\{bampath}", origMapped);
            file = file.replaceAll("\\{pseudoStarIndex}", pseudoStarIndex);
            file = file.replaceAll("\\{genome}", origGenome);
            file = file.replaceAll("\\{pseudogenome}", pseudoGenome);
            if(!pairedEnd) {
                file = file.replaceAll("\\{files}", prefix+"_unmapped_T2C.fastq");
                file = file.replaceAll("\\{samtools}", "samtools view -b -F 4 {prefix}.bam > {prefix}_final.bam");
            } else {
                file = file.replaceAll("\\{files}", prefix+"_unmapped_T2C_1.fastq " +prefix+"_unmapped_T2C_2.fastq");
                file = file.replaceAll("\\{samtools}", "samtools view -b -F 4 {prefix}.bam > {prefix}_mapped.bam\n" +
                        "samtools view -b -f 256 {prefix}_mapped.bam > {prefix}_secondaryMapped.bam\n" +
                        "samtools view -b -F 256 {prefix}_mapped.bam > {prefix}_primaryMapped.bam\n" +
                        "samtools view -b -F 8 {prefix}_primaryMapped.bam > {prefix}_primaryMapped_mappedMate.bam\n" +
                        "samtools merge {prefix}_final.bam {prefix}_primaryMapped_mappedMate.bam {prefix}_secondaryMapped.bam");
            }
            file = file.replaceAll("\\{prefix}", prefix);

            Files.write(Paths.get(pathString), file.getBytes(charset));
            File finalFile = new File(pathString);
            finalFile.setExecutable(true);
        }catch (IOException e){
            e.printStackTrace();
        }

        return pathString;
    }
}
