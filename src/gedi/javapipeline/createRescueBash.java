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
import java.util.ArrayList;

public class createRescueBash {


    public static String createRescueBash(boolean writeAll, String origGenome, String pseudoGenome, String origMapped, String pseudoStarIndex, String tmpDir, String prefix, ArrayList<String> tags) {
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

            if(!writeAll){
                file = file.replace("-all", "");
            }
            file = file.replaceAll("\\{tmp}", tmpDir);
            file = file.replaceAll("\\{bampath}", origMapped);
            file = file.replaceAll("\\{pseudoStarIndex}", pseudoStarIndex);
            file = file.replaceAll("\\{genome}", origGenome);
            file = file.replaceAll("\\{pseudogenome}", pseudoGenome);
            file = file.replaceAll("\\{files}", prefix+"_unmapped_T2C.fastq");
            file = file.replaceAll("\\{prefix}", prefix);
            file = file.replaceAll("\\{tags}", tagnames);



            Files.write(Paths.get(pathString), file.getBytes(charset));
            File finalFile = new File(pathString);
            finalFile.setExecutable(true);
        }catch (IOException e){
            e.printStackTrace();
        }

        return pathString;
    }
}
