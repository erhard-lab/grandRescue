package gedi.javapipeline;

import executables.Pipeline;
import gedi.util.FileUtils;
import gedi.util.io.text.LineIterator;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.net.URL;
import java.nio.charset.Charset;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.ArrayList;

public class createNo4sUBash {

    public static String createNo4sUBash(String origGenome, String origMapped, String tmpDir, String prefix) {
        String pathString = Paths.get(prefix + ".sh").toAbsolutePath().toString();
        try {
            Charset charset = StandardCharsets.UTF_8;
            String name = "no4sUBash.sh";

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
            file = file.replaceAll("\\{bampath}", origMapped);
            file = file.replaceAll("\\{genome}", origGenome);
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
