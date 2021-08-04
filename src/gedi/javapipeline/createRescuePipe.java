package gedi.javapipeline;

import gedi.util.program.GediProgram;
import gedi.util.program.GediProgramContext;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;

import static gedi.javapipeline.createRescueBash.createRescueBash;

public class createRescuePipe extends GediProgram {

    public createRescuePipe(RescueParameterSet params) {
        addInput(params.nthreads);
        addInput(params.genome);
        addInput(params.pseudogenome);
        addInput(params.prefix);
        addInput(params.all);
        addInput(params.pseudoSTAR);
        addInput(params.tmp);
        addInput(params.files);
        addInput(params.slam);
        addInput(params.pairedEnd);

        addOutput(params.bashFile);
    }

    @Override
    public String execute(GediProgramContext context) throws IOException {
        int threads = getIntParameter(0);
        String origGenome = getParameter(1);
        String pseudoGenome = getParameter(2);
        String prefix = getParameter(3);
        boolean writeAll = getBooleanParameter(4);
        String pseudoStarIndex = getParameter(5);
        String tmpDir = getParameter(6);
        ArrayList<String> files = getParameters(7);
        boolean slam = getBooleanParameter(8);
        boolean pairedEnd = getBooleanParameter(9);

        int threadCount = 1;
        int round = 0;
        String mergeCommand = "gedi -t . -e MergeCIT -c "+prefix+"_rescued.cit";
        String wait = "";

        File bash = new File(prefix + "_readRescue.sh");
        BufferedWriter writer = new BufferedWriter(new FileWriter(bash));

        writer.append("#!/bin/bash\n\n");
        writer.append("STAR --genomeLoad LoadAndExit --genomeDir "+pseudoStarIndex+"\n\n");


        for(String file : files){
            String pathString = createRescueBash(writeAll, origGenome, pseudoGenome, file, pseudoStarIndex, tmpDir, getPrefix(file), pairedEnd);
            System.out.println("--" + getPrefix(file));

            writer.append("echo $( date +\"%F %T\" ) Starting " + getPrefix(file)+"\n");
            mergeCommand = mergeCommand + " "+pathString.replace(".sh", "_rescued.cit");
            if(threadCount%threads != 0){
                writer.append(pathString+" & \n");
                writer.append("PIDS["+(threadCount+(round*threads))+"]=$!\n\n");
                threadCount++;
            } else {
                writer.append(pathString+" & \n");
                writer.append("PIDS["+(threadCount+(round*threads))+"]=$!\n\n");
                wait = "wait";
                threadCount = 1;

                for(int i = threadCount; i <= threads; i++){
                    wait = wait + " ${PIDS["+(i+(round*threads))+"]}";
                }
                writer.append(wait+"\n\n");
                round++;
            }
        }
        if(threadCount%threads != 1 && threads != 1){
            wait = "wait";

            for(int i = 1; i < threadCount; i++){
                wait = wait + " ${PIDS["+(i+(round*threads))+"]}";
            }
            writer.append(wait+"\n");
        }

        writer.append("\n\n");
        writer.append("STAR --genomeLoad Remove --genomeDir "+pseudoStarIndex+"\n\n");

        writer.append(mergeCommand+"\n");
        if(slam) {
            writer.append("gedi -t . -e Slam -nthreads " + threads + " -trim5p 15 -genomic " + origGenome + " -prefix grandslam_t15/" + prefix + " -reads " + prefix+"_rescued.cit -plot -D\n\n");
        }
        writer.close();
        bash.setExecutable(true);

        return null;

    }

    public static String getPrefix(String path){
        return path.substring(path.lastIndexOf("/")+1, path.lastIndexOf(".bam"));
    }
}
