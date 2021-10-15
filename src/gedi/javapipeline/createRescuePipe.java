package gedi.javapipeline;

import gedi.core.reference.Strandness;
import gedi.util.program.GediProgram;
import gedi.util.program.GediProgramContext;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;

import static gedi.javapipeline.createRescueBash.createRescueBash;
import static gedi.javapipeline.createNo4sUBash.createNo4sUBash;
import static gedi.RescueReads.getPrefix;

public class createRescuePipe extends GediProgram {

    public createRescuePipe(RescuePipeParameterSet params) {
        addInput(params.nthreads);
        addInput(params.genome);
        addInput(params.pseudogenome);
        addInput(params.prefix);
        addInput(params.pseudoSTAR);
        addInput(params.tmp);
        addInput(params.files);
        addInput(params.slam);
        addInput(params.k);
        addInput(params.tags);
        addInput(params.strandness);

        addOutput(params.bashFile);
    }

    @Override
    public String execute(GediProgramContext context) throws IOException {
        int threads = getIntParameter(0);
        String origGenome = getParameter(1);
        String pseudoGenome = getParameter(2);
        String prefix = getParameter(3);
        String pseudoStarIndex = getParameter(4);
        String tmpDir = getParameter(5);
        ArrayList<String> files = getParameters(6);
        boolean slam = getBooleanParameter(7);
        boolean k = getBooleanParameter(8);
        ArrayList<String> tags = getParameters(9);
        Strandness strandness = getParameter(10);

        int threadCount = 1;
        int round = 0;
        String mergeCommand = "gedi -t . -e MergeCIT -c "+prefix+"_rescued.cit";
        String wait = "";

        File bash = new File(prefix + "_readRescue.sh");
        BufferedWriter writer = new BufferedWriter(new FileWriter(bash));

        writer.append("#!/bin/bash\n\n");
        writer.append("STAR --genomeLoad LoadAndExit --genomeDir "+pseudoStarIndex+"\n\n");


        for(String file : files){
            String pathString = "";
            if(!file.contains("no4sU")) {
                pathString = createRescueBash(origGenome, pseudoGenome, file, pseudoStarIndex, tmpDir, getPrefix(file), tags, strandness);
                mergeCommand = mergeCommand + " " + pathString.replace(".sh", "_rescued.cit");
            } else {
                pathString = createNo4sUBash(origGenome, file, tmpDir, getPrefix(file));
                mergeCommand = mergeCommand + " " + pathString.replace(".sh", ".cit");
            }
            System.out.println("--" + getPrefix(file));
            writer.append("echo $( date +\"%F %T\" ) Starting " + getPrefix(file) + "\n");
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
        if(!k) {
            writer.append("rm -r " + tmpDir + "\n\n");
        }
        writer.append("STAR --genomeLoad Remove --genomeDir "+pseudoStarIndex+"\n\n");

        writer.append(mergeCommand+"\n");
        if(slam) {
            writer.append("gedi -t . -e Slam -nthreads 10" + " -trim5p 15 -genomic " + origGenome + " -prefix grandslam_t15/" + prefix + " -reads " + prefix+"_rescued.cit -modelall -plot -D\n\n");
        }
        writer.close();
        bash.setExecutable(true);

        return null;

    }


}
