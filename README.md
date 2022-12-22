## SLAM-Rescue
SLAM-Rescue is a software to circumvent mappability problems and correct for 4sU-induced changes of RNA metabolism in SLAM-seq data sets. To achieve this, SLAM-Rescue offers the tools to align previously unmappable reads in a T-to-C mismatch independent manner. Additionally, the tools provided by SLAM-Rescue can also be used to apply the mismatch-independent read alignment to similar methods that deal with other kinds of mismatches (e.g. C-to-T mismatches in Bisulfite-seq) 

# Prerequisites
- Java >= 1.8
- samtools >= 1.10
- STAR >= 2.7.3a (or any other read mapping tool you like)

# Download and installation

> To Do

# Preparation: Create a pseudo genome

To align the reads in a T-to-C (or any other) mismatch independent manner, you need to create a "pseudo genome". For this, you need the FASTA and GTF file of the reference genome you want to use and call the following function:

    gedi -e CreatePseudo -fasta reference.fa -gtf reference.gtf

Alternatively, you can specify the mismatches (e.g. A-to-G mismatches) for other methods and a separate output folder:

    gedi -e CreatePseudo -fasta reference.fa -gtf reference.gtf -from A -to G -out pseudoGenomeFolder

This will generate the FASTA and GTF files, which you can then use for the read mapper of your choice (e.g. to create a STAR-index)

# Preparation: Prepare a gedi genome

SLAM-Resuce uses the [GEDI framework](https://github.com/erhard-lab/gedi) and therefore, the reference genome as well as the pseudo genome need to be indexed with gedi as well.

    gedi -e IndexGenome -s reference.fa -a reference.gtf -n h.ens90

    gedi -e IndexGenome -s reference_pseudo_T2C.fa -a reference_pseudo.gtf -n h.ens90_pseudo

With -n you can specify a name to access the indexed gedi genome in other applications like SLAM-Rescue. For more on gedi genomes, see [here](https://github.com/erhard-lab/gedi/wiki/Genomic).


# SLAM-Rescue functions

The complete mapping and rescue procedure consists of the following steps:

1. Map your sequenced reads to a reference genome (keeping the unmapped reads!)
2. Extract the unmapped reads
3. Map the previously unmappable reads to a pseudogenome
4. Revert the pseudo-mapped reads to their original sequence and position on the reference genome
5. Combine your original mapping run and the rescued reads to a final bam file

We will go over all these steps one by one, but under the section "Automated SLAM-Rescue" you can find a function to create a ready-to-use script that will run all the necessary steps for you.

**1. Map your sequenced reads to a reference genome**

You can map your reads with any tool you like. For the following steps, you only need a BAM file which contains the mapped and unmapped reads.
For STAR, you can keep the unmapped reads by using the parameter *--outSAMunmapped Within*

**2. Extract the unmapped reads** 

We want to map the previously unmapped reads in a T-to-C mismatch independent manner, so we first need to extract those reads from the given BAM file and write them to a new FASTQ file (or two, if we are dealing with paired end data):

    gedi -e ExtractReads -f reads.bam

The ExtractReads function has the following [optional] parameters:

    gedi -e ExtractReads
    
    -f: The BAM file of your mapping run. (Needs to contain unmapped reads!)
    [-strandness: Set the strandness of your data set: Sense or Antisense. Default: Sense]
    [-from: Set which nucleotide to convert if you're not dealing with T-to-C mismatches.]
    [-to: Set to which nucleotide to convert if you're not dealing with T-to-C mismatches.]
    
As an output, you get two (or three for paired end) files: reads_unmapped_T2C.fastq and reads.idMap.
The FASTQ file will be used for mapping, the idMap is necessary for a later function.

**3. Map the previously unmappable reads to a pseudogenome**

Now, you want to map the FASTQ file to the pseudo genome created earlier. You can again use the same mapping tool and we suggest reducing the maximum amount of allowed mismatches compared to your original run, since one type of mismatches (e.g. T-to-C) will be ignored here. Let's call this file pseudoGenomeMapped.bam. It will be used again in the next step.

For example, the parameters for mapping with STAR could look like this:
    
    STAR --runMode alignReads --genomeDir /path/to/pseudogenome/STAR-index/ --genomeLoad LoadAndKeep --outFilterMismatchNmax 10 --outFilterScoreMinOverLread 0.4 --outFilterMatchNminOverLread 0.4  --readFilesIn reads_unmapped_T2C.fastq --outSAMmode Full --outSAMtype BAM SortedByCoordinate --outSAMattributes nM MD NH
    mv Aligned.sortedByCoord.out.bam pseudoGenomeMapped.bam

**4. Revert the pseudo-mapped reads to their original sequence and position on the reference genome**

The reads mapped to the pseudo genome now need to be returned to their original position in the reference genome and their original sequence needs to be restored. For this, you will need the idMap file from step 2.

    gedi -e RescuePseudoReads -genome h.ens90 -pseudogenome h.ens90_pseudo -origmaps reads.bam -pseudomaps pseudoGenomeMapped.bam -idMap reads.idMap

The RescuePseudoReads function has the following [optional] parameters:

    gedi -e RescuePseudoReads
    
    -genome: The gedi name of your reference genome (e.g. h.ens90)
    -pseudogenome: The gedi name of your pseudo geonme (e.g. h.ens90_pseudo)
    -origmaps: The BAM file of your mapping run. (Needs to contain unmapped reads!)
    [-strandness: Set the strandness of your data set: Sense or Antisense. Default: Sense]
    [-maxMM: Reads with more mismatches than the specified int will be removed. Default: 75% of read length]
    [-chrPrefix: Prefix/Pattern for chromosome names, * replaces chromosome number. E.g. chr*]
    
As an output, you will get one bam file: reads_reverted.bam.

**5. Combine your original mapping run and the rescued reads to a final bam file**

Finally, you only need to combine your original BAM file with the rescued reads to a final BAM file. You can easily do this with samtools.
First, remove the unmapped reads from your original BAM file and then merge it with the rescued reads:
    
    samtools view -b -F 4 reads.bam > reads_mapped.bam
    samtools merge reads_final.bam reads_mapped.bam reads_reverted.bam

# Automated SLAM-Rescue

The easiest way to use SLAM-Rescue for your data set is the function gedi -e ReadRescue. It creates a ready-to-use bash-file with all intermediate steps to process your mapped BAM file. This function assumes you use STAR. If you use another mapping tool, just edit the bash-file after its creation and replace the STAR call with the mapping call of your choice.


    gedi -e ReadRescue -genome h.ens90 -pseudogenome h.ens90_pseudo -pseudoSTAR /folder/to/pseudoGenome/STAR-index -f reads.bam


The created file will be named like your inserted bam-file. You can then open and edit it to your liking and run the bash-file as follows:

    ./reads.sh


The ReadRescue function has the following [optional] parameters:

    gedi -e ReadRescue
    
    -genome: The gedi name of your reference genome (e.g. h.ens90)
    -pseudogenome: The gedi name of your pseudo geonme (e.g. h.ens90_pseudo)
    -f: The BAM file of your mapping run. (Needs to contain unmapped reads!)
    [-pseudoSTAR: The folder of your STAR-indexed genome. If you dont use STAR, you don't need this parameter, but have to edit the bash-file accordingly.]
    [-pe: Set this flag, if you have paired-end data]
    [-strandness: Set the strandness of your data set: Sense or Antisense. Default: Sense]
    [-from: Set which nucleotide to convert if you're not dealing with T-to-C mismatches.]
    [-to: Set to which nucleotide to convert if you're not dealing with T-to-C mismatches.]
    [-maxMM: Reads with more mismatches than the specified int will be removed. Default: 75% of read length]
    [-chrPrefix: Prefix/Pattern for chromosome names, * replaces chromosome number. E.g. chr*]


