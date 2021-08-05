#!/bin/bash

wd=$(pwd)

mkdir -p {tmp}/{prefix}
cd {tmp}/{prefix}

cp {bampath} .
gedi -e ExtractUnmappedReads -all {tags} -f {prefix}.bam
{samtools}
STAR --runMode alignReads --runThreadN 8   --genomeDir {pseudoStarIndex} --genomeLoad LoadAndKeep --limitBAMsortRAM 8000000000 --outFilterMismatchNmax 20 --outFilterScoreMinOverLread 0.4 --outFilterMatchNminOverLread 0.4 --alignIntronMax 1 --outSAMmode Full --readFilesIn *_unmapped_T2C*fastq --outSAMtype BAM SortedByCoordinate --alignEndsType Extend5pOfReads12 --outSAMattributes nM MD NH
samtools view -b -F 256 Aligned.sortedByCoord.out.bam > {prefix}_pseudoMapped.bam


gedi -e RescuePseudoReads -genome {genome} -pseudogenome {pseudogenome} -origmaps {prefix}_final.bam -pseudomaps {prefix}_pseudoMapped.bam -prefix {prefix} -pairedEnd
samtools merge {prefix}_rescued.bam {prefix}_final.bam {prefix}_pseudoMapped_reverted.bam
samtools index {prefix}_rescued.bam
gedi -e Bam2CIT {prefix}_rescued.cit {prefix}_rescued.bam

#mv {prefix}_pseudoMapped_reverted_merged.cit {prefix}_rescued.cit
#mv {prefix}_pseudoMapped_reverted_merged.cit.metadata.json {prefix}_rescued.cit.metadata.json
gedi -t . -e CorrectCIT {prefix}_rescued.cit
gedi -t . -e ReadCount -g {genome} -m Weight {prefix}_rescued.cit
mv {prefix}_rescued.cit* $wd
echo "Done"
