#!/bin/bash


samtools view -b -f 4 {prefix}.bam > {prefix}_unmapped.bam
{pe_unmapped}
gedi -e ExtractReads -strandness {strandness} {from} {to} -f {prefix}_unmapped.bam
rm {prefix}_unmapped.bam
samtools view -b -F 4 {prefix}.bam > {prefix}_final.bam
STAR --runMode alignReads --runThreadN 8   --genomeDir {pseudoStarIndex} --genomeLoad LoadAndKeep --limitBAMsortRAM 8000000000 --outFilterMismatchNmax 10 --outFilterScoreMinOverLread 0.4 --outFilterMatchNminOverLread 0.4 --alignIntronMax 1 --outSAMmode Full --readFilesIn *_unmapped_T2C*fastq --outSAMtype BAM SortedByCoordinate --alignEndsType Extend5pOfReads12 --outSAMattributes nM MD NH
samtools view -b -F 256 Aligned.sortedByCoord.out.bam > {prefix}_pseudoMapped.bam


gedi -e RescuePseudoReads -genome {genome} -pseudogenome {pseudogenome} -origmaps {prefix}_final.bam -pseudomaps {prefix}_pseudoMapped.bam  -idMap {prefix}_unmapped.idMap -strandness {strandness} {maxMM} {chrPrefix}
samtools merge {prefix}_rescued.bam {prefix}_final.bam {prefix}_pseudoMapped_reverted.bam
samtools index {prefix}_rescued.bam
gedi -e Bam2CIT {prefix}_rescued.cit {prefix}_rescued.bam


gedi -t . -e CorrectCIT {prefix}_rescued.cit
gedi -t . -e ReadCount -g {genome} -m Weight {prefix}_rescued.cit

{grandslam}

echo "Done"
