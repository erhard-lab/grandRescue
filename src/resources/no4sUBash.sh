#!/bin/bash

wd=$(pwd)

mkdir -p {tmp}/{prefix}
cd {tmp}/{prefix}

cp {bampath} .

samtools sort -o {prefix}_sorted.bam {prefix}.bam
samtools index {prefix}_sorted.bam
gedi -e Bam2CIT {prefix}.cit {prefix}_sorted.bam


gedi -t . -e CorrectCIT {prefix}.cit
gedi -t . -e ReadCount -g {genome} -m Weight {prefix}.cit
mv {prefix}.cit* $wd

echo "Done"
