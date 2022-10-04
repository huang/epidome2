#!/bin/bash
#SBATCH -D /srv/data/MPV/ -c 2 --mem=10G -J cutadapt_g216 -p daytime

cd ${1}

for i in ${2}; do SAMPLE=$(echo ${i} | sed "s/_R1_\001\.fastq\.gz//"); echo ${SAMPLE}_R1_001.fastq.gz ${SAMPLE}_R2_001.fastq.gz; cutadapt -e 0.06 -g TGGGTATGRCAATCACTTTACA -G GCATCAAAAGCACTCTCATTACC --pair-filter=any -o ${3}/${SAMPLE}_R1_001.fastq.gz --paired-output ${3}/${SAMPLE}_R2_001.fastq.gz --discard-untrimmed ${SAMPLE}_R1_001.fastq.gz ${SAMPLE}_R2_001.fastq.gz; done