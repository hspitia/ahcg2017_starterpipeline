#!/bin/bash

localExec=$1

if [[ "$localExec" == "" ]]; then
    python3 ahcg_pipeline.test.py \
    -t /data2/AHCG2017FALL/bin/Trimmomatic-0.36/trimmomatic-0.36.jar \
    -g /data2/AHCG2017FALL/bin/GenomeAnalysisTK-3.8-0-ge9d806836/GenomeAnalysisTK.jar \
    -b /data2/AHCG2017FALL/bin/bowtie2-2.2.9/bowtie2 \
    -p /data2/AHCG2017FALL/bin/picard/picard.jar \
    -i reads_1.fastq reads_2.fastq \
    -w /data2/AHCG2017FALL/reference_genome/Bowtie2Index/genome \
    -r /data2/AHCG2017FALL/reference_genome/genome.fa \
    -a /data2/AHCG2017FALL/bin/Trimmomatic-0.36/adapters/NexteraPE-PE.fa \
    -o ./out \
    -d /data2/AHCG2017FALL/reference_genome/GATKResourceBundle/dbsnp_146.hg38.vcf.gz \
    -c /data2/AHCG2017FALL/bin/scripts/getCoverage.sh \
    -s /data2/AHCG2017FALL/guardant360/guardant360.refGene_hg38.genes.bed \
    -m /data2/AHCG2017FALL/bin/samtools-1.5/samtools    
else
    python3 ahcg_pipeline.test.py \
    -t lib/Trimmomatic-0.36/trimmomatic-0.36.jar \
    -g lib/GenomeAnalysisTK.jar \
    -b lib/bowtie2-2.2.9/bowtie2 \
    -p lib/picard.jar \
    -i reads_1.fastq reads_2.fastq \
    -w /data2/AHCG2017FALL/reference_genome/Bowtie2Index/genome \
    -r /data2/AHCG2017FALL/reference_genome/genome.fa \
    -a lib/Trimmomatic-0.36/adapters/NexteraPE-PE.fa \
    -o ./out \
    -d /data2/AHCG2017FALL/reference_genome/GATKResourceBundle/dbsnp_146.hg38.vcf.gz \
    -c lib/getCoverage.sh \
    -s data/guardant360.refGene_hg38.genes.bed \
    -m lib/samtools-1.6/samtools
fi


