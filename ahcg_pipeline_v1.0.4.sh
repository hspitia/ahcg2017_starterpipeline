#!/bin/bash
script=$1
# python3 ahcg_pipeline_v1.0.4.py \
python3 $script \
-t /data2/AHCG2017FALL/bin/Trimmomatic-0.36/trimmomatic-0.36.jar \
-b /data2/AHCG2017FALL/bin/bowtie2-2.2.9/bowtie2 \
-p /data2/AHCG2017FALL/bin/picard/picard.jar \
-g /data2/AHCG2017FALL/bin/GenomeAnalysisTK-3.8-0-ge9d806836/GenomeAnalysisTK.jar \
-i /data2/AHCG2017FALL/data3/MenPa004DNA/Patient4_R1.fastq /data2/AHCG2017FALL/data3/MenPa004DNA/Patient4_R2.fastq \
-w /data2/AHCG2017FALL/reference_genome/Bowtie2Index/genome \
-r /data2/AHCG2017FALL/reference_genome/genome.fa \
-a /data2/AHCG2017FALL/bin/Trimmomatic-0.36/adapters/NexteraPE-PE.fa \
-o /data2/AHCG2017FALL/output3 \
-d /data2/AHCG2017FALL/reference_genome/GATKResourceBundle/dbsnp_146.hg38.vcf.gz \
-c /data2/AHCG2017FALL/bin/scripts/getCoverage.sh \
-s /data2/AHCG2017FALL/guardant360/guardant360.refGene_hg38.genes.bed \
-m /data2/AHCG2017FALL/bin/samtools-1.5/samtools
