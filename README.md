# ahcg_pipeline

Variant calling pipeline for genomic data analysis. This pipeline is part of the class Applied Human Computational Genomics - Fall 2017 of the Georgia Institute of Technology.


## Requirements

1. [Python3 - version 3.4.1](https://www.python.org/download/releases/3.4.1/)
2. [Trimmomatic - version 0.36](http://www.usadellab.org/cms/uploads/supplementary/Trimmomatic/Trimmomatic-0.36.zip)
3. [Bowtie2 - version 2.3.2](https://sourceforge.net/projects/bowtie-bio/files/bowtie2/2.3.2)
4. [Picard tools - version 2.11.0](https://github.com/broadinstitute/picard/releases/download/2.11.0/picard.jar)
5. [GATK - version 3.8](https://software.broadinstitute.org/gatk/download/)

## Reference genome

Reference genomes can be downloaded from [Illumina iGenomes](http://support.illumina.com/sequencing/sequencing_software/igenome.html)

## Test data

Use the following protocol to download and prepare test dataset from NIST sample NA12878

```{sh}
wget ftp://ftp-trace.ncbi.nih.gov/giab/ftp/data/NA12878/Garvan_NA12878_HG001_HiSeq_Exome/NIST7035_TAAGGCGA_L001_R1_001.fastq.gz
wget ftp://ftp-trace.ncbi.nih.gov/giab/ftp/data/NA12878/Garvan_NA12878_HG001_HiSeq_Exome/NIST7035_TAAGGCGA_L001_R2_001.fastq.gz
gunzip NIST7035_TAAGGCGA_L001_R1_001.fastq.gz
gunzip NIST7035_TAAGGCGA_L001_R2_001.fastq.gz
head -100000 NIST7035_TAAGGCGA_L001_R1_001.fastq > test_r1.fastq
head -100000 NIST7035_TAAGGCGA_L001_R2_001.fastq > test_r2.fastq
```

## Help

To access help use the following command:

```{sh}
python3 ahcg_pipeline.py -h
```
##  Execution example

### Building directory structure

```{sh}
mkdir -p data/reads data/reference data/adapters output 
```
### Downloading data
    
#### Read sample

We will download the sample SRR948996 from the SRA databse using the SRA-Toolkit:
    
```{sh}
fastq-dump --split-files SRR948994
```

#### Reference genome

```{sh}
wget ftp://igenome:G3nom3s4u@ussd-ftp.illumina.com/Homo_sapiens/NCBI/GRCh38/Homo_sapiens_NCBI_GRCh38.tar.gz
```
###  Running the pipeline

```{sh}
./ahcg_pipeline_v1.0.1.py \
-t /data2/AHCG2017FALL/bin/Trimmomatic-0.36/trimmomatic-0.36.jar  \
-b /data2/AHCG2017FALL/bin/bowtie2-2.2.9/bowtie2 \
-p /data2/AHCG2017FALL/bin/picard/picard.jar \
-g /data2/AHCG2017FALL/bin/GenomeAnalysisTK-3.8-0-ge9d806836 \
-i /data2/AHCG2017FALL/data \
-w /data2/AHCG2017FALL/reference_genome/Bowtie2Index/genome \
-d /data2/AHCG2017FALL/reference_genome/GATKResourceBundle/dbsnp_146.hg38.vcf.gz \
-r /data2/AHCG2017FALL/reference_genome/genome.fa \
-a /home/hfen3/adapters.fa \
-o output
```


