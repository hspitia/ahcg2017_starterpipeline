# ahcg_pipeline

Variant calling pipeline for genomic data analysis. This pipeline is part of the class Applied Human Computational Genomics - Fall 2017 of the Georgia Institute of Technology.


## Requirements

1. [Python3 - v3.4.1](https://www.python.org/download/releases/3.4.1/)
2. [Trimmomatic - v0.36](http://www.usadellab.org/cms/uploads/supplementary/Trimmomatic/Trimmomatic-0.36.zip)
3. [Bowtie2 - v2.3.2](https://sourceforge.net/projects/bowtie-bio/files/bowtie2/2.3.2)
4. [Picard tools - v2.11.0](https://github.com/broadinstitute/picard/releases/download/2.11.0/picard.jar)
5. [GATK - v3.8](https://software.broadinstitute.org/gatk/download/)
6. [Samtools - v1.6](https://downloads.sourceforge.net/project/samtools/samtools/1.6/samtools-1.6.tar.bz2?r=https%3A%2F%2Fsourceforge.net%2Fprojects%2Fsamtools%2F&ts=1510018121&use_mirror=phoenixnap)
7. [Control-FREEC v11.0](https://github.com/BoevaLab/FREEC/archive/v11.0.tar.gz)
8. [R language - v3.3.2](https://cran.cnr.berkeley.edu/)

## Help

To access help use the following command:

```{sh}
./ahcg_pipeline.py -h
```

## Configuration file

A configuration file is required to run the pipeline. This config file has two 
required sections `[tools]` and `[data]` in which the required software tools and input data are described. An additional optional section `[freec-control]` can be included to specify options for the `[control]` section of the configuration of the Control-FREEC tool (see details [here](http://boevalab.com/FREEC/tutorial.html))

All the options are 



### `[data]` section

| Option       | Description                                                                  |
|--------------|------------------------------------------------------------------------------|
| `inputfiles` | List of paired end read files (comma sparated)                               |
| `geneset`    | Path to the bed file with genes of interest to calculate coverage statistics |
| `outputdir`  | Path to the output directory                                                 |
| `adapters`   | Path to adapters fasta file to perform sequence trimming with Trimmomatic    |
| `chrlenfile` | Path to file with chromosome lengths for Control-FREEC                       |
| `chrfiles`   | Path to the directory with chromosomes fasta files for Control-FREEC         |
| `dbsnp`      | Path to dbSNP vcf file for GATK                                              |
| `index`      | Path to the prefix of the reference Bowtie2 index                            |
| `reference`  | Path to the reference genome fasta file                                      |

### `[tools]` section

| Option        | Description                                                                                                |
|---------------|------------------------------------------------------------------------------------------------------------|
| `bowtie2`     | Path to Bowtie2 executable                                                                                 |
| `freec`       | Path to Control-FREEC executable                                                                           |
| `gatk`        | Path to GATK jar file                                                                                      |
| `makegraph`   | Path to Control-FREEC `makeGraph.R` script (Usually in the folder `scripts` at the Control-FREEC root dir) |
| `picard`      | Path to Picard jar file                                                                                    |
| `samtools`    | Path to Samtools executable                                                                                |
| `trimmomatic` | Path to Trimmomatic jar file                                                                               |

### `[freec-control]` section

Optional section for Control-FREEC's config file `control` section parameters

| Option          | Description                                                                                                                                                           |
| --------------- | --------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| mateFile        | Path to file to act as control of the current sample. See control-FREEC manual for details>                                                                           |
| inputFormat     | Format of mateFile (SAM, BAM, pileup and others. See control-FREEC manual for details)>                                                                               |
| mateOrientation | Orientation of reads in mateFile. 0 - single ends), RF - Illumina mate-pairs, FR - Illumina paired-ends), FF - SOLiD mate-pairs. See control-FREEC manual for details |

### Example of config file

```
[data]
# inputfiles      = /data2/AHCG2017FALL/data4/SRR2530742_1.fastq,/data2/AHCG2017FALL/data4/SRR2530742_2.fastq
sraid           = SRR2530742
geneset         = /data2/AHCG2017FALL/guardant360/guardant360.refGene_hg38.genes.bed
outputdir       = /data2/AHCG2017FALL/output5

adapters        = /data2/AHCG2017FALL/bin/Trimmomatic-0.36/adapters/NexteraPE-PE.fa
chrlenfile      = /data2/AHCG2017FALL/reference_genome/chromosomeSizes.txt
chrfiles        = /data2/AHCG2017FALL/reference_genome/chroms/
dbsnp           = /data2/AHCG2017FALL/reference_genome/GATKResourceBundle/dbsnp_146.hg38.vcf.gz
index           = /data2/AHCG2017FALL/reference_genome/Bowtie2Index/genome
reference       = /data2/AHCG2017FALL/reference_genome/genome.fa

[tools]
assesssig       = /data2/AHCG2017FALL/bin/FREEC/scripts/assess_significance.R
bowtie2         = /data2/AHCG2017FALL/bin/bowtie2-2.2.9/bowtie2
fastq-dump      = /data2/AHCG2017FALL/bin/sratoolkit/bin/fastq-dump
freec           = /data2/AHCG2017FALL/bin/FREEC/src/freec
gatk            = /data2/AHCG2017FALL/bin/GenomeAnalysisTK-3.8-0-ge9d806836/GenomeAnalysisTK.jar
java            = /data2/AHCG2017FALL/bin/java-1.8/bin/java
makegraph       = /data2/AHCG2017FALL/bin/FREEC/scripts/makeGraph.R
picard          = /data2/AHCG2017FALL/bin/picard/picard.jar
samtools        = /data2/AHCG2017FALL/bin/samtools-1.5/samtools
trimmomatic     = /data2/AHCG2017FALL/bin/Trimmomatic-0.36/trimmomatic-0.36.jar

[freec-control]
mateFile        = /data2/AHCG2017FALL/output4/SRR2530741_1_trimmed_final.bam
inputFormat     = BAM
mateOrientation = FR
```

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

##  Execution example

### Building directory structure

```{sh}
mkdir -p data/reads data/reference data/adapters output
```
### Downloading data

#### Read sample

We can download the sample SRR948996 from the SRA databse using the SRA-Toolkit
to specify their path using `inputfiles` option in the config file:

```{sh}
fastq-dump --split-files SRR948994
```

Alternatively, you can specify the id of the sample (SRR948994) using the
`sraid` option in the config file to delegate the download of the sample files
to the pipeline.

#### Reference genome

```{sh}
wget ftp://igenome:G3nom3s4u@ussd-ftp.illumina.com/Homo_sapiens/NCBI/GRCh38/Homo_sapiens_NCBI_GRCh38.tar.gz
```
###  Running the pipeline

```{sh}
./ahcg_pipeline_v1.0.1.py -c config_file.txt
```
