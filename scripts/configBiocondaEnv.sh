#!/bin/bash

# Go to home directory
cd ~
# Get the Anaconda Python distribution installer
wget https://repo.continuum.io/archive/Anaconda3-5.0.1-Linux-x86_64.sh
# Execute the installer and follow the instructions
bash Anaconda3-4.4.0-Linux-x86_64.sh
source ~/.bashrc
# Add conda channels
conda config --add channels r
conda config --add channels defaults
conda config --add channels conda-forge
conda config --add channels bioconda
# Create an environment
conda create -y --name ahcg_pipeline
# Activate the environment
source activate ahcg_pipeline
# Install required applications and packages
conda install -y samtools trimmomatic bowtie2 picard control-freec openjdk r-gplots r-ggplot2 r-gsalib r-reshape gatk
# Download and register GATK into bioconda environment
# Get GATK latest version package
wget "https://software.broadinstitute.org/gatk/download/auth?package=GATK" -O GenomeAnalysisTK-3.8-0.tar.bz2
# Uncompress GATK package
tar xvf GenomeAnalysisTK-3.8-0.tar.bz2
# Register GATK into bioconda environment
gatk-register ~/GenomeAnalysisTK-3.8-0-ge9d806836/GenomeAnalysisTK.jar
# Remove GATK installer files
rm -fr ~/GenomeAnalysisTK-3.8-0-ge9d806836 GenomeAnalysisTK-3.8-0.tar.bz2

