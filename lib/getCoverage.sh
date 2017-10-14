#!/bin/bash

########################################################################
# Author:       Hector Espitia
# Intitution:   Georgia Institute of Technology
# Version:      0.4
# Date:         10/14/2017

# Description:  
#       Given a BAM file, calculates the mean and max coverage of each 
#       region defined in a BED file.
# Requirements:
#   samtools >= 1.5
#   awk >= 3.1.7
#
########################################################################

bamFile=$1;
bedFile=$2;
outFile=$3;

IFS= read -d '' usage << "EOF"
gePerBaseCoverage.sh\\n
\\b====================\\n
\\n
\\bGiven a BAM file, this script calculates the mean and max coverage of the regions defined in a BED file.
\\n\\n
\\busage: \\n\\n

  ./gePerBaseCoverage.sh <BAM_FILE> <BED_FILE> <OUT_FILE>\\n

EOF

# Check for help argument
if [[ $bamFile == "" || $bamFile == "-h" || $bamFile == "--help" ]]; then
    echo -e $usage;
    exit 1;
fi

# Check for required files
if [[ ! -e "$bamFile" ]]; then
    echo -e "ERROR: The file '$bamFile' does not exist."
    exit 1;
fi

if [[ ! -e "$bedFile" ]]; then
    echo -e "ERROR: The file '$bedFile' does not exist."
    exit 1;
fi

# Read the bed file
while IFS='' read -r line || [[ -n "$line" ]]; do
    chr=$(echo -e $line | cut -d " " -f1); 
    start=$(echo -e $line | cut -d " " -f2);
    end=$(echo -e $line | cut -d " " -f3);
    # name=$(echo -e $line | cut -d " " -f4);
    # score=$(echo -e $line | cut -d " " -f5);
    # strand=$(echo -e $line | cut -d " " -f6);
    
    region="${chr}:${start}-${end}"
    
    # extract the region's per-base coverage info and compute average and max values
    echo -n ${line} | sed 's/ /\t/g' | tee -a $outFile
    samtools depth -r ${region} ${bamFile} | awk 'BEGIN{OFS="\t"; max = 0} {sum += $3; if($3 > max) max = $3 } END{mean=0; if(NR>0) mean = sum/NR; print "\t"mean, max}' | tee -a $outFile;
    
done < "$bedFile"

