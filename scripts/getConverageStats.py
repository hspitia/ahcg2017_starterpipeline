#!/usr/bin/env python3

import os
import sys
import glob
import logging
import argparse
import shutil
import subprocess
import statistics as st
import configparser as cp
import textwrap

# ==============================================================================
def calculateCoverageStats(bamFilePath, bedFilePath, outFilePath, samtoolsPath):
    """Calculate mean, median and average coverage of each region
    in bedFilePath by calling samtools depth, and save results to outFilePath."""
    covFile  = open(outFilePath, 'w')
    covFile.write("#chr\tstart\tstop\tname\tscore\tstrand\tmedian_cov\tavg_cov\tmax_cov\n")
    with open(bedFilePath) as fp:  
        line = fp.readline()
        while line:
            tokens     = line.split("\t")
            region     = "{}:{}-{}".format(tokens[0],tokens[1],tokens[2])
            samcmd     = [samtoolsPath, "depth", "-r", region, bamFilePath]
            covLines   = subprocess.check_output(samcmd, universal_newlines=True).splitlines()
            perBaseCov = []
            for l in covLines:
                tokens = l.split("\t")
                perBaseCov.append(int(tokens[2]))
            
            covMedian  = 0
            covAverage = 0
            covMax     = 0

            if len(perBaseCov) > 0:
                covMedian  = st.median(perBaseCov)
                covAverage = st.mean(perBaseCov)
                covMax     = max(perBaseCov)
            
            covFile.write("{0}\t{1:.2f}\t{2:.2f}\t{3}\n".format(line.rstrip(), 
                covMedian, covAverage, covMax))
                
            line = fp.readline()
    
    covFile.close()
# ==============================================================================
def getlist(option, sep=',', chars=None):
    """Return a list from a ConfigParser option. By default, 
       split on a comma and strip whitespaces."""
    return [ chunk.strip(chars) for chunk in option.split(sep) ]
# ==============================================================================
def checkTool(key, name, confOptions):
    """Check if the executable file key defined in confOptions or present the
    system PATH, exists"""
    if key in confOptions:
        confOptions[key] = os.path.abspath(confOptions[key])
        if not os.path.exists(confOptions[key]):
            sys.stderr.write('ERROR: {} not found at {}\n'.format(name, confOptions[key]))
            sys.exit(1)
    else:
        toolPath = shutil.which(key)
        if toolPath != None:
                confOptions[key] = toolPath
        else:
            sys.stderr.write('ERROR: {} not found in system and not provided in the config file.\n'.format(name))
            sys.exit(1)
# ==============================================================================
def checkDataFile(key, name, confOptions):
    """Check if the data file defined by key in confOptions, exists"""
    if key in confOptions:
        confOptions[key] = os.path.abspath(confOptions[key])
        if key == 'index':
            indicies = glob.glob('{0}.*.bt2'.format(confOptions[key]))
            if len(indicies) == 0:
                sys.stderr.write('ERROR: {} not found at {}\n'.format(name, confOptions[key]))
                sys.exit(1)
        else:
            if not os.path.exists(confOptions[key]):
                sys.stderr.write('ERROR: {} not found at {}\n'.format(name, confOptions[key]))
                sys.exit(1)
    else:
        sys.stderr.write('ERROR: {} ({}) not found not provided in the config file\n'.format(name, key))
        sys.exit(1)
# ==============================================================================
def main(bamFilePath, bedFilePath, outputFilename, samtoolsPath):
    # Process config file
    bamFilePath = os.path.abspath(bamFilePath)
    bedFilePath = os.path.abspath(bedFilePath)
    # outputFilename = os.path.abspath(outputFilename)
    # 
    print("Samtools: {}".format(samtoolsPath))
    
    if samtoolsPath == None:
        toolPath = shutil.which(samtoolsPath)
        if toolPath != None:
                samtoolsPath = toolPath
        else:
            sys.stderr.write('ERROR: {} not found in system. You can provide the Samtools binary location using the option -s (--samtools)\n'.format("Samtools"))
            sys.exit(1)
    else:
        if not os.path.exists(samtoolsPath):
            sys.stderr.write('ERROR: Samtools not found at {}\n'.format(samtoolsPath))
            sys.exit(1)
    
    # ==========================================================================
    # Check for files
    if not os.path.exists(bamFilePath):
        sys.stderr.write('ERROR: BAM file not found at {0}'.format(bamFilePath))
        sys.exit(1)
    
    if not os.path.exists(bedFilePath):
        sys.stderr.write('ERROR: BED file not found at {0}'.format(bedFilePath))
        sys.exit(1)
    
    # ==========================================================================
    # Coverage stats calculation from the final BAM file
    # covfile_path = '{0}/{1}_{2}_coverage.tsv'.format(
    #     confData['outputdir'], 
    #     os.path.splitext(os.path.basename(confData['geneset']))[0], 
    #     os.path.splitext(os.path.basename(fbam_path))[0])
    
    covfile_path = outputFilename
    
    print("INFO: Coverage stats calculation - Started")
    calculateCoverageStats(bamFilePath, bedFilePath, covfile_path, 
                           "samtools")
    print("INFO: Coverage stats calculation - Done!\n")
    return 0;
    # ====================================================================

if __name__ == '__main__':

    parser = argparse.ArgumentParser(
        formatter_class=argparse.RawTextHelpFormatter,
        epilog = textwrap.dedent(''' '''))
    parser.add_argument('bamFilePath', metavar = 'BAM_FILE', type = str, 
        help = 'Path to BAM file')
    
    parser.add_argument('bedFilePath', metavar = 'BED_FILE', type = str, 
        help = 'Path to BED file')
    
    parser.add_argument('-s', '--samtools', dest='samtools_path', type=str, 
        help='Path to Samtools')
    
    parser.add_argument('outputFilename', metavar = 'OUTPUT_FILENAME', type = str, 
        help = 'Output filename')
    
    parser.add_argument('-v', '--version', action = 'version', 
        version = '''%(prog)s 0.1
Copyright (c) 2017 Georgia Institute of Technology
Applied Human Computational Genomics - Fall 2017''')
    
    # if len(sys.argv) == 1:
    #     parser.print_help()
    #     print('--------------------------------------------------------------------------------')
    #     print('{0} - ERROR: One or more required arguments are missing.\n'.format(parser.prog))
    #     sys.exit(1)
    
    args = parser.parse_args()
    
    sys.exit(main(args.bamFilePath, args.bedFilePath, args.outputFilename, 
        args.samtools_path))

