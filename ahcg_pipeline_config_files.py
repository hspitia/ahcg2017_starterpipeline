#!/usr/bin/env python3
import os
import sys
import glob
import logging
import argparse
import subprocess
import statistics as st
import ConfigParser
# import io


# ====================================================================
def calculateCoverageStats(bamFilePath, bedFilePath, outFilePath, samtoolsPath):
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
            
            covFile.write("{0}\t{1:.2f}\t{2:.2f}\t{3}\n".format(line.rstrip(), covMedian, covAverage, covMax))
                
            line = fp.readline()
            
# ====================================================================
def processConfigFile(config_path):
    # Load the configuration file
    with open(config_path) as f:
        sample_config = f.read()
        
    config = ConfigParser.RawConfigParser(allow_no_value=True)
    config.readfp(io.BytesIO(sample_config))

    # List all contents
    print("List all contents")
    for section in config.sections():
        print("Section: %s" % section)
        for options in config.options(section):
            print("x %s:::%s:::%s" % (options,
                                      config.get(section, options),
                                      str(type(options))))

    # Print some contents
    print("\nPrint some contents")
    print(config.get('other', 'use_anonymous'))  # Just get the value
    print(config.getboolean('other', 'use_anonymous'))  # You know the datatype?
    

# ====================================================================
def main(config_path):
    
    #Get complete path
    config_path      = os.path.abspath(config_path)
    
    processConfigFile(config_path)
    
    # ====================================================================

    print('Variant call pipeline completed')
    print('VCF file can be found at {0}'.format(vcf_path))
    return

if __name__ == '__main__':

    parser = argparse.ArgumentParser()
    parser.add_argument('-c', '--config', dest='config_path', type=str, help='Path to config file')
    
    # Enrichment kit (Nextera Enrichment Capture kit)
    # Download samples given an SRA accession
    # 
    
    if len(sys.argv) == 1:
        parser.print_help()
        print('----------------------------------------------------------------------------')
        print('{0} - ERROR: One or more required arguments are missing.\n'.format(parser.prog))
        sys.exit(1)
    
    args = parser.parse_args()
    

main(args.config_path)

