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
def checkTool(key, name, confOptions):
    """Check if the executable file key defined in confOptions or present the
    system PATH, exists"""
    if key in confOptions:
        confOptions[key] = os.path.abspath(confOptions[key])
        if not os.path.exists(confOptions[key]):
            raise FileNotFoundError('{} not found at {}'.format(name, confOptions[key]))
            # sys.stderr.write('ERROR: {} not found at {}\n'.format(name, confOptions[key]))
            # sys.exit(1)
    else:
        toolPath = shutil.which(key)
        if toolPath != None:
                confOptions[key] = toolPath
        else:
            raise FileNotFoundError('{} not found in system and not provided in the config file'.format(name))
            # sys.stderr.write('ERROR: {} not found in system and not provided in the config file.\n'.format(name))
            # sys.exit(1)
# ==============================================================================
def getlist(option, sep=',', chars=None):
    """Return a list from a ConfigParser option. By default, 
       split on a comma and strip whitespaces."""
    return [ chunk.strip(chars) for chunk in option.split(sep) ]
# ==============================================================================
def checkDataOption(key, name, confOptions):
    """Check if the data file key defined in confOptions, exists"""
    if key in confOptions:
        if key == 'inputfiles':
            files = getlist(confOptions[key])
            files = [os.path.abspath(fs) for fs in files]
            for f in files:
                if not os.path.exists(f):
                    raise FileNotFoundError('Fastq file not found at {}'.format(f))
                    # sys.stderr.write('ERROR: Fastq file not found at {}\n'.format(f))
                    # sys.exit(1)
        elif key == 'index':
            indicies = glob.glob('{0}.*.bt2'.format(confOptions[key]))
            if len(indicies) == 0:
                raise FileNotFoundError('Bowtie index not found at {0}'.format(confOptions[key]))
        elif key == 'outputdir':
            if not os.path.exists(confOptions['outputdir']):
                os.mkdir(confOptions['outputdir'])
        else:
            confOptions[key] = os.path.abspath(confOptions[key])
            if not os.path.exists(confOptions[key]):
                raise FileNotFoundError('{} not found at {}'.format(name, confOptions[key]))
                # sys.stderr.write('ERROR: {} not found at {}\n'.format(name, confOptions[key]))
                # sys.exit(1)
    elif key == 'outputdir':
        confOptions['outputdir'] = 'out'
        if not os.path.exists(confOptions['outputdir']):
            os.mkdir(confOptions['outputdir'])
    else:
        raise FileNotFoundError('{} ({}) not provided in the config file'.format(name, key))
        # sys.stderr.write('ERROR: {} ({}) not found not provided in the config file\n'.format(name, key))
        # sys.exit(1)
# ==============================================================================
def createControlFreecConfigFile(bamPath, confOptions, freecConfPath):
    '''Creates a config file for Control-FREEC based on files generated in this
    pipeline'''
    freecConf = cp.RawConfigParser()
    freecConf.optionxform = lambda option: option
    freecConf['general'] = {
        'chrLenFile' : confOptions['data']['chrlenfile'],
        'outputDir'  : '{}/freecOut'.format(os.path.abspath(confOptions['data']['outputdir'])),
        'chrFiles'   : confOptions['data']['chrfiles'],
        'samtools'   : confOptions['tools']['samtools'],
        'ploidy'     : 2,
        'window'     : 5000,
    }
    freecConf['sample'] = {
        'mateFile'        : bamPath,
        'inputFormat'     : 'BAM',
        'mateOrientation' : 0
    }
    with open(freecConfPath, 'w') as configfile:
        freecConf.write(configfile)
# ==============================================================================
def main(config_path):
    # Process config file
    config_path = os.path.abspath(config_path)
    config      = cp.ConfigParser(allow_no_value = True)
    config.read(config_path)
    confTools   = config['tools']
    confData    = config['data']
    tools = {
        'bowtie2'    : 'Bowtie2',
        'freec'      : 'Control-FREEC',
        'gatk'       : 'GATK',
        'picard'     : 'Picard',
        'samtools'   : 'Samtools',
        'trimmomatic': 'Trimmomatic'
    }
    data = {
        'adapters'    : 'Adapters file',
         'dbsnp'      : 'dbSNP vcf file',
         'geneset'    : 'Gene set bed file (genes of interest)',
         'index'      : 'Reference Bowtie index prefix',
         'inputfiles' : 'Paired end read file list',
         'outputdir'  : 'Output directory',
         'reference'  : 'Reference genome file'
    }

    # Check for required tools
    for key in tools:
        checkTool(key, tools[key], confTools)
    # Check for required data files
    for key in data:
        checkDataOption(key, data[key], confData)
        
    # ==========================================================================
    # Trim fastq files
    inputFiles = getlist(confData['inputfiles'])
    read1      = inputFiles[0]
    read2      = inputFiles[1]
    tread1     = '{1}_trimmed.fq'.format(confData['outputdir'], os.path.splitext(read1)[0])
    tread2     = '{1}_trimmed.fq'.format(confData['outputdir'], os.path.splitext(read2)[0])
    sread1     = '{1}_unused.fq'.format(confData['outputdir'], os.path.splitext(read1)[0])
    sread2     = '{1}_unused.fq'.format(confData['outputdir'], os.path.splitext(read2)[0])
    
    tcmd       = ['java', '-jar', confTools['trimmomatic'], 'PE', '-phred33', read1, read2, tread1,
                  sread1, tread2, sread2, 'ILLUMINACLIP:{0}:2:30:10'.format(confData['adapters']),
                  'LEADING:0', 'TRAILING:0', 'SLIDINGWINDOW:4:15', 'MINLEN:36']
            
    trun = subprocess.Popen(tcmd, shell=False)
    trun.wait() 
    
    if trun.returncode != 0:
        sys.stderr.write('Fastq trimming failed; Exiting program')
        sys.exit(1)
    
    # ==========================================================================
    # Align the reads using Bowtie2
    sam_path = '{1}.sam'.format(confData['outputdir'], os.path.splitext(tread1)[0])
    bcmd     = [ confTools['bowtie2'], '-x', confData['index'], '-S', 
                 sam_path, '-p', '1' , '-1', tread1, '-2', tread2]
    
    brun = subprocess.Popen(bcmd, shell=False)
    brun.wait()
    
    if brun.returncode != 0:
        sys.stderr.write('Bowtie2 failed; Exiting program')
        sys.exit(1)
    
    # ==========================================================================
    # Add read group information
    add_path = '{0}/{1}_RG.bam'.format(confData['outputdir'], 
        os.path.splitext(os.path.basename(sam_path))[0])
    acmd     = ['java', '-Xmx1g', '-jar', confTools['picard'], 
                'AddOrReplaceReadGroups', 'I='+sam_path , 'O='+add_path, 
                'SORT_ORDER=coordinate', 'RGID=Test', 'RGLB=ExomeSeq', 
                'RGPL=Illumina', 'RGPU=HiSeq2500', 'RGSM=Test', 
                'RGCN=AtlantaGenomeCenter', 'RGDS=ExomeSeq', 
                'RGDT=2016-08-24', 'RGPI=null', 
                'RGPG=Test', 'RGPM=Test', 'CREATE_INDEX=true']
    
    arun = subprocess.Popen(acmd, shell=False)
    arun.wait()
    
    if arun.returncode != 0:
        print('Picard add read groups failed; Exiting program')
        sys.exit(1)
    
    # ==========================================================================
    # Mark PCR duplicates
    dup_path = '{0}/{1}_MD.bam'.format(confData['outputdir'], 
        os.path.splitext(os.path.basename(sam_path))[0])
    met_path = '{0}/{1}_MD.metrics'.format(confData['outputdir'], 
        os.path.splitext(os.path.basename(sam_path))[0])
    mdcmd    = ['java', '-Xmx1g', '-jar', confTools['picard'], 'MarkDuplicates',
                'I='+add_path, 'O='+dup_path, 'METRICS_FILE='+met_path, 
                'REMOVE_DUPLICATES=false', 'ASSUME_SORTED=true', 
                'CREATE_INDEX=true']
    
    mdrun = subprocess.Popen(mdcmd, shell=False)
    mdrun.wait()
    
    if mdrun.returncode != 0:
        print('Picard mark duplicate failed; Exiting program')
        sys.exit(1)

    # ==========================================================================
    # Fix mate information
    fix_path = '{0}/{1}_FM.bam'.format(confData['outputdir'], 
        os.path.splitext(os.path.basename(sam_path))[0])
    fcmd     = ['java', '-Xmx1g', '-jar', confTools['picard'], 
                'FixMateInformation', 'I='+dup_path, 'O='+fix_path, 
                'ASSUME_SORTED=true', 'ADD_MATE_CIGAR=true', 
                'CREATE_INDEX=true']

    frun = subprocess.Popen(fcmd, shell=False)
    frun.wait()
    
    if frun.returncode != 0:
        print('Picard fix mate information failed; Exiting program')
        sys.exit(1)
    
    # ==========================================================================
    # Run realigner target creator
    interval_path = '{0}/{1}.intervals'.format(confData['outputdir'], 
        os.path.splitext(os.path.basename(sam_path))[0]) 
    trcmd         = ['java', '-jar', confTools['gatk'], 
                     '-T', 'RealignerTargetCreator', '-o', interval_path, 
                     '-nt', '1', '-I', fix_path, '-R', confData['reference'], 
                     '-known', confData['dbsnp']]
    
    trrun = subprocess.Popen(trcmd, shell=False)
    trrun.wait()
    
    if trrun.returncode != 0:
        print('Realigner Target creator failed; Exiting program')
        sys.exit(1)

    # =========================================================================
    # Run indel realigner
    ral_path = '{0}/{1}_IR.bam'.format(confData['outputdir'], 
        os.path.splitext(os.path.basename(sam_path))[0])
    recmd    = ['java', '-jar', confTools['gatk'], '-T', 'IndelRealigner',
                '--targetIntervals', interval_path, '-o', ral_path,
                '-I', fix_path, '-R', confData['reference']]

    rerun = subprocess.Popen(recmd, shell=False)
    rerun.wait()

    if rerun.returncode != 0:
        print('Indel realigner creator failed; Exiting program')
        sys.exit(1)

    # ==========================================================================
    # Base quality score recalibration
    bqs_path = '{0}/{1}.table'.format(confData['outputdir'], 
        os.path.splitext(os.path.basename(sam_path))[0])
    bqscmd = ['java', '-jar', confTools['gatk'], '-T', 'BaseRecalibrator', 
              '-R', confData['reference'], '-I', ral_path, '-o', bqs_path, 
              '-nct', '1', '-cov', 'ReadGroupCovariate', 
              '-knownSites', confData['dbsnp']]

    bqsrun = subprocess.Popen(bqscmd, shell=False)
    bqsrun.wait()

    if bqsrun.returncode != 0:
        print('Base quality score recalibrator failed; Exiting program')
        sys.exit(1)
    
    # ==========================================================================
    # Print Reads (generate final BAM)
    fbam_path = '{0}/{1}_final.bam'.format(confData['outputdir'], 
        os.path.splitext(os.path.basename(sam_path))[0])
    prcmd     = ['java', '-jar', confTools['gatk'], '-T', 'PrintReads', 
                 '-R', confData['reference'], '-I', ral_path, 
                 '-o', fbam_path, '-BQSR', bqs_path, '-nct', '1']

    prrun = subprocess.Popen(prcmd, shell=False)
    prrun.wait()

    if prrun.returncode != 0:
        print('Print reads failed; Exiting program')
        sys.exit(1)
    
    # ==========================================================================
    # Coverage stats calculation from the final BAM file
    covfile_path = '{0}/{1}_{2}_coverage.tsv'.format(
        confData['outputdir'], 
        os.path.splitext(os.path.basename(confData['geneset']))[0], 
        os.path.splitext(os.path.basename(fbam_path))[0])
    
    calculateCoverageStats(fbam_path, confData['geneset'], covfile_path, 
                           confTools['samtools'])
    
    # ==========================================================================
    # Haplotype caller
    vcf_path = '{0}/variants.vcf'.format(confData['outputdir'], 
        os.path.splitext(os.path.basename(sam_path))[0])
    hcmd     = ['java', '-jar', confTools['gatk'], '-T', 'HaplotypeCaller', 
                '-R', confData['reference'], '-I', fbam_path, '--dbsnp',
                 confData['dbsnp'], '-o', vcf_path, '-nct', '1', 
                '-gt_mode', 'DISCOVERY']
            
    hrun = subprocess.Popen(hcmd, shell=False)
    hrun.wait()
    
    if hrun.returncode != 0:
        print('Haplotype caller failed; Exiting program')
        sys.exit(1)
    
    # ==========================================================================
    # Variants filtering step
    vcf_filtered_path = '{0}/variants.filtered.vcf'.format(
        confData['outputdir'], os.path.splitext(os.path.basename(sam_path))[0])
    fcmd = ['java', '-jar', confTools['gatk'], '-T', 'SelectVariants', 
            '-R', confData['reference'], '-V', vcf_path, 
            '-select', 'DP>=25 && QUAL>=30', '-o', vcf_filtered_path]
            
    frun = subprocess.Popen(fcmd, shell=False)
    frun.wait()
    
    if frun.returncode != 0:
        print('Variants filtering failed; Exiting program')
        sys.exit(1)
        
    # ==========================================================================
    # CNVs calling
    freecConfPath = os.path.abspath('{0}/freec_conf.txt'.format(confData['outputdir']))
    
    createControlFreecConfigFile(os.path.abspath(fbam_path), config, freecConfPath)
    
    freecConf = cp.ConfigParser(allow_no_value = True)
    freecConf.read(freecConfPath)
    
    if not os.path.exists(freecConf['general']['outputDir']):
        os.mkdir(freecConf['general']['outputDir'])
        
    freecCmd = [confTools['freec'], '-conf', freecConfPath]
    
    freecRun = subprocess.Popen(freecCmd, shell=False)
    freecRun.wait()
    
    if freecRun.returncode != 0:
        print('CNVs calling failed; Exiting program')
        sys.exit(1)
        
    # ====================================================================
    # CNVs plots
    # cat /data2/AHCG2017FALL/bin/FREEC/scripts/makeGraph.R | 
    # R --slave --args 2 outtest/freecOut/reads_1_trimmed_final.bam_ratio.txt
    ratioFilePath = os.path.abspath('{0}/{1}_ratio.txt'.format(freecConf['general']['outputDir'],
        os.path.basename(fbam_path)))
    catCmd  = ['cat', confTools['makegraph']]
    plotCmd = ['R', '--slave', '--args', freecConf['general']['ploidy'], 
               ratioFilePath]
    
    catRun  = subprocess.Popen(catCmd, stdout=subprocess.PIPE, shell=False)
    plotRun = subprocess.Popen(plotCmd, stdin=catRun.stdout, shell=False)
    catRun.stdout.close()
    plotRun.wait()
    catRun.wait()
    
    if plotRun.returncode != 0:
        print('CNVs plotting failed; Exiting program')
        sys.exit(1)
    # ====================================================================

    print('Variant call pipeline completed')
    print('VCF file can be found at {0}'.format(vcf_path))
    return 0;

if __name__ == '__main__':

    parser = argparse.ArgumentParser(
        formatter_class=argparse.RawTextHelpFormatter,
        epilog = textwrap.dedent('''
config file format and options:
  [tools]
  bowtie2     = <Path to Bowtie2 executable>
  freec       = <Path to Control-FREEC executable>
  gatk        = <Path to GATK jar file>
  makegraph   = <Path to Control-FREEC makeGraph.R script (Usually in the folder scripts at the Control-FREEC root dir)>
  picard      = <Path to Picard jar file>
  samtools    = <Path to Samtools executable>
  trimmomatic = <Path to Trimmomatic jar file>
  
  [data]
  adapters    = <Path to adapters fasta file to perform sequence trimming with Trimmomatic>
  chrlenfile  = <Path to file with chromosome lengths for Control-FREEC>
  chrfiles    = <Path to the directory with chromosomes fasta files for Control-FREEC>
  dbsnp       = <Path to dbSNP vcf file for GATK>
  geneset     = <Path to the bed file with genes of interest to calculate coverage statistics>
  index       = <Path to the prefix of the reference Bowtie2 index>
  reference   = <Path to the reference genome fasta file>
  inputfiles  = <List of paired end read files (comma sparated)>
  outputdir   = <Path to the output directory>
  
config file details:
  - Sections [tools] and [data] are required
  - All the options are required except for those that correspond to executable
    files 
  - Options for executable files (bowtie2, freec, and samtools) can be omitted 
    if they are available in the system through the $PATH variable.
'''))
    parser.add_argument('-c', '--config', dest = 'config_path', type = str, 
        help = 'Path to config file')
        
    
    # Enrichment kit (Nextera Enrichment Capture kit)
    # Download samples given an SRA accession
    
    if len(sys.argv) == 1:
        parser.print_help()
        print('--------------------------------------------------------------------------------')
        print('{0} - ERROR: One or more required arguments are missing.\n'.format(parser.prog))
        sys.exit(1)
    
    args = parser.parse_args()
    
    sys.exit(main(args.config_path))

