
"""


    - NanoForenSNP ver0.0.1
    - Written by Ren Zilin; zilin.ren@outlook.com
    - For reads from the forenseq kit and nanopore platform 


"""

import os
import sys
from pathlib import Path
import argparse

## local pack
from myScripts import myGenotypeCaller

def run_mpileup(samtools_path, ref_path, snp_bed_path, bam_file_path, output_pileup):
    cmd = "%s mpileup -f %s -l %s %s -o %s" % (samtools_path, ref_path, snp_bed_path, bam_file_path, output_pileup)
    sys.stdout.write('## Running command: `%s`\n' % cmd)
    os.system(cmd)
    return

def run_isnv(perl_path, output_pileup, output_ntfreq):
    isnv_path = os.path.join(os.path.dirname(os.path.realpath(__file__)), 'myScripts/myISNV.pl')
    cmd = "%s %s -in %s -out %s" % (perl_path, isnv_path, output_pileup, output_ntfreq)
    sys.stdout.write('## Running command: `%s`\n' % cmd)
    os.system(cmd)
    return

def getgeno(args):
    ## tools' path
    samtools_path = args.SAMTOOLS
    perl_path     = args.PERL
    
    ## files' path
    ref_path      = args.ref
    snp_bed_path  = args.snp
    bam_file_path = args.bam
    
    ## parameter set
    homo_cutoff  = args.cutoff
    
    ## out dir
    outdir = Path(args.out)
    Path(outdir.parent).mkdir(parents=True, exist_ok=True)
    
    sampleID = outdir.stem
    ## out files
    output_pileup = '%s.mpileup' % (outdir, )
    output_ntfreq = '%s.ntfreq'  % (outdir, )
    output_genotype = '%s.csv'   % (outdir, )
    
    ## work flow
    run_mpileup(samtools_path, ref_path, snp_bed_path, bam_file_path, output_pileup)
    run_isnv(perl_path, output_pileup, output_ntfreq)
    myGenotypeCaller.get_alleles(sampleID, output_ntfreq, snp_bed_path, output_genotype, homo_cutoff)
    
    return





if __name__ == '__main__':
    
    
    
    parser = argparse.ArgumentParser(description='SNP genotype caller for data from forenseq kit.', 
                                     usage='''python NanoForenSNP.py <command> [<args>]

Available commands are:    
    getgeno    simple genotype caller based on iSNV
    ''')
    subparsers = parser.add_subparsers(help='sub-command help')
    subparsers.required = True

    ############## getgeno
    getgeno_parser = subparsers.add_parser('getgeno', help='simple genotype caller based on iSNV')
    
    ## required arguments
    getgeno_parser.add_argument('--ref', required = True, help='path to genome reference', default = None)
    getgeno_parser.add_argument('--snp', required = True, help='path to the snp information file, see the USAGE for format', default = None)
    getgeno_parser.add_argument('--bam', required = True, help='path to the bam file', default = None)
    # getgeno_parser.add_argument('--id',  required = True, help='sample name for output file', default = 'sample')
    getgeno_parser.add_argument('--out', required = True, help='outdir')
    
    
    ## optional
    # getgeno_parser.add_argument('--subdir', required = False, help='subdir', default = '')
    getgeno_parser.add_argument('--SAMTOOLS', required = False, help='path to samtools', default = 'samtools')
    getgeno_parser.add_argument('--PERL', required = False, help='path to perl', default = 'perl')
    getgeno_parser.add_argument('--cutoff', required = False, 
                                help='the cutoff for homo/heterozygous; recommand range 0.1-0.15 for nanopore data', default = 0.1)
 
    getgeno_parser.set_defaults(func=getgeno)


    ## run
    args = parser.parse_args()
    args.func(args)
