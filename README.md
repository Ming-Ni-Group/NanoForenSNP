# ForenSNP
SNP genotype caller (based on [iSNV](https://github.com/generality/iSNV-calling)) for ForenSeq Kit. It can be applied for the reads processed by ForenSeq Kit and sequenced by NGS platform or Nanopore platform. We've already evaluated the accuracy on 94 identity SNPs from 33 samples. We build it with a snakemake file which will make it easy to install and use. 


## Features
1. Genotype caller for 94 iSNPs from Forenseq Kit;
2. Supporting from NGS/Nanopore platform;


## Workflow



## Usage

```bash
ForenSNP=[PATH_TO_ForenSNP]

python $ForenSNP/ForenSNP.py getgeno --ref $PATH_TO_REFERENCE_GENOME --snp $ForenSNP/config_file/rsID.txt --bam $PATH_TO_BAMFILE --id $SAMPLE_ID
```

## Installation

We recommand you to use conda to build the vitual enviroment to run our SNP caller.

```bash
cd $DIR_TO_SAVE_THE_PROGRAME
conda env create -f environment.yml
git clone XXXX

```

1. create a new conda enviroment
```bash
conda create -n 
```
