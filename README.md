# NanoForenSNP v0.0.1
This script is for Forensic SNP genotyping in the manuscript:  *Nanopore sequencing of forensic STRs and SNPs using Verogen’s ForenSeq DNA Signature Prep Kit and MinION.*

NanoForenSNP (based on [iSNV](https://github.com/generality/iSNV-calling)) is for ForenSeq Kit. It can be applied for the reads processed by ForenSeq Kit and sequenced by NGS platform or Nanopore platform. We've already evaluated the accuracy on 94 identity SNPs from 33 samples. All but one, rs983283, were genotyped correctly.

We build it with python and make it easy to install and use. 



## Installation

We recommend you to use `conda` to build the environment.

1. Download the code from our repo:

```
git clone https://github.com/Ming-Ni-Lab/NanoForenSNP.git
```

2. Build the environment:


```bash
cd NanoForenSNP
conda env create -f env.yml
conda activate nanoforensnp
```

3. Test

We provide a script named `run.sh <fastq-file> <output-dir>` to perform nanoforensnp. Before running, users need to speify the path of some files.
Here is the example.

```bash
#!/bin/bash

fastq=$1
oudir=$2

fullname=${fastq##*/}
stem=${fullname%.*}
bam=$stem.bam

refer=<nanoforensnp-dir>/config_file/ref_472.fa
snpnb=<nanoforensnp-dir>/config_file/ref_472.cfg

mkdir -p $oudir/alignment
minimap2 -ax map-ont ${refer} ${fastq} | samtools view -S -b | samtools sort > $oudir/alignment/${bam}
samtools index $oudir/alignment/${bam}

python <nanoforensnp-dir>/NanoForenSNP.py getgeno --ref ${refer} --snp ${snpnb} --bam $oudir/alignment/${bam} --out $oudir/snp_out/${stem}
```
----------------------------

## python script (old version)

### Quick start

```bash
conda activate nanoforensnp

nanoforensnp="[PATH_TO_NanoForenSNP]"

## required arguments
python $nanoforensnp/nanoforensnp.py getgeno \
--ref $nanoforensnp/config_file/snpRef.fa \
--snp $nanoforensnp/config_file/config.txt \
--bam $BAM_FILE \
--id $SAMPLE_NAME

cd output/$SAMPLE_NAME
```



## Contributing

The NanoForenSNP is based on the EBOV intrahost single nucleotide variation ([iSNV](https://github.com/generality/iSNV-calling])) calling which is created by Dr. Ni.  

