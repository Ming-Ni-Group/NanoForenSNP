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

```
cd test
minimap2 -ax map-ont ../config_file/snpRef.fa alignment/barcode01.part.fastq | samtools view -S -b | samtools sort > ./alignment/p1b01.test.bam

samtools index ./alignment/p1b01.test.bam

python ../nanoforensnp.py getgeno --ref ../config_file/snpRef.fa --snp ../config_file/config.txt --id p1b01test --bam alignment/p1b01.test.bam
```



## Usage

### Quick start

```bash
conda activate nanoforensnp

NanoForenSNP="[PATH_TO_NanoForenSNP]"

## required arguments
REF_GENOME="[PATH_TO_REFERENCE_GENOME]"
BAM_FILE="[PATH_TO_BAM_FILE]"
SAMPLE_NAME="[SAMPLE_NAME]"

python $NanoForenSNP/nanoforensnp.py getgeno \
--ref $NanoForenSNP/config_file/snpRef.fa \
--snp $NanoForenSNP/config_file/config.txt \
--bam $BAM_FILE \
--id $SAMPLE_NAME

cd output/$SAMPLE_NAME
```



## Contributing

The NanoForenSNP is based on the EBOV intrahost single nucleotide variation ([iSNV](https://github.com/generality/iSNV-calling])) calling which is created by Dr. Ni.  

