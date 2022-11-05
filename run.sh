#!/bin/bash

fastq=$1   # dir of output of minknow (such as: barcodeXX, barcodeXX is a directory that stores multi-fastq files)
oudir=$2   

fullname=${fastq##*/}
stem=${fullname%.*}

refer=/home/renz/Projects/nanopore/NanoForenSNP/config_file/ref_472.fa
snpnb=/home/renz/Projects/nanopore/NanoForenSNP/config_file/ref_472.cfg

mkdir -p $oudir/fastqs
mkdir -p $oudir/alignment
mkdir -p $oudir/snp_out

cat ${fastq}/*.fastq > $oudir/fastqs/${stem}.fq

minimap2 -ax map-ont ${refer} $oudir/fastqs/${stem}.fq | samtools view -S -b | samtools sort > $oudir/alignment/$stem.bam

samtools index $oudir/alignment/$stem.bam

python /home/renz/Projects/nanopore/NanoForenSNP/NanoForenSNP.py getgeno --ref ${refer} --snp ${snpnb} --bam $oudir/alignment/$stem.bam --out $oudir/snp_out/${stem}