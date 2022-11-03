#!/bin/bash

fastq=$1
oudir=$2

fullname=${fastq##*/}
stem=${fullname%.*}
bam=$stem.bam

refer=/home/renz/Projects/nanopore/NanoForenSNP/config_file/ref_472.fa
snpnb=/home/renz/Projects/nanopore/NanoForenSNP/config_file/ref_472.cfg

mkdir -p $oudir/alignment
minimap2 -ax map-ont ${refer} ${fastq} | samtools view -S -b | samtools sort > $oudir/alignment/${bam}
samtools index $oudir/alignment/${bam}

python /home/renz/Projects/nanopore/NanoForenSNP/NanoForenSNP.py getgeno --ref ${refer} --snp ${snpnb} --bam $oudir/alignment/${bam} --out $oudir/snp_out/${stem}