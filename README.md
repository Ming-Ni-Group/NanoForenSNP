# ForenSNP
SNP genotype caller (based on [iSNV](https://github.com/generality/iSNV-calling)) for ForenSeq Kit. It can be applied for the reads processed by ForenSeq Kit and sequenced by NGS platform or Nanopore platform. We've already evaluated the accuracy on 94 identity SNPs from 33 samples. We build it with python and make it easy to install and use. 


## Features
1. Genotype caller for 94 iSNPs from Forenseq Kit;
2. Supporting from NGS/Nanopore platform;



## Installation

We recommend you to use `conda` to build the ForenSNP environment.

1. Download the code from our repo:

```
git clone https://github.com/Ming-Ni-Lab/ForenSNP.git
```

2. Build the environment:


```bash
cd ForenSNP
conda env create -f env.yml
```

3. Test

```
conda activate ForenSNP
python ForenSNP.py getgeno -h
conda deactivate
```



## Usage

### Quick start



```bash
conda activate ForenSNP

FORENSNP="[PATH_TO_FORENSNP]"

## required arguments
REF_GENOME="[PATH_TO_REFERENCE_GENOME]"
BAM_FILE="[PATH_TO_BAM_FILE]"
SAMPLE_NAME="[SAMPLE_NAME]"

python $FORENSNP/ForenSNP.py getgeno \
--ref $REF_GENOME \
--snp $FORENSNP/config_file/rsID.txt \
--bam $BAM_FILE \
--id $SAMPLE_NAME

cd output/$SAMPLE_NAME
```



### Learn more







## Contributing

The ForenSNP is based on the EBOV intrahost single nucleotide variation ([iSNV](https://github.com/generality/iSNV-calling])) calling which is created by Ming Ni.  

