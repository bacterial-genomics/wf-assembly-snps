#!/usr/bin/env bash


# Confirm this is being ran with the default directory structure
if [ ! -f main.nf ]; then
  echo "ERROR: expected to run this script in the default dir structure" >&2
  echo "       unable to find the ./main.nf file" >&2
  exit 1
elif [ ! -d "tests/data" ]; then
  echo "ERROR: expected to run this script in the default dir structure" >&2
  echo "       unable to find the ./tests/data directory" >&2
  exit 1
elif [ ! -w "${PWD}" ]; then
  echo "ERROR: unable to write to the ${PWD} directory, which" >&2
  echo "       is the output path of this test script" >&2
  exit 1
fi

# Create an array of assembly files in the test data dir
shopt -s nullglob
ASM=( ./tests/data/*.fna.gz )
shopt -u nullglob

# Grab assembly files if the path lacks them
if [ ${#ASM[@]} -ne 10 ]; then
  if [ ! -w "tests/data" ]; then
    echo "ERROR: unable to write to the ./tests/data directory" >&2
    exit 1
  fi
  wget \
   https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/703/365/GCA_000703365.1_Ec2011C-3609/GCA_000703365.1_Ec2011C-3609_genomic.fna.gz \
   https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/016/766/575/GCA_016766575.1_PDT000040717.5/GCA_016766575.1_PDT000040717.5_genomic.fna.gz \
   https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/003/018/935/GCA_003018935.1_ASM301893v1/GCA_003018935.1_ASM301893v1_genomic.fna.gz \
   https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/012/830/055/GCA_012830055.1_PDT000040719.3/GCA_012830055.1_PDT000040719.3_genomic.fna.gz \
   https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/012/829/335/GCA_012829335.1_PDT000040724.3/GCA_012829335.1_PDT000040724.3_genomic.fna.gz \
   https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/003/018/775/GCA_003018775.1_ASM301877v1/GCA_003018775.1_ASM301877v1_genomic.fna.gz \
   https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/012/829/275/GCA_012829275.1_PDT000040726.3/GCA_012829275.1_PDT000040726.3_genomic.fna.gz \
   https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/016/766/555/GCA_016766555.1_PDT000040728.5/GCA_016766555.1_PDT000040728.5_genomic.fna.gz \
   https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/012/829/195/GCA_012829195.1_PDT000040729.3/GCA_012829195.1_PDT000040729.3_genomic.fna.gz \
   https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/012/829/295/GCA_012829295.1_PDT000040727.3/GCA_012829295.1_PDT000040727.3_genomic.fna.gz \
   -P "${PWD}/tests/data"
fi

# Run the test data set through this workflow
nextflow run main.nf \
 --inpath "${PWD}/tests/data" \
 --outpath "${PWD}/tests/test-output" \
 --reference "${PWD}/tests/data/GCA_000703365.1_Ec2011C-3609_genomic.fna.gz" \
 -resume \
 -ansi-log false
