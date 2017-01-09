#!/bin/bash

###################################################
#       ChIPseq analysis pipeline                 #
#       Michelle Percharde, PhD 2016              #
#                                                 #
#                v vhg19                          #
###################################################

#~~~~~~~~~~EDIT BELOW THIS LINE ~~~~~~~~~~~~~~~~~~#

#Mm10

# Pipeline to take input dir with files "sample.fq" or "sample.fq.gz", outputs sorted bams
# Doesn't work if input files are .fastq!

#usage: ./runChIPseqHS.sh [options] [-i path/to/folder/]

#FOLDERS NEEDED IN ROOT:
  #raw/ (where raw files are)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

flagcheck=0

while getopts ':gchi:' flag; do
    case ${flag} in
      i) dir=$OPTARG
        flagcheck=1 ;;
        g) gz='true' ;;
        h) echo ""
           echo "Usage: $0 [-h] [-g] [-i <path/to/files/>]"
           echo ""
           echo "    -h        Help mode, prints Usage"
           echo "    -g        Input FASTQ is .gz compressed"
           echo "    -i        input file directory/. use "./" for current dir (not recommended)"
           echo ""
           flagcheck=1
           exit 1 ;;
        \?)
          echo ""
          echo "Invalid option, type -h for help"
          echo ""
          flagcheck=1
          exit 1 ;;
    esac
done

if [ "$flagcheck" == 0 ] ; then
  echo ""
  echo "Incorrect or no -i specified, please read help -h and try again!"
  echo ""
  exit 1
fi

if [ "$gz" == "true" ]; then
  echo ""
  echo "your file is compressed, trim and align will be run on .gz files"
fi

mkdir -p HStrimmed/fastqc/
mkdir -p HSaligned/
mkdir -p HSalignment_summaries/
mkdir -p HSsorted_bam/
mkdir -p HSdedup_bam/

for file in "$dir"* ; do
    echo ""
    if [ "$gz" == "true" ]; then
      name=$(basename $file .fq.gz)
      trimfile=${name}_trimmed.fq.gz
    else
      name=$(basename $file .fq)
      trimfile=${name}_trimmed.fq
    fi
    echo "analysing file: $name, trimmed will be $trimfile"
    echo ""
    echo "1. trimming $name"
    echo ""
    trim_galore --fastqc --fastqc_args " --outdir HStrimmed/fastqc/" -a ATCGGAAGAGCAC $file -o HStrimmed/ #picks up more adapters

    echo ""
    echo "2. aligning $name to hg19 with bowtie"
    echo ""
    # mkdir -p ${name}_aligned/

    (bowtie2 -p 4 -x /data/refs/hg19/genome -U HStrimmed/$trimfile | samtools view -Suo - - | \
    samtools sort - -o HSsorted_bam/${name}.sorted.bam) 2> HSalignment_summaries/${name}_alignment.txt

    echo ""
    echo "3. deduplicating $name bam file"
    echo ""
    samtools rmdup -s HSsorted_bam/${name}.sorted.bam HSdedup_bam/${name}.sorted.dedup.bam

    echo ""
    echo "4. Generating bam index for $name"
    echo ""
    samtools index HSsorted_bam/${name}.sorted.bam HSsorted_bam/${name}.sorted.bai
    samtools index HSdedup_bam/${name}.sorted.dedup.bam HSdedup_bam/${name}.sorted.dedup.bai

    echo ""
    echo "$name DONE!"
    echo ""
done
