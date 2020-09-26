#!/bin/bash

# ----------------------------------------------------------------------------------------
# --- Estimate ideal kmer size
# ----------------------------------------------------------------------------------------

IND=$1

# kmergenie conflicts with python 3
PYTHONPATH=
module unload python
module load kmergenie

mkdir -p results/kmer_estimation

echo -e "data/WGS/${IND}_TRIM.R1.fastq.gz\ndata/WGS/${IND}_TRIM.R2.fastq.gz" \
    > tmp_${IND}_file_list.txt

kmergenie tmp_${IND}_file_list.txt \
    --diploid -t 8 \
    -o results/kmer_estimation/${IND} \
    --orig-hist

rm tmp_${IND}_file_list.txt

IDEAL_KMER=`sort -n -k2 -r results/kmer_estimation/${IND}.dat | \
    head -n1 | cut -d' ' -f 1`

echo $IDEAL_KMER > results/kmer_estimation/${IND}.bestkmer.txt
