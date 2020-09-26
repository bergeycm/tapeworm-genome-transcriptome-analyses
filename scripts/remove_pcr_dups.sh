#!/bin/bash

# ----------------------------------------------------------------------------------------
# --- Remove PCR duplicates
# ----------------------------------------------------------------------------------------

IND=$1

READ1=data/WGS/${IND}_TRIM.R1.fastq
READ2=${READ1/.R1/.R2}

# Temporarily unzip
gunzip -c $READ1.gz > $READ1
gunzip -c $READ2.gz > $READ2

wget https://raw.githubusercontent.com/linneas/condetri/master/filterPCRdupl.pl \
     -O scripts/tmp_filterPCRdupl_$IND.pl

perl scripts/tmp_filterPCRdupl_$IND.pl \
    -fastq1=$READ1 \
    -fastq2=$READ2 \
    -prefix=data/WGS/${IND}_TRIM \
    -cmp=30

gzip -c data/WGS/${IND}_TRIM_uniq1.fastq > data/WGS/${IND}_TRIM.R1.rmdup.fastq.gz
gzip -c data/WGS/${IND}_TRIM_uniq2.fastq > data/WGS/${IND}_TRIM.R2.rmdup.fastq.gz

# While here, get length distribution
awk 'NR%4 == 2 {lengths[length($0)]++} END {for (l in lengths) {print l, lengths[l]}}' \
    $READ1 > ${READ1/.fastq/.length}
awk 'NR%4 == 2 {lengths[length($0)]++} END {for (l in lengths) {print l, lengths[l]}}' \
    $READ2 > ${READ2/.fastq/.length}

rm $READ1 $READ2
rm data/WGS/${IND}_TRIM_uniq*.fastq

exit
