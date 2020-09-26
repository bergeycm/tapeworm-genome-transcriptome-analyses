#!/bin/bash

# ----------------------------------------------------------------------------------------
# --- Trim reads
# ----------------------------------------------------------------------------------------

THIS_IND=$1

echo "Trimming for individual $THIS_IND";

module load trimmomatic/0.36

ADAPTERS=`dirname $TRIMMOMATIC`/adapters/TruSeq3-PE.fa

java -jar $TRIMMOMATIC PE -phred33 \
    -threads 8 \
    data/WGS/${THIS_IND}.R1.fastq.gz \
    data/WGS/${THIS_IND}.R2.fastq.gz \
    data/WGS/${THIS_IND}_TRIM.R1.fastq.gz \
    data/WGS/${THIS_IND}_trim_unpaired.R1.fastq.gz \
    data/WGS/${THIS_IND}_TRIM.R2.fastq.gz \
    data/WGS/${THIS_IND}_trim_unpaired.R2.fastq.gz \
    ILLUMINACLIP:$ADAPTERS:2:30:10 \
    LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36

exit;
