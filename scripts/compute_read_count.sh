#!/bin/bash

# ----------------------------------------------------------------------------------------
# --- Compute total number of reads
# ----------------------------------------------------------------------------------------

for FQ in `ls data/WGS/*R1.fastq.gz | grep -iv "trim" | grep "/T"`; do

    IND=`basename $FQ | sed -e "s/\..*//"`

    TOT_LINES=`gunzip -c $FQ | wc -l | cut -d' ' -f1`
    TOT_READS=`echo "$TOT_LINES / 4" | bc`

    TRIM_LINES=`gunzip -c data/WGS/${IND}_TRIM.R1.fastq.gz | wc -l | cut -d' ' -f1`
    TRIM_READS=`echo "$TRIM_LINES / 4" | bc`

    echo -e "$IND\t$TOT_READS\t$TRIM_READS"

done
