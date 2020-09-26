#!/bin/bash

# ----------------------------------------------------------------------------------------
# --- Make aTRAM libraries
# ----------------------------------------------------------------------------------------

module load blast
#module load velvet
ATRAM_PATH=/storage/work/cxb585/bin/aTRAM/  # Turn into module 

mkdir -p results/atram/atram_libraries

IND=$1

R1=data/WGS/${IND}_TRIM.R1.rmdup.fastq
R2=${R1/R1/R2}

gunzip -c $R1.gz > $R1
gunzip -c $R2.gz > $R2

python $ATRAM_PATH/atram_preprocessor.py \
    --cpus 8 \
    -b results/atram/atram_libraries/$IND \
    --end-1 $R1 \
    --end-2 $R2

# Remove temporarily decompressed files
rm $R1 $R2

exit
