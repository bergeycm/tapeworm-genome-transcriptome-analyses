#!/bin/bash

# ----------------------------------------------------------------------------------------
# --- Run aTRAM
# ----------------------------------------------------------------------------------------

module load blast
module load velvet
ATRAM_PATH=/storage/work/cxb585/bin/aTRAM/  # Turn into module 

# pip3 install --upgrade --ignore-installed --user biopython
# etc. for everything in requirements.txt

IND=$1
PROTS=$2          # FASTA file of reference proteome protein sequences
OUTPUT_PREFIX=$3  # For naming output folder

mkdir -p results/atram/atram_output_$OUTPUT_PREFIX/$IND

python $ATRAM_PATH/atram.py \
        -b results/atram/atram_libraries/$IND \
        -Q $PROTS \
        -p \
        -i 10 \
        --cpus 8 \
        --kmer 17 \
        -o results/atram/atram_output_$OUTPUT_PREFIX/$IND/atram \
        --log-file $PROTS.$IND.atram.log \
        -a velvet \
        --min-contig-length 50 \
        --evalue 0.01

# Create marker file to indicate that this proteome FASTA hunk was processed
touch $PROTS.$IND.atram_run_marker.txt

exit
