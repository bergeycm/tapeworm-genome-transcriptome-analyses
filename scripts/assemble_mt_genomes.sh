#!/bin/bash

# ----------------------------------------------------------------------------------------
# --- Assemble mt genomes
# ----------------------------------------------------------------------------------------

IND=$1

MIRA_DIR=/storage/home/cxb585/work/bin/mira_4.0.2_linux-gnu_x86_64_static/bin
export PATH=$PATH:$MIRA_DIR
MITOBIM_DIR=/storage/home/cxb585/work/bin/MITObim

# --- Make directory in which to work

mkdir -p results/MITObim_$IND
cd results/MITObim_$IND

# --- Reference genome downloaded from GenBank in FASTA format

if [[ "$IND" == "T19" || "$IND" == "T23" || "$IND" == "T26" ]]; then

    # Moniezia expansa mitochondrial COX-1 gene for cytochrome oxidase subunit 1, partial cds
    # GenBank: AB099693.1

    REF_FA=../../genomes/mt_AB099693.1.fa

elif [[ "$IND" == "T17" ]]; then

    # Taenia taeniaeformis mitochondrial CO1 gene for cytochrome c oxidase subunit 1, partial cds
    # GenBank: AB221484.1

    REF_FA=../../genomes/mt_AB221484.1.fa

else

    # Taenia solium mitochondrial cox1 gene for cytochrome c oxidase subunit 1, complete cds
    # GenBank: AB271234.1

    REF_FA=../../genomes/mt_AB271234.1.fa

fi

# --- Reads file

R1_FQ=../../data/WGS/${IND}_TRIM.R1.fastq.gz
R2_FQ=../../data/WGS/${IND}_TRIM.R2.fastq.gz

gunzip -c $R1_FQ > tmp.$IND.R1.fq
gunzip -c $R2_FQ > tmp.$IND.R2.fq

# One liner script taken from https://gist.github.com/nathanhaigh/4544979

paste tmp.$IND.R1.fq tmp.$IND.R2.fq | \
    paste - - - - | \
    awk -v OFS="\n" -v FS="\t" '{ print($1,$3,$5,$7,$2,$4,$6,$8) }' > \
    tmp.$IND.interleaved.fastq

rm tmp.$IND.R1.fq tmp.$IND.R2.fq

# --- Call MITObim in quick mode

# Delete old run results
rm -rf iteration*

GENOME_NAME=`basename $REF_FA | sed -e "s/\..*//"`

perl $MITOBIM_DIR/MITObim.pl \
    -start 1 -end 50 \
    -sample MITObim_mt_$IND \
    -ref $GENOME_NAME \
    -readpool tmp.$IND.interleaved.fastq \
    --quick $REF_FA \
    --paired --clean &> \
     MITObim_log_$IND.txt

rm tmp.$IND.interleaved.fastq

# --- Copy final FASTA out of iteration folder

FASTA=`grep "^Final assembly" MITObim_log_$IND.txt | cut -d":" -f2 | sed -e "s/ //"`

sed -e "s/^>.*/>mt_$IND/" $FASTA > MITObim_mt_$IND-mt-final_noIUPAC.fasta

exit
