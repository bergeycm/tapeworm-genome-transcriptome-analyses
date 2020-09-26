#!/bin/bash

# ----------------------------------------------------------------------------------------
# --- Fix FASTAs from aTRAM, add outgroup, filter, infer gene tree
# ----------------------------------------------------------------------------------------

module load samtools
module load emboss/6.6.0
module load gblocks
module load mrbayes

GENE=$1    # E.g. EmuJ_000612800.1

IN_FAS=exon_stitching_echinococcus/final/${GENE/./_}.stitched_exons.fasta

DIR=results/phylogeny/Echin_orthologs_aTRAM/$GENE
mkdir -p $DIR

# --- Make copy of FASTA with duplicate headers removed

awk '!seen[$0]++' $IN_FAS | \
    sed -e "s/^\(>[^\\.]*\)\\..*/\\1/" > \
    $DIR/$GENE.fa

# --- Get outgroup seq

PROTEOME=genomes/echinococcus_multilocularis
PROTEOME=$PROTEOME/echinococcus_multilocularis.PRJEB122.WBPS12.CDS_transcripts.fa

samtools faidx $PROTEOME $GENE > \
    $DIR/Emultilocularis.fa

cat $DIR/Emultilocularis.fa | sed -e "s/^>.*/>Emulti/" >> $DIR/$GENE.fa

# ----------------------------------------------------------------------------------------

muscle -in $DIR/$GENE.fa -out $DIR/$GENE.aln.fa

#Gblocks $DIR/$GENE.aln.fa -t=d    
#mv $DIR/$GENE.aln.fa-gb $DIR/$GENE.flt.fa

python scripts/trim_to_echinococcus_ref.py $DIR/$GENE.aln.fa

seqret $DIR/$GENE.flt.fa $DIR/$GENE.flt.nex -osf nexus
cat data/MrBayes_block.txt >> $DIR/$GENE.flt.nex

mb $DIR/$GENE.flt.nex
