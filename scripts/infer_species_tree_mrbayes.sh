#!/bin/bash

# ----------------------------------------------------------------------------------------
# --- Infer species tree for aTRAM loci with MrBayes
# ----------------------------------------------------------------------------------------

# NEXUS input file is from concatenate_aTRAM_loci.py
# $PHY_DIR/concat/combined.nex

module load mrbayes/3.2.6
module load emboss

ALN_DIR=results/phylogeny/concat/
mkdir -p $ALN_DIR

seqret $ALN_DIR/combined.nex $ALN_DIR/combined.mb.nex -osf nexus

cat data/MrBayes_block.txt >> $ALN_DIR/combined.mb.nex

mb $ALN_DIR/combined.mb.nex

# ----------------------------------------------------------------------------------------

# Fix taxon names

###sed -f scripts/fix_taxon_names.sed \
###    $TREE_DIR/RAxML_bestTree.concat_atram > \
###    $TREE_DIR/RAxML_bestTree.concat_atram.sp.tre  
###
###sed -f scripts/fix_taxon_names.sed \
###    $TREE_DIR/RAxML_bipartitions.concat_atram_boot > \
###    $TREE_DIR/RAxML_bipartitions.concat_atram_boot.sp.tre  

exit
