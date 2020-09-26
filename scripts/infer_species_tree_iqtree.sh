#!/bin/bash

# ----------------------------------------------------------------------------------------
# --- Infer species tree for aTRAM loci with IQ-TREE
# ----------------------------------------------------------------------------------------

export PATH=$PATH:/storage/home/cxb585/work/bin/iqtree-2.0.6-Linux/bin

TREE_DIR=results/phylogeny/concat/wgs_iq_tree/
mkdir -p $TREE_DIR

# --- Copy all loci in

ALN_DIR=tmp_loci
mkdir -p $ALN_DIR
cp results/phylogeny/concat/*nocomment* $ALN_DIR

# --- Infer a concatenation-based species tree with 1000 ultrafast bootstrap and
# --- an edge-linked partition model
iqtree2 -p $ALN_DIR --prefix wgs_concat -B 1000 -T AUTO

# --- Infer the locus trees
iqtree2 -S $ALN_DIR --prefix wgs_loci -T AUTO

# --- Compute concordance factors
iqtree2 -t wgs_concat.treefile --gcf wgs_loci.treefile -p $ALN_DIR --scf 100 --prefix wgs_concord -T 10

# Placeholder until we know output file names      
touch $TREE_DIR/finished_marker.txt                

exit