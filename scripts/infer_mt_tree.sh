#!/bin/bash

# ----------------------------------------------------------------------------------------
# --- Infer mitochondrial tree with IQ-TREE
# ----------------------------------------------------------------------------------------

export PATH=$PATH:/storage/home/cxb585/work/bin/iqtree-2.0.6-Linux/bin

TREE_DIR=results/phylogeny/mt_phylogeny
mkdir -p $TREE_DIR

# --- Infer tree for mt coding genes only

mkdir -p tmp_iqtree_alns
cp results/mt_genes/gene_*.aln.flt.fa tmp_iqtree_alns/

# Infer a concatenation-based species tree with 1000 ultrafast bootstrap and an edge-linked partition model
iqtree2 \
    --redo \
    -p tmp_iqtree_alns \
    --prefix mt_genes_concat \
    -o "NC_002544-Schistosoma-japonicum" \
    -B 1000 \
    -T AUTO

# Infer the locus trees
iqtree2 \
    --redo \
    -S tmp_iqtree_alns \
    --prefix mt_genes_loci \
    -o "NC_002544-Schistosoma-japonicum" \
    -T AUTO

# Compute concordance factors
iqtree2 \
    --redo \
    -t       mt_genes_concat.treefile \
    --gcf    mt_genes_loci.treefile\
    -p       tmp_iqtree_alns \
    --scf    100 \
    --prefix mt_genes_concord \
    -T 8

mv mt_genes_concat.*  $TREE_DIR
mv mt_genes_loci.*    $TREE_DIR
mv mt_genes_concord.* $TREE_DIR

rm -r tmp_iqtree_alns
