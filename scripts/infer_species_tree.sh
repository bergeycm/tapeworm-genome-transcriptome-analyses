#!/bin/bash

# ----------------------------------------------------------------------------------------
# --- Infer species tree for aTRAM loci
# ----------------------------------------------------------------------------------------

module load raxml/2.7.8

TREE_DIR=results/phylogeny/concat/trees/
mkdir -p $TREE_DIR

rm -f RAxML_info.concat_atram

raxmlHPC-SSE3 \
    -s results/phylogeny/concat/combined.phy \
    -n concat_atram \
    -m GTRGAMMA \
    -p 86795 \
    -o T19

mv RAxML_parsimonyTree.concat_atram $TREE_DIR
mv RAxML_log.concat_atram           $TREE_DIR
mv RAxML_result.concat_atram        $TREE_DIR
mv RAxML_bestTree.concat_atram      $TREE_DIR
mv RAxML_info.concat_atram          $TREE_DIR

# Do 100 rapid bootstrap searches, 20 ML searches and
#  return the best ML tree with support values

rm -f RAxML_info.concat_atram_boot

raxmlHPC-SSE3 \
    -s results/phylogeny/concat/combined.phy \
    -n concat_atram_boot \
    -f a \
    -m GTRGAMMA \
    -p 12345 \
    -x 12345 \
    -# 100 \
     -o T19

mv RAxML_bestTree.concat_atram_boot                 $TREE_DIR
mv RAxML_info.concat_atram_boot                     $TREE_DIR
mv RAxML_bipartitionsBranchLabels.concat_atram_boot $TREE_DIR
mv RAxML_bipartitions.concat_atram_boot             $TREE_DIR
mv RAxML_bootstrap.concat_atram_boot                $TREE_DIR

# ----------------------------------------------------------------------------------------

# Fix taxon names

sed -f scripts/fix_taxon_names.sed \
    $TREE_DIR/RAxML_bestTree.concat_atram > \
    $TREE_DIR/RAxML_bestTree.concat_atram.sp.tre

sed -f scripts/fix_taxon_names.sed \
    $TREE_DIR/RAxML_bipartitions.concat_atram_boot > \
    $TREE_DIR/RAxML_bipartitions.concat_atram_boot.sp.tre

exit
