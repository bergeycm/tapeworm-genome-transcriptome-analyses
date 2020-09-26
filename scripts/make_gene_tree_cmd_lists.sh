#!/bin/bash

# ----------------------------------------------------------------------------------------
# --- Infer all gene trees
# ----------------------------------------------------------------------------------------

module load parallel/20170422

EMUS=(`ls exon_stitching_echinococcus/*exons.fasta | \
    sed -e "s/.*echinococcus\/\\(.*\\)\\.exons.*/\\1/" \
        -e "s/_/\\./g" \
        -e "s/EmuJ\\./EmuJ_/"`)

mkdir -p results/phylogeny/gene_tree_commands/
GENE_TREE_CMD_PRE=results/phylogeny/gene_tree_commands/gene_tree_cmds

echo "date;" > $GENE_TREE_CMD_PRE.sh

for EMU in ${EMUS[*]}; do

    DIR=results/phylogeny/Echin_orthologs_aTRAM/$EMU

    # Remove old results
    rm -f $DIR/$EMU*

    echo "sh scripts/align_and_infer_gene_tree.sh $EMU" >> \
        $GENE_TREE_CMD_PRE.sh
done

echo "date" >> $GENE_TREE_CMD_PRE.sh

# --- Split into 24 hunks

NUM_FILES=24
TOT_LINES=`wc -l $GENE_TREE_CMD_PRE.sh | cut -d' ' -f1`
((HUNK_SIZE = (TOT_LINES + $NUM_FILES - 1) / $NUM_FILES))

split -l $HUNK_SIZE -d $GENE_TREE_CMD_PRE.sh ${GENE_TREE_CMD_PRE}_part

exit
