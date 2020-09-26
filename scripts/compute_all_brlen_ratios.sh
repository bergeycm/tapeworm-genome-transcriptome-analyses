#!/bin/bash

# ----------------------------------------------------------------------------------------
# --- Compute human/non-human branch length ratio for all gene trees
# ----------------------------------------------------------------------------------------

ls results/phylogeny/Echin_orthologs_aTRAM/*/*.flt.nex.con.tre | \
    sed -e "s/.*aTRAM\/\([^\/]*\)\/.*/\1/" > tmp.emu.list.txt

for EMU in `cat tmp.emu.list.txt`; do
    Rscript scripts/compute_dist_to_outgroup.R $EMU
done > results/gene_tree_distance_ratios_human_nonhuman.txt

rm tmp.emu.list.txt

exit
