#!/bin/bash

# ----------------------------------------------------------------------------------------
# --- Copy combined aTRAM output
# ----------------------------------------------------------------------------------------

PROTEOME=$1

ls results/atram/atram_output_$PROTEOME/ > tmp.taxon_list.$PROTEOME.txt

COMBINED_DIR=results/atram/atram_output_$PROTEOME/combined
mkdir -p $COMBINED_DIR

for FASTA in `find results/atram/atram_output_$PROTEOME/T*/ -name "*filtered*" -type f`; do
    echo $FASTA
    NEW_FILE=`basename $FASTA | sed -e "s/atram.//"`
    cp $FASTA $COMBINED_DIR/$NEW_FILE
done

touch $COMBINED_DIR/copy_done_marker.txt

exit
