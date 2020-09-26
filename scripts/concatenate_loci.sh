#!/bin/bash

# ========================================================================================
# --- Concatenate aTRAM loci
# ========================================================================================

module load emboss/6.6.0

PHY_DIR=results/phylogeny

EMUS=(`ls $PHY_DIR/Echin_orthologs_aTRAM/*/*aln.fa | cut -d"/" -f4`)

# ----------------------------------------------------------------------------------------
# --- Check taxon counts
# ----------------------------------------------------------------------------------------

echo > $PHY_DIR/concat_tax_count.txt

for EMU in ${EMUS[*]}; do
    TAX_CT=`grep -c "^>" $PHY_DIR/Echin_orthologs_aTRAM/${EMU}/${EMU}.flt.fa`

    echo -e "$EMU\t$TAX_CT" >> $PHY_DIR/concat_tax_count.txt
    echo $TAX_CT

done | sort | uniq -c

TO_INCLUDE=(`awk '{ if ($2 >= 6) print $1 }' $PHY_DIR/concat_tax_count.txt`)

# ----------------------------------------------------------------------------------------
# --- Make NEXUS files ready to concatenate
# ----------------------------------------------------------------------------------------

mkdir -p $PHY_DIR/concat
mkdir -p $PHY_DIR/gtst

for EMU in ${TO_INCLUDE[*]}; do
    ORIG_NEXUS=$PHY_DIR/Echin_orthologs_aTRAM/${EMU}/${EMU}.flt.nex
    NEW_NEXUS=$PHY_DIR/concat/${EMU}.nocomment.nex
    # Already done when inferring gene tree
    # seqret $PHY_DIR/${EMU}.flt.fa $PHY_DIR/concat/${EMU}.aln.nex -osf nexus
    grep -v "\[" $ORIG_NEXUS > $NEW_NEXUS
done

# ----------------------------------------------------------------------------------------
# --- Remove empty NEXUS files
# ----------------------------------------------------------------------------------------

for NEX in `ls $PHY_DIR/concat/*.nocomment.nex`; do
    echo "Checking $NEX..."
    CHAR_CT=`grep -c "nchar" $NEX`
    if [[ $CHAR_CT == 0 ]]; then
        echo "- Removing $NEX"
        rm $NEX
    fi
done

# ----------------------------------------------------------------------------------------
# --- Do concatenation, clean up, and convert to PHYLIP format
# ----------------------------------------------------------------------------------------

module load python/2.7.14-anaconda5.0.1
python scripts/concatenate_aTRAM_loci.py

# Makes $PHY_DIR/concat/combined.fas directly, so can skip next seqret command
###seqret $PHY_DIR/concat/combined.nex $PHY_DIR/concat/combined.fa -osf fasta

# Convert N's to dashes in FASTA alignment
#sed -i"BACK" -e "s/N/-/g" $PHY_DIR/concat/combined.fas

module load gblocks/0.91
# -t=d   # DNA
# -b5=h  # half gaps allowed
Gblocks $PHY_DIR/concat/combined.fas -t=d -b5=h
mv $PHY_DIR/concat/combined.fas-gb $PHY_DIR/concat/combined.flt.fas

seqret $PHY_DIR/concat/combined.flt.fas $PHY_DIR/concat/combined.phy -osf phylip

exit
