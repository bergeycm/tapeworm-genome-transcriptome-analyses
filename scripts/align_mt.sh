#!/bin/bash

# ----------------------------------------------------------------------------------------
# --- Align mt genomes
# ----------------------------------------------------------------------------------------

# --- Concatenate all mt genomes and remove N's

cat data/Guo2017_tapeworm_mt_genomes.fa \
    results/MITObim_all.fasta | \
    sed -e "/^[^>]/ s/N//g" > \
    results/mt_genomes.fa

# --- Split FASTA to separately align those with similar mitochondrial structure

# Align Hymenolepididae and Anoplocephalidae, separately from
#   Taeniidae, Dipylidiidae, and Diphyllobothridae

module load samtools
samtools faidx results/mt_genomes.fa

echo "AF314223-Hymenolepis-diminuta"          >  results/clade_AH_samples.txt
echo "NC_028334-Pseudanoplocephala-crawfordi" >> results/clade_AH_samples.txt
echo "NC_029245-Rodentolepis-nana"            >> results/clade_AH_samples.txt
echo "KU236385-Anoplocephala-magna"           >> results/clade_AH_samples.txt
echo "KR054960-Anoplocephala-perfoliata"      >> results/clade_AH_samples.txt
echo "mt_T19"                                 >> results/clade_AH_samples.txt
echo "KX121040-Moniezia-benedeni"             >> results/clade_AH_samples.txt
echo "KX121041-Moniezia-expansa"              >> results/clade_AH_samples.txt
echo "mt_T26"                                 >> results/clade_AH_samples.txt
echo "NC_028164-Drepanidotaenia-lanceolata"   >> results/clade_AH_samples.txt

grep "^>" results/mt_genomes.fa | \
    sed -e "s/^>//" | \
    grep -v -f results/clade_AH_samples.txt > results/non_AH_samples.txt

xargs samtools faidx results/mt_genomes.fa \
    < results/clade_AH_samples.txt \
    > results/mt_genomes.clade_AH.fa

xargs samtools faidx results/mt_genomes.fa \
    < results/non_AH_samples.txt \
    > results/mt_genomes.non_clade_AH.fa

# --- Do alignment

module load clustalw2/2.1

for FA in results/mt_genomes.clade_AH.fa results/mt_genomes.non_clade_AH.fa; do

    clustalw2 \
        -INFILE=$FA \
        -OUTFILE=${FA/.fa/.aln.fa} \
        -OUTPUT=FASTA \
        -ALIGN
done

# --- Pull out genes

Rscript scripts/extract_mt_gene_alns.R

# --- Align each gene now that they've been grabbed from the two alignments

for FA in `ls results/mt_genes/gene_*.fa | grep -v "aln" | grep -v "fix"`; do

    sed -e "/^[^>]/ s/-//g" $FA > ${FA/.fa/.fix.fa}

    clustalw2 \
        -INFILE=${FA/.fa/.fix.fa} \
        -OUTFILE=${FA/.fa/.aln.fa} \
        -OUTPUT=FASTA \
        -ALIGN &
done

# --- Clean up alignments

module load gblocks/0.91
module load emboss/6.6.0

# --- Clean up - Gene alignments

# Replicate the webserver's less stringent options:
# - Allow smaller final blocks
# - Allow gap positions within the final blocks
# - Allow less strict flanking positions

# Parameters used:
# - Minimum Number Of Sequences For A Conserved Position: 29; 50% + 1 (default)
# - Minimum Number Of Sequences For A Flanking Position: 29; 50% + 1 (reduced)
# - Maximum Number Of Contiguous Nonconserved Positions: 8 (default)
# - Minimum Length Of A Block: 5 (reduced)
# - Allowed Gap Positions: With Half (changed from default of None)

for GENE_FA in `ls results/mt_genes/gene_*.aln.fa`; do
    Gblocks $GENE_FA -t=d -b2=29 -b4=5 -b5=half
    mv $GENE_FA-gb ${GENE_FA/.aln.fa/.aln.flt.fa}
done

# Marker that we're done
ls results/mt_genes/*flt.fa > results/mt_genes/gene_aln_list.txt

exit
