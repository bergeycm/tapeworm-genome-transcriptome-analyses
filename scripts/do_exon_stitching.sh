#!/bin/bash

# ========================================================================================
# --- Do exon stitching for reference proteone hunk
# --- Based on original script by Julie Allen
# ========================================================================================

PROTEOME_NAME=$1    # E.g. schistosoma
PROTEOME_FASTA=$2   # E.g. genomes/schistosoma_mansoni/
                    #      schistosoma_mansoni.PRJEA36577.
                    #      WBPS12.protein.oneline.fa

OVERLAP=10;

TAXON_LIST=tmp.taxon_list.$PROTEOME_NAME.txt
FASTA_DIR=results/atram/atram_output_$PROTEOME_NAME/combined/

module load exonerate

mkdir -p exon_stitching_$PROTEOME_NAME
cd exon_stitching_$PROTEOME_NAME

# ----------------------------------------------------------------------------------------
# --- Edit aTRAM assemblies to remove empty files, lose weird characters, and
# ---  add a unique number to the end of each contig
# --- Output:
# ---  cleaned FASTA files, *.ed.fasta, in same passed $fasta_dir directory
# ---  FASTA files for all reference genes, *.reference.fasta
# ----------------------------------------------------------------------------------------

echo "--- Editing aTRAM assemblies... [`date`]"

# perl scripts/exon_stitching_edit_fastas.pl [fasta_dir] [gene_list_out] [proteome_query]

# Remove old *.ed.fasta files
find ../$FASTA_DIR -name "*.ed.fasta" -type f -delete

GENE_LIST=../results/atram/atram_output_$PROTEOME_NAME/gene_list.txt

# NOTE: Runs in parallel with 20 processors

perl ../scripts/exon_stitching_edit_fastas.pl \
    ../results/atram/atram_output_$PROTEOME_NAME/combined/ \
    $GENE_LIST \
    ../$PROTEOME_FASTA

echo "--- Finished editing aTRAM assemblies. [`date`]"

# ----------------------------------------------------------------------------------------
# --- Create list of all assembly files for each gene
# ----------------------------------------------------------------------------------------

echo "--- Creating lists of assembly files for each gene... [`date`]"

while read -r line  || [[ -n "$line" ]]; do
    gene=${line/./_}
    gene=${gene%".reference.fasta"};
    echo "Processing gene [$gene]"
    find ../$FASTA_DIR -name "*$gene*.ed.fasta" > $gene.list.txt
done < $GENE_LIST

echo "--- Finished creating lists of assembly files for each gene. [`date`]"

# ----------------------------------------------------------------------------------------
# --- For each proteome reference and each individual (taxon), run exonerate
# ----------------------------------------------------------------------------------------

# This creates a results file with the exonerate information for each exon
# The results file is sorted by taxon and by location of the exon
# Exonerate is used to pull out the exons

echo "--- Running exonerate... [`date`]"

# Remove old CSV files (exonerate output)
rm -r *.csv

while read -r proteome_ref_file || [[ -n "$proteome_ref_file" ]]; do

    gene=${proteome_ref_file%".reference.fasta"}
    gene=${gene/./_}

    echo "Processing reference proteome file: [$proteome_ref_file] for gene [$gene]"

    while read -r tax || [[ -n "$tax" ]]; do

        if [[ "$tax" == "combined" ]]; then
            continue
        fi

        echo "- Taxon [$tax]"

        while read -r assembly || [[ -n "$assembly" ]]; do

            if [[ $assembly = *${tax}_* || $assembly = *$tax.* ]]; then

                echo "Match! $assembly contains $tax"

                # Call exonerate and output gene results
                exonerate \
                    --verbose 0 \
                    --model protein2genome \
                    $proteome_ref_file \
                    $assembly \
                    --showvulgar no \
                    --showalignment no \
                    --ryo "$gene,$tax,%ql,%qal,%qab,%qae,%ti\n" >> $gene.results.csv

                # Call exonerate and output exon results
                exonerate \
                    --verbose 0 \
                    --model protein2genome \
                    $proteome_ref_file \
                    $assembly \
                    --showvulgar no \
                    --showalignment no \
                    --ryo ">$tax,%ti,%qab\n%tcs\n" >> $gene.exons.fasta

                # Sort gene results
                LC_ALL=C sort -t, -k 1,1d -k 2,2d -k 5,5d \
                    $gene.results.csv > $gene.results.sorted.csv

            fi
        done < "$gene.list.txt"
    done < ../$TAXON_LIST
done < $GENE_LIST

echo "--- Finished running exonerate. [`date`]"

# ----------------------------------------------------------------------------------------
# --- Figure out which overlapping contigs will be stiched together
# --- Output:
# ---  CSV files with info on which contigs to be stitched together
# ---      {gene}.overlap.{overlap}.contig_list.csv
# ----------------------------------------------------------------------------------------

echo "--- Getting contigs to stitch... [`date`]"

perl ../scripts/exon_stitching_get_contigs.pl ../exon_stitching_$PROTEOME_NAME/ $OVERLAP

echo "--- Finished getting contigs to stitch. [`date`]"

# ----------------------------------------------------------------------------------------
# --- Stitch contigs together
# --- Output:
# ---  Files of stitched together contigs from exonerate, separated by NNN
# ---      {out_dir}/final/{gene}.stitched_exons.fasta
# ---  Summary file (by gene)
# ---      {out_dir}/final/summary_stats.per.gene.csv
# ----------------------------------------------------------------------------------------

echo "--- Stitching contigs... [`date`]"

perl ../scripts/exon_stitching_stitch_contigs.pl \
    ../exon_stitching_$PROTEOME_NAME \
    $OVERLAP \
    &> exon_stitching.out.txt

echo "--- Finished stitching contigs. [`date`]"

# ----------------------------------------------------------------------------------------
# --- Compute summary stats
# --- Output:
# ---  Summary file (by taxon)
# ---      {out_dir}/final/summary_stats.per.taxon.csv
# ----------------------------------------------------------------------------------------

echo "--- Computing summary stats... [`date`]"

perl ../scripts/exon_stitching_summary_stats.pl \
    ../exon_stitching_$PROTEOME_NAME \
    &> summary_stat.out.txt

echo "--- Finished computing summary stats. [`date`]"

# ----------------------------------------------------------------------------------------
# --- Clean up
# ----------------------------------------------------------------------------------------

# In this directory
rm *.reference.fasta
rm *.list.txt
# Keep {gene}.exons.fasta ?
# Keep {gene}.results.csv ?
rm *.contig_list.csv

# In final/
# Keep final/{gene}.stitched_exons.fasta
# Keep final/summary_stats.per.gene.csv
# Keep final/summary_stats.per.taxon.csv
