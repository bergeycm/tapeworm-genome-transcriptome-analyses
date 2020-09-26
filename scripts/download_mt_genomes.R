#!/usr/bin/env Rscript

# ========================================================================================
# --- Download tapeworm mt genomes used in Guo 2017 paper
# ========================================================================================

options(stringsAsFactors=FALSE)

library(ape)

# --- Download Guo 2017 genomes from GenBank

mt.acc = read.table("data/Guo2017_S2_tapeworm_mt_genome_accessions.txt", sep="\t")

mt.gb = read.GenBank(mt.acc$V2, species.names=T)

# Add species to name
names(mt.gb) = paste(names(mt.gb), mt.acc$V1)
names(mt.gb) = gsub(" ", "-", names(mt.gb))

write.dna(mt.gb, "data/Guo2017_tapeworm_mt_genomes.fa", format="fasta")
