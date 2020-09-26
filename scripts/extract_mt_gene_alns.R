#!/usr/bin/env Rscript

# ----------------------------------------------------------------------------------------
# --- Extract protein coding genes from mt genome alignment
# ----------------------------------------------------------------------------------------

library(ape)

# Align Hymenolepididae and Anoplocephalidae, separately from
#   Taeniidae, Dipylidiidae, and Diphyllobothridae

# --- Process alignment of Taeniidae, Dipylidiidae, and Diphyllobothridae

aln = read.FASTA("results/mt_genomes.non_clade_AH.aln.fa")

tsol.idx = which(names(aln) == "NC_004022-Taenia-solium")

tsol.seq = as.character(aln[tsol.idx])
tsol.gaps    = which(tsol.seq[[1]] == "-")
tsol.notgaps = which(tsol.seq[[1]] != "-")

tsol.mapping = data.frame(orig.coord = 1:length(tsol.notgaps),
                          new.coord  = tsol.notgaps)

anno = read.table("data/annotations_Tsolium_NC_004022.txt",
    header=TRUE)

anno = merge(anno, tsol.mapping, by.x="start", by.y="orig.coord")
anno = merge(anno, tsol.mapping, by.x="end",   by.y="orig.coord")

names(anno)[which(grepl("new.coord", names(anno)))] = c("new.start", "new.end")

aln.m = as.matrix(aln)

dir.create("results/mt_genes", showWarnings=FALSE)

tmp = lapply(1:nrow(anno), function (gene.idx) {
    gene      = anno$gene     [gene.idx]
    new.start = anno$new.start[gene.idx]
    new.end   = anno$new.end  [gene.idx]

    gene.aln = aln.m[,new.start:new.end]

    write.FASTA(gene.aln, file=paste0("results/mt_genes/non-AH_gene_", gene, ".aln.fa"))
})

# --- Process alignment of Hymenolepididae and Anoplocephalidae

aln = read.FASTA("results/mt_genomes.clade_AH.aln.fa")

mexp.idx = which(names(aln) == "KX121041-Moniezia-expansa")

mexp.seq = as.character(aln[mexp.idx])
mexp.gaps    = which(mexp.seq[[1]] == "-")
mexp.notgaps = which(mexp.seq[[1]] != "-")

mexp.mapping = data.frame(orig.coord = 1:length(mexp.notgaps),
                          new.coord  = mexp.notgaps)

anno = read.table("data/annotations_Mexpansa_KX121041.txt",
    header=TRUE)

anno = merge(anno, mexp.mapping, by.x="start", by.y="orig.coord")
anno = merge(anno, mexp.mapping, by.x="end",   by.y="orig.coord")

names(anno)[which(grepl("new.coord", names(anno)))] = c("new.start", "new.end")

aln.m = as.matrix(aln)

dir.create("results/mt_genes", showWarnings=FALSE)

tmp = lapply(1:nrow(anno), function (gene.idx) {
    gene      = anno$gene     [gene.idx]
    new.start = anno$new.start[gene.idx]
    new.end   = anno$new.end  [gene.idx]

    gene.aln = aln.m[,new.start:new.end]

    write.FASTA(gene.aln, file=paste0("results/mt_genes/AH_gene_", gene, ".aln.fa"))
})

# --- Combine genes from each assembly

tmp = lapply(anno$gene, function (gene) {

    system(paste0("cat results/mt_genes/non-AH_gene_", gene, ".aln.fa ",
                       "results/mt_genes/AH_gene_", gene, ".aln.fa >",
                       "results/mt_genes/gene_", gene, ".fa"))
})
