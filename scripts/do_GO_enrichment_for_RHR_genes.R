#!/usr/bin/env Rscript

# ----------------------------------------------------------------------------------------
# --- Search for over-represented GO terms in top RHR genes
# ----------------------------------------------------------------------------------------

options(stringsAsFactors=FALSE)

library(gprofiler2)

# ----------------------------------------------------------------------------------------

bl = read.table("results/gene_tree_distance_ratios_human_nonhuman.txt",
    header=FALSE)

names(bl) = c("emu.id", "human.brlen", "nonhuman.brlen",
              "brlen.ratio", "dists.to.out")

# Remove samples with ratio > 5 or < 0.2 (assuming these are sequencing errors)
bl = bl[bl$brlen.ratio <= 5 & bl$brlen.ratio >= 0.2,]

bl = bl[order(-bl$brlen.ratio),]

# ----------------------------------------------------------------------------------------

brlen.cutoff.percentile = 0.99

bl.high.cutoff = quantile(bl$brlen.ratio,     brlen.cutoff.percentile, na.rm=TRUE)
bl.low.cutoff  = quantile(bl$brlen.ratio, 1 - brlen.cutoff.percentile, na.rm=TRUE)

is.outlier.high = bl$brlen.ratio > bl.high.cutoff
is.outlier.low  = bl$brlen.ratio < bl.low.cutoff

bl.high = bl[which(is.outlier.high),]
bl.low  = bl[which(is.outlier.low), ]

bl.to.print = bl[, 1:4]
write.table(bl.to.print, file="results/gene_tree_distance_ratios_human_nonhuman.tbl.txt",
    sep="\t", quote=FALSE, row.names=FALSE)

# Old way of doing things before switching to R package:

#   # Use g:Profiler's tool for functional profiling, g:GOSt
#   # https://biit.cs.ut.ee/gprofiler/gost
#
#   # Input top 1% of RHR genes and bottom 1% of RHR genes
#   # and using list of all genes with branch length info as BG
#
#   # Organism is Echinococcus multilocularis (PRJEB122)
#   # "Statistical domain scope" is "Custom over annotated genes" (paste in BG list)
#   # "Significance threshold" is "Benjamini-Hochberg FDR"
#   # "User threshold" is 0.1
#
#   # ("Show all results" unchecked)
#   # ("Ordered query" and "Run as multiquery" unchecked)
#   # ("Measure underrepresentation" unchecked)
#   # ("Numeric IDs treated as" left blank
#
#   # Save CSV output:
#   # reports/outlier_high_RHR.gProfiler_results.csv
#   # reports/outlier_low_RHR.gProfiler_results.csv

# ----------------------------------------------------------------------------------------

all.genes = bl.to.print$emu.id[!is.na(bl.to.print$emu.id)]

# Top RHR genes
gprofiler.top_RHR.0.01 = gost(bl.high$emu.id,
    organism = "ecmultprjeb122", ordered_query = FALSE,
    multi_query = FALSE, significant = FALSE, exclude_iea = FALSE,
    measure_underrepresentation = FALSE, evcodes = FALSE,
    user_threshold = 0.2, correction_method = "fdr",
    domain_scope = "annotated",
    custom_bg = all.genes,
    numeric_ns = "NO_NAMESPACE", sources = NULL, as_short_link = FALSE)

res.high = gprofiler.top_RHR.0.01$result[order(gprofiler.top_RHR.0.01$result$p_value),]

write.csv(res.high[,c(3,4,5,6,9,10,11)],
    file="reports/outlier_high_RHR.gProfiler_results.csv")

# Bottom RHR genes
gprofiler.low_RHR.0.01 = gost(bl.low$emu.id,
    organism = "ecmultprjeb122", ordered_query = FALSE,
    multi_query = FALSE, significant = FALSE, exclude_iea = FALSE,
    measure_underrepresentation = FALSE, evcodes = FALSE,
    user_threshold = 0.2, correction_method = "fdr",
    domain_scope = "annotated",
    custom_bg = all.genes,
    numeric_ns = "NO_NAMESPACE", sources = NULL, as_short_link = FALSE)

res.low = gprofiler.low_RHR.0.01$result[order(gprofiler.low_RHR.0.01$result$p_value),]

write.csv(res.low[,c(3,4,5,6,9,10,11)],
    file="reports/outlier_low_RHR.gProfiler_results.csv")

# --- Find genes significant in high but not significant in low

res.both = merge(res.high, res.low, by=c("term_id", "term_name"), suffixes=c(".high", ".low"))
res.only_high = res.both[res.both$p_value.low > 0.2 & res.both$p_value.high < 0.2,]

res.only_high = res.only_high[order(res.only_high$p_value.high),]
