#!/usr/bin/env Rscript

library(ape)
library(ggplot2)
library(ggrepel)
library(gprofiler2)

# ========================================================================================
# --- Test if human accelerated ev genes are overrep'd for heat-DE genes
# ========================================================================================

bl = read.table("results/gene_tree_distance_ratios_human_nonhuman.txt",
    header=FALSE)

names(bl) = c("emu.id", "human.brlen", "nonhuman.brlen",
              "brlen.ratio", "dists.to.out")

# Remove samples with ratio > 5 or < 0.2 (assuming these are sequencing errors)
bl = bl[bl$brlen.ratio <= 5 & bl$brlen.ratio >= 0.2,]

de.prefix = paste0("/gpfs/group/ghp3/default/Bergey/",
                           "tapeworm_RNAseq_MIRROR/tapeworm_RNAseq/",
                           "tapeworm-rnaseq/results/")

# All heated vs. control
de.all  = read.table(paste0(de.prefix,
                            "DE_genes.heated-vs-ctl.known.all.txt"),
                     header=TRUE)

# 56C vs. control
de.high = read.table(paste0(de.prefix,
                            "DE_genes.56-vs-ctl.known.all.txt"),
                     header=TRUE)

# Add ortholog info to the DE data
# Orthology info comes from putting all Taenia IDs through g:Profiler's g:Orth tool
ortho.info = read.csv("data/orthology_info_gProfiler_tasoliprjna170813_ecmultprjeb122.csv")

de.all  = merge(de.all,  ortho.info[,c(3,5)], by.x="GeneID", by.y="initial_ensg")
de.high = merge(de.high, ortho.info[,c(3,5)], by.x="GeneID", by.y="initial_ensg")

names(de.all) [names(de.all)  == "ortholog_name"] = "emu.id"
names(de.high)[names(de.high) == "ortholog_name"] = "emu.id"

bl$emu.id = gsub("\\..*", "", bl$emu.id)

# Combine branch length and DE data
bl.de.all  = merge(bl, de.all,  by="emu.id")
bl.de.high = merge(bl, de.high, by="emu.id")

# ----------------------------------------------------------------------------------------
# --- Overrep / enrichment testing function
# ----------------------------------------------------------------------------------------

do.overrep.tests = function (bl.de, is.brlen.outlier, is.heatde.outlier, out.prefix) {

    keep.cols = c(1:4, 6:ncol(bl.de))

    # ("Fast" might also be a slow outlier in the inverted control test.)
    fast.de       = unique(bl.de[  is.brlen.outlier  &   is.heatde.outlier,  keep.cols])
    fast.notde    = unique(bl.de[  is.brlen.outlier  & !(is.heatde.outlier), keep.cols])
    notfast.de    = unique(bl.de[!(is.brlen.outlier) &   is.heatde.outlier,  keep.cols])
    notfast.notde = unique(bl.de[!(is.brlen.outlier) & !(is.heatde.outlier), keep.cols])

    fast.de.ct       = nrow(fast.de      )
    fast.notde.ct    = nrow(fast.notde   )
    notfast.de.ct    = nrow(notfast.de   )
    notfast.notde.ct = nrow(notfast.notde)

    # --- Chi-squared test

    sink(paste0(out.prefix, ".chisq.txt"), split=TRUE)

    cat(paste("Heat DE percentage among branch length outlier human ev proteins:\n"))
    cat(paste(fast.de.ct, "of", (fast.de.ct + fast.notde.ct)))
    cat("\n")
    cat(fast.de.ct / (fast.de.ct + fast.notde.ct))
    cat("\n")

    cat(paste("Heat DE percentage among other proteins:\n"))
    cat(paste(notfast.de.ct, "of", (notfast.de.ct + notfast.notde.ct)))
    cat("\n")
    cat(notfast.de.ct / (notfast.de.ct + notfast.notde.ct))
    cat("\n")

    cat(paste("Chi-squared test output:\n"))
    chi.res = chisq.test(c(fast.de.ct, fast.notde.ct, notfast.de.ct, notfast.notde.ct))
    chi.res.str = capture.output(chi.res)
    cat(chi.res.str)
    cat("\n")

    cat(paste("Chi-squared p-value:\n"))
    cat(chi.res$p.value)
    cat("\n")

    sink()

    # --- Hypergeometric test

    m = fast.de.ct + notfast.de.ct         # Genes heat DE
    n = fast.notde.ct + notfast.notde.ct   # Genes NOT heat DE
    k = fast.de.ct + fast.notde.ct         # Human accelerated genes
    x = 0:(fast.de.ct + fast.notde.ct)     # Genes both heat DE and accelerated

    # Use the dhyper built-in function for hypergeometric density
    probabilities = dhyper(x, m, n, k, log = FALSE)

    sink(paste0(out.prefix, ".hypergeometric.txt"), split=TRUE)

    cat(paste("Hypergeometric test p-value:\n"))
    cat(sum(probabilities[fast.de.ct:(fast.de.ct + fast.notde.ct)]))
    cat("\n")

    cat(paste("Fisher's test output:\n"))

    fish.res = fisher.test(matrix(c(fast.de.ct,    fast.notde.ct,
                                    notfast.de.ct, notfast.notde.ct),
                                  nrow=2))
    fish.res.str = capture.output(fish.res)
    cat(fish.res.str)
    cat("\n")

    cat(paste("Fisher's test p-value:\n"))
    cat(fish.res$p.value)
    cat("\n")

    sink()

    # --- Wilcox test

    sink(paste0(out.prefix, ".wilcox.txt"), split=TRUE)

    cat(paste("Wilcox test output:\n"))

    wilcox.test(bl.de$PValue ~ is.brlen.outlier)

    cat(paste("Median and Mean of heat DE p-values for branch length outliers:\n"))
    cat(paste("Median", median(bl.de[is.brlen.outlier,]$PValue), "\n"))
    cat(paste("Mean",   mean(bl.de[is.brlen.outlier,]$PValue),   "\n"))
    cat("\n")

    cat(paste("Median and Mean of heat DE p-values for other proteins:\n"))
    cat(paste("Median", median(bl.de[!is.brlen.outlier,]$PValue), "\n"))
    cat(paste("Mean",   mean(bl.de[!is.brlen.outlier,]$PValue),   "\n"))
    cat("\n")

    sink()

    return(NA)
}

# ----------------------------------------------------------------------------------------
# --- All heated vs. control
# ----------------------------------------------------------------------------------------

# --- Figure out which genes are upregulated

bl.de.all$direction = "None"
bl.de.all$direction[rowSums(bl.de.all[,7:12] > 0) >= 5] = "Up"
bl.de.all$direction[rowSums(bl.de.all[,7:12] > 0) <= 1] = "Down"

# --- Of the accelerated human ev proteins, are all-heat DE proteins over-represented?

brlen.cutoff.percentile = 0.95

bl.high.cutoff  = quantile(bl$brlen.ratio, brlen.cutoff.percentile, na.rm=TRUE)
de.high.cutoff  = 0.1
de.logFC.cutoff = 0.5

is.brlen.outlier  = bl.de.all$brlen.ratio > bl.high.cutoff
is.heatde.outlier = abs(bl.de.high$logFC) > de.logFC.cutoff &
                    bl.de.all$FDR <= de.high.cutoff &
                    bl.de.all$direction == "Up"

do.overrep.tests (bl.de.all, is.brlen.outlier, is.heatde.outlier,
                  paste0("results/human_accel_heat_DE.allheated.",
                         "brlen_perc", brlen.cutoff.percentile))

# --- INVERSE TEST: Of the *slowest* evolving human ev proteins,
# ---               are all-heat DE genes over-represented?

brlen.cutoff.percentile = 0.05

bl.low.cutoff   = quantile(bl$brlen.ratio, brlen.cutoff.percentile, na.rm=TRUE)
de.high.cutoff  = 0.1
de.logFC.cutoff = 0.5

is.brlen.outlier  = bl.de.all$brlen.ratio < bl.low.cutoff
is.heatde.outlier = abs(bl.de.high$logFC) > de.logFC.cutoff &
                    bl.de.all$FDR <= de.high.cutoff &
                    bl.de.all$direction == "Up"

do.overrep.tests (bl.de.all, is.brlen.outlier, is.heatde.outlier,
                  paste0("results/human_slowed_heat_DE.allheated.",
                         "brlen_perc", brlen.cutoff.percentile))

# ----------------------------------------------------------------------------------------
# --- 56C vs. control
# ----------------------------------------------------------------------------------------

# --- Figure out which genes are upregulated

bl.de.high$direction = "Down"
bl.de.high$direction[bl.de.high$logFC > 0] = "Up"

# --- Of the accelerated human ev proteins, are highest-heat DE proteins over-represented?

brlen.cutoff.percentile = 0.95

bl.high.cutoff  = quantile(bl$brlen.ratio, brlen.cutoff.percentile, na.rm=TRUE)
de.high.cutoff  = 0.1
de.logFC.cutoff = 0.5

is.brlen.outlier  = bl.de.high$brlen.ratio > bl.high.cutoff
is.heatde.outlier = abs(bl.de.high$logFC) > de.logFC.cutoff &
                    bl.de.high$FDR < de.high.cutoff &
                    bl.de.high$direction == "Up"

do.overrep.tests (bl.de.high, is.brlen.outlier, is.heatde.outlier,
                  paste0("results/human_accel_heat_DE.highestheated.",
                         "brlen_perc", brlen.cutoff.percentile))

# --- INVERSE TEST: Of the *slowest* evolving human ev proteins,
# ---               are highest-heat DE genes over-represented?

brlen.cutoff.percentile = 0.05

bl.low.cutoff   = quantile(bl$brlen.ratio, brlen.cutoff.percentile, na.rm=TRUE)
de.high.cutoff  = 0.1
de.logFC.cutoff = 0.5

is.brlen.outlier  = bl.de.high$brlen.ratio < bl.low.cutoff
is.heatde.outlier = abs(bl.de.high$logFC) > de.logFC.cutoff &
                    bl.de.high$FDR <= de.high.cutoff &
                    bl.de.high$direction == "Up"

do.overrep.tests (bl.de.high, is.brlen.outlier, is.heatde.outlier,
                  paste0("results/human_slowed_heat_DE.highestheated.",
                         "brlen_perc", brlen.cutoff.percentile))

# ----------------------------------------------------------------------------------------
# --- Make list of interesting genes for GO overrep
# ----------------------------------------------------------------------------------------

bl.high.cutoff = quantile(bl$brlen.ratio, 0.95, na.rm=TRUE)
bl.low.cutoff  = quantile(bl$brlen.ratio, 0.05, na.rm=TRUE)

de.cutoff = 0.1
de.logFC.cutoff = 0.5

# High RHR, up-regulated - all heat vs. control
dbl.outliers.high.all = bl.de.all[bl.de.all$brlen.ratio > bl.high.cutoff &
                                  abs(bl.de.all$logFC.group56) > de.logFC.cutoff &
                                  bl.de.all$FDR <= de.high.cutoff &
                                  bl.de.all$direction == "Up",]

# Low RHR, up-regulated - all heat vs. control
dbl.outliers.low.all = bl.de.all[bl.de.all$brlen.ratio < bl.low.cutoff &
                                  abs(bl.de.all$logFC.group56) > de.logFC.cutoff &
                                  bl.de.all$FDR <= de.high.cutoff &
                                  bl.de.all$direction == "Up",]

# High RHR, up-regulated - 56 vs. control
dbl.outliers.high.high = bl.de.high[bl.de.high$brlen.ratio > bl.high.cutoff &
                                    abs(bl.de.high$logFC) > de.logFC.cutoff &
                                    bl.de.high$FDR <= de.high.cutoff &
                                    bl.de.high$direction == "Up",]

# Low RHR, up-regulated - 56 vs. control
dbl.outliers.low.high = bl.de.high[bl.de.high$brlen.ratio < bl.low.cutoff &
                                   abs(bl.de.high$logFC) > de.logFC.cutoff &
                                   bl.de.high$FDR <= de.high.cutoff &
                                   bl.de.high$direction == "Up",]

# Sort by branch length ratio
dbl.outliers.high.all  = dbl.outliers.high.all[order(dbl.outliers.high.all$brlen.ratio,
    decreasing=TRUE),]
dbl.outliers.high.high = dbl.outliers.high.high[order(dbl.outliers.high.high$brlen.ratio,
    decreasing=TRUE),]

dbl.outliers.low.all  = dbl.outliers.low.all[order(dbl.outliers.low.all$brlen.ratio,
    decreasing=FALSE),]
dbl.outliers.low.high = dbl.outliers.low.high[order(dbl.outliers.low.high$brlen.ratio,
    decreasing=FALSE),]

cols.to.keep.all = c(1,6,2:4,12,14,15)
dbl.outliers.high.all  = dbl.outliers.high.all[,cols.to.keep.all]
dbl.outliers.low.all   = dbl.outliers.low.all [,cols.to.keep.all]

cols.to.keep.high = c(1,6,2:4,7,9,10)
dbl.outliers.high.high = dbl.outliers.high.high[,cols.to.keep.high]
dbl.outliers.low.high  = dbl.outliers.low.high [,cols.to.keep.high]

names(dbl.outliers.high.all)  = names(dbl.outliers.low.all) =
names(dbl.outliers.high.high) = names(dbl.outliers.low.high) = c(
    "ID.Echinococcus_multilocularis", "ID.Taenia_solium",
    "Branch_length.human_infective", "Branch_length.non_human_infective",
    "Branch_length.ratio",
    "logFC", "PValue", "FDR")

# All heated vs. control
write.table(dbl.outliers.high.all,
    file="reports/outlier_high_RHR_DE_all_heated_vs_control.txt",
    sep="\t", row.names=FALSE, col.names=TRUE)
write.table(dbl.outliers.high.all[,1],
    file="reports/outlier_high_RHR_DE_all_heated_vs_control.justEmuID.txt",
    sep="\t", row.names=FALSE, col.names=FALSE, quote=FALSE)

write.table(dbl.outliers.low.all,
    file="reports/outlier_low_RHR_DE_all_heated_vs_control.txt",
    sep="\t", row.names=FALSE, col.names=TRUE)
write.table(dbl.outliers.low.all[,1],
    file="reports/outlier_low_RHR_DE_all_heated_vs_control.justEmuID.txt",
    sep="\t", row.names=FALSE, col.names=FALSE, quote=FALSE)

# 56 vs. control
write.table(dbl.outliers.high.high,
    file="reports/outlier_high_RHR_DE_56_vs_control.txt",
    sep="\t", row.names=FALSE, col.names=TRUE)
write.table(dbl.outliers.high.high[,1],
    file="reports/outlier_high_RHR_DE_56_vs_control.justEmuID.txt",
    sep="\t", row.names=FALSE, col.names=FALSE, quote=FALSE)

write.table(dbl.outliers.low.high,
    file="reports/outlier_low_RHR_DE_56_vs_control.txt",
    sep="\t", row.names=FALSE, col.names=TRUE)
write.table(dbl.outliers.low.high[,1],
    file="reports/outlier_low_RHR_DE_56_vs_control.justEmuID.txt",
    sep="\t", row.names=FALSE, col.names=FALSE, quote=FALSE)

# --- Also make BG list
write.table(unique(de.all[,11]), file="reports/BG.justEmuID.DE_all_heated.txt",
    sep="\t", row.names=FALSE, col.names=FALSE, quote=FALSE)
write.table(unique(de.high[,6]), file="reports/BG.justEmuID.DE_56.txt",
    sep="\t", row.names=FALSE, col.names=FALSE, quote=FALSE)

# ----------------------------------------------------------------------------------------
# --- Then, do GO overrepresentation analysis
# ----------------------------------------------------------------------------------------

# --- RHR+heat

# Old way of doing things before switching to R package:

#   # Use g:Profiler's tool for functional profiling, g:GOSt
#   # https://biit.cs.ut.ee/gprofiler/gost
#
#   # All heated vs. control:
#   #   Using high RHR + heat DE (all) input:
#   #     reports/outlier_high_RHR_DE_all_heated_vs_control.justEmuID.txt
#   #   and low RHR + heat DE (all) input:
#   #     reports/outlier_high_RHR_DE_all_heated_vs_control.justEmuID.txt
#   #   and using list of all genes with branch length info as BG:
#   #     reports/BG.justEmuID.DE_all_heated.txt
#
#   # Highest heated vs. control:
#   #   Using high RHR + heat DE (56) input:
#   #     reports/outlier_high_RHR_DE_56_vs_control.justEmuID.txt
#   #   and low RHR + heat DE (56) input:
#   #     reports/outlier_high_RHR_DE_56_vs_control.justEmuID.txt
#   #   and using list of all genes with branch length info as BG:
#   #     reports/BG.justEmuID.DE_56.txt
#
#   # Organism is Echinococcus multilocularis (PRJEB122)
#   # "Statistical domain scope" is "Custom over annotated genes" (paste in BG list)
#   # "Significance threshold" is "Benjamini-Hochberg FDR"
#   # "User threshold" is 0.2
#
#   # ("Show all results" unchecked)
#   # ("Ordered query" and "Run as multiquery" unchecked)
#   # ("Measure underrepresentation" unchecked)
#   # ("Numeric IDs treated as" left blank
#
#   # Save CSV output:
#   # reports/outlier_high_RHR_all_heated_vs_control.gProfiler_results.csv
#   # reports/outlier_low_RHR_all_heated_vs_control.gProfiler_results.csv
#   # reports/outlier_high_RHR_56_vs_control.gProfiler_results.csv
#   # reports/outlier_low_RHR_56_vs_control.gProfiler_results.csv

# --- All heated vs. control

all.genes = as.character(unique(de.all[,11]))

# high RHR + heat DE
gprofiler.top_RHR_and_heat_DE = gost(dbl.outliers.high.all$ID.Echinococcus_multilocularis,
    organism = "ecmultprjeb122", ordered_query = FALSE,
    multi_query = FALSE, significant = FALSE, exclude_iea = FALSE,
    measure_underrepresentation = FALSE, evcodes = FALSE,
    user_threshold = 0.2, correction_method = "fdr",
    domain_scope = "annotated",
    custom_bg = all.genes,
    numeric_ns = "NO_NAMESPACE", sources = NULL, as_short_link = FALSE)

res.high = gprofiler.top_RHR_and_heat_DE$result[order(gprofiler.top_RHR_and_heat_DE$result$p_value),]

write.csv(res.high[,c(3,4,5,6,9,10,11)],
    file="reports/outlier_high_RHR_all_heated_vs_control.gProfiler_results.csv")

# low RHR + predicted/known drug targets genes
gprofiler.low_RHR_and_heat_DE = gost(dbl.outliers.low.all$ID.Echinococcus_multilocularis,
    organism = "ecmultprjeb122", ordered_query = FALSE,
    multi_query = FALSE, significant = FALSE, exclude_iea = FALSE,
    measure_underrepresentation = FALSE, evcodes = FALSE,
    user_threshold = 0.2, correction_method = "fdr",
    domain_scope = "annotated",
    custom_bg = all.genes,
    numeric_ns = "NO_NAMESPACE", sources = NULL, as_short_link = FALSE)

res.low = gprofiler.low_RHR_and_heat_DE$result[order(gprofiler.low_RHR_and_heat_DE$result$p_value),]

write.csv(res.low[,c(3,4,5,6,9,10,11)],
    file="reports/outlier_low_RHR_all_heated_vs_control.gProfiler_results.csv")

# --- 56 vs. control

all.genes = as.character(unique(de.high[,6]))

# high RHR + heat DE
gprofiler.top_RHR_and_heat_DE = gost(dbl.outliers.high.high$ID.Echinococcus_multilocularis,
    organism = "ecmultprjeb122", ordered_query = FALSE,
    multi_query = FALSE, significant = FALSE, exclude_iea = FALSE,
    measure_underrepresentation = FALSE, evcodes = FALSE,
    user_threshold = 0.2, correction_method = "fdr",
    domain_scope = "annotated",
    custom_bg = all.genes,
    numeric_ns = "NO_NAMESPACE", sources = NULL, as_short_link = FALSE)

res.high = gprofiler.top_RHR_and_heat_DE$result[order(gprofiler.top_RHR_and_heat_DE$result$p_value),]

write.csv(res.high[,c(3,4,5,6,9,10,11)],
    file="reports/outlier_high_RHR_56_vs_control.gProfiler_results.csv")

# low RHR + predicted/known drug targets genes
gprofiler.low_RHR_and_heat_DE = gost(dbl.outliers.low.high$ID.Echinococcus_multilocularis,
    organism = "ecmultprjeb122", ordered_query = FALSE,
    multi_query = FALSE, significant = FALSE, exclude_iea = FALSE,
    measure_underrepresentation = FALSE, evcodes = FALSE,
    user_threshold = 0.2, correction_method = "fdr",
    domain_scope = "annotated",
    custom_bg = all.genes,
    numeric_ns = "NO_NAMESPACE", sources = NULL, as_short_link = FALSE)

res.low = gprofiler.low_RHR_and_heat_DE$result[order(gprofiler.low_RHR_and_heat_DE$result$p_value),]

write.csv(res.low[,c(3,4,5,6,9,10,11)],
    file="reports/outlier_low_RHR_56_vs_control.gProfiler_results.csv")

# ----------------------------------------------------------------------------------------

# --- Plot DE p-value against branch length ratio       

bl.de.high$to.label = FALSE
bl.de.high$to.label[-log(bl.de.high$FDR, base=10) > 15 & log(bl.de.high$brlen.ratio, base=2) > 0.5] = TRUE

bl.de.high$label = ""
bl.de.high$label[bl.de.high$emu.id == "EmuJ_000417100"] = "Cytosolic malate dehydrogenase"
bl.de.high$label[bl.de.high$emu.id == "EmuJ_000065300"] = "Universal stress protein"
bl.de.high$label[bl.de.high$emu.id == "EmuJ_000413500"] = "Sodium:potassium dependent ATPase beta subunit"
bl.de.high$label[bl.de.high$emu.id == "EmuJ_000634800"] = "L-lactate dehydrogenase"
bl.de.high$label[bl.de.high$emu.id == "EmuJ_000859400"] = "Ras protein Rap 1b"
bl.de.high$label[bl.de.high$emu.id == "EmuJ_000873800"] = "Universal stress protein"
bl.de.high$label[bl.de.high$emu.id == "EmuJ_000886400"] = "Tubulin alpha chain"

p = ggplot(bl.de.high[bl.de.high$direction == "Up",],
        aes(-1 * log(FDR, base=10), log(brlen.ratio, base=2), col=-logFC)) +
    geom_point() +
    geom_label_repel(data=bl.de.high[bl.de.high$direction == "Up" & bl.de.high$to.label,],
        aes(label=label), nudge_y = 1, max.iter = 200000) +
    ylim(c(-2.5, 2)) +
    xlab("-log10(adjusted p-value) for heat DE test") +
    ylab("log2(human vs. non-human infective branch length ratio)") +
    theme_bw() +
    #scale_size(range = c(0.1, 3)) +
    theme(legend.position = "none")

ggsave(p, file="reports/branch_vs_heat_DE_pval_56-vs-ctl.pdf")

quit()
