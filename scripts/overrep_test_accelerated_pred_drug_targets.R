#!/usr/bin/env Rscript

library(ape)
library(gprofiler2)

# ----------------------------------------------------------------------------------------
# --- Test if human accelerated ev genes are overrep'd for likely drug targets
# ----------------------------------------------------------------------------------------

bl = read.table("results/gene_tree_distance_ratios_human_nonhuman.txt",
    header=FALSE)

names(bl) = c("emu.id", "human.brlen", "nonhuman.brlen",
              "brlen.ratio", "dists.to.out")

# Remove samples with ratio > 5 or < 0.2 (assuming these are sequencing errors)
bl = bl[bl$brlen.ratio <= 5 & bl$brlen.ratio >= 0.2,]

bl$emu.id = gsub("\\..*", "", bl$emu.id)

drug = read.csv("data/Tsai_etal_2013/Tsai_etal_2013_Suppl13.10.csv",
    sep=",", comment.char="#")

names(drug)[names(drug) == "ID"] = "emu.id"

dg = drug[,c("emu.id", "Product",
             "S..mansoni", "T..solium",
             "Drugbank.approved", "Total.score")]

bl.high.cutoff = quantile(bl$brlen.ratio, 0.95, na.rm=TRUE)
dg.high.cutoff = quantile(dg$Total.score, 0.95, na.rm=TRUE)

bl.low.cutoff  = quantile(bl$brlen.ratio, 0.05, na.rm=TRUE)
dg.low.cutoff  = quantile(dg$Total.score, 0.05, na.rm=TRUE)

bl.dg = merge(bl, dg, by="emu.id")

#head(bl.dg[order(-bl.dg$brlen.ratio),])

# ----------------------------------------------------------------------------------------
# --- Overrep / enrichment testing function
# ----------------------------------------------------------------------------------------

do.overrep.tests = function (bl.dg, is.brlen.outlier, is.drug.target, out.prefix) {

    keep.cols = c(1:4, 6:ncol(bl.dg))

    # ("Fast" might also be a slow outlier in the inverted control test.)
    fast.dg       = unique(bl.dg[  is.brlen.outlier  &   is.drug.target,  keep.cols])
    fast.notdg    = unique(bl.dg[  is.brlen.outlier  & !(is.drug.target), keep.cols])
    notfast.dg    = unique(bl.dg[!(is.brlen.outlier) &   is.drug.target,  keep.cols])
    notfast.notdg = unique(bl.dg[!(is.brlen.outlier) & !(is.drug.target), keep.cols])

    fast.dg.ct       = nrow(fast.dg      )
    fast.notdg.ct    = nrow(fast.notdg   )
    notfast.dg.ct    = nrow(notfast.dg   )
    notfast.notdg.ct = nrow(notfast.notdg)

    # --- Chi-squared test

    sink(paste0(out.prefix, ".chisq.txt"), split=TRUE)

    cat(paste("Drug target percentage among branch length outlier human ev proteins:\n"))
    cat(paste(fast.dg.ct, "of", (fast.dg.ct + fast.notdg.ct)))
    cat("\n")
    cat(fast.dg.ct / (fast.dg.ct + fast.notdg.ct))
    cat("\n")

    cat(paste("Drug target percentage among other proteins:\n"))
    cat(paste(notfast.dg.ct, "of", (notfast.dg.ct + notfast.notdg.ct)))
    cat("\n")
    cat(notfast.dg.ct / (notfast.dg.ct + notfast.notdg.ct))
    cat("\n")

    cat(paste("Chi-squared test output:\n"))
    chi.res = chisq.test(c(fast.dg.ct, fast.notdg.ct, notfast.dg.ct, notfast.notdg.ct))
    chi.res.str = capture.output(chi.res)
    cat(chi.res.str)
    cat("\n")

    cat(paste("Chi-squared p-value:\n"))
    cat(chi.res$p.value)
    cat("\n")

    sink()

    # --- Hypergeometric test

    m = fast.dg.ct + notfast.dg.ct         # Genes drug target
    n = fast.notdg.ct + notfast.notdg.ct   # Genes NOT drug target
    k = fast.dg.ct + fast.notdg.ct         # Human accelerated genes
    x = 0:(fast.dg.ct + fast.notdg.ct)     # Genes both drug target and accelerated

    # Use the dhyper built-in function for hypergeometric density
    probabilities = dhyper(x, m, n, k, log = FALSE)

    sink(paste0(out.prefix, ".hypergeometric.txt"), split=TRUE)

    cat(paste("Hypergeometric test p-value:\n"))
    cat(sum(probabilities[fast.dg.ct:(fast.dg.ct + fast.notdg.ct)]))
    cat("\n")

    cat(paste("Fisher's test output:\n"))

    fish.res = fisher.test(matrix(c(fast.dg.ct,    fast.notdg.ct,
                                    notfast.dg.ct, notfast.notdg.ct),
                                  nrow=2))
    fish.res.str = capture.output(fish.res)
    cat(fish.res.str)
    cat("\n")

    cat(paste("Fisher's test p-value:\n"))
    cat(fish.res$p.value)
    cat("\n")

    sink()

    return(NA)
}

# ----------------------------------------------------------------------------------------
# --- All heated vs. control
# ----------------------------------------------------------------------------------------

# --- Of the accelerated human ev proteins, are likely drug targets over-represented?

is.brlen.outlier = bl.dg$brlen.ratio >  bl.high.cutoff
is.drug.target   = bl.dg$Total.score >= dg.high.cutoff

do.overrep.tests (bl.dg, is.brlen.outlier, is.drug.target,
                  "results/human_accel_drug_targets")

# --- Of the accelerated human ev proteins, are known drug targets over-represented?

is.brlen.outlier = bl.dg$brlen.ratio > bl.high.cutoff
is.drug.target   = bl.dg$Drugbank.approved != ""

do.overrep.tests (bl.dg, is.brlen.outlier, is.drug.target,
                  "results/human_accel_drug_targets_known")

# --- INVERSE TEST: Of the slowest evolving human ev proteins,
# ---               are predicted drug targets over-represented?

is.brlen.outlier.low = bl.dg$brlen.ratio < bl.low.cutoff
is.drug.target   = bl.dg$Total.score >= dg.high.cutoff

do.overrep.tests (bl.dg, is.brlen.outlier.low, is.drug.target,
                  "results/human_slowed_drug_targets")

# --- INVERSE TEST: Of the slowest evolving human ev proteins,
# ---               are known drug targets over-represented?

is.brlen.outlier.low = bl.dg$brlen.ratio < bl.low.cutoff
is.drug.target   = bl.dg$Drugbank.approved != ""

do.overrep.tests (bl.dg, is.brlen.outlier.low, is.drug.target,
                  "results/human_slowed_drug_targets_known")

# --- WEIRD INVERSE TEST: Of the accelerated human ev proteins,
# ---                     are definitely-not-predicted drug targets over-represented?

is.brlen.outlier  = bl.dg$brlen.ratio >  bl.high.cutoff
is.NOTdrug.target = bl.dg$Total.score <= dg.low.cutoff

do.overrep.tests (bl.dg, is.brlen.outlier, is.NOTdrug.target,
                  "results/human_accel_NOTdrug_targets")

# --- WEIRD INVERSE TEST: Of the slowest evolving human ev proteins,
# ---                     are definitely-not-predicted drug targets over-represented?

is.brlen.outlier.low = bl.dg$brlen.ratio < bl.low.cutoff
is.NOTdrug.target = bl.dg$Total.score <= dg.low.cutoff

do.overrep.tests (bl.dg, is.brlen.outlier.low, is.NOTdrug.target,
                  "results/human_slowed_NOTdrug_targets")

# ----------------------------------------------------------------------------------------
# --- Make list of interesting genes for GO overrep
# ----------------------------------------------------------------------------------------

bl.high.cutoff = quantile(bl$brlen.ratio, 0.95, na.rm=TRUE)
dg.high.cutoff = quantile(dg$Total.score, 0.95, na.rm=TRUE)

bl.low.cutoff  = quantile(bl$brlen.ratio, 0.05, na.rm=TRUE)
dg.low.cutoff  = quantile(dg$Total.score, 0.05, na.rm=TRUE)


dbl.outliers.high = bl.dg[bl.dg$brlen.ratio >  bl.high.cutoff &
                         (bl.dg$Total.score >= dg.high.cutoff |
                          bl.dg$Drugbank.approved != ""),]

dbl.outliers.low  = bl.dg[bl.dg$brlen.ratio <  bl.low.cutoff &
                         (bl.dg$Total.score >= dg.high.cutoff |
                          bl.dg$Drugbank.approved != ""),]

# Sort by branch length ratio
dbl.outliers.high = dbl.outliers.high[order(dbl.outliers.high$brlen.ratio,
    decreasing=TRUE),]
dbl.outliers.low = dbl.outliers.low  [order(dbl.outliers.low$brlen.ratio,
    decreasing=FALSE),]

cols.to.keep = c(1,7,8,6,2:4,10,9)
dbl.outliers.high = dbl.outliers.high[,cols.to.keep]
dbl.outliers.low  = dbl.outliers.low [,cols.to.keep]

names(dbl.outliers.high) = names(dbl.outliers.low) = c(
    "ID.Echinococcus_multilocularis", "ID.Schistosoma_mansoni", "ID.Taenia_solium",
    "Product",
    "Branch_length.human_infective", "Branch_length.non_human_infective",
    "Branch_length.ratio",
    "Tsai.drug_pred.Total_score", "Tsai.Drugbank_info")

write.table(dbl.outliers.high,
    file="reports/outlier_high_RHR_drug_targets.txt",
    sep="\t", row.names=FALSE, col.names=TRUE)
write.table(dbl.outliers.high[,1],
    file="reports/outlier_high_RHR_drug_targets.justEmuID.txt",
    sep="\t", row.names=FALSE, col.names=FALSE, quote=FALSE)

write.table(dbl.outliers.low,
    file="reports/outlier_low_RHR_drug_targets.txt",
    sep="\t", row.names=FALSE, col.names=TRUE)
write.table(dbl.outliers.low[,1],
    file="reports/outlier_low_RHR_drug_targets.justEmuID.txt",
    sep="\t", row.names=FALSE, col.names=FALSE, quote=FALSE)

# --- Also make BG list
write.table(unique(bl.dg[,1]), file="reports/BG.justEmuID.txt",
    sep="\t", row.names=FALSE, col.names=FALSE, quote=FALSE)

# ----------------------------------------------------------------------------------------
# --- Then, do GO overrepresentation analysis
# ----------------------------------------------------------------------------------------

# --- RHR+drug

# Old way of doing things before switching to R package:

#   # Use g:Profiler's tool for functional profiling, g:GOSt
#   # https://biit.cs.ut.ee/gprofiler/gost
#
#   # Using high RHR + predicted drug targets input:
#   #   outlier_high_RHR_drug_targets.justEmuID.txt
#   # and low RHR + predicted drug targets input:
#   #   outlier_low_RHR_drug_targets.justEmuID.txt
#   # and using list of all genes with branch length info as BG:
#   #   BG.justEmuID.txt
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
#   # reports/outlier_high_RHR_drug_targets.gProfiler_results.csv
#   # reports/outlier_low_RHR_drug_targets.gProfiler_results.csv

all.genes = unique(bl.dg[,1])

# high RHR + predicted/known drug targets genes
gprofiler.top_RHR_and_pred_drug = gost(dbl.outliers.high$ID.Echinococcus_multilocularis,
    organism = "ecmultprjeb122", ordered_query = FALSE,
    multi_query = FALSE, significant = FALSE, exclude_iea = FALSE,
    measure_underrepresentation = FALSE, evcodes = FALSE,
    user_threshold = 0.2, correction_method = "fdr",
    domain_scope = "annotated",
    custom_bg = all.genes,
    numeric_ns = "NO_NAMESPACE", sources = NULL, as_short_link = FALSE)

res.high = gprofiler.top_RHR_and_pred_drug$result[order(gprofiler.top_RHR_and_pred_drug$result$p_value),]

write.csv(res.high[,c(3,4,5,6,9,10,11)],
    file="reports/outlier_high_RHR_drug_targets.gProfiler_results.csv")

# low RHR + predicted/known drug targets genes
gprofiler.low_RHR_and_pred_drug = gost(dbl.outliers.low$ID.Echinococcus_multilocularis,
    organism = "ecmultprjeb122", ordered_query = FALSE,
    multi_query = FALSE, significant = FALSE, exclude_iea = FALSE,
    measure_underrepresentation = FALSE, evcodes = FALSE,
    user_threshold = 0.2, correction_method = "fdr",
    domain_scope = "annotated",
    custom_bg = all.genes,
    numeric_ns = "NO_NAMESPACE", sources = NULL, as_short_link = FALSE)

res.low = gprofiler.low_RHR_and_pred_drug$result[order(gprofiler.low_RHR_and_pred_drug$result$p_value),]

write.csv(res.low[,c(3,4,5,6,9,10,11)],
    file="reports/outlier_low_RHR_drug_targets.gProfiler_results.csv")

# --- Find terms significant in high but not significant in low

res.both = merge(res.high, res.low, by=c("term_id", "term_name"), suffixes=c(".high", ".low"))
res.only_high = res.both[res.both$p_value.low > 0.2 & res.both$p_value.high < 0.2,]

res.only_high = res.only_high[order(res.only_high$p_value.high),]

quit()
