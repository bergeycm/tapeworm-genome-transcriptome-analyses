#!/usr/bin/env Rscript

library(ape)

# ----------------------------------------------------------------------------------------
# --- Compute distance of human and non-human worms to outgroup (T19 - Cittotaenia)
# ----------------------------------------------------------------------------------------

args = commandArgs(trailingOnly=TRUE)

emu.id = args[1]    # E.g. EmuJ_000619800.1

# Ignoring: T07 - other T. solium
#           T11 - T. multiceps
#           T23 - Hymenolepis diminuta
#           T26 - Unknown worm
#           Emu - Echinococcus

tx.human    = c("T06", "T07", "T08", "T09")
tx.nonhuman = c("T10", "T13", "T16", "T17")

gt = read.nexus(file = paste0("results/phylogeny/Echin_orthologs_aTRAM/",
                              emu.id, "/", emu.id, ".flt.nex.con.tre"))

dist.to.out = cophenetic(gt)["T19",]

avg.dist.human    = mean(dist.to.out[tx.human],    na.rm=TRUE)
avg.dist.nonhuman = mean(dist.to.out[tx.nonhuman], na.rm=TRUE)

ratio.dist = avg.dist.human / avg.dist.nonhuman

dist.str = paste(paste(names(dist.to.out), dist.to.out, sep=":"), collapse=";")

write(paste(emu.id,
            avg.dist.human, avg.dist.nonhuman,
            ratio.dist, dist.str, sep="\t"), stdout())

quit()
