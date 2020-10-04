#!/usr/bin/env Rscript

# ----------------------------------------------------------------------------------------
# --- Plot WGS tree
# ----------------------------------------------------------------------------------------

library(treeio)
library(ape)
library(ggplot2)
library(ggtree)
library(plyr)

options(stringsAsFactors=FALSE)

tree = read.iqtree("results/phylogeny/concat/wgs_iq_tree/wgs_concord.cf.tree")

# This reads data goofily:
#  UFboot is actually:
#    sCF: Site concordance factor averaged over 100 quartets (=sCF_N/sN %)
#  and SH_aLRT is actually:
#    Label: Existing branch label
#  The label contains three values separated by slashes, like this:
#    100/69/43.4
#  These correspond to:
#    Label: Existing branch label (bootstrap values, for us)
#    gCF: Gene concordance factor (=gCF_N/gN %)
#    sCF: Site concordance factor averaged over 100 quartets (=sCF_N/sN %)

# --- Store stats

tree.stats = do.call(rbind, strsplit(tree@phylo$node.label, split="/"))
class(tree.stats) = "numeric"

tree@data$bootstrap = c(NA, tree.stats[,1])
tree@data$gCF       = c(NA, tree.stats[,2])
tree@data$sCF       = c(NA, tree.stats[,3])

tree@phylo = root.phylo(tree@phylo,
    which(tree@phylo$tip.label %in% c("T19", "T23", "T26")),
    resolve.root=TRUE)

# tree@phylo = root.phylo(tree@phylo,
#     which(tree@phylo$tip.label %in% c("T26")),
#     resolve.root=TRUE)

#root = rootnode(tree)
#root = 999  # Ignore for now since unrooted
root = which(tree@phylo$node.label == "Root")

# --- Set human infective Taenia to red

df.meta = data.frame(node=1:Nnode2(tree), category = 'Non-human infective Taenia')
hum.infect = which(tree@phylo$tip.label %in% c("T06", "T07", "T08", "T09"))
df.meta[hum.infect, 2] = "Human infective Taenia"

df.meta = df.meta[-c(which(tree@phylo$tip.label == "T07")),]

# --- Remove one of the T. solium samples

solium.node = getMRCA(tree@phylo, c("T06", "T07"))

tree@phylo = ape::drop.tip(tree@phylo, which(tree@phylo$tip.label == "T07"),
    trim.internal=FALSE)

# --- Rename tips with species

tree@phylo$tip.label = revalue(tree@phylo$tip.label,
    c("Emulti" = paste0('italic(', 'Echinococcus~multilocularis', ')'),
      "T06" = paste0('italic(', 'Taenia~solium', ')'),
      "T07" = paste0('italic(', 'Taenia~solium', ')'),
      "T08" = paste0('italic(', 'Taenia~saginata', ')'),
      "T09" = paste0('italic(', 'Taenia~asiatica', ')'),
      "T10" = paste0('italic(', 'Taenia~hydatigena', ')'),
      "T11" = paste0('italic(', 'Taenia~multiceps', ')'),
      "T13" = paste0('italic(', 'Taenia~pisiformis', ')'),
      "T16" = paste0('italic(', 'Taenia~crassiceps', ')'),
      "T17" = paste0('italic(', 'Hydatigera~taeniaeformis', ')'),
      "T19" = paste0('italic(', 'Cittotaenia', ')', "~sp."),
      "T23" = paste0('Unknown~sp.~', 'italic(', 'ex', ')', '~goat'),
      "T26" = paste0('Unknown~sp.~', 'italic(', 'ex', ')', '~chicken')))

# --- Plot tree

p = ggtree(tree, color="black", size=1, linetype=1,  right=TRUE) +
    ggplot2::xlim(-1.5, 6) +
    geom_tiplab(aes(color=factor(df.meta$category)),
                parse=TRUE, size=4.5, hjust = -0.060, fontface="bold",
                align=TRUE, linetype='dotted', linesize=.3) +
    geom_point2(aes(subset = !isTip & node != root,
                    fill = cut(gCF, c(0, 50, 75, 100))),
                    shape = 21, size = 3) +
    scale_fill_manual(values=c("black", "grey", "white"), guide='legend',
                    name='Gene concordance factor (gCF)',
                    breaks=c('(75,100]', '(50,75]', '(0,50]'),
                    labels=expression(gCF >= 75,50 <= gCF * " < 75", gCF < 50)) +
    scale_color_manual(values=c("black", "firebrick"), guide='legend',
                       name='',
                       breaks=c("Non-human infective Taenia", "Human infective Taenia")) +
    theme(legend.position = c(0.15, 0.25))

ggsave(p, file="reports/WGS_tree.pdf", height=5, width=7)
