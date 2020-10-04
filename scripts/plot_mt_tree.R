#!/usr/bin/env Rscript

# ----------------------------------------------------------------------------------------
# --- Plot mtDNA tree
# ----------------------------------------------------------------------------------------

library(treeio)
library(ape)
library(ggplot2)
library(ggtree)
library(plyr)

options(stringsAsFactors=FALSE)

tree = read.newick("results/phylogeny/mt_phylogeny/mt_genes_concat.contree")

#  The label contains bootstrap values

tree = root.phylo(tree,
    which(tree$tip.label == "NC_002544-Schistosoma-japonicum"),
    resolve.root=TRUE)

root = which(tree$node.label == "Root")

# --- Set our samples to blue

df.meta = data.frame(node=1:Nnode2(tree), category = 'Published')
new.data = which(grepl("mt_", tree$tip.label))
df.meta[new.data, 2] = "Newly generated"

# --- Rename tips with species

pub.samps = data.frame(do.call(rbind,
    strsplit(tree$tip.label[!grepl("mt_", tree$tip.label)], split="-")))

pub.samps.str = do.call(rbind, lapply(1:nrow(pub.samps), function (x) {
    paste0('italic(', pub.samps[x,2], "~", pub.samps[x,3], ')')
}))

tree$tip.label[!grepl("mt_", tree$tip.label)] = pub.samps.str

tree$tip.label = revalue(tree$tip.label,
    c("mt_T06" = paste0('italic(', 'Taenia~solium', ')'),
      "mt_T07" = paste0('italic(', 'Taenia~solium', ')'),
      "mt_T08" = paste0('italic(', 'Taenia~saginata', ')'),
      "mt_T09" = paste0('italic(', 'Taenia~asiatica', ')'),
      "mt_T10" = paste0('italic(', 'Taenia~hydatigena', ')'),
      "mt_T11" = paste0('italic(', 'Taenia~multiceps', ')'),
      "mt_T13" = paste0('italic(', 'Taenia~pisiformis', ')'),
      "mt_T16" = paste0('italic(', 'Taenia~crassiceps', ')'),
      "mt_T17" = paste0('italic(', 'Hydatigera~taeniaeformis', ')'),
      "mt_T19" = paste0('italic(', 'Cittotaenia', ')', "~sp."),
      "mt_T23" = paste0('Unknown~sp.~', 'italic(', 'ex', ')', '~goat'),
      "mt_T26" = paste0('Unknown~sp.~', 'italic(', 'ex', ')', '~chicken')))

# --- Figure out which clades to label

node.diphyllobothriidae = getMRCA(tree, c(
    tree$tip.label[grepl("nihonkaiense", tree$tip.label)],
    tree$tip.label[grepl("decipiens", tree$tip.label)]))

node.anoplocephalidae = getMRCA(tree, c(
    tree$tip.label[grepl("expansa", tree$tip.label)],
    tree$tip.label[grepl("magna", tree$tip.label)]))

node.hymenolepididae = getMRCA(tree, c(
    tree$tip.label[grepl("lanceolata", tree$tip.label)],
    tree$tip.label[grepl("crawfordi", tree$tip.label)]))

node.dipylidiidae = grep("caninum", tree$tip.label)

node.taeniidae = getMRCA(tree, c(
    tree$tip.label[grepl("oligarthrus", tree$tip.label)],
    tree$tip.label[grepl("saginata", tree$tip.label)]))

node.versteria = grep("mustelae", tree$tip.label)

node.echinococcus = getMRCA(tree, c(
    tree$tip.label[grepl("oligarthrus", tree$tip.label)],
    tree$tip.label[grepl("canadensis", tree$tip.label)]))

node.hydatigera = getMRCA(tree, c(
    tree$tip.label[grepl("parva", tree$tip.label)],
    tree$tip.label[grepl("taeniaeformis", tree$tip.label)]))

node.taenia = getMRCA(tree, c(
    tree$tip.label[grepl("crassiceps", tree$tip.label)],
    tree$tip.label[grepl("saginata", tree$tip.label)]))

# --- Plot tree

boot.cats = cut(c(rep(NA, length(tree$tip.label)), as.numeric(tree$node.label)),
    c(0, 90, 99, 100))

p = ggtree(tree, color="black", linetype=1, right=TRUE) +
    geom_point2(aes(subset = !isTip & node != root, fill = boot.cats),
        shape = 21, size = 2, stroke=0) +
    geom_cladelabel(node=node.diphyllobothriidae, label='bold(Diphyllobothriidae)',
        align=TRUE, offset=1, offset.text = 0.05, extend=-0.2, parse=TRUE) +
    geom_cladelabel(node=node.anoplocephalidae, label='bold(Anoplocephalidae)',
        align=TRUE, offset=1, offset.text = 0.05, extend=-0.2, parse=TRUE) +
    geom_cladelabel(node=node.hymenolepididae, label='bold(Hymenolepididae)',
        align=TRUE, offset=1, offset.text = 0.05, extend=-0.2, parse=TRUE) +
    geom_cladelabel(node=node.dipylidiidae, label="bold(Dipylidiidae)",
        align=TRUE, offset=1, offset.text = 0.05, extend=-0.2, parse=TRUE) +
    geom_cladelabel(node=node.taeniidae, label="bold(Taeniidae)",
        align=TRUE, offset=1, offset.text = 0.05, extend=-0.2, parse=TRUE) +
    geom_cladelabel(node=node.versteria, label='italic(Versteria)',
        align=TRUE, offset=-0.5, offset.text = 0.05, extend=-0.2, parse=TRUE) +
    geom_cladelabel(node=node.echinococcus, label='italic(Echinococcus)',
        align=TRUE, offset=-0.5, offset.text = 0.05, extend=-0.2, parse=TRUE) +
    geom_cladelabel(node=node.hydatigera, label='italic(Hydatigera)',
        align=TRUE, offset=-0.5, offset.text = 0.05, extend=-0.2, parse=TRUE) +
    geom_cladelabel(node=node.taenia, label='italic(Taenia)',
        align=TRUE, offset=-0.5, offset.text = 0.05, extend=-0.2, parse=TRUE) +
    geom_treescale(x=4, y=5, width=0.5, color='black') +
    ggplot2::xlim(0, 5) +
    geom_tiplab(aes(color=factor(df.meta$category)),
                parse=TRUE, size=3, hjust = -0.060, fontface="bold",
                align=FALSE) +
    scale_color_manual(values=c("black", "deepskyblue4"), guide='legend',
                       name='',
                       breaks=c("Published", "Newly generated")) +
    scale_fill_manual(values=c(NA, "#EECF22", "#EA530B"), guide='legend',
                    name='Bootstrap support',
                    breaks=c('(99,100]', '(90,99]', '(0,90]'),
                    labels=c("100", "90-99", "0-90")) +
    theme(legend.position="None",
          plot.margin = unit(c(0,0,0,0), "inch"))

ggsave(p, file="reports/mt_tree.pdf", height=10, width=10)
