#!/usr/bin/env python

from Bio.Nexus import Nexus
from os import listdir
from os.path import isfile, join
import re

file_list = [join("results/phylogeny/concat/", f)
    for f in listdir("results/phylogeny/concat/")
    if isfile(join("results/phylogeny/concat/", f))]

r = re.compile(".*\.nocomment\.nex")
file_list_emu = filter(r.match, file_list)

nexi = [(fname, Nexus.Nexus(fname)) for fname in file_list_emu]

r = re.compile(".*copy")

nexi_clean = [nex for nex in nexi if len(filter(r.match, nex[1].taxlabels)) == 0]

combined = Nexus.combine(nexi_clean)
combined.taxlabels
combined.write_nexus_data(filename=open('results/phylogeny/concat/combined.nex', 'w'))

combined.export_fasta(filename='results/phylogeny/concat/combined.fas')
