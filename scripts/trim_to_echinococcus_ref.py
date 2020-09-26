#!/usr/bin/env python

# ----------------------------------------------------------------------------------------
# --- Remove gaps in sequence that don't match Echinococcus reference
# ----------------------------------------------------------------------------------------

from Bio import AlignIO
import sys

# e.g. "results/phylogeny/Echin_orthologs_aTRAM/EmuJ_000001200.1/EmuJ_000001200.1.aln.fa"
in_aln = sys.argv[1]

aln = AlignIO.read(in_aln, "fasta")

for record in aln:
    if record.id == "Emulti":
        to_cut = []
        for pos, char in enumerate(record.seq):
            if char == "-":
                to_cut.append(pos)

n = int(len(aln[0]))
i = n

while i >= 0:
    if i in to_cut:
        if i == 0:
            aln = aln[:, 1:]
        elif i+1 == n:
            aln = aln[:, :i]
        else:
            aln = aln[:, :i] + aln[:, i+1:]
        n -= 1
    i -= 1

print (aln)

# --- Write trimmed alignment

AlignIO.write(aln, in_aln.replace("aln.fa", "flt.fa"), "fasta")
