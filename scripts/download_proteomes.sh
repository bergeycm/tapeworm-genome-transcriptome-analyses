#!/bin/bash

# ========================================================================================
# --- Download reference proteomes
# ========================================================================================

taxon=$1

if [[ "$taxon" == "echinococcus_multilocularis" ]]; then

    mkdir -p genomes/echinococcus_multilocularis

    FASTA=genomes/echinococcus_multilocularis/echinococcus_multilocularis.PRJEB122.WBPS12.protein.fa

    URL=ftp://ftp.ebi.ac.uk/pub/databases/wormbase/parasite/releases/WBPS12/species
    URL=$URL/echinococcus_multilocularis/PRJEB122/echinococcus_multilocularis.PRJEB122.WBPS12.protein.fa.gz

    wget -O $FASTA.gz \
        $URL

    gunzip -f $FASTA.gz

    # Remove line breaks...
    # ...and stuff in header after space
    awk '!/^>/ { printf "%s", $0; n = "\n" }
         /^>/  { print n $0; n = "" }
         END   { printf "%s", n }' $FASTA | \
        sed -e "/>/ s/ .*//" > ${FASTA/.fa/.oneline.fa}

    # Also download DNA seqs for all genes

    FASTA=genomes/echinococcus_multilocularis/echinococcus_multilocularis.PRJEB122.WBPS12.CDS_transcripts.fa

    URL=ftp://ftp.ebi.ac.uk/pub/databases/wormbase/parasite/releases/WBPS12/species
    URL=$URL/echinococcus_multilocularis/PRJEB122/echinococcus_multilocularis.PRJEB122.WBPS12.CDS_transcripts.fa.gz

    wget -O $FASTA.gz \
        $URL

    gunzip -f $FASTA.gz

fi

exit
