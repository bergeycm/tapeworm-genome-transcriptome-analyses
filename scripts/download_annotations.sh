#!/bin/sh

# ----------------------------------------------------------------------------------------
# --- Download reference genomes
# ----------------------------------------------------------------------------------------

mkdir -p genomes/annotations
cd genomes/annotations

# --- Download Taenia solium annotations, Tsolium_Mexico_v1

# http://parasite.wormbase.org/Taenia_solium_prjna170813/Info/Index/
GFF=taenia_solium.PRJNA170813.WBPS9.annotations.gff3

TSOLIUM_URL=ftp://ftp.ebi.ac.uk/pub/databases/wormbase/parasite/releases/WBPS9/species
TSOLIUM_URL=${TSOLIUM_URL}/taenia_solium/PRJNA170813
TSOLIUM_URL=${TSOLIUM_URL}/$GFF.gz

wget $TSOLIUM_URL \
    -O $GFF.gz
gunzip -c $GFF.gz > $GFF

# --- Download TGD project annotations (Only asiatica and saginata, solium not available)

# --- T. asiatica

GFF=Tas.v1.gff3

TASIATICA_URL=http://taenia.big.ac.cn/data/TAS/Tas.v1.gff3.gz

wget $TASIATICA_URL \
    -O $GFF.gz
gunzip -c $GFF.gz > $GFF

# --- T. saginata

GFF=Tsa.v1.gff3

TSAGINATA_URL=http://taenia.big.ac.cn/data/TSA/Tsa.v1.gff3.gz

wget $TSAGINATA_URL \
    -O $GFF.gz
gunzip -c $GFF.gz > $GFF

cd ../..

exit
