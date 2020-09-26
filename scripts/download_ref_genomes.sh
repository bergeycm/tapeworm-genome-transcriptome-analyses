#!/bin/sh

# ----------------------------------------------------------------------------------------
# --- Download reference genomes
# ----------------------------------------------------------------------------------------

mkdir -p genomes
cd genomes

# --- Download Taenia solium genome, Tsolium_Mexico_v1

# http://parasite.wormbase.org/Taenia_solium_prjna170813/Info/Index/

mkdir Tsolium_Mexico_v1
cd Tsolium_Mexico_v1

GENOME_FA=taenia_solium.PRJNA170813.WBPS8.genomic_softmasked.fa

TSOLIUM_URL=ftp://ftp.ebi.ac.uk/pub/databases/wormbase/parasite/releases/WBPS8/species
TSOLIUM_URL=${TSOLIUM_URL}/taenia_solium/PRJNA170813
TSOLIUM_URL=${TSOLIUM_URL}/$GENOME_FA.gz

wget $TSOLIUM_URL \
    -O ${GENOME_FA}.gz
gunzip ${GENOME_FA}.gz

mv ${GENOME_FA} Tsolium_Mexico_v1.fa

cd ..

# --- Download genomes from Taenia Genome Database (TGD)

# http://taenia.big.ac.cn/

# T. solium is not actually tarred, despite its name
wget http://taenia.big.ac.cn/data/T.solium.ch.genome.fasta.tar.gz \
    -O TGD_Tsolium.fa.gz
wget http://taenia.big.ac.cn/data/TSA/Tsa.v1.genome.fa.gz \
    -O TGD_Tsaginata.fa.gz
wget http://taenia.big.ac.cn/data/TAS/Tas.v1.genome.fa.gz \
    -O TGD_Tasiatica.fa.gz

for TGD_ROOT in TGD_Tsolium TGD_Tsaginata TGD_Tasiatica; do

    mkdir $TGD_ROOT
    cd $TGD_ROOT

    mv ../$TGD_ROOT.fa.gz .

    gunzip $TGD_ROOT.fa.gz

    cd ..
done

cd ..

exit
