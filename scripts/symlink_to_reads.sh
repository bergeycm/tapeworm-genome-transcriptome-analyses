#!/bin/bash

# ----------------------------------------------------------------------------------------
# --- Link to data in shared directory
# ----------------------------------------------------------------------------------------

# --- WGS data

mkdir -p data/WGS/

for FA in $GRP/tapeworm_WGS/tapeworms_hiseq_2017-02/GP01272017/*.fastq.gz; do
    NEW_NAME=`basename $FA | sed -e "s/_S[0-9]\+_L[0-9]\+_R\([0-9]\)_[0-9]\+/.R\1/"`
    ln -s $FA data/WGS/$NEW_NAME
done
