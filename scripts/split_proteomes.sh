#!/bin/bash

# ----------------------------------------------------------------------------------------
# --- Split proteomes into hunks for aTRAM
# ----------------------------------------------------------------------------------------

PROTEOME=$1

split -l 100 \
    --suffix-length=3 \
    -d $PROTEOME \
    ${PROTEOME}_PART

touch ${PROTEOME}-split_marker

exit
