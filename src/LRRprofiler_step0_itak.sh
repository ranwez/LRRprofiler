#!/bin/bash
#========================================================
# PROJET : LRRprofiler
# SCRIPT : LRRprofiler_step0_itak.sh
# AUTHOR : Celine Gottin
#========================================================
# DESCRIPTION :
# ARGUMENTS : o $1 : Proteome Path (fasta)
#             o $2 : NAME
# DEPENDENCIES : o iTAK v 1.7
#========================================================


#========================================================
#                Environment, variables, files
#========================================================
set -euo pipefail

# Variables
MAIN=$(pwd)
NAME=$2
RESDIR=$LRRPROFILER_RESDIR/Res_step0_itak

WD=$LRRPROFILER_TMP/wd_itak_${NAME}
mkdir $WD; cd $WD

# files
ln -s $1 .
PROTEOME=$(basename $1)

#========================================================
#                       SCRIPT
#========================================================

# running iTAK (-m p mean running itak only for kinase search, not TF)
perl $LG_ITAK -m p $PROTEOME

# Saving results
if [[ ! -e $RESDIR ]];then
    mkdir -p $RESDIR
fi

ITAK_KINASE_FILE="${PROTEOME}_output/shiu_alignment.txt"
OUTPUT_KINASE_FILE="$RESDIR/${NAME}_shiu_alignment.txt"

if [[ -e "$ITAK_KINASE_FILE" ]]; then
    cp "$ITAK_KINASE_FILE" "$OUTPUT_KINASE_FILE"
else
    echo "No protein kinase identified by iTAK."
    : > "$OUTPUT_KINASE_FILE"
fi

# Cleaning
cd $MAIN ; rm -r $WD

echo "END STEP 0"

#========================================================
#                    END of SCRIPT
#========================================================
