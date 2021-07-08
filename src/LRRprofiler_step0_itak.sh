#!/bin/bash
#========================================================
# PROJET : LRRannotation
# SCRIPT : LRRannotation_step0_itak.sh
# AUTHOR : Celine Gottin
# CREATION : 2020.01.08
#========================================================
# DESCRIPTION : 
# ARGUMENTS : o $1 : Proteome Path (fasta)
#             o $2 : NAME
# DEPENDENCIES : o iTAK v 1.7               
#========================================================


#========================================================
#                Environment, variables, files
#========================================================


# Variables
MAIN=$(pwd)
NAME=$2
RESDIR=$LRRPROFILER_RESDIR/Res_step0_itak
WD=$MAIN/wd_itak_${NAME}_$(date +'%H%M%S')

# working dir
if [[ -e $WD ]];then
    rm -r $WD
fi

mkdir $WD; cd $WD

# files
ln -s $1 .
PROTEOME=$(basename $1)

#========================================================
#                       SCRIPT
#========================================================

# running iTAK (-m p mean running itak only for kinase search, not TF)
#perl /usr/local/bioinfo/iTAK/1.7/iTAK.pl -m p $PROTEOME
perl $LG_ITAK -m p $PROTEOME

# Saving results
if [[ ! -e $RESDIR ]];then
    mkdir -p $RESDIR
fi

cp ${PROTEOME}_output/shiu_alignment.txt $RESDIR/${NAME}_shiu_alignment.txt

# Cleaning
cd $MAIN ; rm -r $WD

echo "END STEP 0"

#========================================================
#                    END of SCRIPT
#========================================================
