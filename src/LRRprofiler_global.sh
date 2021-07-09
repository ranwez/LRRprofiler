#!/bin/bash
#========================================================
# PROJET : LRRprofiler
# SCRIPT : LRRprofiler_global.sh
# AUTHOR : Celine Gottin
#========================================================
# DESCRIPTION : 
# ARGUMENTS : o --in_proteome : Path of proteome file in fasta format
#             o --name : project name for outputfile identification
#             o --dev : if given, run in dev mode (conserve all temporary files for debug)
#             o --nobuild : if provided, the pipeline will not rebuild specific LRR HMM profils
#                           and use HMM from the provided HMM_lib exclusively.
# DEPENDENCIES : o iTAK v 1.7
#                o HMMER v 3.1b2
#                o MAFFT v7.271
#                o TMHMM 2.0a
#========================================================


#========================================================
#                Environment & variables
#========================================================

# on error exit flag : set -e
set -o errexit

# error if a var is unset : set -u
set -o nounset

# raise error in pipe
set -o pipefail


function quit_pb_option() {
    printf "\nOptions : --in_proteome <fastafile> ; --name <string>; --dev ; --nobuild\n"
    exit 1
}

if (( $# == 0)); then
  quit_pb_option
fi


echo "======================================================================================"
echo "||                                      LRRprofiler                                 ||"
echo "======================================================================================"

PROTEOME=""
NAME=""
devopt=0
nobuild=0



printf "\n\n============== PARSING OPTIONS\n"
while (( $# > 0 )); do
    echo "parse option $1"
    case "$1" in
    --in_proteome)
        PROTEOME=$(readlink -f "$2"); shift 2;
        if [[ ! -f $PROTEOME ]];then
            echo "File $PROTEOME does not exist"
            quit_pb_option
        fi
        echo "    proteome file is $PROTEOME"
        ;;
    --name)
        NAME="$2" ; shift 2
        echo "    job name is $NAME"
        ;;
    --dev)
      ## development option; will conserve running folder with every temporary files
        devopt=1; shift 1
        echo "    dev option set in"
        ;;
    --nobuild)
      ## will skip HMM refining step (step1)
        nobuild=1; shift 1
        echo "    nobuild option set in"
        ;;
    *)
      echo "Option $1 is unknown please ckeck your command line"
      quit_pb_option
      ;;
    esac
done


# Variables
MAIN=$(pwd)
SCRIPT=${LG_SCRIPT}

LRRPROFILER_TMP=$MAIN/tmp_LRRprofiler_$(date +'%Y%m%d_%H%M%S')
export LRRPROFILER_TMP
mkdir $LRRPROFILER_TMP 

LRRPROFILER_RESDIR=$MAIN/Res_$NAME
export LRRPROFILER_RESDIR
if [[ ! -e $LRRPROFILER_RESDIR ]];then
    mkdir $LRRPROFILER_RESDIR
fi

#========================================================
#                       SCRIPT
#========================================================

echo -e "\n\n"

## STEP 0 : iTAK
##--------------------------------------------

echo "-----------------------"
echo "  STEP 0 : running iTAK"
echo -e "-----------------------\n"

if [[ ! -e $LRRPROFILER_RESDIR/Res_step0_itak ]];then
  $SCRIPT/LRRprofiler_step0_itak.sh $PROTEOME $NAME
else
  echo -e "Found iTAK results; skipping step 0\n"
fi

## STEP 1 : Refining Profiles
##--------------------------------------------

echo -e "\n----------------------------"
echo "  STEP 1 : refining profiles"
echo -e "----------------------------\n"

## list Kinase and NB-ARC file
gawk '{print($1,$3,$4)}' $LRRPROFILER_RESDIR/Res_step0_itak/${NAME}_shiu_alignment.txt | sort -k1,1 -Vk2,3 | gawk 'BEGIN{prot="";start=0;end=0}{if($1==prot && $2<end){end=$3}else{if(prot!=""){print(prot,"Kinase",start,end)};prot=$1;start=$2;end=$3}}END{print(prot,"Kinase",start,end)}' > $LRRPROFILER_TMP/ListeKinase.txt

hmmsearch -o /dev/null -Z 35000 --nobias --noali --domtblout $LRRPROFILER_TMP/NBARC.tbl ${LG_HMMlib}/NB-ARC.hmm $PROTEOME
gawk '$1!~/#/{print($1,$18,$19)}' $LRRPROFILER_TMP/NBARC.tbl | sort -k1,1 -Vk2,3 | gawk 'BEGIN{prot="";start=0;end=0}{if($1==prot && $2<end){end=$3}else{if(prot!=""){print(prot,"NBARC",start,end)};prot=$1;start=$2;end=$3}}END{print(prot,"NBARC",start,end)}' > $LRRPROFILER_TMP/ListeNBARC.txt


if [[ $nobuild -eq 1 ]];then
    echo -e "nobuild option; skipping step 1\n"
else
    # A. LRR-RLK
    if [[ ! -e $LRRPROFILER_RESDIR/Res_step1/LRR_kinase_$NAME.hmm ]]; then
        echo -e "Running step 1 for LRR-RLK motifs"
        gawk '{if(NR==FNR){P[">"$1]=1}else{if($1~/>/){if(P[$1]==1){pr=1}else{pr=0}};if(pr==1){print}}}' $LRRPROFILER_TMP/ListeKinase.txt $PROTEOME > $LRRPROFILER_TMP/subproteome_kinase.fasta
        $SCRIPT/LRRprofiler_step1_AmelioProfil.sh --in_proteins $LRRPROFILER_TMP/subproteome_kinase.fasta --in_profile $LG_HMMlib/SMART_LRR.hmm --out_profile_name LRR_kinase_${NAME}.hmm
        cp $LRRPROFILER_TMP/ListeKinase.txt $LRRPROFILER_RESDIR/Res_step1/.
        echo "moving list ListeKinase.txt to $LRRPROFILER_RESDIR/Res_step1/ "
    else
        echo -e "Found profile for LRR-RLK motifs; skipping step 1 for kinase protein\n"
    fi

    # B. NLR
    if [[ ! -e $LRRPROFILER_RESDIR/Res_step1/LRR_NLR_${NAME}.hmm ]]; then
        echo -e "Running step 1 for LRR-NLR motifs"
        gawk '{if(NR==FNR){P[">"$1]=1}else{if($1~/>/){if(P[$1]==1){pr=1}else{pr=0}};if(pr==1){print}}}' $LRRPROFILER_TMP/ListeNBARC.txt $PROTEOME > $LRRPROFILER_TMP/subproteome_NBARC.fasta
        $SCRIPT/LRRprofiler_step1_AmelioProfil.sh --in_proteins $LRRPROFILER_TMP/subproteome_NBARC.fasta --in_profile ${LG_HMMlib}/SMART_LRR.hmm --out_profile_name LRR_NLR_$NAME.hmm
        cp $LRRPROFILER_TMP/ListeNBARC.txt
        echo "moving list ListeNBARC.txt to $LRRPROFILER_RESDIR/Res_step1/ "
    else
        echo -e "Found profile for LRR-NLR motifs; skipping step 1 for NB-ARC protein\n"
    fi
fi

## STEP 2 : LRR detection
##--------------------------------------------

echo -e "\n--------------------------------"
echo "  STEP 2 : Search for LRR motifs"
echo -e "--------------------------------\n"


if [[ ! -e $LRRPROFILER_RESDIR/Res_step2 ]];then
    $SCRIPT/LRRprofiler_step2_RechercheLRR.sh --in_proteome $PROTEOME --name $NAME --rlk_profile $LRRPROFILER_RESDIR/Res_step1/LRR_kinase_$NAME.hmm --nlr_profile $LRRPROFILER_RESDIR/Res_step1/LRR_NLR_$NAME.hmm --out_dir $LRRPROFILER_RESDIR/Res_step2
else
    echo -e "Found results for step 2; skipping step 2\n"
fi

## STEP 3 : Sequence annotation
##--------------------------------------------

echo -e "\n-----------------------------------------------"
echo "  STEP 3 : Sequence annotation & classification"
echo -e "-----------------------------------------------\n"

$SCRIPT/LRRprofiler_step3_classification.sh --in_proteome $PROTEOME --name $NAME --out_dir $LRRPROFILER_RESDIR/Res_step3


## cleaning
##-------------------------------

if [[ $devopt -eq 0 ]];then
    cd $MAIN ; rm -r $LRRPROFILER_TMP
fi

#========================================================
#                      END
#========================================================





