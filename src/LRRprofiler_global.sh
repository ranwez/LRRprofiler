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


# Variables
function quit_pb_option() {
    printf "\nOptions : --in_proteome <fastafile> ; --name <string>; --dev\n"
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



printf "\n\n============== PARSING OPTIONS\n"
while (( $# > 0 )); do
    echo "parse option $1"
    case "$1" in
	--in_proteome)
	    PROTEOME=$(readlink -f "$2"); shift 2;
	    if [[ ! -f $PROTEOME ]];then
		    echo "File $Proteome does not exist"
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

	*)
      echo "Option $1 is unknown please ckeck your command line"
      quit_pb_option
      ;;
    esac
done


# Variables
MAIN=$(pwd)
SCRIPT=${LG_SCRIPT}

TMP=tmp_$(date +'%Y%m%d_%H%M%S')
RESDIR=$MAIN/Res_$NAME

mkdir $TMP 

if [[ ! -e $RESDIR ]];then
    mkdir $RESDIR
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

if [[ ! -e $RESDIR/Res_step0_itak ]];then

  $SCRIPT/LRRprofiler_step0_itak.sh $PROTEOME $NAME $RESDIR/Res_step0_itak ${NAME}_shiu_alignment.txt
else
  echo -e "Found iTAK results; skipping step 0\n\n"
fi

## STEP 1 : Refining Profiles
##--------------------------------------------

echo "----------------------------"
echo "  STEP 1 : refining profiles"
echo -e "----------------------------\n"


# A. LRR-RLK

if [[ ! -e $RESDIR/Res_step1/LRR_kinase_$NAME.hmm ]]; then

    echo -e "Running step 1 for LRR-RLK motifs"

    gawk '{print($1,$3,$4)}' $RESDIR/Res_step0_itak/${NAME}_shiu_alignment.txt | sort -k1,1 -Vk2,3 | gawk 'BEGIN{prot="";start=0;end=0}{if($1==prot && $2<end){end=$3}else{if(prot!=""){print(prot,"Kinase",start,end)};prot=$1;start=$2;end=$3}}END{print(prot,"Kinase",start,end)}' > $MAIN/$TMP/ListeKinase.txt
    
    if [[ $devopt -eq 1 ]];then
        $SCRIPT/LRRprofiler_step1_AmelioProfil.sh --in_proteome $PROTEOME --list_proteins $MAIN/$TMP/ListeKinase.txt --in_profile ${LG_HMMlib}/SMART_LRR.hmm --out_dir $RESDIR/Res_step1 --out_profile_name LRR_kinase_$NAME.hmm --dev
    else
        $SCRIPT/LRRprofiler_step1_AmelioProfil.sh --in_proteome $PROTEOME --list_proteins $MAIN/$TMP/ListeKinase.txt --in_profile ${LG_HMM_lib}/SMART_LRR.hmm --out_dir $RESDIR/Res_step1 --out_profile_name LRR_kinase_$NAME.hmm
    fi

else
    echo -e "Found profile for LRR-RLK motifs; skipping step 1 for kinase protein\n"
fi

# B. NLR

if [[ ! -e $RESDIR/Res_step1/LRR_NLR_$NAME.hmm ]]; then

    echo -e "Running step 1 for LRR-NLR motifs"
    
    hmmsearch -o $TMP/del.tmp -Z 35000 --nobias --noali --domtblout $MAIN/$TMP/NBARC.tbl ${LG_HMMlib}/NB-ARC.hmm $PROTEOME

    gawk '$1!~/#/{print($1,$18,$19)}' $MAIN/$TMP/NBARC.tbl | sort -k1,1 -Vk2,3 | gawk 'BEGIN{prot="";start=0;end=0}{if($1==prot && $2<end){end=$3}else{if(prot!=""){print(prot,"NBARC",start,end)};prot=$1;start=$2;end=$3}}END{print(prot,"NBARC",start,end)}' > $MAIN/$TMP/ListeNBARC.txt

    if [[ $devopt -eq 1 ]];then
        $SCRIPT/LRRprofiler_step1_AmelioProfil.sh --in_proteome $PROTEOME --list_proteins $MAIN/$TMP/ListeNBARC.txt --in_profile ${LG_HMMlib}/SMART_LRR.hmm --out_dir $RESDIR/Res_step1 --out_profile_name LRR_NLR_$NAME.hmm --dev
    else
        $SCRIPT/LRRprofiler_step1_AmelioProfil.sh --in_proteome $PROTEOME --list_proteins $MAIN/$TMP/ListeNBARC.txt --in_profile ${LG_HMMlib}/SMART_LRR.hmm --out_dir $RESDIR/Res_step1 --out_profile_name LRR_NLR_$NAME.hmm
    fi

else
    echo -e "Found profile for LRR-NLR motifs; skipping step 1 for NB-ARC protein\n"
fi


## STEP 2 : Recherche LRR
##--------------------------------------------

echo -e "\n--------------------------------"
echo "  STEP 2 : Search for LRR motifs"
echo -e "--------------------------------\n"


if [[ ! -e $RESDIR/Res_step2 ]];then

    if [[ $devopt -eq 1 ]];then
        $SCRIPT/LRRprofiler_step2_RechercheLRR.sh --in_proteome $PROTEOME --name $NAME --rlk_profile $MAIN/Res_$NAME/Res_step1/LRR_kinase_$NAME.hmm --nlr_profile $MAIN/Res_$NAME/Res_step1/LRR_NLR_$NAME.hmm --out_dir $RESDIR/Res_step2 --dev
    else
        $SCRIPT/LRRprofiler_step2_RechercheLRR.sh --in_proteome $PROTEOME --name $NAME --rlk_profile $MAIN/Res_$NAME/Res_step1/LRR_kinase_$NAME.hmm --nlr_profile $MAIN/Res_$NAME/Res_step1/LRR_NLR_$NAME.hmm --out_dir $RESDIR/Res_step2
    fi
    
else
    echo -e "Found results for step 2; skipping step 2\n"
fi

## STEP 3 : Annotation sequences
##--------------------------------------------

echo -e "\n-----------------------------------------------"
echo "  STEP 3 : Sequence annotation & classification"
echo -e "-----------------------------------------------\n"

if [[ $devopt -eq 1 ]];then
    $SCRIPT/LRRprofiler_step3_classification.sh --in_proteome $PROTEOME --name $NAME --out_dir $RESDIR/Res_step3 --dev
else
    $SCRIPT/LRRprofiler_step3_classification.sh --in_proteome $PROTEOME --name $NAME --out_dir $RESDIR/Res_step3
fi



## cleaning
##-------------------------------

if [[ $devopt -eq 0 ]];then
    cd $MAIN ; rm -r $TMP
fi

#========================================================
#                      END
#========================================================





