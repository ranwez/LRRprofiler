#!/bin/bash
#========================================================
# PROJET : LRRprofiler
# SCRIPT : LRRprofiler_step1_AmelioProfil.sh
# AUTHOR : Celine Gottin
#========================================================
# DESCRIPTION : 
# ARGUMENTS : o --in_proteins : Proteome subset (fasta)
#             o --in_profile : Initial Profile Path
#             o --out_profile_name : Final Profile name
# DEPENDENCIES : o HMMER v 3.1b2
#                o MAFFT v 7.313
#========================================================


#========================================================
#             Environment, variables, files
#========================================================


# on error exit flag : set -e
set -o errexit

# error if a var is unset : set -u
set -o nounset

# raise error in pipe
set -o pipefail


# Variables
function quit_pb_option() {
    printf "\nOptions : --in_proteins ; --in_profile ; --out_profile_name\n"
    exit 1
}

if (( $# == 0)); then
  quit_pb_option
fi

Proteins=""
initProfile=""
OUT_DIR="RES_step1_AmelioProfile_$(date +'%Y%m%d_%H%M%S')"
outname="new_LRR_profile.hmm"
itMax=15

while (( $# > 0 )); do
    case "$1" in
    --in_proteins)
        Proteins=$(readlink -f "$2"); shift 2;
        if [[ ! -f $Proteins ]];then
            echo "File $Proteins does not exist"
            quit_pb_option
        fi
        ;;
    --in_profile)
        initProfile=$(readlink -f "$2"); shift 2;
        if [[ ! -f $initProfile ]];then
            echo "File $initProfile does not exist"
            quit_pb_option
        fi
        ;;
    --out_profile_name)
        outname="$2"; shift 2
        ;;
    *)
      echo "Option $1 is unknown please ckeck your command line"
      quit_pb_option
      ;;
    esac
done

if [[ -e $LRRPROFILER_RESDIR ]];then
    OUT_DIR=$LRRPROFILER_RESDIR/Res_step1
fi

MAIN=$(pwd)

# Working dir
WD=$LRRPROFILER_TMP/wd_amelioProfil_${outname%.hmm}

mkdir $WD ; cd $WD


#========================================================
#                       Script
#========================================================

echo "============================================="
echo "Proteins : $(basename $Proteins)"
echo "Initial Profile : $(basename $initProfile)"
echo "Max number of iteration : $itMax "
echo "============================================="

#--------------------------------------------------------------------
# 1. Extract protein sample from protein List

nbPf=$(grep -c ">" $Proteins)

echo -e "\nFound $nbPf proteins in $Proteins file."

gawk 'BEGIN{OFS=";";protName=""}{if($1~/^>/){if(length(protName)>1){print(protName,sequence)};gsub(">","");protName=$1;sequence=""}else{sequence=sequence""$1}}END{print(protName,sequence)}' $Proteins > proteome_tmp.csv

#--------------------------------------------------------------------
# 2. Iterative search of LRR repeat

profile=$initProfile ;
size=$(gawk 'NR==3{print($2)}' $initProfile) ; ## size of initial profile
nbMotifs=0 ;
nbMotifsOld=0 ;
nbProt=0 ;
nbProtOld=0 ;
nbChar=0 ;
nbCharOld=0 ;

penalty=0 ;

echo ""
echo "---------------------------"

#res file header
echo "profile\tnbProt\tnb_5%\tnb_1%\tnb_0.1%\tnb_motif\tmean_size\tsum_length" > RES_construction_profile.txt


for it in `seq 1 $itMax`
do   
    echo ""
    echo "iteration $it ..."


# a) execute hmmsearch the output file contain comment lines starting with # that should be ignored
    
    hmmsearch -E 1000 --domE 1000 -o /dev/null --nobias --noali --domtblout HMMres${it}.tmp ${profile} $Proteins 

    nbProtOld=$nbProt
    nbProt=$(grep -v "^#"  HMMres${it}.tmp | cut -f1 -d" "| sort -u | wc -l)



# b) Profile Results : Profile; nbProteine from list found; threshold 0.5,0.1,0.01 ; nbMotifs; mean length of motifs

    gawk -v prof=$(basename ${profile}) 'BEGIN{OFS="\t"; nb=0; sizeLRR=0; nbProt=0; nbProt5=0; nbProt1=0; nbProt01=0;}
            $1!~/^#/{P[$1]++;
                    if($7<=0.001){nb++;sizeLRR=sizeLRR+($19-$18+1)};
                    if(P[$1]==1){ nbProt++;
                      if($7<=0.05){nbProt5++};
                      if($7<=0.01){ nbProt1++};
                      if($7<=0.001){ nbProt01++;}
                    }}END{print(prof,nbProt,nbProt5,nbProt1,nbProt01,nb,sizeLRR/nb,sizeLRR)}' HMMres${it}.tmp >> RES_construction_profile.txt




# c) identify motifs with length consistant with the expected familly pattern length 
    
    gawk -v Psize=$size 'BEGIN{OFS=";"} 
        $1!~/^#/{alignMatchLg=$(19)-$(18)+1;envMatchLg=$(21)-$(20)+1; 
          if(alignMatchLg>(0.8*Psize) && envMatchLg<(1.1*Psize))
            {print($1,$10,$11,$18,$19,$20,$21)}}' HMMres${it}.tmp > motifs.tmp
    
    nbMotifsOld=$nbMotifs ;
    nbMotifs=$(grep -c . motifs.tmp) ;

 # d) do we have a better HMM or not ? If not continue a few more steps to avoid plateus then stop

    nbCharOld=$nbChar ;
    nbChar=$(tail -1 RES_construction_profile.txt | cut -f8) ;

    if (( nbChar <= nbCharOld ));then
        let penalty=${penalty}+1;
    fi

    if (( $penalty >= 3 ));then
        break ;
    fi


# e) extract LRR motifs from the protein sequence

    gawk -F";" 'BEGIN{OFS=";"}{
                    if(NR==FNR){
                        SEQ[$1]=$2}
                    else{
                        len=$7-$6+1;
                        L[$1]=len;
                        seq=substr(SEQ[$1],$6,L[$1]);
                        print(">"$1"."$2);
                        print(seq)}
                }' proteome_tmp.csv motifs.tmp > motifsIt${it}.fasta


# f) Align sequence using MAFFT and build new profile
    
    echo "Building new HMM ..."
    mafft --quiet motifsIt${it}.fasta > motifsIt${it}_align.fasta

    hmmbuild -o /dev/null --amino Profile_${it}.hmm motifsIt${it}_align.fasta 


# g) update variables, nettoyage fichiers temporaires

    profile=Profile_${it}.hmm
    size=$(gawk 'NR==3{print($2)}' $profile)
    rm *tmp

done


#--------------------------------------------------------------------
## 3. Export Best Profile


##we select the profile with the best results (max number of amino acids identified within repeats)
best=$(gawk 'BEGIN{nb=0}{if($8>nb){nb=$8;profile=$1}}END{print(profile)}' RES_construction_profile.txt)

echo "Retrieve $best as best profile"


if [[ ! -e $OUT_DIR ]];then
    mkdir -p $OUT_DIR ;
fi


gawk -v name=${outname%.hmm} '{if(NR==2){$2=name};print}' $best > $OUT_DIR/${outname}
echo "moving HMM file ${outname} to $OUT_DIR/ "
cp RES_construction_profile.txt $OUT_DIR/${outname%.hmm}.log
echo "moving log file ${outname%.hmm}.log to $OUT_DIR/ "

cd $MAIN ;

echo "END STEP 1"

#
#========================================================
#                  Fin du Script
#========================================================
