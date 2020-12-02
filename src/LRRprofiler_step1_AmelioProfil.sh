#!/bin/bash
#========================================================
# PROJET : LRRannotation
# SCRIPT : LRRannotation_step1_AmelioProfil.sh
# AUTHOR : Celine Gottin
# CREATION : 2020.01.08
#========================================================
# DESCRIPTION : 
# ARGUMENTS : o $1 : Proteome (fasta)
#             o $2 : Liste proteines
#             o $3 : Profil initial Path
#             o $4 : Profil final Path
#             o $5 : Profil final name
# DEPENDENCIES : o HMMER v 3.1b2
#                o MAFFT v 7.313
#========================================================


#========================================================
#             Environment, variables, files
#========================================================

# Modules
module load bioinfo/hmmer/3.1b2
module load bioinfo/mafft/7.313

# on error exit flag : set -e
set -o errexit

# error if a var is unset : set -u
set -o nounset

# raise error in pipe
set -o pipefail


# Variables
function quit_pb_option() {
    printf "\nOptions : --in_proteome ; --list_proteins ; --in_profile ; --out_dir ; --out_profile_name ; --dev\n"
    exit 1
}

if (( $# == 0)); then
  quit_pb_option
fi

Proteome=""
ListProt=""
initProfile=""
OUT_DIR="RES_step1_AmelioProfile_$(date +'%Y%m%d_%H%M%S')"
outname="new_profile.hmm"
itMax=15

devopt=false

while (( $# > 0 )); do
    case "$1" in
	--in_proteome)
	    Proteome=$(readlink -f "$2"); shift 2;
	    if [[ ! -f $Proteome ]];then
		    echo "File $Proteome does not exist"
		    quit_pb_option	
	    fi
	    ;;
	--list_proteins)
	    ListProt=$(readlink -f "$2"); shift 2;
	    if [[ ! -f $ListProt ]];then
		    echo "File $ListProt does not exist"
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
  --out_dir)
      OUT_DIR="$2"; shift 2;
      if [ ! -e  $OUT_DIR ]; then
        mkdir "$OUT_DIR"
      fi
      ;;
	--out_profile_name)
	    outname="$2"; shift 2
	    ;;
  --dev)
      devopt=true; shift 1
      ;;
	*)
      echo "Option $1 is unknown please ckeck your command line"
      quit_pb_option
      ;;
    esac
done

#for param in $Proteome $initProfile $ListProt
#do
#    if [[ $param -eq "" ]];then
#	echo "error : mandatory parameter is missing"
#    fi
#done



MAIN=$(pwd)

# Working dir
WD=wd_amelioProfil_$(date +'%H%M%S')

if [[ -e $WD ]];then
    rm -r $WD
fi

mkdir $WD ; cd $WD


#========================================================
#                       Script
#========================================================

echo "============================================="
echo "Proteome : $(basename $Proteome)"
echo "Initial Profile : $(basename $initProfile)"
echo "List of proteins : $(basename $ListProt)"
echo "Max number of iteration : $itMax "
echo "============================================="

#--------------------------------------------------------------------
# 1. Extract protein sample from protein List

gawk '{if(NR==FNR){P[">"$1]=1}else{if($1~/>/){if(P[$1]==1){pr=1}else{pr=0}};if(pr==1){print}}}' $ListProt $Proteome > proteome_tmp.fasta

nbPf=$(grep -c ">" proteome_tmp.fasta)
nbPe=$(grep -c . $ListProt)

echo "found $nbPf proteins over $nbPe expected in the list"

gawk 'BEGIN{OFS=";";protName=""}{if($1~/^>/){if(length(protName)>1){print(protName,sequence)};gsub(">","");protName=$1;sequence=""}else{sequence=sequence""$1}}END{print(protName,sequence)}' proteome_tmp.fasta > proteome_tmp.csv

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
    
    hmmsearch -E 1000 --domE 1000 -o del.tmp --nobias --noali --domtblout HMMres${it}.tmp ${profile} proteome_tmp.fasta 

    nbProtOld=$nbProt
    nbProt=$(grep -v "^#"  HMMres${it}.tmp | cut -f1 -d" "| sort -u | wc -l)



# b) Resultats des profils : Profil; nbProteine liste trouvees; idem seuil 0.5,0.1,0.01 ; nbMotifs; taille moyenne motifs

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

#    if (( $nbProt > $nbProtOld || (( $nbProt == $nbProtOld && $nbMotifs >  $nbMotifsOld )) )) ; then # we have a better HMM
#	penalty=0;
#    elif ((  $nbProt == $nbProtOld &&  $nbMotifs <=  $nbMotifsOld )) ;  then # we have an equivalent HMM, use a moderate penalty (stop after 3 such iterations) 
#	let penalty=${penalty}+2;
 #   else # we got a worse HMM use a higher penalty (stop after 2 such iterations)
#	les penalty=${penalty}+3;
#    fi

    if (( $penalty >= 3 ));then
	break ;
    fi


# e) extract LRR motifs from the protein sequence

    gawk -F";" 'BEGIN{OFS=";"}{if(NR==FNR){SEQ[$1]=$2}else{len=$7-$6+1;L[$1]=len;seq=substr(SEQ[$1],$6,L[$1]);print(">"$1"."$2);print(seq)}}' proteome_tmp.csv motifs.tmp > motifsIt${it}.fasta


# f) Align sequence using MAFFT and build new profile
    
    echo "Building new HMM ..."
    #~/bin/famsa-1.2.5-linux-static -go 9 motifsIt${it}.fasta motifsIt${it}_align.fasta
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

echo "Conserve $best"


if [[ ! -e $OUT_DIR ]];then
    mkdir -p $OUT_DIR ;
fi


gawk -v name=${outname%.hmm} '{if(NR==2){$2=name};print}' $best > $OUT_DIR/${outname}
cp RES_construction_profile.txt $OUT_DIR/${outname%.hmm}.log
echo "moving logfile ${outname%.hmm}.log to $OUT_DIR/ "
cp $ListProt $OUT_DIR/.
echo "moving list $ListProt to $OUT_DIR/ "


cd $MAIN ;

if [[ $devopt==false ]];then
   rm -r $WD
fi

echo "END STEP 1"

#
#========================================================
#                  Fin du Script
#========================================================
