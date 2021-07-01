#!/bin/bash
#========================================================
# PROJET : LRRprofiler
# SCRIPT : LRRannotation_step2_RechercheLRR.sh
# AUTHOR : Celine Gottin
# CREATION : 2020.01.09
#========================================================
# DESCRIPTION : 
# ARGUMENTS : o $1 : Proteome (fasta)
#             o $2 : Name
#             o $3 : New Profile RLK
#             o $4 : New Profile NLR
#             o $5 : output directory
# DEPENDENCIES : o HMMER v 3.1b2
#                o MAFFT v 7.313

#========================================================


#========================================================
#                Environment & variables
#========================================================

# variables
function quit_pb_option() {
    printf "\nOptions : --in_proteome ; --name ; --rlk_profile ; --nlr_profile ; --out_dir ; --dev\n"
    exit 1
}

if (( $# == 0)); then
  quit_pb_option
fi

MAIN=$(pwd)
SCRIPT=${LG_SCRIPT}

PROTEOME=""
NAME="LRRprofiler"
OUT_DIR="" 
LRR_RLK=""
LRR_NLR=""

devopt=0

while (( $# > 0 )); do
    case "$1" in
	--in_proteome)
	    PROTEOME=$(readlink -f "$2"); shift 2;
	    if [[ ! -f $PROTEOME ]];then
		    echo "File $PROTEOME does not exist"
		    quit_pb_option	
	    fi
	    ;;
	--name)
	    NAME="$2"; shift 2
	    ;;
	--rlk_profile)
	    LRR_RLK=$(readlink -f "$2"); shift 2;
	    if [[ ! -f $LRR_RLK ]];then
		    echo "File $LRR_RLK does not exist"
		    quit_pb_option
	    fi
	    ;;
	--nlr_profile)
      LRR_NLR=$(readlink -f "$2"); shift 2;
      if [ ! -f  $LRR_NLR ]; then
        echo "File $LRR_NLR does not exist"
		    quit_pb_option
      fi
      ;;
	--out_dir)
      OUT_DIR="$2"; shift 2;
      if [ ! -e  $OUT_DIR ]; then
        mkdir "$OUT_DIR"
      fi
	    ;;
	--dev)
      devopt=1; shift 1
      ;;
	*)
            echo "Option $1 is unknown please ckeck your command line"
            quit_pb_option
            ;;
    esac
done


# working dir
WD=$MAIN/wd_rechercheLRR_${NAME}_$(date +'%H%M%S')

if [[ -e $WD ]];then
    rm -r $WD
fi

mkdir $WD; cd $WD

#========================================================
#                       SCRIPT
#========================================================

cat $MAIN/Res_$NAME/Res_step1/Liste* > ${NAME}_PriorClassif.txt

ln -s ${LG_HMMlib}/SMART_LRR_*.hmm .
ln -s ${LG_HMMlib}/LRR_NLR_ORYSJ_canonic_refined.hmm .
ln -s ${LRR_RLK} .
ln -s ${LRR_NLR} .



## A . SEARCH FOR LRR MOTIFS WITH HMM
##------------------------------------

## a. Extract LRR motifs from several hmmsearch results

for profile in *.hmm
do
    outfile=${profile%.hmm}.tbl
    
    hmmsearch -o del.tmp -E 1000 --domE 1000 -Z 500 --nobias --noali --domtblout ${outfile} $profile $PROTEOME

    if [[ $(wc -l ${outfile} | gawk '{print $1}') -le 13 ]];then
	    rm ${outfile} #empty file
    else
	    python3.6 $SCRIPT/Extract_LRR_motifs.py -s $PROTEOME -t ${outfile} -o ${outfile%.tbl}.csv
    fi
    
done

## b. Concatenate all LRR motifs
python3 $SCRIPT/Concat_all_motifs.py -s $PROTEOME -d $(pwd) -o $(pwd)/LRR_${NAME}_concat.tmp

## c. Look for LRR motifs in interMotifs regions
ln -s $MAIN/HMM_lib/SMART* .

grep "interLRR" LRR_${NAME}_concat.tmp | gawk -F";" '($5-$4>8){print(">"$1";"$4";"$5"\n"$8)}' - > LRR_${NAME}_interLRR.fasta

for profile in *.hmm
do
    outfile=LRR_${NAME}_${profile%.hmm}_interLRR_search.tbl
    
    hmmsearch -o del.tmp -E 1000 --domE 1000 --noali --nobias --domtblout ${outfile} $profile LRR_${NAME}_interLRR.fasta

    if [[ $(wc -l ${outfile} | gawk '{print $1}') -le 13 ]];then
	    rm ${outfile} #enpty file
    else
	    python3.6 $SCRIPT/Extract_interLRR_motifs.py -s LRR_${NAME}_interLRR.fasta -t ${outfile} -o $(pwd)/${outfile%.tbl}.csv
    fi
    
done

## e. Concatenate all motifs
python3.6 $SCRIPT/Concat_all_motifs.py -s $PROTEOME -d $(pwd) -o $(pwd)/LRR_${NAME}_ALLMOTIFS_init.csv 



## B . SEARCH FOR LRR MOTIFS WITH BLAST
##---------------------------------------

# Extract noLRR regions with length over 4 aa : interLRR extracted from LRR_${NAME}_ALLMOTIFS.csv and regions 0-LRR et LRR-END
# if NBARC or f-box --> after LRR end domain position, if kinase --> before kinase strat position

gawk -F";" 'BEGIN{OFS=";"}{
     if(NR==FNR){
          gsub(" ",";");
          CLASS[$1]=$2;
          if($2~/Kinase/){
               lim1[$1]=$3}
          else{
               lim2[$1]=$4}}
     else{
          if($2!~/interLRR/){
               if($1!=prot){
                    if(CLASS[$1]!~/NBARC/ || (CLASS[$1]~/NBARC/ && $4>lim2[$1])){prot=$1;S[$1]=$(4)-1;E[$1]=$(5)+1};}
               else{
                    if(CLASS[$1]!~/Kinase/ || (CLASS[$1]~/Kinase/ && $5<lim1[$1])){E[$1]=$(5)+1}}}}
}END{for(i in S){
     if(lim1[i]!=""){
          print(i,1,S[i]);print(i,E[i],lim1[i]-1)}
     else{
          if(lim2[i]){
               print(i,lim2[i]+1,S[i]);print(i,E[i],"end")}
          else{print(i,1,S[i]);print(i,E[i],"end")}}}
}' ${NAME}_PriorClassif.txt LRR_${NAME}_ALLMOTIFS_init.csv > extraLRR_dom.tmp

gawk 'BEGIN{OFS=";";prot="";seq=""}{if($1~/>/){if(prot!=""){print(prot,seq);prot="";seq=""};gsub(">","");prot=$1}else{gsub("*","");seq=seq""$1}}END{print(prot,seq)}' $PROTEOME | gawk -F";" 'BEGIN{OFS=";"}{if(NR==FNR){SEQ[$1]=$2}else{if($3=="end"){end=length(SEQ[$1]);string=substr(SEQ[$1],$2);}else{string=substr(SEQ[$1],$2,$3-$2+1);end=$3};if(end-$2>10){print(">"$1";"$2";"end"\n"string)}}}' - extraLRR_dom.tmp > extraLRR_regions.fasta

# Create fasta for blast database
cat LRR_${NAME}_interLRR.fasta extraLRR_regions.fasta > ${NAME}_noLRR_regions.fasta
#gawk -F";" 'BEGIN{OFS=";"}{if(NR==FNR){if($2~/interLRR/ && $5-$4>=8){print(">"$1,$4,$5);print($8)}}else{if($3-$2>=8){print(">"$1,$2,$3);print($4)}}}' LRR_${NAME}_ALLMOTIFS_init.csv extraLRR_regions.txt > ${NAME}_noLRR_regions.fasta

makeblastdb -in ${NAME}_noLRR_regions.fasta -dbtype prot -parse_seqids

# Select LRR motifs with eval<5 for Query file
gawk -F";" 'BEGIN{OFS=";"}NR>1{if($2!~/interLRR/ && $7<=5 && $6>=20 && $6<=26){print(">"$1,$4,$5);print($8)}}' LRR_${NAME}_ALLMOTIFS_init.csv > ${NAME}_LRR_seq.fasta

time blastp -db ${NAME}_noLRR_regions.fasta -query ${NAME}_LRR_seq.fasta -max_hsps 2 -evalue 1.0 -out testblast.out -outfmt "6 qseqid sseqid length nident positive gaps bitscore evalue qstart qend sstart send qlen sseq"

#filtering results
sort -k2,2 -Vk11,11 testblast.out | gawk 'BEGIN{OFS="\t"}{if($3>7 && $5/$3>0.6){print($0)}}' - > Resblast.tmp

#overlopping of more than 80%, we keep the hit with best eval
gawk -F"\t" 'BEGIN{OFS="\t";current="";line=""}{
      if($2==current && ((oldE-$11+1)/($12-$11))>0.8){
           if($8<evalue){
                line=$0;
                oldS=$11;
                oldE=$12;
                evalue=$8;}}
      else{
           gsub("\t",";",line);
           if(line!=""){print(line)};
           current=$2;
           oldS=$11;
           oldE=$12;
           evalue=$8;
           line=$0;}
}' Resblast.tmp > Resblast.txt

#processed
python3.6 $SCRIPT/Extract_blast_motifs.py -s $PROTEOME -b Resblast.txt -o motifsblast.csv

##suppr csv to avoid redumdancies with next step
rm *${NAME}.csv
rm *search.csv
rm SMART_LRR_CC.csv

#concatenate
python3.6 $SCRIPT/Concat_all_motifs.py -s $PROTEOME -d $(pwd) -o $(pwd)/LRR_${NAME}_ALLMOTIFS.csv

#SAVE results
cp LRR_${NAME}_ALLMOTIFS.csv $OUT_DIR/.
cp ${NAME}_PriorClassif.txt $OUT_DIR/.

cd $MAIN 

if [[ $devopt -eq 0 ]];then
  rm -r $WD
fi

echo "END STEP 2"



