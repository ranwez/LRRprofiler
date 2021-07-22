#!/bin/bash
#========================================================
# PROJET : LRRprofiler
# SCRIPT : LRRprofiler_step3_classification.sh
# AUTHOR : Celine Gottin
#========================================================
# DESCRIPTION : 
# ARGUMENTS : o $1 : Proteome (fasta)
#             o $2 : Name
# DEPENDENCIES : o HMMER v 3.1b2
#                o TMHMM 2.0c
#========================================================


#========================================================
#                Environment & variables
#========================================================


# variables
function quit_pb_option() {
    printf "\nOptions : --in_proteome <fastafile> ; --name <string>; --dev\n"
    exit 1
}

if (( $# == 0)); then
  quit_pb_option
fi

MAIN=$(pwd)
SCRIPT=${LG_SCRIPT}

PROTEOME=""
NAME=""
OUT_DIR="RES_step3_classif_$(date +'%Y%m%d_%H%M%S')"
devopt=0

while (( $# > 0 )); do
    case "$1" in
    --in_proteome)
        PROTEOME=$(readlink -f "$2"); shift 2;
        if [[ ! -f $PROTEOME ]];then
            echo "File $Proteome does not exist"
            quit_pb_option	
        fi
        ;;
    --name)
        NAME="$2" ; shift 2
        ;;
    *)
      echo "Option $1 is unknown please ckeck your command line"
      quit_pb_option
      ;;
    esac
done


if [[ -e $LRRPROFILER_RESDIR ]];then
    OUT_DIR=$LRRPROFILER_RESDIR/Res_step3
fi

# working dir
WD=$LRRPROFILER_TMP/wd_classificationLRR_${NAME}
mkdir $WD; cd $WD

#========================================================
#                       SCRIPT
#========================================================

#ln -s $WD/Output_PriorClassif_${NAME}/${NAME}_PriorClassif.txt .
ln -s $LRRPROFILER_RESDIR/Res_step2/LRR_${NAME}_ALLMOTIFS.csv .
ln -s $LRRPROFILER_RESDIR/Res_step0_itak/${NAME}_shiu_alignment.txt .

cat $LRRPROFILER_TMP/Liste* > ${NAME}_PriorClassif.txt

## a. Extract LRR sequences in fasta format
gawk '{if(NR==FNR){split($0,T,";");P[">"T[1]]=1}else{if($1~/>/){if(P[$1]==1){p=1}else{p=0}};if(p==1){print}}}' LRR_${NAME}_ALLMOTIFS.csv $PROTEOME > LRR_${NAME}.fasta

# prot size
gawk 'BEGIN{OFS=";"}{if($1~/>/){if(length(prot)>0){print(prot,L)};gsub(">","");prot=$1;L=0}else{L=L+length($0)}}END{if(length(prot)>0){print(prot,L)}}' LRR_${NAME}.fasta > LRR_prot_size.csv

## b. Search for functionnal domains
for hmm in F-box FBD NB-ARC Malectin Malectin_like TIR Cys-typeA Cys-typeB Cys-typeC Cys-typeD Cys-typeE
do
    hmmsearch -o del.tmp --nobias --noali --domtblout ${hmm}_$NAME.tbl ${LG_HMMlib}/${hmm}.hmm LRR_${NAME}.fasta
done

## c. Extract Domains
gawk 'BEGIN{OFS=";"}{if(NR==FNR){gsub(";"," ");P[$1]=1}else{if(P[$1]==1){print($1,$2,$3,$4)}}}' LRR_prot_size.csv ${NAME}_PriorClassif.txt > LRR_domains.tmp ##Kinase, NBARC

for file in *.tbl
do
    gawk 'BEGIN{OFS=";"}$1!~/#/{print($1,$4,$18,$19)}' $file >> LRR_domains.tmp
done

## tmhmm and keep peptide signal
cat LRR_${NAME}.fasta | ${LG_TMHMM} -noplot > TMHMM_out.txt

gawk 'BEGIN{OFS=";"}{if($3~/TMhelix/){if($4>30){print($1,"TM",$4,$5)}else{print($1,"PS",$4,$5)}}}' TMHMM_out.txt >> LRR_domains.tmp

##LRR
gawk -F";" 'BEGIN{OFS=";"}NR>1{if($2!="interLRR"){print($1,$2,$4,$5)}}' LRR_${NAME}_ALLMOTIFS.csv >> LRR_domains.tmp

sort -V -t";" -k1,1 -k3,3 -k4,4r LRR_domains.tmp | gawk -F";" 'BEGIN{OFS=";"}{if($0!=""){print}}' - > LRR_domains.csv

## d. Structure filtering
# rename LRR motifs according to observed consensus

gawk -F";" 'BEGIN{OFS=";"}{
       gsub("SMART_LRR_align","LRR_TYP");
       gsub("SMART_LRR_TYP","LRR_RI");
       gsub("SMART_LRR_RI","LRR_CC");
       gsub("SMART_LRR_CC","LRR_PS");
       gsub("SMART_LRR_BAC","LRR_BAC");
       if($2~/LRR_kinase/){$2="LRR_PS"};
       if($2~/LRR_NLR/){$2="LRR_NLR"};
       if($2~/LRR_Fbox/){$2="LRR_FBOX"};
       if($2~/Cys/){$2="Cys-Pair"};print}' LRR_domains.csv > LRR_domains_filtered.tmp

##Correct blast resultats : if inside another domain --> suppr
gawk -F";" 'BEGIN{OFS=";";prot=""}{if(prot!=$1){prot=$1;lim=0};if($2!~/BLAST/){if($4>lim){print;lim=$4}}else{if($3>lim){print}}}' LRR_domains_filtered.tmp > tmp

## if less than 24aa from TM --> suppr; if less than 60 aa from start --> suppr; if less than 40 aa after NBARC --> suppr
gawk -F";" 'BEGIN{OFS=";"}{if(NR==FNR){if($2~/TM/){Lim[$1]=$(3)-24};if($2~/NBARC/){Lim2[$1]=$(4)+40}}else{if($2~/BLAST/ && $3<60){NEXT}else{if(Lim[$1]!="" && $2~/BLAST/ && $4>Lim[$1]){NEXT}else{if(Lim2[$1]!="" && $2~/BLAST/ && $3<Lim2[$1]){NEXT}else{print}}}}}' tmp tmp > LRR_domains_filtered.csv

rm *.tmp

## e. Protein classification (NLR, RLK, RLP, F-box and other)

## Functionnal Domain list per proteins
gawk -F";" 'BEGIN{OFS=";"}{if(P[$1]){P[$1]=P[$1]"-"$2}else{P[$1]=$2}}END{for(i in P){print(i,P[i])}}' LRR_domains_filtered.csv > tmp
#gawk -F";" 'BEGIN{OFS=";"}{if($2!~/LRR/){if(P[$1]){P[$1]=P[$1]"-"$2}else{P[$1]=$2}}}END{for(i in P){print(i,P[i])}}' LRR_domains_filtered.csv > tmp

## Classif RLK; NLR and f-box
gawk -F";" 'BEGIN{OFS=";"}{if($2~/NBARC/){print($1,"NLR")}else{if($2~/Fbox/ || $2~/F-box/ || $2~/FBD/){print($1,"F-box")}else{if($2~/Kinase/){print($1,"RLK")}else{if($2~/Malectin/ || $2~/TM/ || $2~/Cys-Pair/){print($1,"RLP")}else{print($1)>"putativeRLP.tmp"}}}}}' tmp > LRR_classification.tmp

## For putative RLP : 
#   - pct = % of PS type motifs
#   - pct2 = % of NLR type motifs

#  if pct1 > 65% and more than 13 motifs in protein => RLP
#  if pct1> 10% and pct2 > 40% and pct1+pct2>90% => PIRL
#gawk -F";" 'BEGIN{OFS=";"}{if(NR==FNR){if($2~/LRR/){LRR[$1]++;if($2~/LRR_PS/){PS[$1]++}else{if($2~/LRR_NLR/){NLR[$1]++}}}}else{pct=PS[$1]/LRR[$1];if(pct>0.65 && PS[$1]>13){print($1,"RLP")}else{pct2=NLR[$1]/LRR[$1];if(pct2>0.4 && pct>0.1 && pct+pct2>0.9){print($1,"PIRL")}else{print($1,"other")}}}}' LRR_domains_filtered.csv putativeRLP.tmp >> LRR_classification.tmp
gawk -F";" 'BEGIN{OFS=";"}{if(NR==FNR){if($2~/LRR/){LRR[$1]++;if($2~/LRR_PS/){PS[$1]++}else{if($2~/LRR_NLR/){NLR[$1]++}}}}else{pct=PS[$1]/LRR[$1];if(pct>0.65 && PS[$1]>13){print($1,"RLP")}else{print($1,"other")}}}' LRR_domains_filtered.csv putativeRLP.tmp >> LRR_classification.tmp


## Add number of identified LRR

gawk -F"[; ]" 'BEGIN{OFS=";"}{if(NR==FNR){if($2~/LRR/ || $2~/BLAST/){P[$1]++}}else{print($1,$2,P[$1])}}' LRR_domains_filtered.csv LRR_classification.tmp > LRR_classification.csv

## f. structure 
cp $SCRIPT/LRR_structure.Rmd .
cp $SCRIPT/render.R .
#cat > render.R <<EOF
#library(rmarkdown)
#render("LRR_structure.Rmd", output_format = "html_document", output_file = "$(pwd)/LRR_structure_${NAME}.html", params = list(domains="${WD}/LRR_domains_filtered.csv",sizes="${WD}/LRR_prot_size.csv",class="${WD}/LRR_classification.csv"))
#EOF

R CMD BATCH '--args $(pwd)' render.R

# SAVE results
mkdir -p $OUT_DIR

if [[ ! -e LRR_classification.csv ]];then
   echo "File LRR_classification does not exist. Process has failed"
else
    cp LRR_classification.csv $OUT_DIR/.
fi

if [[ ! -e LRR_domains_filtered.csv ]];then
   echo "File LRR_domains_filtered does not exist. Process has failed"
else
    cp LRR_domains_filtered.csv $OUT_DIR/.
fi

if [[ ! -e LRR_structure.html ]];then
    echo "File LRR_structure.html was not created. This could be due to memory error."
else
    cp LRR_structure.html $OUT_DIR/.
fi

cd $MAIN 

echo "END STEP 3"





