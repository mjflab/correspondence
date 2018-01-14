a=pos; 
b=pos; 
thresh=${1}; 
for cell in K562 GM12878 HUVEC HeLa-S3 NHEK IMR90; 
do  
    intersectBed -a ${cell}/output-eep/enhancers.${a}.bed -b ${cell}/output-eep/enhancers.${b}.bed -wo -f $thresh -F $thresh | cut -f 4,8 | awk 'BEGIN{FS=OFS="\t"}{if($1 > $2){temp=$1; $1=$2; $2=temp; print} else if($1 != $2){print}}' | sort -u > ${cell}/output-eep/overlap_${a}_${b}_${thresh}_pairs.txt; 
done
