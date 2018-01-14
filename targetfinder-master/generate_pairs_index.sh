a=pos; 
b=pos; 
thresh=${1}; 
for cell in K562 GM12878 HUVEC HeLa-S3 NHEK IMR90; 
do  
    intersectBed -a ${cell}/output-epw/windows.${a}.bed -b ${cell}/output-epw/windows.${b}.bed -wo -f $thresh -F $thresh | cut -f 4,11 | awk 'BEGIN{FS=OFS="\t"}{if($1 > $2){temp=$1; $1=$2; $2=temp; print} else if($1 != $2){print}}' | sort -u > ${cell}/output-epw/overlap_${a}_${b}_${thresh}_pairs.txt; 
done
