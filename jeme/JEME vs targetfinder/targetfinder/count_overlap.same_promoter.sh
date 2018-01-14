a=$1
b=$2
thresh=$3
filename=overlap_${a}_${b}.same_promoter.txt
> ${filename}
for cell in  K562 GM12878
do
    wc -l ${cell}/windows.promoter_name.${a}.bed | awk '{printf("%s ", $1)}' >> ${filename}
    wc -l ${cell}/windows.promoter_name.${b}.bed | awk '{printf("%s ", $1)}' >> ${filename}
    echo -n "$cell " >> ${filename}
    for i in 1e-9 0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9 0.99
    do 
        intersectBed  -a ${cell}/windows.promoter_name.${a}.bed -b ${cell}/windows.promoter_name.${b}.bed -f $i -F $i -wo | awk '$1==$6 && ($2!=$7 || $3!=$8) && ($4==$9 || $5==$10)' | cut -f 1-5 | sort -u | wc -l
    done | xargs echo -n >> ${filename}
    echo >> ${filename}
done
