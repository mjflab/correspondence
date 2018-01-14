a=$1
b=$2
filename=overlap_${a}_${b}.same_promoter.txt
> ${filename}
for cell in  K562 GM12878 HeLa-S3 IMR90 NHEK HUVEC
do
    wc -l ${cell}/output-epw/windows.promoter_name.${a}.bed | awk '{printf("%s ", $1)}' >> ${filename}
    wc -l ${cell}/output-epw/windows.promoter_name.${b}.bed | awk '{printf("%s ", $1)}' >> ${filename}
    echo -n "$cell " >> ${filename}
    for i in 1e-9 0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9 0.99
    do 
        intersectBed -a ${cell}/output-epw/windows.promoter_name.${a}.bed -b ${cell}/output-epw/windows.promoter_name.${b}.bed -wo -f $i -F $i  | awk '!($1==$6 && $2==$7 && $3==$8) && ($4==$9 || $5==$10)' | cut -f 1-5 | sort -u | wc -l
    done | xargs echo -n >> ${filename}
    echo >> ${filename}
done
