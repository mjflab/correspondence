a=$1
b=$2
thresh=$3
filename=eep_overlap_${a}_${b}.txt
> ${filename}
for cell in  K562 GM12878 HeLa-S3 IMR90 NHEK HUVEC
do
    wc -l ${cell}/output-eep/enhancers.${a}.bed | awk '{printf("%s ", $1)}' >> ${filename}
    wc -l ${cell}/output-eep/enhancers.${b}.bed | awk '{printf("%s ", $1)}' >> ${filename}
    echo -n "$cell " >> ${filename}
    for i in 1e-9 0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9 0.99
    do 
        intersectBed  -a ${cell}/output-eep/enhancers.${a}.bed -b ${cell}/output-eep/enhancers.${b}.bed -f $i -F $i -c | awk -v thresh=$thresh '$5 > thresh' | wc -l
    done | xargs echo -n >> ${filename}
    echo >> ${filename}
done
