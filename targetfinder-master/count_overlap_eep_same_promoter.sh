a=$1
b=$2
filename=eep_overlap_${a}_${b}_same_promoter.txt
> ${filename}
for cell in  K562 GM12878 HeLa-S3 IMR90 NHEK HUVEC
do
    pairs=${cell}/output-eep/pairs.bedpe;  
    cat ${pairs} | awk -v a=$a '$7==a' | cut -f 1-3 | sort -u | wc -l | awk '{printf("%s ", $1)}' >> ${filename}
    cat ${pairs} | awk -v b=$b '$7==b' | cut -f 1-3 | sort -u | wc -l | awk '{printf("%s ", $1)}' >> ${filename}
    echo -n "$cell " >> ${filename}
    for thresh in 1e-9 0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9 0.99
    do 
        pairToPair -a  <(cat ${pairs} | awk -v a=$a '$7==a') -b <(cat ${pairs} | awk -v b=$b '$7==b') -type both -f $thresh | awk -v thresh=$thresh '!($1==$8 && $2==$9 && $3==$10) && ($4==$11 && $5==$12 && $6==$13){start=$2>$9?$2:$9; end=$3<$10?$3:$10; if((end-start)/($3-$2) >= thresh && (end-start)/($10-$9)>=thresh){print}}' | cut -f 1-3 | sort -u | wc -l 
    done | xargs echo -n >> ${filename}
    echo >> ${filename}
done
