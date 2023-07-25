#!/bin/sh
#1) Prepare output directories
mkdir ../star_output
mkdir ../count_result_star
cd ../raw_data

#2) run STAR and FeatureCounts for every sample
for var1 in *_1.fq.gz
do
  var2="$(echo $var1 | sed -e "s/_1.fq/_2.fq/g")"
  outcount="$(echo $var1 | sed -e "s/_1.fq.gz/.count/g")"
  dirname="$(echo $var1 | sed -e "s/_1.fq.gz//g")"
  mkdir ../star_output/$dirname
  STAR --runThreadN 35 --genomeDir /Ramen12TB/mitchell/genome_celegans --readFilesIn $var1 $var2 --readFilesCommand zcat --outSAMtype BAM SortedByCoordinate --outFileNamePrefix ../star_output/$dirname/
  featureCounts -F GTF -a /Ramen12TB/mitchell/genome_celegans/ce11_RefSeq.gtf -g gene_id -t exon -o ../count_result_star/$outcount -s 0 -p -C -B -T 32 ../star_output/$dirname/Aligned.sortedByCoord.out.bam
done