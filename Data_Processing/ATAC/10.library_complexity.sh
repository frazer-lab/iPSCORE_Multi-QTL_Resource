#!/bin/bash

#$ -N complexity
#$ -V -cwd 
#$ -pe smp 4
#$ -o logs/complexity
#$ -e logs/complexity

source /frazer01/home/jennifer/.bash_profile
source activate encode-chip

out_dir=$1
log=$2

echo "Processing ${out_dir}" >& 2

bam=${out_dir}/Aligned.sorted.filt.nodup.bam

cmd="samtools sort -@ 4 -n -o ${out_dir}/tmp.srt.tmp ${bam}"
echo $cmd >& 2; eval $cmd

echo "TotalReadPairs    DistinctReadPairs   OneReadPair TwoReadPairs    NRF=Distinct/Total  PBC1=OnePair/Distinct   PBC2=OnePair/TwoPair" > ${out_dir}/qc/library_complexity.txt

bedtools bamtobed -bedpe -i ${out_dir}/tmp.srt.tmp | awk 'BEGIN{OFS="\t"}{print $1,$2,$4,$6,$9,$10}' | grep -v 'chrM' | sort | uniq -c | awk 'BEGIN{mt=0;m0=0;m1=0;m2=0} ($1==1){m1=m1+1} ($1==2){m2=m2+1} {m0=m0+1} {mt=mt+$1} END{printf "%d\t%d\t%d\t%d\t%f\t%f\t%f\n",mt,m0,m1,m2,m0/mt,m1/m0,m1/m2}' >> ${out_dir}/qc/library_complexity.txt

cmd="rm ${out_dir}/tmp.srt.tmp"
echo $cmd >& 2; eval $cmd
