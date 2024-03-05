#!/bin/bash

#$ -pe smp 8
#$ -V -cwd

set -e 

out_dir=$1
log=$2

bam=${out_dir}/Aligned.sorted.bam
nomito=${out_dir}/Aligned.sorted.nomito.bam
tmp=${out_dir}/Aligned.sorted.nomito.sorted.bam
mdup=${out_dir}/Aligned.sorted.nomito.mdup.bam
mdup_qc=${out_dir}/qc/Aligned.sorted.nomito.mdup.dup.qc
nodup=${out_dir}/Aligned.sorted.nomito.mdup.nodup.bam

# 1. Make mito-free bam
date >& 2

cmd="samtools idxstats $bam | cut -f 1 | grep -v -P "^chrM" | xargs samtools view $bam -@ 8 -b > ${nomito}"
echo $cmd >& 2; echo $cmd >> $log; eval $cmd

cmd="samtools sort -@ 8 -o ${tmp} ${nomito}"
echo $cmd >& 2; echo $cmd >> $log; eval $cmd

cmd="samtools index -@ 8 ${tmp} ${tmp}.bai"
echo $cmd >& 2; echo $cmd >> $log; eval $cmd

# 2. Mark duplicates
cmd="picard MarkDuplicates INPUT=${tmp} OUTPUT=${mdup} METRICS_FILE=${mdup_qc} VALIDATION_STRINGENCY=LENIENT ASSUME_SORTED=true REMOVE_DUPLICATES=false"
echo $cmd >& 2; echo $cmd >> $log; eval $cmd

cmd="rm ${tmp}"

# 3. Remove duplicates
cmd="samtools view -F 1024 -@ 8 -o ${nodup} ${mdup}"
echo $cmd >& 2; echo $cmd >> $log; eval $cmd

# 4. flagstat
cmd="samtools flagstat -@ 8 ${nodup} > ${out_dir}/qc/Aligned.sorted.nomito.mdup.nodup.flagstat.qc"
echo $cmd >& 2; echo $cmd >> $log; eval $cmd

# 5. clean
cmd="rm ${nomito} ${tmp} ${mdup} ${nodup}"
echo $cmd >& 2; echo $cmd >> $log; eval $cmd