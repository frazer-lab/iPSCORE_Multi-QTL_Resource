#!/bin/bash

#$ -N peaks
#$ -V -cwd
#$ -o logs/peaks
#$ -e logs/peaks
#$ -pe smp 12

source /frazer01/home/jennifer/.bash_profile
source activate encode-chip

set -e

bam=$1
name=${bam::-4}
out_dir=$2

log=${out_dir}/pipeline.log

date >& 2

echo "Create tagalign file" >> $log

date >& 2

cmd="samtools sort -@ 12 -n -o ${name}.nsort.bam ${bam}"
echo $cmd >> $log; echo $cmd >& 2; eval $cmd

date >& 2

cmd="samtools fixmate -@ 12 ${name}.nsort.bam ${name}.nsort.fixed.bam"
echo $cmd >> $log; echo $cmd >& 2; eval $cmd

date >& 2

cmd="bedtools bamtobed -bedpe -mate1 -i ${name}.nsort.fixed.bam | gzip -nc > ${name}.tmp"
echo $cmd >> $log; echo $cmd >& 2; eval $cmd
zcat ${name}.tmp | awk 'BEGIN{OFS="\t"}{printf "%s\t%s\t%s\tN\t1000\t%s\n%s\t%s\t%s\tN\t1000\t%s\n",$1,$2,$3,$9,$4,$5,$6,$10}' | gzip -nc > ${name}.tagAlign.gz

date >& 2

cmd="rm ${name}.nsort.bam ${name}.nsort.fixed.bam ${name}.tmp"
echo $cmd >> $log; echo $cmd >& 2; eval $cmd

date >& 2

echo "Saved: ${name}.tagAlign.gz" >> $log
