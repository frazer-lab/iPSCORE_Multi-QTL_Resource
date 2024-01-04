#!/bin/bash

#$ -N peaks
#$ -V -cwd
#$ -pe smp 4

date >& 2

set -e

out_dir=$1
log=$2
prefix=${out_dir}/peaks/narrow
pval_thresh=0.01
smooth_window=$3 # 150
shiftsize=$4 # 75
gensz=hs

if [ ! -d ${out_dir}/peaks ]; then mkdir ${out_dir}/peaks; fi

source /frazer01/home/jennifer/.bash_profile
source activate encode-atac

date >& 2; date >> $log
tag=${out_dir}/Aligned.sorted.filt.nodup.nomito.tn5.tagAlign.gz
cmd="macs2 callpeak -t $tag -f BEDPE -n ${prefix}_tn5_tagAlign -g $gensz -q $pval_thresh --shift $shiftsize --extsize $smooth_window --nomodel -B --SPMR --keep-dup all --call-summits"
echo $cmd >& 2; echo $cmd >> $log; eval $cmd; date >& 2

