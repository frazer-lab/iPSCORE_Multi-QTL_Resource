#!/bin/bash

#$ -N counts
#$ -V -cwd
#$ -pe smp 4

date >& 2 

set -e

source /frazer01/home/jennifer/.bash_profile
source activate encode-atac

script_dir=/projects/PPC/pipeline/ATAC-Seq/work/2023_0720/scripts

out_dir=$1
ref_saf=$2
name=`basename $saf`
name=${name::-4}

if [ ! -d ${out_dir}/ref_peaks ]; then mkdir ${out_dir}/ref_peaks; fi

bam=${out_dir}/Aligned.sorted.filt.nodup.nomito.bam

date >& 2

# 2. run feature counts
# -B = only count read pairs that have both ends aligned
# -C = do not cound read pairs that have their two ends mapping to diff chr
# -T = threads
# -F = input format
# -p = paired-end reads
# --countReadPairs = fragments will be counted instead of reads

cmd="featureCounts -p --countReadPairs -B -C -T 4 -F SAF -a ${ref_saf} -o ${out_dir}/ref_peaks/ref_peaks.counts ${bam}"
echo $cmd >& 2; eval $cmd

date >& 2
