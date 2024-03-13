#!/bin/bash

source /frazer01/home/jennifer/.bash_profile
source activate encode-atac

set -e

script_dir=/frazer01/home/jennifer/software/encode-atac-seq-pipeline

dir=`mktemp -d -p scratch`
if [ ! -d ${dir} ]; then exit 1; fi
fq_dir=$1
merge_uuid=$2
map_file=$3 #premerge to merge map file (column 1 premerge uuids, column 2 merge uuids)
out_dir=$4
log=$5

echo fq_dir=$fq_dir
echo merge_uuid=$merge_uuid
echo map_file=$map_file
echo out_dir=$out_dir
echo log=$log

# 1. concatenate fastqs from multiple sequencing runs
date >& 2

premerge_uuids=(`grep $merge_uuid $map_file | cut -f1`)

for premerge in ${premerge_uuids[@]}; do
    cmd="rsync ${fq_dir}/${premerge}/*.fastq.gz ${dir}"
    echo $cmd >& 2; echo $cmd >> $log; eval $cmd
done

date >& 2

cmd="zcat ${dir}/*R1*.fastq.gz > ${dir}/R1.fastq"
echo $cmd >& 2; echo $cmd >> $log; eval $cmd

date >& 2

cmd="zcat ${dir}/*R2*.fastq.gz > ${dir}/R2.fastq"
echo $cmd >& 2; echo $cmd >> $log; eval $cmd

cmd="rm ${dir}/*.fastq.gz"
echo $cmd >& 2; echo $cmd >> $log; eval $cmd

in_files=( ${dir}/R1.fastq ${dir}/R2.fastq )
echo fastq input.. ${in_files[@]} >& 2

# 2. Detect and trim adapters
#fqs=($fastqs)
for fq in ${in_files[@]}; do
    prefix=`basename $fq .fastq.gz`
	out=${out_dir}/detect_adapters/${prefix}.adapter.txt
	cmd="python3 ${script_dir}/src/detect_adapter.py $fq > $out" # get script from encode github
	echo $cmd >& 2; echo $cmd >> $log; eval $cmd
done

date >& 2


ADAPTER_FWD=`cat ${out_dir}/detect_adapters/R1*.adapter.txt`
ADAPTER_REV=`cat ${out_dir}/detect_adapters/R2*.adapter.txt`
echo "ADAPTER_FWD" $ADAPTER_FWD
echo "ADAPTER_REV" $ADAPTER_REV
cmd="cutadapt --cores 4 -a $ADAPTER_FWD -A $ADAPTER_REV -o ${out_dir}/detect_adapters/reads.1.fastq -p ${out_dir}/detect_adapters/reads.2.fastq ${in_files[0]} ${in_files[1]}"
echo $cmd >& 2; echo $cmd >> $log; eval $cmd

date >& 2


