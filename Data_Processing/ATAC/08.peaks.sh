#!/bin/bash

#$ -N peaks
#$ -V -cwd
#$ -pe smp 4

date >& 2 

set -e

out_dir=$1
log=$2
prefix=${out_dir}/peaks/narrow
pval_thresh=0.01 # default in ENCODE
gensz=hs # human genome size
smooth_window=$3 #150
shiftsize=$4 #75

#smooth_window=150 # default in ENCODE
#shiftsize=$(( -$smooth_window/2 ))

script_dir=/projects/PPC/pipeline/ATAC-Seq/work/2023_0720/scripts

source /frazer01/home/jennifer/.bash_profile
source activate encode-atac

if [ ! -d ${out_dir}/peaks ]; then mkdir ${out_dir}/peaks; fi

# 1. Call peaks
if [ ! -f ${out_dir}/peaks/narrow_tn5_tagAlign_peaks.narrowPeak ]; then
    tag=${out_dir}/Aligned.sorted.filt.nodup.nomito.tn5.tagAlign.gz
    cmd="macs2 callpeak -t $tag -f BEDPE -n ${prefix}_tn5_tagAlign -g $gensz -q $pval_thresh --shift $shiftsize  --extsize $smooth_window --nomodel -B --SPMR --keep-dup all --call-summits"
    echo $cmd >& 2; echo $cmd >> $log; eval $cmd
fi
   
# 2. Remove blacklist
BLACKLIST=/reference/public/ENCODE/hg38-blacklist.v2.bed
arr=( ${out_dir}/peaks/narrow_bam ${out_dir}/peaks/narrow_tn5_tagAlign )
for prefix in ${arr[@]}; do
    FILTERED_PEAK=${prefix}_peaks_noblacklist.narrowPeak
    bedtools intersect -v -a ${prefix}_peaks.narrowPeak -b ${BLACKLIST} | awk 'BEGIN{OFS="\t"} {if ($5>1000) $5=1000; print $0}' | grep -P 'chr[\dXY]+[ \t]'   > ${FILTERED_PEAK}
    cut -f1,2,3 ${FILTERED_PEAK} | sort -u | wc -l > ${prefix}_peaks_noblacklist.narrowPeak.npeaks.txt
done

# 3. Feature counts
bam=${out_dir}/Aligned.sorted.filt.nodup.nomito.bam
peak=${out_dir}/peaks/narrow_tn5_tagAlign_peaks_noblacklist.narrowPeak

cmd="cat $peak | ${script_dir}/peaks2saf.pl > ${out_dir}/saf"
echo $cmd >& 2; echo $cmd >> $log; eval $cmd

cmd="featureCounts -p --countReadPairs -B -C -T 4 -F SAF -a ${out_dir}/saf -o ${peak}.counts ${bam}"
echo $cmd >& 2; echo $cmd >> $log; eval $cmd

# 4. Clean
cmd="rm ${out_dir}/saf"
echo $cmd >& 2; echo $cmd >> $log; eval $cmd
