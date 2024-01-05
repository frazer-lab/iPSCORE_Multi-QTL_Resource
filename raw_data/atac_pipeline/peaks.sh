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

# 1. Call peaks
tag=${out_dir}/Aligned.sorted.filt.nodup.nomito.tn5.tagAlign.gz
cmd="macs2 callpeak \
-t $tag \
-f BEDPE \
-n ${prefix}_tn5_tagAlign \
-g $gensz \
-q $pval_thresh \
--shift $shiftsize \
--extsize $smooth_window \
--nomodel -B --SPMR --keep-dup all --call-summits"
echo $cmd >& 2; echo $cmd >> $log; eval $cmd; date >& 2

# 2. Remove blacklist
blacklist=/reference/public/ENCODE/hg38-blacklist.v2.bed
prefix=${out_dir}/peaks/narrow_tn5_tagAlign
filtered_peak=${prefix}_peaks_noblacklist.narrowPeak
bedtools intersect -v -a ${prefix}_peaks.narrowPeak -b ${blacklist} | awk 'BEGIN{OFS="\t"} {if ($5>1000) $5=1000; print $0}' | grep -P 'chr[\dXY]+[ \t]' > ${filtered_peak}
cut -f1,2,3 ${filtered_peak} | sort -u | wc -l > ${prefix}_peaks_noblacklist.narrowPeak.npeaks.txt

# 3. Call featureCounts
date >& 2; date >> $log

peak=${out_dir}/peaks/narrow_tn5_tagAlign_peaks_noblacklist.narrowPeak
cmd="cat $peak | ${script_dir}/peaks2saf.pl > ${out_dir}/saf"
echo $cmd >& 2; echo $cmd >> $log; eval $cmd; date >& 2

cmd="featureCounts -p --countReadPairs -B -C -T 4 -F SAF -a ${out_dir}/saf -o ${peak}.counts ${bam}"
echo $cmd >& 2; echo $cmd >> $log; eval $cmd; date >& 2

# 4. Calculate FRIP (using ENCODE method)
annot_bed=${out_dir}/peaks/narrow_tn5_tagAlign_peaks_noblacklist.narrowPeak
prefix=${out_dir}/peaks/narrow_tn5_tagAlign_peaks_noblacklist
ta_file=${out_dir}/Aligned.sorted.filt.nodup.tn5.tagAlign.gz
cmd="bedtools sort -i ${annot_bed} | bedtools merge -i stdin | bedtools intersect -u -nonamecheck -a ${ta_file} -b stdin | wc -l > ${prefix}.encode_frip.txt"
echo $cmd >& 2; echo $cmd >> $log; eval $cmd

# 5. Clean
cmd="rm ${out_dir}/saf"
echo $cmd >& 2; echo $cmd >> $log; eval $cmd; date >& 2



