#!/bin/bash

#$ -N peaks
#$ -V -cwd
#$ -o logs/peaks
#$ -e logs/peaks
#$ -pe smp 4

source /frazer01/home/jennifer/.bash_profile
source activate encode-chip

set -e

sample_dir=$1
input_dir=$2
log_file=${sample_dir}/macs2.out

if [ -f ${sample_dir}/peaks_q0.01/broad_tagAlign_peaks.broadPeak.counts ]
then 
    echo "${sample_dir}/peaks_q0.01/broad_tagAlign_peaks.broadPeak.counts already exists"
    exit 1
fi 

if [ ! -d ${sample_dir}/peaks_q0.01 ]; then mkdir ${sample_dir}/peaks_q0.01; fi

# 0. Get fragment length
cc_qc_log=${sample_dir}/Aligned.merged.subsampled.tagAlign.cc.qc
echo $cc_qc_log
fraglen=`cut -f3 $cc_qc_log`

# 0. Setup input
if [ ! -f ${sample_dir}/Aligned.merged.bam ]
then
    bam_chip=${sample_dir}/Aligned.filt.srt.nodup.bam
    bed_chip=${sample_dir}/Aligned.filt.srt.nodup.tagAlign
else
    bam_chip=${sample_dir}/Aligned.merged.bam
    bed_chip=${sample_dir}/Aligned.merged.tagAlign
fi

bed_input=${input_dir}/Aligned.filt.srt.nodup.tagAlign

BLACKLIST=/reference/public/ENCODE/hg38-blacklist.v2.bed

echo "Fragment length: $fraglen" >& 2
echo "Sample directory:" $sample_dir >& 2
echo "Input directory:" $input_dir >& 2 
echo "Sample bed:" $bed_chip >& 2
echo "Input bed:" $bed_input >& 2 
echo "Sample bam:" $bam_chip >& 2

# 1. Narrow peaks
prefix=${sample_dir}/peaks_q0.01/narrow_tagAlign
cmd="macs2 callpeak -t ${bed_chip}.gz -c ${bed_input}.gz -g 3.0e9 -n $prefix -q 1e-2 --nomodel --shift 0 --extsize $fraglen --keep-dup all -B --SPMR --call-summits"
echo $cmd >> $log; echo $cmd >& 2; eval $cmd

# Remove blacklist peaks
FILTERED_PEAK=${prefix}_peaks_noblacklist.narrowPeak
bedtools intersect -v -a ${prefix}_peaks.narrowPeak -b ${BLACKLIST} | grep -P 'chr[\dXY]+[ \t]'   > ${FILTERED_PEAK}

# Run featureCounts
cmd="awk -F \"\t\" '{print \$4,\$1,\$2,\$3,\$6}' ${prefix}_peaks_noblacklist.narrowPeak | tr ' ' '\\t' > ${sample_dir}/narrow_saf"
echo $cmd >> $log; echo $cmd >& 2; eval $cmd

cmd="featureCounts -p --countReadPairs -B -C -T 4 -F SAF -a ${sample_dir}/narrow_saf -o ${prefix}_peaks_noblacklist.narrowPeak.counts $bam_chip"
echo $cmd >> $log; echo $cmd >& 2; eval $cmd

# 2. Broad peaks
prefix=${sample_dir}/peaks_q0.01/broad_tagAlign
cmd="macs2 callpeak -t ${bed_chip}.gz -c ${bed_input}.gz -g 3.0e9 -n $prefix -q 1e-2 --nomodel --shift 0 --extsize $fraglen --keep-dup all -B --SPMR --broad"
echo $cmd >> $log; echo $cmd >& 2; eval $cmd

# Remove blacklist peaks (ended up not being used)
FILTERED_PEAK=${prefix}_peaks_noblacklist.broadPeak
bedtools intersect -v -a ${prefix}_peaks.broadPeak -b ${BLACKLIST} | grep -P 'chr[\dXY]+[ \t]'   > ${FILTERED_PEAK}

# Run featureCounts
cmd="awk -F \"\t\" '{print \$4,\$1,\$2,\$3,\$6}' ${prefix}_peaks_noblacklist.broadPeak | tr ' ' '\\t' > ${sample_dir}/broad_saf"
echo $cmd >> $log; echo $cmd >& 2; eval $cmd

cmd="featureCounts -p --countReadPairs -B -C -T 4 -F SAF -a ${sample_dir}/broad_saf -o ${prefix}_peaks_noblacklist.broadPeak.counts $bam_chip"
echo $cmd >> $log; echo $cmd >& 2; eval $cmd
