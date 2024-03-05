#!/bin/bash

source /frazer01/home/jennifer/.bash_profile
source activate encode-atac

OUT_DIR=$1
log=$2
BAM_FILE=${OUT_DIR}/Aligned.sorted.bam
NON_MITO_BAM=${OUT_DIR}/Aligned.sorted.non_mito.bam
NON_MITO_FLAGSTAT=${OUT_DIR}/qc/Aligned.sorted.non_mito.flagstat.qc
MITO_FLAGSTAT=${OUT_DIR}/qc/Aligned.sorted.mito.flagstat.qc
MITO_BAM=${OUT_DIR}/Aligned.sorted.mito.bam
NTHREADS=4

# 1. Make mito-free bam
date >& 2

cmd="samtools idxstats ${BAM_FILE} | cut -f 1 | grep -v -P "^chrM" | xargs samtools view ${BAM_FILE} -@ ${NTHREADS} -b> ${NON_MITO_BAM}"
echo $cmd >& 2; echo $cmd >> $log; eval $cmd

date >& 2

# 2. Flagstat on mito-free bam
cmd="samtools sort -n --threads 8 ${NON_MITO_BAM} -O SAM  | SAMstats --sorted_sam_file -  --outf ${NON_MITO_FLAGSTAT}"
echo $cmd >& 2; echo $cmd >> $log; eval $cmd

cmd="rm $NON_MITO_BAM"
echo $cmd >& 2; echo $cmd >> $log; eval $cmd

date >& 2

# 3. Make mito bam
cmd="samtools view -b ${BAM_FILE} chrM > ${MITO_BAM}"
echo $cmd >& 2; echo $cmd >> $log; eval $cmd

date >& 2

# 4. Flagstat on mito bam
cmd="samtools sort -n --threads 8 ${MITO_BAM} -O SAM  | SAMstats --sorted_sam_file -  --outf ${MITO_FLAGSTAT}"
echo $cmd >& 2; echo $cmd >> $log; eval $cmd

date >& 2

cmd="rm $MITO_BAM"
echo $cmd >& 2; echo $cmd >> $log; eval $cmd

date >& 2

