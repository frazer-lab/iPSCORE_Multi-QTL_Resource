#!/bin/bash

#$ -N qc
#$ -V -cwd
#$ -pe smp 8

source /frazer01/home/jennifer/.bash_profile
source activate encode-atac

# ATAC-seq pipeline
# Documentation: https://docs.google.com/document/d/1f0Cm4vRyDQDu0bMehHD7P7KOMxTOP-HiNoIvL1VcBt8/edit
# Github: https://github.com/ENCODE-DCC/atac-seq-pipeline/tree/master

# ===================
# Create tagAlign file
# ===================
OUT_DIR=$1
log=$2
PREFIX=${OUT_DIR}/Aligned.sorted
FINAL_BAM_PREFIX=${OUT_DIR}/Aligned.sorted.filt.nodup.nomito
FINAL_BAM_FILE=${FINAL_BAM_PREFIX}.bam
FINAL_BEDPE_FILE=${FINAL_BAM_PREFIX}.bedpe.gz
FINAL_TA_FILE=${FINAL_BAM_PREFIX}.tagAlign.gz
SCRIPT_DIR=/projects/PPC/pipeline/ATAC-Seq/work/2023_0720/scripts
PBC_FILE_QC=${FINAL_BAM_PREFIX}.pbc.qc

# =============================
# Compute library complexity
# =============================

# TotalReadPairs [tab] DistinctReadPairs [tab] OneReadPair [tab] TwoReadPairs [tab] NRF=Distinct/Total [tab] PBC1=OnePair/Distinct [tab] PBC2=OnePair/TwoPair

cmd="samtools sort -@ 8 -n ${FINAL_BAM_FILE} -o ${FINAL_BAM_PREFIX}.srt.tmp.bam"
echo $cmd >& 2; echo $cmd >> $log; eval $cmd

bedtools bamtobed -bedpe -i ${FINAL_BAM_PREFIX}.srt.tmp.bam | awk 'BEGIN{OFS="\t"}{print $1,$2,$4,$6,$9,$10}' | grep -v 'chrM' | sort | uniq -c | awk 'BEGIN{mt=0;m0=0;m1=0;m2=0} ($1==1){m1=m1+1} ($1==2){m2=m2+1} {m0=m0+1} {mt=mt+$1} END{printf "%d\t%d\t%d\t%d\t%f\t%f\t%f\n",mt,m0,m1,m2,m0/mt,m1/m0,m1/m2}' > ${PBC_FILE_QC}

cmd="mv ${PBC_FILE_QC} ${OUT_DIR}/qc"
echo $cmd >& 2; echo $cmd >> $log; eval $cmd

# =============================
# Make tagAlign file for downstream QC
# =============================

cmd="bedtools bamtobed -bedpe -mate1 -i ${FINAL_BAM_PREFIX}.srt.tmp.bam | gzip -nc > ${FINAL_BEDPE_FILE}"
echo $cmd >& 2; echo $cmd >> $log; eval $cmd

# Create final TA file
zcat ${FINAL_BEDPE_FILE} | awk 'BEGIN{OFS="\t"}{printf "%s\t%s\t%s\tN\t1000\t%s\n%s\t%s\t%s\tN\t1000\t%s\n",$1,$2,$3,$9,$4,$5,$6,$10}' | gzip -nc > ${FINAL_TA_FILE}

cmd="rm ${FINAL_BAM_PREFIX}.srt.tmp.bam"
echo $cmd >& 2; echo $cmd >> $log; eval $cmd

# =================================
# Tn5 shifting
# ================================

SHIFTED_TAG_FILE="${FINAL_BAM_PREFIX}.tn5.tagAlign.gz"

zcat -f ${FINAL_TA_FILE} | awk 'BEGIN {OFS = "\t"} { if ($6 == "+") {$2 = $2 + 4} else if ($6 == "-") {$3 = $3 - 5} if ($2 >= $3) { if ($6 == "+") {$2 = $3 - 1} else {$3 = $2 + 1} } print $0}' | gzip -nc > ${SHIFTED_TAG_FILE}

# =================================
# Fragment length statistics
# ================================
INSERT_DATA=${FINAL_BAM_PREFIX}.insertsize.metrics
INSERT_PLOT=${FINAL_BAM_PREFIX}.insertsize.metrics.plot

cmd="picard CollectInsertSizeMetrics \
INPUT=${FINAL_BAM_FILE} OUTPUT=${INSERT_DATA} H=${INSERT_PLOT} \
VERBOSITY=ERROR QUIET=TRUE \
USE_JDK_DEFLATER=TRUE USE_JDK_INFLATER=TRUE \
W=1000 STOP_AFTER=5000000"
echo $cmd >& 2; echo $cmd >> $log; eval $cmd

cmd="python3 ${SCRIPT_DIR}/plot_fragment_length.py $INSERT_DATA $FINAL_BAM_PREFIX"
echo $cmd >& 2; echo $cmd >> $log; eval $cmd

cmd="mv $INSERT_DATA ${OUT_DIR}/*.plot ${OUT_DIR}/*.png ${OUT_DIR}/*.qc ${OUT_DIR}/qc"
echo $cmd >& 2; echo $cmd >> $log; eval $cmd

#cmd="rm $FINAL_BEDPE_FILE"
#echo $cmd >& 2; echo $cmd >> $log; eval $cmd
