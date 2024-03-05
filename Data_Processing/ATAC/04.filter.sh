#!/bin/bash

source /frazer01/home/jennifer/.bash_profile
source activate encode-atac

OUT_DIR=$1
log=$2
PREFIX=${OUT_DIR}/Aligned.sorted
RAW_BAM_FILE=${PREFIX}.bam
FILT_BAM_PREFIX="${PREFIX}.filt"
FILT_BAM_FILE="${FILT_BAM_PREFIX}.bam"
TMP_FILT_BAM_PREFIX="${FILT_BAM_PREFIX}.nmsrt"
TMP_FILT_BAM_FILE="${TMP_FILT_BAM_PREFIX}.bam"
TMP_FILT_FIXMATE_BAM_FILE="${TMP_FILT_BAM_PREFIX}.fixmate.bam"
MULTIMAPPING=4

TMP_FILT_BAM_FILE="${FILT_BAM_PREFIX}.dupmark.bam"
MARKDUP="picard MarkDuplicates"
DUP_FILE_QC="${FILT_BAM_PREFIX}.dup.qc"
TMP_FILT_BAM_FILE_QC=${FILT_BAM_PREFIX}.dupmark.flagstat.qc

FINAL_BAM_PREFIX="${FILT_BAM_PREFIX}.nodup"
FINAL_BAM_FILE="${FINAL_BAM_PREFIX}.bam" # To be stored
FINAL_BAM_INDEX_FILE="${FINAL_BAM_FILE}.bai"
FINAL_BAM_FILE_MAPSTATS="${FINAL_BAM_PREFIX}.flagstat.qc" # QC file

# 1. Remove  unmapped, mate unmapped, not primary alignment, reads failing platform
cmd="samtools view -@ 8 -F 524 -f 2 -q 30 -u ${RAW_BAM_FILE} | samtools sort -@ 8 -n -o ${TMP_FILT_BAM_FILE} -"
echo $cmd >& 2; echo $cmd >> $log; eval $cmd

# 2. Take multimappers and randomly assign
ENCODE_DIR=/frazer01/home/jennifer/software/encode-atac-seq-pipeline/src
cmd="samtools view -h ${TMP_FILT_BAM_FILE} | ${ENCODE_DIR}/assign_multimappers.py -k $MULTIMAPPING --paired-end | samtools fixmate -@ 8 -r - ${TMP_FILT_FIXMATE_BAM_FILE}"
echo $cmd >& 2; echo $cmd >> $log; eval $cmd

# 3. Remove reads unmapped, mate unmapped, not primary alignment, failed platform QC, PCR duplicates. Keep reads in proper pair
cmd="samtools view -@ 8 -F 1804 -f 2 -q 30 -u ${TMP_FILT_FIXMATE_BAM_FILE} | samtools sort -@ 8 -o ${FILT_BAM_FILE} -"
echo $cmd >& 2; echo $cmd >> $log; eval $cmd

cmd="rm ${TMP_FILT_FIXMATE_BAM_FILE} ${TMP_FILT_BAM_FILE}"
echo $cmd >& 2; echo $cmd >> $log; eval $cmd

# 4. Mark duplicates
cmd="${MARKDUP} INPUT=${FILT_BAM_FILE} OUTPUT=${TMP_FILT_BAM_FILE} METRICS_FILE=${DUP_FILE_QC} VALIDATION_STRINGENCY=LENIENT ASSUME_SORTED=true REMOVE_DUPLICATES=false"
echo $cmd >& 2; echo $cmd >> $log; eval $cmd

cmd="samtools sort -@ 8 -n --threads 10 ${TMP_FILT_BAM_FILE} -O SAM  | SAMstats --sorted_sam_file -  --outf ${TMP_FILT_BAM_FILE_QC}"
echo $cmd >& 2; echo $cmd >> $log; eval $cmd

cmd="mv ${TMP_FILT_BAM_FILE} ${FILT_BAM_FILE}"
echo $cmd >& 2; echo $cmd >> $log; eval $cmd

cmd="mv ${DUP_FILE_QC} ${OUT_DIR}/qc"
echo $cmd >& 2; echo $cmd >> $log; eval $cmd

# 5. Remove duplicates
cmd="samtools view -@ 8 -F 1804 -f 2 -b ${FILT_BAM_FILE} > ${FINAL_BAM_FILE}"
echo $cmd >& 2; echo $cmd >> $log; eval $cmd

# Index Final BAM file
cmd="samtools index -@ 8 ${FINAL_BAM_FILE}"
echo $cmd >& 2; echo $cmd >> $log; eval $cmd

cmd="samtools sort -@ 8 -n --threads 10 ${FINAL_BAM_FILE} -O SAM  | SAMstats --sorted_sam_file -  --outf ${FINAL_BAM_FILE_MAPSTATS}"
echo $cmd >& 2; echo $cmd >> $log; eval $cmd

cmd="mv $FINAL_BAM_FILE_MAPSTATS ${OUT_DIR}/qc"
echo $cmd >& 2; echo $cmd >> $log; eval $cmd

cmd="rm ${FILT_BAM_FILE}"
echo $cmd >& 2; echo $cmd >> $log; eval $cmd

# 6. Remove mitochondrial chromosomes
nomito=${FILT_BAM_PREFIX}.nodup.nomito.bam
cmd="samtools idxstats ${FINAL_BAM_FILE} | cut -f 1 | grep -v -P "^chrM" | xargs samtools view ${FINAL_BAM_FILE} -@ 8 -b > ${nomito}"
echo $cmd >& 2; echo $cmd >> $log; eval $cmd

cmd="samtools index -@ 8 ${nomito} ${nomito}.bai"
echo $cmd >& 2; echo $cmd >> $log; eval $cmd


