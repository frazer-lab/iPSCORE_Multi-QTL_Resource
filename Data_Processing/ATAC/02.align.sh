#!/bin/bash

source /frazer01/home/jennifer/.bash_profile
source activate encode-atac

ENCODE_DIR=/frazer01/home/jennifer/software/encode-atac-seq-pipeline

OUT_DIR=$1
log=$2
FASTQ1=${OUT_DIR}/detect_adapters/reads.1.fastq
FASTQ2=${OUT_DIR}/detect_adapters/reads.2.fastq
BAM_FILE=${OUT_DIR}/Aligned.bam
BAM_SORTED_FILE=${OUT_DIR}/Aligned.sorted.bam
FLAGSTAT_FILE=${OUT_DIR}/qc/Aligned.sorted.flagstat.qc
ref=/reference/public/ucsc/hg38/hg38.fa.gz

#BWT_IDX=/reference/private/Gencode.v34lift37/bowtie2_index/bowtie2_index

# 1. Align with bowtie
date >& 2

#cmd="bowtie2 -k 4 -X 2000 --mm --threads 8 -x $BWT_IDX -1 $FASTQ1 -2 $FASTQ2 > ${SAM_FILE}"
#echo $cmd; eval $cmd
cmd="bwa mem -t 8 $ref $FASTQ1 $FASTQ2 | samtools sort -@ 12 -O bam -o $BAM_FILE -"
echo $cmd >& 2; echo $cmd >> $log; eval $cmd

# 2. Sort 
cmd="samtools sort -@ 8 -o ${BAM_SORTED_FILE} ${BAM_FILE}"
echo $cmd >& 2; echo $cmd >> $log; eval $cmd

date >& 2

cmd="rm $FASTQ1 $FASTQ2 ${BAM_FILE}"
echo $cmd >& 2; echo $cmd >> $log; eval $cmd

date >& 2

cmd="samtools sort -n -@ 8 ${BAM_SORTED_FILE} -O SAM | SAMstats --sorted_sam_file -  --outf ${FLAGSTAT_FILE}"
echo $cmd >& 2; echo $cmd >> $log; eval $cmd

date >& 2

cmd="samtools index -@ 8 $BAM_SORTED_FILE"
echo $cmd >& 2; echo $cmd >> $log; eval $cmd

date >& 2




