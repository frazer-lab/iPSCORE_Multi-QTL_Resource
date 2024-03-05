#!/bin/bash

#$ -N peaks
#$ -V -cwd
#$ -o logs/cc
#$ -e logs/cc
#$ -pe smp 8

source /frazer01/home/jennifer/.bash_profile
source activate encode-chip

out_dir=$1

dir=`mktemp -d -p scratch`

cmd="samtools view -@ 12 -F 1804 -f 2 -q 30 -o ${dir}/Aligned.filt.bam ${out_dir}/Aligned.bam"
echo $cmd >& 2; eval $cmd

# =================================
# make tagAlign for filtered (but not deduped) BAM
# and subsample it for cross-correlation analysis 
# ================================
NREADS=10000000
TA_FILE=${dir}/Aligned.merged.tagAlign
SUBSAMPLED_TA_FILE=${dir}/Aligned.merged.subsampled.tagAlign

bedtools bamtobed -i ${dir}/Aligned.filt.bam | awk 'BEGIN{OFS="\t"}{$4="N";$5="1000";print $0}' > ${TA_FILE}

grep -v “chrM” ${TA_FILE} | shuf -n ${NREADS} > ${SUBSAMPLED_TA_FILE}

# Estimate read length from first 100 reads.
READ_LEN=$(head -n 100 ${TA_FILE} | awk 'function abs(v) {{return v < 0 ? -v : v}} BEGIN{{sum=0}} {{sum+=abs($3-$2)}} END{{print int(sum/NR)}}')

# Determine exclusion range for fragment length estimation.
# Use a fixed lowerbound at -500.
# Upperbound EXCLUSION_RANGE_MAX is 
#   TF ChIP-seq:  max(read_len + 10, 50)
#   Histone ChIP-seq:  max(read_len + 10, 100)

# lowerbound is fixed at 500 for both
EXCLUSION_RANGE_MIN=-500
EXCLUSION_RANGE_MAX=$(python -c "read_len = $READ_LEN; print(max(read_len + 10, 50))")
#rm -f ${TA_FILE}.tmp

### cross-correlation analysis
NTHREADS=8
CC_SCORES_FILE="${SUBSAMPLED_TA_FILE}.cc.qc"
CC_PLOT_FILE="${SUBSAMPLED_TA_FILE}.cc.plot.pdf"

# CC_SCORE FILE format
# Filename <tab> numReads <tab> estFragLen <tab> corr_estFragLen <tab> PhantomPeak <tab> corr_phantomPeak <tab> argmin_corr <tab> min_corr <tab> phantomPeakCoef <tab> relPhantomPeakCoef <tab> QualityTag

source /frazer01/home/jennifer/.bash_profile

cmd="Rscript /frazer01/home/jennifer/software/phantompeakqualtools/run_spp.R -c=${SUBSAMPLED_TA_FILE} -p=${NTHREADS} -filtchr=chrM -savp=${CC_PLOT_FILE} -out=${CC_SCORES_FILE} -x=${EXCLUSION_RANGE_MIN}:${EXCLUSION_RANGE_MAX} -rf"
echo $cmd >& 2; eval $cmd

sed -r 's/,[^\t]+//g' ${CC_SCORES_FILE} > temp
mv temp ${CC_SCORES_FILE}

cmd="mv ${CC_PLOT_FILE} ${CC_SCORES_FILE} ${out_dir}"
echo $cmd >& 2; eval $cmd

cmd="rm -r ${dir}"
echo $cmd >& 2; eval $cmd
