
pipe_dir=$1
map_file=$2
fq_dir=$3
SGE_TASK_ID=$4
line=(`tail -n +2 $map_file | tail -n +$SGE_TASK_ID | head -1`)
npremerge=${line[4]}
merge=${line[3]}
#fq=`ls ${fq_dir}/${premerge} | head -1`
#fq=${fq_dir}/${premerge}/${fq}
#bam=${pipe_dir}/${premerge}/Aligned.bam

echo "Processing $merge"
echo "Premerge: ${premerge}"
echo "Merge: ${merge}"
echo "FASTQ: ${fq}"
echo "BAM: ${bam}"

dir=`mktemp -d -p scratch`

# 1. Merge bams
if [ $npremerge == 1 ]
then
    fq=${fq_dir}/${line[0]}/`ls ${fq_dir}/${line[0]} | head -1`

    cmd="samtools merge -@ 12 -o ${dir}/Aligned.merged.bam ${pipe_dir}/${line[0]}/Aligned.bam"
    echo $cmd; eval $cmd
    
elif [ $npremerge == 2 ]
then
    fq1=${fq_dir}/${line[0]}/`ls ${fq_dir}/${line[0]} | head -1`
    fq2=${fq_dir}/${line[0]}/`ls ${fq_dir}/${line[0]} | head -1`
    cmd="zcat $fq1 $fq2 > ${dir}/input.fq"
    echo $cmd; eval $cmd
    
    cmd="samtools merge -@ 12 -o ${dir}/Aligned.merged.bam ${pipe_dir}/${line[0]}/Aligned.bam ${pipe_dir}/${line[1]}/Aligned.bam"
    echo $cmd; eval $cmd
    
elif [ $npremerge == 3 ]
then
    fq1=${fq_dir}/${line[0]}/`ls ${fq_dir}/${line[0]} | head -1`
    fq2=${fq_dir}/${line[0]}/`ls ${fq_dir}/${line[1]} | head -1`
    fq3=${fq_dir}/${line[0]}/`ls ${fq_dir}/${line[2]} | head -1`
    cmd="zcat $fq1 $fq2 $fq3 > ${dir}/input.fq"
    echo $cmd; eval $cmd
    
    cmd="samtools merge -@ 12 -o ${dir}/Aligned.merged.bam ${pipe_dir}/${line[0]}/Aligned.bam ${pipe_dir}/${line[1]}/Aligned.bam ${pipe_dir}/${line[2]}/Aligned.bam"
    echo $cmd; eval $cmd
    
fi

cmd="samtools view -@ 12 -F 1804 -f 2 -q 30 -o ${dir}/Aligned.filt.bam ${dir}/Aligned.merged.bam"
echo $cmd; eval $cmd

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
NTHREADS=12
CC_SCORES_FILE="${SUBSAMPLED_TA_FILE}.cc.qc"
CC_PLOT_FILE="${SUBSAMPLED_TA_FILE}.cc.plot.pdf"

# CC_SCORE FILE format
# Filename <tab> numReads <tab> estFragLen <tab> corr_estFragLen <tab> PhantomPeak <tab> corr_phantomPeak <tab> argmin_corr <tab> min_corr <tab> phantomPeakCoef <tab> relPhantomPeakCoef <tab> QualityTag

cmd="Rscript /frazer01/home/jennifer/software/phantompeakqualtools/run_spp.R -c=${SUBSAMPLED_TA_FILE} -p=${NTHREADS} -filtchr=chrM -savp=${CC_PLOT_FILE} -out=${CC_SCORES_FILE} -x=${EXCLUSION_RANGE_MIN}:${EXCLUSION_RANGE_MAX} -rf"
echo $cmd; eval $cmd

sed -r 's/,[^\t]+//g' ${CC_SCORES_FILE} > temp
mv temp ${CC_SCORES_FILE}

cmd="mv ${CC_PLOT_FILE} ${CC_SCORES_FILE} ${pipe_dir}/${merge}"
echo $cmd; eval $cmd

cmd="rm -r ${dir}"
echo $cmd; eval $cmd

