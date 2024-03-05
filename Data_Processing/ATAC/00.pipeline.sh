#!/bin/bash

#$ -pe smp 8
#$ -V -cwd

source /frazer01/home/jennifer/.bash_profile
source activate encode-atac

encode_dir=/frazer01/home/jennifer/software/encode-atac-seq-pipeline
script_dir=/projects/PPC/pipeline/ATAC-Seq/work/2023_0720/scripts

fq_dir=$1
merge_uuids=$2
merge_uuid=`tail -n +$SGE_TASK_ID $merge_uuids | head -1`
map_file=$3
out_dir=$4
out_dir=${out_dir}/${merge_uuid}
ref_saf=$5

mkdir $out_dir $out_dir/detect_adapters $out_dir/logs $out_dir/qc

log=${out_dir}/pipeline.log

# 1. Trim adapters 
# Note: ENCODE uses the same adapter seq for both fq1 and fq1; see "detect_most_likely_adapter" function; https://github.com/ENCODE-DCC/atac-seq-pipeline/blob/master/src/encode_task_trim_adapter.py
sh ${script_dir}/cutadapt.sh $fq_dir $merge_uuid $map_file $out_dir $log

# 1.1. Index reference genome for alignment (already done! do not run again. takes a very long time.)
cmd="source activate encode-atac; bowtie2-build ~/references/Gencode.v34lift37/GRCh37.primary_assembly.genome.fa.gz /reference/private/Gencode.v34lift37/bowtie2_index"
echo $cmd | qsub -N index -V -cwd -o logs/index.out -e logs/index.err -pe smp 8

# 2. Align
sh $script_dir/align.sh $out_dir $log

# 3. Calculate % mito reads
sh $script_dir/calculate_mitochondrial_reads.sh $out_dir $log

# 4. Filter reads
sh $script_dir/filter.sh $out_dir $log

# 5. Sample id
log_out=${out_dir}/logs/id.out
log_err=${out_dir}/logs/id.err
qsub -N id -o $log_out -e $log_err -V -cwd ${script_dir}/identity.sh $out_dir $log

# 6. Remove mitochondrial and duplicates
if [ ! -f ${out_dir}/Aligned.sorted.filt.nodup.nomito.bam ] || [ ! -f ${out_dir}/Aligned.sorted.filt.nodup.nomito.bam.bai ]
then
    log_out=${out_dir}/logs/mito_nondup.out
    log_err=${out_dir}/logs/mito_nondup.err
    qsub -hold_jid qc -N mito -o $log_out -e $log_err -V -cwd -pe smp 8 ${script_dir}/remove_mitochondrial_dup.sh $out_dir $log
fi

# 7. QC metrics (GC bias and fragment length)
if [ ! -f ${out_dir}/qc/Aligned.sorted.filt.nodup.nomito.pbc.qc ]
then
    log_out=${out_dir}/logs/qc.out
    log_err=${out_dir}/logs/qc.err
    qsub -N qc -o $log_out -e $log_err -V -cwd $script_dir/get_qc_metrics.sh $out_dir $log
fi

# 8. Peaks
if [ ! -f ${out_dir}/peaks/narrow_tn5_tagAlign_peaks_noblacklist.narrowPeak.counts ]
then
    log_out=${out_dir}/logs/peaks.out
    log_err=${out_dir}/logs/peaks.err
    qsub -hold_jid qc -N peaks -o $log_out -e $log_err -V -cwd ${script_dir}/peaks.sh $out_dir $log
fi

# 9. Tsse
if [ ! -f ${out_dir}/qc/TSSEscore.robj ]
then
    cmd="Rscript ${script_dir}/tsse_atacqc.R --pipe_dir $out_dir"
    echo $cmd >& 2; echo $cmd >> $log; eval $cmd
fi

# 10. Library Complexity
if [ ! -f ${out_dir}/qc/library_complexity.txt ]
then
    log_out=${out_dir}/logs/library_complexity.out
    log_err=${out_dir}/logs/library_complexity.err
    qsub -hold_jid qc -N peaks -o $log_out -e $log_err -V -cwd ${script_dir}/library_complexity.sh $out_dir $log
fi

# 11. Get counts for reference peaks
# Only run after identifying reference peaks
if [ ! -f ${out_dir}/ref_peaks/ref_peaks.counts ]
then
   if [ ! -d ${out_dir}/ref_peaks ]; then mkdir ${out_dir}/ref_peaks; fi
   
   # ALWAYS keep the log from featureCounts to get FRIP
   log_out=${out_dir}/ref_peaks/log_counts.out
   log_err=${out_dir}/ref_peaks/log_counts.err 
   
   qsub -N counts -o $log_out -e $log_err -V -cwd -pe smp 4 ${script_dir}/counts_for_reference_peaks.sh $out_dir $ref_saf
fi
