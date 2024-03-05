#!/bin/bash/

#$ -N chip
#$ -t 292-292:1
#$ -tc 5
#$ -pe smp 12
#$ -V -cwd
#$ -o logs
#$ -e logs

source /frazer01/home/jennifer/.bash_profile
source activate encode-chip

script_dir=/projects/CARDIPS/pipeline/ChIP-Seq/encode_script
pipe_dir=/projects/CARDIPS/pipeline/ChIP-Seq/sample_hg38
data_dir=/projects/CARDIPS/data/ChIP-Seq/sample

sample_list=$1 # List of samples
call_peaks=$2 # True /False to call peaks
sample_input_list=$3 # first column = input id, second column = chip id
cell=$4 # ipsc, input, or cvpc

if [ $call_peaks == F ]
then
    id=`tail -n +$SGE_TASK_ID $sample_list | head -1`
    out_dir=${pipe_dir}/${id}
    fq_dir=${data_dir}/${id}
    input="skipping input"
else
    id=`tail -n +$SGE_TASK_ID $sample_list | head -1`
    input=`grep $id $sample_input_list | cut -f1`
fi

echo "Sample ID: $id" >& 2
echo "Input ID: $input" >& 2
echo "Calling peaks: $call_peaks" >& 2
echo "Output directory: $out_dir" >& 2
echo "Fastq directory: $fq_dir" >& 2

# get genome reference
# (already done, do not re-run)
# wget https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz
# wget https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.chrom.sizes
# bwa=/software/bwa-0.7.17/bwa
# ${bwa} index hg38.fa.gz                    

if [ $call_peaks == F ]
then

    mkdir ${out_dir} ${out_dir}/detect_adapters ${out_dir}/qc
    
    # 1. Concatenate fastqs. Trim adapters. Run fastqc. 
    cmd="sh ${script_dir}/pre-align.sh $id $fq_dir $out_dir"
    echo $cmd >& 2; eval $cmd

    # 2. Align
    cmd="sh ${script_dir}/align.sh $id $out_dir"
    echo $cmd >& 2; eval $cmd

    # 3. Filter 
    cmd="sh ${script_dir}/filter.sh $id $out_dir"
    echo $cmd >& 2; eval $cmd
    
    # 4. Create tagAlign
    cmd="sh ${script_dir}/tagAlign.sh ${out_dir}/Aligned.filt.srt.nodup.bam"
    echo $cmd >& 2; eval $cmd
    
    # 5. Run cross-correlation
    cmd="sh ${script_dir}/cross_correlation.sh sample_hg38/${id}"
    echo $cmd >& 2; eval $cmd

    # 6. Run plink
    cmd="sh ${script_dir}/identity.sh sample_hg38/${id}"
    echo $cmd >& 2; eval $cmd

    # 7. Run Library Complexity
    cmd="sh ${script_dir}/library_complexity.sh sample_hg38/${id}"
    echo $cmd >& 2; eval $cmd
    
else

    # 8. Call peaks
    sample_dir=${pipe_dir}/${id}
    input_dir=${pipe_dir}/${input}
    
    if [ ! -d ${sample_dir}/logs ]; then mkdir ${sample_dir}/logs; fi
    
    log_out=${sample_dir}/logs/peaks_counts.out
    log_err=${sample_dir}/logs/peaks_counts.err
    
    cmd="qsub -N peaks -V -cwd -pe smp 4 -o $log_out -e $log_err ${script_dir}/peaks.sh ${sample_dir} ${input_dir}"
    echo $cmd >& 2; eval $cmd
fi
