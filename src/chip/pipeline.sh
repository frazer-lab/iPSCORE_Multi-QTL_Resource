#!/bin/bash

#$ -N chip
#$ -t 1-292:1
#$ -tc 3
#$ -pe smp 16
#$ -V -cwd
#$ -o logs
#$ -e logs

source /frazer01/home/jennifer/.bash_profile
source activate encode-chip
source /projects/CARDIPS/pipeline/ChIP-Seq/work/2023_0731_encode/encode_script/functions.sh

script_dir=/projects/CARDIPS/pipeline/ChIP-Seq/encode_script
pipe_dir=/projects/CARDIPS/pipeline/ChIP-Seq/encode_sample
data_dir=/projects/CARDIPS/data/ChIP-Seq/sample

if [ $# != 2 ]
then 
    
    id=`tail -n +$SGE_TASK_ID samples_uniq.txt | head -1`
    out_dir=${pipe_dir}/${id}
    fq_dir=${data_dir}/${id}
    call_peaks=$1
    input=$2
    
    echo "Sample ID: $id"
    echo "Input ID: $input"
    echo "Calling peaks: $call_peaks"
    echo "Output directory: $out_dir"
    echo "Fastq directory: $fq_dir"

else
    echo "Error. Provide input: sh pipeline.sh <TRUE/FALSE-call_peaks> <input_id>"
    exit
fi

mkdir ${out_dir} ${out_dir}/detect_adapters ${out_dir}/qc

if [ $call_peaks == F ]
then

    echo "Trim and qc.."
    if [ ! -f ${out_dir}/detect_adapters/reads.1.fastq ] && [ ! -f ${out_dir}/detect_adapters/reads.2.fastq ]
    then
        cmd="sh ${script_dir}/pre-align.sh $id $fq_dir $out_dir"
        echo $cmd >& 2; eval $cmd
    else
        echo Pre-align already done. Skipping.
    fi

    echo "Align.."
    cmd="sh ${script_dir}/align.sh $id $out_dir"
    echo $cmd >& 2; eval $cmd

    echo "Filter reads.."
    cmd="sh ${script_dir}/filter.sh $id $out_dir"
    echo $cmd >& 2; eval $cmd
else
    
    echo "Call peaks.."
    cmd="sh ${script_dir}/peaks.sh $id $input"
    echo $cmd >& 2; eval $cmd
fi
