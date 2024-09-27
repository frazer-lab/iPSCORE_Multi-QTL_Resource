#!/bin/bash

# Author: Jennifer Nguyen

#$ -N gwas
#$ -V -cwd
#$ -o logs
#$ -e logs

#source /frazer01/home/jennifer/.bash_profile
#export PATH=/home/tarthur/software/R-4.2.1/bin:$PATH

input_file=$1 # contains columns type, element_id, egene (T/F), qtl_id, tissue, analysis (eqtls/caqtls/haqtls), taskid
manifest=$2
outdir=$3

Rscript 05.1c.run_qtl_gwas_coloc.R --taskid $SGE_TASK_ID --input_file $input_file --manifest $manifest --outdir $outdir
