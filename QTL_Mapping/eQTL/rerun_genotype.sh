#!/usr/bin/sh

#$ -V -cwd 
#$ -pe smp 2
#$ -t 1-1130:1
#$ -tc 200
#$ -o logs/redo_genotypes.out
#$ -e logs/redo_genotypes.err

source /frazer01/home/jennifer/.bash_profile

SGE_TASK_ID=`tail -n +$SGE_TASK_ID /projects/CARDIPS/analysis/epigenome_resource/eqtls/CVPC/notebooks/geno_taskids.txt | head -1`

Rscript /projects/CARDIPS/analysis/epigenome_resource/eqtls/CVPC/notebooks/prepare_genotype.R --config_file /projects/CARDIPS/analysis/epigenome_resource/eqtls/CVPC/notebooks/qtl.config.sh --taskid $SGE_TASK_ID --functions_file /projects/CARDIPS/analysis/epigenome_resource/eqtls/CVPC/notebooks/functions.R

