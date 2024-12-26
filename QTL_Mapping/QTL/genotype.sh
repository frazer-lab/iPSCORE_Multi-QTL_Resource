#!/usr/bin/sh

source /frazer01/home/jennifer/.bash_profile

Rscript /projects/CARDIPS/analysis/epigenome_resource/eqtls/CVPC/notebooks/prepare_genotype.R --config_file /projects/CARDIPS/analysis/epigenome_resource/eqtls/CVPC/notebooks/qtl.config.sh --taskid $SGE_TASK_ID --functions_file /projects/CARDIPS/analysis/epigenome_resource/eqtls/CVPC/notebooks/functions.R

