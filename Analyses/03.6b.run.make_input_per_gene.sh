#!/bin/bash

# Author: Jennifer Nguyen

#$ -N make_input
#$ -V -cwd
#$ -pe smp 2

genelist=$1
date=$2

gene=`tail -n +$SGE_TASK_ID $genelist | head -1`

Rscript /projects/CARDIPS/analysis/epigenome_resource/analyses/jennifer/ipscore_unique_qtls/mashr/scripts/make_input_per_gene.R --gene_id $gene --date $date