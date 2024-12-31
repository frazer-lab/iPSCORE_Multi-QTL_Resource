#!/bin/bash

# Author: Jennifer Nguyen

export PATH=~/software/R-4.3.3/bin:$PATH

genelist=$1
date=$2

gene=`tail -n +$SGE_TASK_ID $genelist | head -1`

Rscript /projects/CARDIPS/analysis/epigenome_resource/analyses/jennifer/ipscore_unique_qtls/mashr/scripts/mashr_run.R --gene_id $gene --date $date