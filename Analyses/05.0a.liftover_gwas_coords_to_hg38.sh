#!/bin/bash

# Author: Jennifer Nguyen

#$ -N lift
#$ -V -cwd
#$ -pe smp 4
#$ -o logs
#$ -e logs

id=$1

work_dir=/projects/CARDIPS/analysis/epigenome_resource/analyses/jennifer/gwas_liftover
chain_file="/reference/public/ucsc/hg19ToHg38.over.chain"
dir=`mktemp -d -p ${work_dir}/scratch`
         
# Convert summary statistics to BED format (assuming CHR is column 1 and POS is column 2)
awk 'NR > 1 {print "chr"$2"\t"$3-1"\t"$3"\t"$0}' ${work_dir}/hg19_summary_statistics/${id}.hg19.tsv > ${dir}/tmp.bed

# Run liftOver
cmd="/software/ucsc.linux.x86_64.20151103/liftOver -bedPlus=3 -tab ${dir}/tmp.bed $chain_file ${work_dir}/hg38_summary_statistics/${id}.hg38.tsv ${dir}/unmapped.txt"
echo $cmd >& 2; eval $cmd

# Reformat statistics
cmd="Rscript /projects/CARDIPS/analysis/epigenome_resource/analyses/jennifer/notebooks/08.00.process_gwas_liftover.R --id $id"
echo $cmd >& 2; eval $cmd

cmd="rm -r ${dir}"
echo $cmd >& 2; eval $cmd

# Tabix
cmd="sh /projects/CARDIPS/analysis/epigenome_resource/analyses/jennifer/notebooks/08.00.tabix.sh ${work_dir}/hg38_summary_statistics/${id}.hg38.tsv ${work_dir}/hg38_summary_statistics"
echo $cmd >& 2; eval $cmd
