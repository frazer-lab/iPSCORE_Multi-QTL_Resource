#!/bin/bash

#$ -N prep_1kg
#$ -V -cwd
#$ -o logs
#$ -e logs
#$ -t 1-22:1
#$ -tc 22
#$ -pe smp 8

source /frazer01/home/jennifer/.bash_profile
source activate frazer-rna

pca_dir=/projects/CARDIPS/analysis/epigenome_resource/eqtls/prepare_wgs/pca
chr=$SGE_TASK_ID

ref_vcf=/reference/public/1000Genomes_hg38/ALL.chr${chr}.shapeit2_integrated_snvindels_v2a_27022019.GRCh38.phased.vcf.gz

# 1. Separate multi-allelic snps
cmd="bcftools norm --threads 8 -m - -Oz -o ${pca_dir}/input/chr${chr}.vcf.gz $ref_vcf"
echo $cmd; eval $cmd

# 2. Extract SNPs only, rename variants, output as PLINK bed files
cmd="plink2 --snps-only --vcf ${pca_dir}/input/chr${chr}.vcf.gz  --make-bed  --memory 3000 --threads 8 --out ${pca_dir}/input/chr${chr}.rename --set-all-var-ids @:#:\\\$1:\\\$2"
echo $cmd; eval $cmd

# 3. Write filename to merge later
cmd="echo ${pca_dir}/input/chr${chr}.rename >> ${pca_dir}/merge_list.txt"
echo $cmd; eval $cmd