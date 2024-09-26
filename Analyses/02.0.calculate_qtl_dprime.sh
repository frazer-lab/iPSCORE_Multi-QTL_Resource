#!/bin/bash
#$ -cwd
#$ -V 
#$ -t 1-23478

# Variables
plink_dir="/projects/CARDIPS/analysis/epigenome_resource/analyses/jennifer/gwas_independent/reference"
pair_file="/projects/CARDIPS/analysis/epigenome_resource/analyses/tim/ld_modules/scripts/conditional_pairs_v3.txt"
outdir="/projects/CARDIPS/analysis/epigenome_resource/analyses/tim/ld_modules/conditional_dprime_v3"

line=$(sed -n "${SGE_TASK_ID}p" $pair_file)
read var1 var2 element chrom tissue <<< "$line"

plink_file="${plink_dir}/chr${chrom}.renamed"
# Define output file name based on variant IDs
output_name="${outdir}/${tissue}_${element}_${var1}_${var2}"

echo $var1 > ${output_name}.snps
echo $var2 >> ${output_name}.snps

plink --bfile ${plink_file} --extract ${output_name}.snps --make-bed --memory 30000 --threads 8 --out ${output_name}_filtered

plink \
  --bfile ${output_name}_filtered \
  --ld $var1 $var2 \
  --r2 \
  --ld-window-r2 0 \
  --memory 30000 --threads 8 \
  --out $output_name

grep R-sq ${output_name}.log > ${output_name}_dprime.txt
rm ${output_name}_filtered* ${output_name}.*
