#!/bin/bash
#$ -cwd
#$ -V 
#$ -t 1-22  # By chromosome

# Variables
plink_dir="path_to_1000G_plink_files"
outdir="path_to_out_directory"
chr_num="chr${SGE_TASK_ID}"

# Step 1: Filter the VCF for the specific chromosome using the SGE_TASK_ID
plink --bfile ${plink_dir}/${chr_num}.renamed --extract ${outdir}/${chr_num}.snps --make-bed --out ${outdir}/${chr_num}

# Step 2: Calculate linkage disequilibrium (LD) for the specific chromosome
plink --bfile ${outdir}/${chr_num} --r2 --ld-window-r2 0 --ld-window 99999 --ld-window-kb 10000 --out ${outdir}/${chr_num}_ld

echo "LD calculation is complete for ${chr_num}. Results are saved in ${outdir}/${chr_num}_ld.ld"

