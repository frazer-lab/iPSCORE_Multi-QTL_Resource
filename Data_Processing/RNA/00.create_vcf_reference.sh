# Get VCF input for sample identity
vcf=/reference/private/CARDIPS/CARDIPS.GATK.autosome.common_snp.vcf.gz
ref_fa=/reference/public/UCSC/hg38/hg38.fa
chain=/reference/public/UCSC/hg38/hg19ToHg38.over.chain.gz
out_prefix=plink/cardips_common_hg38
region_bed=/reference/private/Gencode.v44lift38/exon.bed

source activate frazer-rna

# 1. Retain common variants (MAF between 45-55%)
cmd="vcftools --gzvcf $vcf --maf 0.45 --max-maf 0.55 --recode --bed --out ${out_dir}/filtered"
echo $cmd >& 2; eval $cmd

# 2. Lift to hg38
cmd="CrossMap.py vcf ${chain} ${out_dir}/filtered.recode.vcf ${ref_fa} ${out_dir}/filtered_hg38.vcf"
echo $cmd >& 2; eval $cmd

# 3. Clean
cmd="rm ${out_dir}/filtered.recode.vcf"
echo $cmd >& 2; eval $cmd

# 4. Select only variants in exon regions (for RNA only)
cmd="vcftools --vcf ${out_dir}/filtered_hg38.vcf --bed ${region_bed} --recode --out ${out_dir}/filtered_region"
echo $cmd >& 2; eval $cmd

# 5. Sort VCF
cmd="bcftools sort -Oz -o ${out_dir}/filtered_region.vcf.gz ${out_dir}/filtered_region.recode.vcf"
echo $cmd >& 2; eval $cmd

# 6. Index VCF
cmd="tabix -p vcf ${out_dir}/filtered_region.vcf.gz"
echo $cmd >& 2; eval $cmd

# 7. Clean
cmd="rm ${out_dir}/filtered_region.recode.vcf"
echo $cmd >& 2; eval $cmd
