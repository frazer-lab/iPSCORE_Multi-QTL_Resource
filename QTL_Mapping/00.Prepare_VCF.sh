
outdir=/projects/CARDIPS/analysis/epigenome_resource/eqtls/prepare_wgs
vcf=${outdir}/iPSCORE_hg38.vcf.gz
log=${outdir}/vcf_kinship_pca.log
    
# 1. Separate multi-allelic snps
cmd="bcftools norm --threads 8 -m - -Oz -o ${outdir}/iPSCORE_hg38.norm.vcf.gz $vcf "
echo $cmd >> $log; echo $cmd; eval $cmd

# 2. Filter for snps that passed instrument QC
cmd="bcftools view -f PASS -o ${outdir}/iPSCORE_hg38.norm.pass.vcf.gz -Oz ${outdir}/iPSCORE_hg38.norm.vcf.gz"
echo $cmd >> $log; echo $cmd; eval $cmd

# 3. Clean
cmd="rm ${outdir}/iPSCORE_hg38.norm.vcf.gz"
echo $cmd >> $log; echo $cmd; eval $cmd

# 4. Filter with MAF > 5%, HWE > 1e-06, present in >99% samples
cmd="vcftools --gzvcf ${outdir}/iPSCORE_hg38.norm.pass.vcf.gz --maf 0.05 --hwe 0.000001 --max-missing 0.99 --recode --out ${outdir}/iPSCORE_hg38.norm.pass.MAF0.05.HWE.GT99"
echo $cmd >> $log; echo $cmd; eval $cmd

# ended up not using 0.01 MAF threshold because most of the rare variants only have two GTs. We want to have a good number of samples for all 3 GTs. 
#cmd="vcftools --gzvcf ${outdir}/iPSCORE_hg38.norm.pass.vcf.gz --maf 0.01 --hwe 0.000001 --max-missing 0.99 --recode --out ${outdir}/iPSCORE_hg38.norm.pass.MAF0.01.HWE.GT99"
# echo $cmd >> $log; echo $cmd; eval $cmd

# 5. Bgzip VCF
cmd="bgzip -c ${outdir}/iPSCORE_hg38.norm.pass.MAF0.01.HWE.GT99.recode.vcf > ${outdir}/iPSCORE_hg38.norm.pass.MAF0.05.HWE.GT99.recode.vcf.gz"
echo $cmd >> $log; echo $cmd; eval $cmd

# 6. Index VCF
cmd="tabix -p vcf ${outdir}/iPSCORE_hg38.norm.pass.MAF0.05.HWE.GT99.recode.vcf.gz"
echo $cmd >> $log; echo $cmd; eval $cmd

# 7. Re-name SNP IDs usig CHROM_POS_REF_ALT
cmd="bcftools annotate --threads 8 --set-id '%CHROM\_%POS\_%REF\_%ALT' ${outdir}/iPSCORE_hg38.norm.pass.MAF0.05.HWE.GT99.recode.vcf.gz -Oz -o ${outdir}/iPSCORE_hg38.norm.pass.MAF0.05.HWE.GT99.recode.renamed.vcf.gz"
echo $cmd >> $log; echo $cmd; eval $cmd
