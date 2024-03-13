
source /frazer01/home/jennifer
source activate frazer-rna

outdir=/projects/CARDIPS/analysis/epigenome_resource/eqtls/prepare_wgs
log=${outdir}/vcf_kinship_pca.log

pca_dir=${outdir}/pca
mkdir ${pca_dir}/input

# 1. Extract independent variants
cmd="plink --threads 8 --bfile ${outdir}/kinship/to_kinship --extract ${outdir}/kinship/to_kinship.prune.in --make-bed --out ${outdir}/pca/input/ipscore.to_pca"
echo $cmd >> $log; echo $cmd; eval $cmd

# 2. Make PLINK bed files
cmd="plink2 --threads 8 --bfile ${outdir}/pca/input/ipscore.to_pca --make-bed --out ${outdir}/pca/input/ipscore.to_pca.renamed --set-all-var-ids @:#:\\\$1:\\\$2"
echo $cmd >> $log; echo $cmd; eval $cmd

# 3. Prepare reference VCF files (1000 Genomes)
qsub /projects/CARDIPS/analysis/epigenome_resource/analyses/jennifer/notebooks/prep_1000genomes.sh

# 4. Merge VCF files
cmd="plink --merge-list ${pca_dir}/merge_list.txt --make-bed --out ${pca_dir}/input/reference.to_pca --memory 30000 --threads 8"
echo $cmd >> $log; echo $cmd; eval $cmd

# 5. Clean
cmd="rm ${pca_dir}/merge_list.txt ${pca_dir}/input/*prune* ${pca_dir}/input/*.bim ${pca_dir}/input/*.bed ${pca_dir}/input/*.fam ${pca_dir}/input/*.nosex"
echo $cmd >> $log; echo $cmd; eval $cmd

# 6. Make input for "within" parameter (i.e., a list of all samples. 3 columns - FID, ID, Cluster name)
awk -F "\t" '{print 0,$1,$6}' /reference/public/1000Genomes_hg38/sample_info.tsv | tail -n +2 > ${pca_dir}/within.txt
project=ipscore
bcftools query -l ${outdir}/iPSCORE_hg38.norm.pass.MAF0.05.HWE.GT99.recode.vcf |  awk -F "\t" -v a="$project" '{print $1, $1, a}' >> ${pca_dir}/within.txt
awk '$3 == "EUR" || $3 == "AMR" || $3 == "AFR" || $3 == "EAS" || $3 == "SAS" || $3 == "ipscore"' ${pca_dir}/within.txt > ${pca_dir}/within_filt.txt

rm ${pca_dir}/within.txt

# 7. Intersect snps between 1000 genomes and ipscore
cut -f2 ${outdir}/pca/input/ipscore.to_pca.renamed.bim | sort -u > snps1.txt
cut -f2 ${outdir}/pca/input/reference.to_pca.bim | sort -u > snps2.txt
comm -12 snps1.txt snps2.txt > ${pca_dir}/input/intersect.snps
rm snps1.txt snps2.txt ${pca_dir}/merge_list.txt

# 4. Filter for those overlapping variants
cmd="plink --bfile ${pca_dir}/input/reference.to_pca --make-bed --extract ${pca_dir}/input/intersect.snps --out ${pca_dir}/input/reference.to_pca.filt"
echo $cmd >> $log; echo $cmd; eval $cmd

echo "${pca_dir}/input/reference.to_pca.filt" >> ${pca_dir}/merge_list.txt
echo "${pca_dir}/input/ipscore.to_pca.renamed" >> ${pca_dir}/merge_list.txt

# 5. Merge and run PCA on those SNPs
plink --merge-list ${pca_dir}/merge_list.txt --make-bed --keep ${outdir}/pca/within_filt.txt --extract ${pca_dir}/input/intersect.snps --pca-cluster-names AFR EUR AMR EAS SAS --pca --within ${outdir}/pca/within_filt.txt --out ${pca_dir}/pca




