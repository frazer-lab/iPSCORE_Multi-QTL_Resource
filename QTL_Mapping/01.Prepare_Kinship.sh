#!/bin/bash/

# From Chris's paper: We constructed an empirical kinship matrix for all subjects with WGS variant calls by intersecting biallelic SNVs with 1000 Genomes phase 3 variants and LD pruning the resulting variants using plink 1.90b3x (–biallelic-only –indep-pairwise 50 5 0.2) for unrelated EUR 1000 Genomes subjects (1000 Genomes Project Consortium, 2015; Chang et al., 2015). We used the remaining LD-pruned variants to construct the kinship matrix using EPACTS 3.2.6 (epacts make-kin –min-maf 0.01 –min-callrate 0.95) keeping variants whose frequency was above 1% in our cohort and that were called in at least 95% of our samples.

### Kinship
source /frazer01/home/jennifer
source activate frazer-rna

outdir=/projects/CARDIPS/analysis/epigenome_resource/eqtls/prepare_wgs
log=${outdir}/vcf_kinship_pca.log

cmd="plink --indep-pairwise 50 5 0.2 --vcf ${outdir}/iPSCORE_hg38.norm.pass.MAF0.05.HWE.GT99.recode.renamed.vcf.gz --make-bed --out ${outdir}/kinship/to_kinship"
echo $cmd >> $log; echo $cmd; eval $cmd

cmd="plink --bfile ${outdir}/kinship/to_kinship --make-rel square --extract ${outdir}/kinship/to_kinship.prune.in --out ${outdir}/kinship/kinship"
echo $cmd >> $log; echo $cmd; eval $cmd

