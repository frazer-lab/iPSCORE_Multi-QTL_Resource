##############################################################
# Configuration file for QTL calling
##############################################################

# General instructions
# TODO

##############################################################
# General parameters
# Required parameters:

out_folder=/projects/CARDIPS/analysis/epigenome_resource/eqtls/CVPC # output folder
script_dir=/projects/CARDIPS/analysis/epigenome_resource/eqtls/scripts # scripts folder
functions_file=/projects/CARDIPS/analysis/epigenome_resource/eqtls/scripts/functions.R # path for functions file

##############################################################
# General software
# Required software (in parenthesis is the corresponding variable name in the config file). If not provided, it is assumed that the software is specified in .bashrc:

bcftools= # if empty, it is assumed to be called as "bcftools"
python= # if empty, it is assumed to be called as "python"
qsub_queue= # if empty, qsub is run with no specific queue

##############################################################
# 1- Prepare input data
# Generate phenotype and genotype information, used as input for QTL analysis

config_file=/projects/CARDIPS/analysis/epigenome_resource/eqtls/CVPC/notebooks/qtl.config.sh

# Sample metadata: lookup table with M rows and at least 3 columns: phenotype_id (ID used in the phenotype table); genotype_id (ID in the header of the VCF files); and subject_id. All the other columns can contain covariates that will be normalized and used for QTL analysis
metadata_sample=/projects/CARDIPS/analysis/epigenome_resource/eqtls/CVPC/input/metadata_sample.txt

# (optional) Subject metadata: subject-associated metadata. If present, must contain at least four columns: genotype_id, subject_id, phenotype_id and at least one covariate. Can be used to provide genotype PCs
metadata_subject=/projects/CARDIPS/analysis/epigenome_resource/eqtls/CVPC/input/metadata_subject.txt

# Phenotype data, TMM or TPM: a matrix with N rows (number of elements to test: genes/peaks) and M columns (number of samples to test) - un-normalized data
input_phenotype_tpm_matrix=/projects/CARDIPS/analysis/epigenome_resource/eqtls/CVPC/input/phenotype_tpm_matrix.txt

# Phenotype data, read counts: a matrix with N rows (number of elements to test: genes/peaks) and M columns (number of samples to test) - un-normalized data
input_phenotype_count_matrix=/projects/CARDIPS/analysis/epigenome_resource/eqtls/CVPC/input/phenotype_count_matrix.txt

# Phenotype information: a BED file with N rows. Required are the following (first four columns): chromosome, start, end, element_id; where element_id corresponds to the elements in input_phenotype_matrix
input_phenotype_info=/projects/CARDIPS/analysis/epigenome_resource/eqtls/CVPC/input/phenotype_info.txt

# VCF files with genotype data, comma separated
# All vcfs need to be bgzip'd and have index files
input_vcfs=/projects/CARDIPS/analysis/epigenome_resource/eqtls/prepare_wgs/iPSCORE_hg38.norm.pass.MAF0.05.HWE.GT99.recode.vcf.gz

# SNP allele frequency file
# Has columns: chr, pos, af.1, af.2, allele.1, allele.2, id.1, id.2, where id.1 = chr_pos_allele.1_allele.2 and id.2 = chr_pos_allele.2_allele.1
freq_file=/projects/CARDIPS/analysis/epigenome_resource/eqtls/prepare_wgs/iPSCORE_hg38.norm.pass.MAF0.05.HWE.GT99.recode.frq

# max number of PEER factors to calculate: default = M-30
# original developer recommends 25% of individuals
# setting a higher number makes the model harder to converge
peer_factors_n_peer=50

# number of elements used to calculate PEER factors: default = 5000
peer_factors_n_elements=2000 

# Kinship matrix .rel file: row and column names are genotype_id. Generated using plink --bfile \[infile\] --make-rel square. It is assumed that a file kinship.rel.id exists in the same folder
kinship_file=/projects/CARDIPS/analysis/epigenome_resource/eqtls/prepare_wgs/kinship/kinship.rel 

# distance threshold to test QTLs: default = 2000 bp
qtl_distance=1000000 

# MAF threshold: minimum allele frequency to test: default = 5%
maf_threshold=0.05

# Chromosome sizes (chromsizes_file): obtained from http://hgdownload.cse.ucsc.edu/goldenpath/hg19/bigZips/hg19.chrom.sizes: default = download from UCSC (hg19)
chromsizes_file=/reference/public/ucsc/hg38/hg38.chrom.sizes

# Filter genes
expressed_pct = 0.2 # percentage of samples to have expressed genes
expressed_tpm = 0.01 # min TPM for gene to be expressed
expressed_counts = 6 # min read counts for gene to be expressed




