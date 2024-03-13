#$ -V
#$ -cwd
#$ -e logs
#$ -o logs
#$ -pe smp 4

hostname >& 2

source /frazer01/home/jennifer/.bash_profile
source activate encode-atac

out_dir=$1
log=$2
plink_dir=${out_dir}/plink
id=`basename $out_dir`
bam=${out_dir}/Aligned.sorted.filt.nodup.bam
reference_vcf=/projects/PPC/pipeline/RNA-Seq/work/2023_0804/plink/cardips_common_hg38.vcf.gz 
reference=/reference/public/ucsc/hg38/hg38.fa

if [ ! -d ${plink_dir} ]; then mkdir ${plink_dir}; fi

if [ ! -f ${bam}.bai ]; then samtools index ${bam}; fi


# 1. Prepare file to rename the sample as "target"
echo "${bam} target" > ${plink_dir}/sample.txt

# 2. Call genotypes
cmd="bcftools mpileup --threads 4 -Ou -f $reference -R $reference_vcf $bam |\
bcftools call --threads 4 -Ou -mv |\
bcftools filter -e 'DP<10'  | bcftools reheader -s ${plink_dir}/sample.txt |\
awk 'BEGIN{OFS=\"\t\";} {if(\$1 !~ /^#/) {\$3=\$1\":\"\$2;} print;}' |\
bgzip -c > ${plink_dir}/call.vcf.gz"
date >> $log; echo $cmd >& 2; echo $cmd >> $log; eval $cmd

# 3. Index VCF
cmd="tabix -p vcf ${plink_dir}/call.vcf.gz"
date >> $log; echo $cmd >& 2; echo $cmd >> $log; eval $cmd

# 4. Merge with reference VCF
cmd="bcftools merge ${plink_dir}/call.vcf.gz $reference_vcf |\
bcftools view --threads 4 -m2 -M2 -v snps -o ${plink_dir}/merged.vcf.gz -O z"
date >> $log; echo $cmd >& 2; echo $cmd >> $log; eval $cmd

# 5. Convert VCF to plink bed files
cmd="plink --threads 4 --vcf ${plink_dir}/merged.vcf.gz --make-bed --out ${plink_dir}/merged --allow-extra-chr"
date >> $log; echo $cmd >& 2; echo $cmd >> $log; eval $cmd

# 6. Run plink genome
cmd="plink --threads 4 --bfile ${plink_dir}/merged --genome full --out ${plink_dir}/plink --allow-extra-chr"
date >> $log; echo $cmd >& 2; echo $cmd >> $log; eval $cmd

# 7. Clean
cmd="rm ${plink_dir}/call.vcf.gz ${plink_dir}/call.vcf.gz.tbi ${plink_dir}/sample.txt"
date >> $log; echo $cmd >& 2; echo $cmd >> $log; eval $cmd
