#$ -N chip-identity
#$ -V
#$ -cwd
#$ -e logs/identity
#$ -o logs/identity
#$ -pe smp 8
#$ -t 1-6:1
#$ -tc 20

set -e

source /frazer01/home/jennifer/.bash_profile
source activate frazer-rna

hostname >& 2

reference_vcf=/projects/PPC/pipeline/RNA-Seq/work/2023_0804/plink/cardips_common_hg38.vcf.gz 
reference=/reference/public/ucsc/hg38/hg38.fa
out_dir=$1 # output directory to save
#bam_list=$2 # text file containing list of bams
#bam=`tail -n +$SGE_TASK_ID $bam_list | head -1` # get bam file to test
#uuid=$(echo "$bam" | cut -d '/' -f 2)

log=${out_dir}/pipeline.log

if [ ! -f ${out_dir}/Aligned.merged.bam ]
then
    bam=${out_dir}/Aligned.filt.srt.nodup.bam # for samples that were sequenced once 
else
    bam=${out_dir}/Aligned.merged.bam # for samples that were sequenced multiple times - needed to merge bam files
fi

if [ -f ${out_dir}/plink.genome ]
then
    echo "${out_dir}/plink.genome already exists" >& 2
    exit 1
fi

if [ ! -f ${bam} ]
then
    echo "Error. ${bam} does not exist." >& 2
    exit 1
fi

if [ ! -d scratch ]
then
    cmd="mkdir scratch"
    echo $cmd >& 2; eval $cmd
fi

dir=`mktemp -d -p scratch` # make scratch directory

date >& 2

echo "$bam target" > ${dir}/sample.txt # make file to rename bamfile with "target" in VCF header

cat ${dir}/sample.txt >& 2

cmd="bcftools mpileup -Ou -f $reference -R $reference_vcf --threads 12 $bam |\
bcftools call --threads 12 -Ou -mv |\
bcftools filter -e 'DP<10'  | bcftools reheader -s ${dir}/sample.txt |\
awk 'BEGIN{OFS=\"\t\";} {if(\$1 !~ /^#/) {\$3=\$1\":\"\$2;} print;}' |\
bgzip -c > ${dir}/call.vcf.gz"
echo $cmd >> $log; echo $cmd >& 2; eval $cmd

date >& 2

cmd="tabix -p vcf ${dir}/call.vcf.gz"
echo $cmd >> $log; echo $cmd >& 2; eval $cmd

date >& 2

cmd="bcftools merge ${dir}/call.vcf.gz $reference_vcf |\
bcftools view --threads 8 -m2 -M2 -v snps -o ${dir}/merged.vcf.gz -O z"
echo $cmd >> $log; echo $cmd >& 2; eval $cmd

date >& 2

cmd="plink --threads 12 --vcf ${dir}/merged.vcf.gz --make-bed --out ${dir}/merged --allow-extra-chr"
echo $cmd >> $log; echo $cmd >& 2; eval $cmd

date >& 2

cmd="plink --threads 12 --bfile ${dir}/merged --genome full --out ${dir}/plink --allow-extra-chr"
echo $cmd >> $log; echo $cmd >& 2; eval $cmd

date >& 2

cmd="rsync ${dir}/plink.genome ${out_dir}/."
echo $cmd >> $log; echo $cmd >& 2; eval $cmd

date >& 2

cmd="rm -rf $dir"
echo $cmd >> $log; echo $cmd >& 2; eval $cmd

date >& 2
