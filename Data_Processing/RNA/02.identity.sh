#$ -N identity
#$ -V
#$ -cwd
#$ -e logs/identity
#$ -o logs/identity
#$ -pe smp 8

set -e

source /frazer01/home/jennifer/.bash_profile
source activate frazer-rna

hostname >& 2

reference=/reference/public/ucsc/hg38/hg38.fa
reference_vcf=$1
in_dir=$2
out_dir=$3
uuid=$4

in_dir=${in_dir}/${uuid}
out_dir=${out_dir}/${uuid}
bam=${in_dir}/Aligned.out.sorted.mdup.bam

log=${out_dir}/pipeline.log

if [ ! -f ${bam} ]
then
    echo "Error. ${bam} does not exist." >& 2
    exit 1
fi

if [ ! -f ${bam}.bai ]
then
    cmd="samtools index -@ 8 ${bam} ${bam}.bai" 
    echo $cmd >& 2
    eval $cmd
fi

if [ ! -d scratch ]
then
    cmd="mkdir scratch"
    echo $cmd >& 2; eval $cmd
fi

if [ -f ${out_dir}/plink.genome ]
then
    echo "${out_dir}/plink.genome already exists!" >& 2
    exit 1
fi

dir=`mktemp -d -p scratch`

# 1. Update header for bcftools
cmd="samtools view -H ${bam} > ${dir}/header.sam; awk '{gsub(\"SM:\", \"SM:target\"); print}' ${dir}/header.sam > ${dir}/header.txt"
#cmd="java -jar /software/picard-2.20.1/picard.jar AddOrReplaceReadGroups I=$bam O=${dir}/out.bam RGSM=target RGLB=lib1 RGPL=illumina RGPU=unit1"
echo $cmd >> $log; echo $cmd >& 2; eval $cmd

cmd="samtools reheader ${dir}/header.txt ${bam} > ${dir}/out.bam"
echo $cmd >> $log; echo $cmd >& 2; eval $cmd

cmd="samtools index ${dir}/out.bam"
echo $cmd >> $log; echo $cmd >& 2; eval $cmd

cmd="bam=${dir}/out.bam"
echo $cmd >> $log; echo $cmd >& 2; eval $cmd

# 2. Call genotypes on bam
echo "$bam target" > ${dir}/sample.txt
cat ${dir}/sample.txt >& 2
cmd="bcftools mpileup -Ou -f $reference -R $reference_vcf --threads 8 $bam |\
bcftools call --threads 8 -Ou -mv |\
bcftools filter -e 'DP<10'  | bcftools reheader -s ${dir}/sample.txt |\
awk 'BEGIN{OFS=\"\t\";} {if(\$1 !~ /^#/) {\$3=\$1\":\"\$2;} print;}' |\
bgzip -c > ${dir}/call.vcf.gz"
echo $cmd >> $log; echo $cmd >& 2; eval $cmd

cmd="tabix -p vcf ${dir}/call.vcf.gz"
echo $cmd >> $log; echo $cmd >& 2; eval $cmd

# 3. Merge with WGS vcf
cmd="bcftools merge ${dir}/call.vcf.gz $reference_vcf |\
bcftools view --threads 8 -m2 -M2 -v snps -o ${dir}/merged.vcf.gz -O z"
echo $cmd >> $log; echo $cmd >& 2; eval $cmd

# 4. Convert VCF to plink bed files
cmd="plink --threads 8 --vcf ${dir}/merged.vcf.gz --make-bed --out ${dir}/merged --allow-extra-chr"
echo $cmd >> $log; echo $cmd >& 2; eval $cmd

# 5. Run plink genome
cmd="plink --threads 8 --bfile ${dir}/merged --genome full --out ${dir}/plink --allow-extra-chr"
echo $cmd >> $log; echo $cmd >& 2; eval $cmd

# 6. Move 
cmd="rsync ${dir}/plink.genome ${out_dir}/."
echo $cmd >> $log; echo $cmd >& 2; eval $cmd

# 7. Clean
cmd="rm -rf $dir"
echo $cmd >> $log; echo $cmd >& 2; eval $cmd

