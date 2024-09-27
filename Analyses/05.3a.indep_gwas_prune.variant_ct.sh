#!/bin/bash

# Author: Jennifer Nguyen

#$ -N ld
#$ -V -cwd
#$ -o /projects/CARDIPS/analysis/epigenome_resource/analyses/jennifer/logs
#$ -e /projects/CARDIPS/analysis/epigenome_resource/analyses/jennifer/logs
#$ -pe smp 8
#$ -t 1-22:1
#$ -tc 22

# original was 500kb r2 0.1

chr=$SGE_TASK_ID
kb=$1
step=$2
r2=$3
outdir=/projects/CARDIPS/analysis/epigenome_resource/analyses/jennifer/gwas_independent
manifest=/projects/CARDIPS/analysis/epigenome_resource/analyses/jennifer/gwas_independent/subset_manifest.txt
dir=`mktemp -d -p /projects/CARDIPS/analysis/epigenome_resource/analyses/jennifer/gwas_independent/scratch`

# filter vcf
#vcf=/reference/public/1000Genomes_hg38/ALL.chr${chr}.shapeit2_integrated_snvindels_v2a_27022019.GRCh38.phased.vcf.gz

# rename position (already done!)
#bcftools annotate --set-id '%CHROM\_%POS' $vcf > ${outdir}/reference/chr${chr}.vcf
#plink --memory 30000 --threads 8 --vcf ${outdir}/reference/chr${chr}.vcf --keep #/projects/CARDIPS/analysis/epigenome_resource/analyses/jennifer/gwas_independent/indep/EUR.txt --make-bed --out ${outdir}/reference/chr${chr}

tail -n +2 $manifest | cut -f5 | while read name;
do
    #name=Mahajan.NatGenet2018b.T2D.European_sorted
    file=${outdir}/indep/${name}.hg38_5e-08.txt

    # get trait description
    name=`basename $file`
    name=${name::-4}
    if [ ! -d ${outdir}/indep/${name} ]
    then
        mkdir ${outdir}/indep/${name}
    fi

    # print list of snps to prune
    
    # 1. using either chr_pos as id names
    # awk -F"\t" '{print $1"_"$2}' $file  | sort -u > ${dir}/${name}_${chr}_snps.txt
    # sed -i "s/chr//g" ${dir}/${name}_${chr}_snps.txt
    
    # or 2. using chr_pos_ref_alt as id names
    awk '{print $NF}' $file | tail -n +2 > ${dir}/${name}_${chr}_snps.txt

    # print number of snps to prune
    # nsnps=`grep -w "${chr}_" ${dir}/${name}_${chr}_snps.txt | wc -l`
    # echo "SNPs detected: ${nsnps}" 

    # run plink pruning
    cmd="plink \
    --memory 30000 \
    --threads 8 \
    --bfile ${outdir}/reference/chr${chr}.renamed \
    --extract ${dir}/${name}_${chr}_snps.txt \
    --indep-pairwise ${kb} ${step} ${r2} \
    --out ${outdir}/indep/${name}/chr${chr}.${kb}ct_${step}step_r${r2}"

    echo $cmd >& 2; eval $cmd

    # move to outdir
    rsync ${dir}/${name}_${chr}_snps.txt ${outdir}/indep/${name}/chr${chr}.snps

done

#cmd="rm -r ${dir}"
#echo $cmd >& 2; eval $cmd
