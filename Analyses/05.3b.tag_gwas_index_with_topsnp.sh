#!/bin/bash

# Author: Jennifer Nguyen

#$ -N tag
#$ -V -cwd
#$ -o logs
#$ -e logs
#$ -pe smp 8
#$ -t 1-22:1
#$ -tc 22

chr=$SGE_TASK_ID
kb=$1
wind=$2
r2=$3
type=$4
r2_tag=$5
kb_tag=$6
datadir=/projects/CARDIPS/analysis/epigenome_resource/analyses/jennifer/gwas_independent
outdir=/projects/CARDIPS/analysis/epigenome_resource/analyses/tim/gwas_coloc2
manifest=/projects/CARDIPS/analysis/epigenome_resource/analyses/tim/gwas_coloc/scripts/manifest_subset.txt
dir=`mktemp -d -p /projects/CARDIPS/analysis/epigenome_resource/analyses/tim/gwas_coloc2/scratch`

tail -n +2 $manifest| cut -f5 | while read name;
do
    file=${datadir}/indep/${name}.hg38_5e-08.txt
    echo $file
    name=`basename $file`
    name=${name::-4}
    
    if [ ! -d ${outdir}/qtl_tag2/${name} ]
    then
        mkdir ${outdir}/qtl_tag2/${name}
    fi
    
    cat ${datadir}/indep/${name}/chr${chr}.${kb}${type}_${wind}wind_r${r2}.prune.in ${outdir}/qtl_tag2/gwas_topsnp2.txt | sort -u >> ${dir}/${name}_${chr}_snps.txt

    cmd="plink --memory 30000 --threads 8 --bfile ${datadir}/reference/chr${chr}.renamed --extract ${dir}/${name}_${chr}_snps.txt --show-tags all --tag-kb ${kb_tag} --tag-r2 ${r2_tag} --out ${outdir}/qtl_tag2/${name}/chr${chr}.${kb}${type}_${wind}step_r${r2}"
    echo $cmd >& 2; eval $cmd
    
done

cmd="rm -r ${dir}"
echo $cmd >& 2; eval $cmd
