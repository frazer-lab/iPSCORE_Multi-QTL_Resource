#$ -N rna
#$ -pe smp 16
#$ -V -cwd

set -e

hostname >& 2

source /frazer01/home/jennifer/.bash_profile
source activate frazer-rna
export PATH=/software/biobambam2-2.0.95/bin:/software/sambamba_v0.6.7:$PATH
reference_star=/reference/private/STAR/hg38_gencode44/reference
reference_rsem=/reference/private/RSEM/hg38_gencode44/reference/hg38_gencode44

args=($@)

if [ ${#args[@]} -ne 4 ]
then
    echo "Exiting. Need input: sh merge_pipeline.sh <fastq_dir> <merge_uuid_list> <map_file> <out_dir>" >& 2
    exit 1
fi

dir=`mktemp -d -p scratch`
if [ ! -d ${dir} ]; then exit 1; fi
fq_dir=$1
merge_uuid=`tail -n +$SGE_TASK_ID $2 | head -1`
map_file=$3 #premerge to merge map file (column 1 premerge uuids, column 2 merge uuids)
out_dir=$4
out_dir=${out_dir}/${merge_uuid}

if [ -f ${out_dir}/rsem.genes.results ]
then
    echo "Sample ${merge_uuid} is already processed!" >& 2
    cmd="rm -r ${dir}"
    echo $cmd >& 2; eval $cmd
    exit 1
fi

premerge_uuids=(`grep $merge_uuid $map_file | cut -f1`)

echo Processing .. ${merge_uuid} >& 2
echo Premerge uuids.. ${premerge_uuids[@]} >& 2
echo Fastq directory.. ${fq_dir} >& 2
echo Output directory.. ${out_dir} >& 2
echo Working directory.. ${dir} >& 2

if [ ! -d $out_dir ]; then mkdir $out_dir; fi

# 1. Combine fastqs
date >& 2

for premerge in ${premerge_uuids[@]}; do
    cmd="rsync ${fq_dir}/${premerge}/*.fastq.gz ${dir}"
    echo $cmd >& 2; eval $cmd
done

date >& 2

cmd="zcat ${dir}/*R1*.fastq.gz > ${dir}/R1.fastq"
echo $cmd >& 2; eval $cmd

date >& 2

cmd="zcat ${dir}/*R2*.fastq.gz > ${dir}/R2.fastq"
echo $cmd >& 2; eval $cmd

cmd="rm ${dir}/*.fastq.gz"
echo $cmd >& 2; eval $cmd

in_files=( ${dir}/R1.fastq ${dir}/R2.fastq )
echo fastq input.. ${in_files[@]} >& 2

date >& 2

# 2. Run STAR Alignment
# good STAR manual: https://physiology.med.cornell.edu/faculty/skrabanek/lab/angsd/lecture_notes/STARmanual.pdf
cmd="STAR \
--runThreadN 16 \
--genomeDir ${reference_star} \
--genomeLoad NoSharedMemory \
--readFilesIn ${in_files[@]} \
--outSAMattributes All \
--outSAMunmapped Within \
--outSAMattrRGline ID:1 PL:ILLUMINA PU:CARDIPS LB:${name} SM:${name} \
--outFilterMultimapNmax 20 \
--outFilterMismatchNmax 999 \
--alignIntronMin 20 \
--alignIntronMax 1000000 \
--alignMatesGapMax 1000000 \
--outSAMtype BAM Unsorted \
--outFileNamePrefix ${dir}/ \
--quantMode TranscriptomeSAM"
echo $cmd >& 2; eval $cmd

cmd="rm ${dir}/R1.fastq ${dir}/R2.fastq"
echo $cmd >& 2; eval $cmd

# 3. Sort bam 
date >& 2
#cmd="sambamba sort -m 32GB -t 16 --tmpdir ${dir} ${dir}/Aligned.out.bam"
cmd="samtools sort -m 2G -n -o ${dir}/Aligned.out.namesorted.bam -@ 16 ${dir}/Aligned.out.bam"
echo $cmd >& 2; eval $cmd

# 4. Fill in mate coordinates
date >& 2
cmd="samtools fixmate -@ 16 -m ${dir}/Aligned.out.namesorted.bam ${dir}/Aligned.out.namesorted.fixmate.bam"
echo $cmd >& 2; eval $cmd

# 5. Sort bam
date >& 2
cmd="samtools sort -m 2G -o ${dir}/Aligned.out.sorted.bam -@ 16 ${dir}/Aligned.out.namesorted.fixmate.bam"
echo $cmd >& 2; eval $cmd

cmd="rm ${dir}/Aligned.out.namesorted.bam ${dir}/Aligned.out.namesorted.fixmate.bam"
echo $cmd >& 2; eval $cmd

# 6. Index bam
date >& 2
cmd="samtools index -@ 16 ${dir}/Aligned.out.sorted.bam ${dir}/Aligned.out.sorted.bam.bai"
echo $cmd >& 2; eval $cmd

# 7. Mark duplicates
date >& 2
# Why we don't remove duplicates: https://dnatech.genomecenter.ucdavis.edu/faqs/should-i-remove-pcr-duplicates-from-my-rna-seq-data/
cmd="samtools markdup -@ 16 -s -T ${dir}/temp ${dir}/Aligned.out.sorted.bam ${dir}/Aligned.out.sorted.mdup.bam"
echo $cmd >& 2; eval $cmd

cmd="rm ${dir}/Aligned.out.sorted.bam ${dir}/Aligned.out.sorted.bam.bai"
echo $cmd >& 2; eval $cmd

# 8. Index bam
date >& 2
cmd="samtools index -@ 16 ${dir}/Aligned.out.sorted.mdup.bam ${dir}/Aligned.out.sorted.mdup.bam.bai"
echo $cmd >& 2; eval $cmd

# 9. Run flagstat
date >& 2
cmd="samtools flagstat -@ 16 ${dir}/Aligned.out.sorted.mdup.bam > ${dir}/Aligned.out.sorted.mdup.bam.flagstat"
echo $cmd >& 2; eval $cmd

# 10. Run idxstats
date >& 2
cmd="samtools idxstats -@ 16 ${dir}/Aligned.out.sorted.mdup.bam > ${dir}/Aligned.out.sorted.mdup.bam.idxstats"
echo $cmd >& 2; eval $cmd

# 10. Run samtools stats
date >& 2
cmd="samtools stats -@ 16 ${dir}/Aligned.out.sorted.mdup.bam > ${dir}/Aligned.out.sorted.mdup.bam.stats"
echo $cmd >& 2; eval $cmd

# 11. Run RSEM
date >& 2
cmd="rsem-calculate-expression \
--bam \
--num-threads 16 \
--no-bam-output \
--seed 3272015 \
--estimate-rspd \
--forward-prob 0 \
--paired-end"
cmd=$cmd" ${dir}/Aligned.toTranscriptome.out.bam ${reference_rsem} ${dir}/rsem"
echo $cmd >& 2; eval $cmd

# 12. Reformat header for picard
date >& 2
cmd="samtools view -H ${dir}/Aligned.out.sorted.mdup.bam | sed 's,^@RG.*,@RG\tID:None\tSM:None\tLB:None\tPL:Illumina,g' | samtools reheader - ${dir}/Aligned.out.sorted.mdup.bam > ${dir}/Aligned.out.sorted.mdup.reheader.bam"
echo $cmd >& 2; eval $cmd

# 14. Run Picard's CollectRnaSeqMetrics
date >& 2
refFlat=/reference/public/UCSC/hg38/refFlat.txt.gz

cmd="picard CollectRnaSeqMetrics \
I=${dir}/Aligned.out.sorted.mdup.reheader.bam \
O=${dir}/picard.RNA_Metrics \
VALIDATION_STRINGENCY=SILENT \
REF_FLAT=$refFlat \
STRAND=NONE"
echo $cmd >& 2; eval $cmd

# 15. Clean
date >& 2
cmd="rm -r ${dir}/Aligned.out.bam ${dir}/Aligned.out.sorted.mdup.reheader.bam ${dir}/Aligned.toTranscriptome.out.bam"
echo $cmd >& 2; eval $cmd

# Move
cmd="mv ${dir}/* ${out_dir}/"
echo $cmd >& 2; eval $cmd

cmd="rm -r ${dir}"
echo $cmd >& 2; eval $cmd

date >& 2
