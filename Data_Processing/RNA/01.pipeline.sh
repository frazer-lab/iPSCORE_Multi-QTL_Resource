#$ -pe smp 12
#$ -V
#$ -cwd
#$ -e logs
#$ -o logs

set -e 

hostname >& 2

source /frazer01/home/jennifer/.bash_profile
source activate frazer-rna
export PATH=/software/biobambam2-2.0.95/bin:/software/sambamba_v0.6.7:$PATH

fq_dir=$1 # directory containing fastqs. ex: fastq/SRRXXXX
out_dir=$2
reference_vcf=$3
id=`basename $fq_dir`
out_dir=${out_dir}/${id}

dir=`mktemp -d -p scratch`
reference_star=/reference/private/STAR/hg38_gencode44/reference
reference_rsem=/reference/private/RSEM/hg38_gencode44/reference/hg38_gencode44

log=${out_dir}/pipeline.log

if [ ! -d $out_dir ]; then mkdir $out_dir; fi

if [ -f ${out_dir}/picard.RNA_Metrics ]
then
    echo "Sample $id is already processed!" >& 2
    cmd="rm -r ${dir}"
    echo $cmd >& 2; eval $cmd
    exit 1
fi

# 1. Combine fastqs
cmd="zcat ${fq_dir}/*R1*.fastq.gz > ${dir}/R1.fastq"
echo $cmd >> $log; echo $cmd >& 2; eval $cmd

cmd="zcat ${fq_dir}/*R2*.fastq.gz > ${dir}/R2.fastq"
echo $cmd >> $log; echo $cmd >& 2; eval $cmd

in_files=( ${dir}/R1.fastq ${dir}/R2.fastq )
echo fastq input.. ${in_files[@]} >& 2

# 2. Run STAR Alignment
# good STAR manual: https://physiology.med.cornell.edu/faculty/skrabanek/lab/angsd/lecture_notes/STARmanual.pdf
cmd="STAR \
--runThreadN 12 \
--genomeDir ${reference_star} \
--genomeLoad NoSharedMemory \
--readFilesIn ${in_files[@]} \
--outSAMattributes All \
--outSAMunmapped Within \
--outSAMattrRGline ID:1 PL:ILLUMINA PU:CARDIPS LB:${id} SM:${id} \
--outFilterMultimapNmax 20 \
--outFilterMismatchNmax 999 \
--alignIntronMin 20 \
--alignIntronMax 1000000 \
--alignMatesGapMax 1000000 \
--outSAMtype BAM Unsorted \
--outFileNamePrefix ${dir}/ \
--quantMode TranscriptomeSAM"
echo $cmd >> $log; echo $cmd >& 2; eval $cmd

cmd="rm ${dir}/R1.fastq ${dir}/R2.fastq"
echo $cmd >> $log; echo $cmd >& 2; eval $cmd

# 3. Sort bam 
#cmd="sambamba sort -m 32GB -t 12 --tmpdir ${dir} ${dir}/Aligned.out.bam"
cmd="samtools sort -m 2G -n -o ${dir}/Aligned.out.namesorted.bam -@ 12 ${dir}/Aligned.out.bam"
echo $cmd >> $log; echo $cmd >& 2; eval $cmd

# 4. Fill in mate coordinates
cmd="samtools fixmate -@ 12 -m ${dir}/Aligned.out.namesorted.bam ${dir}/Aligned.out.namesorted.fixmate.bam"
echo $cmd >> $log; echo $cmd >& 2; eval $cmd

# 5. Sort bam
cmd="samtools sort -m 2G -o ${dir}/Aligned.out.sorted.bam -@ 12 ${dir}/Aligned.out.namesorted.fixmate.bam"
echo $cmd >> $log; echo $cmd >& 2; eval $cmd

cmd="rm ${dir}/Aligned.out.namesorted.bam ${dir}/Aligned.out.namesorted.fixmate.bam"
echo $cmd >> $log; echo $cmd >& 2; eval $cmd

# 6. Index bam
cmd="samtools index -@ 12 ${dir}/Aligned.out.sorted.bam ${dir}/Aligned.out.sorted.bam.bai"
echo $cmd >> $log; echo $cmd >& 2; eval $cmd

# 7. Mark duplicates
# Why we don't remove duplicates: https://dnatech.genomecenter.ucdavis.edu/faqs/should-i-remove-pcr-duplicates-from-my-rna-seq-data/
cmd="samtools markdup -@ 12 -s -T ${dir}/temp ${dir}/Aligned.out.sorted.bam ${dir}/Aligned.out.sorted.mdup.bam"
echo $cmd >> $log; echo $cmd >& 2; eval $cmd

cmd="rm ${dir}/Aligned.out.sorted.bam ${dir}/Aligned.out.sorted.bam.bai"
echo $cmd >> $log; echo $cmd >& 2; eval $cmd

# 8. Index bam
cmd="samtools index -@ 12 ${dir}/Aligned.out.sorted.mdup.bam ${dir}/Aligned.out.sorted.mdup.bam.bai"
echo $cmd >> $log; echo $cmd >& 2; eval $cmd

# 9. Run flagstat
cmd="samtools flagstat -@ 12 ${dir}/Aligned.out.sorted.mdup.bam > ${dir}/Aligned.out.sorted.mdup.bam.flagstat"
echo $cmd >> $log; echo $cmd >& 2; eval $cmd

# 10. Run idxstats
cmd="samtools idxstats -@ 12 ${dir}/Aligned.out.sorted.mdup.bam > ${dir}/Aligned.out.sorted.mdup.bam.idxstats"
echo $cmd >> $log; echo $cmd >& 2; eval $cmd

# 11. Run samtools stats
cmd="samtools stats -@ 12 ${dir}/Aligned.out.sorted.mdup.bam > ${dir}/Aligned.out.sorted.mdup.bam.stats"
echo $cmd >> $log; echo $cmd >& 2; eval $cmd

# 12. Reformat header for picard
cmd="samtools view -H ${dir}/Aligned.out.sorted.mdup.bam | sed 's,^@RG.*,@RG\tID:None\tSM:None\tLB:None\tPL:Illumina,g' | samtools reheader - ${dir}/Aligned.out.sorted.mdup.bam > ${dir}/Aligned.out.sorted.mdup.reheader.bam"
echo $cmd >> $log; echo $cmd >& 2; eval $cmd

# 13. Run Picard's CollectRnaSeqMetrics
refFlat=/reference/public/UCSC/hg38/refFlat.txt.gz

cmd="picard CollectRnaSeqMetrics \
I=${dir}/Aligned.out.sorted.mdup.reheader.bam \
O=${dir}/picard.RNA_Metrics \
VALIDATION_STRINGENCY=SILENT \
REF_FLAT=$refFlat \
STRAND=NONE"
echo $cmd >> $log; echo $cmd >& 2; eval $cmd

# 14. Clean
cmd="rm -r ${dir}/Aligned.out.bam ${dir}/Aligned.out.sorted.mdup.reheader.bam ${dir}/Aligned.toTranscriptome.out.bam"
echo $cmd >> $log; echo $cmd >& 2; eval $cmd

# 15. Move
cmd="mv ${dir}/* ${out_dir}"
echo $cmd >> $log; echo $cmd >& 2; eval $cmd

# 16. Identity
cmd="bash identity.sh $reference_vcf $out_dir $out_dir $id"
echo $cmd >> $log; echo $cmd >& 2; eval $cmd

# 17. Delete scratch
cmd="rm -r ${dir}"
echo $cmd >> $log; echo $cmd >& 2; eval $cmd
