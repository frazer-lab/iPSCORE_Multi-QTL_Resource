#!/bin/bash

#$ -N atac
#$ -V -cwd
#$ -o logs
#$ -e logs
#$ -pe smp 12

source /frazer01/home/jennifer/.bash_profile
source activate encode-atac

set -e

script_dir=/frazer01/home/jennifer/software/encode-atac-seq-pipeline

fq_dir=$1
merge_uuid=$2
map_file=$3 #premerge to merge map file (column 1 premerge uuids, column 2 merge uuids)
out_dir=$4
log=$5

echo fq_dir=$fq_dir
echo merge_uuid=$merge_uuid
echo map_file=$map_file
echo out_dir=$out_dir
echo log=$log

date >& 2

# 0. Make scratch
dir=`mktemp -d -p scratch`

###### Prepare FASTQs

# 1. Move all fastqs into scratch
premerge_uuids=(`grep $merge_uuid $map_file | cut -f1`)
for premerge in ${premerge_uuids[@]}; do
    cmd="rsync ${fq_dir}/${premerge}/*.fastq.gz ${dir}"
    echo $cmd >& 2; echo $cmd >> $log; eval $cmd; date >& 2
done

# 2. Concatenate FQ1
cmd="zcat ${dir}/*R1*.fastq.gz > ${dir}/R1.fastq"
echo $cmd >& 2; echo $cmd >> $log; eval $cmd; date >& 2

date >& 2

# 3. Concatenate FQ2
cmd="zcat ${dir}/*R2*.fastq.gz > ${dir}/R2.fastq"
echo $cmd >& 2; echo $cmd >> $log; eval $cmd; date >& 2

# 4. Clean
cmd="rm ${dir}/*.fastq.gz"
echo $cmd >& 2; echo $cmd >> $log; eval $cmd; date >& 2

# 2. Get adapter sequence
in_files=( ${dir}/R1.fastq ${dir}/R2.fastq )
echo fastq input.. ${in_files[@]} >& 2

for fq in ${in_files[@]}; do
    prefix=`basename $fq .fastq.gz`
        out=${out_dir}/detect_adapters/${prefix}.adapter.txt
        cmd="python3 ${script_dir}/src/detect_adapter.py $fq > $out"
        echo $cmd >& 2; echo $cmd >> $log; eval $cmd
done

date >& 2

adapter_fwd=`cat ${out_dir}/detect_adapters/R1*.adapter.txt`
adapter_rev=`cat ${out_dir}/detect_adapters/R2*.adapter.txt`

# 3. Trim adapters
cmd="cutadapt --cores 4 -a $adapter_fwd -A $adapter_rev -o ${out_dir}/detect_adapters/reads.1.fastq -p ${out_dir}/detect_adapters/reads.2.fastq ${in_files[0]} ${in_files[1]}"
echo $cmd >& 2; echo $cmd >> $log; eval $cmd; date >& 2

date >& 2

###### Align reads

# 4. Align reads
ref=/reference/public/ucsc/hg38/hg38.fa.gz
fastq1=${out_dir}/detect_adapters/reads.1.fastq
fastq2=${out_dir}/detect_adapters/reads.2.fastq
cmd="bwa mem -t 8 $ref $fastq1 $fastq2 | samtools sort -@ 12 -O bam -o ${out_dir}/Aligned.bam -"
echo $cmd >& 2; echo $cmd >> $log; eval $cmd; date >& 2

# 5. Sort 
cmd="samtools sort -@ 8 -o ${out_dir}/Aligned.sorted.bam ${out_dir}/Aligned.bam"
echo $cmd >& 2; echo $cmd >> $log; eval $cmd; date >& 2

# 6. Clean
cmd="rm $fastq1 $fastq2 ${out_dir}/Aligned.bam"
echo $cmd >& 2; echo $cmd >> $log; eval $cmd; date >& 2

# 7. Flagstat
#cmd="samtools sort -n -@ 8 ${BAM_SORTED_FILE} -O SAM | SAMstats --sorted_sam_file -  --outf ${FLAGSTAT_FILE}" # this was used in the paper but recommend regular flagstat below
cmd="samtools flagstat -@ 8 ${out_dir}/Aligned.sorted.bam > ${out_dir}/qc/Aligned.sorted.flagstat"
echo $cmd >& 2; echo $cmd >> $log; eval $cmd; date >& 2

# 8. Index
cmd="samtools index -@ 8 ${out_dir}/Aligned.sorted.bam"
echo $cmd >& 2; echo $cmd >> $log; eval $cmd; date >& 2

###### Calculate number of mitochondrial reads

# 9. Make mito-free bam
cmd="samtools idxstats ${out_dir}/Aligned.sorted.bam | cut -f 1 | grep -v -P "^chrM" | xargs samtools view ${out_dir}/Aligned.sorted.bam -@ 12 -b> ${out_dir}/Aligned.sorted.non_mito.bam"
echo $cmd >& 2; echo $cmd >> $log; eval $cmd; date >& 2

# 10. Flagstat on mito-free bam
cmd="samtools flagstat -@ 8 ${out_dir}/Aligned.sorted.non_mito.bam > ${out_dir}/qc/Aligned.sorted.non_mito.flagstat"
echo $cmd >& 2; echo $cmd >> $log; eval $cmd; date >& 2

# 11. Make mito bam
cmd="samtools view -b ${out_dir}/Aligned.sorted.bam chrM > ${out_dir}/Aligned.sorted.mito.bam"
echo $cmd >& 2; echo $cmd >> $log; eval $cmd; date >& 2

# 12. Flagstat on mito bam
cmd="samtools flagstat -@ 8 ${out_dir}/Aligned.sorted.mito.bam > ${out_dir}/qc/Aligned.sorted.mito.flagstat"
echo $cmd >& 2; echo $cmd >> $log; eval $cmd; date >& 2

# 13. Clean
cmd="rm ${out_dir}/Aligned.sorted.mito.bam"
echo $cmd >& 2; echo $cmd >> $log; eval $cmd; date >& 2

###### Read filtering

# 14. Remove unmapped, mat unmapped, not primary alignment, and reads failing platform
cmd="samtools view -@ 8 -F 524 -f 2 -q 30 -u ${out_dir}/Aligned.sorted.bam | samtools sort -@ 8 -n -o ${out_dir}/Aligned.sorted.filt.nmsrt.bam -"
echo $cmd >& 2; echo $cmd >> $log; eval $cmd; date >& 2

# 15. Randomly assign multimappers
ENCODE_DIR=/frazer01/home/jennifer/software/encode-atac-seq-pipeline/src
MULTIMAPPING=4
cmd="samtools view -h ${out_dir}/Aligned.sorted.filt.nmsrt.bam | ${ENCODE_DIR}/assign_multimappers.py -k $MULTIMAPPING --paired-end | samtools fixmate -@ 8 -r - ${out_dir}/Aligned.sorted.filt.nmsrt.fixmate.bam"
echo $cmd >& 2; echo $cmd >> $log; eval $cmd; date >& 2

# 16. Again, remove reads
cmd="samtools view -@ 8 -F 1804 -f 2 -q 30 -u ${out_dir}/Aligned.sorted.filt.nmsrt.fixmate.bam | samtools sort -@ 8 -o ${out_dir}/Aligned.sorted.filt.bam -"
echo $cmd >& 2; echo $cmd >> $log; eval $cmd; date >& 2

# 17. Clean
cmd="rm ${out_dir}/Aligned.sorted.filt.nmsrt.fixmate.bam ${out_dir}/Aligned.sorted.filt.nmsrt.bam"
echo $cmd >& 2; echo $cmd >> $log; eval $cmd; date >& 2

# 18. Mark duplicates
cmd="picard MarkDuplicates \
INPUT=${out_dir}/Aligned.sorted.filt.bam \
OUTPUT=${out_dir}/Aligned.sorted.filt.dupmark.bam \
METRICS_FILE=${out_dir}/qc/Aligned.sorted.filt.dup.qc \
VALIDATION_STRINGENCY=LENIENT ASSUME_SORTED=true REMOVE_DUPLICATES=false"
echo $cmd >& 2; echo $cmd >> $log; eval $cmd; date >& 2

# 19. Flagstat on dupmarked bam
cmd="samtools flagstat -@ 8 ${out_dir}/Aligned.sorted.filt.dupmark.bam > ${out_dir}/qc/${out_dir}/Aligned.sorted.filt.dupmark.flagstat"
echo $cmd >& 2; echo $cmd >> $log; eval $cmd; date >& 2

# 20. Clean
cmd="mv ${out_dir}/Aligned.sorted.filt.dupmark.bam ${out_dir}/Aligned.sorted.filt.bam"
echo $cmd >& 2; echo $cmd >> $log; eval $cmd; date >& 2

# 21. Remove duplicates
cmd="samtools view -@ 8 -F 1804 -f 2 -b ${out_dir}/Aligned.sorted.filt.bam > ${out_dir}/Aligned.sorted.filt.nodup.bam"
echo $cmd >& 2; echo $cmd >> $log; eval $cmd; date >& 2

# 22. Index
cmd="samtools index -@ 8 ${out_dir}/Aligned.sorted.filt.nodup.bam"
echo $cmd >& 2; echo $cmd >> $log; eval $cmd; date >& 2

# 23. Flagstat
cmd="samtools flagstat -@ 8 ${out_dir}/Aligned.sorted.filt.nodup.bam > ${out_dir}/qc/Aligned.sorted.filt.nodup.flagstat"
echo $cmd >& 2; echo $cmd >> $log; eval $cmd; date >& 2

# 24. Clean
cmd="rm ${out_dir}/Aligned.sorted.filt.bam"
echo $cmd >& 2; echo $cmd >> $log; eval $cmd; date >& 2

# 25. Remove mitochondrial reads
cmd="samtools idxstats ${out_dir}/Aligned.sorted.filt.nodup.bam | cut -f 1 | grep -v -P "^chrM" | xargs samtools view ${out_dir}/Aligned.sorted.filt.nodup.bam -@ 8 -b > ${out_dir}/Aligned.sorted.filt.nodup.nomito.bam"
echo $cmd >& 2; echo $cmd >> $log; eval $cmd; date >& 2

# 26. Index
cmd="samtools index -@ 8 ${out_dir}/Aligned.sorted.filt.nodup.nomito.bam"
echo $cmd >& 2; echo $cmd >> $log; eval $cmd; date >& 2
















