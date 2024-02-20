source /frazer01/home/jennifer/.bash_profile
source activate encode-chip

set -e 

id=$1
out_dir=$2
fastq1=${out_dir}/detect_adapters/reads.1.fastq
fastq2=${out_dir}/detect_adapters/reads.2.fastq
sam=${out_dir}/Aligned.sam
bam=${out_dir}/Aligned.bam
flagstat_qc=${out_dir}/qc/Aligned.flagstat.qc
ref=/reference/public/ucsc/hg38/hg38.fa.gz
log=${out_dir}/pipeline.log

# 1. Align
cmd="bwa mem -t 12 $ref $fastq1 $fastq2 | samtools sort -@ 12 -O bam -o $bam -"
echo $cmd >> $log; echo $cmd >& 2; eval $cmd
       
# 2. Clean
cmd="rm $fastq1 $fastq2"
echo $cmd >> $log; echo $cmd >& 2; eval $cmd

# 3. Flagstat
cmd="samtools sort -n --threads 12 ${bam} -O SAM  | SAMstats  --sorted_sam_file - --outf ${flagstat_qc}"
echo $cmd >> $log; echo $cmd >& 2; eval $cmd
