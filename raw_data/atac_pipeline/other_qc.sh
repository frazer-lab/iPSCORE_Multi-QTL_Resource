#!/bin/bash

#$ -N qc
#$ -V -cwd
#$ -pe smp 8

source /frazer01/home/jennifer/.bash_profile
source activate encode-atac

# ATAC-seq pipeline
# Documentation: https://docs.google.com/document/d/1f0Cm4vRyDQDu0bMehHD7P7KOMxTOP-HiNoIvL1VcBt8/edit
# Github: https://github.com/ENCODE-DCC/atac-seq-pipeline/tree/master

# 1. Library Complexity
pbc_file_qc=${out_dir}/qc/Aligned.sorted.filt.nodup.nomito.pbc.qc

cmd="samtools sort -@ 8 -n ${out_dir}/Aligned.sorted.filt.nodup.nomito.bam -o ${out_dir}/Aligned.sorted.filt.nodup.nomito.srt.tmp.bam"
echo $cmd >& 2; echo $cmd >> $log; eval $cmd; date >& 2

bedtools bamtobed -bedpe -i ${out_dir}/Aligned.sorted.filt.nodup.nomito.srt.tmp.bam | awk 'BEGIN{OFS="\t"}{print $1,$2,$4,$6,$9,$10}' | grep -v 'chrM' | sort | uniq -c | awk 'BEGIN{mt=0;m0=0;m1=0;m2=0} ($1==1){m1=m1+1} ($1==2){m2=m2+1} {m0=m0+1} {mt=mt+$1} END{printf "%d\t%d\t%d\t%d\t%f\t%f\t%f\n",mt,m0,m1,m2,m0/mt,m1/m0,m1/m2}' > ${pbc_file_qc}

# 2. Create tagAlign file
cmd="bedtools bamtobed -bedpe -mate1 -i ${out_dir}/Aligned.sorted.filt.nodup.nomito.srt.tmp.bam | gzip -nc > ${out_dir}/Aligned.sorted.filt.nodup.nomito.bedpe.gz"
echo $cmd >& 2; echo $cmd >> $log; eval $cmd; date >& 2

zcat ${out_dir}/Aligned.sorted.filt.nodup.nomito.bedpe.gz | awk 'BEGIN{OFS="\t"}{printf "%s\t%s\t%s\tN\t1000\t%s\n%s\t%s\t%s\tN\t1000\t%s\n",$1,$2,$3,$9,$4,$5,$6,$10}' | gzip -nc > ${out_dir}/Aligned.sorted.filt.nodup.nomito.tagAlign.gz

# 3. Clean
cmd="rm ${out_dir}/Aligned.sorted.filt.nodup.nomito.srt.tmp.bam"
echo $cmd >& 2; echo $cmd >> $log; eval $cmd; date >& 2

# 4. Tn5 shifting
zcat -f ${out_dir}/Aligned.sorted.filt.nodup.nomito.tagAlign.gz | awk 'BEGIN {OFS = "\t"} { if ($6 == "+") {$2 = $2 + 4} else if ($6 == "-") {$3 = $3 - 5} if ($2 >= $3) { if ($6 == "+") {$2 = $3 - 1} else {$3 = $2 + 1} } print $0}' | gzip -nc > ${out_dir}/Aligned.sorted.filt.nodup.nomito.tn5.tagAlign.gz

# 5. Fragment Length
cmd="picard CollectInsertSizeMetrics \
INPUT=${out_dir}/Aligned.sorted.filt.nodup.nomito.bam \ 
OUTPUT=${out_dir}/qc/Aligned.sorted.filt.nodup.nomito.insertsize.metrics \
H=${out_dir}/qc/Aligned.sorted.filt.nodup.nomito.insertsize.metric.plot \
VERBOSITY=ERROR QUIET=TRUE \
USE_JDK_DEFLATER=TRUE USE_JDK_INFLATER=TRUE \
W=1000 STOP_AFTER=5000000"
echo $cmd >& 2; echo $cmd >> $log; eval $cmd; date >& 2

cmd="python3 ${script_dir}/plot_fragment_length.py ${out_dir}/qc/Aligned.sorted.filt.nodup.nomito.insertsize.metrics ${out_dir}/qc/Aligned.sorted.filt.nodup.nomito"
echo $cmd >& 2; echo $cmd >> $log; eval $cmd; date >& 2
