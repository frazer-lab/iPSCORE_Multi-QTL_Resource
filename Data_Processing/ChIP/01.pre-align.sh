source /frazer01/home/jennifer/.bash_profile
source activate encode-chip

encode_dir=/frazer01/home/jennifer/software/encode-atac-seq-pipeline

id=$1
fq_dir=$2
out_dir=$3

echo Sample ID: $id >& 2
echo Fastq directory: $fq_dir >& 2
echo Output directory: $out_dir >& 2
    
log=${out_dir}/pipeline.log
fastqs="${out_dir}/detect_adapters/R1.fastq.gz ${out_dir}/detect_adapters/R2.fastq.gz"

date >& 2
echo "Concatenate fastqs" >& 2
cmd="cat ${fq_dir}/*R1* > ${out_dir}/detect_adapters/R1.fastq.gz"
echo $cmd >> $log; echo $cmd >& 2; eval $cmd

date >& 2
cmd="cat ${fq_dir}/*R2* > ${out_dir}/detect_adapters/R2.fastq.gz"
echo $cmd >> $log; echo $cmd >& 2; eval $cmd

date >& 2
echo "Detect adapters.." >& 2
fqs=($fastqs)
for fq in ${fqs[@]}; do
    date >& 2
    prefix=`basename $fq .fastq.gz`
    out=${out_dir}/detect_adapters/${prefix}.adapter.txt
    cmd="python3 ${encode_dir}/src/detect_adapter.py $fq > $out"
    echo $cmd >> $log; echo $cmd >& 2; eval $cmd
done

date >& 2
echo "Trim adapters.." >& 2
ADAPTER_FWD=`cat ${out_dir}/detect_adapters/*R1*.adapter.txt`
ADAPTER_REV=`cat ${out_dir}/detect_adapters/*R2*.adapter.txt`
echo "ADAPTER_FWD:" $ADAPTER_FWD >& 2
echo "ADAPTER_REV:" $ADAPTER_REV >& 2

date >& 2
cmd="cutadapt --cores 12 -a $ADAPTER_FWD -A $ADAPTER_REV -o ${out_dir}/detect_adapters/reads.1.fastq -p ${out_dir}/detect_adapters/reads.2.fastq ${fqs[0]} ${fqs[1]}"
echo $cmd >> $log; echo $cmd >& 2; eval $cmd

date >& 2
echo "Run fastqc.." >& 2
if [ ! -f ${out_dir}/detect_adapters/reads.1_fastqc.zip ]
then 
    date >& 2
    cmd="fastqc -o ${out_dir}/detect_adapters ${out_dir}/detect_adapters/reads.*.fastq"
    echo $cmd >> $log; echo $cmd >& 2; eval $cmd
    
    date >& 2
    cmd="unzip ${out_dir}/detect_adapters/reads.1_fastqc.zip -d ${out_dir}/detect_adapters; unzip ${out_dir}/detect_adapters/reads.2_fastqc.zip -d ${out_dir}/detect_adapters"
    echo $cmd >> $log; echo $cmd >& 2; eval $cmd
else
    echo Fastqc already done!
fi

date >& 2
echo "Clean up" >& 2
cmd="rm ${out_dir}/detect_adapters/R1.fastq.gz ${out_dir}/detect_adapters/R2.fastq.gz"
echo $cmd >> $log; echo $cmd >& 2; eval $cmd