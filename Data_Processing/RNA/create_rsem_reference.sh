#$ -V
#$ -cwd
#$ -pe smp 32
#$ -o logs
#$ -e logs

export PATH=/software/rsem-1.2.20:$PATH

# 0. Make scratch directory 
dir=`mktemp -d -p /scratch`

date >& 2

# 1. Copy GTF reference to scratch (server was unstable so need to copy everything into one directory)
cmd="rsync /reference/public/Gencode.v44lift38/gencode.v44.annotation.gtf.gz ${dir}/gencode.v44.annotation.gtf.gz"
echo $cmd >& 2
$cmd &

date >& 2

# 2. Copy fasta to scratch
cmd="rsync /reference/public/UCSC/hg38/hg38.fa.gz ${dir}/hg38.fa.gz"
echo $cmd >& 2
$cmd &

date >& 2

cmd="wait"
echo $cmd >& 2
$cmd

date >& 2

# 3. Gzip fasta and gtf reference
cmd="gunzip ${dir}/hg38.fa.gz"
echo $cmd >& 2
$cmd &

date >& 2

cmd="gunzip ${dir}/gencode.v44.annotation.gtf.gz"
echo $cmd >& 2
$cmd &

date >& 2

# 4. Make output directory
cmd="mkdir ${dir}/reference"
echo $cmd >& 2
$cmd &

date >& 2

cmd="wait"
echo $cmd >& 2
$cmd

date >& 2

# 5. Prepare RSEM reference
cmd="rsem-prepare-reference --gtf ${dir}/gencode.v44.annotation.gtf \
${dir}/hg38.fa ${dir}/reference/hg38_gencode44"
echo $cmd >& 2
$cmd

date >& 2

# 6. Move output to final directory
cmd="mv ${dir}/reference ."
echo $cmd >& 2
$cmd

date >& 2

# 7. Remove scratch
cmd="rm -rf ${dir}"
echo $cmd >& 2
$cmd

date >& 2
