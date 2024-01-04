#$ -o logs/create.out
#$ -e logs/create.err
#$ -pe smp 32
#$ -V
#$ -cwd

#tutorial: https://www.reneshbedre.com/blog/star-aligner.html
export PATH=/software/STAR-2.7.10b/bin/Linux_x86_64:$PATH
export PATH=/software/STAR-2.7.10b/source/:$PATH

# 0. Make scratch
dir=`mktemp -d -p /scratch`

cmd="rm -r Log.out _STARtmp"
echo $cmd >& 2
$cmd

date >& 2

# 1. Copy GTF and fasta reference
cmd="rsync /reference/public/Gencode.v44lift38/gencode.v44.annotation.gtf.gz ${dir}/gencode.v44.annotation.gtf.gz"
echo $cmd >& 2
$cmd &

date >& 2

cmd="rsync /reference/public/UCSC/hg38/hg38.fa.gz ${dir}/hg38.fa.gz"
echo $cmd >& 2
$cmd &

date >& 2

# 2. Make output directory
cmd="mkdir ${dir}/reference"
echo $cmd >& 2
$cmd &

date >& 2

cmd="wait"
echo $cmd >& 2
$cmd

# 3. Gzip GTF and fasta 
cmd="gunzip ${dir}/gencode.v44.annotation.gtf.gz"
echo $cmd >& 2
$cmd &

cmd="gunzip ${dir}/hg38.fa.gz"
echo $cmd >& 2
$cmd &

cmd="wait"
echo $cmd >& 2
$cmd

date >& 2

# 4. Generate STAR reference
# references on sjdbOverhang: https://groups.google.com/g/rna-star/c/h9oh10UlvhI/m/BfSPGivUHmsJ; https://groups.google.com/g/rna-star/c/pHKU0vGGvGk
cmd="STAR --runThreadN 32 --runMode genomeGenerate --genomeDir ${dir}/reference \
--genomeFastaFiles ${dir}/hg38.fa --sjdbGTFfile ${dir}/gencode.v44.annotation.gtf \
--sjdbOverhang 99"
echo $cmd >& 2
$cmd

date >& 2

# 5. Move output to final directory
cmd="mv ${dir}/reference ./"
echo $cmd >& 2
$cmd

date >& 2

# 6. Delete scratch
cmd="rm -rf ${dir}"
echo $cmd >& 2
$cmd

date >& 2
