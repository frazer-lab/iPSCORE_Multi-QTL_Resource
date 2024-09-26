# Author: Jennifer Nguyen

input_file=$1
date=$2

header=$(cat $input_file | head -n 1)

# Column number to split by (1-based index)
column_number=1
outdir=/projects/CARDIPS/analysis/epigenome_resource/analyses/jennifer/ipscore_unique_qtls/mashr/${date}/input/by_gene
mkdir $outdir
cd $outdir

# Using zcat and awk to split the file
cat $input_file | tail -n +2 | awk -v col=$column_number -v header="$header" '
{
    split($col, arr, ".")
    file = arr[1] ".txt"
    print "Writing to file: " file
    if (!(file in seen)) {
        print header > file
        seen[file] = 1
    }
    print >> file
}'
