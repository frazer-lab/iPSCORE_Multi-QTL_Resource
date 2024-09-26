#!/bin/bash
#$ -cwd
#$ -V
#$ -j y
#$ -S /bin/bash
#$ -t 1-XX # XXX number of egenes tested from mashr


SNP_DIR="path_to_tmpdir_for_plink_snps"
PLINK_DIR="path_to_1000G_plink_files"
OUTPUT_DIR="path_to_out_directory"

FILE_PATH=$(ls $SNP_DIR | sed -n "${SGE_TASK_ID}p")
FILE_NAME=$(basename "$FILE_PATH" "_snps.txt")

# Extract the chromosome number from the first SNP in the file
CHROM=$(awk 'NR==1 {split($1, a, "_"); print a[1]}' $SNP_DIR/$FILE_PATH)

# Path to the corresponding PLINK chromosome file
PLINK_FILE="${PLINK_DIR}/chr${CHROM}.renamed"

# Output file for this gene
OUTPUT_FILE="${OUTPUT_DIR}/${FILE_NAME}"

# Execute PLINK to calculate LD
#plink --bfile $PLINK_FILE --ld $SNP_DIR/$FILE_PATH --r2 --ld-window-r2 0 --out $OUTPUT_FILE



plink --bfile ${PLINK_DIR}/chr${CHROM}.renamed --extract $SNP_DIR/$FILE_PATH --make-bed --out $OUTPUT_FILE

# Step 2: Calculate linkage disequilibrium (LD) for the specific chromosome
plink --bfile ${OUTPUT_FILE} --r2 --ld-window-r2 0 --ld-window 99999 --ld-window-kb 10000 --out ${OUTPUT_FILE}

rm ${OUTPUT_FILE}.bim ${OUTPUT_FILE}.bed ${OUTPUT_FILE}.fam ${OUTPUT_FILE}.log ${OUTPUT_FILE}.nosex
