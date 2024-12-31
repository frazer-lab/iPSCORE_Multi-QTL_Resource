# Author: Jennifer Nguyen

setwd("/projects/CARDIPS/analysis/epigenome_resource")
source("analyses/jennifer/notebooks/functions.R")

set.seed(5366)

option_list = list(make_option("--tissue", type = "character", default = NA, help = "tissue", metavar = "character"),
                   make_option("--qtl_list", type = "character", default = NA, help = "table of SNP-eGenes to test", metavar = "character"),
                   make_option("--date", type = "character", default = NA, help = "date (used as analysis id)", metavar = "character"))

# parse arguments
opt_parser    = OptionParser(option_list = option_list)
opt           = parse_args(opt_parser)

# set arguments
tissue  = opt$tissue
qtl_list = opt$qtl_list
date = opt$date

# variants to get qtl stats for
variants = fread(qtl_list, data.table = F)

# read all gene-variant pairs
allpairs = paste("/reference/public/GTEX_v8/all_pairs/GTEx_Analysis_v8_eQTL_all_associations", paste(tissue, "allpairs.txt.gz", sep = "."), sep = "/")

message(tissue)

# filter for same genes and variants
message("filter for variants..")
qtls = fread(allpairs, data.table = F) %>% filter(variant_id %in% variants$gtex_variant) 

message("change gene ids..")
qtls$gene_id = simplify_id(qtls$gene_id)

message("filter for genes..")
qtls = qtls %>% filter(gene_id %in% variants$gene_id)

qtls$gene_variant = paste(qtls$gene_id, gsub("_b38", "", gsub("chr", "VAR_", qtls$variant_id)))

message("filter for gene-variant pairs...")
qtls = qtls %>% filter(gene_variant %in% variants$gene_variant)

# Save
dir.create(paste("analyses/jennifer/ipscore_unique_qtls/mashr", date, "input", "all_ipscore_qtlstats", sep = "/"))
outfile = paste("analyses/jennifer/ipscore_unique_qtls/mashr", date, "input",  "all_ipscore_qtlstats",  paste(tissue, "txt", sep = "."), sep = "/")

fwrite(qtls, outfile, row.names = F, sep = "\t")
message(paste("Saved:", outfile))
