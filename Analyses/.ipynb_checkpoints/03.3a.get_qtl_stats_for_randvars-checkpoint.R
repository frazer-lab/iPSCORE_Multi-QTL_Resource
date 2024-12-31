# Author: Jennifer Nguyen

setwd("/projects/CARDIPS/analysis/epigenome_resource")
source("analyses/jennifer/notebooks/functions.R")

set.seed(5366)

option_list = list(make_option("--tissue", type = "character", default = NA, help = "tissue", metavar = "character"),
                   make_option("--qtl_list", type = "character", default = NA, help = "table of SNP-eGenes to test", metavar = "character"),
                   make_option("--date", type = "character", default = NA, help = "date", metavar = "character"))

# parse arguments
opt_parser    = OptionParser(option_list = option_list)
opt           = parse_args(opt_parser)

# set arguments
tissue  = opt$tissue
qtl_list = opt$qtl_list
date = opt$date

# variants to get qtl stats for
variants = readLines(qtl_list)

# read all gene-variant pairs
file = paste("analyses/jennifer/ipscore_unique_qtls/mashr", date, "input",  "all_ipscore_qtlstats", paste(tissue, "txt", sep = "."), sep = "/")

message(tissue)

# filter for same genes and variants
message("filter for variants..")
qtls = fread(file, data.table = F) %>% 
        filter(gene_variant %in% variants) %>% 
        select(gene_variant, variant_id, slope, slope_se, pval_nominal) %>% 
        mutate(tissue = tissue) %>% 
        dplyr::rename(beta = slope, se = slope_se, pval = pval_nominal) 

# Save
outfile = paste("analyses/jennifer/ipscore_unique_qtls/mashr", date, "input", "randvars_200k_qtlstats", paste(tissue, "txt", sep = "."), sep = "/")

fwrite(qtls, outfile, row.names = F, sep = "\t")
message(paste("Saved:", outfile))