# Author: Jennifer Nguyen

setwd("/projects/CARDIPS/analysis/epigenome_resource")
source("analyses/jennifer/notebooks/functions.R")

set.seed(5366)

option_list = list(make_option("--tissue", type = "character", default = NA, help = "tissue", metavar = "character"),
                   make_option("--qtl_list", type = "character", default = NA, help = "table of SNP-eGenes to test", metavar = "character"),
                   make_option("--date", type = "character", default = NA, help = "date of analysis", metavar = "character"))

# parse arguments
opt_parser    = OptionParser(option_list = option_list)
opt           = parse_args(opt_parser)

# set arguments
tissue   = opt$tissue
qtl_list = opt$qtl_list
date     = opt$date

# variants to get qtl stats for
snps2test = fread(qtl_list, data.table = F)

# read all gene-variant pairs
message("opening allpairs..")
allpairs = paste("/reference/public/GTEX_v8/all_pairs/GTEx_Analysis_v8_eQTL_all_associations", paste(tissue, "allpairs.txt.gz", sep = "."), sep = "/")

message(tissue)

# filter for same genes and variants
message("filter for variants..")
qtls = fread(allpairs, data.table = F) %>% 
    filter(variant_id %in% c(snps2test$gtex_snp1, snps2test$gtex_snp2)) %>% 
    tidyr::separate(variant_id, into = c("chr", "pos", "ref", "alt", "genome")) %>%
    mutate(id1 = paste("VAR", gsub("chr", "", chr), pos, ref, alt, sep = "_"),
           id2 = paste("VAR", gsub("chr", "", chr), pos, alt, ref, sep = "_")) %>%
    mutate(gene_id = simplify_id(gene_id)) %>% 
    mutate(gene_variant1 = paste(gene_id, id1),
           gene_variant2 = paste(gene_id, id2)) %>% 
    filter(gene_id %in% snps2test$gene_id)

dir.create(paste(getwd(), "analyses/jennifer/ipscore_unique_qtls/mashr", date, sep = "/"))
dir.create(paste(getwd(), "analyses/jennifer/ipscore_unique_qtls/mashr", date, "leadvars", sep = "/"))
outfile = paste("analyses/jennifer/ipscore_unique_qtls/mashr", date, "input", "leadvars", paste(tissue,  "txt", sep = "."), sep = "/")
  
message("writing..")
fwrite(qtls, outfile, row.names = F, sep = "\t")
message(paste("Saved:", outfile))
