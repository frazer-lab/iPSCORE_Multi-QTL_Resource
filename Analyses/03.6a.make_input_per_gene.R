# Author: Jennifer Nguyen

setwd("/projects/CARDIPS/analysis/epigenome_resource")
source("analyses/jennifer/notebooks/functions.R")

set.seed(5366)

option_list = list(make_option("--gene_id", type = "character", default = NA, help = "gene_id", metavar = "character"),
                   make_option("--date", type = "character", default = NA, help = "date", metavar = "character"))

# parse arguments
opt_parser    = OptionParser(option_list = option_list)
opt           = parse_args(opt_parser)

# set arguments
gene_id = opt$gene_id
date = opt$date
filename = paste("/projects/CARDIPS/analysis/epigenome_resource/analyses/jennifer/ipscore_unique_qtls/mashr", date, "input", "by_gene", paste(gene_id, "txt", sep = "."), sep = "/")

message(paste("Reading", filename))
this_stats = fread(filename, data.table = F)

# Convert to matrix
Bhat       = suppressMessages(reshape2::dcast(data = this_stats, gene_variant ~ tissue, value.var = "beta", fun.aggregate = mean))
Shat       = suppressMessages(reshape2::dcast(data = this_stats, gene_variant ~ tissue, value.var = "se"  , fun.aggregate = mean))

# Fill in missing
Bhat[is.na(Bhat)] = 1e-6
Shat[is.na(Shat)] = 1

# Count number of conditions
message(paste("Detected", ncol(Bhat), ncol(Shat), "conditions"))
message(paste("Detected", nrow(Bhat), nrow(Shat), "gene-variant pairs"))

# Write out
datalist = list("bhat" = Bhat, "shat" = Shat)
outfile = paste(getwd(), "analyses/jennifer/ipscore_unique_qtls/mashr", date, "input", "by_gene", paste(gene_id, "robj", sep = "."), sep = "/")
save(datalist, file = outfile)
message(paste("Saved:", outfile))    




