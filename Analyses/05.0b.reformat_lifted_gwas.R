# Author: Jennifer Nguyen

setwd("/projects/CARDIPS/analysis/epigenome_resource")
source("analyses/jennifer/notebooks/functions.R")

set.seed(5366)

option_list = list(make_option("--id", type = "character", default = NA, help = "trait id", metavar = "character"))
opt_parser = OptionParser(option_list = option_list)
opt = parse_args(opt_parser)
id = opt$id

message(paste("Input:", id))

work_dir = "/projects/CARDIPS/analysis/epigenome_resource/analyses/jennifer/gwas_liftover"
outfile = paste(work_dir, "hg38_summary_statistics", paste(id, "hg38.tsv", sep = "."), sep = "/")

message(paste(Sys.time(), "Reading", outfile))

data = fread(outfile, data.table = F, header = F)

out = data  %>% 
    select(-V1, -V2) %>% 
    mutate(V6 = V3) %>% 
    select(-V3)

cols = colnames(fread(cmd = paste("head -1", paste(work_dir, "hg19_summary_statistics", paste(id, "hg19.tsv", sep = "."), sep = "/")), data.table = F))
colnames(out) = cols

out = out %>% dplyr::relocate(id, .after = a2)
out$chr = paste0("chr", out$chr)
out$chrpos = paste(out$chr, out$pos, sep = "_")
if ("n_case" %in% colnames(out) & "n_control" %in% colnames(out) & !"total" %in% colnames(out))
{
    out$total = out$n_case + out$n_control
    out$cases_fr = out$n_case / out$total
}

fwrite(out, outfile, row.names = F, sep = "\t")

message(paste(Sys.time(), "Saved:", outfile))
