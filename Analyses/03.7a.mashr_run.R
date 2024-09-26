# Author: Jennifer Nguyen

setwd("/projects/CARDIPS/analysis/epigenome_resource")

set.seed(5366)

library(mashr)
library(ash)
library(data.table)
library(dplyr)
library(optparse)

option_list = list(make_option("--testname", type = "character", default = NA, help = "testname", metavar = "character"),
                   make_option("--date", type = "character", default = NA, help = "modelname", metavar = "character"))

# parse arguments
opt_parser    = OptionParser(option_list = option_list)
opt           = parse_args(opt_parser)

# set arguments
testname  = opt$testname
date      = opt$date

# Add rownames function
add_rownames = function(x)
{
    rownames(x) = x[,1]
    x[,1] = NULL
    return(x)
}

dir = "/projects/CARDIPS/analysis/epigenome_resource/analyses/jennifer/ipscore_unique_qtls/mashr"

dir.create(paste(dir, date, "results", "by_gene", sep = "/"))

input_datafile    = paste (dir, date, "input", "by_gene", paste(testname , "robj", sep = "."), sep = "/")
results_outfile   = paste (dir, date, "results", "by_gene", paste(testname , "robj", sep = "."), sep = "/")
model_outfile     = paste (dir, date, "results/mashr_model.robj", sep = "/")
conditions_file   = paste (dir, date, "input", "conditions.txt", sep = "/")

# Read conditions
conditions        = readLines(conditions_file)

# Open input datafile
message(paste("Opening", input_datafile))
load(input_datafile, verbose = T)

print(dim(datalist[["bhat"]]))
print(dim(datalist[["shat"]]))

# Open model
message(paste("Opening", model_outfile))
load(model_outfile, verbose = T)
m    = model_list[["m"]]
Vhat = model_list[["Vhat"]]

# Check if all conditions are present
missing = setdiff(conditions, colnames(datalist[["bhat"]]))
message(paste("Missing:", length(missing), "conditions"))
message(paste(missing, collapse = "\n"))

for (cond in missing)
{
    datalist[["bhat"]][cond] = 1e-6
    datalist[["shat"]][cond] = 1
}

datalist[["bhat"]] = as.matrix(add_rownames(datalist[["bhat"]])) 
datalist[["shat"]] = as.matrix(add_rownames(datalist[["shat"]])) 

if (nrow(datalist[["bhat"]] ) > 1)
{
    datalist[["shat"]] = datalist[["shat"]][rownames(datalist[["bhat"]]),] 
}

# Prepare input for mash using the same null correlation matrix
data.input = mash_set_data(datalist[["bhat"]],datalist[["shat"]], V=Vhat)

# Run mash
m2 = mash(data.input, g=get_fitted_g(m), fixg=TRUE)

# Get LFSR
lfsr = as.data.frame(get_lfsr(m2))

outlist = list("m2" = m2, "lfsr" = lfsr)
save(outlist, file = results_outfile)
message(paste("Saved:", results_outfile))


