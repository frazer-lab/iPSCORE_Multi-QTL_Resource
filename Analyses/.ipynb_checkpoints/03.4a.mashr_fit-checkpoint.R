# Author: Jennifer Nguyen

setwd("/projects/CARDIPS/analysis/epigenome_resource")

set.seed(5366)

library(mashr)
library(ash)
library(data.table)
library(dplyr)
library(optparse)

option_list = list(make_option("--date", type = "character", default = NA, help = "date", metavar = "character"))

# parse arguments
opt_parser    = OptionParser(option_list = option_list)
opt           = parse_args(opt_parser)

# set arguments
date = opt$date

npcs = 3

dir = paste("/projects/CARDIPS/analysis/epigenome_resource/analyses/jennifer/ipscore_unique_qtls/mashr", date, sep = "/")

random_datafile    = paste(dir, "input/mashr_input_random.robj", sep = "/")
strong_datafile    = paste(dir, "input/mashr_input_strong.robj", sep = "/")
model_outfile      = paste(dir, "results/mashr_model.robj", sep = "/")
results_outfile    = paste(dir, "results/mashr_results.robj", sep = "/")

dir.create(paste(dir, "results", sep = "/"))

# Open random dataset (used to fit the model)
message(paste("Opening", random_datafile))
load(random_datafile, verbose = T)
random_datalist = datalist

# Open random dataset (used to fit the model)
message(paste("Opening", strong_datafile))
load(strong_datafile, verbose = T)
strong_datalist = datalist

# Function to add rownames 
add_rownames = function(df)
{
  rownames(df) = df[,1]
  df[,1] = NULL
  return(df)
}

# Convert to matrix
random_datalist[["bhat"]] = as.matrix(add_rownames(random_datalist[["bhat"]]))
random_datalist[["shat"]] = as.matrix(add_rownames(random_datalist[["shat"]]))
random_datalist[["shat"]] = random_datalist[["shat"]][rownames(random_datalist[["bhat"]]),]

# Set missing values in random (already done for strong)
random_datalist[["bhat"]][is.na(random_datalist[["bhat"]])] = 1e-6
random_datalist[["shat"]][is.na(random_datalist[["shat"]])] = 1

# Estimate null correlation from random dataset
data.temp = mash_set_data(random_datalist[["bhat"]],random_datalist[["shat"]])
Vhat = estimate_null_correlation_simple(data.temp)
rm(data.temp)

# Prepare random dataset using the null correlation matrix
data.random = mash_set_data(random_datalist[["bhat"]],random_datalist[["shat"]], V=Vhat)

# Prepare strong dataset using the null correlation matrix
strong_datalist[["bhat"]] = as.matrix(add_rownames(strong_datalist[["bhat"]]))
strong_datalist[["shat"]] = as.matrix(add_rownames(strong_datalist[["shat"]]))
strong_datalist[["shat"]] = strong_datalist[["shat"]][rownames(strong_datalist[["bhat"]]),]
data.strong = mash_set_data(strong_datalist[["bhat"]],strong_datalist[["shat"]], V=Vhat)

# Estimate data-driven covariance from strong data
U.pca = cov_pca(data.strong, npcs)
U.ed = cov_ed(data.strong, U.pca)

# Estimate canonical covariance from random data
U.c = cov_canonical(data.random)

# Fit the model
# The outputlevel=1 option means that it will not compute posterior summaries for these tests (which saves time).
message(Sys.time())
m = mash(data.random, Ulist = c(U.ed,U.c), outputlevel = 1)
message(Sys.time())

# Save the model for future runs on more variants
model_list = list("m" = m, "Vhat" = Vhat)
save(model_list, file = model_outfile)
message(paste("Saved:", model_outfile))

## Messages outputted from running this script (from an old run)
#Opening #/projects/CARDIPS/analysis/epigenome_resource/analyses/jennifer/ipscore_unique_qtls/mashr/input/mashr_input_random.robj
#Loading objects:
#  datalist
#Opening #projects/CARDIPS/analysis/epigenome_resource/analyses/jennifer/ipscore_unique_qtls/mashr/input/mash_input_strong.obj
#Loading objects:
# datalist
#2024-06-25 17:22:12.815945
# - Computing 200000 x 2597 likelihood matrix.
# - Likelihood calculations took 25281.18 seconds.
# - Fitting model with 2597 mixture components.
# - Model fitting took 85417.81 seconds.
#2024-06-27 00:07:46.404344
#Saved: #/projects/CARDIPS/analysis/epigenome_resource/analyses/jennifer/ipscore_unique_qtls/mashr/results/mashr_model.robj
