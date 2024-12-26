library(optparse)
library(peer)

option_list = list(make_option("--config_file", type="character", default=0, help="config file", metavar="character"),
                   make_option("--out_folder" , type="character", default=0, help="out folder", metavar="character"),
                   make_option("--functions_file" , type="character", default=0, help="functions file", metavar="character")
				  ) 

opt_parser   = OptionParser(option_list=option_list)
opt          = parse_args(opt_parser)
config_file  = opt$config_file
out_folder   = opt$out_folder
functions_file   = opt$functions_file

message(Sys.time())
source(functions_file)

opt = list(config = config_file, step = 1) # change with optparse
config_file   = opt$config
pipeline_step = opt$step

config = parse_config(config_file)

load(paste(config$out_folder, "step_1", "expdata.robj", sep = "/"), verbose = T)

peerdata = calculate_peer_factors(config, expdata)

