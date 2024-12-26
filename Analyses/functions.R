##############################################################
# General functions and required packages
##############################################################

##############################################################
# Packages
suppressPackageStartupMessages(library(data.table    ))
#suppressPackageStartupMessages(library(peer          ))
suppressPackageStartupMessages(library(optparse      ))
suppressPackageStartupMessages(library(plyr          ))
suppressPackageStartupMessages(library(preprocessCore))
suppressPackageStartupMessages(library(rhdf5         ))
suppressPackageStartupMessages(library(caret         ))
suppressPackageStartupMessages(library(dplyr         ))

set.seed(5366)

jn_theme = theme_classic() + theme(axis.text = element_text(size = 20), plot.title = element_text(hjust = 0.5, size = 22), axis.title = element_text(size = 20))

psize = function(w, h)
{
    options(repr.plot.width = w, repr.plot.height = h) 
}

col2 = rev(c("#67001F", "#B2182B", "#D6604D", "#F4A582", "#FDDBC7", "#FFFFFF", "#D1E5F0", "#92C5DE", "#4393C3", "#2166AC", "#053061"))


##############################################################
# Functions
add_rownames = function(x) # add rownames to fread
{
	rownames(x) = x[,1]
	x[,1]       = NULL
	return(x)
}

# Configuration
read_config_params = function(x)
{
    x = gsub(" ", "", unlist(strsplit(unlist(strsplit(x, "#"))[[1]], "=")))
    
    return(x)
}

parse_config = function(config_file)
{
    indata = readLines(config_file)
    indata = indata[ grepl("^#", indata) == FALSE & indata != ""]
    params = lapply(indata, read_config_params)
    
    param_list = list()
    
    for(x in params)
    {
        param_list[[x[[1]]]] = x[(2:length(x))]
    }
    
    if(param_list$python                  == ""){param_list$python                  = "python"  }
    if(param_list$bcftools                == ""){param_list$bcftools                = "bcftools"}
    if(param_list$out_folder              == ""){param_list$out_folder              = getwd()}
    if(param_list$peer_factors_n_peer     == ""){param_list$peer_factors_n_peer     = nrow(fread(param_list$metadata_sample, select = 1L))}
    if(param_list$peer_factors_n_elements == ""){param_list$peer_factors_n_elements = 5000   }
    if(param_list$qtl_distance            == ""){param_list$qtl_distance            = 2000   }
    if(param_list$maf_threshold           == ""){param_list$maf_threshold           =    0.05}
	if(param_list$chromsizes_file         == "")
	{
		command = paste("wget", "http://hgdownload.cse.ucsc.edu/goldenpath/hg19/bigZips/hg19.chrom.sizes")

		system(command)

		chromsizes_file = paste(config$out_folder, "input", "chrom.sizes", sep = "/")
		invisible(file.copy ("hg19.chrom.sizes", chromsizes_file))
		invisible(file.remove("hg19.chrom.sizes"))
		
		param_list$chromsizes_file = chromsizes_file
	}
	
	param_list$peer_factors_n_peer     = as.numeric(param_list$peer_factors_n_peer    )
	param_list$peer_factors_n_elements = as.numeric(param_list$peer_factors_n_elements)
	param_list$qtl_distance            = as.numeric(param_list$qtl_distance           )
	param_list$maf_threshold           = as.numeric(param_list$maf_threshold          )
	param_list$input_vcfs              = unlist(strsplit(param_list$input_vcfs, ","))
    
    return(param_list)
}

create_workspace = function(config)
{
    message(paste("Creating input directories for:", config$out_folder))
    dir.create(config$out_folder, showWarnings = FALSE)
    dir.create(paste(config$out_folder, "input", sep = "/"), showWarnings = FALSE)
    
    invisible(lapply(1:5, function(x){dir.create(paste(config$out_folder, paste("step", x, sep = "_"), sep = "/"), showWarnings = FALSE)}))
	
    dir.create(paste(config$out_folder, "tmp"                              , sep = "/"), showWarnings = FALSE)
    dir.create(paste(config$out_folder, "logs"                             , sep = "/"), showWarnings = FALSE)
    dir.create(paste(config$out_folder, "script"                           , sep = "/"), showWarnings = FALSE)
    dir.create(paste(config$out_folder, "step_1", "phenotype"              , sep = "/"), showWarnings = FALSE)
    dir.create(paste(config$out_folder, "step_1", "genotype"               , sep = "/"), showWarnings = FALSE)
    dir.create(paste(config$out_folder, "step_1", "phenotype", "by_element", sep = "/"), showWarnings = FALSE)
    dir.create(paste(config$out_folder, "step_1", "genotype" , "by_element", sep = "/"), showWarnings = FALSE)
    dir.create(paste(config$out_folder, "step_2", "test_by_peer"           , sep = "/"), showWarnings = FALSE)
	dir.create(paste(config$out_folder, "step_4", "qtl_by_element"         , sep = "/"), showWarnings = FALSE)
	dir.create(paste(config$out_folder, "step_4", "input"                  , sep = "/"), showWarnings = FALSE)
}

# Run steps

# Prepare input and calculate PEER factor
run_step_1 = function(config)
{
	# Phenotype
	expdata  = normalize_tpm(config)
    
    # Save phenotype data
    save(expdata, file = paste(config$out_folder, "step_1", "expdata.robj", sep = "/"))
    writeLines(expdata$element_ids  , paste(config$out_folder, "step_1", "phenotype", "elements_ids.txt" , sep = "/"), sep = "\n")
    writeLines(expdata$phenotype_ids, paste(config$out_folder, "step_1", "phenotype", "phenotype_ids.txt", sep = "/"), sep = "\n")
    
    # Calculate peer factors
    qsub_peer_factors(config)
    
    # Covariates and metadata 
	qsub_process_covariates(config)

    # Divide phenotype data
	#qsub_divide_phenotypes_by_element(config, expdata, 50)

	# Genotype
	#n_elements = length(expdata$element_ids)
	#run_qsub_genotype(n_elements, config)
}

# Pilot eQTLs using different # of PEERs
run_step_2 = function(config, nelements)
{
	element_ids_file = paste(config$out_folder, "step_2", "elements_ids.random.txt" , sep = "/")
    element_ids      = select_random_elements(readLines(paste(config$out_folder, "step_1", "phenotype", "elements_ids.txt" , sep = "/")), n = 500, seed = 5366)

    element_ids_file = paste(config$out_folder, "step_2", "element_ids.txt", sep = "/")
    writeLines(element_ids, element_ids_file, sep = "\n")
    
	peer_to_test = ceiling((0:(config$peer_factors_n_peer/10)) * 10)
    message(paste("Peers to test:", paste(peer_to_test, collapse = " ")))

    for(peer_n in peer_to_test)
    {
        prefix = paste("test_by_peer/peer", peer_n, sep = ".")
        prepare_qtl_input(config, "step_2", prefix, element_ids, peer_n)
        run_qsub_qtl     (config, "step_2", prefix, element_ids, peer_n)
    }
}

run_step_2.filter_gt = function(config, nelements)
{
	element_ids_file = paste(config$out_folder, "step_2", "elements_ids.random.txt" , sep = "/")
    element_ids      = select_random_elements(readLines(paste(config$out_folder, "step_1", "phenotype", "elements_ids.txt" , sep = "/")), n = 500, seed = 5366)

    element_ids_file = paste(config$out_folder, "step_2", "element_ids.txt", sep = "/")
    #writeLines(element_ids, element_ids_file, sep = "\n")
    
	peer_to_test = ceiling((0:(config$peer_factors_n_peer/10)) * 10)
    message(paste("Peers to test:", paste(peer_to_test, collapse = " ")))

    for(peer_n in peer_to_test)
    {
        prefix = paste("test_by_peer/peer", peer_n, sep = ".")
        prepare_qtl_input(config, "step_2", prefix, element_ids, peer_n)
        run_qsub_qtl2    (config, "step_2", prefix, element_ids, peer_n)
    }
}

run_step_3 = function(config, peer_to_test, type, qval_threshold = 0.05)
{
    #peer_to_test = ceiling((0:(config$peer_factors_n_peer/10)) * 10)
    #peer_to_test = peer_to_test[which(peer_to_test != config$peer_factors_n_peer)]
    pp_step      = "step_2"
    peer         = data.frame(name = paste(type, peer_to_test, sep = "_"))
    peer[,type]  = peer_to_test
    
    qtl2peer = lapply(peer[,type], function(peer_n)
    {
        prefix = paste("test_by_peer", paste(type, peer_n, sep = "."), sep = "/")
        qtl_folder       = paste(config$out_folder, pp_step, prefix, "qtl", sep = "/")
        outfile = paste(qtl_folder, "txt", sep = ".")

        if (file.exists(outfile) == F)
        {
            if (length(list.files(qtl_folder)) != 0)
            {
                message(paste(peer_n, "processing"), appendLF = F)
                process_qtls(config = config, 
                             prefix = prefix, 
                             pp_step = pp_step, 
                             to_return = TRUE, 
                             qval_threshold = qval_threshold, 
                             nconds_to_use = 0)
            } else
            {
                message(paste(peer_n, "empty"), appendLF = F)
            }
        } else
        {
            message(paste(peer_n, "done"), appendLF = F)
            return(fread(outfile, data.table = F))
        }
    })
    
    names(qtl2peer) = peer$name
    
    process_by_peer = function(qtl)
    {
        out = data.frame(n_elements = nrow(qtl[qtl$type == 0,]), n_egenes = nrow(qtl[ qtl$egene == TRUE & qtl$type == 0,]))
    }

    peer = cbind(peer, as.data.frame(rbindlist(lapply(peer$name, function(x){process_by_peer(qtl2peer[[x]])}))))
    peer$pct_egenes = peer$n_egenes / peer$n_elements * 100

    fwrite(peer, paste(config$out_folder, "step_3", paste0(type, "_test.txt"), sep = "/"), sep = "\t", col.names = TRUE, row.names = FALSE)

    top                   = peer[ which.max(peer$n_egenes), type]
    peer_n_optimized_file = paste(config$out_folder, "step_3", paste0(type, "_n_optimized.txt"), sep = "/")

    message(paste("Optimized", type, ":", top))

    writeLines(as.character(top), peer_n_optimized_file, sep = "\n")
    
    message(paste("Saved:", paste(config$out_folder, "step_3", paste0(type, "_test.txt"), sep = "/")))
    message(paste("Saved:", peer_n_optimized_file))

}

run_step_4 = function(config)
{
	element_ids_file = paste(config$out_folder, "step_1", "phenotype", "elements_ids.txt" , sep = "/")
	element_ids      = readLines(element_ids_file)
	prefix           = "qtl_by_element"
	peer_n_optimized = readLines(paste(config$out_folder, "step_3", "peer_n_optimized.txt", sep = "/"))
	
	prepare_qtl_input(config, "step_4", prefix, element_ids, peer_n_optimized)
	run_qsub_qtl2    (config, "step_4", prefix, element_ids, peer_n_optimized)
}

run_step_5 = function(config)
{

}

# Prepare phenotype data
my.invnorm = function(x)
{
    res = rank(x)
    res = qnorm(res/(length(res)+0.5))
    return(res)
}

transform_standard_normal = function(df)
{
    message("Quantile normalizing across genes per sample")
    df                           = as.matrix(df)
    data_valid_expressed_full_qn = normalize.quantiles(df, copy=FALSE)
    rownames(data_valid_expressed_full_qn) = rownames(df)
    colnames(data_valid_expressed_full_qn) = colnames(df)
    
    message("Inverse normalizing across samples per gene")
    input_mat                    = as.data.frame(t(apply(t(data_valid_expressed_full_qn), 2, my.invnorm)))
    
    return(input_mat)
}

normalize_exp = function(config)
{
    library(edgeR)
    counts = add_rownames(fread(config$input_phenotype_count_matrix, data.table = F))
    tpm = add_rownames(fread(config$input_phenotype_tpm_matrix, data.table = F))

    message("TMM normalize across all genes")
    # https://www.biostars.org/p/317701/
    dge = DGEList(counts=as.matrix(counts), group=colnames(counts))
    dge = calcNormFactors(dge, method = "TMM")
    tmm = cpm(dge)

    message(paste("Filtering for expressed genes"))
    message(paste("Require:", ">=", config$expressed_counts, "counts (unnormalized)", "in >=", paste0(as.double(config$expressed_pct) * 100, "%"), "samples"))
    message(paste("Require:", ">=", config$expressed_tpm, "TPM", "in >=", paste0(as.double(config$expressed_pct) * 100, "%"), "samples"))
    message(paste("Starting:", nrow(counts), nrow(tpm), "elements"))

    expressed = as.matrix(tpm)
    expressed[as.matrix(tpm) <  as.double(config$expressed_tpm)] = 0
    expressed[as.matrix(tpm) >= as.double(config$expressed_tpm)] = 1

    tpm_f = tpm[rowSums(expressed) >= (as.double(config$expressed_pct) * ncol(tpm)),]

    expressed = as.matrix(counts)
    expressed[as.matrix(counts) <  as.double(config$expressed_counts)] = 0
    expressed[as.matrix(counts) >= as.double(config$expressed_counts)] = 1

    counts_f = counts[rowSums(expressed) >= (as.double(config$expressed_pct) * ncol(counts)),]

    expressed = intersect(rownames(tpm_f), rownames(counts_f))

    tpm_f = tpm_f[expressed,]
    counts_f = counts_f[expressed,]

    message("Removing sex chromosomal genes")
    length(unique(expressed))
    auto_genes = fread("/reference/private/Gencode.v44lift38/gene_info.txt", data.table = F) %>% filter(chrom %in% paste0("chr", c(1:22)))
    expressed = expressed[which(expressed %in% auto_genes$gene_id)]
    length(unique(expressed))

    # filtered based on TPM and counts
    tmm_f = as.matrix(tmm[expressed,]) 

    # quantile normalize and inverse norm transform across samples
    tmm_f_norm = transform_standard_normal(tmm_f) 

    message(paste("Expressed:", length(unique(expressed)), nrow(counts_f[expressed,]), nrow(tpm_f[expressed,]), nrow(tmm_f[expressed,]), "elements"))
    message(paste("Samples:", ncol(tmm_f_norm)))
    
    return(list("all_tmm"          = tmm,
                "all_counts"       = counts,
                "all_tpm"          = tpm,
                "expressed_tmm"    = tmm_f[expressed,], 
                "expressed_tpm"    = tpm_f[expressed,], 
                "expressed_counts" = counts_f[expressed,], 
                "normalized_tmm"   = tmm_f_norm[expressed,], 
                "element_ids"      = expressed, 
                "phenotype_ids"    = colnames(tmm_f_norm)))


}

normalize_tpm = function(config, input_phenotype_matrix, element_ids)
{
    message(paste("Reading phenotype matrix:", config$input_phenotype_matrix))
    input_phenotype_matrix         = add_rownames(suppressWarnings(fread(config$input_phenotype_matrix, sep = "\t", data.table = FALSE)))
    message(paste("Read", nrow(input_phenotype_matrix), "elements"))
    
    message(paste("Reading phenotype info:", config$input_phenotype_info))
	input_phenotype_info           =                               fread(config$input_phenotype_info  , sep = "\t", data.table = FALSE)[,1:4]
	colnames(input_phenotype_info) = c("chrom", "start", "end", "element_id")
	rownames(input_phenotype_info) = input_phenotype_info$element_id
	element_ids                    = input_phenotype_info$element_id
    
	if(is.null(element_ids) == FALSE){input_phenotype_matrix = input_phenotype_matrix[element_ids,]}
    message(paste("Found info for", nrow(input_phenotype_matrix), "elements"))
    
    expressed  = as.matrix(input_phenotype_matrix)
    normalized = transform_standard_normal(expressed)
    
    message(paste("# elements normalized:", nrow(normalized)))
    message(paste("# samples total:" , ncol(normalized)))
    
    return(list(expressed = expressed, normalized = normalized, element_ids = rownames(normalized), phenotype_ids = colnames(normalized)))
}

qsub_peer_factors = function(config)
{
    logout = paste(config$out_folder, "logs", "calculate_peer_factors.out", sep = "/")
    logerr = paste(config$out_folder, "logs", "calculate_peer_factors.err", sep = "/")
    
    message("Executing..")
    cmd  = paste("Rscript", paste(config$script_dir, "calculate_peer_factors.R", sep = "/"), 
                 "--config_file", config$config_file, "--out_folder", config$out_folder, "--functions_file", config$functions_file)
    qsub = paste("echo", paste0("\"", cmd, "\""), "| qsub -N peer -o", logout, "-e", logerr, "-pe smp 8 -V -cwd -l short")
    message(qsub)
    system(qsub)
}

calculate_peer_factors = function(config, expdata)
{
    expressed  = expdata$expressed_tmm
    normalized = expdata$normalized_tmm
    
    # Select the most variable features for PEER factor calculation
    totest     = data.frame(gene_id = rownames(expressed), sd = unlist(apply(expressed, 1, sd)))
    totest     = totest[order(totest$sd, decreasing = TRUE), "gene_id"]
    normalized = normalized[ totest[1:config$peer_factors_n_elements],]
    
    message(paste(Sys.time(), paste("calculating",  config$peer_factors_n_peer, "peer factors using", nrow(normalized), "elements")))
   
    model = PEER()
    PEER_setPhenoMean (model, t(as.matrix(normalized))  )
    PEER_setNk        (model, config$peer_factors_n_peer)
    PEER_update       (model                            )

    factors   = PEER_getX        (model)
    weights   = PEER_getW        (model)
    precision = PEER_getAlpha    (model)
    residuals = PEER_getResiduals(model)
    precision = data.frame(peer = paste("peer", 1:config$peer_factors_n_peer, sep = "") , precision = precision, stringsAsFactors = FALSE)
    residuals = t(residuals)

    rownames(factors  ) = colnames(normalized)
    rownames(weights  ) = rownames(normalized)
    rownames(residuals) = rownames(normalized)
    colnames(residuals) = colnames(normalized)
    colnames(factors  ) = paste("peer", 1:config$peer_factors_n_peer, sep = "")
    colnames(weights  ) = paste("peer", 1:config$peer_factors_n_peer, sep = "")

    factors          = as.data.frame(factors  )
    weights          = as.data.frame(weights  )
    residuals        = as.data.frame(residuals)
    
    fwrite(factors  , file = paste(config$out_folder, "step_1", "phenotype", paste("peer", "factors"  , "txt", sep = "."), sep = "/"), sep = "\t", row.names = TRUE , col.names = TRUE)
    fwrite(weights  , file = paste(config$out_folder, "step_1", "phenotype", paste("peer", "weights"  , "txt", sep = "."), sep = "/"), sep = "\t", row.names = TRUE , col.names = TRUE)
    fwrite(precision, file = paste(config$out_folder, "step_1", "phenotype", paste("peer", "precision", "txt", sep = "."), sep = "/"), sep = "\t", row.names = TRUE , col.names = TRUE)
    fwrite(residuals, file = paste(config$out_folder, "step_1", "phenotype", paste("peer", "residuals", "txt", sep = "."), sep = "/"), sep = "\t", row.names = TRUE , col.names = TRUE)
    
    save(model, file = paste(config$out_folder, "step_1", "phenotype", paste("peer", "model", "robj", sep = "."), sep = "/"))
    
    message(paste(Sys.time(), "saving"))
                    
    message(paste("Saved:", paste(config$out_folder, "step_1", "phenotype", paste("peer", "factors"  , "txt", sep = "."), sep = "/")))
    message(paste("Saved:", paste(config$out_folder, "step_1", "phenotype", paste("peer", "weights"  , "txt", sep = "."), sep = "/")))
    message(paste("Saved:", paste(config$out_folder, "step_1", "phenotype", paste("peer", "precision", "txt", sep = "."), sep = "/")))
    message(paste("Saved:", paste(config$out_folder, "step_1", "phenotype", paste("peer", "residuals", "txt", sep = "."), sep = "/")))
    message(paste("Saved:", paste(config$out_folder, "step_1", "phenotype", paste("peer", "model"    , "robj", sep = "."), sep = "/")))
    
    return(factors)
}
            
qsub_divide_phenotypes_by_element = function(config, expdata, n)
{
    n_elements = length(expdata$element_ids)
    logout = paste(config$out_folder, "logs", "divide_phenotypes_by_element.out", sep = "/")
    logerr = paste(config$out_folder, "logs", "divide_phenotypes_by_element.err", sep = "/")
    qsub = paste("qsub -N divide_pheno -o", logout, "-e", logerr, "-V -cwd -pe smp 1", "-tc", n, paste0("-t 1-", n_elements, ":1"),
                 paste(config$script_dir, "run_divide_phenotypes_by_element.sh", sep = "/"), config$config_file, config$out_folder, config$functions_file)
    message(paste("Executing:", qsub))
    system(qsub)
}

divide_phenotypes_by_element = function(config, expdata)
{
    invisible(lapply(expdata$element_ids, function(element)
    {
        outdata = data.frame(sample_id = expdata$phenotype_ids, raw = as.numeric(expdata$expressed[ element, expdata$phenotype_ids]), norm = as.numeric(expdata$normalized[ element, expdata$phenotype_ids]))
        fwrite(outdata, paste(config$out_folder, "step_1", "phenotype", "by_element", paste(element, "txt", sep = "."), sep = "/"), sep = "\t", row.names = FALSE, col.names = TRUE)
    }))
}

# Process covariates
scale_metadata = function(metadata_sample)
{
    to_scale                   = colnames(metadata_sample)[ !colnames(metadata_sample) %in% c("phenotype_id", "genotype_id", "subject_id")]
    metadata_sample[,to_scale] = scale(metadata_sample[,to_scale])
    
    return(metadata_sample)
}

qsub_process_covariates = function(config)
{
    logout = paste(config$out_folder, "logs", "process_covariates.out", sep = "/")
    logerr = paste(config$out_folder, "logs", "process_covariates.err", sep = "/")
    
    cmd  = paste("Rscript", paste(config$script_dir, "process_covariates.R", sep = "/"), "--config_file", config$config_file, "--out_folder", config$out_folder, "--functions_file", config$functions_file)
    qsub = paste0("echo \"", cmd, "\" | qsub -hold_jid peer -N covariates -V -cwd -o ", logout, " -e ", logerr, " -pe smp 1")
    message(qsub)
    system(qsub)
}

process_covariates = function(config, peerdata)
{
    metadata_sample      = fread(config$metadata_sample , sep = "\t", header = TRUE, data.table = FALSE)
    metadata_subject     = fread(config$metadata_subject, sep = "\t", header = TRUE, data.table = FALSE)
    metadata             = metadata_sample[,c("phenotype_id", "genotype_id", "subject_id")]
    covariates           = merge(scale_metadata(metadata_sample), peerdata, by.x = "phenotype_id", by.y = "row.names")
    covariates           = merge(covariates                     , metadata_subject, by = c("phenotype_id", "genotype_id", "subject_id"))
    rownames(covariates) = covariates$phenotype_id
    rownames(metadata  ) = metadata  $phenotype_id

    message(paste("Number of rows in covariates:", nrow(covariates)))
    message(paste("Number of rows in metadata:", nrow(metadata)))

    covariates[,c("subject_id", "genotype_id")] = NULL

    covariates = covariates[ metadata$phenotype_id,]    

    fwrite(covariates, file = paste(config$out_folder, "step_1", "covariates.txt", sep = "/"), sep = "\t", row.names = FALSE, col.names = TRUE)
    fwrite(metadata  , file = paste(config$out_folder, "step_1", "metadata.txt"  , sep = "/"), sep = "\t", row.names = FALSE, col.names = TRUE)

    message(paste("Saved:", paste(config$out_folder, "step_1", "covariates.txt", sep = "/")))
    message(paste("Saved:", paste(config$out_folder, "step_1", "metadata.txt"  , sep = "/")))
    
    #return(sort(unique(metadata$genotype_id)))
}

# Prepare genotype data
run_qsub_genotype = function(n_elements, config, tc)
{
    logs_out = paste(config$out_folder, "logs"  , "genotype.err", sep = "/")
    logs_err = paste(config$out_folder, "logs"  , "genotype.out", sep = "/")
    sh_file  = paste(config$out_folder, "script", "genotype.sh" , sep = "/")
    
    my_queue = ""
    to_sh    = c("#!/usr/bin/sh",
                 "source /frazer01/home/jennifer/.bash_profile",
                 "export PATH=/home/tarthur/software/R-4.2.1/bin:$PATH",
                 paste("Rscript", 
                       paste(config$script_dir, "prepare_genotype.R", sep = "/"),
                       "--config_file", config$config_file,
                       "--taskid"     , "$SGE_TASK_ID",
                       "--functions_file", config$functions_file
                      )
                )
    
    writeLines(text = to_sh, con = sh_file, sep = "\n\n")
    message(paste(sh_file, "written!"))
    
    if(config$qsub_queue != ""){my_queue = paste("-l", config$qsub_queue)}
    
    qsub_command = paste("qsub",
                         "-V", "-cwd",
                         "-t", paste(1, "-", n_elements, ":1", sep = ""),
                         "-tc", tc, 
                         "-o" , logs_out,
                         "-e" , logs_err,
                         my_queue,
                         sh_file
                        )
    
    message(qsub_command)
    system(qsub_command)
}

run_qsub_genotype_rerun = function(n_elements, config, tc, missing_file)
{
    logs_out = paste(config$out_folder, "logs", "genotype_rerun.err", sep = "/")
    logs_err = paste(config$out_folder, "logs", "genotype_rerun.out", sep = "/")
    sh_file = paste(config$out_folder, "script", "genotype_rerun.sh", sep = "/")
    my_queue = ""
    
    to_sh = c("#!/usr/bin/sh", 
              "source /frazer01/home/jennifer/.bash_profile", 
              "export PATH=/home/tarthur/software/R-4.2.1/bin:$PATH",
              "SGE_TASK_ID=`tail -n +${SGE_TASK_ID} $1 | head -1`",
        paste("Rscript", paste(config$script_dir, "prepare_genotype.R", sep = "/"), 
              "--config_file", config$config_file, 
              "--taskid", "$SGE_TASK_ID",
              "--functions_file", config$functions_file))
    
    writeLines(text = to_sh, con = sh_file, sep = "\n\n")
    message(paste(sh_file, "written!"))
    
    if (config$qsub_queue != "") {
        my_queue = paste("-l", config$qsub_queue)
    }
    
    qsub_command = paste("qsub", "-V", "-cwd", "-t", paste(1, 
        "-", n_elements, ":1", sep = ""), "-tc", tc, "-o", logs_out, 
        "-e", logs_err, my_queue, sh_file, missing_file)
    
    message(qsub_command)
    system(qsub_command)
}

test_chr = function(region, vcf, bcftools)
{
    command = paste("bcftools", "view", "-H", vcf, "|", "head", "-n", 1)
    message(command)
    chrom   = as.character(fread(cmd = command, sep = "\t", header = FALSE))[[1]]
    
    if(grepl("chr", chrom) == TRUE){region = paste0("chr", region)}
    
    message(paste("Region:", region))
    return(region)
}

get_gtdata = function(config, element_id)
{
    out_folder   = paste(config$out_folder, "step_1", "genotype", "by_element", sep = "/")
    chromsizes   = read.table(config$chromsizes_file                                         , sep = "\t", header = FALSE, col.names = c("chrom", "size"))[,1:2]
    metadata     = fread     (paste(config$out_folder, "step_1", "metadata.txt"  , sep = "/"), sep = "\t", header = TRUE , data.table = FALSE)
    pheno_info   = fread     (config$input_phenotype_info                                    , sep = "\t", header = FALSE, data.table = FALSE)[,1:4]
    pheno_info   = pheno_info[ pheno_info[,4] == element_id,]
    chrom        = pheno_info[,1]
    element_from = as.numeric(pheno_info[,2])
    element_to   = as.numeric(pheno_info[,3])

    region_from = max(0                                            , element_from - config$qtl_distance)
    region_to   = min(chromsizes[chromsizes$chrom == chrom, "size"], element_to   + config$qtl_distance)
    region      = paste(sub("chr", "", chrom), ":", region_from, "-", region_to, sep = "")
    geno_ids    = sort(unique(metadata$genotype_id))

    gttable            = as.data.frame(rbindlist(lapply(config$input_vcfs, function(vcf){filter_vcf_by_region(geno_ids, region, config$maf_threshold, vcf, config$bcftools)})), stringsAsFactors = FALSE)
    gt_info            = gttable[, c("chrom", "pos", "ref", "alt", "rsid", "id", "genotyped")]
    gtmatrix           = as.matrix(gttable[, geno_ids])
    gtmatrix           = mapvalues(gtmatrix, from=c("0|0", "0|1", "1|0", "1|1", ".|."), to=c(0, 0.5, 0.5, 1, NA), warn_missing = FALSE)
    rownames(gtmatrix) = gt_info$id

    if (nrow(gt_info) > 0)
    {
        class(gtmatrix) = "numeric" 
        gtmatrix        = as.data.frame(gtmatrix[unlist(apply(gtmatrix, 1, FUN = function(x){sd(x, na.rm = TRUE)})) > 0,])

        if(ncol(gtmatrix) == 1)
        {
            gtmatrix = as.data.frame(t(gtmatrix))
            rownames(gtmatrix) = gt_info$id
        }
        
        #fwrite(gt_info , paste(out_folder, paste("gt_info", element_id, "txt", sep = "."), sep = "/"), sep = "\t", row.names = FALSE, col.names = TRUE)
        #fwrite(gtmatrix, paste(out_folder, paste("gt_data", element_id, "txt", sep = "."), sep = "/"), sep = "\t", row.names = TRUE , col.names = TRUE)

        #message(paste("Saved:", paste(out_folder, paste("gt_info", element_id, "txt", sep = "."), sep = "/")))
        #message(paste("Saved:", paste(out_folder, paste("gt_data", element_id, "txt", sep = "."), sep = "/")))

        return(list("gt_info" = gt_info, "gtmatrix" = gtmatrix))
    }
}

get_gtdata2 = function(config, element_id)
{
    out_folder   = paste(config$out_folder, "step_1", "genotype", "by_element", sep = "/")
    chromsizes   = read.table(config$chromsizes_file                                         , sep = "\t", header = FALSE, col.names = c("chrom", "size"))[,1:2]
    metadata     = fread     (paste(config$out_folder, "step_1", "metadata.txt"  , sep = "/"), sep = "\t", header = TRUE , data.table = FALSE)
    pheno_info   = fread     (config$input_phenotype_info                                    , sep = "\t", header = FALSE, data.table = FALSE)[,1:4]
    pheno_info   = pheno_info[ pheno_info[,4] == element_id,]
    chrom        = pheno_info[,1]
    element_from = pheno_info[,2]
    element_to   = pheno_info[,3]

    region_from = max(0                                            , element_from - config$qtl_distance)
    region_to   = min(chromsizes[chromsizes$chrom == chrom, "size"], element_to   + config$qtl_distance)
    region      = paste(sub("chr", "", chrom), ":", region_from, "-", region_to, sep = "")
    geno_ids    = sort(unique(metadata$genotype_id))

    gttable            = as.data.frame(rbindlist(lapply(config$input_vcfs, function(vcf){filter_vcf_by_region2(geno_ids, region, config$maf_threshold, vcf, config$bcftools)})), stringsAsFactors = FALSE)
    gt_info            = gttable[, c("chrom", "pos", "ref", "alt", "rsid", "id", "af", "genotyped")]
    gtmatrix           = as.matrix(gttable[, geno_ids])
    gtmatrix           = mapvalues(gtmatrix, from=c("0/0", "0/1", "1/0", "1/1", "./."), to=c(0, 0.5, 0.5, 1, NA), warn_missing = FALSE)
    rownames(gtmatrix) = gt_info$id

    if (nrow(gt_info) > 0)
    {
        class(gtmatrix) = "numeric" 
        gtmatrix        = as.data.frame(gtmatrix[unlist(apply(gtmatrix, 1, FUN = function(x){sd(x, na.rm = TRUE)})) > 0,])

        if(ncol(gtmatrix) == 1)
        {
            gtmatrix = as.data.frame(t(gtmatrix))
            rownames(gtmatrix) = gt_info$id
        }

        #fwrite(gt_info , paste(out_folder, paste("gt_info", element_id, "txt", sep = "."), sep = "/"), sep = "\t", row.names = FALSE, col.names = TRUE)
        #fwrite(gtmatrix, paste(out_folder, paste("gt_data", element_id, "txt", sep = "."), sep = "/"), sep = "\t", row.names = TRUE , col.names = TRUE)

        #message(paste("Saved:", paste(out_folder, paste("gt_info", element_id, "txt", sep = "."), sep = "/")))
        #message(paste("Saved:", paste(out_folder, paste("gt_data", element_id, "txt", sep = "."), sep = "/")))

        return(list("gt_info" = gt_info, "gtmatrix" = gtmatrix))
    }

}

filter_vcf_by_region = function(geno_ids, region, maf_threshold, vcf, bcftools)
{
    region   = test_chr(region, vcf, bcftools)
    command  = paste(bcftools, "view", 
                     "-r", region, 
                     vcf, 
                     "|",
                     bcftools, "query", 
                     "-H", 
                     "-f", '"%CHROM\\t%POS\\t%REF\\t%ALT{0}\\t%ID[\\t%GT]\\n"'
                    )
    
    message(command)
    indata           = fread(cmd = command, sep = "\t", header = TRUE, data.table = FALSE)
    
    message(paste("# SNPs extracted:", nrow(indata)))
    
    colnames(indata) = unlist(lapply(colnames(indata), function(x)
    {
        x = unlist(strsplit(unlist(strsplit(x, "]"))[[2]], ":"))[[1]]

        if(x %in% c("CHROM", "POS", "REF", "ALT", "ID", "AF")){x = tolower(x)}
        return(x)
    }))

    to_add      = setdiff(geno_ids, colnames(indata))
    out_columns = c("chrom", "pos", "ref", "alt", "id", "rsid", "genotyped", geno_ids)

    if(nrow(indata) >  0)
    {
        indata[,to_add] = NA

        indata$chrom = paste0("chr", sub("chr", "", indata$chrom))
        indata$rsid  = indata$id
        indata$id    = sub("chr", "", indata$id)

        if(nrow(indata[ grepl("^rs", indata$id) == TRUE | indata$id == ".", ]) == nrow(indata))
        {
            indata$id = paste("VAR", sub("chr", "", indata$chrom), indata$pos, indata$ref, indata$alt, sep = "_")
        }

        if(length(to_add) > 0)
        {
            if(nrow(indata) >  0)
            {
                message(paste("adding missing genotypes for", length(to_add), "subjects"))
                indata[,to_add] = ".|."
            }
        }

        indata$genotyped = apply(indata[, geno_ids], 1, function(x){length(x[ grepl("\\.", x) == FALSE]) / length(x)})
        message(paste("filtering snps with 99% genotyped"))
        message(paste("starting:", length(unique(indata$id)), length(indata$id)))
        indata           = indata[ indata$genotyped >= 0.99,]
        message(paste("ending:", length(unique(indata$id)), length(indata$id)))
        
    }else
    {
        indata = data.frame(matrix(vector(), nrow = 0, ncol = length(out_columns), dimnames=list(c(), out_columns)), stringsAsFactors = FALSE, check.names = FALSE)
    }

    indata = indata[,out_columns]
    
    #region = ifelse(!region %like% "chr", paste0("chr", region), region)
    #command  = paste("bcftools", "view", 
    #             "-r", region, 
    #             "/reference/private/CARDIPS/CARDIPS.GATK.autosome.vcf.gz", 
    #             "|",
    #             "bcftools", "query", 
    #             "-H", 
    #             "-f", '"%CHROM\\t%POS\\t%REF\\t%ALT{0}\\t%ID\\t%AF\\n"'
    #            )

    #message(command)

    #af           = fread(cmd = command, sep = "\t", header = TRUE, data.table = FALSE)
    #colnames(af) = unlist(lapply(colnames(af), function(x)
    #{
    #    x = unlist(strsplit(unlist(strsplit(x, "]"))[[2]], ":"))[[1]]
    #    if(x %in% c("CHROM", "POS", "REF", "ALT", "ID", "AF")){x = tolower(x)}
    #    return(x)
    #})) 

    #af$id = paste("VAR", gsub("chr", "", af$chrom), af$pos, af$ref, af$alt, sep = "_")

    #indata = merge(indata, af[,c("id", "af")], all.x = T, by = "id") %>% relocate(af, .before = genotyped)
    
    return(indata)
}

filter_vcf_by_region2 = function(geno_ids, region, maf_threshold, vcf, bcftools)
{
    region   = test_chr(region, vcf, bcftools)
    command  = paste(bcftools, "view", 
                     "-r", region, 
                     vcf, 
                     "|",
                     bcftools, "query", 
                     "-H", 
                     "-f", '"%CHROM\\t%POS\\t%REF\\t%ALT{0}\\t%ID[\\t%GT]\\n"'
                    )
    
    message(command)
    
    indata           = fread(cmd = command, sep = "\t", header = TRUE, data.table = FALSE)
    colnames(indata) = unlist(lapply(colnames(indata), function(x)
    {
        x = unlist(strsplit(unlist(strsplit(x, "]"))[[2]], ":"))[[1]]

        if(x %in% c("CHROM", "POS", "REF", "ALT", "ID")){x = tolower(x)}
        return(x)
    }))

    to_add      = setdiff(geno_ids, colnames(indata))
    out_columns = c("chrom", "pos", "ref", "alt", "id", "rsid", "genotyped", geno_ids)

    if(nrow(indata) >  0)
    {
        indata[,to_add] = NA

        indata$chrom = paste0("chr", sub("chr", "", indata$chrom))
        indata$rsid  = indata$id
        indata$id    = sub("chr", "", indata$id)

        if(nrow(indata[ grepl("^rs", indata$id) == TRUE | indata$id == ".", ]) == nrow(indata))
        {
            indata$id = paste("VAR", sub("chr", "", indata$chrom), indata$pos, indata$ref, indata$alt, sep = "_")
        }

        if(length(to_add) > 0)
        {
            if(nrow(indata) >  0)
            {
                indata[,to_add] = "./."
            }
        }

        a = data.frame(t(indata[1, geno_ids]))
        colnames(a) = "gt"
        print(table(a$gt))
        
        indata$genotyped = apply(indata[, geno_ids], 1, function(x){length(x[ grepl("\\.", x) == FALSE]) / length(x)}) 
        indata$hom_ref = apply(indata[, geno_ids], 1, function(x){length(x[ grepl("0/0", x) == T]) / length(x)})
        indata$het = apply(indata[, geno_ids], 1, function(x){length(x[ grepl("1/0", x) == T | grepl("0/1",x) == T]) / length(x)})
        indata$hom_alt = apply(indata[, geno_ids], 1, function(x){length(x[ grepl("1/1", x) == T]) / length(x)})
        
        print(summary(indata$hom_ref))
        print(summary(indata$hom_ref))
        
        print(head(indata[1,c("hom_ref", "het", "hom_alt")]))
        
        message(paste("Starting SNPs:", nrow(indata)))
        indata           = indata[ indata$genotyped >= 0.95 & indata$hom_ref >= 0.07 & indata$hom_alt >= 0.07 & indata$het >= 0.07,]
        message(paste("SNPs after filtering by GT:", nrow(indata)))
        
    }else
    {
        indata = data.frame(matrix(vector(), nrow = 0, ncol = length(out_columns), dimnames=list(c(), out_columns)), stringsAsFactors = FALSE, check.names = FALSE)
    }

    indata = indata[,out_columns]
    
    region = ifelse(!region %like% "chr", paste0("chr", region), region)
    command  = paste("bcftools", "view", 
                 "-r", region, 
                 "/frazer01/projects/CARDIPS/analysis/epigenome_qtls_all_ipscore/pipeline/QTL/general_inputs/snv.vcf.gz", 
                 "|",
                 "bcftools", "query", 
                 "-H", 
                 "-f", '"%CHROM\\t%POS\\t%REF\\t%ALT{0}\\t%ID\\t%AF\\n"'
                )

    message(command)

    af           = fread(cmd = command, sep = "\t", header = TRUE, data.table = FALSE)
    colnames(af) = unlist(lapply(colnames(af), function(x)
    {
        x = unlist(strsplit(unlist(strsplit(x, "]"))[[2]], ":"))[[1]]

        if(x %in% c("CHROM", "POS", "REF", "ALT", "ID", "AF")){x = tolower(x)}
        return(x)
    })) 

    af$id = paste("VAR", gsub("chr", "", af$chrom), af$pos, af$ref, af$alt, sep = "_")

    indata = merge(indata, af[,c("id", "af")], all.x = T, by = "id") %>% relocate(af, .before = genotyped)
    
    return(indata)
}

# Test PEER factors
select_random_elements = function(element_ids, n = 1000, seed = 123)
{
    set.seed(seed)
    element_ids = sample(element_ids, size = n, replace = FALSE)
    message(paste("Selected", length(element_ids), "random genes"))
    
    return(element_ids)
}

qsub_prepare_qtl_input = function(config, my_step, prefix, peer_n, element_ids_file)
{
    logout = paste(config$out_folder, "logs", paste("prepare_qtl_input", my_step, prefix, "out", sep = "."), sep = "/")
    logerr = paste(config$out_folder, "logs", paste("prepare_qtl_input", my_step, prefix, "err", sep = "."), sep = "/")
    
    cmd = paste("Rscript prepare_qtl_input.R", 
                "--config_file", config$config_file, 
                "--step", my_step, 
                "--prefix", prefix, 
                "--peer_n", peer_n,
                "--element_ids_file", element_ids_file)
    
    qsub = paste("echo", paste0("\"", cmd, "\""), "|", "qsub -hold_jid peer,covariates,prepare_vcf.sh -o", logout, "-e", logerr, "-pe smp 1", "-N prepare_qtl_input")
    message(qsub)
    system(qsub)
}

# Run QTLs
prepare_qtl_input_pca = function(config, my_step, prefix, element_ids, peer_n, covs)
{
    message(paste("Testing with", pc_n, "pcs factors:", prefix, my_step))

    phenotype_ids         =              readLines(paste(config$out_folder, "step_1", "phenotype", "phenotype_ids.txt", sep = "/"))
    covariates            =              fread    (paste(config$out_folder, "step_1", "covariates.txt"               , sep = "/"), sep = "\t", header = TRUE , data.table = FALSE)
    metadata              =              fread    (paste(config$out_folder, "step_1", "metadata.txt"                 , sep = "/"), sep = "\t", header = TRUE , data.table = FALSE)
    phenotype_info        =              fread    (config$input_phenotype_info                                                   , sep = "\t", header = FALSE, data.table = FALSE)[,1:4]
    kinship               =              fread    (config$kinship_file                                                           , sep = "\t", header = FALSE, data.table = FALSE)
    kinship_ids           = as.character(fread    (paste(config$kinship_file, "id", sep = ".")                                   , sep = "\t", header = FALSE, data.table = FALSE)[,1])
    phenotype_pcs         =              fread    (paste(config$out_folder, "input", "phenotype_pcs.txt", sep = "/"), data.table = F)

    if(my_step == "step_4")
    {
        pc_n_optimized_file = paste(config$out_folder, "step_3", "pc_n_optimized.txt", sep = "/")

        if(file.exists(pc_n_optimized_file) == TRUE)
        {
            pc_n = as.numeric(readLines(pc_n_optimized_file))[[1]]
        }
    }

    if (pc_n == 0)
    {
        covariates = covariates[,covs]
    } else
    {
        covariates = merge(covariates[,c("phenotype_id", covs)], phenotype_pcs[,c("phenotype_id", paste0("PC", c(1:pc_n)))], by = "phenotype_id")
    }

    covariates = add_rownames(merge(metadata, covariates, by = "phenotype_id"))

    colnames(phenotype_info) = c("chrom", "start", "end", "element_id")
    rownames(phenotype_info) = phenotype_info$element_id
    rownames(kinship       ) = kinship_ids
    colnames(kinship       ) = kinship_ids
    kinship                  = kinship[ covariates$genotype_id, covariates$genotype_id]
    rownames(kinship       ) = rownames(covariates)
    colnames(kinship       ) = rownames(covariates)
    rownames(metadata      ) = metadata$phenotype_id
    metadata                 = metadata[ rownames(covariates),]
    phenotype_info           = phenotype_info[ element_ids,]
    phenotype_info           = phenotype_info[ order(phenotype_info$chrom, phenotype_info$start),]

    covariates[,c("genotype_id", "subject_id")] = NULL

    dir.create(paste(config$out_folder, my_step, prefix       , sep = "/"), showWarnings = FALSE)
    dir.create(paste(config$out_folder, my_step, prefix, "qtl", sep = "/"), showWarnings = FALSE)
    dir.create(paste(config$out_folder, my_step, prefix, "tmp", sep = "/"), showWarnings = FALSE)

    saveRDS(list(covariates = covariates, metadata = metadata, element_ids = element_ids, phenotype_ids = rownames(covariates), phenotype_info = phenotype_info), paste(config$out_folder, my_step, prefix, "qtl_input.rds", sep = "/"))

    fwrite(covariates, paste(config$out_folder, my_step, prefix, "covariates.csv", sep = "/"), sep = ",", col.names = TRUE, row.names = TRUE)
    fwrite(kinship   , paste(config$out_folder, my_step, prefix, "kinship.csv"   , sep = "/"), sep = ",", col.names = TRUE, row.names = TRUE)

    message(paste("Saved:", paste(config$out_folder, my_step, prefix, "covariates.csv", sep = "/")), appendLF = F)
    message(paste("Saved:", paste(config$out_folder, my_step, prefix, "kinship.csv"   , sep = "/")), appendLF = F)
}

prepare_qtl_input = function(config, my_step, prefix, element_ids, peer_n)
{
    message(paste("Testing with", peer_n, "peer factors:", prefix, my_step))

    phenotype_ids         =              readLines(paste(config$out_folder, "step_1", "phenotype", "phenotype_ids.txt", sep = "/"))
    covariates            =              fread    (paste(config$out_folder, "step_1", "covariates.txt"               , sep = "/"), sep = "\t", header = TRUE , data.table = FALSE)
    metadata              =              fread    (paste(config$out_folder, "step_1", "metadata.txt"                 , sep = "/"), sep = "\t", header = TRUE , data.table = FALSE)
    phenotype_info        =              fread    (config$input_phenotype_info                                                   , sep = "\t", header = FALSE, data.table = FALSE)[,1:4]
    kinship               =              fread    (config$kinship_file                                                           , sep = "\t", header = FALSE, data.table = FALSE)
    kinship_ids           = as.character(fread    (paste(config$kinship_file, "id", sep = ".")                                   , sep = "\t", header = FALSE, data.table = FALSE)[,1])

    if(my_step == "step_4")
    {
        peer_n_optimized_file = paste(config$out_folder, "step_3", "peer_n_optimized.txt", sep = "/")

        if(file.exists(peer_n_optimized_file) == TRUE)
        {
            peer_n = as.numeric(readLines(peer_n_optimized_file))[[1]]
        }
    }

    if(peer_n < config$peer_factors_n_peer){covariates = covariates[, !colnames(covariates) %in% paste0("peer", (peer_n + 1):config$peer_factors_n_peer)]}

    colnames(phenotype_info) = c("chrom", "start", "end", "element_id")
    rownames(phenotype_info) = phenotype_info$element_id
    rownames(kinship       ) = kinship_ids
    colnames(kinship       ) = kinship_ids
    covariates               = add_rownames(merge(covariates, metadata))[ phenotype_ids,]
    kinship                  = kinship[ covariates$genotype_id, covariates$genotype_id]
    rownames(kinship       ) = rownames(covariates)
    colnames(kinship       ) = rownames(covariates)
    rownames(metadata      ) = metadata$phenotype_id
    metadata                 = metadata[ rownames(covariates),]
    phenotype_info           = phenotype_info[ element_ids,]
    phenotype_info           = phenotype_info[ order(phenotype_info$chrom, phenotype_info$start),]

    covariates[,c("genotype_id", "subject_id")] = NULL

    dir.create(paste(config$out_folder, my_step, prefix       , sep = "/"), showWarnings = FALSE)
    dir.create(paste(config$out_folder, my_step, prefix, "qtl", sep = "/"), showWarnings = FALSE)
    dir.create(paste(config$out_folder, my_step, prefix, "tmp", sep = "/"), showWarnings = FALSE)

    saveRDS(list(covariates = covariates, metadata = metadata, element_ids = element_ids, phenotype_ids = rownames(covariates), phenotype_info = phenotype_info), paste(config$out_folder, my_step, prefix, "qtl_input.rds", sep = "/"))

    fwrite(covariates, paste(config$out_folder, my_step, prefix, "covariates.csv", sep = "/"), sep = ",", col.names = TRUE, row.names = TRUE)
    fwrite(kinship   , paste(config$out_folder, my_step, prefix, "kinship.csv"   , sep = "/"), sep = ",", col.names = TRUE, row.names = TRUE)

    message(paste("Saved:", paste(config$out_folder, my_step, prefix, "covariates.csv", sep = "/")), appendLF = F)
    message(paste("Saved:", paste(config$out_folder, my_step, prefix, "kinship.csv"   , sep = "/")), appendLF = F)
}


run_qsub_qtl = function(config, my_step, prefix, element_ids, peer_n, type, run_conditionals, run_qsub)
{
    logs_out = paste(config$out_folder, "logs"  , paste("qtl", my_step, type, peer_n, "err", sep = "."), sep = "/")
    logs_err = paste(config$out_folder, "logs"  , paste("qtl", my_step, type, peer_n, "out", sep = "."), sep = "/")
    sh_file  = paste(config$out_folder, "script", paste("qtl", my_step, type, peer_n, "sh" , sep = "."), sep = "/")
    nohup_file  = paste(config$out_folder, "script", paste("qtl", my_step, type, peer_n, "nohup.sh" , sep = "."), sep = "/")
    my_queue = ""
    to_sh    = c("#!/usr/bin/sh",
                 "source /frazer01/home/jennifer/.bash_profile",
                 "export PATH=/home/tarthur/software/R-4.2.1/bin:$PATH",
                 paste("Rscript", 
                       paste(config$script_dir, "run_qtl.by_element.R", sep = "/"),
                       "--config_file", config$config_file,
                       "--step"       , my_step, 
                       "--prefix"     , prefix,
                       "--taskid"     , "$SGE_TASK_ID",
                       "--run_conditionals", run_conditionals,
                       "--functions_file", config$functions_file
                      )
                )
    
    writeLines(text = to_sh, con = sh_file, sep = "\n\n")
    
    to_sh    = c("#!/usr/bin/sh",
                 "source /frazer01/home/jennifer/.bash_profile",
                 "for SGE_TASK_ID in {1..1000}; do",
                 paste("Rscript", 
                       paste(config$script_dir, "run_qtl.by_element.R", sep = "/"),
                       "--config_file", config$config_file,
                       "--step"       , my_step, 
                       "--prefix"     , prefix,
                       "--taskid"     , "$SGE_TASK_ID",
                       "--run_conditionals", run_conditionals,
                       "--functions_file", config$functions_file
                      ),
                 "done"
                )
    
    writeLines(text = to_sh, con = nohup_file, sep = "\n\n")
    
    message(paste("Written:", sh_file))
    message(paste("Written:", nohup_file))
    
    if(config$qsub_queue != ""){my_queue = paste("-l", config$qsub_queue)}
    
    n_elements   = length(element_ids)
    qsub_command = paste("qsub",
	                     "-hold_jid", "genotype.sh,genotype_redo.sh,genotype_redo,genotype_rerun.sh,genotype_rerun",
                         "-V", "-cwd",
                         "-t", paste(1, "-", n_elements, ":1", sep = ""),
                         "-tc", 300, 
                         "-o" , logs_out,
                         "-e" , logs_err,
                         "-pe smp 2",
                         my_queue,
                         sh_file
                        )
    
    message(qsub_command)
    
    if (run_qsub == T) { system(qsub_command) } else { message("Not running qsub") }
}

write_h5_file = function(element_id, phenotype_info, expdata, gtinfo, gtdata, h5_file)
{
    invisible(suppressWarnings(file.remove(h5_file)))

    invisible(h5createFile (h5_file))
    invisible(h5createGroup(h5_file, "genotype"))
    invisible(h5createGroup(h5_file, "genotype/col_header"))
    invisible(h5createGroup(h5_file, "genotype/row_header"))
    invisible(h5createGroup(h5_file, "phenotype"))
    invisible(h5createGroup(h5_file, "phenotype/col_header"))
    invisible(h5createGroup(h5_file, "phenotype/row_header"))

    h5write(gtinfo$chrom                  , file = h5_file, name="genotype/col_header/chrom"        )
    h5write(gtinfo$pos                    , file = h5_file, name="genotype/col_header/pos"          )
    h5write(gtinfo$pos                    , file = h5_file, name="genotype/col_header/pos_cum"      )
    h5write(as.matrix(gtdata)             , file = h5_file, name="genotype/matrix"                  )
    h5write(colnames (gtdata)             , file = h5_file, name="genotype/row_header/sample_ID"    )
    h5write(c(element_id)                 , file = h5_file, name="phenotype/col_header/gene_ID"     )
    h5write(phenotype_info[, "chrom"     ], file = h5_file, name="phenotype/col_header/gene_chrom"  )
    h5write(phenotype_info[, "end"       ], file = h5_file, name="phenotype/col_header/gene_end"    )
    h5write(phenotype_info[, "start"     ], file = h5_file, name="phenotype/col_header/gene_start"  )
    h5write(c("+")                        , file = h5_file, name="phenotype/col_header/gene_strand" )
    h5write(element_id                    , file = h5_file, name="phenotype/col_header/phenotype_ID")
    h5write(as.matrix(expdata)            , file = h5_file, name="phenotype/matrix"                 )
    h5write(colnames( expdata)            , file = h5_file, name="phenotype/row_header/sample_ID"   )
}

find_redundancy = function(gtdata, cor_threshold = 0.8)
{
    filteredDescr = t(gtdata)
    descrCor      = cor(filteredDescr)
    
    descrCor[ is.na(descrCor) == TRUE] = 1
    
    highlyCorDescr = findCorrelation(descrCor, cutoff = cor_threshold)
    #if(length(highlyCorDescr) == 0){highlyCorDescr = 0}
    
    n_independent  = nrow(gtdata) - length(highlyCorDescr)
    
    return(n_independent)
}

run_eigenmt = function(config, gt_file, geneloc_file, snploc_file, qtl_file, fdr_file, chrom)
{
    command = paste("python", paste(config$script_dir,"eigenMT.py", sep = "/"), 
                    "--CHROM"   , chrom,
                     "--QTL"    , qtl_file,
                     "--GEN"    , gt_file,
                     "--GENPOS" , snploc_file,
                     "--PHEPOS" , geneloc_file,
                     "--OUT"    , fdr_file
                   )
    
	message(paste("running eigenMT:", command))
    system(command)
}

run_qtl = function(config, element_id, tmp_folder, phenotype_info, expdata, gtinfo, gtdata, h5_file, cov_file, kin_file, qtl_file_tmp, geneloc_file, snploc_file, gt_file, lead_vars)
{
    suppressWarnings(write_h5_file(element_id, phenotype_info, expdata, gtinfo, gtdata, h5_file))
	
    command_py = paste(config$python, paste(config$script_dir, "run_limix.py", sep = "/"), tmp_folder, element_id, "scan", "normal")
    
    message(paste("Running limix:", command_py))
    system(command_py)
    
    message(paste("Reading limix output:", qtl_file_tmp))
    
    indata                   = add_rownames(fread(qtl_file_tmp, sep = ",", header = TRUE, data.table = FALSE))
	colnames(indata)         = c("beta", "se", "pval")
    indata                   = cbind(gtinfo, indata)
    indata$element_id        = element_id
    if (length(lead_vars) != 0) { indata = indata[!indata$id %in% unlist(lead_vars),] }
    indata$bonferroni        = p.adjust(indata$pval, method = "bonferroni")
    
    chrom = paste0("chr", unique(unlist(lapply(indata$id, function(x) { unlist(strsplit(x, "_"))[2] }))))
	
	#if(nrow(indata[indata$se > 100, ]) > 0){indata[indata$se > 100, "se"] = 100}
    if(nrow(indata[indata$se > 100 | is.na(indata$se), ]) > 0){indata[indata$se > 100 | is.na(indata$se), "se"] = 100}
	
    if(min(indata$bonferroni) == 1)
    {
        message("Minimum bonferroni is 1")
        fdrdata       = indata[which.min(indata$pval), ]
        fdrdata$fdr   = 1
        fdrdata$tests = nrow(indata)
    }else
    {
        eigenmt_input            = indata[,c("id", "element_id", "beta", "se", "pval", "bonferroni")]
        colnames(eigenmt_input)  = c("SNP", "gene", "beta", "t-stat", "p-value", "FDR")
        qtl_file_emt_in          = sub("qtl.csv", "eigenmt_input.txt" , qtl_file_tmp)
        qtl_file_emt_out         = sub("qtl.csv", "eigenmt_output.txt", qtl_file_tmp)
	
        fwrite(eigenmt_input, qtl_file_emt_in , sep = "\t", row.names = FALSE, col.names = TRUE)
	
        run_eigenmt(config, gt_file, geneloc_file, snploc_file, qtl_file_emt_in, qtl_file_emt_out, chrom)
	
        message(paste("Reading eigenMT output:", qtl_file_emt_out))
        fdrdata           = fread(qtl_file_emt_out, sep = "\t", header = TRUE, data.table = FALSE)
        colnames(fdrdata) = c("id", "element_id", "beta", "se", "pval", "bonferroni", "fdr", "tests")
        
        cmd = paste("rm", qtl_file_emt_in, qtl_file_emt_out, qtl_file_tmp)
        message(cmd)
        system(cmd)
    }
    
    fdrdata           = merge(gtinfo, fdrdata)
    fdrdata           = fdrdata[,c(colnames(indata), "fdr", "tests")]
    
    return(list(lead = fdrdata, qtl = indata))
}

remove_mhc = function(config, mhc_bedfile)
{
    a = mhc_bedfile

    # 1Mb window from genes
    b = paste(config$out_folder, "input", "phenotype_info.bed", sep = "/")
    b = fread(b, data.table = F) %>% mutate(V2 = V2 - config$qtl_distance, V3 = V3 + config$qtl_distance)
    b[b$V2 < 0,]$V2 = 0
    fwrite(b, "analyses/jennifer/scratch/tmp.bed", row.names = F, sep = "\t", col.names = F)
    b = paste(getwd(), "analyses/jennifer/scratch/tmp.bed", sep = "/")

    cmd = paste("bedtools intersect", "-a", a, "-b", b, "-wa -wb")
    message(cmd)
    int = fread(cmd = cmd, data.table = F)
    return(int)
}

recorrect_after_mhc = function(config, prefix, pp_step, to_return, qval_threshold, nconds_to_use, indata)
{
    outdata          = list()
    element_ids      = unique(indata$element_id)

    message(paste("Recorrecting", length(unique(element_ids)), "elements"))

    for(type in sort(unique(indata$type)))
    {
        if (type <= nconds_to_use)
        {
            message(paste("Processing conditional", type), appendLF = F)
            if(length(element_ids) > 0)
            {
                this = indata[indata$type == type & indata$element_id %in% element_ids,]
                this$filt_qval  = p.adjust(this$fdr, method = "BH")
                this$new_egene = FALSE

                this[this$filt_qval <= qval_threshold, "new_egene"] = TRUE

                element_ids         = this[this$new_egene == TRUE, "element_id"]
                outdata[[type + 1]] = this
            }else
            {
                break
            }
        } else
        {
            message(paste("Skipping conditional", type), appendLF = F)
        }
    }

    indata  = as.data.frame(rbindlist(outdata), stringsAsFactors = FALSE)

    prefix  = "qtl_by_element"
    pp_step = "step_4"
    outfile = paste(config$out_folder, pp_step, prefix, "qtl.no_mhc.txt", sep = "/")
    fwrite(indata, outfile, row.names = F, sep = "\t")
    
    message(paste("# qElements:", length(unique(indata[indata$new_egene == T,]$element_id))))

    message(paste("Saved:", outfile))
    
    if(to_return == TRUE){return(indata)}
}

process_qtls = function(config, prefix, pp_step, to_return = FALSE, qval_threshold = 0.05, nconds_to_use)
{
    qtl_folder       = paste(config$out_folder, pp_step, prefix, "qtl", sep = "/")
    infiles          = list.files(qtl_folder, pattern = "^fdr", full.names = TRUE)
    if (length(infiles) != 0)
    {
        indata           = suppressWarnings(as.data.frame(rbindlist(lapply(infiles, function(x){fread(x, sep = "\t", header = TRUE, data.table = FALSE)})), stringsAsFactors = FALSE))
        
        outdata          = list()
        element_ids      = unique(indata$element_id)
        
        message(paste("Processing", length(unique(element_ids)), "elements"), appendLF = F)

        for(type in sort(unique(indata$type)))
        {
            if (type <= nconds_to_use)
            {
                message(paste("Processing conditional", type), appendLF = F)
                if(length(element_ids) > 0)
                {
                    this = indata[indata$type == type & indata$element_id %in% element_ids,]
                    this$qval  = p.adjust(this$fdr, method = "BH")
                    this$egene = FALSE

                    this[this$qval <= qval_threshold, "egene"] = TRUE

                    element_ids         = this[this$egene == TRUE, "element_id"]
                    outdata[[type + 1]] = this
                }else
                {
                    break
                }
            } else
            {
                message(paste("Skipping conditional", type), appendLF = F)
            }
        }

        indata  = as.data.frame(rbindlist(outdata), stringsAsFactors = FALSE)
        outfile = paste(qtl_folder, "txt", sep = ".")
        fwrite(indata, outfile, sep = "\t", col.names = TRUE, row.names = FALSE)
        message(paste("Saved:", outfile), appendLF = T)

        if(to_return == TRUE){return(indata)}
    } else
    {
        message(paste(prefix, "hasn't run yet"))
    }    
}

run_regression_qtl = function(config, element_id, tmp_folder, phenotype_info, expdata, gtinfo, gtdata, h5_file, cov_file, kin_file, qtl_file_tmp, geneloc_file, snploc_file, gt_file, lead_vars)
{
    message("Running QTL regression limix")

    command_py = paste(config$python, paste(config$script_dir, "run_limix.py", sep = "/"), tmp_folder, element_id, "scan", "normal")

    message(command_py)
    system(command_py)

    message(paste("Reading limix output:", qtl_file_tmp))

    indata                   = add_rownames(fread(qtl_file_tmp, sep = ",", header = TRUE, data.table = FALSE))
    colnames(indata)         = c("beta", "se", "pval")
    indata                   = cbind(gtinfo, indata)
    indata$element_id        = element_id
    if (length(lead_vars) != 0) { indata = indata[!indata$id %in% unlist(lead_vars),] }
    indata$bonferroni        = p.adjust(indata$pval, method = "bonferroni")

    chrom = paste0("chr", unique(unlist(lapply(indata$id, function(x) { unlist(strsplit(x, "_"))[2] }))))

    if(nrow(indata[indata$se > 100, ]) > 0){indata[indata$se > 100, "se"] = 100}

    message("Running eigenMT")
    eigenmt_input            = indata[,c("id", "element_id", "beta", "se", "pval", "bonferroni")]
    
    colnames(eigenmt_input)  = c("SNP", "gene", "beta", "t-stat", "p-value", "FDR")
    qtl_file_emt_in          = sub("qtl.csv", "eigenmt_input.txt" , qtl_file_tmp)
    qtl_file_emt_out         = sub("qtl.csv", "eigenmt_output.txt", qtl_file_tmp)

    fwrite(eigenmt_input, qtl_file_emt_in , sep = "\t", row.names = FALSE, col.names = TRUE)
    
    run_eigenmt(config, gt_file, geneloc_file, snploc_file, qtl_file_emt_in, qtl_file_emt_out, chrom)

    message(paste("Reading eigenMT output:", qtl_file_emt_out))
    
    fdrdata           = fread(qtl_file_emt_out, sep = "\t", header = TRUE, data.table = FALSE)

    colnames(fdrdata) = c("id", "element_id", "beta", "se", "pval", "bonferroni", "fdr", "tests")
    
    fdrdata           = merge(gtinfo, fdrdata)
    fdrdata           = fdrdata[,c(colnames(indata), "fdr", "tests")]


    cmd = paste("rm", qtl_file_emt_in, qtl_file_emt_out, qtl_file_tmp)
    message(cmd)
    system(cmd)

    return(list(lead = fdrdata, qtl = indata))
}


run_qsub_regression_qtl = function(config, my_step, prefix, element_ids, peer_n, run_qsub, threads)
{
    logs_out    = paste(config$out_folder, "logs"  , paste("qtl", my_step, "regression_qtl", "peer", peer_n, "err", sep = "."), sep = "/")
    logs_err    = paste(config$out_folder, "logs"  , paste("qtl", my_step, "regression_qtl", "peer", peer_n, "out", sep = "."), sep = "/")
    sh_file     = paste(config$out_folder, "script", paste("qtl", my_step, "regression_qtl", "peer", peer_n, "sh" , sep = "."), sep = "/")
    nohup_file  = paste(config$out_folder, "script", paste("qtl", my_step, "regression_qtl", "peer", peer_n, "nohup.sh" , sep = "."), sep = "/")


    my_queue = ""
    to_sh    = c("#!/usr/bin/sh",
                 "source /frazer01/home/jennifer/.bash_profile",
                 "export PATH=/home/tarthur/software/R-4.2.1/bin:$PATH",
                 paste("Rscript", 
                       paste(config$script_dir, "run_dprime_regression_qtl.by_element.R", sep = "/"),
                       "--config_file", config$config_file,
                       "--step"       , my_step, 
                       "--prefix"     , prefix,
                       "--taskid"     , "$SGE_TASK_ID",
                       "--functions_file", config$functions_file
                      )
                )

    writeLines(text = to_sh, con = sh_file, sep = "\n\n")

    to_sh    = c("#!/usr/bin/sh",
                 "source /frazer01/home/jennifer/.bash_profile",
                 "for SGE_TASK_ID in {1..1000}; do",
                 paste("Rscript", 
                       paste(config$script_dir, "run_dprime_regression_qtl.by_element.R", sep = "/"),
                       "--config_file", config$config_file,
                       "--step"       , my_step, 
                       "--prefix"     , prefix,
                       "--taskid"     , "$SGE_TASK_ID",
                       "--functions_file", config$functions_file
                      ),
                 "done"
                )

    writeLines(text = to_sh, con = nohup_file, sep = "\n\n")

    message(paste("Written:", sh_file))
    message(paste("Written:", nohup_file))

    if(config$qsub_queue != ""){my_queue = paste("-l", config$qsub_queue)}

    n_elements   = length(element_ids)
    qsub_command = paste("qsub",
                         "-hold_jid", "genotype.sh,genotype_redo.sh,genotype_redo,genotype_rerun.sh,genotype_rerun",
                         "-V", "-cwd",
                         "-t", paste(1, "-", n_elements, ":1", sep = ""),
#                          "-tc", 300, 
                         "-o" , logs_out,
                         "-e" , logs_err,
                         "-pe smp", threads,
                         my_queue,
                         sh_file
                        )

    message(qsub_command)

    if (run_qsub == T) { system(qsub_command) } else { message("Not running qsub") }
}


# ------------------- Plotting functions
plot_box_qtl = function(config, element_id, snp, qtls)
{
    meta_subject = fread(config$metadata_subject, data.table = F)

    load(paste(config$out_folder, "step_1", "expdata.robj", sep = "/"))
    expdata = add_rownames(data.frame(sample_id = expdata$phenotype_ids, 
                                      raw = as.numeric(expdata$expressed_tmm[ element_id, expdata$phenotype_ids]), 
                                      norm = as.numeric(expdata$normalized_tmm[ element_id, expdata$phenotype_ids])))
    gt_data = get_gtdata(config, element_id)

    expdata$phenotype_id = rownames(expdata)

    expdata = merge(expdata, meta_subject[,c("phenotype_id", "genotype_id")] %>% distinct(), all.x = T)

    gtinfo = gt_data$gt_info
    gtmatrix = gt_data$gtmatrix

    gtmatrix = melt(gtmatrix[snp,]) %>% dplyr::rename(genotype_id = variable, gt = value)
    
    toplot = merge(expdata, gtmatrix, by = "genotype_id")

    peak = element_id
    beta = signif(qtls[qtls$element_id == peak,]$beta, 2)
    qval = signif(qtls[qtls$element_id == peak,]$qval, 2)

    title = paste(element_id, snp, paste(paste("beta =", beta), paste("q =", qval), sep = "; "), sep = "\n")
    
    toplot = toplot %>% filter(!is.na(gt))

    psize(w = 6, h = 7)
    ggplot(toplot, aes(x = factor(gt, levels = c("0", "0.5", "1")), y = norm)) + 
        geom_boxplot() + geom_jitter(width = 0.15) + jn_theme + 
        geom_smooth(aes(x = factor(gt), y = norm, group = 1), method = "lm", color = "blue") +
        xlab("") + ylab("Normalized TMM") + 
        ggtitle(title)
}
