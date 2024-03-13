##############################################################
# General functions and required packages
##############################################################

##############################################################
# Packages
suppressPackageStartupMessages(library(data.table    ))
suppressPackageStartupMessages(library(peer          ))
suppressPackageStartupMessages(library(optparse      ))
suppressPackageStartupMessages(library(plyr          ))
suppressPackageStartupMessages(library(preprocessCore))
suppressPackageStartupMessages(library(rhdf5         ))
suppressPackageStartupMessages(library(caret         ))
suppressPackageStartupMessages(library(dplyr         ))
suppressPackageStartupMessages(library(corrplot      ))

set.seed(5366)

##############################################################
# Functions

# add rownames to fread
add_rownames = function(x) 
{
	rownames(x) = x[,1]
	x[,1]       = NULL
	return(x)
}

# add personal theme to ggplot
jn_theme = theme_classic() + theme(axis.text = element_text(size = 20), plot.title = element_text(hjust = 0.5, size = 22), axis.title = element_text(size = 20))

# adjust size of plot
psize = function(w, h) { options(repr.plot.width = w, repr.plot.height = h) }

# blue-white-red color gradient
col2 = c("#67001F", "#B2182B", "#D6604D", "#F4A582", "#FDDBC7", "#FFFFFF", "#D1E5F0", "#92C5DE", "#4393C3", "#2166AC", "#053061")

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

# Step 1 Functions
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
                 "--config_file", config$config_file, 
                 "--out_folder", config$out_folder, 
                 "--functions_file", config$functions_file)
    qsub = paste("echo", paste0("\"", cmd, "\""), "| qsub -N peer -o", logout, "-e", logerr, "-pe smp 8 -V -cwd -l short")
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

process_qtls = function(config, prefix, pp_step, to_return = FALSE, qval_threshold)
{
    qtl_folder       = paste(config$out_folder, pp_step, prefix, "qtl", sep = "/")
    infiles          = list.files(qtl_folder, pattern = "^fdr", full.names = TRUE)
    if (length(infiles) != 0)
    {
        indata           = suppressWarnings(as.data.frame(rbindlist(lapply(infiles, function(x){fread(x, sep = "\t", header = TRUE, data.table = FALSE)})), stringsAsFactors = FALSE))
        outdata          = list()
        element_ids      = unique(indata$element_id)

        for(type in sort(unique(indata$type)))
        {
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

run_step_3 = function(config, peer_to_test, type)
{
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
                process_qtls(config, prefix, pp_step, TRUE, 0.05)
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

select_random_elements = function(element_ids, n = 1000, seed = 123)
{
    set.seed(seed)
    element_ids = sample(element_ids, size = n, replace = FALSE)
    message(paste("Selected", length(element_ids), "random genes"))
    
    return(element_ids)
}
