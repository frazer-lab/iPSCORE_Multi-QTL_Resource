setwd("/projects/CARDIPS/analysis/epigenome_resource")
suppressMessages(source("analyses/jennifer/notebooks/functions.R"))
set.seed(5366)
library(coloc)

option_list = list(make_option("--taskid", type = "integer", default = NA, help = "taskid", metavar = "integer"),
                   make_option("--input_file", type = "character", default = NA, help = "analysis", metavar = "character"),
                   make_option("--manifest", type = "character", default = NA, help = "manifest file", metavar = "character"),
                   make_option("--outdir", type = "character", default = NA, help = "${out_dir}/${analysis}", metavar = "character"))

# parse arguments
opt_parser    = OptionParser(option_list = option_list)
opt           = parse_args(opt_parser)

# set arguments
task_id       = opt$taskid
input_file    = opt$input_file # list of qtls to colocalization GWAS with (require columns: type, element_id, qtl_id, tissue, analysis)
manifest_file = opt$manifest
outdir        = opt$outdir

message(paste("Input file:"   , input_file   ))
message(paste("Manifest file:", manifest_file))
message(paste("Taskid:"       , task_id      ))

# read gwas manifest
manifest = fread(manifest_file, data.table = F)

# read list of qtls to process
# filter for the qtl to process
input    = fread(input_file   , data.table = F)
input    = input[input$taskid == task_id,]

# get info about the qtl
element_id      = input$element_id
discovery_order = input$type
analysis        = input$datatype
tissue          = input$tissue
qtl_file        = input$path

message(paste("Element ID:"     , element_id     ))
message(paste("Discovery Order:", discovery_order))
message(paste("Analysis:"       , analysis       ))
message(paste("Tissue:"         , tissue         ))
message(paste("QTL file:"       , qtl_file       ))

suppressWarnings(dir.create(paste(outdir, sep = "/")))
suppressWarnings(dir.create(paste(outdir, analysis, sep = "/")))
suppressWarnings(dir.create(paste(outdir, analysis, tissue, sep = "/")))

# set file path for output
outdir = paste(outdir, analysis, tissue, sep = "/")
outfile = paste(outdir, paste(paste(discovery_order, element_id, sep = "-"), "robj", sep = "."), sep = "/")

# read phenotype info
pheninfo = fread(paste(analysis, tissue, "input/phenotype_info.txt", sep = "/"), data.table = F)
pheninfo = pheninfo[pheninfo[,4] == element_id,]

# set sample size for each analysis type
# number corresponds to the number of samples
ipscore_sample_size = list("iPSC eqtls"  = 220, 
                           "iPSC caqtls" = 142, 
                           "iPSC haqtls" = 43, 
                           "CVPC eqtls"  = 178, 
                           "CVPC caqtls" = 140, 
                           "CVPC haqtls" = 101, 
                           "PPC eqtls"   = 107, 
                           "PPC caqtls"  = 109)

# read qtldata
qtldata        = fread(qtl_file, data.table = F)
qtldata$chrpos = paste(qtldata$chrom, qtldata$pos, sep = "_") # set chrom and pos
qtldata        = qtldata %>% dplyr::rename(a1 = ref, a2 = alt, p = pval) # rename columns to match those in gwas
qtldata$maf    = ifelse(qtldata$af > 0.5, 1-qtldata$af, qtldata$af) # set maf
qtldata        = qtldata[qtldata$type == discovery_order, c("chrpos", "a1", "a2", "beta", "se", "maf", "p", "bonferroni")] # filter for the type

#fdrdata        = fread(paste(analysis, tissue, "step_4/qtl_by_element/qtl", paste("fdr", element_id, "txt", sep = "."), sep = "/"), data.table = F)
#qtldata$tests  = fdrdata[fdrdata$type == discovery_order,]$tests

#if (file.exists(outfile) == T)
#{
#    load(outfile, verbose = T)
#} else
#{
outlist = list()
#}

# iterate through each row in GWAS manifest
for (gwas_row in c(1:nrow(manifest)))
{
    gwas_file = manifest$filename[gwas_row] # path to GWAS sumstats
    gwas_type = manifest$trait_type[gwas_row] # trait type (continuous or case-control)
    trait_id = manifest$full_trait_id[gwas_row] # trait-id
    description = manifest$description[gwas_row] # trait-id
    
    # check if the trait has already been colocalized
    if (trait_id %in% names(outlist) & length(outlist[[trait_id]]) == 4)
    {
        message(paste("Skipping", gwas_row, trait_id, description, ". Already exists."))
    } else
    {
        message(paste(Sys.time(), gwas_row, trait_id, description))

        # get association statistics for the region of interest (1 MB from element coordinates)
        coord = paste0(pheninfo$chrom, ":", pheninfo$start - 1e6, "-", pheninfo$end + 1e6)
        cmd1 = paste("tabix", gwas_file, coord)
        cmd2 = paste("zcat", gwas_file, "| head -1")

        message(cmd1)
        message(cmd2)

        gwasdata = suppressWarnings(fread(cmd = cmd1, data.table = F, header = F))
        header   = suppressWarnings(fread(cmd = cmd2, data.table = F))

        # if there are GWAS variants in the region of interest, continue
        if (nrow(gwasdata) > 0)
        {
            colnames(gwasdata) = colnames(header)

            # check if there are GWAS variants that are above genome-wide significance, otherwise, do not continue
            if (nrow(gwasdata[gwasdata$p < 5e-05,]) > 0 & "maf" %in% colnames(header))
            {
                if (gwas_type == "case_control")
                {
                    cols = c("chrpos", "a1", "a2", "beta", "se", "p", "cases_fr", "maf", "total")
                    
                    print(colnames(gwasdata))
                    
                    missing_cols = cols[which(!cols %in% colnames(gwasdata))]
                    
                    if (length(missing_cols) > 0)
                    {
                        stop(paste("Missing", paste(cols, collase = ", ")))
                    } 
                    
                    if (!"cases_fr" %in% colnames(gwasdata)) { gwasdata$cases_fr = gwasdata$n_case / gwasdata$total }
                    
                    gwasdata = gwasdata[,c("chrpos", "a1", "a2", "beta", "se", "p", "cases_fr", "maf", "total")]
                } else
                {
                    if (!"n" %in% colnames(gwasdata))
                    {
                        gwasdata$n = gwasdata$n_case 
                        gwasdata = gwasdata[,c("chrpos", "a1", "a2", "beta", "se", "p", "n", "maf")]
                    } else
                    {
                        gwasdata = gwasdata[,c("chrpos", "a1", "a2", "beta", "se", "p", "n", "maf")]
                    }
                }

                merge = merge(qtldata, gwasdata, by = "chrpos")

                # fix opp. alleles
                tmp1 = merge[merge$a1.x == merge$a1.y & merge$a2.x == merge$a2.y,]
                tmp2 = merge[merge$a1.x == merge$a2.y & merge$a2.x == merge$a1.y,]
                tmp2$beta.y = -1 * tmp2$beta.y
                merge = rbind(tmp1, tmp2) 

                # remove multi-allelic snps
                fq = data.frame(table(merge$chrpos)) %>% filter(Freq != 1)
                merge = merge[!merge$chrpos %in% fq$Var1,]
                message(paste("Multi-allelic snps:", nrow(fq)))

                # update id
                merge = merge %>% mutate(chrpos = gsub("chr", "VAR_", paste(chrpos, a1.x, a2.x, sep = "_"))) %>% dplyr::rename(id = chrpos)

                # remove empty entries
                merge = merge[!is.na(merge$beta.x) & 
                              !is.na(merge$beta.y) & 
                              !is.na(merge$maf.x) & 
                              !is.na(merge$maf.y) & 
                              complete.cases(merge) & 
                              merge$maf.x > 0 & 
                              merge$maf.x < 1 & 
                              merge$maf.y > 0 & 
                              merge$maf.y < 1,]
                merge$maf.x  = as.double(merge$maf.x)
                merge$maf.y  = as.double(merge$maf.y)
                merge$beta.x = as.double(merge$beta.x)
                merge$beta.y = as.double(merge$beta.y)
                merge$se.x   = as.double(merge$se.x)
                merge$se.y   = as.double(merge$se.y)

                if (nrow(merge) > 50)
                {
                    dataset1 = list(type = "quant", snp = merge$id, beta = merge$beta.x, varbeta = merge$se.x^2, MAF = merge$maf.x, N = ipscore_sample_size[[paste(tissue, analysis)]])

                    if (gwas_type == "case_control")
                    {
                        dataset2 = list(type = "cc", snp = merge$id, beta = merge$beta.y, varbeta = merge$se.y^2, MAF = merge$maf.y, s = merge$cases_fr, N = merge$total)
                    } else
                    {
                        dataset2 = list(type = "quant", snp = merge$id, beta = merge$beta.y, varbeta = merge$se.y^2, MAF = merge$maf.y, N = merge$n)
                    }

                    coloc = suppressWarnings(coloc.abf(dataset1 = dataset1, dataset2 = dataset2))
                    coloc = process_coloc(coloc)
                    coloc$input = merge
                    outlist[[trait_id]] = coloc

                } else
                {
                    message("Overlapping SNPs < 50")
                }
            } else
            {
                if (!"maf" %in% colnames(header))
                {
                    message("No MAF")
                } else if (nrow(gwasdata[gwasdata$p < 5e-05,]) == 0)
                {
                    message("No GWAS variants passed significance (p < 5e-05)")
                } else
                {
                    message("Unknown error")
                }
            }
        } else
        {
            message("No GWAS variants in window")
        }
    }  
}

if (length(outlist) > 0)
{
    suppressWarnings(dir.create(outdir))
    save(outlist, file = outfile)
    message(paste("Saved:", outfile))
} else
{
    message("Empty results. No overlap with GWAS")
}



