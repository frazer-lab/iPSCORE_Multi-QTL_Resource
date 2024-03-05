.libPaths(c(.libPaths(), "/frazer01/home/tarthur/software/R-4.1.0/lib64/R/library"))

suppressMessages(library(optparse))
suppressMessages(library(ATACseqQC))
suppressMessages(library(TxDb.Hsapiens.UCSC.hg38.knownGene))
suppressMessages(library(Rsamtools))

option_list = list(make_option("--pipe_dir", type = "character", default = 0, help = "pipe dir", metavar = "character"))

opt_parser = OptionParser(option_list = option_list)
opt = parse_args(opt_parser)
pipe_dir = opt$pipe_dir

message(Sys.time())
message(paste("Out directory:", pipe_dir))

#seqlev <- "chr1" ## subsample data for quick run
seqlev <- "chr21"
seqinformation <- seqinfo(TxDb.Hsapiens.UCSC.hg38.knownGene)
which <- as(seqinformation[seqlev], "GRanges")

## bamfile tags to be read in
possibleTag <- list("integer"=c("AM", "AS", "CM", "CP", "FI", "H0", "H1", "H2", 
                                "HI", "IH", "MQ", "NH", "NM", "OP", "PQ", "SM",
                                "TC", "UQ"), 
                 "character"=c("BC", "BQ", "BZ", "CB", "CC", "CO", "CQ", "CR",
                               "CS", "CT", "CY", "E2", "FS", "LB", "MC", "MD",
                               "MI", "OA", "OC", "OQ", "OX", "PG", "PT", "PU",
                               "Q2", "QT", "QX", "R2", "RG", "RX", "SA", "TS",
                               "U2"))

calc_tsse = function(pipe_dir, tags, which, possibleTag)
{
    bamfile = paste(pipe_dir, "Aligned.sorted.filt.nodup.nomito.bam", sep = "/")
    bamTop100 <- scanBam(BamFile(bamfile, yieldSize = 100),param = ScanBamParam(tag=unlist(possibleTag)))[[1]]$tag
    tags <- names(bamTop100)[lengths(bamTop100)>0]

    gal <- readBamFile(bamfile, tag=tags, which=which, asMates=TRUE, bigFile=TRUE)

    shiftedBamfile <- file.path(paste(pipe_dir, "Aligned.sorted.filt.nodup.shifted.bam", sep = "/"))
    gal1 <- shiftGAlignmentsList(gal)

    txs <- transcripts(TxDb.Hsapiens.UCSC.hg38.knownGene)

    tsse <- TSSEscore(gal1, txs)
    save(tsse, file = paste(pipe_dir, "qc", "TSSEscore.robj", sep = "/"))
    message(paste("Saved:", paste(pipe_dir, "qc", "TSSEscore.robj", sep = "/")), appendLF = F)
    return(tsse) 
}

calc_tsse(pipe_dir, tags, which, possibleTag) 
