# splice analysis
# work directory: /data/seed/Chordoma/Pita/SpliceAnalysis

###
# this scripts will take a list of entrez gene id as input and print their corresponding splice
# variation expression levels in Chordoma/Pita/mapping data. if you specify "verbose" in command
# line arguments, some checkpoints will be activated
# it can only use R/3.3.0, due to version of some packages in library
# For example: `Rscript get_splice_result.R 9821 9823 verbose`

library("rlang", lib.loc="/PHShome/mj137/R/x86_64-pc-linux-gnu-library/3.3/")
library("SGSeq")
library("TxDb.Hsapiens.UCSC.hg19.knownGene")
library("org.Hs.eg.db")
setwd("/data/seed/Chordoma/Pita/SpliceAnalysis")

#---------------------------------------------------------------------------------------------------------

# interpret command line arguemnt into two parts: verbose and a list of entrez ids
# output an R list containing all this information
argument_interpreter <- function() {
  args <- commandArgs(trailingOnly = TRUE)
  verbose <- FALSE
  name_of_args <- grepl("-", args)
  
  if (length(args) == 1) {
    stop("need to type at least one entrez ID")
  }
  if ("-verbose" %in% args[name_of_args]) {
    verbose <- TRUE
  }
  
  # load entrez ID
  if ('-entrez_id' %in% args){
    entrez_id <- args[!name_of_args]
  } else if ('-gene_name' %in% args) {
    entrez_id <- geneName2Entrez(args[!name_of_args])
    print(entrez_id)
  } else if ('-input_file' %in% args) {
    entrez_id <- geneName2Entrez(readLines(args[!name_of_args]))
  }
  
  return(list(entrez_id=entrez_id, verbose=verbose))
}

#---------------------------------------------------------------------------------------------------------

geneName2Entrez <- function(gene_name) {
  symbols <- toupper(gene_name)
  converter <- mapIds(org.Hs.eg.db, symbols, 'ENTREZID', 'SYMBOL')
  return(as.character(converter))
}

#---------------------------------------------------------------------------------------------------------

# the sample file include name of sample, path to their bam files, whether it's paired end or not,
# the average length of reads, the average length of fragments, total library size
# the average length of reads, the average length of fragments can be fetched using following command:
# `awk 'function abs(v) {return v <0 ? -v : v} !/^@/ {len_count+=length($10);count+=1;frag_length+=abs($9)} END {print len_count/count;print frag_length/count}' Aligned.out.sam`
get_sample_file <- function() {
  path <- vector()
  for (i in 1:length(list.files("../Mapping/"))) {
    path[i] <- paste("../Mapping/",list.files("../Mapping/")[i], "/Aligned.sorted.bam", sep="")
  }
  
  sample_data <- data.frame(sample_name = list.files("../Mapping/"),
                            file_bam = path,
                            paired_end = c(TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE),
                            read_length = c(50.8559, 50.8629, 50.8534, 50.8538, 50.8629, 50.8653, 50.8096, 50.8256, 50.8264),
                            frag_length = c(1777.92, 1829.46, 1885.85, 1658.52, 1601.32, 1614.83, 2306.78, 2010.58, 1860.93),
                            lib_size = c(116910295, 345523260, 11576825, 112959626, 142482562, 125112631, 66309442, 96109950, 85481880),
                            stringsAsFactors = FALSE)
  
  return(sample_data)
}
#---------------------------------------------------------------------------------------------------------

# this function will generate a gene_range object of certain gene
get_gene_IRange <- function(entrezID) {
  anno_db <- select(org.Hs.eg.db, keys=as.character(entrezID),
                    keytype="ENTREZID",
                    columns=c("SYMBOL","CHRLOC","CHRLOCEND"))
  
  if (anno_db$CHRLOC > 0) {
    anno_sign <- "+"
    anno_start <- anno_db$CHRLOC
    anno_end <- anno_db$CHRLOCEND
  } else {
    anno_sign <- "-"
    anno_start <- -1 * anno_db$CHRLOC
    anno_end <- -1 * anno_db$CHRLOCEND
  }
  gr_obj <- GRanges(seqnames=anno_db$CHRLOCCHR,
                    ranges = anno_start:anno_end,
                    strand = anno_sign)
  return(gr_obj)
}

#---------------------------------------------------------------------------------------------------------

# get transcript information from TxDb.Hsapiens.UCSC.hg19.knownGene database
# onle select the region where the certain gene is located
get_trancript_features <- function(entrez_id, verbose) {
  gr <- get_gene_IRange(entrez_id)
  chr_num <- paste("chr",as.character(seqnames(gr)),sep="")
  
  txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene
  restoreSeqlevels(txdb)
  txdb_keep <- keepSeqlevels(txdb, unique(chr_num))
  seqlevelsStyle(txdb) <- "NCBI"
  txf_ucsc <- convertToTxFeatures(txdb)
  txf_ucsc <- txf_ucsc[txf_ucsc %over% gr]
  
  if (verbose) {
    cat("\n\n+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n         print transcript features in gene region\n++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n")
    print(head(txf_ucsc))
    print(type(txf_ucsc))
    print(head(txName(txf_ucsc)))
    print(head(geneName(txf_ucsc)))
  }
  return(txf_ucsc)
}

#---------------------------------------------------------------------------------------------------------

# count the expression level of each part of transcripts and return a transcript object
analyze_transcript <- function(transcript_features, sample_data, verbose) {
  sgf_ucsc <- convertToSGFeatures(transcript_features)
  head(sgf_ucsc)
  
  sgfc_ucsc <- analyzeFeatures(sample_data, features = transcript_features)
  if (verbose) {
    cat("\n\n+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n         print the number of alternative splices\n++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n")
    print(colData(sgfc_ucsc))
    print(rowRanges(sgfc_ucsc))
    print(head(counts(sgfc_ucsc)))
    print(head(FPKM(sgfc_ucsc)))
  }
  return(sgfc_ucsc)
}

#---------------------------------------------------------------------------------------------------------

main <- function() {
  arguments <- argument_interpreter()
  sample_data <- get_sample_file()
  # loop through each entrez ID
  for (i in 1:length(arguments[['entrez_id']])) {
    entrez_id <- arguments[['entrez_id']][i]
    transcript_features <- get_trancript_features(entrez_id, verbose=arguments[["verbose"]])
    print("===========")
    transcripts_expressions <- analyze_transcript(transcript_features, sample_data, verbose=arguments[["verbose"]])
    # plot the expression levels of splice alternative in a pdf file/
    pdf(paste(entrez_id, ".pdf", sep = ""), width = 14, height = 10)
    df <- plotFeatures(transcripts_expressions, geneID = 1, include = "both")
    dev.off()
    cat("\n\n===================================================\n         a new pdf has created\n===================================================\n")
  }
}

options(warn=-1)
main()