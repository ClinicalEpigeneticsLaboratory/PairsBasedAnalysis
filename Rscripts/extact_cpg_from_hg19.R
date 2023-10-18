args = commandArgs(trailingOnly=TRUE)

if (length(args)==0) {
  stop("No arguments!", call.=FALSE)
}

hg19_cpgs = methrix::extract_CPGs(ref_genome = 'BSgenome.Hsapiens.UCSC.hg19')
hg19_cpgs <- hg19_cpgs$cpgs

arrow::write_parquet(hg19_cpgs, args[1])
