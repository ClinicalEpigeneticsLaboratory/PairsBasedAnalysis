if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install(c("LOLA", "methrix", "BSgenome.Hsapiens.UCSC.hg19"))
