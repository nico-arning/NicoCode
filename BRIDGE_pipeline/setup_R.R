options(repos = c(CRAN = "https://cran.r-project.org"))

install.packages(c("BiocManager", 
                 "gplots",
                 "argparse",
                 "pheatmap",
                 "pcaPP"))
#install.packages("pcaPP")

library(BiocManager)

packages = c(
    "DESeq2",
    "ggalt",
    "EnhancedVolcano",
    "S4Vectors",
    "BiocGenerics",
    "IRanges",
    "GenomicRanges",
    "GenomeInfoDb",
    "SummarizedExperiment",
    "MatrixGenerics",
    "matrixStats",
    "Biobase",
    "ggrepel",
    "limma",
    "sva", 
    "NOISeq",
    "edgeR",
    #"bigPint",
    "latticeExtra",
    "enrichplot",
    "AnnotationHub",
    "clusterProfiler",
    "pathview",
    "GOplot"
    #"FactoMineR"
)

BiocManager::install(packages)
#, lib = "/mnt/datapond/Dependencies/R_libs", dependencies=TRUE)
BiocManager::install("ggpubr")

BiocManager::install("S4Vectors")
devtools::install_github("valentint/rrcov")
