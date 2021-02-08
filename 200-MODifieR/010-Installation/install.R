need = c("AnnotationDbi",
                       "MODA",
                       "STRINGdb",
                       "limma",
                       "org.Hs.eg.db",
                       "foreach",
                       "doParallel",
                       "Rcpp",
                       "dynamicTreeCut",
                       "reticulate",
                       "plyr",
                       "parallel",
                       "igraph",
                       "WGCNA",
                       "RSQLite",
                       "devtools",
                       "stackoverflow",
                       "preprocessCore",
                       "DESeq2",
                       "edgeR",
                       "openxlsx",
                       "ggplot2",
                       "ggdendro",
                       "ggrepel");

have = sapply(need, require, character.only = TRUE, quietly = TRUE)

BiocManager::install(pkgs = need[!have], lib = .Library, ask = FALSE)

library(devtools)

install_git(url = "https://gitlab.com/Gustafsson-lab/MODifieR.git")

