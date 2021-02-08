# Code and data for the article Core liver homeostatic co-expression networks are preserved but respond to perturbations in an organism and disease specific manner

This repository contains R analysis code and (nearly) all data necessary to reproduce the analyses in the article *Core
liver homeostatic co-expression networks are preserved but respond to perturbations in an organism and disease specific
manner.* 

The code runs in reasonably recent versions of R (R 3.5.0 or newer should work). The code requires various R packages
available from CRAN and Bioconductor. At a minimum, the code needs the `WGCNA` package and its dependencies and packages
[`anRichment` and `anRichmentMethods`](https://horvath.genetics.ucla.edu/html/CoexpressionNetwork/GeneAnnotation/). Some of
the plotting code needs the `vioplot` package. 

Each script requires that it be run in with the directory where it resides to be set as the working directory in R (use
`setwd` to set the working directory. Each script creates, as needed, up to 3 subdirectories named `Plots` (for PDF plots), 
`Results`
(for tables with result, usually in csv format and compressed using gzip) 
and `RData` (for intermediate data that can be re-used when re-running the analysis to save time). 

Folders in the main repository are named with a numerical prefix that roughly represents the order in which the scripts or
analysis steps should be performed. 

1. `010-Preprocessing`: Preprocessing of the mouse liver RNA-seq data.
2. `020-IndividualAnalysis`: Individual gene DE analysis of the mouse liver RNA-seq data.
3. `030-NetworkAnalysis`: The first part of the analysis should be run before the network stability analysis in 
   `029-StabilityAnalysis` is carried out. The second part needs the results of the stability analysis.
4. `029-StabilityAnalysis`: Network stability analysis of the mouse liver RNA-seq data. This should be run after the first
   step in the network analysis in folder `030-NetworkAnalysis`.
5. `039-DownloadGEOData`: This is an optional step of downloading the human NAFLD (GSE126848) data from GEO. The data is
   already included in the `Data` folder.
6. `040-Preprocessing-Suppli2019-HumanNAFLD-GSE126848`: Preprocessing of the human NAFLD data.
7. `050-IndividualAnalysis-Suppli2019-HumanNAFLD-GSE126848`: Individual gene DE analysis of the human NAFLD data.
8. `060-NetworkAnalysis-Suppli2019-HumanNAFLD-GSE126848` and `059-StabilityAnalysis-Suppli2019-HumanNAFLD-GSE126848`: as
   with the mouse network analysis, the first part of the network analysis shoudl be carried out first to set up the
   network stability analysis whose results are used in the second part of the main network analysis script.
9. `070-ModulePreservation-MouseAndHuman`: Network module preservation analysis between mouse and human NAFLD data.
10. `110-Preprocessing-TCGA`: Preprocessing of TCGA data.
11. `120-IndividualAnalysis-TCGA`: Individual gene DE analysis of TCGA data using limma-voom.
12. `120-IndividualAnalysis-TCGA-V2-DESeq2`: Individual gene DE analysis of TCGA data using DESeq2.
13. `120-IndividualAnalysis-TCGA-CompareLimmaAndDESeq2`: A short script comparing limma-voom and DESeq2 results from TCGA.
14. `130-NetworkAnalysis-TCGAPhase` and `129-StabilityAnalysis-TCGA`: The first part of the network analysis shoudl
    be carried out first to set up the network stability analysis whose results are used in the second part of the main 
    network analysis script.
15. `140-ModulePreservation-TCGA`: Module preservation calculations between TCGA and the mouse and human NAFLD data.
16. `200-MODifieR`: Network analysis using MODA, ARACNE and MRNet and comparison with WGCNA.
17. `210-EnrichmentInPPINetworks`: Enrichment analysis of WGCNA modules in PPI networks.

Folder `Data` contains the raw and preprocessed data used in this study as well as annotation data needed for the
analysis. Please note that not all annotation files used by the analysis scripts are included here. For licensing reasons,
we are unable to distribute [Enrichr libraries](https://maayanlab.cloud/Enrichr/#stats) and the [Molecular Signatures
Database (MSigDB)](https://www.gsea-msigdb.org/gsea/msigdb/). Analysts who wish to re-run and/or adapt our code should
(complying with all licensing requirements) download the necessary files from these web sites, or adapt the code to avoid
using these files. 

Folder `Functions` contains supporting functions scattered over multiple files.


