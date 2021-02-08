README for code for project eqtl-tcga
=====================================

Daniel Dvorkin <daniel@bioinformaticscro.com> February 8, 2021

This document describes the R code for the `eqtl-tcga` project.  The following sections describe the directory structure, required packages, and execution of the code.  Please contact the author with any questions.


Directory structure
===================

The expected structure of the project directory is as follows.  Note that the "WGCNA-TCGA-liver" and "tcga_survival_files" directories contain protected TCGA data and are not included here.

├── code
├── processed_data
├── raw_data
│   ├── WGCNA-TCGA-liver
│   │   ├── 110-Preprocessing-TCGA
│   │   │   ├── Plots
│   │   │   └── Results
│   │   ├── 120-IndividualAnalysis-TCGA
│   │   │   ├── Plots
│   │   │   └── Results
│   │   │       ├── 01-Association
│   │   │       └── 02-Enrichment
│   │   ├── 129-StabilityAnalysis-TCGA
│   │   ├── 130-NetworkAnalysis-TCGAPhase
│   │   │   ├── Plots
│   │   │   │   ├── LICA-FR.main
│   │   │   │   ├── LICA-FR.sub
│   │   │   │   ├── LIHC-US.main
│   │   │   │   ├── LIHC-US.sub
│   │   │   │   ├── LIRI-JP.main
│   │   │   │   └── LIRI-JP.sub
│   │   │   └── Results
│   │   │       └── EnrichmentAnalysis
│   │   │           ├── LICA-FR.main
│   │   │           ├── LICA-FR.sub
│   │   │           ├── LIHC-US.main
│   │   │           ├── LIHC-US.sub
│   │   │           ├── LIRI-JP.main
│   │   │           └── LIRI-JP.sub
│   │   ├── 140-ModulePreservation-TCGA
│   │   │   ├── Plots
│   │   │   └── Results
│   │   ├── Data
│   │   │   ├── Annotation
│   │   │   │   ├── Enrichr
│   │   │   │   ├── Ensembl
│   │   │   │   ├── MSigDB
│   │   │   │   └── NCBI
│   │   │   ├── TCGA
│   │   │   │   ├── Expression
│   │   │   │   │   ├── 010-RawData
│   │   │   │   │   ├── 011-RawData-Reformatted
│   │   │   │   │   ├── 015-XDataWithValidEntrez
│   │   │   │   │   └── 030-OutlierSamplesRemoved
│   │   │   │   └── SampleAnnotation
│   │   │   │       ├── 010-AsSupplied
│   │   │   │       ├── 020-AsSupplied
│   │   │   │       └── 030-OutlierSamplesRemoved
│   │   │   └── anRichmentData
│   │   ├── Functions
│   │   └── Report
│   └── tcga_survival_files
└── reports
    ├── main_KM_plots
    └── sub_KM_plots


Required packages
=================

Execution of the code requires the following CRAN packages to be installed:

- survival 
- survminer 
- matrixStats 
- stringr 
- readxl 
- reshape2 
- data.table 
- knitr
- seqminer
- MatrixEQTL

and the Bioconductor library `biomaRt`.

It also requires two packages written by the author, `ARC.utils` and `lcmix`.  These packages are included in the project directory and can be installed in the usual way from the command line, i.e. `R CMD INSTALL ARC.utils_0.1-2.tar.gz` and `R CMD INSTALL lcmix_0.3.tar.gz`.


Execution
=========

The analysis can be executed from the R console, inside the "code" directory, with the usual command `source("analyze_data.R")`.  Note that the code is separated into blocks of the form `if(TRUE) { ... }` or `if(1) { ... }` to allow execution of one block at a time.  Change `FALSE` to `TRUE` or `0` to `1` in these outer `if` statements to execute.
