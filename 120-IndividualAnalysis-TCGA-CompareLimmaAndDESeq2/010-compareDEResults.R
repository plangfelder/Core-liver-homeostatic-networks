# Standard analysis; only difference is that here we use DESeq2 to run the analysis.

source("../Functions/networkFunctions-extras-20.R")
source("../Functions/labelPoints2-01.R");
source("../Functions/outlierRemovalFunctions.R")
source("../Functions/preprocessing-General-013-02.R")
source("../Functions/GNVFunctions-016-02.R")
source("../Functions/individualAnalysis-General-008-03.R");

dir.create("Plots", recursive = TRUE);

library(anRichment)
library(DESeq2)

# Load 

files = c("associationWithPhenotypes-LICA.FR.candidateCovars-all.csv.gz",
          "associationWithPhenotypes-LICA.FR.other-all.csv.gz",
          "associationWithPhenotypes-LICA.FR.selected-all.csv.gz",
          "associationWithPhenotypes-LIHC.US.candidateCovars-all.csv.gz",
          "associationWithPhenotypes-LIHC.US.other-all.csv.gz",
          "associationWithPhenotypes-LIHC.US.selected-all.csv.gz",
          "associationWithPhenotypes-LIRI.JP.candidateCovars-all.csv.gz",
          "associationWithPhenotypes-LIRI.JP.other-all.csv.gz",
          "associationWithPhenotypes-LIRI.JP.selected-all.csv.gz");

res.limma = lapply(files, function(f) 
  read.csv(gzfile(file.path("../120-IndividualAnalysis-TCGA/Results/01-Association", f)), check.names = FALSE));

res.DESeq = lapply(files, function(f)
  read.csv(gzfile(file.path("../120-IndividualAnalysis-TCGA-V2-DESeq2/Results/01-Association", f)), check.names = FALSE));

nFiles = length(files);
correlations = lapply(1:nFiles, function(f)
{
  res1 = res.limma[[f]];
  res2 = res.DESeq[[f]];
  ZCols1 = grepv("Z for ", names(res1));
  ZCols2 = grepv("Z for ", names(res2));
  #if (!isTRUE(all.equal(make.names(ZCols1), make.names(ZCols2)))) browser();
  setNames(diag(cor(res1[, ZCols1], res2[, ZCols2], use = 'p')), ZCols2);
});

cor.flat = unlist(correlations)

pdf(file = "Plots/histogramOfDESeq2-LimmaCorrelations.pdf", wi = 6, he = 4);
scpp(1.5);
hist(cor.flat, breaks = 20, main = "Correlations between Z statistics from limma-voom and DESeq2",
   xlab = "Correlation", cex.main = 1);

dev.off();
  

f = 5

res1 = res.limma[[f]];
  res2 = res.DESeq[[f]];
  ZCols1 = grepv("Z for ", names(res1));
  ZCols2 = grepv("Z for ", names(res2));

thinnedScatterplot(res1[[ZCols1[[1]] ]], res2[[ZCols2[[1]] ]], verbose = TRUE)
