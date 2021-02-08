# Standard preprocessing

source("../Functions/networkFunctions-extras-17.R")
source("../Functions/labelPoints2-01.R");
source("../Functions/outlierRemovalFunctions.R")
source("../Functions/preprocessing-General-013.R")
source("../Functions/GNVFunctions-015.R")

library(anRichment)

dir.create("RData", recursive = TRUE);
dir.create("Results", recursive = TRUE);
dir.create("Plots", recursive = TRUE);

# Load data

exprFiles = file.path("../Data/Expression/010-RawData",
    c("Diet 2m Raw Read Counts(for WGCNA).csv.bz2"));

setNames0 = c("Liver");

data0 = lapply(exprFiles, function(f) 
{
  printFlush(f);
  read.csv(bzfile(f));
})

expr0 = list();
geneAnnot0 = list();

expr0[1] = lapply(data0[1], function(.data) 
{
  out = t(.data[, -1]);
  setColnames(out, .data$Gene);
});
  
geneAnnot0[1] = lapply(data0[1], function(.data) cbind(.data[1], 
        geneAnnotationFromEntrez(convert2entrez(symbol = .data[, 1], organism = "mouse"), organism = "mouse",
                                 includePosition = TRUE)));

table(is.na(geneAnnot0[[1]]$Entrez))
#FALSE  TRUE 
#28974  5462 

expr0 = list2multiData(expr0);
geneAnnot0 = list2multiData(geneAnnot0);

nSets = length(expr0);

#===============================================================================================================
#
# Load sample data
#
#===============================================================================================================

sampleData = read.csv("../Data/SampleAnnotation/010-AsSupplied/Phenotypic data.csv");

sampleAnnot1 = list();
expr1 = list();
for (set in 1:nSets)
{
  exprSamples = rownames(expr0[[set]]$data);
  common = intersect(exprSamples, sampleData$Samples);
  common2expr = match(common, exprSamples);
  expr1[[set]] = expr0[[set]]$data [ common2expr, ];
  sampleAnnot1[[set]] = cbind(sampleData[ match(common, sampleData$Samples), ], exprSample = rownames(expr1[[set]]));
  sampleAnnot1[[set]]$Diet.machine = gsub("2\\%", "2percent", sampleAnnot1[[set]]$Diets);
  sampleAnnot1[[set]]$Diet.machine = multiGSub(c("\\%",  "\\.+", "_+", "_$"), c(".perc", "_", "_", ""), 
                                         make.names(sampleAnnot1[[set]]$Diet.machine));
}

expr1 = list2multiData(expr1);
sampleAnnot1 = list2multiData(sampleAnnot1);

# Binarize diets

dietLevelTranslation.paper = cbind(c("NC", "HS", "Cholate (CA) 2m", "ChR 2% (2m)", "Chol 2%_HS (2m)", "Chol2%-CA"),
            c("Normal", "HS", "Cholate (CA)", "HS-Chol2%-CA", "HS-Chol2%", "Chol2%-CA"));

dietOrder.paper = c("Normal", "HS", "HS-Chol2%", "HS-Chol2%-CA", "Chol2%-CA", "Cholate (CA)");
dietOrder = translateUsingTable(dietOrder.paper, dietLevelTranslation.paper[, c(2,1)]);
dietLevels = dietOrder;
dietLevels.machine = translateUsingTable(dietLevels, mtd.rbindSelf(mtd.subset(sampleAnnot1, , c("Diets", "Diet.machine"))))
sampleAnnot2 = mtd.apply(sampleAnnot1, function(sa)
{
  bin = binarizeCategoricalColumns.pairwise(sa["Diet.machine"], convertColumns = "Diet.machine",
            levelOrder = list(dietLevels.machine), includePrefix = FALSE, prefixSep = "");
  cbind(sa, bin);
});


#===============================================================================================================
#
# prettifyList
#
#===============================================================================================================


prettifyList = list(
  from = c(dietLevels.machine, mtd.colnames(sampleAnnot1)),
  to = c(dietLevelTranslation.paper[, 2], gsub("\\.", " ", mtd.colnames(sampleAnnot1))));

as.data.frame(prettifyList)

sampleAnnot = mtd.apply(sampleAnnot2, function(sa) 
{
   sa$Diet.paper = translateUsingTable(sa$Diets, dietLevelTranslation.paper);
   sa$Diet.paper.factor = factor(sa$Diet.paper, levels = dietOrder.paper);
   sa;
});
#================================================================================================================
#
# Preprocessing
#
#================================================================================================================

nSamples = checkSets(expr1, checkStructure = TRUE)$nSamples
setNames = setNames0;
phenoColors = mtd.apply(mtd.subset(sampleAnnot2, , c("Diets", "Weight", "Liver.weight", "Liver.body.weight.ratio",
                   "Spleen.weight", "Epi.Fat", "Cecum.weight", "Liver.total.cholesterol", "Liver.free.cholesterol",
                   "Lobular.inflammation", "Steatosis", "balooning", "Portal.inflammation")),
                dataFrame2colors, maxOrdinalLevels = 2)
library(sva)
library(DESeq2)

groupProp = mtd.apply(sampleAnnot2, function(sa) table(sa$Diets)/nrow(sa));

prepr = mtd.mapply( preprocessGeneral,
   # Main input: xdata
   xdata = expr1,
   analysisID = setNames,

   sampleAnnotation = sampleAnnot,
   phenoColors = phenoColors,
   geneAnnotation = geneAnnot0,
   stp.mainBase = setNames,

   minProportion = mtd.apply(groupProp, min),

   bw.groupsForMinWeightRestriction = mtd.subset(sampleAnnot, , "Diets", drop = TRUE),

 MoreArgs = list(
   organism = "mouse",
   # Flow control, for some measure of interactivity
   # intermediateResults = prepr
   #stopAt = "30-VST",
   addToAnalysisManager = FALSE,
   idsAreSymbols = FALSE,
   minValue = 1,

   vst.design = "~Diets", 
   adj.calculateNPCs = 4,
   adj.removedCovariates = NULL,
   adj.calculateSVAFactors = TRUE,
   adj.sva.useNFactors = 1,
   adj.svaModel = "~Diets",
   bw.groupBy = NULL,
   bw.otherArgs = list(maxPOutliers = 0.10,
                       outlierReferenceWeight = 0.1,
                       minWeightInGroups = 0.9,
                       maxPropUnderMinWeight = 0.4,
                       defaultWeight = 1),
   # Outlier removal options
   outlierRemovalZ = 6,
   iorFromHEGenes = TRUE,
   ior.replace = FALSE,
   ior.remove = FALSE,

   restrictIORExprToHE = TRUE,
   restrictIORXDataToHE = FALSE,
   stp.width = 8,
   stp.height = 5,
   stp.mainSep = ", ",
   stp.marAll = c(1, 15, 3, 1),
   stp.maxCharPerLine = 47,
   plotDir = "Plots",

   pca.correctForColumn = NULL,
   pca.colorColumn = "Diet.paper.factor",
   pca.colorPrefix = "",
   pca.mainSep = ", ",
   indent = 2,
   saveDir.base = "../Data"), mdmaVerbose = TRUE)

dir.create("RData", recursive = TRUE)
save(prepr, file = "RData/prepr.RData");



multiExpr.VST = mtd.apply(prepr, lastHE.VST.OR, "OSR", "expr.OSR");
multiExpr.all = mtd.apply(prepr, lastHE.VST.OR, "OSR", "xdata.OSR");
multiExpr = mtd.mapply(function(x, ref) x[, colnames(ref)], multiExpr.all, multiExpr.VST);

multiWeights = mtd.apply(prepr, function(p) p$IOR$weightsForIOR);
weightFactors = mtd.apply(prepr, function(p) p$IOR$weightFactorsForIOR);

mtd.apply(multiWeights, function(x) table(x >0.1))
mtd.apply(weightFactors, function(x) table(x >1))

multiSampleAnnot = mtd.apply(prepr, lastHE.VST.OR, "OSR", "sampleAnnotation.OSR");
multiGeneAnnot = mtd.apply(prepr, lastHE.VST.OR, "HEFilter", "geneAnnotation.he");

names(multiExpr) = names(multiExpr.VST) = names(multiSampleAnnot) = names(multiGeneAnnot) = setNames;
names(multiWeights) = names(weightFactors) = setNames;

numericCols = mtd.apply(multiSampleAnnot, sapply, is.numeric);

multiNumericPheno = mtd.subset(multiSampleAnnot, , numericCols[[1]]$data);
multiNumericPheno = mtd.subset(multiNumericPheno, , 
             multiGrep(c("SVAFactor.[2-9]", "PrincipalComponent.[2-9]"), mtd.colnames(multiNumericPheno), invert = TRUE));

save(multiExpr, multiExpr.VST, multiSampleAnnot, multiGeneAnnot, multiWeights, weightFactors,
     prettifyList, multiNumericPheno,
     setNames, file = "RData/preprocessedData.RData");

sink(file = "Results/sessionInfo.txt");
sessionInfo();
sink(NULL);




#=============================================================================================================
#
# Sanity check for weights/weight factors
#
#=============================================================================================================

nSamples.OR = checkSets(multiExpr, checkStructure = TRUE)$nSamples;
nPlots = 24;

pdf(file = "Plots/GenesWithOutlierWeights.pdf", wi = 15, he = 10)
for (set in 1:nSets)
{
  ord = order(-weightFactors[[set]]$data);
  plotIndex = ord[1:nPlots];
  plotGenes = unique(floor((plotIndex-1)/nSamples.OR[set]) + 1);
  par(mfrow = c(4,6));
  scpp(1.4);
  par(cex = 1);
  for (g in plotGenes)
  {
    plot(multiExpr.VST[[set]]$data[, g], weightFactors[[set]]$data[, g],
         pch = 21,
         bg = labels2colors(as.numeric(factor(multiSampleAnnot[[set]]$data$Diets))), cex = 2,
         main = spaste(multiGeneAnnot[[set]]$data$Gene[g], ", chr ", multiGeneAnnot[[set]]$data$Chr[g]),
         xlab = "VST expression", ylab = "Outlier statistic");
    abline(h = 1, col = "grey");
  }
}

for (set in 1:nSets)
{
  ord = order(-weightFactors[[set]]$data);
  plotIndex = ord[1:nPlots];
  plotGenes = unique(floor((plotIndex-1)/nSamples.OR[set]) + 1);
  x= as.numeric(factor(multiSampleAnnot[[set]]$data$Diets));
  xl = levels(factor(multiSampleAnnot[[set]]$data$Diets));
  par(mfrow = c(4,2));
  scpp(1.4);
  par(cex = 1);
  for (g in plotGenes)
  {
    plot(x, multiExpr.VST[[set]]$data[, g], 
         pch = 21,
         bg = numbers2colors(weightFactors[[set]]$data[, g], colors = grey.colors(100, start = 0, end = 1), lim = c(0,1)), 
         cex = 2,
         main = spaste(multiGeneAnnot[[set]]$data$Gene[g], ", chr ", multiGeneAnnot[[set]]$data$Chr[g]),
         ylab = "VST expression", xlab = "Age", xaxt = "n");
    axis(1, at = sort(unique(x)), labels = xl, cex.axis = 0.5);
    abline(h = 1, col = "grey");
  }
}


dev.off()

#=============================================================================================================
#
# Check if the removed 1st PC/1st SVA factor correlates with any of the traits.
#
#=============================================================================================================
cp = corAndPvalue(multiNumericPheno[[1]]$data);
diag(cp$p) = NA


pdf(file = "Plots/correlationHeatmapOfTraitsAndSurrogateVariables.pdf", wi = 17, he = 10);
corHeatmapWithDendro(cp$cor, pValues = cp$p, main = "Correlations of phenotypes and surrogate variables",
      mar.main = c(14, 18, 2,1), cex.text = 0.6, dendroWidth = 0.17);

dev.off()







   





