# Standard preprocessing

source("../Functions/networkFunctions-extras-18.R")
source("../Functions/labelPoints2-01.R");
source("../Functions/outlierRemovalFunctions.R")
source("../Functions/preprocessing-General-013.R")
source("../Functions/GNVFunctions-016.R")

library(anRichment)

dir.create("RData", recursive = TRUE);
dir.create("Results", recursive = TRUE);
dir.create("Plots", recursive = TRUE);

# Load data

dataDir = "../Data/Suppli2019-HumanNAFLD-GSE126848";
exprFiles = file.path(dataDir, "Expression/010-RawData", c("GSE126848_Gene_counts_raw.txt.gz"));

setNames0 = c("HumanNAFLDLiver");

data0 = lapply(exprFiles, function(f) 
{
  printFlush(f);
  read.table(gzfile(f), sep = "\t", quote = "", comment.char = "", header = TRUE);
})

expr0 = list();
geneAnnot0 = list();

expr0[1] = lapply(data0[1], function(.data) 
{
  out = t(.data[, -1]);
  setColnames(out, .data[, 1]);
});

ensAnnot = read.csv(gzfile("../Data/Annotation/Ensembl/EnsembleGeneAnnotation-Homo_sapiens-GRCh38.89.csv.gz"))
  
geneAnnot0 = lapply(data0, function(.data) cbind(.data[1], 
                   ensAnnot[ match(.data[, 1], ensAnnot[, 1]), -1]));
        
table(is.na(geneAnnot0[[1]]$EnsSymbol))

#FALSE  TRUE 
#19747    39 

expr0 = list2multiData(expr0);
geneAnnot0 = list2multiData(geneAnnot0);

nSets = length(expr0);

#===============================================================================================================
#
# Load sample data
#
#===============================================================================================================

sampleData = read.csv(file.path("../Data/Suppli2019-HumanNAFLD-GSE126848/SampleAnnotation/010-AsSupplied",
                                "GSE126848-sampleAnnotation-processed.csv"));

sampleAnnot1 = list();
expr1 = list();
for (set in 1:nSets)
{
  exprSamples = as.numeric(sub("X", "", rownames(expr0[[set]]$data)))
  common = intersect(exprSamples, sampleData$description);
  common2expr = match(common, exprSamples);
  expr1[[set]] = expr0[[set]]$data [ common2expr, ];
  sampleAnnot1[[set]] = cbind(sampleData[ match(common, sampleData$description), ], exprSample = rownames(expr1[[set]]));
}

expr1 = list2multiData(expr1);
sampleAnnot1 = list2multiData(sampleAnnot1);

table(sampleAnnot1[[1]]$data$disease)
# Binarize traits

diseaseLevels = c("healthy", "obese", "NAFLD", "NASH");
sampleAnnot2 = mtd.apply(sampleAnnot1, function(sa)
{
  bin = binarizeCategoricalColumns.pairwise(sa[c("gender", "disease")], convertColumns = c("gender", "disease"),
            levelOrder = list(c("Female", "Male"), diseaseLevels), includePrefix = FALSE, prefixSep = "");
  cbind(sa, bin);
});


#===============================================================================================================
#
# prettifyList
#
#===============================================================================================================

# Here the initial prettifyList is trivial

prettifyList = list(c(".vs."), c(" vs. "));
as.data.frame(prettifyList)
sampleAnnot = sampleAnnot2;

#===============================================================================================================
#
# Check gender/sex
#
#===============================================================================================================

if (FALSE)
{
gene = "DBY"
entrez1 = convert2entrez(symbol = gene, organism = "human")


  # Create a quick nearest neighbor predictor of sex from Xist

  for (set in 1:nSets)
  {
     ds1 = DESeqDataSetFromMatrix.PL(t(expr1[[set]]$data),
                                     colData = sampleAnnot[[set]]$data,
                                     design = ~1);
     ds1 = estimateSizeFactors(ds1);
     expr11 = log2(t(counts(ds1, normalized = TRUE)) + 1);
     index = match(entrez1, geneAnnot0[[set]]$data$Entrez);
     cn1 = expr11[, index]

     group = apply(sampleAnnot[[set]]$data[, c("gender", "disease")], 1, base::paste, collapse = ".");
     groupOrder = sort(unique(group));

     pdf(file = spaste("Plots/", setNames[set], "-", gene, "ExpressionVsGroup.pdf"), wi = 10, he = 4)
     #sizeGrWindow(10, 4);
     groupOrder = sort(unique(group));
     boxplot(cn1~group)

     dev.off();

     normals = tapply(cn1, sampleAnnot[[set]]$data$gender, median)
     predictedSex = names(normals)[ apply(   do.call(cbind, lapply(normals, function(x) abs(cn1-x))),
                                           1, which.min)]
     printFlush(spaste("=========== ", setNames[set], " ======================="));
     print(table(predictedSex, sampleAnnot[[set]]$data$gender));
     # perfect agreement.
     multiSampleAnnot[[set]]$data$predictedGender = predictedSex;
  }
}


#================================================================================================================
#
# Preprocessing
#
#================================================================================================================

nSamples = checkSets(expr1, checkStructure = TRUE)$nSamples
setNames = setNames0;
phenoColors = mtd.apply(mtd.subset(sampleAnnot2, , c("gender", "disease")),
                dataFrame2colors, maxOrdinalLevels = 2)
library(sva)
library(DESeq2)

groupProp = mtd.apply(sampleAnnot2, function(sa) table(sa$disease)/nrow(sa));

prepr = mtd.mapply( preprocessGeneral,
   # Main input: xdata
   xdata = expr1,
   analysisID = setNames,

   sampleAnnotation = sampleAnnot,
   phenoColors = phenoColors,
   geneAnnotation = geneAnnot0,
   stp.mainBase = setNames,

   minProportion = mtd.apply(groupProp, min),

   bw.groupsForMinWeightRestriction = mtd.subset(sampleAnnot, , "disease", drop = TRUE),

 MoreArgs = list(
   organism = "human",
   # Flow control, for some measure of interactivity
   # intermediateResults = prepr
   #stopAt = "30-VST",
   addToAnalysisManager = FALSE,
   idsAreSymbols = FALSE,
   minValue = 1,

   vst.design = "~gender + disease", 
   adj.calculateNPCs = 4,
   adj.removedCovariates = NULL,
   adj.calculateSVAFactors = TRUE,
   adj.sva.useNFactors = 2,
   adj.svaModel = "~gender + disease",
   bw.groupBy = "gender",
   bw.otherArgs = list(maxPOutliers = 0.10,
                       outlierReferenceWeight = 0.1,
                       minWeightInGroups = 0.9,
                       maxPropUnderMinWeight = 0.4,
                       defaultWeight = 1),
   # Outlier removal options
   outlierRemovalZ = 4,
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
   pca.colorColumn = "disease",
   pca.colorPrefix = "",
   pca.shapeColumn = "gender",
   pca.mainSep = ", ",
   indent = 2,
   saveDir.base = dataDir), mdmaVerbose = TRUE)

dir.create("RData", recursive = TRUE)
save(prepr, file = "RData/prepr.RData");


# Check the correspondence between original outlier Z and total sample expression

Z1 = prepr[[1]]$data$HE.VST.OR[[1]]$OSR$OSRDetails$Step.1$Z;
te = rowSums(expr1[[1]]$data);
plot(te, Z1, log ="x")


multiExpr.VST = mtd.apply(prepr, lastHE.VST.OR, "OSR", "expr.OSR");
multiCounts.all = mtd.apply(prepr, lastHE.VST.OR, "OSR", "xdata.OSR");
multiCounts = mtd.mapply(function(x, ref) x[, colnames(ref)], multiCounts.all, multiExpr.VST);

Z2 = prepr[[1]]$data$HE.VST.OR[[2]]$OSR$OSRDetails$Step.1$Z;
te2 = rowSums(multiCounts[[1]]$data);
plot(te2, Z2, log ="x")

multiWeights = mtd.apply(prepr, function(p) p$IOR$weightsForIOR);
weightFactors = mtd.apply(prepr, function(p) p$IOR$weightFactorsForIOR);

mtd.apply(multiWeights, function(x) table(x >0.1))
mtd.apply(weightFactors, function(x) table(x >1))

multiSampleAnnot = mtd.apply(prepr, lastHE.VST.OR, "OSR", "sampleAnnotation.OSR");
multiGeneAnnot = mtd.apply(prepr, lastHE.VST.OR, "HEFilter", "geneAnnotation.he");

names(multiCounts) = names(multiExpr.VST) = names(multiSampleAnnot) = names(multiGeneAnnot) = setNames;
names(multiWeights) = names(weightFactors) = setNames;

numericCols = mtd.apply(multiSampleAnnot, sapply, is.numeric);

multiNumericPheno = mtd.subset(multiSampleAnnot, , numericCols[[1]]$data);
multiNumericPheno = mtd.subset(multiNumericPheno, , 
             multiGrep(c("SVAFactor.[2-9]", "PrincipalComponent.[2-9]"), mtd.colnames(multiNumericPheno), invert = TRUE));

save(multiCounts, multiExpr.VST, multiSampleAnnot, multiGeneAnnot, multiWeights, weightFactors,
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
         bg = labels2colors(as.numeric(factor(multiSampleAnnot[[set]]$data$disease))), cex = 2,
         main = spaste(multiGeneAnnot[[set]]$data$Symbol[g], ", chr ", multiGeneAnnot[[set]]$data$EnsChromosome[g]),
         xlab = "VST expression", ylab = "Outlier statistic");
    abline(h = 1, col = "grey");
  }
}

for (set in 1:nSets)
{
  ord = order(-weightFactors[[set]]$data);
  plotIndex = ord[1:nPlots];
  plotGenes = unique(floor((plotIndex-1)/nSamples.OR[set]) + 1);
  x= as.numeric(factor(multiSampleAnnot[[set]]$data$disease));
  xl = levels(factor(multiSampleAnnot[[set]]$data$disease));
  par(mfrow = c(4,2));
  scpp(1.4);
  par(cex = 1);
  for (g in plotGenes)
  {
    plot(x, multiExpr.VST[[set]]$data[, g], 
         pch = 21,
         bg = numbers2colors(weightFactors[[set]]$data[, g], colors = grey.colors(100, start = 0, end = 1), lim = c(0,1)), 
         cex = 2,
         main = spaste(multiGeneAnnot[[set]]$data$Symbol[g], ", chr ", multiGeneAnnot[[set]]$data$EnsChromosome[g]),
         ylab = "VST expression", xlab = "disease", xaxt = "n");
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

multiNumericPheno = mtd.subset(multiSampleAnnot, , multiGrep(c("\\.vs\\.", "SVAFactor"), mtd.colnames(multiSampleAnnot)));
cp = corAndPvalue(multiNumericPheno[[1]]$data);
diag(cp$p) = NA


pdf(file = "Plots/correlationHeatmapOfTraitsAndSurrogateVariables.pdf", wi = 11, he = 7);
corHeatmapWithDendro(cp$cor, pValues = cp$p, main = "Correlations of phenotypes and surrogate variables",
      mar.main = c(6, 8, 2,1), cex.text = 0.6, dendroWidth = 0.17);

dev.off()









   





