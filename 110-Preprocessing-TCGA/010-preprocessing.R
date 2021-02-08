# Standard preprocessing

source("../Functions/networkFunctions-extras-19.R")
source("../Functions/labelPoints2-01.R");
source("../Functions/outlierRemovalFunctions.R")
source("../Functions/preprocessing-General-013.R")
source("../Functions/GNVFunctions-016.R")

library(anRichment)

dir.create("RData", recursive = TRUE);
dir.create("Results", recursive = TRUE);
dir.create("Plots", recursive = TRUE);

source("../../../../RLibs/inHouseGeneAnnotation/anRichmentMethods/R/enrichmentAnalysis.R");

# Load data

baseDir = "../Data/TCGA/";
setNames = c("LICA-FR", "LIHC-US", "LIRI-JP");
exprDir.11 = file.path(baseDir, "Expression/011-RawData-Reformatted");
nSets = length(setNames)
counts0 = list();
for (set in 1:nSets)
{
  counts0[[set]] = loadTable(file = gzfile(spaste(exprDir.11, "/countMatrix-", setNames[set], ".csv.gz")),
        transpose = TRUE, convertToMatrix = TRUE, sep = ",", header = TRUE);
}

names(counts0) = setNames;

sampleDir.11 = file.path(baseDir, "SampleAnnotation/011-Combined");
sampleAnnot1 = lapply(setNames, function(.name)
  read.csv( file = gzfile(spaste(exprDir.11, "/sampleData-ICGCcombinedAnnotation-", .name, ".csv.gz"))));

names(sampleAnnot1) = setNames;

clinicalData1 = read.csv(file.path("../Data/TCGA/SampleAnnotation/010-AsSupplied",
                   "List_of_HCC_samples_in_TCGA_and_ICGC-ICGCLiverAndNormalTissue-Bioportal.csv"),
                   check.names = FALSE)

sampleData2 = read.csv(file.path("../Data/TCGA/SampleAnnotation/010-AsSupplied",
                   "List_of_HCC_samples_in_TCGA_and_ICGC-ICGCLiverAndNormalTissue.csv"),
                   check.names = FALSE)


ensemblAnnot = read.csv(gzfile("../Data/Annotation/Ensembl/Homo_sapiens_GRCh37p7.txt.gz"));
NCBIAnnot = read.csv(gzfile("../Data/Annotation/NCBI/geneAnnotation-NCBIHuman37.3.csv.gz"));



lapply(sampleAnnot1, function(sa) table(sa$icgc_specimen_id %in% sampleData2$icgc_specimen_id))

sampleID1 = lapply(sampleAnnot1, getElement, "icgc_specimen_id");
sampleID2 = lapply(sampleAnnot1, getElement, "submitted_specimen_id");
sampleID3 = sampleID2;
sampleID3[[2]] = sub("[A-Z]$", "", sampleID3[[2]]);

rows.clinicalData1 = lapply(sampleID3, match, clinicalData1$`Sample ID`)

lapply(rows.clinicalData1, function(x) table(is.na(x)))

clinicalData1.IDs = clinicalData1[ , grep("ID", names(clinicalData1))];

#============================================================================================================
#
# Fish out clinical data. 
#
#============================================================================================================
# We're missing a lot of the samples, but at least we have some.

sampleAnnot2 = mymapply(function(sa, rows) 
     data.frame(sa, clinicalData1[rows, !names(clinicalData1) %in% names(sa)], check.names = FALSE),
     sampleAnnot1, rows.clinicalData1);   ### This one contains lots of missing samples

ensemblAnnot2 = setNames(ensemblAnnot, multiSub(c("EntrezGene.ID", "Gene.Start..bp.", "Chromosome.Name", "HGNC.symbol"),
                                                 c("Entrez", "Loc", "Chr", "Symbol"), names(ensemblAnnot)));

geneAnnot0 = list();

geneAnnot0[[1]] = data.frame(GeneID = colnames(counts0[[1]]),
                      ensemblAnnot2[ match(colnames(counts0[[1]]), ensemblAnnot2$Ensembl.Gene.ID), ])
table(is.na(geneAnnot0[[1]]$Ensembl.Gene.ID))
# FALSE  TRUE 
# 46812  6268 



for (set in 2:3)
{
  geneAnnot0[[set]] = data.frame(GeneID = colnames(counts0[[set]]),
                       geneAnnotationFromEntrez(convert2entrez(symbol = colnames(counts0[[set]]), organism = "human"), 
                                                organism = "human", includePosition = TRUE));
}

sapply(geneAnnot0, function(ga) sum(!is.na(ga$Entrez)))
#[1] 21817 19903 21721


#===============================================================================================================
#
# prettifyList
#
#===============================================================================================================

prettifyList = list(from = "A", to = "A");

as.data.frame(prettifyList)

sampleAnnot3 = list2multiData(sampleAnnot2);
names(counts0) = names(geneAnnot0) = names(sampleAnnot2) = setNames;
maxInt = 2^31-1;
counts0 = lapply(counts0, function(x) {x[x > maxInt] = maxInt; x});
counts1 = list2multiData(counts0);

#===============================================================================================================
#
# predict sex from expression data
#
#===============================================================================================================

library(DESeq2)
genes = c("XIST")
for (set in 1:nSets)
{
     sa = sampleAnnot3[[set]]$data;
     rownames(sa) = rownames(counts1[[set]]$data);
     ds1 = DESeqDataSetFromMatrix.PL(t(counts1[[set]]$data), colData = data.frame(a = rnorm(nrow(sa))), design = ~1);
     ds1 = estimateSizeFactors(ds1);
     expr1 = log2(t(counts(ds1, normalized = TRUE)) + 1);
     for (gene in genes)
     {
        index = match(gene, geneAnnot0[[set]]$Symbol);
        cn1 = expr1[, index]

        pdf(file = spaste("Plots/expressionVsSex-", gene, "-", setNames[set], ".pdf"), wi = 5, he = 5);
        scpp(1.3);
        labeledBoxplot(tapply(cn1, sampleAnnot3[[set]]$data$donor_sex, identity), 
                names = sort(unique(sampleAnnot3[[set]]$data$donor_sex)),
                main = setNames[set], ylab = spaste(gene, " expression"),
                xlab = "Sex", addScatterplot = TRUE, notch = TRUE);

        dev.off();
        normals = tapply(cn1, sampleAnnot3[[set]]$data$donor_sex, median, na.rm = TRUE)
        predictedSex = names(normals)[ apply(   do.call(cbind, lapply(normals, function(x) abs(cn1-x))),
                                              1, which.min)]
        printFlush(spaste("=========== ", setNames[set], " ======================="));
        print(table(predictedSex, sampleAnnot3[[set]]$data$donor_sex));
        sampleAnnot3[[set]]$data$predictedSex = predictedSex;
     }
}


#================================================================================================================
#
# Preprocessing
#
#================================================================================================================

nSamples = checkSets(counts1, checkStructure = TRUE)$nSamples


traitsOfInterest = c("Disease Stage", "Liver fibrosis ishak score category", "Tumor Size", "Vascular Invasion",
          "Subtype");

nWithData = do.call(cbind, mtd.apply(sampleAnnot3, function(x) sapply(traitsOfInterest, function(tr)
{
  x1 = x[[tr]];
  printFlush("================================================================")
  printFlush(tr);
  print(table(!is.na(x1)));
  print(table(x1));
  sum(!is.na(x1));
}), returnList = TRUE))

colnames(nWithData) = setNames;

nWithData = rbind(nSamples = t(as.matrix(c(mtd.apply(sampleAnnot3, nrow, mdaSimplify = TRUE)))), nWithData)
rownames(nWithData)[1] = "Total samples";

write.csv.nr(dataForTable(nWithData, transpose = FALSE, IDcolName = "Trait"), "Results/numbersOfSamplesWithClinicalData.csv");



traits1 = c("donor_sex", "donor_age_at_diagnosis", "specimen_type",
               "specimen_donor_treatment_type", 
               "tumour_histological_type", traitsOfInterest);

mtd.apply(mtd.subset(sampleAnnot3, , traits1), lapply, table);

# For pheno colors: drop uniformative columns

sampleAnnot4 = mtd.apply(mtd.subset(sampleAnnot3, , traits1), dropConstantColumns);

phenoColors = mtd.apply(sampleAnnot4, dataFrame2colors, maxOrdinalLevels = 5, maxLegendLength = 120)
library(sva)
library(DESeq2)

groupProp = rep(1/4, nSets);

prepr = list();
for (set in 1:nSets)
prepr[[set]] = list(data = preprocessGeneral(
   # Main input: xdata
   xdata = counts1[[set]]$data,
   analysisID = setNames[set],
   sampleAnnotation = sampleAnnot3[[set]]$data,
   phenoColors = phenoColors[[set]]$data,
   geneAnnotation = geneAnnot0[[set]],
   stp.mainBase = setNames[set],
   minProportion = groupProp[set],
   bw.groupsForMinWeightRestriction = mtd.subset(sampleAnnot3, , "donor_sex", drop = TRUE)[[set]]$data,

   organism = "human",
   # Flow control, for some measure of interactivity
   # intermediateResults = prepr[[set]]$data,
   #stopAt = "30-VST",
   addToAnalysisManager = FALSE,
   idsAreSymbols = FALSE,
   minValue = 1,

   vst.design = "~donor_sex + specimen_type", 
   adj.calculateNPCs = 4,
   adj.removedCovariates = NULL,
   adj.calculateSVAFactors = FALSE,
   adj.sva.useNFactors = 1,
   adj.svaModel = "~ ",
   bw.groupBy = c("donor_sex", "specimen_type"),
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
   stp.width = 15,
   stp.height = 5,
   stp.mainSep = ", ",
   stp.marAll = c(1, 20, 3, 1),
   stp.maxCharPerLine = 60,
   plotDir = "Plots",

   pca.correctForColumn = NULL,
   pca.colorColumn = "specimen_type",
   pca.colorPrefix = "",
   pca.shapeColumn = "donor_sex",
   pca.shapePrefix = "",
   pca.mainSep = ", ",
   pca.legendWidth = 3.5,
   indent = 2, verbose = 5,
   saveDir.base = "../Data/TCGA"))

dir.create("RData", recursive = TRUE)
save(prepr, file = "RData/prepr.RData");
#load(file = "RData/prepr.RData");


multiExpr = mtd.apply(prepr, lastHE.VST.OR, "OSR", "expr.OSR");
multiCounts.all = mtd.apply(prepr, lastHE.VST.OR, "OSR", "xdata.OSR");
multiCounts = mtd.mapply(function(x, ref) x[, colnames(ref)], multiCounts.all, multiExpr);

multiWeights = mtd.apply(prepr, function(p) p$IOR$weightsForIOR);
weightFactors = mtd.apply(prepr, function(p) p$IOR$weightFactorsForIOR);

mtd.apply(multiWeights, function(x) table(x >0.1)/length(x))
mtd.apply(weightFactors, function(x) table(x >1)/length(x))

multiSampleAnnot = mtd.apply(prepr, lastHE.VST.OR, "OSR", "sampleAnnotation.OSR");
multiGeneAnnot = mtd.apply(prepr, lastHE.VST.OR, "HEFilter", "geneAnnotation.he");

names(multiCounts) = names(multiExpr) = names(multiSampleAnnot) = names(multiGeneAnnot) = setNames;
names(multiWeights) = names(weightFactors) = setNames;

numericCols = mtd.apply(multiSampleAnnot, sapply, is.numeric);

multiNumericPheno = mtd.subset(multiSampleAnnot, , numericCols[[1]]$data);
multiNumericPheno = mtd.subset(multiNumericPheno, , 
             multiGrep(c("SVAFactor.[2-9]", "PrincipalComponent.[2-9]"), mtd.colnames(multiNumericPheno), invert = TRUE));

save(multiCounts, multiExpr, multiSampleAnnot, multiGeneAnnot, multiWeights, #weightFactors,
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
    plot(multiExpr[[set]]$data[, g], weightFactors[[set]]$data[, g],
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
    plot(x, multiExpr[[set]]$data[, g], 
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







   





