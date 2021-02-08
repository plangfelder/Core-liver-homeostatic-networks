source("../Functions/networkFunctions-extras-19.R")
source("../Functions/labelPoints2-01.R");
source("../Functions/outlierRemovalFunctions.R")
source("../Functions/preprocessing-General-013.R")
source("../Functions/GNVFunctions-016.R")
source("../Functions/individualAnalysis-General-008.R");

library(anRichment)

dir.create("RData", recursive = TRUE);
dir.create("Results", recursive = TRUE);
dir.create("Plots", recursive = TRUE);

# Load data

x = load("../110-Preprocessing-TCGA/RData/preprocessedData.RData");

#prettifyList = loadAsList("../120-IndividualAnalysis-Suppli2019-HumanNAFLD-GSE126848/RData/prettifyList2.RData")[[1]]

traitNames = mtd.colnames(multiSampleAnnot)
traitNames2 = traitNames[traitNames!=make.names(traitNames)];
order = order(-nchar(traitNames2));
prettifyList1 = data.frame(from = make.names(traitNames2[order]), to = traitNames2[order])

prettifyList2 = data.frame(from = c(".for.", ".in.", ".vs."), to = c(" for ", " in ", " vs. "));

pretifyList = rbind(prettifyList1, prettifyList2);

nSets = length(multiExpr);
tissues = rep("Liver", nSets)
ltissues = tolower(tissues);

library(Cairo);
library(doParallel)

# Add one more trait to multiSampleAnnot

multiSampleAnnot2 = mtd.apply(multiSampleAnnot, function(sa)
{
  if ("Etiology" %in% names(sa))
  {
    unclear = is.na(sa$Etiology) | replaceMissing(sa$Etiology=="W/O Etiology") | grepl(".+ METAB", replaceMissing(sa$Etiology));
    metab = replaceMissing(sa$Etiology=="METAB");
    metab[unclear] = NA;
    sa$MetabEtiology = 0 + metab;
  }
  sa;
})


expr.flat = multiExpr;
sampleAnnot.flat = multiSampleAnnot2;
weights.flat = multiWeights;

setNames.flat = setNames;
nSets.flat = length(expr.flat);
names(expr.flat) = names(weights.flat) = names(sampleAnnot.flat) = setNames

cor = WGCNA::cor;

setsByTissue = tapply(setNames, tissues, identity)

setNames.pretty.flat = setNames;
lsetNames.pretty.flat = setNames;

#========================================================================================================
#
# Scale-free topology
#
#========================================================================================================

enableWGCNAThreads(4);

powers = c(2:10);
if (!checkRDataFile("RData/sft.RData"))
{
   system.time({sft = mtd.mapply(pickSoftThreshold, expr.flat, weights = weights.flat, 
                          MoreArgs = list(dataIsExpr = TRUE, corFnc = cor,
                          networkType = "signed hybrid", blockSize = 5000,
                         powerVector = powers , verbose = 2),
                         mdmaVerbose = 1 )});
   save(sft, file = "RData/sft.RData");
   names(sft) = setNames.flat;
   indices = mtd.apply(sft, getElement, "fitIndices");
   R.sq = mtd.apply(indices, getElement, "SFT.R.sq", mdaSimplify = TRUE);
   mean.k = mtd.apply(indices, getElement, "mean.k.", mdaSimplify = TRUE);
   median.k = mtd.apply(indices, getElement, "median.k.", mdaSimplify = TRUE);

   pdf(file = "Plots/scaleFreeTopologyAnalysis.pdf", wi = 12, he = 4);
   #sizeGrWindow(12, 9);
   par(mfrow = c(1,3))
   par(mar = c(3.5, 3.5, 1.5, 1))
   par(mgp = c(2.1, 0.8, 0))
   par(cex = 1);

     colors = 1:nSets;
     multiPlot(powers, R.sq, pch = 21, col = 1, bg = colors, cex = 1.7,
               main = "Scale-free topology fit index", xlab = "Soft-thresholding power",
               ylab = "Fit R squared");
     addGrid(v = TRUE, linesPerTick = 2)
     legendClean("bottomright", legend = setNames.pretty.flat, pch = 21, col = 1,  pt.bg = colors,
                 pt.cex = 1, cex = 0.8)
     multiPlot(powers, log10(mean.k), pch = 21, col = 1, bg = colors, cex = 1.7,
               main = "Mean connectivity", xlab = "Soft-thresholding power",
               ylab = "log10(Mean connectivity)");
     addGrid(v = TRUE, linesPerTick = 2)
     legendClean("topright", legend = setNames.pretty.flat, pch = 21, col = 1,  pt.bg = colors,
                 pt.cex = 1, cex = 0.8)
     multiPlot(powers, log10(median.k), pch = 21, col = 1, bg = colors, cex = 1.7,
               main = "Median connectivity", xlab = "Soft-thresholding power",
               ylab = "Median connectivity");
     addGrid(v = TRUE, linesPerTick = 2)
     legendClean("topright", legend = setNames.pretty.flat, pch = 21, col = 1,  pt.bg = colors,
                 pt.cex = 1, cex = 0.8)
   dev.off();
}

setNames.short = setNames.flat;

#========================================================================================================
#
# Network analysis
#
#========================================================================================================

# Will use power 5
softPowers = rep(5, nSets.flat);

networkOptions = list2multiData(lapply(softPowers, function(power)
   newNetworkOptions(
      correlationOptions = newCorrelationOptions(
          corType = "pearson", maxPOutliers = 0.05,
          corOptions = list(use = 'p', nThreads = 0)),
      power = power,
      networkType = "signed hybrid",
      TOMType = "signed",
      suppressNegativeTOM = FALSE,
      TOMDenom = "mean")));

names(networkOptions) = setNames.flat;


#====================================================================================================
#
# Define the analyses. Must calculate TOMs separately because the genes in each data set are different.
#
#====================================================================================================

consensusOptions = newConsensusOptions(
   calibration = "none");

analysisTrees = lapply(setNames.flat, function(.set)
  humanLiver = newConsensusTree(
      consensusOptions = consensusOptions,
      inputs = .set, 
      analysisName = .set)
);

analysisNames = sapply(analysisTrees, getElement, "analysisName");

analysisNames.pretty = analysisNames;
analysisSets = lapply(analysisTrees, consensusTreeInputs);

nAnalyses = length(analysisSets)

expr.ana = lapply(analysisSets, function(sets) mtd.restrictToCommonColumns(expr.flat[sets]));

weights.ana = mymapply(function(sets, .multiExpr) 
{
  mtd.subset(weights.flat[sets], , mtd.colnames(.multiExpr), permissive = TRUE)
}, analysisSets, expr.ana);


sampleAnnot.ana = lapply(analysisSets, function(sets) sampleAnnot.flat[sets]);

geneAnnotation.ana = lapply(expr.ana, function(.multiExpr) 
{
   ga = multiGeneAnnot[[ names(.multiExpr)[1] ]]$data;
   ga[ match(mtd.colnames(.multiExpr), ga$GeneID), ];
});

networkOptions.ana = lapply(analysisSets, function(sets) networkOptions[sets]);

TOMDir.indiv= lapply(analysisNames, function(name) spaste("RData/TOMs/Individual-", name));
lapply(TOMDir.indiv, dir.create, recursive = TRUE);

gc();
if (!checkRDataFile("RData/individualTOMinfo.RData"))
{
  individualTOMInfo = mymapply( individualTOMs,
        expr.ana, multiWeights = weights.ana, networkOptions = networkOptions.ana,
        individualTOMFileNames = spaste(TOMDir.indiv, "/individualTOM-Set%s-Block%b.RData"),
        MoreArgs = list(
          maxBlockSize = 45000,
          saveTOMs = TRUE,
          verbose = 5));

  save(individualTOMInfo, file = spaste("RData/individualTOMinfo.RData"));
}
#=======================================================================================================
#
# Read in results of standard analysis to get gene significance;
#
#=======================================================================================================

# Read in the meta-analysis results
individualAssoc1 = lapply(setNames.flat, function(name)
  read.csv(gzfile(spaste("../120-IndividualAnalysis-TCGA/Results/01-Association/",
                         "associationWithPhenotypes-", make.names(name), ".candidateCovars-all.csv.gz")), check.names = FALSE));

individualAssoc2 = lapply(setNames.flat, function(name)
{
  out = try(read.csv(gzfile(spaste("../120-IndividualAnalysis-TCGA/Results/01-Association/",
                         "associationWithPhenotypes-", make.names(name), ".selected-all.csv.gz")), check.names = FALSE));
  if (inherits(out, "try-error")) return(NULL);
  out;
})

individualAssoc = mymapply (function(x, y) if (!is.null(y)) data.frame(x, y, check.names = FALSE) else x, 
    individualAssoc1, individualAssoc2)


names(individualAssoc) = setNames;

# For each analysis, put together the Q association from the relevant tissues

individualResults.ana = mymapply(function(.ga, .sets)
{
  IDs = .ga$GeneID;
  out = do.call(cbind, removeListNames(lapply(individualAssoc[.sets], function(ia)
     ia[ match(IDs, ia$GeneID, incomparables = NA), grep("Z for ", names(ia))])))
  as.data.frame(lapply(out, clipData, max = 8), check.names = FALSE);
}, geneAnnotation.ana, analysisSets);

individualSignificance.ana = mymapply(function(.ga, .sets)
{
  IDs = .ga$Entrez;
  out = do.call(cbind, removeListNames(lapply(individualAssoc[.sets], function(ia)
     ia[ match(IDs, ia$Entrez), grep("FDR for ", names(ia))])))
  0 + (out < 0.05)
}, geneAnnotation.ana, analysisSets);

individualColors.ana = lapply(individualResults.ana, numbers2colors, commonLim = TRUE, lim = c(-7, 7));

plotTraitNames = lapply(individualResults.ana, function(.res)
{
  multiSub(c("Z for "), c(""), colnames(.res));
})

#==================================================================================================
#
# Consensus calculation
#
#==================================================================================================

if (!checkRDataFile("RData/consensusTOMInfo.RData"))
{
  dirs = spaste("RData/TOMs/noCalibration-", analysisNames);
  lapply(dirs, dir.create, recursive = TRUE);

  system.time({ consensusTOMInfo = mymapply(hierarchicalConsensusTOM,
     consensusTree = analysisTrees,
     individualTOMInfo = individualTOMInfo,
     consensusTOMFilePattern = spaste(dirs, "consensusTOM-%a-Block%b.RData"),
     MoreArgs = list(
       useDiskCache = FALSE,
       keepIntermediateResults = FALSE,
       verbose = 5))});

  save(consensusTOMInfo, file = "RData/consensusTOMInfo.RData");
}

#===============================================================================================================
#
# Load the mouse network labels to use as reference labels.
#
#===============================================================================================================

kmeOld = read.csv(gzfile(
  "../060-NetworkAnalysis-Suppli2019-HumanNAFLD-GSE126848/Results/networkAnalysisResults-humanLiver.main.csv.gz"));

oldLabels0 = data.frame(Entrez = kmeOld$Entrez, Label = kmeOld$module);

oldLabels = oldLabels0[ !is.na(oldLabels0$Entrez) & oldLabels0$Label!=0, ];

nSets.ana = sapply(analysisSets, length)
detectCutHeight = 0.998;
deepSplit = 2.5;
stop(10)
#===============================================================================================================
#
# Module analysis
#
#===============================================================================================================

mods.ns = list();
for (ana in 1:nAnalyses)
 mods.ns[[ana]] = hierarchicalConsensusModules(
     consensusTOMInfo = consensusTOMInfo[[ana]],
     multiExpr = expr.ana[[ana]],
     multiWeights = weights.ana[[ana]],
     checkMissingData = FALSE,
     detectCutHeight = detectCutHeight,
     deepSplit = deepSplit,
     minModuleSize = 30,
     useBranchEigennodeDissim = TRUE,
     mergeCutHeight = 0.2,
     minCoreKME = 0.4, minKMEtoStay = 0.2,
     verbose = 5)

save(mods.ns, file = "RData/mods.ns.RData");

#load(file = "RData/mods.ns.RData");

labels.ns.0 = lapply(mods.ns, getElement, "labels");

labels.original.ns = mymapply(function(l, .geneAnnot) data.frame(Entrez = .geneAnnot$Entrez, Label = l),
                              labels.ns.0, geneAnnotation.ana);


# Save the original labels 
save(labels.original.ns, file = "RData/originalLabels-labels.original.ns.RData");


labels.ns = lapply(labels.original.ns, matchLabelsToReference.2, oldLabels);

save(labels.ns, file = "RData/labels.ns.RData");
#load(file = "RData/labels.ns.RData");

labels.ns = list2multiData(labels.ns);

colors.ns = mtd.apply(labels.ns, labels2colors)

#========================================================================================================
#
# First look at dendrograms
#
#========================================================================================================

plotColors.ns = mtd.mapply(cbind, colors.ns, individualColors.ana)

plotTrees.ns = lapply(mods.ns, function(md) clipTree(md$dendrograms[[1]], p = 0.01));

nPlotPheno = lapply(plotTraitNames, length)

pdf(file = "Plots/geneDendrogramsAndGS-withoutStabilityAnalysis.pdf", wi = 14, he = 9);
#sizeGrWindow(14, 9);
par(lheight = 0.85)
for (ana in 1:nAnalyses)
plotDendroAndColors(plotTrees.ns[[ana]],
                    plotColors.ns[[ana]]$data,
                    c("Modules", plotTraitNames[[ana]]),
                    autoColorHeight = FALSE,
                    colorHeight = 0.5,
                    rowText = labels.ns[[ana]]$data, 
                    textPositions = 1,
                    cex.rowText = 0.8,
                    rowWidths = c(1, 10, rep(1, nPlotPheno[[ana]])),
                    rowTextAlignment = "center",
                    dendroLabels = FALSE,
                    main = analysisNames.pretty[ana],
                    marAll = c(0, 7, 2, 0), hang = 0.02, guideHang = 0.08,
                    cex.colorLabels = 0.8,
                    addGuide = TRUE);

dev.off();

#========================================================================================================
#
# save data for stability analysis
#
#========================================================================================================

dir.create("../129-StabilityAnalysis-TCGA/RData", recursive = TRUE);

save(expr.ana, weights.ana, analysisTrees, networkOptions.ana, detectCutHeight, deepSplit,
     file = "../129-StabilityAnalysis-TCGA/RData/dataForStabilityAnalysis.RData");


#========================================================================================================
#
# Re-calculate modules with the input of sampled labels.
#
#========================================================================================================

nRuns = 50;
sampledLabels.lst = lapply(1:nAnalyses, function(ana)
  lapply(1:nRuns, function(run)
  {
     out = try(loadAsList(spaste("../129-StabilityAnalysis-TCGA/RData/sampledLabels-run",
               prependZeros((ana-1)*nRuns + run-1, 3),
                ".RData"))[[1]]);
     if (inherits(out, "try-error")) NULL else out;
  })
);

sampledLabels = lapply(sampledLabels.lst, function(x) do.call(cbind, x))

lapply(sampledLabels, dim)
mods.lst = list();

minStabilityDissim = list(
    main = rep(0.6, nAnalyses),
    sub = c(0.3, 0.3, 0.3));

for (size in names(minStabilityDissim))
{
  mods.lst[[size]] = list();
  for (ana in 1:nAnalyses)
  {
    gc();
    mods.lst[[size]] [[ana]] = hierarchicalConsensusModules(
        consensusTOMInfo = consensusTOMInfo[[ana]],
        multiExpr = expr.ana[[ana]],
        multiWeights = weights.ana[[ana]],
        checkMissingData = FALSE,
        detectCutHeight = detectCutHeight,
        deepSplit = deepSplit,
        minModuleSize = 30,
        useBranchEigennodeDissim = TRUE,
        stabilityLabels = sampledLabels[[ana]],
        minStabilityDissim = minStabilityDissim[[size]] [ana],
        mergeCutHeight = 0.2,
        minCoreKME = 0.4, minKMEtoStay = 0.2,
        verbose = 5)
   
  }
}
gc();

save(mods.lst, file = "RData/mods.lst.RData");
#load(file = "RData/mods.lst.RData");

mods = mods.lst[["main"]];

names(mods) = analysisNames;
save(mods, file = "RData/mods.RData");

#load(file = "RData/mods.RData");

labels.0 = lapply(mods, getElement, "labels");

labels.original = mymapply(function(l, .geneAnnot) data.frame(Entrez = .geneAnnot$Entrez, Label = l),
                              labels.0, geneAnnotation.ana);

# Save the original labels for calculation of reference labels
save(labels.original, file = "RData/originalLabels-labels.original.RData");

labels = lapply(labels.original, matchLabelsToReference.2, oldLabels);

labels = list2multiData(labels);
names(labels) = analysisNames;

save(labels, file = "RData/labels.RData");
#load(file = "RData/labels.RData");

colors = mtd.apply(labels, labels2colors)

#========================================================================================================
#
# Match the labels of sub-modules to their main modules
#
#========================================================================================================

subLabels0 = lapply(mods.lst[["sub"]], getElement, "labels");
overlaps = mtd.mapply(overlapTable, labels, subLabels0, MoreArgs = list(ignore = 0), mdmaVerbose = TRUE);

subLabels = mtd.mapply(function(.overlap, .subLabels)
{
  maxRow = apply(.overlap$pTable, 2, which.min);
  modSizes = table(.subLabels);
  oldLabels = colnames(.overlap$pTable);
  modSizes = modSizes[match(oldLabels, names(modSizes))];
  newLabels = rownames(.overlap$pTable)[maxRow];
  for (r in 1:nrow(.overlap$pTable))
  {
    mods1 = which(maxRow==r);
    if (length(mods1) > 1)
    {
      sizes1 = modSizes[ mods1 ];
      order = order(-sizes1);
      newLabels[mods1] = spaste(newLabels[mods1], ".", prependZeros(1:length(mods1)));
    }
  }
  translateUsingTable(.subLabels, cbind(as.numeric(oldLabels), as.numeric(newLabels)), keepUntranslated = TRUE)
}, overlaps, subLabels0);

maxLabel = max(range(labels))

labels.x = c(labels, subLabels);
labelNames = c(outer(analysisNames, c(".main", ".sub"), FUN = spaste));
labelNames.pretty = c(outer(analysisNames.pretty, c(", main modules", ", sub-modules"), FUN = spaste));
nLabels = length(labelNames);
names(labels.x) = labelNames;

subAnalyses = grep(".sub$", labelNames);

labelTranslationForColors = mtd.apply(labels.x, function(.labels)
{
  fractional1 = .labels-0.1 == round(.labels-0.1)
  fractional2 = .labels!=round(.labels) & !fractional1
  fractionalLevels1 = sort(unique(.labels[fractional1]));
  fractionalLevels2 = sort(unique(.labels[fractional2]));
  integralLevels = sort(unique(.labels[!fractional1 & !fractional2]));
  out = cbind(fractional = integralLevels, integral = integralLevels);
  if (length(fractionalLevels1) > 0)
    out = rbind(out, cbind(fractional = fractionalLevels1, integral = fractionalLevels1-0.1));
  if (length(fractionalLevels2)>0)
    out = rbind(out, cbind(fractional = fractionalLevels2, integral = (maxLabel + 1):(maxLabel + length(fractionalLevels2))));
  out;
})

labelsForColors = mtd.mapply(translateUsingTable, labels.x, labelTranslationForColors, MoreArgs = list(keepUntranslated = TRUE));
subLabelsForColors = labelsForColors[subAnalyses];

xxxx = mtd.mapply(function(x, y) print(table(x, y)), subLabels[1], subLabelsForColors[1])

subColors = mtd.apply(subLabelsForColors, labels2colors);

save(labels.x, labelNames, labelNames.pretty, nLabels, labelsForColors, subLabelsForColors, subColors,
   labelTranslationForColors,
   file = "RData/labels.x.etc.RData");

#load("RData/labels.x.etc.RData")

#========================================================================================================
#
# Save data for module preservation
#
#========================================================================================================

save(expr.ana, weights.ana, labels.x, analysisNames, labelNames, 
     analysisNames.pretty, labelNames.pretty,
     geneAnnotation.ana, sampleAnnot.ana,
     file = "RData/dataForModulePreservation-TCGA.RData");
    
#========================================================================================================
#
# First look at dendrograms
#
#========================================================================================================

plotColors = mtd.mapply(cbind, subColors, colors, individualColors.ana)
#plotColors = mtd.mapply(cbind, subColors, colors);
plotTrees = lapply(mods, function(md) clipTree(md$dendrograms[[1]], p = 0.01));
nPlotPheno = lapply(plotTraitNames, length)

pdf(file = "Plots/geneDendrogramsAndGS.pdf", wi = 14, he = 8);
#sizeGrWindow(14, 9);
par(lheight = 0.85)
for (ana in 1:nAnalyses)
plotDendroAndColors(plotTrees[[ana]],
                    plotColors[[ana]]$data,
                    c("Sub-modules", "Modules", formatLabels(plotTraitNames[[ana]], maxCharPerLine = 120)),
                    autoColorHeight = FALSE,
                    colorHeight = 0.5,
                    rowText = labels[[ana]]$data, 
                    textPositions = 2,
                    cex.rowText = 0.8,
                    rowWidths = c(1, 1, 4, rep(1, nPlotPheno[[ana]])),
                    #rowWidths = c(1, 1, 4),
                    rowTextAlignment = "center",
                    dendroLabels = FALSE,
                    main = spaste("Gene clustering and modules, ", analysisNames.pretty[ana]),
                    marAll = c(0, 20, 1.5, 0), hang = 0.02, guideHang = 0.08,
                    cex.colorLabels = 0.8,
                    addGuide = TRUE);

dev.off();

# Plot gene dendrograms and indicator of significance (FDR<0.1)
plotColors2 = mtd.mapply(cbind, subColors, colors, individualSignificance.ana)

pdf(file = "Plots/geneDendrogramsAndGS-significanceIndicators.pdf", wi = 14, he = 8);
#sizeGrWindow(14, 9);
par(lheight = 0.85)
for (ana in 1:nAnalyses)
plotDendroAndColors(plotTrees[[ana]],
                    plotColors2[[ana]]$data,
                    c("Sub-modules", "Modules", plotTraitNames[[ana]]),
                    autoColorHeight = FALSE,
                    colorHeight = 0.5,
                    rowText = labels[[ana]]$data,
                    textPositions = 2,
                    cex.rowText = 0.8,
                    rowWidths = c(1, 1, 4, rep(1, nPlotPheno[[ana]])),
                    rowTextAlignment = "center",
                    dendroLabels = FALSE,
                    main = analysisNames.pretty[ana],
                    marAll = c(0, 20, 1.5, 0), hang = 0.02, guideHang = 0.08,
                    cex.colorLabels = 0.8,
                    addGuide = TRUE);

dev.off();

#====================================================================================================
#
# Enrichment analysis
#
#====================================================================================================

x = load("../050-IndividualAnalysis-Suppli2019-HumanNAFLD-GSE126848/RData/stdCollections.RData");

mouseColl = loadAsList("../020-IndividualAnalysis/RData/collectionOfDEGenes-Liver.RData")[[1]];

allCollections = c(stdCollections, list(mouseLiver = convertCollectionToOrganism(mouseColl, organism = "human")))


allCollections.combined = do.call(mergeCollections, allCollections);

identifiers = lapply(geneAnnotation.ana, getElement, "Entrez");
identifiers.x = c(identifiers, identifiers);


enr = mtd.mapply(function(lab1, ids)
     enrichmentAnalysis(classLabels= lab1, identifiers = ids,
                     refCollection = allCollections.combined, ignoreLabels = 0, threshold = 5e-2,
                     nBestDataSets = 5,
                     getOverlapSymbols = TRUE, getOverlapEntrez = FALSE,
                     getDataSetDetails = FALSE, useBackground = "given",
                     maxReportedOverlapGenes = 10000),
    labels.x, ids = identifiers.x, mdmaVerbose = TRUE)

save(enr, file = "RData/enrichment-enr.RData");

load(file = "RData/enrichment-enr.RData");

names(enr) = labelNames;

enrichmentTables.all = mtd.apply(enr, getElement, "enrichmentTable")

enrichmentTables.specific = mtd.apply(enrichmentTables.all, 
   splitEnrichmentTableBySubcollections,
      allCollections.combined, allCollections,
      dropColumns = "dataSetID");

nc = length(enrichmentTables.specific[[1]]$data);

nCollSets = length(allCollections)
for (ana in 1:nLabels) 
{
  dir.create(spaste("Results/EnrichmentAnalysis/", labelNames[ana]), recursive = TRUE);
  for (ic in 1:nc)
    write.csv(signifNumeric(enrichmentTables.specific[[ana]]$data[[ic]], 2),
      file = gzfile(spaste("Results/EnrichmentAnalysis/", labelNames[ana], "/enrichmentOfModules-", labelNames[ana], "-",
                 prependZeros(ic, 2), "-",
                 names(enrichmentTables.specific[[ana]]$data)[ic], ".csv.gz")),
      row.names = FALSE, quote = TRUE);

  write.csv(signifNumeric(enrichmentTables.all[[ana]]$data, 2),
        file = gzfile(spaste("Results/EnrichmentAnalysis/", labelNames[ana], "/enrichmentOfModules-", labelNames[ana],
                 "-00-CombinedCollection.csv.gz")),
        row.names = FALSE, quote = TRUE);
}

# For eLabels: create an enrichment table that excludes WGCNA and DE for mHtt sets

eLabels = mtd.apply(enrichmentTables.all, enrichmentLabels,
  focusOnGroups = c("DE genes in liver of mice on various diets", 
                    "KEGG", "GO.BP", "GO.MF", "GO.CC", "REACTOME", "Genomic position", "MSigDB", 
                    "ChEA", "ENCODE Histone Modification", "mirTarBase", "all"),
  groupShortNames = c("mouse liver", "KEGG", "GO.BP", "GO.MF", "GO.CC", "Reactome", "GenPos", "MSigDB", 
                    "ChEA", "HistMod", "mirTar", "all"))

eLabels = mtd.apply(eLabels, function(x) setColnames(x, sub("^class", "module", colnames(x))));
eLabels = mtd.apply(eLabels, shortenEnrichmentLabels);

eLabels2 = mtd.mapply(function(.eLabels, .transTable)
  cbind(.eLabels, labelForColor = translateUsingTable(.eLabels$module, .transTable, keepUntranslated = TRUE)),
  eLabels, labelTranslationForColors);

save(eLabels, eLabels2, file = "RData/eLabels.RData");
load(file = "RData/eLabels.RData")

#========================================================================================================
#
# Get KME values
#
#========================================================================================================

# Consensus KME analysis
expr.ana.x = c(expr.ana, expr.ana);
weights.ana.x = c(weights.ana, weights.ana);
analysisSets.x = c(analysisSets, analysisSets);
analysisTrees.x = c(analysisTrees, analysisTrees);
geneAnnotation.ana.x = c(geneAnnotation.ana, geneAnnotation.ana);
sampleAnnot.ana.x = c(sampleAnnot.ana, sampleAnnot.ana)

names(expr.ana.x) = names(weights.ana.x) = names(analysisSets.x) = labelNames;
names(analysisTrees.x) = names(geneAnnotation.ana.x) = names(sampleAnnot.ana.x) = labelNames;


conKME = mtd.mapply(hierarchicalConsensusKME,
                multiExpr = expr.ana.x,
                multiWeights = weights.ana.x,
                moduleLabels = labels.x,
                setNames = analysisSets.x,
                consensusTree = analysisTrees.x,
                additionalGeneInfo = geneAnnotation.ana.x,
            MoreArgs = list(
                corAndPvalueFnc = "corAndPvalue", corComponent = "cor",
                corOptions = list(use = 'p'),
                getSetZ = TRUE, getSetP = TRUE, 
                useRankPvalue = FALSE,
                reportWeightType = "rootDoF",
                includeWeightTypeInColnames = FALSE,
                includeID = FALSE));

conKME = mtd.apply(conKME, function(.conKME) setNames(.conKME, multiSub(c("[Cc]onsKME", ".LI[HCR][ICA].*", "consensus."), 
                                                                        c("kME", "", ""),
                                    names(.conKME))));

for (ana in 1:nLabels)
{
  write.csv(signifNumeric(conKME[[ana]]$data, 3),
          file = gzfile(spaste("Results/networkAnalysisResults-", labelNames[ana], ".csv.gz")),
          row.names = FALSE, quote = TRUE);
}


3
# Get module hubs

conKME2 = mtd.apply(conKME, function(.conKME) 
   .conKME[, multiGrep(c("Symbol", "module$", "^Z.kME[0-9]"),
                       names(.conKME), value = TRUE, sort = FALSE)])

moduleHubs = mtd.apply(conKME2, 
                    consensusHubs,
                    nHubs = 20, selectBy = "Z.kME",
                    select = "positive", returnColumns = "Symbol",
                    simpleOutput = TRUE);

moduleHubs.ps = mtd.apply(moduleHubs, apply, 2, base::paste, collapse = ", ");

save(moduleHubs, moduleHubs.ps, file = "RData/moduleHubs.RData");

load(file = "RData/moduleHubs.RData");

#========================================================================================================
#
# Relationship to traits
#
#========================================================================================================

meInfo = mtd.mapply(cbind, eLabels, topConsensusHubsByZ.KME = moduleHubs.ps);

MEs = mtd.mapply(multiSetMEs,  expr.ana.x, universalColors = labels.x, 
                    MoreArgs = list(excludeGrey = TRUE));

for (ana in 1:nLabels)
  write.csv.nr(dataForTable(MEs[[ana]]$data [[1]]$data, transpose = FALSE, IDcolName = "SampleID"),
     file = spaste("Results/eigengenes-", labelNames[ana], ".csv"));


 

modelInfo = loadAsList("../120-IndividualAnalysis-TCGA/RData/sampleAnnotation-designInfo-prettifyList.RData");


names(modelInfo$multiSampleAnnot.present) = setNames.flat

sampleAnnot.model.flat = mtd.mapply(cbind, 
   modelInfo$multiSampleAnnot.covars, modelInfo$multiSampleAnnot.selected, modelInfo$multiSampleAnnot.otherTraits);

prettifyList = modelInfo$prettifyList5;
designInfo.covars = mtd.apply(modelInfo$designInfo.covars, function(di) di[1:2, ]);

designInfo.flat = mtd.mapply(function(...)
{
  args = list(...);
  args = args[sapply(args, length) > 0];
  args = lapply(args, function(x) 
  {
     out = x[names(x) != "reduced"];
     out$design = sub("donor_sex", "M.vs.F", out$design);
     out$coefficientName = sub("donor_sex.*", "M.vs.F", out$coefficientName);
     out;
   });
 
  do.call(rbind, args);
}, designInfo.covars, modelInfo$designInfo.selected, modelInfo$designInfo.otherTraits)

names(sampleAnnot.model.flat) = names(designInfo.flat) = setNames.flat;

# Shortcut
sampleAnnot.model.ana = lapply(analysisSets, function(sets) mtd.restrictToCommonColumns(sampleAnnot.model.flat[sets]));

sampleAnnot.model.ana.x = c(sampleAnnot.model.ana, sampleAnnot.model.ana);
designInfo.ana.x = c(designInfo.flat, designInfo.flat)

names(sampleAnnot.model.ana.x) = names(designInfo.ana.x) = analysisNames;

ss.ME = mtd.mapply(function(.MEs, .preds, .meInfo, .designInfo)
{
  ss = prettifyNames(associationScreening(
          response = .MEs[[1]]$data,
          predictors = .preds[[1]]$data,
          full = .designInfo$design,
          termsOfInterest = .designInfo$coefficientName,
          reduced = NULL,
          testName = .designInfo$testName,
          additionalResponseInfo = .meInfo,
          addID = TRUE, getFDR = TRUE,
          FDRFromCombinedPValues = TRUE,
          getCLS = FALSE), prettifyList);
  ss;
}, MEs, sampleAnnot.model.ana.x, meInfo, designInfo.ana.x, mdmaVerbose = TRUE);



for (ana in 1:nLabels)
  write.csv(signifNumeric(ss.ME[[ana]]$data, 3),
            file = gzfile(spaste("Results/moduleSummaryAndAssociationWithTraits-", labelNames[ana], ".csv.gz")),
            quote = TRUE, row.names = FALSE);

save(ss.ME, file = "RData/ss.ME.RData");

#================================================================================================
#
# Plot the relationships to traits
#
#================================================================================================

# Start with relationship in most detail.

zPattern = c("Z for ");
matPattern = c("Z for ");
pPattern = "FDR for";

labelPlotDirs = file.path("Plots", labelNames);
lapply(labelPlotDirs, dir.create, recursive = TRUE)

maxRowsPerPage = 30;
maxColsPerPage = 20;
for (ana in 1:nLabels)
{
  matCols = multiGrep(matPattern, colnames(ss.ME[[ana]]$data));
  pCols = multiGrep(pPattern, colnames(ss.ME[[ana]]$data));
  mat = as.matrix(ss.ME[[ana]] $ data[, matCols, drop = FALSE]);
  ZCols = multiGrep(zPattern, colnames(ss.ME[[ana]]$data));
  mat.plot = mat;
  max = max(abs(mat.plot), na.rm = TRUE);
  mat.text = round(as.matrix(ss.ME[[ana]]$data[, ZCols]), 1);

  p = ss.ME[[ana]]$data[, pCols];

  textMat = spaste(c(as.matrix(mat.text)), "\n", signif(as.matrix(p), 1) );
  dim(textMat) = dim(mat);
  textMat[ is.na(p) | replaceMissing(p > 0.05)] = "";
  
  xlabs = multiSub(c("Z for "), c(""), colnames(mat));

  nc = ncol(mat);
  colsPerPage = allocateJobs(nc, ceil(nc/maxColsPerPage));
  nc1 = max(sapply(colsPerPage, length));
  nr = nrow(mat)
  rowsPerPage = allocateJobs(nr, ceil(nr/maxRowsPerPage));
  nr1 = max(sapply(rowsPerPage, length));
  wi = 7 + 0.5 * nc1; he = 3 + 0.4 * nr1
  pdf(file = spaste(labelPlotDirs[ana], "/Association-MEs-vs-traits-", labelNames[ana], "-all.pdf"), 
      wi = wi, he = he)
  #sizeGrWindow(wi, he);
  par(mar = c(22, 35, 2.8, 1));
  par(lheight = 0.85);
  labeledHeatmap.multiPage(mat.plot,
             xLabels = xlabs,
             yLabels = spaste("ME", labels2colors(
                translateUsingTable(as.numeric(eLabels[[ana]]$data$module), labelTranslationForColors[[ana]]$data))),
             ySymbols = formatLabels(eLabels[[ana]]$data$enrichmentLabel, maxCharPerLine = 130, maxLines = 2,
                                     split = "[ _]", fixed = FALSE, newsplit=" ", keepSplitAtEOL = FALSE),
             textMatrix = textMat,
             cex.text = 0.8,
             colors = blueWhiteRed(100, gamma = 1),
             zlim = c( -max, max),
             setStdMargins = FALSE,
             main = spaste("Module association with disease contrasts\n", labelNames.pretty[ana]),
             cex.lab.y = 0.7, cex.lab.x = 0.8, colsPerPage = colsPerPage, rowsPerPage = rowsPerPage);

  dev.off();
}

#================================================================================================
#
# Show boxplots of module eigengenes vs. disease
#
#================================================================================================

cex.main = 0.8;
minWidth = 5;
for (ana in 1:nLabels)
{
  sa1 = modelInfo$multiSampleAnnot.present[[analysisSets.x[[ana]] ]]$data;
  plotTraits = c("donor_age_at_diagnosis", multiGrep(c("[^r]_sex$", "_id$", "donor_vital_status",
                   "disease_status_last_followup", "_age_", "PrincipalComponent", "Sex", "donor_interval_of_last_followup",
                   "donor_age_at_enrollment", "Form.completion.date",
                    "Number of Samples Per Patient", "Birth from Initial Pathologic Diagnosis Date",
                 "ID$", "Age$", "Code$", "Sample Type", "Race.Category", "Center.of.sequencing",
                 "donor_tumour_stage_at_diagnosis", "donor_survival_time",
                 "American.Joint.Committee.on.Cancer.Publication.Version.Type",
                 "New.Neoplasm.Event.Post.Initial.Therapy.Indicator", "Informed.consent.verified",
                 "Ethnicity.Category", "Person.Gender", "Tissue.Source.Site_9",
                 "Tissue Prospective Collection Indicator_7", "Tissue Retrospective Collection Indicator_8",
                 "In PanCan Pathway Analysis"), names(sa1), value = TRUE, invert = TRUE));

  for (tr in plotTraits)
  {
    trv = sa1[[tr]];
    num = is.numeric(trv);
    if (num) 
    {
       width = minWidth
       height = 5; mar1 = 3.1; mar2 = 3.1;
    } else {
       width = max(minWidth, length(levels(factor(trv)))/2);
       fdt1 = factor(trv); levels1 = levels(fdt1); mw = marginWidth(levels1);
       height = 4 + 0.7 * mw$inch;
       mar1 = 0.7 * mw$lines + 0.5;
       mar2 = 3.1 + max(0, mw$lines - 10);
    }
    pdf(file = spaste(labelPlotDirs[ana], "/eigengenesVsTraits-", labelNames[ana], "-", make.names(tr), ".pdf"),
      wi = width, he = height);
    par(mar = c(mar1, mar2, if (num) 6 else 5, 0.5));
    par(mgp = c(2, 0.7, 0));
    mes1 = MEs[[ana]]$data [[1]]$data;
    nm = ncol(mes1);
    for (m in 1:nm)
    {
      mm = eLabels[[ana]]$data$module[m];
      x = MEs[[ana]]$data [[1]] $ data [[m]];
      if (num)
      {
         plot(trv, x, xlab = tr, ylab = spaste("ME", mm),
              pch = 21, col = "grey25", bg = "grey75",
              cex.main = 1, cex.lab = 1, cex.axis = 1, main = "");
         mainEnd = spaste("\ncor = ", round(cor(trv, x, use = 'p'), 2), 
                    ", p = ", signif(corAndPvalue(trv, x)$p, 1));
      } else {
         fdt1 = factor(trv);
         levels1 = levels(fdt1);
         xl = tapply(x, fdt1, identity);
         labeledBoxplot(xl, g = NULL, names = levels1, xlab = "", ylab = spaste("ME", mm), main = "",
                        addScatterplot = TRUE, notch = TRUE, pars = list(boxwex = 0.4, staplewex = 0.4, outwex = 0.4))
         mainEnd = "";
      }
      box = par("usr");
      width = box[2] - box[1];
      main = spaste(
        formatLabels(eLabels[[ana]]$data$enrichmentLabel[m], maxWidth = 0.92 * width, maxLines = 4, cex = cex.main, font = 2), 
          "\n",
        formatLabels(moduleHubs.ps[[ana]]$data[[m]], maxWidth = 0.92 * width, maxLines = 1, cex = cex.main, font = 2))
      title(spaste(main, mainEnd), cex.main = cex.main);
    }
    dev.off();
  }
}

#=====================================================================================================================
#
# Export sub-networks to cytoscape
#
#=====================================================================================================================

# For the top connections, keep 1% of all connections (hopefully about 10% of module genes)
# For smaller modules, need to keep more.

effectiveThreshold = function(topThreshold, moduleSize)
{
  (max(sqrt(topThreshold), (1000-moduleSize)/1000))^2
}

topThreshold = 0.01;

for (ana in 2:nLabels)
{
  dirs = spaste("CytoscapeFiles-", labelNames[ana], "/", c("TopConnections", "WholeModules", "AllModules-TopConnections"));
  lapply(dirs, dir.create, recursive = TRUE);

  ana.TOM = (ana-1) %% nAnalyses + 1;
  TOM = as.matrix(BD.getData(consensusTOMInfo[[ana.TOM]]$consensusData));
  gc()
  indivAssoc1 = individualAssoc[[analysisSets[[ ana.TOM]] ]]
  nodeAttr = data.frame(conKME[[ana]]$data[, c("GeneID", "Entrez", "Symbol", "Chr", "Loc", 
                                               "module", "kME.inOwnModule")],
                     check.names = FALSE);

  # Export top connections
  modules = setdiff(sort(unique(labels.x[[ana]]$data)), 0)
  topGenes = character(0);
  stopifnot(all.equal(conKME[[ana]]$data$GeneID, geneAnnotation.ana.x[[ana]]$GeneID))
  stopifnot(all.equal(conKME[[ana]]$data$GeneID, mtd.colnames(expr.ana.x[[ana]])))
  mxChar = max(nchar(modules))
  for (m in modules)
  {
    printFlush("Working on module", m);
    moduleGenes = which(labels.x[[ana]]$data==m);
    modTOM = TOM[moduleGenes, moduleGenes];
    colnames(modTOM) = geneAnnotation.ana.x[[ana]]$Gene[ moduleGenes];

    net = exportNetworkToCytoscape(modTOM,
      threshold = 0, nodeNames =  colnames(modTOM),
      nodeAttr = nodeAttr[moduleGenes, ]);

    write.table(net$edgeData, file = gzfile(spaste(dirs[2], "/M", prependZeros(m, mxChar), "-edgeFile.txt.gz")),
      col.names = TRUE, row.names = FALSE, quote = FALSE, sep = "\t");
    write.table(net$nodeData, file = gzfile(spaste(dirs[2], "/M", prependZeros(m, mxChar), "-nodeFile.txt.gz")),
      col.names = TRUE, row.names = FALSE, quote = FALSE, sep = "\t");
 
    threshold = quantile( as.dist(modTOM), prob = 1-effectiveThreshold(topThreshold, length(moduleGenes)));

    net.top = exportNetworkToCytoscape(modTOM,
      threshold = threshold, nodeNames =  colnames(modTOM),
      nodeAttr = nodeAttr[moduleGenes, ]);

    write.table(net.top$edgeData, file = gzfile(spaste(dirs[1], "/M", prependZeros(m, mxChar), 
                                                "-topConnections-edgeFile.txt.gz")),
      col.names = TRUE, row.names = FALSE, quote = FALSE, sep = "\t");
    write.table(net.top$nodeData, file = gzfile(spaste(dirs[1], "/M", prependZeros(m, mxChar), 
                                                "-topConnections-nodeFile.txt.gz")),
      col.names = TRUE, row.names = FALSE, quote = FALSE, sep = "\t");

    topGenes = c(topGenes, net.top$nodeData$Gene);
  }
  topIndex = match(topGenes, geneAnnotation.ana.x[[ana]]$Gene);
  topTOM = TOM[ topIndex, topIndex];
  net.top.all = exportNetworkToCytoscape(topTOM,
      threshold = 0, nodeNames =  topGenes,
      nodeAttr = nodeAttr[topIndex, ]);

  write.table(net.top.all$edgeData, file = gzfile(spaste(dirs[3], "/topGeneConnections-edgeFile.txt.gz")),
    col.names = TRUE, row.names = FALSE, quote = FALSE, sep = "\t");
  write.table(net.top.all$nodeData, file = gzfile(spaste(dirs[3], "/topGeneConnections-nodeFile.txt.gz")),
    col.names = TRUE, row.names = FALSE, quote = FALSE, sep = "\t");
}

    


#===================================================================================================================
#
# Add a quick prediction analysis for metabolic etiology in LICA-FR data. 
#
#===================================================================================================================

set = 1
keepSamples = !is.na(sampleAnnot.flat[[set]]$data$MetabEtiology);

expr1 = expr.flat[[1]]$data[keepSamples, ];
trait1 = sampleAnnot.flat[[set]]$data$MetabEtiology[keepSamples];
library(randomForest)
rf1 = randomForest(x = expr1, y = factor(trait1))

table(rf1$predicted, trait1);
am.rf = accuracyMeasures(rf1$predicted, trait1, negativeLevel = 0, positiveLevel = 1);

library(randomGLM)
rg1 = randomGLM(x = expr1, y = trait1);

table(rg1$predictedOOB, trait1)


# Try the same on the module eigengenes

ana = 1;
expr1 = MEs[[ana]]$data[[1]]$data[keepSamples, ];
trait1 = sampleAnnot.flat[[set]]$data$MetabEtiology[keepSamples];
library(randomForest)
rf1 = randomForest(x = expr1, y = factor(trait1))

table(rf1$predicted, trait1);
am.rf = accuracyMeasures(rf1$predicted, trait1, negativeLevel = 0, positiveLevel = 1);

library(randomGLM)
rg1 = randomGLM(x = expr1, y = trait1);

table(rg1$predictedOOB, trait1)

# The results aren't any better

#===================================================================================================================
#
# Module overlaps of TCGA modules with human and mouse modules
#
#===================================================================================================================

netFiles1 = spaste("../030-NetworkAnalysis/Results/networkAnalysisResults-Liver.", c("main", "sub"), ".csv.gz");
netNames1 = spaste("mouseLiver.", c("main", "sub"));
nets1 = lapply(netFiles1, read.csv);
names(nets1) = netNames1;

nets12 = lapply(nets1, function(.net)
{
  .net$Entrez = mapEntrez(.net$Entrez, orgFrom = "mouse", orgTo = "human");
  .net[!is.na(.net$Entrez), c("Entrez", "module")]
});

netFiles2 = spaste("../060-NetworkAnalysis-Suppli2019-HumanNAFLD-GSE126848/Results/networkAnalysisResults-humanLiver.", 
                     c("main", "sub"), ".csv.gz");
netNames2 = spaste("humanLiver.", c("main", "sub"));
nets2 = lapply(netFiles2, read.csv);
names(nets2) = netNames2;

nets22 = lapply(nets2, function(.net) .net[!is.na(.net$Entrez), c("Entrez", "module")]);

targetNets = c(nets12, nets22);

nTargets = length(targetNets);

targetNames.pretty = c("mouse liver main modules", "mouse liver sub-modules", "human NAFLD liver main modules",
   "human NAFLD liver sub-modules");

targetNames.short = c("mouse", "mouse", "NAFLD", "NAFLD");

for (ana in 1:nLabels) for (tg in 1:nTargets)
{
  targetLabels = translateUsingTable(geneAnnotation.ana.x[[ana]]$Entrez, targetNets[[tg]]);
  ot = overlapTable(targetLabels, labels.x[[ana]]$data, ignore = 0)

  t1 = table(targetLabels);
  t1 = t1[match(rownames(ot$pTable), names(t1))];

  t2 = table(labels.x[[ana]]);
  t2 = t2[match(colnames(ot$pTable), names(t2))];
  
  plotOverlapHeatmap(ot$countTable, ot$pTable,
     rowLabels = spaste("ME", labels2colors(as.numeric(rownames(ot$countTable)))),
     rowSymbols = spaste(targetNames.short[tg], " M", rownames(ot$countTable), " (", t1, ")"),
     colLabels = spaste("ME", labels2colors(as.numeric(colnames(ot$countTable)))),
     colSymbols = spaste(sub("\\..*", "", labelNames[ana]), " M", colnames(ot$countTable), " (", t2, ")"),
     main = spaste("Overlap of all ", labelNames.pretty[ana], " with ", targetNames.pretty[tg]),
     mar = c(9, 11, 1.5, 0.2),
     plotFile = spaste(labelPlotDirs[ana], "/overlapWith.", names(targetNets)[tg], ".pdf"),
     separatorInterval = 3);

  plotOverlapHeatmap(ot$countTable, ot$pTable,
     rowLabels = spaste("ME", labels2colors(as.numeric(rownames(ot$countTable)))),
     rowSymbols = spaste(targetNames.short[tg], " M", rownames(ot$countTable), " (", t1, ")"),
     colLabels = spaste("ME", labels2colors(as.numeric(colnames(ot$countTable)))),
     colSymbols = spaste(sub("\\..*", "", labelNames[ana]), " M", colnames(ot$countTable), " (", t2, ")"),
     main = spaste("Overlap of selected ", labelNames.pretty[ana], " with ", targetNames.pretty[tg]),
     threshold = 1e-4,
     mar = c(9, 11, 1.5, 0.2),
     plotFile = spaste(labelPlotDirs[ana], "/overlapWith.", names(targetNets)[tg], "-significant.pdf"),
     separatorInterval = 3);

}
     

if (FALSE)
{
  # This code plots the heatmaps transposed, caution, overwrites the files created above.
  for (ana in 1:nLabels) for (tg in 1:nTargets)
  {
    targetLabels = translateUsingTable(geneAnnotation.ana.x[[ana]]$Entrez, targetNets[[tg]]);
    ot = overlapTable(labels.x[[ana]]$data, targetLabels, ignore = 0)

    t1 = table(targetLabels);
    t1 = t1[match(colnames(ot$pTable), names(t1))];

    t2 = table(labels.x[[ana]]);
    t2 = t2[match(rownames(ot$pTable), names(t2))];

    plotOverlapHeatmap(ot$countTable, ot$pTable,
       colLabels = spaste("ME", labels2colors(as.numeric(colnames(ot$countTable)))),
       colSymbols = spaste(targetNames.short[tg], " M", colnames(ot$countTable), " (", t1, ")"),
       rowLabels = spaste("ME", labels2colors(as.numeric(rownames(ot$countTable)))),
       rowSymbols = spaste(sub("\\..*", "", labelNames[ana]), " M", rownames(ot$countTable), " (", t2, ")"),
       main = spaste("Overlap of all ", labelNames.pretty[ana], " with ", targetNames.pretty[tg]),
       mar = c(9, 11, 1.5, 0.2),
       plotFile = spaste(labelPlotDirs[ana], "/overlapWith.", names(targetNets)[tg], ".pdf"),
       separatorInterval = 3);

    plotOverlapHeatmap(ot$countTable, ot$pTable,
       colLabels = spaste("ME", labels2colors(as.numeric(colnames(ot$countTable)))),
       colSymbols = spaste(targetNames.short[tg], " M", colnames(ot$countTable), " (", t1, ")"),
       rowLabels = spaste("ME", labels2colors(as.numeric(rownames(ot$countTable)))),
       rowSymbols = spaste(sub("\\..*", "", labelNames[ana]), " M", rownames(ot$countTable), " (", t2, ")"),
       main = spaste("Overlap of selected ", labelNames.pretty[ana], " with ", targetNames.pretty[tg]),
       threshold = 1e-4,
       mar = c(9, 11, 1.5, 0.2),
       plotFile = spaste(labelPlotDirs[ana], "/overlapWith.", names(targetNets)[tg], "-significant.pdf"),
       separatorInterval = 3);
  }
}



