source("../Functions/networkFunctions-extras-18.R")
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

x = load("../040-Preprocessing-Suppli2019-HumanNAFLD-GSE126848/RData/preprocessedData.RData");

prettifyList = loadAsList("../050-IndividualAnalysis-Suppli2019-HumanNAFLD-GSE126848/RData/prettifyList2.RData")[[1]]

tissues = "Liver"
ltissues = tolower(tissues);

nSets = length(setNames);

library(Cairo);
library(doParallel)

expr.flat = multiExpr.VST;
sampleAnnot.flat = multiSampleAnnot;
weights.flat = multiWeights;

setNames.flat = setNames;
nSets.flat = length(expr.flat);
names(expr.flat) = names(weights.flat) = names(sampleAnnot.flat) = setNames

cor = WGCNA::cor;

setsByTissue = tapply(setNames, tissues, identity)

setNames.pretty.flat = "Human NAFLD liver";
lsetNames.pretty.flat = "human NAFLD liver"

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

# Will use power 6

softPowers = rep(6, nSets.flat);

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

analysisTrees = list(
  humanLiver = newConsensusTree(
      consensusOptions = consensusOptions,
      inputs = setNames.flat, 
      analysisName = "humanLiver")
);

analysisNames = sapply(analysisTrees, getElement, "analysisName");

analysisNames.pretty = c("human liver");
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
   ga[ match(mtd.colnames(.multiExpr), ga$key), ];
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
individualAssoc = lapply(setNames.flat, function(name)
  read.csv(gzfile(spaste("../050-IndividualAnalysis-Suppli2019-HumanNAFLD-GSE126848/Results/01-Association/",
                         "associationWithPhenotypes-", name, "-all.csv.gz")), check.names = FALSE));

names(individualAssoc) = setNames;

# For each analysis, put together the Q association from the relevant tissues

individualResults.ana = mymapply(function(.expr, .sets)
{
  IDs = mtd.colnames(.expr);
  do.call(cbind, removeListNames(lapply(individualAssoc[.sets], function(ia)
     ia[ match(IDs, ia$key), grep("Z for ", names(ia))])))
}, expr.ana, analysisSets);

individualSignificance.ana = mymapply(function(.expr, .sets)
{
  IDs = mtd.colnames(.expr);
  out = do.call(cbind, removeListNames(lapply(individualAssoc[.sets], function(ia)
     ia[ match(IDs, ia$key), grep("FDR for ", names(ia))])))
  0 + (out < 0.05)
}, expr.ana, analysisSets);

individualColors.ana = lapply(individualResults.ana, numbers2colors, commonLim = TRUE, lim = c(-10, 10));

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

kmeMouse = read.csv(gzfile("../030-NetworkAnalysis/Results/networkAnalysisResults-Liver.main.csv.gz"));

mouseLabels0 = data.frame(Entrez = mapEntrez(kmeMouse$Entrez, orgFrom = "mouse", orgTo = "human"),
                          Label = kmeMouse$module);

mouseLabels = mouseLabels0[ !is.na(mouseLabels0$Entrez) & mouseLabels0$Label!=0, ];


#===============================================================================================================
#
# Module analysis
#
#===============================================================================================================

mods.ns = list();
nSets.ana = sapply(analysisSets, length)
detectCutHeight = 0.998;
deepSplit = 3;
for (ana in 1:nAnalyses)
 mods.ns[[ana]] = hierarchicalConsensusModules(
     consensusTOMInfo = consensusTOMInfo[[ana]],
     multiExpr = expr.ana[[ana]],
     multiWeights = weights.ana[[ana]],
     checkMissingData = FALSE,
     detectCutHeight = detectCutHeight,
     deepSplit = deepSplit,
     minModuleSize = 20,
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


labels.ns = lapply(labels.original.ns, matchLabelsToReference.2, mouseLabels);

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

dir.create("../059-StabilityAnalysis-Suppli2019-HumanNAFLD-GSE126848/RData", recursive = TRUE);

save(expr.ana, weights.ana, analysisTrees, networkOptions.ana, detectCutHeight, deepSplit,
     file = "../059-StabilityAnalysis-Suppli2019-HumanNAFLD-GSE126848/RData/dataForStabilityAnalysis.RData");


#========================================================================================================
#
# Re-calculate modules with the input of sampled labels.
#
#========================================================================================================

nRuns = 50;
sampledLabels.lst = lapply(1:nAnalyses, function(ana)
  lapply(1:nRuns, function(run)
  {
     out = try(loadAsList(spaste("../059-StabilityAnalysis-Suppli2019-HumanNAFLD-GSE126848/RData/sampledLabels-run",
               prependZeros((ana-1)*nRuns + run-1, 3),
                ".RData"))[[1]]);
     if (inherits(out, "try-error")) NULL else out;
  })
);

sampledLabels = lapply(sampledLabels.lst, function(x) do.call(cbind, x))

lapply(sampledLabels, dim)
gc()
mods.lst = list();

minStabilityDissim = list(
    main = c(0.53),
    sub = c(0.3));

for (size in names(minStabilityDissim))
{
  mods.lst[[size]] = list();
  for (ana in 1:nAnalyses)
    mods.lst[[size]] [[ana]] = hierarchicalConsensusModules(
        consensusTOMInfo = consensusTOMInfo[[ana]],
        multiExpr = expr.ana[[ana]],
        multiWeights = weights.ana[[ana]],
        checkMissingData = FALSE,
        detectCutHeight = detectCutHeight,
        deepSplit = deepSplit,
        minModuleSize = 20,
        useBranchEigennodeDissim = TRUE,
        stabilityLabels = sampledLabels[[ana]],
        minStabilityDissim = minStabilityDissim[[size]] [ana],
        mergeCutHeight = 0.2,
        minCoreKME = 0.4, minKMEtoStay = 0.2,
        verbose = 5)
}

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

labels = lapply(labels.original, matchLabelsToReference.2, mouseLabels);

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
     file = "RData/dataForModulePreservation-humanLiver.RData");
    
#========================================================================================================
#
# First look at dendrograms
#
#========================================================================================================

plotColors = mtd.mapply(cbind, subColors, colors, individualColors.ana)
plotTrees = lapply(mods, function(md) clipTree(md$dendrograms[[1]], p = 0.01));
nPlotPheno = lapply(plotTraitNames, length)

pdf(file = "Plots/geneDendrogramsAndGS.pdf", wi = 14, he = 9);
#sizeGrWindow(14, 9);
par(lheight = 0.85)
for (ana in 1:nAnalyses)
plotDendroAndColors(plotTrees[[ana]],
                    plotColors[[ana]]$data,
                    c("Sub-modules", "Modules", plotTraitNames[[ana]]),
                    autoColorHeight = FALSE,
                    colorHeight = 0.5,
                    rowText = labels[[ana]]$data, 
                    textPositions = 2,
                    cex.rowText = 0.8,
                    rowWidths = c(1, 1, 4, rep(1, nPlotPheno[[ana]])),
                    rowTextAlignment = "center",
                    dendroLabels = FALSE,
                    main = spaste("Gene clustering, modules and association with disease, ", analysisNames.pretty[ana]),
                    marAll = c(0, 9, 1.5, 0), hang = 0.02, guideHang = 0.08,
                    cex.colorLabels = 0.8,
                    addGuide = TRUE);

dev.off();

# Plot gene dendrograms and indicator of significance (FDR<0.1)
plotColors2 = mtd.mapply(cbind, subColors, colors, individualSignificance.ana)

pdf(file = "Plots/geneDendrogramsAndGS-significanceIndicators.pdf", wi = 14, he = 9);
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
                    marAll = c(0, 9, 1.5, 0), hang = 0.02, guideHang = 0.08,
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
  focusOnGroups = c("KEGG", "GO.BP", "GO.MF", "GO.CC", "REACTOME", "Genomic position", "MSigDB", 
                    "ChEA", "ENCODE Histone Modification", "mirTarBase", "all"),
  groupShortNames = c("KEGG", "GO.BP", "GO.MF", "GO.CC", "Reactome", "GenPos", "MSigDB", 
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

conKME = mtd.apply(conKME, function(.conKME) setNames(.conKME, multiSub(c("[Cc]onsKME", ".HumanNAFLDLiver", "consensus."), 
                                                                        c("kME", "", ""),
                                    names(.conKME))));

for (ana in 1:nLabels)
{
  write.csv(signifNumeric(conKME[[ana]]$data, 3),
          file = gzfile(spaste("Results/networkAnalysisResults-", labelNames[ana], ".csv.gz")),
          row.names = FALSE, quote = TRUE);
}



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

designInfo0 = loadAsList("../050-IndividualAnalysis-Suppli2019-HumanNAFLD-GSE126848/RData/designInfo.RData")[[1]];

di = designInfo0;
  di$design = mapply(function(ds, ts) multiSub(c("gender", "disease"), c("Male.vs.Female", ts), ds), di$design, di$testName);
  di$contrast.name = NA;
  di$contrast.numerator = NA;
  di$contrast.denominator = NA;
  di$coefficientName = di$testName;
designInfo = di;

ss.ME = mtd.mapply(function(.MEs, .preds, .meInfo)
{
  ss = prettifyNames(associationScreening(
          response = .MEs[[1]]$data,
          predictors = .preds[[1]]$data,
          full = designInfo$design,
          termsOfInterest = designInfo$coefficientName,
          reduced = NULL,
          testName = designInfo$testName,
          additionalResponseInfo = .meInfo,
          addID = TRUE, getFDR = TRUE,
          FDRFromCombinedPValues = TRUE,
          getCLS = FALSE), prettifyList);
  # Add test of overall dependence with  F statistic
  ss1 = lm (as.matrix(.MEs[[1]]$data)~disease + SVAFactor.1 + SVAFactor.2 + gender, data = .preds[[1]]$data);
  sum = summary(ss1);
  fdof = do.call(rbind, lapply(sum, getElement, "fstatistic"));
  colnames(fdof) = c("F for disease", "DoF1 for disease", "DoF2 for disease");
  p = pf(fdof[,1], fdof[,2], fdof[,3], lower.tail = FALSE);
  pb = p*length(p);
  pb[pb > 1] = 1;
  fdr = p.adjust(p, method = "fdr");
  Z = ZfromP(p, p)
  cbind(ss, fdof, `p for disease` = p, `FDR for disease` = fdr, `Z for disease` = Z);
}, MEs, sampleAnnot.ana.x, meInfo, mdmaVerbose = TRUE);


for (ana in 1:nLabels)
  write.csv(signifNumeric(ss.ME[[ana]]$data, 3),
            file = gzfile(spaste("Results/moduleSummaryAndAssociationWithDisease-", labelNames[ana], ".csv.gz")),
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
  wi = 10.5; he = 3 + 0.33 * nrow(mat);
  pdf(file = spaste(labelPlotDirs[ana], "/Association-MEs-vs-disease-", labelNames[ana], ".pdf"), 
      wi = wi, he = he)
  #sizeGrWindow(wi, he);
  par(mar = c(6, 30, 2.8, 1));
  par(lheight = 0.85);
  labeledHeatmap.multiPage(mat.plot,
             xLabels = xlabs,
             yLabels = spaste("ME", labels2colors(
                translateUsingTable(as.numeric(eLabels[[ana]]$data$module), labelTranslationForColors[[ana]]$data))),
             ySymbols = formatLabels(eLabels[[ana]]$data$enrichmentLabel, maxCharPerLine = 110, maxLines = 3,
                                     split = "[ _]", fixed = FALSE, newsplit=" ", keepSplitAtEOL = FALSE),
             textMatrix = textMat,
             cex.text = 0.8,
             colors = blueWhiteRed(100, gamma = 1),
             zlim = c( -max, max),
             setStdMargins = FALSE,
             main = spaste("Module association with disease contrasts\n", labelNames.pretty[ana]),
             cex.lab.y = 0.7, cex.lab.x = 0.8, maxRowsPerPage = 200, maxColsPerPage = 25);

  dev.off();
}

#================================================================================================
#
# Show boxplots of module eigengenes vs. disease
#
#================================================================================================

diseaseOrder = c("healthy", "obese", "NAFLD", "NASH");

cex.main = 0.8;
for (ana in 1:nLabels)
{
  pdf(file = spaste(labelPlotDirs[ana], "/eigengenesVsDisease-", labelNames[ana], ".pdf"),
      wi = 7, he = 5);

  mes1 = MEs[[ana]]$data [[1]]$data;
  nm = ncol(mes1);

  for (m in 1:nm)
  {
    mm = eLabels[[ana]]$data$module[m];
    x = MEs[[ana]]$data [[1]] $ data [[m]];
    dt1 = sampleAnnot.ana.x[[ana]] [[1]] $ data $ disease;
    fdt1 = factor(dt1, levels = diseaseOrder);
    xl = tapply(x, fdt1, identity);
    par(mar = c(4, 3.1, 5, 0.5));
    par(mgp = c(2, 0.7, 0));
    labeledBoxplot(xl, g = NULL, names = diseaseOrder, xlab = "", ylab = spaste("ME", mm), main = "",
                   addScatterplot = TRUE, notch = TRUE, pars = list(boxwex = 0.5, staplewex = 0.5, outwex = 0.5))
    box = par("usr");
    width = box[2] - box[1];
    main = spaste(
      formatLabels(eLabels[[ana]]$data$enrichmentLabel[m], maxWidth = 0.92 * width, maxLines = 4, cex = cex.main, font = 2), 
        "\n",
      formatLabels(moduleHubs.ps[[ana]]$data[[m]], maxWidth = 0.92 * width, maxLines = 1, cex = cex.main, font = 2))

    title(main, cex.main = cex.main);
  }

  dev.off();
}

#=====================================================================================================================
#
# Scatterplots of module eigengenes
#
#=====================================================================================================================

for (ana in 1:nLabels)
{

  col = as.numeric(factor(sampleAnnot.ana.x[[ana]] [[1]] $ data $ disease, levels = diseaseOrder));
  pch = 21 + (col-1)%%5
  nMEs = ncol(MEs[[ana]]$data[[1]]$data)
  pdf(file = spaste(labelPlotDirs[ana], "/eigengeneScatterplots.pdf"), wi = 5, he = 5);
  scpp(2.6);
  for (m1 in 1:(nMEs-1)) for (m2 in (m1+1):nMEs)
  {
     verboseScatterplot(MEs[[ana]]$data[[1]]$data[, m1], MEs[[ana]]$data[[1]]$data[, m2],
            pch = pch, bg = col, col = "grey", cex = 2,
            cex.lab = 1, cex.axis = 1, cex.main = 1.2,
            main = spaste(colnames(MEs[[ana]]$data[[1]]$data)[m2], " vs. ", colnames(MEs[[ana]]$data[[1]]$data)[m1], "\n"),
            xlab = colnames(MEs[[ana]]$data[[1]]$data)[m1],
            ylab = colnames(MEs[[ana]]$data[[1]]$data)[m2]);
     legendClean("auto", points.x = MEs[[ana]]$data[[1]]$data[, m1], points.y = MEs[[ana]]$data[[1]]$data[, m2],
         legend = diseaseOrder, pch = 21 + (c(1:length(diseaseOrder))-1)%%5,
          pt.bg = 1:length(diseaseOrder), col = "grey", cex = 1, pt.cex = 1.7);
  }
  dev.off()
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

for (ana in 1:nLabels)
{
  dirs = spaste("CytoscapeFiles-", labelNames[ana], "/", c("TopConnections", "WholeModules", "AllModules-TopConnections"));
  lapply(dirs, dir.create, recursive = TRUE);

  ana.TOM = (ana-1) %% nAnalyses + 1;
  TOM = as.matrix(BD.getData(consensusTOMInfo[[ana.TOM]]$consensusData));
  gc()
  indivAssoc1 = individualAssoc[[analysisSets[[ ana.TOM]] ]]
  nodeAttr = data.frame(conKME[[ana]]$data[, c("key", "Entrez", "Symbol", "EnsChromosome", "EnsStart", 
                                               "module", "kME.inOwnModule")],
                     indivAssoc1[match(conKME[[ana]]$data$key , indivAssoc1$key),
                           multiGrep(c("meanExpr for", "Z for", "log2FC for", "p for", "FDR for"), names(indivAssoc1))],
                     check.names = FALSE);

  # Export top connections
  modules = setdiff(sort(unique(labels.x[[ana]]$data)), 0)
  topGenes = character(0);
  stopifnot(all.equal(conKME[[ana]]$data$key, geneAnnotation.ana.x[[ana]]$key))
  stopifnot(all.equal(conKME[[ana]]$data$key, mtd.colnames(expr.ana.x[[ana]])))
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

    
