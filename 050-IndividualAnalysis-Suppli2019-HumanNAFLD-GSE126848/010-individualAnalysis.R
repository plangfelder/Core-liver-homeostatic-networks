# Standard preprocessing

source("../Functions/networkFunctions-extras-18.R")
source("../Functions/labelPoints2-01.R");
source("../Functions/outlierRemovalFunctions.R")
source("../Functions/preprocessing-General-013.R")
source("../Functions/GNVFunctions-016.R")
source("../Functions/individualAnalysis-General-008.R");

library(anRichment)

dir.create("RData", recursive = TRUE);
#dir.create("Results", recursive = TRUE);
dir.create("Plots", recursive = TRUE);

# Load data

x = load("../040-Preprocessing-Suppli2019-HumanNAFLD-GSE126848/RData/preprocessedData.RData");

tissues = "Liver"
ltissues = tolower(tissues);

nSets = length(setNames);

#=============================================================================================================
#
# Create designs and contrasts
#
#=============================================================================================================


# Compare all to all

diseaseOrder = c("healthy", "obese", "NAFLD", "NASH");
nDiseases = length(diseaseOrder)

numerators = denominators = character(0);

for (d in 1:(nDiseases-1)) for (n in (d + 1):nDiseases)
{
  denominators = c(denominators, diseaseOrder[d]);
  numerators = c(numerators, diseaseOrder[n]);
}

numerators.pretty = prettifyStrings(numerators, prettifyList);
denominators.pretty = prettifyStrings(denominators, prettifyList);

designInfo.diseases = data.frame(
  design = spaste("~ disease + SVAFactor.1 + SVAFactor.2 + gender"),
  contrast.name = "disease",
  contrast.numerator = numerators,
  contrast.denominator = denominators,
  coefficientName = NA,
  testName = spaste(numerators, ".vs.", denominators),
  testName.pretty = spaste(numerators.pretty, " vs. ", denominators.pretty),
  combineSignificantHits = FALSE,
  splitSignificantHitsBy = NA);

designInfo = designInfo.diseases;

save(designInfo, file = "RData/designInfo.RData");


#=============================================================================================================
#
# Significance threhsolds for reporting genes
#
#=============================================================================================================

topThresholds = list(._FDR.lessThan.0.05_. = list(pAdj = 0.05, p = 1));

nTopThresholds = length(topThresholds);

prettifyList[[1]] = c(prettifyList[[1]], 
            "Differentially.expressed", ".and.", "._FDR.lessThan.", "._p.lessThan.", "_.", "inter.", ".with.");
prettifyList[[2]] = c(prettifyList[[2]],
            "Differentially expressed", " and ", " (FDR < ", " (p < ", ")", "interaction of ", " with ");


#=============================================================================================================
#
# Collections for enrichment analysis
#
#=============================================================================================================

if (!checkRDataFile("RData/stdCollections.RData"))
{
  coll.GO = buildGOcollection(organism = "human");
  coll.MSigDB = MSigDBCollection(file = "../Data/Annotation/MSigDB/msigdb_v6.2.xml",
                  organism = "human", MSDBVersion = "6.2");
  coll.BioSystems = BioSystemsCollection(organism = "human");
  coll.genomicPosition = genomicPositionCollection(organism = "human", 
                  spacings = 5e6, namePrefix = "Genomic position", 
                  dataSource = "Entrez Gene, via Bioconductor annotation packages", 
       individualSignificance.ana[[1]]$data[, 1]           overlapFactor = 2, membershipBy = "start", useUnits = "Mb", 
                  unit = 1e+06);

  enrichrMeta = data.frame(file = c("ChEA_2016.bz2", "ENCODE_Histone_Modifications_2015.bz2", "ENCODE_TF_ChIP-seq_2015.txt.bz2",
                                     "miRTarBase_2017.bz2"),
     IDBase = spaste("Enrichr.", c("ChEA", "HistoneModif", "TF-ChIP-seq", "mirTarBase"), "."),
     nameBase = c("ChEA ", "Histone Modification ", "TF-ChIP-seq ", "mirTarBase "),
     descriptionBase = spaste("Enrichr ", c("ChEA ", "ENCODE Histone Modification ", "ENCODE TF-ChIP-seq ", "mirTarBase "), 
                              " library"),
     enrichrClass = c("ChEA", "ENCODE Histone Modification", "ENCODE TF-ChIP-seq", "mirTarBase"));

  coll.Enrichr = mymapply(genericEnrichrCollection, 
         enrichrFile = file.path("../Data/Annotation/Enrichr", enrichrMeta$file),
         IDBase = enrichrMeta$IDBase, nameBase = enrichrMeta$nameBase,
         descriptionBase = enrichrMeta$descriptionBase,
         enrichrClass = enrichrMeta$enrichrClass,
         outputOrganism = "human");  

  names(coll.Enrichr) = c("ChEA", "HistoneModif", "TFChIPSeq", "mirTarBase");

  stdCollections = c(list(GO = coll.GO, MSigDB = coll.MSigDB, BioSystems = coll.BioSystems, 
                          genomicPosition = coll.genomicPosition),
                     coll.Enrichr);

  save(stdCollections, file = "RData/stdCollections.RData");
}

sum(sapply(stdCollections, nDataSets))
# [1] 36652

#==================================================================================================
#
# Add the collection from the mouse analysis
#
#==================================================================================================

mouseColl = loadAsList("../020-IndividualAnalysis/RData/collectionOfDEGenes-Liver.RData")[[1]];

allCollections = c(stdCollections, list(mouseLiver = convertCollectionToOrganism(mouseColl, organism = "human")))

#==================================================================================================
#
# prepare groups for collection of DE genes
#
#==================================================================================================

baseGroupNames = spaste("DE genes in human ", ltissues, " between normal, obese, NALFD and NASH samples"); 
baseGroups = mymapply(function(name, ts)
   list(newGroup(name,
                  spaste("Differentially expressed genes (FDR < 0.05) between various disease status pairs ",
                         "in human liver. "),
                  "Analysis by Peter Langfelder for Bioinformatics CRO",
                  parents = c(ts, "mRNA differential expression"))),
  baseGroupNames, tissues);

#==================================================================================================
#
# Individual analysis
#
#==================================================================================================

dir.create("RData", recursive = TRUE);
if (!checkRDataFile("RData/analyses.RData"))
{
  analyses = listRep(list(data = list()), nSets)

  analyses = mtd.mapply(
    individualAssociation.general,
    data = multiCounts,
    dataWeights = multiWeights,
    pheno = multiSampleAnnot,
    intermediateResults = analyses,

    analysisName = setNames,
    geneAnnotation.start = multiGeneAnnot,

    IDBase = setNames,
    # character string with tags %t and %d, %p for test, direction and significance
    geneSetNamePattern = spaste("%d%p for %t in human ", tolower(setNames)),
    geneSetDescriptionPattern = spaste("Genes %d for %t at the significance threshold %p in human ",
                              tolower(setNames), "."),
    geneSetGroups = baseGroupNames,
    collectionFile = spaste("RData/collectionOfDEGenes-", setNames, ".RData"),
    groupList = baseGroups,

    MoreArgs = list(
      addToAnalysisManager = FALSE,
      designInfo = designInfo,
      combineResults = list(all = designInfo$testName),
      combineNames = "all",
      combineColNameExtensions = "",
      geneSource = "Analysis by Peter Langfelder for Bioinformatics CRO.",
      resultDir = file.path("Results"),
      topThresholds = topThresholds,
      topList.minN = 0,
      enrichment.minN = 0,
      prettifyList = prettifyList,
      extendPrettifyList = TRUE,
      enrichment.collections = allCollections,
      enrichment.compress = TRUE,
      keepFullDESeqObjects = FALSE,
      topList.machineColName = "key",
      topList.humanColName = "Symbol",
      topList.backupColName = "key",
      forceExport = TRUE,
      organism = "human",
      internalClassification = "Human liver (normal, obese, NAFLD, NASH)",
      indent = 2, verbose = 2),
    mdmaVerbose = TRUE);
  save(analyses, file = "RData/analyses.RData")
}

prettifyList2 = analyses[[1]]$data$prettifyList;

save(prettifyList2, file = "RData/prettifyList2.RData");

enrichment.dropCols = "dataSetID";

combCollection = do.call(mergeCollections, allCollections);
analyses = mtd.mapply(function(ana, collections)
{
  printFlush(ana$analysisName)

  for (cri in 1:length(ana$enrichment)) for (thr in 1:length(ana$enrichment[[cri]]))
  {
    enr.subColl = splitEnrichmentResultsBySubcollections(ana$enrichment[[cri]][[thr]]$CombinedCollection, combCollection,
             collections, dropColumns = enrichment.dropCols);
    ana$enrichment[[cri]][[thr]] = c(ana$enrichment[[cri]][[thr]], enr.subColl);
  }
  ana;

}, analyses, MoreArgs = list(collections = allCollections));


#==================================================================================================
#
# Plot a christmas tree of numbers of significant genes
#
#==================================================================================================

nSignif = mtd.apply(analyses, getElement, "nSignif.all")

col.right = "#FF9070"; col.left = "skyblue";
border.right = "brown"; border.left = "blue";
plotDir = "Plots";
dir.create(plotDir, recursive = TRUE);

for (set in 1:nSets)
{
  ns1 = nSignif[[set]]$data[[1]];
  labels = prettifyStrings(rownames(ns1), prettifyList2);
  labWidth = marginWidth(labels);
  pdf(file = spaste(plotDir, "/numbersOfSignifGenes-", setNames[set], "-allTests.pdf"), 
                    wi = 4 + labWidth$inch, he = 0.5 + 0.25 * nrow(ns1));
  par(mar = c(0.2, labWidth$lines, 2, 0.2));
  for (thr in 1:length(topThresholds))
  {
    ns1 = nSignif[[set]]$data[[thr]];
    twoSideBarplot(ns1[, 1], ns1[, 2],
       col.left = col.left, col.right = col.right,
       border.left = border.left, border.right = border.right,
       yLabels = labels, barGap = 0.2, 
       main = spaste("Number of significant genes ",  prettifyStrings(names(topThresholds)[thr], prettifyList),
                     "   "), cex.main = 1,
       separatorPositions = c(5, 9, 12), sep.ext = TRUE, sep.col = "grey80");
  }

  dev.off();
}


#==================================================================================================
#
# Scatterplots of Z statistics and overlap heatmaps
#
#==================================================================================================

results = mtd.apply(analyses, function(ana) ana$combinedDESeqResults$all);
results = mtd.apply(results, prettifyNames, prettifyList2)

Zs = mtd.subset(results, , grep("Z for ", mtd.colnames(results)))

setNames.pretty = "human NAFLD liver";

pdf(file = "Plots/correlationHeatmap-Zs.pdf", wi = 7, he = 4);
for (set in 1:nSets)
  corHeatmapWithDendro(cor(Zs[[set]]$data),
     main = spaste("Correlation of DE Z statistics, ", setNames.pretty[set], "                      "),
     xLabels = sub("Z for ", "", colnames(Zs[[set]]$data)),
     yLabels = sub("Z for ", "", colnames(Zs[[set]]$data)),
     mar.main = c(6, 8, 2, 1));

dev.off();


# Load mouse results

mouseResults0 = read.csv(
      gzfile("../020-IndividualAnalysis/Results/01-Association/associationWithPhenotypes-Liver-all.csv.gz"),
      check.names = FALSE);

mouseEntrez.human = mapEntrez(mouseResults0$Entrez, orgFrom = "mouse", orgTo = "human");

fin = !is.na(mouseEntrez.human);
table(fin)
#FALSE  TRUE 
# 4586 13485 

mouseResults1 = data.frame(mouseResults0[ fin, ], humanEntrez = mouseEntrez.human[fin], check.names = FALSE);

mouseResults = mouseResults1[ match(results[[1]]$data$Entrez, mouseResults1$humanEntrez), ];

table(is.na(mouseResults$Entrez))
#FALSE  TRUE 
#12811  2680 

Z.mouse = mouseResults[, grep("Z for ", names(mouseResults))]

pdf(file = "Plots/correlationHeatmaps-Zs-humanVsMouse.pdf", wi = 8, he = 6);
for (set in 1:nSets)
{
  cor1 = cor(Z.mouse, Zs[[set]]$data, use = 'p');
  colTree = bestTree(1-cor(Zs[[set]]$data));
  rowTree = bestTree(1-cor(Z.mouse, use = 'p'));

  heatmapWithDendro(cor1, rowTree = rowTree, colTree = colTree,
          mar.main = c(6, 13, 0, 1),
          xLabels = sub("Z for ", "", colnames(Zs[[set]]$data)),
          yLabels = sub("Z for ", "", colnames(Z.mouse)),
          main = spaste("Correlations of DE statistics in human and mouse data                                  "))
}

dev.off()

par(mfrow = c(1,1))
annotatedScatterplot(Zs[[set]]$data$`Z for NAFLD vs. obese`, 
     Z.mouse$`Z for HS-Chol2%-CA vs. Normal`,
     pointLabels = results[[set]]$data$Symbol, 
     addOverlapPValues = TRUE);


#==================================================================================================
#
# Enrichment barplots
#
#==================================================================================================

existingCollectionNames = names(allCollections);
plotCollections = c(existingCollectionNames);
existingNames2 = existingCollectionNames;
fileIDs = c(existingNames2);

collectionName.pretty = prettifyStrings(fileIDs,
   list(c("BioSystems", "GO", "MSigDB", 
          "genomicPosition", "ChEA", "HistoneModif", "TFChIPSeq", "miTarBase", "mouseLiver"),
        c("NCBI BioSystems pathays", "GO terms", "MSigDB sets", "genomic position sets", 
          "ChEA", "ENCODE histone modification sets", "ENCODE TF ChIP-seq sets", "mirTarBase miRNA targets",
          "Mouse liver")));
          

enrichmentThreshold = 0.05;

keepClasses = list(DEGenes = c("Downregulated.for.[^i]", "Upregulated.for.[^i]"));

plotDir.enrichBarp = file.path(plotDir, "EnrichmentBarplots", setNames);
lapply(unique(plotDir.enrichBarp), dir.create, recursive = TRUE);

library(Cairo)
useSets = 1:nSets
for (th in 1:nTopThresholds) for (cl in 1:length(keepClasses))
{
  clName = names(keepClasses)[cl];
  thName = names(topThresholds)[th];
  dropTh = names(topThresholds)[-th];
  for (pcoll in 1:length(plotCollections))
    mtd.mapply( enrichmentBarplot.standardAnalysis,
     analyses[useSets],
     analysisName.pretty = setNames[useSets],
     main = spaste("Enrichment of top genes ",
                   "in ", collectionName.pretty[pcoll]),
     plotDir = plotDir.enrichBarp,
    MoreArgs = list(
       #plotDev = NULL,
       plotFileID = spaste(clName, "-", thName, "-", prependZeros(pcoll, 2), "-", fileIDs[pcoll]),
       collectionName = plotCollections[pcoll],
       plotFileBase = "enrichment",
       width = 12,
       height = NULL,
       heightPerRow = 0.23,
       maxRowsPerPage = 50,
       keepClassOnOnePage = TRUE,
       mar = c(3.3, 30, 2, 0.8),
       plotComponent = "all",
       prettifyList = prettifyList2,
       dropModuleEnrichmentLabels = TRUE,
       frame = FALSE,
       rescueOnly = FALSE,
       classBesideBars = FALSE,
       splitClassOnAnd = FALSE,
       enrichmentThreshold = enrichmentThreshold,
       thresholdColumn = "Bonferroni",
       maxRank = 10,
       separatorSpace = 2,
       dropClassPattern = c(dropTh),
       keepClassPattern = keepClasses[[cl]],
       maxLines.genes = 1, maxLines.geneSets = 1,
       printGenes = TRUE,
       cex.classLabels = 0.8, cex.genes = 0.8, cex.terms = 1,
       pageNumberEnd = ")                                           "),

    mdmaVerbose = TRUE
  )
}

# One special plot

classPattern.special = list( NASHvsNAFLD = c("NASH.vs.NAFLD"));

for (pcoll in match("mouseLiver", plotCollections)) for (cps in 1:length(classPattern.special))
 mtd.mapply( enrichmentBarplot.standardAnalysis,
     analyses[useSets],
     analysisName.pretty = setNames[useSets],
     main = spaste("Enrichment of top genes ",
                   "in ", collectionName.pretty[pcoll]),
     plotDir = plotDir.enrichBarp,
    MoreArgs = list(
       #plotDev = NULL,
       plotFileID = spaste(names(classPattern.special)[cps], "-", thName, "-", prependZeros(pcoll, 2), "-", fileIDs[pcoll]),
       collectionName = plotCollections[pcoll],
       plotFileBase = "enrichment",
       width = 12,
       height = NULL,
       heightPerRow = 0.23,
       maxRowsPerPage = 50,
       keepClassOnOnePage = TRUE,
       mar = c(3.3, 30, 2, 0.8),
       plotComponent = "all",
       prettifyList = prettifyList2,
       dropModuleEnrichmentLabels = TRUE,
       frame = FALSE,
       rescueOnly = FALSE,
       classBesideBars = FALSE,
       splitClassOnAnd = FALSE,
       enrichmentThreshold = enrichmentThreshold,
       thresholdColumn = "Bonferroni",
       maxRank = 10,
       separatorSpace = 2,
       dropClassPattern = c(dropTh),
       keepClassPattern = classPattern.special[[cps]],
       maxLines.genes = 1, maxLines.geneSets = 1,
       printGenes = TRUE,
       cex.classLabels = 0.8, cex.genes = 0.8, cex.terms = 1,
       pageNumberEnd = ")                                           "))



#==================================================================================================
#
# Enrichment matrix plot for human vs. mouse DE genes
#
#==================================================================================================


for (set in 1:nSets)
{
  enr1 = analyses[[set]]$data$enrichment$all$._FDR.lessThan.0.05_.$mouseLiver
  n1 = enr1$countsInDataSet;
  p1 = enr1$pValues;
  setNames1 = multiSubr(c("regulated", " (FDR < 0.05)", " in liver, study of WT mice on various diets"), 
                        as.character(dataSetNames(allCollections[["mouseLiver"]])), fixed = TRUE);
  deNames1 = multiSubr(c("regulated", " (FDR < 0.05)"), prettifyStrings(colnames(n1), prettifyList2), fixed = TRUE);

  colTree = bestTree(1-cor(Zs[[set]]$data));
  rowTree = bestTree(1-cor(Z.mouse, use = 'p'));

  colOrder = order(substring(deNames1, 1, 2))[c(colTree$order, colTree$order + ncol(Zs[[set]]$data))];
  rowOrder = order(substring(setNames1, 1, 2))[c(rowTree$order, rowTree$order + ncol(Z.mouse))];

  plotOverlapHeatmap(n1[rowOrder, colOrder], p1[rowOrder, colOrder],
       rowLabels = setNames1[rowOrder],
       colLabels = deNames1[colOrder],
       plotFile = spaste("Plots/enrichmentInMouseDEGenes-", setNames[set], ".pdf"),
       baseWidth = 4.5, baseHeight = 3,
       mar = c(8.5, 17, 2, 1),
       separatorInterval = 3,
       main = spaste("Enrichment of DE genes in ", setNames.pretty[set], " in mouse DE genes                     "));
}
       

