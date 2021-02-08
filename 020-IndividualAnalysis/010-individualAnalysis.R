# Individual analysis

source("../Functions/networkFunctions-extras-17.R")
source("../Functions/labelPoints2-01.R");
source("../Functions/outlierRemovalFunctions.R")
source("../Functions/preprocessing-General-013.R")
source("../Functions/GNVFunctions-015.R")
source("../Functions/individualAnalysis-General-007-02.R");

library(anRichment)

dir.create("RData", recursive = TRUE);
#dir.create("Results", recursive = TRUE);
dir.create("Plots", recursive = TRUE);

# Load data

x = load("../010-Preprocessing/RData/preprocessedData.RData");

multiCounts = multiExpr;

tissues = setNames
ltissues = tolower(tissues);

nSets = length(setNames);

#=============================================================================================================
#
# Check if the removed 1st PC/1st SVA factor correlates with any of the traits.
#
#=============================================================================================================
numericCols = mtd.apply(multiSampleAnnot, sapply, is.numeric);

multiPheno = mtd.subset(multiSampleAnnot, , numericCols[[1]]$data);

cp = corAndPvalue(multiPheno[[1]]$data);
diag(cp$p) = NA

dir.create("Plots", recursive = TRUE);
pdf(file = "Plots/correlationHeatmapOfTraitsAndSurrogateVariables.pdf", wi = 12, he = 10);
corHeatmapWithDendro(cp$cor, pValues = cp$p, main = "Correlations of phenotypes and surrogate variables",
      mar.main = c(10, 14, 2,1), cex.text = 0.6);

dev.off()


#=============================================================================================================
#
# Create designs and contrasts
#
#=============================================================================================================

dietOrder.paper = c("Normal", "HS", "HS-Chol2%", "HS-Chol2%-CA", "Chol2%-CA", "Cholate (CA)");
dietOrder = translateUsingTable(dietOrder.paper, as.data.frame(prettifyList[c(2,1)]));


## Note diet name mappings

#                       NC                  Normal
#                       HS                      HS
#            Cholate_CA_2m            Cholate (CA)
#          ChR_2percent_2m            HS-Chol2%-CA
#      Chol_2percent_HS_2m               HS-Chol2%
#          Chol2percent_CA               Chol2%-CA


#                       NC                      NC
#                       HS                      HS
#            Cholate_CA_2m         Cholate (CA) 2m
#          ChR_2percent_2m             ChR 2% (2m)
#      Chol_2percent_HS_2m         Chol 2%_HS (2m)
#          Chol2percent_CA               Chol2%-CA

# Compare all to all

nDiets = length(dietOrder)

numerators = denominators = character(0);

for (d in 1:(nDiets-1)) for (n in (d + 1):nDiets)
{
  denominators = c(denominators, dietOrder[d]);
  numerators = c(numerators, dietOrder[n]);
}

numerators.pretty = prettifyStrings(numerators, prettifyList);
denominators.pretty = prettifyStrings(denominators, prettifyList);

designInfo.diets = data.frame(
  design = spaste("~ Diet.machine + SVAFactor.1"),
  contrast.name = "Diet.machine",
  contrast.numerator = numerators,
  contrast.denominator = denominators,
  coefficientName = NA,
  testName = spaste(numerators, ".vs.", denominators),
  testName.pretty = spaste(numerators.pretty, " vs. ", denominators.pretty),
  combineSignificantHits = FALSE,
  splitSignificantHitsBy = NA);

designInfo = designInfo.diets;

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
  coll.GO = buildGOcollection(organism = "mouse");
  coll.MSigDB = MSigDBCollection(file = "../Data/Annotation/MSigDB/msigdb_v6.2.xml",
                  organism = "mouse", MSDBVersion = "6.2");
  coll.BioSystems = BioSystemsCollection(organism = "mouse");
  coll.genomicPosition = genomicPositionCollection(organism = "mouse", 
                  spacings = 5e6, namePrefix = "Genomic position", 
                  dataSource = "Entrez Gene, via Bioconductor annotation packages", 
                  overlapFactor = 2, membershipBy = "start", useUnits = "Mb", 
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
         outputOrganism = "mouse");  

  names(coll.Enrichr) = c("ChEA", "HistoneModif", "TFChIPSeq", "mirTarBase");

  stdCollections = c(list(GO = coll.GO, MSigDB = coll.MSigDB, BioSystems = coll.BioSystems, 
                          genomicPosition = coll.genomicPosition),
                     coll.Enrichr);

  save(stdCollections, file = "RData/stdCollections.RData");
}

sum(sapply(stdCollections, nDataSets))
# [1] 35232

#==================================================================================================
#
# prepare groups for collection of DE genes
#
#==================================================================================================

baseGroupNames = spaste("DE genes in ", ltissues, " of mice on various diets");
baseGroups = mymapply(function(name, ts)
   list(newGroup(name,
                  spaste("Differentially expressed genes between various diet combinations ",
                         "in wild-type mice."),
                  "Analysis by Peter Langfelder for CRO Bioinformatics",
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
    geneSetNamePattern = spaste("%d%p for %t in ", tolower(setNames), ", study of WT mice on various diets"),
    geneSetDescriptionPattern = spaste("Genes %d for %t at the significance threshold %p in ",
                              tolower(setNames),
                               ", study of WT mice on various diets."),
    geneSetGroups = baseGroupNames,
    collectionFile = spaste("RData/collectionOfDEGenes-", setNames, ".RData"),
    groupList = baseGroups,

    MoreArgs = list(
      addToAnalysisManager = FALSE,
      designInfo = designInfo,
      combineResults = list(all = designInfo$testName),
      combineNames = "all",
      combineColNameExtensions = "",
      geneSource = "Analysis by Peter Langfelder for CRO Bioinformatics.",
      resultDir = file.path("Results"),
      #organism = "mouse",
      #stopAt = "20-individual results",
      topThresholds = topThresholds,
      topList.minN = 0,
      enrichment.minN = 0,
      prettifyList = prettifyList,
      extendPrettifyList = TRUE,
      enrichment.collections = stdCollections,
      enrichment.compress = TRUE,
      keepFullDESeqObjects = FALSE,
      topList.machineColName = "Gene",
      topList.humanColName = "Gene",
      topList.backupColName = "Gene",
      forceExport = TRUE,
      organism = "mouse",
      internalClassification = "WT mice on various diets",
      indent = 2, verbose = 2),
    mdmaVerbose = TRUE);
  save(analyses, file = "RData/analyses.RData")
}

prettifyList2 = analyses[[1]]$data$prettifyList;

save(prettifyList2, file = "RData/prettifyList2.RData");

enrichment.dropCols = "dataSetID";

combCollection = do.call(mergeCollections, stdCollections);
analyses = mtd.mapply(function(ana, collections)
{
  printFlush(ana$analysisName)
  for (cri in 1:length(ana$enrichment))
  {
    enr.subColl = splitEnrichmentResultsBySubcollections(ana$enrichment[[cri]]$CombinedCollection, combCollection,
             collections, dropColumns = enrichment.dropCols);
    ana$enrichment[[cri]] = c(ana$enrichment[[cri]], enr.subColl);
  }
  ana;
}, analyses, MoreArgs = list(collections = stdCollections));


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
# Enrichment barplots
#
#==================================================================================================

existingCollectionNames = names(stdCollections);
plotCollections = c(existingCollectionNames);
existingNames2 = existingCollectionNames;
fileIDs = c(existingNames2);

collectionName.pretty = prettifyStrings(fileIDs,
   list(c("BioSystems", "GO", "MSigDB", 
          "genomicPosition", "ChEA", "HistoneModif", "TFChIPSeq", "miTarBase"),
        c("NCBI BioSystems pathays", "GO terms", "MSigDB sets", "genomic position sets", 
          "ChEA", "ENCODE histone modification sets", "ENCODE TF ChIP-seq sets", "mirTarBase miRNA targets")));
          

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






