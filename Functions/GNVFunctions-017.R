# This file contains functions specific to data sets that are generated within the GNV project

# Version 002: 
#   . order statistically fully rescued genes by rescue score as well. This means that rescue gene
#     lists with "top N genes" rather than genes actually passing the criterion will be the same as the
#     corresponding gene lists for "significant rescue". 
#   . Important: order of genotype level has changed; perturbed background is now highest. This switches the
#     interpretation of the sign of rescue scatterplot and makes many other plots a bit more logically 
#     consistent.

# Version 003: naming convention for combineNames changed: the name is not assumed to contain the analysis
# name. This affects enrichmentPlot.matrix and possibly other functions.

# Version 004: GNV.createDesignsAndContrasts also creates an (NA-filled) column coefficientNames.

# Version 005: Rescue manhattan plot and rescue colors/status now have only 3 values: no significant change,
# exacerbation and rescue. Statistically full rescue has been removed. 
# Plots of numbers of rescued plot by default only numbers of rescued and exacerbated genes, and the default
# limit has been raised to 10.

# Version 006: Changing input to addRescueColumns. The input can be general data (not just counts) but only one
# matrix (instead of a list of which only first component was used); plus options for exponentiating
# log-transformed data.

# Version 007: Adding the CHDI/Rancho collection to standard collections.

# Version 007-01: to avoid conflicts with commonFunctions-extras, legendForNumberPlots is renamed to
# legenedForNumberPlots.fromInfo

# Version 008: New standard collections, based on a function in networkFunctions.

# Version 009: adapting to new names of gene lists in output of individualAnalysis. This version is
# compatible with individualAnalysis-General-006.R

# 009-01: Re-tuning colors of "no significant move" points in GNV.rescuePlotAttributes
#  . Adding p-value histograms to GNV.standardPlots

# 010: Changing input to addRescueColumns. Instead of using baseline change from the same data set, use an
# externally supplied fold change and Z statistic. The output will not contain statististically full rescue
# since we don't use it anymore. 

# 011: Changing default AM analysis name in GNV.newPreprocessingInfo.

# 012: Adding extra output to GNV.rescueStatistics

# 012-01:
#   Chaning order of genotype and covariates in GNV.createDesignsAndContrasts

# 013: GNV.newPreprocessingInfo now distinguishes input and output base directories; default is for output to be
# /CommonReanalysis/ under GNVProject

# 014: Replacing aging collection for anRichment: instead of Allelic Series controls add the newer aging analysis of 2, 6,
# 10, 14, 18 month striatum. For now just the DE analysis, the network analysis is not added since at this time the
# collection has not been generated yet.

# 015: Adapting to networkFunctions-extras-16, specifically function marginWidth now has output with different names

# 016: by default, GNV.addAgingWGCNAcollection will also add the collection of aging WGCNA modules;
#   . the argument removeInteractionSets in GNV.addAgingWGCNAcollection is now TRUE by default

# 016-02: certain functions will check for presence of combinedDESeqResults.ext in analysis and use that if present since
# this data frame will contain columns from extended gene annotation
#   . Adapting to new gene sets in cell type markers, especially Gokce and Saunders

# 17: changing the meaning of "all" in GNV.createDesignsAndContrasts to keep the first 4 rather than first 3 designs.
#     . DIFF: the meaning of highKME in GNV.addCircadianCollection was switched; that is fixed now.
#     . GNV.addAgingWGCNAcollection, GNV.addCircadianWGCNAcollection and GNV.addWTGNVModuleCollection have changes in
#       arguments: instead of rootDir, they take directly file names whose defaults reproduce default behavior.
#     . Bugfix in GNV.addAgingWGCNAcollection: by default, the function did not seem to include aging WGCNA modules
#     . GNV.addWTGNVModuleCollection is changed to GNV.addStriatumCoExMap and its output is also named StriatumCoExMap

#=======================================================================================================
#
# Create designs and contrasts
#
#=======================================================================================================

GNV.createDesignsAndContrasts = function(wt = "WT", background="Q175", perturbation, 
                          pertBack = spaste(background, "/", perturbation), 
                          genotypeName = "Genotype", covariates = NULL,
                          testNameBase = "", 
                          testNameBase.pretty = "", 
                          levelSep = ".vs.",
                          levelSep.pretty = " vs. ",
                          all = TRUE,
                          testSubsets = FALSE,
                          # should significant hits up and down be combined for enirchment analysis (in
                          # addition to direction-specific list)? Only affects first two tests.
                          combineSignificantHits = FALSE)
{
   levels = c(wt, perturbation, background, pertBack);
   nLevels = length(levels);
   levelMatrix = matrix(levels, nrow = length(levels), ncol = length(levels));
   # manual ordering:
   contrastOrder = c(1,6,2,3,5,4);
   #contrastOrder = c(1:6);
   lev1 = t(levelMatrix)[lower.tri(levelMatrix)][contrastOrder];  # fold change denominators
   lev2 = levelMatrix[lower.tri(levelMatrix)][contrastOrder];  # fold change numerators
   lev1.n = make.names(lev1)
   lev2.n = make.names(lev2)
   nContrasts = length(lev1)
   contrasts = do.call(rbind, mymapply(c, rep(genotypeName, nContrasts), lev2.n, lev1.n))
   colnames(contrasts) = c("contrast.name", "contrast.numerator", "contrast.denominator");
   testNames.pretty = spaste( testNameBase.pretty, lev2, levelSep.pretty, lev1);
   testNames = spaste( testNameBase, lev2.n, levelSep, lev1.n);
   if (testSubsets) designs0 = testNames else designs0 = rep(genotypeName, nContrasts);

   designs = paste("~", designs0, prependPrefix("+ ", base::paste(covariates, collapse = " + ")));
   if (testSubsets)
   {
     out = data.frame(
        design = designs, contrast.name = rep(NA, nContrasts),
        contrast.numerator = rep(NA, nContrasts),
        contrast.denominator = rep(NA, nContrasts), 
        coefficientName = testNames,
        testName = testNames,
        testName.pretty = testNames.pretty,
        combineSignificantHits = c(rep(combineSignificantHits, 2), rep(FALSE, nContrasts-2)),
        splitSignificantHitsBy = rep(NA, nContrasts)
        );
   } else
     out = data.frame(design = designs, contrasts, 
        coefficientName = rep(NA, nContrasts),
        testName = testNames,
        testName.pretty = testNames.pretty,
        combineSignificantHits = c(rep(combineSignificantHits, 2), rep(FALSE, nContrasts-2)),
        splitSignificantHitsBy = rep(NA, nContrasts));
   if (!all)
   {
     keepIndex = c(1:4);
     out = out[keepIndex, ];
   }
   out;
}






# Alternative matrix for when phenotypes are pairwise binarized 
GNV.createDesignsForPairwisePheno = function(wt = "WT", background="Q175", perturbation,
                          pertBack = spaste(background, "/", perturbation),
                          genotypeName = "Genotype", covariates = NULL,
                          testNameBase = "",
                          testNameBase.pretty = "",
                          levelSep = ".vs.",
                          levelSep.pretty = " vs. ",
                          all = TRUE,
                          # should significant hits up and down be combined for enirchment analysis (in
                          # addition to direction-specific list)? Only affects first two tests.
                          combineSignificantHits = FALSE)
{
   levels = c(wt, perturbation, background, pertBack);
   nLevels = length(levels);
   levelMatrix = matrix(levels, nrow = length(levels), ncol = length(levels));
   # manual ordering:
   contrastOrder = c(1,6,2,3,5,4);
   #contrastOrder = c(1:6);
   lev1 = t(levelMatrix)[lower.tri(levelMatrix)][contrastOrder];  # fold change denominators
   lev2 = levelMatrix[lower.tri(levelMatrix)][contrastOrder];  # fold change numerators
   lev1.n = make.names(lev1)
   lev2.n = make.names(lev2)
   nContrasts = length(lev1)
   tests = make.names(spaste(genotypeName, ".", lev2, levelSep, lev1));
   contrasts = matrix(NA_character_, nContrasts, 3);
   colnames(contrasts) = c("contrast.name", "contrast.numerator", "contrast.denominator");
   designs = paste("~", sapply(tests, base::paste, covariates, collapse = " + "))
   designs = sub("~ *\\+", "~", designs);
   testNames.pretty = spaste( testNameBase.pretty, lev2, levelSep.pretty, lev1);
   testNames = spaste( testNameBase, lev2.n, levelSep, lev1.n);
   out = data.frame(
       design = designs, contrasts,
       coefficientName = tests,
       testName = testNames,
       testName.pretty = testNames.pretty,
       combineSignificantHits = c(rep(combineSignificantHits, 2), rep(FALSE, nContrasts-2)),
       splitSignificantHitsBy = rep(NA, length(designs)));
   if (!all)
   {
     keepIndex = c(1:4);
     out = out[keepIndex, ];
   }
   out;
}

extendedResults = function(analysis)
{
  if (is.null(analysis$combinedDESeqResults.ext)) analysis$combinedDESeqResults else analysis$combinedDESeqResults.ext;
}

#=======================================================================================================
#
# addRescueColumns
#
#=======================================================================================================

rescueScoreFnc = function(baselineStat, perturbationStat, inverse = TRUE)
{
  c(1, -1)[inverse+1] * replaceMissing(pmin(abs(baselineStat), abs(perturbationStat)) *
                                sign(baselineStat* perturbationStat))
}
  

# This adds rescue information columns. It also removes redundant columns, typically baseMean which under
# default settings will be the same for all comparisons.
addRescueColumns = function(
    normalizedData,
    dataIsLogTransformed = FALSE,
    logBase = 2,
    logShift = 1,
    pheno,
    res,
    #testNames,
    #combName,
    wt = "WT", background="Q175", perturbation,
    pertBack = spaste(background, "/", perturbation),
    levelSep = ".vs.",
    testNameSep = ".for.",
    pattern.Z = spaste("Z", testNameSep),
    pattern.p = spaste("p", testNameSep),
    pattern.FDR = spaste("FDR", testNameSep),
    pattern.mean = spaste("baseMean", testNameSep),
    pattern.LFC = spaste("log2FC", testNameSep),
    stopOnMultipleMatches = TRUE,
    genotypeColumnName = "Genotype",
    topThresholds,
    topThresholdNames,
    baselineLog2FC = NULL,
    baselineZ = NULL,
    baselineP = NULL,
    baselineFDR = NULL,
    rescueColNameExtension = "",
    ...)
{
  if (!genotypeColumnName %in% colnames(pheno))
     stop("addRescueColumns: 'genotypeColumnName' ", genotypeColumnName, 
          " was not found among column names of pheno.");
  group = make.names(pheno[, genotypeColumnName]);
  if (dataIsLogTransformed) normalizedData = logBase^normalizedData - logShift;

  pertBack.n = make.names(pertBack);
  bg.n = make.names(background);
  wt.n = make.names(wt);
  bg.vs.wt = spaste(bg.n, levelSep, wt.n);
  pb.vs.bg = spaste(pertBack.n, levelSep, bg.n);
  pb.vs.wt = spaste(pertBack.n, levelSep, wt.n);

  if (length(baselineZ) > 0 || length(baselineLog2FC) > 0)
  {
    if (length(baselineZ) != length(baselineLog2FC))
       stop("When 'baselineZ' or 'baselineLog2FC' are given, they must be vectors of the same length");
    printFlush("   addRescueColumns: using external baseline data.");
    Z.bg.vs.wt = baselineZ;
    if (is.null(baselineP)) baselineP = 2*pnorm(abs(baselineZ), lower.tail = FALSE);
    if (is.null(baselineFDR)) baselineFDR = p.adjust(baselineP);
    p.bg.vs.wt = baselineP;
    FDR.bg.vs.wt = baselineFDR;
    LFC.pb.vs.bg = res[[grepSingleOrError(spaste(pattern.LFC, pb.vs.bg), names(res), fixed = TRUE,
                                          stopOnMultiple = stopOnMultipleMatches)]];
    # Note that the rescue fraction may be biased here
    rescueFraction = (2^LFC.pb.vs.bg - 1) / (2^(-baselineLog2FC) - 1)
  } else {
    p.bg.vs.wt = res[[grepSingleOrError(spaste(pattern.p, bg.vs.wt), names(res), fixed = TRUE,
                                          stopOnMultiple = stopOnMultipleMatches)]];
    FDR.bg.vs.wt = res[[grepSingleOrError(spaste(pattern.FDR, bg.vs.wt), names(res), fixed = TRUE,
                                          stopOnMultiple = stopOnMultipleMatches)]];
    Z.bg.vs.wt = res[[grepSingleOrError(spaste(pattern.Z, bg.vs.wt), names(res), fixed = TRUE,
                                          stopOnMultiple = stopOnMultipleMatches)]];
    baseMeans = t(colStatsByGroup(normalizedData, group)$mean);
    rescueFraction = (baseMeans[, pertBack.n] - baseMeans[, bg.n])/
                      (baseMeans[, wt.n] - baseMeans[, bg.n]);
  }
    
  FDR.pb.vs.bg = res[[grepSingleOrError(spaste(pattern.FDR, pb.vs.bg), names(res), fixed = TRUE,
                                          stopOnMultiple = stopOnMultipleMatches)]];
  p.pb.vs.bg = res[[grepSingleOrError(spaste(pattern.p, pb.vs.bg), names(res), fixed = TRUE,
                                          stopOnMultiple = stopOnMultipleMatches)]];

  # Note: the move could also be in the wrong direction
  significantMove = mymapply(function(tt)
    replaceMissing(FDR.bg.vs.wt < tt$pAdj & FDR.pb.vs.bg < tt$pAdj & 
                   p.bg.vs.wt < tt$p & p.pb.vs.bg < tt$p),
    topThresholds)

  significantRescue = do.call(cbind, lapply(significantMove, function(sm) 
             replaceMissing(as.numeric(sm & rescueFraction > 0))));

  significantExacerbation = do.call(cbind, lapply(significantMove, function(sm) 
             replaceMissing(as.numeric(sm & rescueFraction < 0))));

  significantMove2 = do.call(cbind, lapply(significantMove, as.numeric));

  if (length(rescueColNameExtension) == 0) rescueColNameExtension = "";

  colnames(significantMove2) = spaste("significantMove", rescueColNameExtension, topThresholdNames);
  colnames(significantRescue) = spaste("significantRescue", rescueColNameExtension, topThresholdNames);
  colnames(significantExacerbation) = spaste("significantExacerbation", rescueColNameExtension, topThresholdNames);

  rescueInd = interleave(list(significantMove2, significantRescue, significantExacerbation),
                          nameBase = rep("", 3), sep = "");
  # Define a rescue score. Can define it as the minimum of abs(Z.bg.vs.wt), abs(Z.pb.vs.bg) multiplied by a
  # sign such that a move in the rescue direction is positive and in the exacerbation is negative.
  # Rescue direction: Z.bg.vs.wt * Z.pb.vs.bg < 0

  Z.pb.vs.bg = res[[grepSingleOrError(spaste(pattern.Z, pb.vs.bg), names(res), 
                                          stopOnMultiple = stopOnMultipleMatches)]];
  rescueScore = rescueScoreFnc(Z.bg.vs.wt, Z.pb.vs.bg, inverse = TRUE);
  rescueInd2 = cbind(rescueScore = rescueScore, rescueFraction = rescueFraction);
  cbind(rescueInd, setColnames(rescueInd2, spaste(colnames(rescueInd2), rescueColNameExtension)), res); 
}


#=======================================================================================================
#
# rescueGeneLists
#
#=======================================================================================================

# This gets slightly tricky. We want to split the genes by the sign of differential expression but select top
# genes by rescue score except statistically full rescue where I don't have good ordering statistic. 

rescueGeneLists = function(
   results,
   topList.minN = 0,
   geneIDs,
   backupGeneIDs = NULL,
   wt = "WT", background="Q175", perturbation,
   pertBack = spaste(background, "/", perturbation),
   pattern.Z = spaste("Z", testNameSep),
   levelSep = ".vs.",
   signColName = spaste(pattern.Z, make.names(background), levelSep, make.names(wt)),
   signCol = results[[signColName]],
   testNameSep = ".for.",
   andSep = ".and.",
   nameModifier = spaste(make.names(background), levelSep, make.names(wt)),
   sourceCols = c("significantMove", "significantRescue", "significantExacerbation"),
   orderCols = c("rescueScore", "rescueScore", "rescueScore"),
   orderSign = c(0, -1, 1),
   topThreshold1 = NULL,
   ...)
{
  #bg.n = make.names(background);
  #wt.n = make.names(wt);

  #if (is.null(geneIDs)) geneIDs = rownames(results);
  # geneIDs must be explicitly given

  cols1 = grep(sourceCols[1], names(results));
  suffixes = sub(sourceCols[1], "", names(results)[cols1]);  ## these are typically the significance thresholds
  if (!is.null(topThreshold1)) 
  {
    if (!topThreshold1 %in% suffixes) 
      stop("rescueGeneLists: given topThreshold1 is not among the suffixes.");
    suffixes = topThreshold1;
  }
  out = unlist(lapply(suffixes, function(suf)
  {
    unlist(removeNames(mymapply(function(colName, orderColName, orderSign)
    {
      colName1 = spaste(colName, suf);
      colInd1 = match(colName1, names(results));
      if (length(colInd1)==0) stop("Could not find column ", colName1, " among columns of 'results'.");
      ind = results[, colInd1];
      if (!is.na(orderColName))
      {
        orderColInd = match(orderColName, names(results));
        if (is.na(orderColInd)) stop("Could not find column ", orderColName, " among columns of 'results'.");
        orderStat = if (orderSign==0) -abs(results[, orderColInd]) else orderSign * results[, orderColInd];
        out = lapply(c(-1, 1), function(dir)
        {
          selData = cbind(p = 1-results[, colInd1],
                          Z = orderStat);
          selData[dir * signCol < 0, ] = NA;
          significantOrTop(IDs = geneIDs, 
                           backupIDs = backupGeneIDs,
                           Z = selData[, 2], 
                           p = selData[, 1], pThreshold = 0.5, 
                           p.adjusted = selData[, 1], pAdjThreshold = 0.5,
                           nTop = topList.minN);
        })
      } else {
        out = lapply(c(-1, 1), function(dir)
          geneIDs[dir * signCol > 0 & ind > 0]);
      }
      names(out) = spaste(names(results)[colInd1], andSep, c("down", "up"), testNameSep, nameModifier);
      out;
    }, sourceCols, orderCols, orderSign)), recursive = FALSE)
  }), recursive = FALSE);
  out;
}


#=======================================================================================================
#
# prepare standard collections for enrichment analysis
#
#=======================================================================================================

GNV.addWGCNAcollection = function(pathToRLibs, organism = "mouse",
                               tissue = "striatum", WGCNA = c("Allelic Series"),
                               collections = NULL,
                               addCombined = TRUE, combinedName = "AllCombined",
                               WGCNAName = "WGCNA")
{
  require(HuntingtonDiseaseWGCNACollection);
  WGCNACollection = NULL;
  if ("Allelic Series" %in% WGCNA)
  {
    WGCNACollection = 
        mergeCollections(WGCNACollection, 
               subsetCollection(HuntingtonDiseaseWGCNACollection(organism = organism),
               tags = spaste("C-WGCNA of 2-, 6-, 10-month CHDI AS ", tolower(tissue)),
               exactMatch = FALSE));
  } 
  if ("BACHD.dN17" %in% WGCNA)
  {
    dn17 = loadAsList(file.path(pathToRLibs, "../HuntingtonsDisease",
             "IndividualAnalyses/180-Cantle-BACHD-dN17/030-CreateAnRichmentCollection-fromV2",
             "020-Collection-RObject/modules.BACHD.dN17.rda"))
    WGCNACollection = 
        mergeCollections(WGCNACollection,
               subsetCollection(dn17$modules.BACHD.dN17,
               tags = spaste("WGCNA of BACHD-dN17 ", tolower(tissue)),
               exactMatch = FALSE));
    if (organismLabels(organism)[1]!="mouse") 
      WGCNACollection = convertCollectionToOrganism(WGCNACollection, organism = "mouse");
  }
  collections = c(collections, list(xxx123 = WGCNACollection))
  names(collections) = sub("xxx123", WGCNAName, names(collections));
  if (addCombined)
  {
    n = length(collections)+1;
    collections[[n]] = do.call(mergeCollections, collections);
    names(collections)[n] = combinedName;
  }
  collections;
}


GNV.addAgingWGCNAcollection = function(organism = "mouse",
                               getDECollection = TRUE,
                               tissue = "Striatum", 
                               highKME = FALSE,
                               collections = NULL, 
                               removeInteractionSets = TRUE,
                               WGCNAName = "AgingWGCNA",
                               addCombined = TRUE, combinedName = "AllCombined",
                               file.DE = NULL,
                               file.WGCNA = NULL)
{
  if (tolower(tissue)=="striatum")
  {
    if (is.null(file.WGCNA))
      file.WGCNA = file.path( getRootDir(), "Aging/AgingInWTMice-NanWang-2018-9096/030-NetworkAnalysis/RData",
                    if (highKME) "moduleCollection-WTAging.Striatum-highKME.RData" else
                                  "moduleCollection-WTAging.Striatum.RData");
    if (is.null(file.DE))
      file.DE = file.path(getRootDir(), "Aging/AgingInWTMice-NanWang-2018-9096/020-IndividualAnalysis/RData",
                          "collectionOfDEGenes-WTAgingStriatum.RData");
  } else {
    if (is.null(file.WGCNA))
      file.WGCNA = file.path(getRootDir(), 
                     "HuntingtonsDisease/OutputResults/200-CreateWGCNAAgingModuleCollection/Collection-RObject", 
                     "WGCNA.Age.PL.collection.RData");
    if (is.null(file.DE)) file.DE = file.path(getRootDir(), "HuntingtonsDisease/IndividualAnalyses/070-CHDIAllelicSeries",
          "CommonAnalysis-040-mRNA-aging/015-StandardAnalysis/RData/AgeAssociationCollection.RData");
  }

  x = loadAsList(file.WGCNA)[[1]];
  x2 = loadAsList(file.DE)$assocCollection;
  if (removeInteractionSets)
  {
    # Get rid of interaction sets...
    x2 = subsetCollection(x2, tags = "interaction", exactMatch = FALSE, invertSearch = TRUE);
  }

  if (getDECollection)
  {
    coll= mergeCollections(x, x2);
  } else coll = x;

  if (!is.null(tissue))
    coll = subsetCollection(coll, tags = tissue, matchComponents = "groups", exactMatch = FALSE, ignore.case = TRUE);
     
  collections = c(collections,
        list(x128372874 = convertCollectionToOrganism(coll, organism = organism)));
  names(collections) = sub("^x128372874$", WGCNAName, names(collections));

  if (addCombined)
  {
    n = length(collections)+1;
    collections[[n]] = do.call(mergeCollections, collections);
    names(collections)[n] = combinedName;
  }
  collections;
}

GNV.addStriatumCoExMap= function(
   organism = "mouse",
   highKME = FALSE,
   file = NULL,
   collections = NULL,
   addCombined = TRUE, combinedName = "AllCombined")
{
  if (length(file)==0) file = file.path(getRootDir(),
          "HuntingtonsDisease/IndividualAnalyses/GNVProject/CommonAnalyses/",
          "160-WGNCA-perGenotype/030-NetworkAnalysis/RData/",
          if (highKME) "moduleCollection-WT-stringentKME.RData" else
                 "moduleCollection-WT-stringentKME.RData");

  coll = loadAsList(file)[[1]];
  collections = c(collections,
        list(StriatumCoExMap = convertCollectionToOrganism(coll, organism = organism)));

  if (addCombined)
  {
    n = length(collections)+1;
    collections[[n]] = do.call(mergeCollections, collections);
    names(collections)[n] = combinedName;
  }
  collections;
}

GNV.addCircadianWGCNAcollection = function(
     organism = "mouse",
     getDECollection = TRUE,
     tissue = "striatum", 
     highKME = FALSE,
     defaultTissue = NULL,
     file.WGCNA = NULL,
     file.DE = NULL,
     collections = NULL, 
     addCombined = TRUE, combinedName = "AllCombined")
{
  availableTissues = c("Striatum", "Cortex", "Hypothalamus", "Liver");

  keepSV = c("with 6 SV", "with 6 SV", "with 6 SV", "with 0 SV");

  tsIndex = match(tolower(tissue), tolower(availableTissues));
  if (is.na(tsIndex) && !is.null(defaultTissue)) 
  {  
      tissue = defaultTissue;
      tsIndex = match(tolower(defaultTissue), tolower(availableTissues));
  }

  if (is.na(tsIndex))
     stop("Neither 'tissue' (", tissue, ") nor 'defaultTissue' are among the available tissues.\n",
          "Available tissues: ", paste(availableTissues, collapse = ", "));
  
  if (is.null(file.WGCNA))
    file.WGCNA = spaste(getRootDir(),
          "/HuntingtonsDisease/IndividualAnalyses/500-CHDI-CircadianGeneExpression/030-NetworkAnalysis/RData/",
          "moduleCollection-selected-", availableTissues[tsIndex], 
          if (highKME) "-highKME" else "", ".RData");

  if (is.null(file.DE))
    file.DE = spaste(getRootDir(),
          "/HuntingtonsDisease/IndividualAnalyses/500-CHDI-CircadianGeneExpression/020-IndividualAnalysis/RData/",
          "collectionOfDEGenes-CircadianExpression.", availableTissues[tsIndex], ".RData");

  x = loadAsList(file.WGCNA)[[1]];
  x2 = loadAsList(file.DE)$assocCollection;

  # Select subset
  x2 = subsetCollection(x2, tags = keepSV[tsIndex], exactMatch = FALSE, ignore.case = TRUE);

  if (getDECollection)
  {
    coll= mergeCollections(x, x2);
  } else coll = x;

  if (!is.null(tissue))
    coll = subsetCollection(coll, tags = tolower(tissue), matchComponents = "groups", exactMatch = FALSE);
     
  collections = c(collections,
        list(CircadianWGCNA = convertCollectionToOrganism(coll, organism = organism)));

  if (addCombined)
  {
    n = length(collections)+1;
    collections[[n]] = do.call(mergeCollections, collections);
    names(collections)[n] = combinedName;
  }
  collections;
}

GNV.standardCollections = function(pathToRLibs, organism = "mouse",
    removeAllelicSeriesWGCNA = FALSE, removeRedundantKelleySets = TRUE,
    keepKelleyRegions = c("STR", "all brain"),
    removeGokcePattern = c("Top 50 genes in striatum.*Gokce SCS 2016"),
    keepSaundersTissue = "striatum",
    separateCircadianSets = FALSE, circadianPattern = "Circadian oscillations",
    genomicSpacings = 5e6,
    MSigDBxml = MSigDBFile(), getMiRNATargets = TRUE)
{   
  require(anRichment)
  #source(file.path(pathToRLibs, "inHouseGeneAnnotation/anRichmentMethods/R/enrichmentAnalysis.R"))

  collections1 = rearrangeInternalCollections(
     allCollections(organism = organism, merge = FALSE, buildExternal = TRUE,
                  MSigDBxml = MSigDBxml, genomicSpacings = genomicSpacings),
     separateCircadianSets = separateCircadianSets, circadianPattern = circadianPattern);

  removeSets = character(0);
  if (removeRedundantKelleySets)
  {
    removeSets = c(removeSets, multiGrep(keepKelleyRegions,
                  multiGrep("Kelley 2018", dataSetNames(collections1[["CellAndRegionMarkers"]]), value = TRUE),
                  value = TRUE, invert = TRUE));
  }

  if (length(removeGokcePattern) > 0)
  {
    removeSets = c(removeSets, multiGrepv(removeGokcePattern, dataSetNames(collections1[["CellAndRegionMarkers"]])));
  }

  if (length(keepSaundersTissue) > 0)
  {
    sets1 = grepv("\\(Saunders 2018\\)", dataSetNames(collections1[["CellAndRegionMarkers"]]));
    removeSets = c(removeSets, multiGrepv(keepSaundersTissue, sets1, invert = TRUE));
  }

  if (length(removeSets) > 0)
    collections1[["CellAndRegionMarkers"]] = subsetCollection(collections1[["CellAndRegionMarkers"]],
                                                       tags = removeSets, invertSearch = TRUE);

  if (removeAllelicSeriesWGCNA)
  {
    collections1$mRNA.WGCNA = subsetCollection(collections1$mRNA.WGCNA, 
      tags = c("RNA M[0-9]* module .*Langfelder, "),
      exactMatch = FALSE, fixed = FALSE, invertSearch = TRUE);
  }

  org.std = organismLabels(organism)
  if (org.std=="mouse" && getMiRNATargets) {
    mirTargets = loadMirnaTargets(file.path(pathToRLibs, "../Data/miRNA/Mouse"));
  } else mirTargets = NULL;

  if (!is.null(mirTargets))
  {
    mirColl = collectionOfMirnaTargets(mirTargets)
    c(collections1, list(miRNATargets = mirColl));
  } else collections1;
}


#=======================================================================================================
#
# plotPCA
#
#=======================================================================================================

plotPCA.2panels = function(
  data1,
  data2 = NULL,

  colorFactor = NULL,
  shapeFactor = NULL,
  borderFactor = NULL,
  pointLabels = NULL,
  
  colorPrefix = "",
  shapePrefix = "",
  scale = TRUE,
  device = "pdf",
  dir = "",
  file = NULL,

  main1 = "PCA",
  main2 = "",

  includePVE = TRUE,

  plotLegend = TRUE,
  separateLegend = TRUE,
  legendWithinPlot = !separateLegend,
  legendWidth = 0.4,
  setLayout = TRUE,
  height = 4,
  width = height*(1 + !is.null(data2) + legendWidth * separateLegend), 
  pt.cex = 2.5,
  pt.bg = "grey70",
  thinned = FALSE,
  firstShape = 21,
  colorSequence = standardColors(),
  cex.pointLabels = 1,
  fig.cex = 1,
  lwd = 1,
  mar = c(3.2, 3.2, 3, 1),
  mgp = c(2, 0.7, 0),
  legendOnly = FALSE,
  prettifyList = NULL,
  ...)
{
  if (!is.null(device) && !is.null(file))
  {
    device = match.fun(device);
    dir.create(dir, recursive = TRUE, showWarnings = FALSE);
    device(file = file.path(dir, file), wi = width, he = height);
    on.exit(dev.off());
  }
  do2 = !is.null(data2);
  if (separateLegend && plotLegend && !is.null(colorFactor))
  {
    if (setLayout)
      layout(matrix(c(1:(2 + do2)), 1, 2 + do2),
                    widths = c(legendWidth, rep(1, 1+do2)));
    cfLevels = levels(as.factor(colorFactor));
    nCFLevels = length(cfLevels);
    if (!is.null(shapeFactor))
    {
      shape = firstShape - 1 + as.numeric(as.factor(shapeFactor));
      shapeLevels = levels(as.factor(shapeFactor));
    } else {
      shape = firstShape;
      shapeLevels = NULL;
    }
    if (!is.null(borderFactor))
    {
      borderLevels = levels(as.factor(borderFactor));
      nBorderLevels = length(borderLevels);
      if (any(is.na(suppressWarnings(as.numeric(borderLevels)))))
      {
        borderColors = labels2colors(as.numeric(as.factor(borderFactor)), colorSeq = colorSequence);
        borderColorLevels = labels2colors(1:nBorderLevels, colorSeq = colorSequence);
      } else {
        borderColors = as.numeric(as.character(borderFactor));
        borderColorLevels = sort(unique(borderColors));
      }
    } else {
      borderColors = 1
      borderColorLevels = 1;
    }
    if (!is.null(prettifyList))
    {
      shapePrefix = prettifyStrings(shapePrefix, prettifyList);
      shapeLevels = prettifyStrings(shapeLevels, prettifyList);
      cfLevels = prettifyStrings(cfLevels, prettifyList);
    }
    nShapeLevels = length(shapeLevels)
    par(mar = c(mar[1], 0.2, mar[3], 0.2));
    par(cex = fig.cex);
    plot(c(0,1), type = "n", axes = FALSE, frame = FALSE, xlab = "", ylab = "")
    out = legendClean("top",
              legend = c(cfLevels, spaste(rep(shapePrefix, nShapeLevels), shapeLevels)),
              pt.cex = min(pt.cex, 1.5),
              pch = c(rep(firstShape, nCFLevels), firstShape - 1+seq_len(nShapeLevels)),
              pt.bg = labels2colors(c(1:nCFLevels, rep(0, nShapeLevels)), colorSeq = colorSequence),
              col = c(.extend(borderColorLevels, nCFLevels), rep(borderColorLevels[1], nShapeLevels)),
              pt.lwd = lwd);
    if (legendOnly) return(out);
  } else if (setLayout) par(mfrow = c(1,1 + do2))

  plotPCA(data1, genotype = colorFactor,
          shapeFactor = shapeFactor,
          borderFactor = borderFactor,
          pointLabels = pointLabels,
          colorPrefix = colorPrefix,
          shapePrefix = shapePrefix,
          device = NULL,
          analysisID = NULL,
          file = NULL,
          plotLegend = legendWithinPlot,
          main = main1,
          dataStep = "",
          scale = scale,
          pt.cex = pt.cex,
          pt.bg = pt.bg,
          thinned = thinned,
          firstShape = firstShape,
          lwd = lwd,
          colorSequence = colorSequence,
          cex.pointLabels = cex.pointLabels,
          fig.cex = fig.cex,
          mar = mar,
          mgp = mgp,
          prettifyList = prettifyList,
          includePVE = includePVE,
          ...);

  if (do2) plotPCA(data2, genotype = colorFactor,
          shapeFactor = shapeFactor,
          borderFactor = borderFactor,
          pointLabels = pointLabels,
          shapePrefix = shapePrefix,
          device = NULL,
          analysisID = NULL,
          file = NULL,
          plotLegend = legendWithinPlot,
          main = main2,
          dataStep = "",
          scale = scale,
          pt.cex = pt.cex,
          pt.bg = pt.bg,
          thinned = thinned,
          firstShape = firstShape,
          lwd = lwd,
          colorSequence = colorSequence,
          cex.pointLabels = cex.pointLabels,
          fig.cex = fig.cex,
          mar = mar,
          mgp = mgp,
          prettifyList = prettifyList,
          includePVE = includePVE,
          ...);
}

# A wrapper for plotPCA.2panels  with frequently used presets

plotPCA.asIsAndBackgroundCorrection = function(
   expr, genotype,
   setName,
   setName.pretty,
   background,
   colorFactor = fix.dN17(genotype),
   shapeFactor = NULL,
   analysisID = setName,
   dataStep = "",
   scale = TRUE,
   device = "pdf",
   plotDir = "Plots",
   file = spaste(analysisID, "-PCA", prependDash(dataStep), ".pdf"),
   mainBase = "PCA after preprocessing\n", 
   correctedMain = spaste("Corrected for Htt genotype (", fix.dN17(background), ")"),
   plotLegend = TRUE,
   separateLegend = TRUE,
   legendWithinPlot = !separateLegend,
   legendWidth = 0.4,
   plotCorrected = TRUE,
   setLayout = TRUE,
   width = 4*(1 + plotCorrected + legendWidth * separateLegend), height = 4,
   pt.cex = 2.5,
   pt.bg = "grey70",
   thinned = FALSE,
   fig.cex = 1,
   mar = c(3.2, 3.2, 3, 1),
   mgp = c(2, 0.7, 0),
   shapePrefix = "",
   legendOnly = FALSE,
   ...)
{
  if (plotCorrected)
  {
    bgGenotype = grepl(background, genotype) + 0;
    genoCorrectedData = residuals(lm(expr~bgGenotype));
  }

  plotPCA.2panels(
    data1 = expr,
    data2 = if(plotCorrected) genoCorrectedData else NULL,

    colorFactor = colorFactor,
    colorPrefix = "",
    shapeFactor = shapeFactor,
    shapePrefix = shapePrefix,

    scale = scale,
    device = device,
    dir = plotDir,
    file = file,

    main1 = spaste(mainBase, setName.pretty),
    main2 = correctedMain,

    plotLegend = plotLegend,
    separateLegend = separateLegend,
    legendWithinPlot = legendWithinPlot,
    legendWidth = legendWidth,
    setLayout = setLayout,
    height = height,
    width = width,
    pt.cex = pt.cex,
    pt.bg = pt.bg,
    thinned= thinned,
    fig.cex = fig.cex,
    mar = mar,
    mgp = mgp,
    legendOnly = legendOnly,
    ...);

}

plotPCA = function(
  expr, genotype, 
  shapeFactor = NULL,
  borderFactor = NULL,
  pointLabels = NULL,
  colorPrefix = "",
  shapePrefix = "",
  
  analysisID = "",
  dataStep = "", 
  scale = TRUE, 
  plotDir = "Plots",
  file = spaste(analysisID, "-PCA", prependDash(dataStep), ".pdf"),
  width = 5, height = 5, 
  pt.cex = 2.5,
  pt.bg = "grey70",
  thinned = FALSE,
  firstShape = 21,
  colorSequence = standardColors(),
  lwd = 1,
  fig.cex = 1,
  cex.pointLabels = 1,
  mar = c(3.2, 3.2, 3, 1),
  mgp = c(2, 0.7, 0),
  device = "pdf",
  plotLegend = TRUE,
  prettifyList = NULL,
  includePVE = TRUE,
  ...)
{
  if (scale) expr = scale(expr);
  if (any(is.na(expr)))
  {
    gsg = goodSamplesGenes(expr);
    if (!gsg$allOK) expr = expr[, gsg$goodGenes];
    if (any(is.na(expr)))
      expr = t(impute.knn(t(expr))$data);
  }
  svd = svd(expr, nu = 2, nv = 0);

  nSamples = nrow(expr);
  if (length(genotype)>0 && length(genotype)!=nSamples) 
       browser()
       #stop("Number of samples in 'expr' and in 'genotype' disagree.");
  if (length(shapeFactor)>0 && length(shapeFactor)!=nSamples) 
       stop("Number of samples in 'expr' and in 'shapeFactor' disagree.");

  if (!is.null(genotype))
  {
    geno.num = as.numeric(as.factor(genotype));
    genoLevels = levels(as.factor(genotype));
    bgColors = labels2colors(geno.num, colorSeq = colorSequence);
  } else {
    bgColors = pt.bg;
    genoLevels = NULL;
  }
  if (!is.null(borderFactor))
  {
      if (any(is.na(suppressWarnings(as.numeric(as.character(borderFactor))))))
      {
        borderColors = labels2colors(as.numeric(as.factor(borderFactor)), colorSeq = colorSequence);
      } else {
        borderColors = as.numeric(as.character(borderFactor));
      }
  } else
    borderColors = 1;
  nGenoLevels = length(genoLevels)
  if (!is.null(shapeFactor))
  {
    shape = firstShape - 1 + as.numeric(as.factor(shapeFactor));
    shapeLevels = levels(as.factor(shapeFactor));
  } else {
    shape = firstShape;
    shapeLevels = NULL;
  }
  nShapeLevels = length(shapeLevels)

  if (!is.null(file))
  {
    devFnc = match.fun(device);
    dir.create(plotDir, recursive = TRUE, showWarnings = FALSE);
    devFnc(file = file.path(plotDir, file), wi = width, he = height);
    on.exit(dev.off());
  }

  par(cex = fig.cex);
  par(mar = mar);
  par(mgp = mgp);

  if (includePVE) {
    PVE = svd$d^2/sum(svd$d^2);
    PVE1 = spaste(": PVE = ", round(PVE[1] * 100)/100);
    PVE2 = spaste(": PVE = ", round(PVE[2] * 100)/100);
  } else {
    PVE1 = PVE2 = "";
  }
  plotFnc = if (thinned) match.fun("thinnedScatterplot") else match.fun("plot")
  plotFnc(svd$u[, 1], svd$u[, 2],
       xlab = spaste("PC1", PVE1), ylab = spaste("PC2", PVE2),
       cex = pt.cex, pch = shape,
       bg = bgColors,
       col = borderColors, lwd = lwd,
       ...)

  if (!is.null(pointLabels))
  {
    aa = try(labelPoints2(svd$u[, 1], svd$u[, 2], pointLabels, 
                 cex = cex.pointLabels));
    if (inherits(aa, "try-error")) browser();
  }

  if (!is.null(prettifyList))
  {
    shapePrefix = prettifyStrings(shapePrefix, prettifyList);
    shapeLevels = prettifyStrings(shapeLevels, prettifyList);
    genoLevels = prettifyStrings(genoLevels, prettifyList);
  }
  if (plotLegend && nGenoLevels > 0)
    legendClean("auto", points.x = svd$u[, 1], points.y = svd$u[, 2],
              tryNCol = c(1:4),
              nPositions = 10,
              legend = c(spaste(colorPrefix, genoLevels), spaste(shapePrefix, shapeLevels)),
              pt.cex = min(pt.cex, 1.5),
              pch = c(rep(firstShape, nGenoLevels), firstShape - 1 +seq_len(nShapeLevels)), 
              pt.bg = labels2colors(c(1:nGenoLevels, rep(27, nShapeLevels)), colorSeq = colorSequence),
              col = 1);

  #if (!is.null(file)) dev.off();
  NULL;
}


#=========================================================================================================
#
# Add a table of the number of rescued/moved genes to analysis results
#
#=========================================================================================================

GNV.addNumbersOfRescuedGenes = function(
   analysis,
   prettifyList = analysis$prettifyList,
   testNameSep = ".for.", 
   andSep = ".and.", resultDir = NULL,
   addTopGenes = 15,
   maxAddTopGenes = 20)
{
  lst = analysis$topIDs.extraLists[[1]];
  if (is.null(lst))
    stop("'analysis' does not contain element 'topIDs.extraLists'.")

  nAdd = function(n) if (n<=maxAddTopGenes) n else addTopGenes;

  char.top = lapply(analysis$topSymbols.extraLists[[1]], 
         function(x) if (length(x) > 0) spaste(" (",
              paste(x[1:nAdd(length(x))], collapse = ", "), ")") else "");
      
  lengths = sapply(lst, length);

  lengths.up = lengths[grep(spaste(andSep, "up", testNameSep), names(lst))];
  lengths.down = lengths[grep(spaste(andSep, "down", testNameSep), names(lst))];

  topGenes.up = char.top[grep(spaste(andSep, "up", testNameSep), names(lst))];
  topGenes.down = char.top[grep(spaste(andSep, "down", testNameSep), names(lst))];

  nMoved = cbind(lengths.down, lengths.up);

  colnames(nMoved) = spaste(c("Downregulated", "Upregulated"), testNameSep, 
                            sub(spaste(".+", spaste(andSep, "up", testNameSep)), "", names(lengths.up[1])));
  #rownames(nMoved) = prettifyStrings(
  #                      sub(spaste(spaste(andSep, "up", testNameSep), ".+"), "", names(lengths.up)),
  #                      prettifyList);
  rownames(nMoved) = sub(spaste(spaste(andSep, "up", testNameSep), ".+"), "", names(lengths.up));

  analysis$nMoved = data.frame(nMoved, check.names = FALSE);

  nMoved.ext = nMoved;
  nMoved.ext[, 1] = spaste(nMoved[, 1], topGenes.down);
  nMoved.ext[, 2] = spaste(nMoved[, 2], topGenes.up);

  analysis$nMoved.ext = nMoved.ext;

  if (!is.null(resultDir))
  {
    # dir.create(resultDir, recursive = TRUE, showWarnings = FALSE)
    write.csv(prettifyNames(
                 prettifyColumns(dataForTable(nMoved, transpose = FALSE, IDcolName = "Category"),
                                 "Category", prettifyList = analysis$prettifyList),
                 prettifyList = analysis$prettifyList),
              file = spaste(resultDir, "/numbersOfSignificantlyMovedGenes",
                                prependDash(analysis$analysisName), ".csv"),
              row.names = FALSE);
    write.csv(prettifyNames(
                prettifyColumns(dataForTable(nMoved.ext, transpose = FALSE, IDcolName = "Category"),
                           "Category", prettifyList = analysis$prettifyList),
                prettifyList = analysis$prettifyList),
              file = spaste(resultDir, "/numbersOfSignificantlyMovedGenes-WithTopGenes",
                                prependDash(analysis$analysisName), ".csv"),
              row.names = FALSE);
  }
  analysis;
}
  


#=========================================================================================================
#
# Generate a standard set of comparison plots for GNV analysis
#
#=========================================================================================================

# Plots to generate:
#  . Numbers of significantly associated genes, rescued etc. genes; display top genes in the plot if possible
#  . Enrichment barplots, possibly also with top overlapping genes
#  . Heatmap of correlation of Z statistics for all tests
#  . Scatterplot of Z for genotype between this data and old data
#  . Scatterplot of Z for bg vs wt and Z for bg vs pert/bg

#=========================================================================================================
#
# numbers of significant genes
#
#=========================================================================================================

dummy = function(...) {NULL}

legendForNumberPlots.fromInfo = function(
   plotInfo,
   ...)
{
  n = length(plotInfo$ytop)
  nPerRow = 100;
  xmin = if (is.null(plotInfo$textLeft.x)) rep(0, n) else plotInfo$textLeft.box["xmin", ]
  xmax = if (is.null(plotInfo$textRight.x)) rep(0, n) else plotInfo$textRight.box["xmax", ];
  x = y = numeric(0);
  for (row in 1:n)
  {
    x = c(x, rep(seq(from = xmin[row], to = xmax[row], length.out = nPerRow), 3));
    y = c(y, rep(c(plotInfo$ytop[row], plotInfo$ybottom[row], plotInfo$yMid[row]), each = nPerRow));
  }

  legendClean("auto", points.x = x, points.y = y, ...);
}

  

GNV.plotNumbersOfSignifGenes = function(
   analysis,

   analysisName = "",
   plotDir = "Plots",
   device = "CairoPDF",
   width = 7,
   height = 3,
   baseHeight = 1,
   heightPerRow = 0.25, 
   baseWidth = 5,
   prettifyList = NULL,

   mar = c(0.5, 8, 2, 0.5),
   cex = 1,
   main = spaste("Numbers of significantly associated genes %s", 
            prepComma(prettifyStrings.multi(analysisName, list(analysis$prettifyList.plot,
                                                                                prettifyList)))),
   col.right = "#FF9070", col.left = "skyblue",
   border.right = "brown", border.left = "blue",
   barGap = 0.4,
   plotLegend = FALSE,
   #legendPos = "bottomright",
   legend.cex = 1,
   minLim = 1,
   cex.lab = 1,
   ...
)
{
  yLabels1 = prettifyStrings.multi(rownames(analysis$nSignif.all[[1]]),
                                    list(analysis$prettifyList.plot, prettifyList))
  # First plot: just the nSignif table. There is one nSignif table per significance theshold.
  if (!is.null(device))
  {
    if (is.null(height)) height = baseHeight + heightPerRow * nrow(analysis$nSignif.all[[1]]);
    if (is.null(width))
    {
      mw = marginWidth(text = yLabels1, cex = cex.lab, device = device);
      width = baseWidth + mw$inch;
      mar[2] = mw$lines;
    } 
    device = match.fun(device);
    device(file = file.path(plotDir,
        spaste("numberOfSignificantGenes-genotypeComparisons", prependDash(analysisName),
               ".pdf")),
        width = width, height = height);

    on.exit(dev.off());
  }
  st = prettifyStrings(names(analysis$nSignif.all), prettifyList = analysis$prettifyList);
  for (tt in 1:length(analysis$nSignif.all))
  {
    par(mar = mar);
    par(cex = cex);
    tsb = twoSideBarplot(analysis$nSignif.all[[tt]][, 1], analysis$nSignif.all[[tt]][, 2],
                 col.right = col.right, col.left = col.left,
                 border.right = border.right, border.left = border.left,
                 yLabels = prettifyStrings.multi(rownames(analysis$nSignif.all[[tt]]), 
                                    list(analysis$prettifyList.plot, prettifyList)),
                 barGap = barGap,
                 main = WGCNA:::.substituteTags(main, "%s", st[tt]),
                 minLim = minLim,
                 cex.lab = cex.lab,
                 ...);

    if (plotLegend) legendForNumberPlots.fromInfo(plotInfo = tsb, 
           legend = c("Down", "Up"),
           pch = 22, pt.cex = 1.5 * legend.cex,
           cex = legend.cex, opacity = 0.5,
           pt.bg = c(col.left, col.right),
           col = c(border.left, border.right));
  }

  #dev.off();
}

#=========================================================================================================
#
# Numbers of "moved" genes: significant move, rescue, exacerbation
#
#=========================================================================================================


GNV.plotNumbersOfMovedGenes = function(
   analysis,

   testName.pretty,
   analysisName = "",
   plotDir = "Plots",
   device = "CairoPDF",
   width = 7,
   height = 3,
   prettifyList = NULL,

   splitByThreshold = TRUE,

   addThresholdToMain = splitByThreshold,
   thresholdSep = "\n",

   mar = c(0.5, 12, 3, 0.5),
   cex = 1,
   main = spaste("Numbers of rescued or exacerbated genes %s\n", 
                  prettifyStrings.multi(analysisName, 
                      list(analysis$prettifyList.plot, prettifyList)),
                  "\nsplit by DE direction in ", 
                  prettifyStrings.multi(testName.pretty, list(analysis$prettifyList.plot, prettifyList))),
   col.right = "#FF9070", col.left = "skyblue",
   border.right = "brown", border.left = "blue",
   barGap = 0.4,
   maxCharPerLine = 30,
   plotLegend = TRUE,
   legend.cex = 1,
   minLim = 10,
   dropRowPattern = c("significantMove", "statisticallyFullRescue"),
   ...
)
{
  if (!any(analysis$nMoved > 0)) return(NULL);

  if (splitByThreshold)
  {
    nMoved = lapply(analysis$topThresholdNames, function(tn)
    {
      out = analysis$nMoved[ grep(tn, rownames(analysis$nMoved), fixed = TRUE), , drop = FALSE];
      if (length(dropRowPattern) > 0)
        out = out[ !multiGrepl(dropRowPattern, rownames(out)), , drop = FALSE];
      rownames(out) = gsub(" +$", "", sub(tn, "", rownames(out), fixed = TRUE));
      out;
    });
    names(nMoved) = analysis$topThresholdNames;
  } else
    nMoved = list(analysis$nMoved);

  device = match.fun(device);
  device(file = file.path(plotDir, 
    spaste("numberOfSignificantGenes-rescuedGenes", prependDash(analysisName),
           ".pdf")),
    width = width, height = height);

  on.exit(dev.off());
  par(mar = mar);
  par(cex = cex);
  st = prettifyStrings(names(nMoved), prettifyList = analysis$prettifyList);
  for (tt in 1:length(nMoved))
  {
    tsb = twoSideBarplot(nMoved[[tt]][, 1], nMoved[[tt]][, 2],
               col.right = col.right, col.left = col.left,
               border.right = border.right, border.left = border.left,
               yLabels = formatLabels(
                            prettifyStrings.multi(rownames(nMoved[[tt]]), 
                                         list(analysis$prettifyList.plot, prettifyList)),
                            maxCharPerLine = maxCharPerLine),
               barGap = barGap,
               main = WGCNA:::.substituteTags(main, "%s", st[tt]),
               minLim = minLim,
               ...);

    if (plotLegend) legendForNumberPlots.fromInfo(
         plotInfo = tsb, 
         legend = c("Down", "Up"),
         pch = 22, pt.cex = 1.5 * legend.cex,
         cex = legend.cex, opacity = 0.5,
         pt.bg = c(col.left, col.right),
         col = c(border.left, border.right));
  }

}


#=========================================================================================================
#
# Plot histogranms of p-values for all tests
#
#=========================================================================================================

GNV.plotPValueHistograms = function(
  analysis,
  analysisName = analysis$analysisName,
  analysisName.pretty = NULL,
  plotDir = "Plots",
  device = "CairoPDF",
  width = 7,
  height = 3,
  prettifyList = NULL,
  onePerPage = TRUE,
  mfrow = c(2,3),
  pPattern = "^p.for",
  breaks = 200,
  cex = 1)
{
  device = match.fun(device);
  device(file = file.path(plotDir,
    spaste("pValueHistograms", prependDash(analysisName),
           ".pdf")),
    width = width, height = height);

  if (is.null(prettifyList)) prettifyList = analysis$prettifyList else 
    prettifyList = lapply( mymapply(c, analysis$prettifyList, prettifyList),
                           unique);

  if (is.null(analysisName.pretty)) analysisName.pretty = prettifyStrings(analysisName.pretty, prettifyList);

  on.exit(dev.off());
  if (!onePerPage) par(mfrow = mfrow);
  par(mar = c(3.2, 3.2, 3, 1));
  par(mgp = c(2, 0.7, 0));
  par(cex = cex);

  res = analysis$combinedDESeqResults$all;
  pValues = res[, grep(pPattern, colnames(res)), drop = FALSE];
  prettyColNames = prettifyStrings(colnames(pValues), prettifyList);
  np = ncol(pValues);
  out = list();
  for (i in 1:np)
    out[[i]] = hist(pValues[, i], breaks = breaks,
         ylab = "Frequency",
         xlab = "p value",
         main = spaste(prettyColNames[i], "\n", analysisName.pretty));

  invisible(out);
}

#=========================================================================================================
#
# Scatterplot of Z for bg vs wt and Z for pert/bg vs bg
#
#=========================================================================================================

GNV.rescuePlotAttributes = function(analysis,
   sigThreshold, # the significant threshold name to be retained
   perturbation,
   wt = "WT",
   bg = "Q175",
   pertBg = spaste(bg, ".", perturbation),
   statName.bg = spaste("Z.for.", bg, ".vs.", wt),
   thresholdType = c("FDR", "p"),
   thresholdValue = if (thresholdType=="FDR") 0.1 else 0.05,

   results = extendedResults(analysis)[[1]],
   sigName.bg = sub("^Z", thresholdType, statName.bg),
   refStat = NULL, # if not given, Z for bg vs. wt wil be used.
   color.exacerbation = "red",
   color.noSigMove = "grey40",
   color.partialRecue = "blue",
   color.nonSignif = "grey70",

   statusPlotOrder = c("No sig. move", "Rescue", "Exacerbation"),
   bg.exacerbation = "pink",
   bg.noSigMove = "grey75",
   bg.nonSignif = "grey85",
   bg.partialRescue = "royalblue",
   rescueScoreCol = "rescueScore",
   sigRescueColBase = "significantRescue",
   sigExacColBase = "significantExacerbation",
   sigRescName = spaste(sigRescueColBase, sigThreshold),
   sigExacName = spaste(sigExacColBase, sigThreshold))
{
  statName = spaste("Z.for.", pertBg, ".vs.", bg)

  thresholdType = match.arg(thresholdType);
  sigName = spaste(thresholdType, ".for.", pertBg, ".vs.", bg)

  rescueStatusNames = c("Exacerbation", "No sig. move", "Rescue");

  res = results;
  stat.bg = res[[statName.bg]]
  nonSignif = replaceMissing ( res[[sigName]] > thresholdValue & res[[sigName.bg]] > thresholdValue, TRUE)
  if (is.null(refStat)) refStat = stat.bg;
  if (is.null(rescueScoreCol))
  {
    sigMove = replaceMissing(res[[sigName]] < thresholdValue & res[[sigName.bg]] < thresholdValue);
    rescueInd = 0 + replaceMissing(res[[statName]] * stat.bg < 0 & sigMove);
    exacInd = 0 + replaceMissing(res[[statName]] * stat.bg > 0 & sigMove);
  } else {
    rescueInd = replaceMissing(res[[sigRescName]]);
    exacInd = replaceMissing(res[[sigExacName]]);
  }
  #stat = res[[rescueScoreCol]];
  rescueStatus = rescueStatusNames[rescueInd - exacInd + 2];
  rescueColors = c(color.noSigMove, color.partialRecue, color.exacerbation)
  rescueBgColors = c(bg.noSigMove, bg.partialRescue, bg.exacerbation)
  names(rescueColors) = names(rescueBgColors) = statusPlotOrder;
  pointColor = rescueColors[rescueStatus];
  pointBg = rescueBgColors[rescueStatus];

  pointColor[nonSignif] = color.nonSignif;
  pointBg[nonSignif] = bg.nonSignif;

  class.num = match(rescueStatus, statusPlotOrder);
  order = order(class.num, refStat)

  list(pointColor = pointColor, pointBg = pointBg, pointOrder = order,
       pointStatus = rescueStatus,
       rescueColors = rescueColors,
       rescueBgColors = rescueBgColors,
       showPoints = class.num > 1)
}

GNV.addRescueLegend = function(plotAttr, ...)
{

  presentTypes = unique(plotAttr$pointStatus);
  keepTypes = names(plotAttr$rescueColors) %in% presentTypes;
  legendClean(..., pch = 21, col = plotAttr$rescueColors[keepTypes],
              pt.bg = plotAttr$rescueBgColors[keepTypes],
              legend = names(plotAttr$rescueColors)[keepTypes]);
}

GNV.rescueScatterplot = function(
  analysis,
  analysisName = analysis$analysisName,
  wt = "WT",
  background = "Q175",
  perturbation,
  pertBack = spaste(background, "/", perturbation),
  pertBack.plot = pertBack,
  levelSep = ".vs.",
  testNameSep = ".for.",
  Zpattern = "Z",
  FDRpattern = "FDR", 
  reverseSign = TRUE,
  simpleYLab = FALSE,

  baselineZCol = NULL,
  baselineFDRCol = NULL,

  limitZTo = Inf,
  limitZTo.x = limitZTo,
  limitZTo.y = limitZTo,

  plotDir = "Plots",
  plotDevice = "CairoPDF",
  plotFileExt = "pdf",
  width = 6, height = 6,
  cex = 1,
  mar = c(3.4, if (reverseSign)4.6 else 3.4, 3, 1),
  mgp = c(2, 0.7, 0),
  main,

  nLabel = 8,
  nConsider = 1000,
  labelDirections = c("0+", "++", "+0", "+-", "0-", "--", "-0", "-+"),
  cex.labels = 0.9,
  forceLabel = NULL,

  cex.lab = 1,
  cex.axis = 1,
  col = 1,
  bg = 0,
  pch = 1,
 
  pt.cex = 0.5,

  addLegend = TRUE,
  legendPos = "auto",
  leg.pt.cex = 2*pt.cex,
  leg.cex = 0.7*cex.lab,
  addCounts = TRUE,
  addOverlapPValues = TRUE,
  stretchLims = 0.12,
  addAblines = TRUE,
  ab.thresholdType = c("FDR", "p"),
  ab.threshold = if (ab.thresholdType=="p") 0.025 else 0.1,
  addAxisLines = FALSE,
  prettifyList = list(character(0), character(0)),
  showPoints = NULL,
  plotAttr,
  rescuePlotType = "Zstatistics",
  excludeIndex = NULL,
  xlim = NULL,
  ylim = NULL,
  ...)
{

  ab.thresholdType = match.arg(ab.thresholdType);
  ab.thresholdType = ab.thresholdType[1];

  if (!is.null(plotDevice))
  {
    device = match.fun(plotDevice);
    device(file = file.path(plotDir,
    spaste("rescueScatterplot", prependDash(rescuePlotType), prependDash(analysisName),
           ".", plotFileExt)),
    width = width, height = height);

    on.exit(dev.off());
  }

  par(mar = mar);
  par(cex = cex);
  par(mgp = mgp);

  data = extendedResults(analysis)[[1]];

  bg.n = make.names(background);
  pb.n = make.names(pertBack);
  wt.n = make.names(wt);
  if (is.null(baselineZCol)) baselineZCol = spaste(Zpattern, testNameSep, bg.n, levelSep, wt.n);
  col1 = baselineZCol;
  x = data[[col1]];
  col2 = spaste(Zpattern, testNameSep, pb.n, levelSep, bg.n);
  if (is.null(baselineFDRCol))
  {
    fdr.x = fdrFromZ(x, alternative = "two.sided");
  } else {
    fdr.x = data[[baselineFDRCol]];
  }
  fdrCol2 = spaste(FDRpattern, testNameSep, pb.n, levelSep, bg.n);
  col2.inverse = spaste(Zpattern, testNameSep, bg.n, levelSep, pb.n);
  col2.legend = spaste(Zpattern, testNameSep, pb.n, levelSep, bg.n);
  col2.inverse.legend = spaste(Zpattern, testNameSep, bg.n, levelSep, pb.n);

  sign = 1-2*reverseSign;
  y = sign * data[[col2]];

  x[replaceMissing(x>limitZTo.x)] = limitZTo.x;
  x[replaceMissing(x< -limitZTo.x)] = -limitZTo.x;

  y[replaceMissing(y>limitZTo.y)] = limitZTo.y;
  y[replaceMissing(y< -limitZTo.y)] = -limitZTo.y;

  x0 = x;
  y0 = y;
  fdr.y = data[[fdrCol2]];

  if (simpleYLab)
  {
    if (reverseSign) ylab = col2.inverse.legend else ylab = col2.legend;
    ylab = prettifyStrings(ylab, mymapply("c", analysis$prettifyList.plots, prettifyList));
  } else
    ylab = spaste( if(reverseSign)
          spaste(prettifyStrings(col2.inverse.legend, mymapply("c", analysis$prettifyList.plots, prettifyList)),
                     "\n(-") else "",
          prettifyStrings(col2.legend, mymapply("c", analysis$prettifyList.plots, prettifyList)),
          if(reverseSign) ")" else "");

  symbol = data$Symbol;
  if (length(excludeIndex)>0)
  {
    if (is.logical(excludeIndex)) excludeIndex = which(excludeIndex);
    if (length(col)==length(x)) col = col[-excludeIndex];
    if (length(pch)==length(x)) pch = pch[-excludeIndex];
    if (length(bg)==length(x)) bg = bg[-excludeIndex];
    if (length(pt.cex)==length(x)) pt.cex = pt.cex[-excludeIndex];
    if (is.logical(showPoints)) showPoints = showPoints[-excludeIndex] else
       showPoints = setdiff(showPoints, excludeIndex);
    x = x[-excludeIndex];
    y = y[-excludeIndex];
    symbol = symbol[-excludeIndex];
  }

  if (is.null(xlim)) xlim = stretchLim(x, rx = addCounts * stretchLims);
  if (is.null(ylim)) ylim = stretchLim(y, rx = addCounts * stretchLims);

  thinnedScatterplot(x, y,
                     xlab = prettifyStrings(col1, mymapply("c", analysis$prettifyList.plots, prettifyList)),
                     ylab = ylab,
                     verbose = TRUE,
                     main = spaste(main, "\n"),
                     cex.lab = cex.lab, cex.axis = cex.axis,
                     col = col, bg = bg, pch = pch, cex = pt.cex,
                     showPoints = showPoints,
                     xlim = xlim, ylim = ylim,
                     ...);

  if (ab.thresholdType=="p") {
    at.x = at.y = -qnorm(ab.threshold/2);
  } else {
    at.x = ablinePositionFromFDR(Z = x0, fdr = fdr.x, threshold = ab.threshold);
    at.y = ablinePositionFromFDR(Z = y0, fdr = fdr.y, threshold = ab.threshold);
  }
  
  if (addAblines) addAblines(at.x = at.x, at.y = at.y);
  if (addAxisLines) addAblines(0);
  if (nLabel > 0)
  {
    labInfo = labelExtremePoints2(x, y, replaceMissing(symbol), cex = cex.labels,
                       ratio.pointToChar = 0.3 * pt.cex, nLabel = nLabel, nConsider = nConsider, forceLabel = forceLabel,
                       directions = labelDirections, scaleForSelection = TRUE);
    avoid2 = fillLabelSpaceWithPoints(labInfo$x, labInfo$y, labInfo$label, cex = cex.labels);
  } else avoid2 = list(x = numeric(0), y = numeric(0));

  if (addLegend)
  {
    leg = GNV.addRescueLegend(plotAttr, legendPos, points.x = c(x, avoid2$x), points.y = c(y, avoid2$y), 
                        pt.cex = leg.pt.cex, cex = leg.cex);
    pfl = pointsFillingLegend(leg);
  } else pfl = list(x = numeric(0), y = numeric(0));

  if (addCounts && addAblines)
  {
    addScatterplotCounts(x = x, y = y, 
                         avoid.x = c(x, pfl$x, avoid2$x), avoid.y = c(y, pfl$y, avoid2$y), 
                         cuts.x = c(-at.x, at.x), cuts.y = c(-at.y, at.y),
                         cex = leg.cex, addFisherTestPValues = addOverlapPValues);
  }
  
}

GNV.rescueScatterplot.FoldChanges = function(
  analysis,
  commonLim = TRUE,
  wt = "WT",
  background = "Q175",
  perturbation,
  pertBack = spaste(background, "/", perturbation),
  levelSep = ".vs.",
  testNameSep = ".for.",
  LFCpattern = "Log2FC",
  addAxisLines = TRUE,
  reverseSign = TRUE,
  rescuePlotType = "LogFoldChanges",
  excludeIndex = NULL,
  baselineFCCol = NULL,
  baselineFDRCol = NULL,
 
   ...)
{
  data = extendedResults(analysis)[[1]];
  bg.n = make.names(background);
  pb.n = make.names(pertBack);
  wt.n = make.names(wt);

  col1 = if (is.null(baselineFCCol)) spaste(LFCpattern, testNameSep, bg.n, levelSep, wt.n) else baselineFCCol;
  col2 = spaste(LFCpattern, testNameSep, pb.n, levelSep, bg.n);
  sign = 1-2*reverseSign;
  x = data[[col1]];
  y = sign * data[[col2]];
  if (length(excludeIndex)>0)
  {
    if (is.logical(excludeIndex)) excludeIndex = which(excludeIndex);
    x = x[-excludeIndex];
    y = y[-excludeIndex];
  }
  if (commonLim)
  {
    xlim = ylim = range(c(x, y), na.rm = TRUE);
  } else {
    xlim = range(x, na.rm = TRUE);
    ylim = range(y, na.rm = TRUE);
  }

  GNV.rescueScatterplot(
    analysis = analysis,
    wt = wt,
    background = background,
    perturbation = perturbation,
    pertBack  = pertBack,
    levelSep =  levelSep,
    testNameSep = testNameSep,
    Zpattern = LFCpattern,
    baselineZCol = baselineFCCol,
    addAblines = FALSE,
    addAxisLines = addAxisLines,
    xlim = xlim, ylim = ylim,
    rescuePlotType = rescuePlotType,
    excludeIndex = excludeIndex,
    ...)
}

  

#===============================================================================================
#
# Rescue manhattan plot
#
#===============================================================================================

GNV.rescueManhattanPlot = function(
  analysis,
  results = extendedResults(analysis)[[1]],
  analysisName = analysis$analysisName,
  analysisName.pretty,
  statName = "rescueScore",
  baselineZCol = "Z.for.Q140.vs.WT.in.allGNV",
  baselineSignifCol = "pBonferroni.for.Q140.vs.WT.in.allGNV",
  sigThreshold = "._FDR.lessThan.0.10._", 
  background = "Q140",
  perturbation,
  plotAttr = NULL,
  refNetworkModules = results$module,
  refStatName = "Z.for.Q140.vs.Q20.in.AllelicSeries",
  refStat = results [[refStatName]],
  sigRescueColBase = "significantRescue",
  sigExacColBase = "significantExacerbation",
  sigRescName = spaste(sigRescueColBase, sigThreshold),
  sigExacName = spaste(sigExacColBase, sigThreshold),
  rescueScoreCol = statName,

  moduleOrder,
  plotIndex = NULL,
  excludeSymbols = NULL,
  excludeIndex = NULL,

  summaryModuleRescue = NULL, 
  summaryModuleRescueP = NULL,

  plotDir = "Plots",
  plotDevice = "CairoPDF",
  width = 9, height = 6,
  setLayout = TRUE,
  cex = 1,
  mar = c(0.5, 3.3, 3, 0.5),
  mgp = c(2, 0.7, 0),

  summaryRescueHeight = 0.5,

  plotDataDir = file.path(plotDir, "PlotData"),
  main = spaste(analysisName.pretty, "\n"), cex.main = 1,

  summaryModuleRescue.main = "",
  summaryModuleRescue.cex.main = cex.main,
  summaryModuleRescue.ylab = "",

  nLabel = 8,
  nConsider = 1000,
  cex.labels = 0.9,

  cex.lab = 1,
  cex.axis = 1,
  col = 1,
  bg = 0,
  pch = 1,
  pt.cex = 0.5,

  addLegend = TRUE,
  leg.pt.cex = 1.4,
  leg.cex = 0.8,

  ...)
{
  
  if (is.null(plotAttr))
    plotAttr = GNV.rescuePlotAttributes(
      results = results,
      bg = background,
      sigThreshold = sigThreshold,
      perturbation = perturbation,
      statName.bg = baselineZCol,
      sigName.bg = baselineSignifCol, 
      sigRescueColBase = sigRescueColBase,
      sigExacColBase = sigExacColBase,
      sigRescName = sigRescName,
      sigExacName = sigExacName,
      rescueScoreCol = rescueScoreCol);

  res = results;
  stat = res[[statName]]

  if (!is.null(excludeSymbols))
  {
    excludeIndex = match(excludeSymbols, res$Symbol);
    if (any(is.na(excludeIndex)))
    {
      warning("Could not find the following entries in 'excludeSymbols' in 'res$Symbol':\n  ",
               paste(excludeSymbols[is.na(excludeIndex)], collapse = ", "));
      excludeIndex = excludeIndex[is.finite(excludeIndex)];
    }
  }

  if (!is.null(plotDevice))
  {
    plotDevice = match.fun(plotDevice);
    plotDevice(spaste(plotDir, "/rescueDotPlots", sigThreshold, analysisName, ".pdf"), 
                      wi = width, he = height);
    on.exit(dev.off());
  }
  
  nRows = 1 + addLegend + !is.null(summaryModuleRescue);
  heights = c(1, list(numeric(0), summaryRescueHeight)[[1+!is.null(summaryModuleRescue)]],
                 list(numeric(0), 0.14)[[addLegend + 1]])
  if (setLayout) 
     layout(matrix(1:nRows, nRows, 1), heights = heights);

  par(mar = mar); par(mgp = mgp); par(cex = cex);
  plotData = rescuePlot.base(
       plotStat = stat,
       rescueStatus = plotAttr$pointStatus,
       rescueColors = plotAttr$rescueColors,
       rescueBgColors = plotAttr$rescueBgColors,
       rescueLevels = names(plotAttr$rescueColors),
       statusPlotOrder = names(plotAttr$rescueColors),
       class = refNetworkModules,
       referenceStat = refStat,
       plotIndex = plotIndex,
       excludeIndex = excludeIndex,
       classOrder = moduleOrder,
       classPValue = if (is.null(summaryModuleRescue)) summaryModuleRescueP else NULL,
       separator.col = "grey",
       separator.lwd = 2,
       pt.cex = 1,
       pch = 21,
       classPrefix = "M",
       classSuffix = "",
       classLabels.cex = 1,
       classLabels.col = "black",
       ylab = "Rescue score",
       main = main, cex.main = cex.main,
       nodeNames = res$Symbol, plotLegend = FALSE,
       nConsider = nConsider, nLabel = nLabel, ...);
  if (!is.null(summaryModuleRescue))
  {
    x = plotData$moduleLabelData$x;
    if (length(x)!=length(summaryModuleRescue))
      stop("Length of 'summaryModuleRescue' must equal the number of modules.");

    lim = range(summaryModuleRescue);
    if (lim[1] > 0) lim[1] = 0;
    if (lim[2] < 0) lim[2] = 0;
    lim[2] = lim[2] + 0.1 * (lim[2] - lim[1]);

    col = numbers2colors(summaryModuleRescue, color = rev(blueWhiteRed(100)), signed = TRUE);
    n0 = table(refNetworkModules [ refNetworkModules %in% moduleOrder ]);
    nGenesInModules = n0[ match(moduleOrder, names(n0))];
    width = 30 + nGenesInModules/5;
    bottom = ifelse(summaryModuleRescue > 0, 0, summaryModuleRescue);
    top = ifelse(summaryModuleRescue > 0, summaryModuleRescue, 0);
    par(mar = c(0, mar[2], if (summaryModuleRescue.main=="") 0 else 1.2, mar[4]));
    plot(plotData$pointData$x, plotData$pointData$x, ylim = lim,
         type = "n", xaxt = "none", frame = TRUE, main = summaryModuleRescue.main, xlab = "", 
         ylab = summaryModuleRescue.ylab,
         cex.main = summaryModuleRescue.cex.main);
    lapply(plotData$verticalLineData$x, function(x) 
           lines(c(x,x), par("usr")[3:4], col = plotData$verticalLineData$lineColor[1],
                 lwd = 2));
    abline(h=0, col = "grey20");
    rect(xleft = x - width/2, ybottom = bottom, ytop = top, xright = x+width/2,
         border = "grey20", col = col);

    if (!is.null(summaryModuleRescueP))
      text(x = x, y = top, labels = pValueStars2(summaryModuleRescueP),
           adj = c(0.5, 0)); 
  }
    
  if (addLegend)
  {
    par(mar = c(mar[1], mar[2], 0.1, mar[4]));
    plot(c(0,1), type = "n", axes = FALSE, xlab = "", ylab = "");
    GNV.addRescueLegend(plotAttr, "top", pt.cex = leg.pt.cex, ncol = 4, cex = 0.8);
  }

  if (!is.null(plotDataDir))
  {
    dir.create(plotDataDir, recursive = TRUE, showWarnings = FALSE);
    mymapply(function(data, name)
         write.csv(data, file = gzfile(spaste(plotDataDir, "/rescueDotPlots",
               sub("^X", "", make.names(sigThreshold)), analysisName, "-", name, ".csv.gz")),
               quote = TRUE, row.names = FALSE),
        plotData, names(plotData));
  }
  invisible(plotData);
}

#=========================================================================================================
#
# module rescue statistics
#
#=========================================================================================================

moduleRescueStatistics = function(
  rescueScore = NULL,
  #rescueFraction,
  moduleLabels,
  moduleMembershipZ,
  stat.perturbation,
  stat.reference,
  useModules = NULL,
  ignoreModules = 0,
  get = c("mean", "weighted mean", "correlation of RS and MM"),
  simplify = TRUE
)
{
  options.get = c("mean", "weighted mean", "correlation of RS and MM");
  get = options.get[ charmatch(get, options.get)];
  if (any(is.na(get))) stop("'get' must be one or more of ", paste(options.get, collapse = ", "));
 
  if (is.null(rescueScore)) rescueScore = rescueScoreFnc(stat.perturbation, stat.reference, inverse = TRUE);

  if (!is.null(useModules)) 
     moduleLabels[!replaceMissing(moduleLabels, replaceWith = -1000) %in% useModules] = NA;
  if (!is.null(ignoreModules)) 
    moduleLabels[ replaceMissing(moduleLabels, replaceWith = -1000) %in% ignoreModules] = NA;
  modLevels = sort(unique(moduleLabels[!is.na(moduleLabels)]));
  mean = tapply(rescueScore, moduleLabels, mean, na.rm = TRUE);
  weightedMean = sapply(modLevels, function(ml)
  {
    weighted.mean(x = rescueScore[replaceMissing(moduleLabels==ml)], 
                  w = moduleMembershipZ[replaceMissing(moduleLabels==ml)],
                  na.rm = TRUE)
  });
  cor.RS.MM = sapply(modLevels, function(ml)
  {
    stats::cor(rescueScore[replaceMissing(moduleLabels==ml)], 
               moduleMembershipZ[replaceMissing(moduleLabels==ml)],
               use = 'p');
  });
  names(weightedMean) = names(cor.RS.MM) = modLevels;
  order = match(useModules, modLevels);
  stats = list(mean[order], weightedMean[order], cor.RS.MM[order]);
  names(stats) = options.get;
  stats = stats[get]
  if (length(get)==1 && simplify) stats = stats[[1]];
  stats;
}


#=========================================================================================================
#
# Scatterplot of Z for bg vs wt and Z for pert/bg vs bg
#
#=========================================================================================================

## Default addAblines assumes that the externalData also represents Z statistics.
GNV.compareWithExternalData = function(
  analysis,
  analysisName = analysis$analysisName,

  externalDataCol = NULL,
  externalData,

  wt = "WT",
  background = "Q175",
  #perturbation,
  #pertBack = spaste(background, "/", perturbation),
  levelSep = ".vs.",
  testNameSep = ".for.",

  Zpattern = "Z",
  plotDir = "Plots",
  plotDevice = "CairoPDF",
  width = 6, height = 6,
  cex = 1,
  mar = c(3.4, 3.4, 4, 1),
  mgp = c(2, 0.7, 0),

  main,
  externalLabel = NULL,
  externalFileID,

  nLabel = 8,
  nConsider = 1000,
  labelDirections = c("0+", "++", "+0", "+-", "0-", "--", "-0", "-+"),
  cex.labels = 0.9,

  cex.lab = 1,
  cex.axis = 1,
  col = 1,
  bg = 0,
  pch = 1,
 
  pt.cex = 0.5,

  addLegend = TRUE,
  leg.pt.cex = 2.4*pt.cex,
  leg.cex = 0.7*cex.lab,
  
  addCounts = TRUE,
  addOverlapPValues = FALSE,
  cornerOnly = addOverlapPValues,
  stretchLims = 0.12,

  addAblines = TRUE,
  ab.thresholdType = c("FDR", "p"),
  ab.threshold = if (ab.thresholdType=="p") 0.025 else 0.1,
  prettifyList = NULL,

  savePlotData = TRUE,
  plotDataDir = file.path(plotDir, "PlotData"),
  showPoints = NULL,
  plotAttr,
  xlim = NULL,
  ylim = NULL,
  ...)
{

  if (!is.null(plotDevice))
  {
    device = match.fun(plotDevice);
    device(file = file.path(plotDir,
        spaste("Scaterplot", prependDash(analysisName), "-vs-", externalFileID, ".pdf")),
    width = width, height = height);

    on.exit(dev.off());
  }
  par(mar = mar);
  par(cex = cex);
  par(mgp = mgp);

  bg.n = make.names(background);
  wt.n = make.names(wt);

  data = extendedResults(analysis)[[1]];
  col1 = spaste(Zpattern, testNameSep, bg.n, levelSep, wt.n);

  x = data[[col1]];
  y = if (is.null(externalDataCol)) externalData else data[[externalDataCol]];

  if (is.null(externalLabel))
  {
    if (is.null(externalDataCol)) stop("'externalLabel' must be supplied when 'externalDataCol' is not.");
    externalLabel = externalDataCol;
  }

  if (is.null(xlim)) xlim = stretchLim(x, rx = addCounts * stretchLims);
  if (is.null(ylim)) ylim = stretchLim(y, rx = addCounts * stretchLims);

  thinnedScatterplot(x, y,
                     xlab = prettifyStrings.multi(col1, list(analysis$prettifyList.plots, prettifyList)),
                     ylab = prettifyStrings.multi(externalLabel, list(analysis$prettifyList.plots, prettifyList)),
                     verbose = TRUE,
                     main = spaste(prettifyStrings.multi(main, list(analysis$prettifyList.plots, prettifyList)), 
                                   "\n"),
                     cex.lab = cex.lab, cex.axis = cex.axis,
                     col = col, bg = bg, pch = pch, cex = pt.cex,
                     showPoints = showPoints, xlim = xlim, ylim = ylim,
                     ...);

  
  if (addAblines) 
  {
    ab.thresholdType = match.arg(ab.thresholdType);
    if (ab.thresholdType=="p") {
      at.x = at.y = -qnorm(ab.threshold/2);
    } else {
      at.x = ablinePositionFromFDR(Z = x, threshold = ab.threshold);
      at.y = ablinePositionFromFDR(Z = y, threshold = ab.threshold);
    }
    addAblines(at.x = at.x, at.y = at.y);
  }
  if (nLabel > 0)
  {
    labInfo = labelExtremePoints2(x, y, replaceMissing(data$Symbol), cex = cex.labels,
                       ratio.pointToChar = 0.3 * pt.cex, nLabel = nLabel, nConsider = nConsider,
                       directions = labelDirections, scaleForSelection = TRUE);
    avoid2 = fillLabelSpaceWithPoints(labInfo$x, labInfo$y, labInfo$label, cex = cex.labels);
  } else avoid2 = list(x = numeric(0), y = numeric(0));
  if (addLegend)
  {
    leg = GNV.addRescueLegend(plotAttr, "auto", points.x = c(x, avoid2$x), points.y = c(y, avoid2$y), 
                        pt.cex = leg.pt.cex, cex = leg.cex);
    pfl = pointsFillingLegend(leg);
  } else pfl = list(x = numeric(0), y = numeric(0));

  if (addCounts && addAblines)
  {
    addScatterplotCounts(x = x, y = y, 
                         avoid.x = c(x, pfl$x, avoid2$x), avoid.y = c(y, pfl$y, avoid2$y), 
                         cuts.x = c(-at.x, at.x), cuts.y = c(-at.y, at.y),
                         cex = leg.cex, addFisherTestPValues = addOverlapPValues,
                         cornerOnly = cornerOnly)
  }
  if (savePlotData)
  {
    suppressWarnings(dir.create(plotDataDir, recursive = TRUE, showWarnings = FALSE));
    n = length(x);
    plotData = data.frame(x = x, y = y, Symbol = data$Symbol, BorderColor = RColors2hex(checkOrExtend(col, n, "col")),
                          FillColors = RColors2hex(checkOrExtend(bg, n, "bg")));

    write.csv(plotData, file = gzfile(file.path(plotDataDir, 
                            spaste("Scaterplot", prependDash(analysisName), "-vs-", externalFileID, ".csv.gz"))),
              quote = TRUE, row.names = FALSE);
  }
}


#=========================================================================================================
#
# Correlation heatmap of Z statistics
#
#=========================================================================================================

GNV.ZcorrelationHeatmap = function(
  analysis,
  analysisName = analysis$analysisName,
  testNameSep = ".for.",
  Zpattern = spaste("Z",testNameSep),

  plotDir = "Plots",
  plotDevice = "CairoPDF",
  width = 8, height = 4,
  cex.text = 1,
  mar.main = c(6, 8.5, 2, 1.5),
  main,
  prettifyList = NULL,
  ...)
{
  if (!is.null(plotDevice))
  {
    device = match.fun(plotDevice);
    device(file = file.path(plotDir,
    spaste("heatmapOfZCorrelations", prependDash(analysisName),
           ".pdf")),
    width = width, height = height);

    on.exit(dev.off());
  }

  data = extendedResults(analysis)[[1]];
  cols = grep(Zpattern, colnames(data));
  Zs = data[, cols];
  labels1 = sub(Zpattern, "", colnames(Zs));
  labels = prettifyStrings.multi(labels1, list(analysis$prettifyList.plots, prettifyList));
  cor.Zs = cor(Zs, use = 'p');
  corHeatmapWithDendro(cor.Zs, mar.main = mar.main,
         xLabels = labels, yLabels = labels, 
         cex.text = cex.text, 
         main = prettifyStrings.multi(main, list(analysis$prettifyList.plots, prettifyList)))
}  

#=========================================================================================================
#
# Main wrapper
#
#=========================================================================================================

GNV.standardPlots = function(
   analysis,

   analysisName = analysis$analysisName,
   analysisName.pretty = analysis$analysisName,
   wt = "WT",
   background = "Q175",
   perturbation,
   pertBack = spaste(background, "/", perturbation),

   prettifyList.plots = NULL,
   prettifyList.legends = prettifyList.plots,

   levelSep = ".vs.",
   testNameSep = ".for.",
   andSep = ".and.",
   pattern.Z = spaste("Z", testNameSep),

   plotDir = "Plots",
   plotDevice = "CairoPDF",

   nSignif.width = 7,
   nSignif.height = 3,
   nSignif.mar = c(0.5, 8, 2, 0.5),
   nSignif.addThresholdToMain = TRUE,
   nSignif.thresholdSep = "\n",

   nSignif.cex = 1,
   nSignif.main = spaste("Numbers of significantly associated genes, ", analysisName),
   nSignif.col.right = "#FF9070", nSignif.col.left = "skyblue",
   nSignif.border.right = "brown", nSignif.border.left = "blue",
   nSignif.barGap = 0.4,
   nSignif.plotLegend = TRUE,
   nSignif.legendPos = "bottomright",
   nSignif.legend.cex = 1,

   nMoved.main = spaste("Numbers of rescued or exacerbated genes, ", analysisName, "\n",
                 "split by DE direction in ", prettifyStrings(spaste(background, testNameSep, wt),
                                                              analysis$prettifyList.plots)),
   nMoved.maxCharPerLine = 34,
   nMoved.width = 7,
   nMoved.height = 5,
   nMoved.mar = c(0.5, 10, 2, 0.5),
   nMoved.splitByThreshold = TRUE,
   nMoved.addThresholdToMain = TRUE,
   nMoved.thresholdSep = "\n",
   nMoved.plotLegend = FALSE,
   nMoved.legendPos = "bottomright",
   nMoved.legend.cex = 1,

   colorSignificanceThreshold = analysis$topThresholdNames[1],

   plotAttr = NULL,
   sigRescueColBase = "significantRescue",
   sigExacColBase = "significantExacerbation",
   sigRescName = spaste(sigRescueColBase, colorSignificanceThreshold),
   sigExacName = spaste(sigExacColBase,  colorSignificanceThreshold),
   rescueScoreCol = "rescueScore",

   # Recue scatterplot
   rsc.baselineZCol = NULL,
   rsc.baselineFDRCol = NULL,
   rsc.baselineFCCol = NULL,
   rsc.reverseSign = TRUE,
   rsc.width = 6,
   rsc.height = 6,
   rsc.cex = 1, 
   rsc.mar = c(3.4, if (rsc.reverseSign)4.6 else 3.4, 3, 1),
   rsc.mgp = c(2, 0.7, 0),
   rsc.main = "Rescue scatterplot",
   rsc.nLabel = 8, 
   rsc.labelDirections = c("0+", "++", "+0", "+-", "0-", "--", "-0", "-+"),
   rsc.nConsider = 1000,
   rsc.cex.labels = 0.8,
   rsc.cex.lab = 1,
   rsc.cex.axis = 1,
   rsc.col = NULL,
   rsc.bg = NULL,
   rsc.pch = 1,
   rsc.pt.cex = 0.5,
   rsc.ylim = NULL,
   rsc.pointsThatMustBeShown = NULL,
   rsc.addLegend = TRUE,
   rsc.fc.commonLim = TRUE,
   rsc.excludeIndex = NULL,

   rsc.addAblines = TRUE,
   rsc.ab.thresholdType = "FDR",
   rsc.ab.threshold = 0.1,
   rsc.addAxisLines = FALSE,
   rsc.prettifyList = list(character(0), character(0)),
   rsc.rescuePlotType = "Zstatistics",
   rsc.addCounts = TRUE,
   rsc.addOverlapPValues = FALSE,
   rsc.stretchLims = 0.12,

   # Scatterplot against external data
   # most arguments will be copied from the rescue scatterplot
   # Data and description arguments can be lists/vectors; a separate plot will be shown for each component
   eds.externalDataCol = NULL,
   eds.externalData = NULL,
   eds.externalLabel,
   eds.externalFileID,
   eds.main,
   eds.mar = c(3.4, 3.4, 4, 1),
   eds.mgp = c(2, 0.7, 0),
   eds.saveData = TRUE,
   eds.plotDataDir = file.path(plotDir, "PlotData"),
   eds.col = rsc.col,
   eds.bg = rsc.bg,
   eds.pch = rsc.pch,
   eds.pointsThatMustBeShown = NULL,
   eds.addLegend = TRUE,
   eds.addCounts = TRUE,
   eds.addOverlapPValues = FALSE,
   eds.cornerOnly = eds.addOverlapPValues,
   eds.stretchLims = 0.12,
   eds.labelDirections = c("0+", "++", "+0", "+-", "0-", "--", "-0", "-+"),

   # Z correlation heatmap
   zch.mar.main = c(6, 8.5, 2, 1),
   zch.main = spaste("Correlations of DE significance, ", analysisName.pretty, 
                     "                         "),
   zch.cex.text = 1,
   zch.width = 8, zch.height = 4,
   ...
   )
{
  baseTestName.pretty = prettifyStrings(spaste(background, levelSep, wt), analysis$prettifyList.plots);

  dir.create(plotDir, recursive = TRUE, showWarnings = FALSE);

  GNV.ZcorrelationHeatmap(
      analysis,
      analysisName = analysis$analysisName,
      Zpattern = pattern.Z,

      plotDir = plotDir,
      plotDevice = plotDevice,
      width = zch.width, height = zch.height,
      cex.text = zch.cex.text,
      mar.main = zch.mar.main,
      main = zch.main,
      prettifyList = prettifyList.plots,
      ...)

  if (is.null(plotAttr))
    plotAttr = GNV.rescuePlotAttributes(analysis,
      sigThreshold = colorSignificanceThreshold,
      perturbation = make.names(perturbation),
      wt = make.names(wt),
      bg = make.names(background),
      pertBg = make.names(pertBack),
      statName.bg = rsc.baselineZCol,
      sigName.bg = rsc.baselineFDRCol,
      sigRescueColBase = sigRescueColBase,
      sigExacColBase = sigExacColBase,
      sigRescName = sigRescName,
      sigExacName = sigExacName,
      rescueScoreCol = rescueScoreCol);

  if (is.null(rsc.col)) 
  {
    rsc.col = plotAttr$pointColor;
    rsc.col2 = plotAttr$rescueColors[1];
    plotNoColor = TRUE;
  } else {
    plotNoColor = FALSE;
    rsc.col2 = rsc.col;
  }
  if (is.null(rsc.bg)) 
  {
    rsc.bg = plotAttr$pointBg;
    rsc.bg2 = plotAttr$rescueBgColors[1];
    plotNoColor = TRUE;
  } else {
    rsc.bg2 = rsc.bg;
    plotNoColor = FALSE;
  }
  if (is.null(rsc.pointsThatMustBeShown)) rsc.pointsThatMustBeShown = plotAttr$showPoints;

  if (is.null(eds.col)) eds.col = plotAttr$pointColor;
  if (is.null(eds.bg)) eds.bg = plotAttr$pointBg;
  if (is.null(eds.pointsThatMustBeShown)) eds.pointsThatMustBeShown = plotAttr$showPoints;

  if (!is.null(eds.externalData) || !is.null(eds.externalDataCol))
  {
     if (is.null(eds.externalData))
     {
       nPlots = 1;
       eds.externalData = list(NULL);
       if (is.null(eds.externalLabel)) eds.externalLabel = eds.externalDataCol;
     } else {
       if (is.atomic(eds.externalData)) eds.externalData = list(eds.externalData);
       nPlots = length(eds.externalData);
     }
     eds.main = checkOrExtend(eds.main, nPlots, "eds.main");
     eds.externalLabel = checkOrExtend(eds.externalLabel, nPlots);
     if (length(eds.externalFileID)!=nPlots) 
        stop("Lenghts of 'eds.externalData' and 'eds.externalFileID' must equal.");
     for (pl in 1:nPlots)
       GNV.compareWithExternalData(
          analysis, analysisName = analysis$analysisName,
          externalData = eds.externalData[[pl]],
          wt = wt, background = background, #perturbation = perturbation, pertBack = pertBack,
          levelSep = ".vs.", testNameSep = ".for.", Zpattern = "Z",
          plotDir = plotDir, plotDevice = plotDevice,

          width = rsc.width, height = rsc.height,
          cex = rsc.cex,
          mar = eds.mar,
          mgp = eds.mgp,
          main = eds.main[pl],
          nLabel = rsc.nLabel,
          nConsider = rsc.nConsider,
          labelDirections = eds.labelDirections,
          cex.labels = rsc.cex.labels,
          cex.lab = rsc.cex.lab,
          cex.axis = rsc.cex.axis,
          col = eds.col,
          bg = eds.bg,
          pch = rsc.pch,
          pt.cex = rsc.pt.cex,

          addLegend = eds.addLegend,

          externalLabel = eds.externalLabel[pl],
          externalFileID = eds.externalFileID[pl],
          addAblines = TRUE,
          prettifyList = prettifyList.plots,
          savePlotData = eds.saveData,
          plotDataDir = eds.plotDataDir,
          showPoints = eds.pointsThatMustBeShown,
          plotAttr = plotAttr,

          addCounts = eds.addCounts,
          addOverlapPValues = eds.addOverlapPValues,
          cornerOnly = eds.cornerOnly,
          stretchLims = eds.stretchLims,
          ...)
  }  


  GNV.rescueScatterplot(
     analysis, analysisName = analysis$analysisName,
     wt = wt, background = background, perturbation = perturbation, pertBack = pertBack,
     levelSep = ".vs.", testNameSep = ".for.", Zpattern = "Z",
     baselineZCol = rsc.baselineZCol, baselineFDRCol = rsc.baselineFDRCol,
     reverseSign = rsc.reverseSign,
     plotDir = plotDir, plotDevice = plotDevice,
     width = rsc.width, height = rsc.height,
     cex = rsc.cex,
     mar = rsc.mar,
     mgp = rsc.mgp,
     main = rsc.main,
     nLabel = rsc.nLabel,
     nConsider = rsc.nConsider,
     labelDirections = rsc.labelDirections,
     cex.labels = rsc.cex.labels,
     cex.lab = rsc.cex.lab,
     cex.axis = rsc.cex.axis,
     col = rsc.col,
     bg = rsc.bg,
     pch = rsc.pch,
     pt.cex = rsc.pt.cex,
     ylim = rsc.ylim,
     addLegend = rsc.addLegend,

     showPoints = rsc.pointsThatMustBeShown,
     plotAttr = plotAttr,
     excludeIndex = rsc.excludeIndex,
     addAblines = rsc.addAblines,
     ab.thresholdType = rsc.ab.thresholdType,
     ab.threshold = rsc.ab.threshold,
     addAxisLines = rsc.addAxisLines,
     prettifyList = rsc.prettifyList,
     rescuePlotType = rsc.rescuePlotType,
     addCounts = rsc.addCounts,
     addOverlapPValues = rsc.addOverlapPValues,
     stretchLims = rsc.stretchLims,
 
     ...)

  GNV.rescueScatterplot.FoldChanges(analysis, analysisName = analysis$analysisName,
     wt = wt, background = background, perturbation = perturbation, pertBack = pertBack,
     LFCpattern = "log2FC",
     baselineFCCol = rsc.baselineFCCol,
     reverseSign = rsc.reverseSign,
     commonLim = rsc.fc.commonLim,
     plotDir = plotDir, plotDevice = plotDevice,

     width = rsc.width, height = rsc.height,
     cex = rsc.cex,
     mar = rsc.mar,
     mgp = rsc.mgp,
     main = rsc.main,
     nLabel = rsc.nLabel,
     nConsider = rsc.nConsider,
     labelDirections = rsc.labelDirections,
     cex.labels = rsc.cex.labels,
     cex.lab = rsc.cex.lab,
     cex.axis = rsc.cex.axis,
     col = rsc.col,
     bg = rsc.bg,
     pch = rsc.pch,
     pt.cex = rsc.pt.cex,
     addLegend = rsc.addLegend,

     addAxisLines = TRUE,
     prettifyList = prettifyList.plots,
     showPoints = rsc.pointsThatMustBeShown,
     plotAttr = plotAttr,
     excludeIndex = rsc.excludeIndex,
     rescuePlotType = "logFoldChanges",
     addCounts = rsc.addCounts,
     addOverlapPValues = rsc.addOverlapPValues,
     stretchLims = rsc.stretchLims,
     ...)

  if (plotNoColor) GNV.rescueScatterplot.FoldChanges(analysis, analysisName = analysis$analysisName,
     wt = wt, background = background, perturbation = perturbation, pertBack = pertBack,
     LFCpattern = "log2FC",
     baselineFCCol = rsc.baselineFCCol,
     reverseSign = rsc.reverseSign,
     commonLim = rsc.fc.commonLim,
     plotDir = plotDir, plotDevice = plotDevice,

     width = rsc.width, height = rsc.height,
     cex = rsc.cex,
     mar = rsc.mar,
     mgp = rsc.mgp,
     main = rsc.main,
     nLabel = rsc.nLabel,
     nConsider = rsc.nConsider,
     labelDirections = rsc.labelDirections,
     cex.labels = rsc.cex.labels,
     cex.lab = rsc.cex.lab,
     cex.axis = rsc.cex.axis,
     col = rsc.col2,
     bg = rsc.bg2,
     pch = rsc.pch,
     pt.cex = rsc.pt.cex,
     addLegend = FALSE,

     addAxisLines = TRUE,
     prettifyList = prettifyList.plots,
     showPoints = rsc.pointsThatMustBeShown,
     plotAttr = plotAttr,
     excludeIndex = rsc.excludeIndex,
     rescuePlotType = "logFoldChanges-noColor",
     addCounts = rsc.addCounts,
     addOverlapPValues = rsc.addOverlapPValues,
     stretchLims = rsc.stretchLims,
     ...)




  # Numbers of sig
  GNV.plotNumbersOfSignifGenes(
   analysis,
   analysisName = analysisName,
   plotDir = plotDir,
   device = plotDevice,
   width = nSignif.width,
   height = nSignif.height,

   mar = nSignif.mar,
   cex = nSignif.cex,
   main = nSignif.main,
   col.right = nSignif.col.right, col.left = nSignif.col.left,
   border.right = nSignif.border.right, border.left = nSignif.border.left,
   barGap = nSignif.barGap,
   plotLegend = nSignif.plotLegend,
   legendPos = nSignif.legendPos,
   legend.cex = nSignif.legend.cex,
   prettifyList = prettifyList.legends,
   ...)

  GNV.plotNumbersOfMovedGenes(
   analysis,
   analysisName = analysis$analysisName,
   plotDir = plotDir,
   device = plotDevice,
   width = nMoved.width,
   height = nMoved.height,

   mar = nMoved.mar,
   cex = nSignif.cex,
   main = nMoved.main,
   col.right = nSignif.col.right, col.left = nSignif.col.left,
   border.right = nSignif.border.right, border.left = nSignif.border.left,
   barGap = nSignif.barGap,
   plotLegend = nMoved.plotLegend,
   legendPos = nMoved.legendPos,
   legend.cex = nMoved.legend.cex,
   maxCharPerLine = nMoved.maxCharPerLine,
   prettifyList = prettifyList.plots,
   ...)

}


#========================================================================================================
#
# Enrichment plot functions
#
#========================================================================================================

enrichmentPlot.matrix = function(
   analysis,

   analysisName = analysis$analysisName,
   analysisName.pretty = analysis$analysisName,
   prettifyList = NULL,

   plotDir = "Plots",
   plotDevice = "CairoPDF",
   plotFileID = "",

   width = NULL,
   height = NULL,
   baseWidth = 6.6,
   widthPerCol = 0.3,
   baseHeight = 2.2,
   heightPerRow = 0.20,
   mar = c(6, 16, 2, 0.5),

   WGCNACollection, 
   collectionName = "WGCNA",

   plotComponent = "core GNV",

   # useModules can be a vector of module labels or a p value theshold.
   useModules = NULL,
   maxCharPerLine = 35,

   maxColsPerPage = 20,

   rowPThreshold = 1,
   colPThreshold = 1,

   dropRowPattern = NULL,

   ...

#   wt = "WT",
#   background = "Q175",
#   perturbation,
#   pertBack = spaste(background, "/", perturbation),
#
#   levelSep = ".vs.",
#   testNameSep = ".for.",
#   andSep = ".and.",
#   pattern.Z = spaste("Z", testNameSep),
)
{

  enrComp = make.names(plotComponent);
  counts1 = t(analysis$enrichment[[ enrComp ]] [[collectionName]]$countsInDataSet)
  stopifnot(isTRUE(all.equal(dataSetIDs(WGCNACollection), as.character( colnames(counts1)))))

  p1 = t(analysis$enrichment[[ enrComp ]] [[collectionName]]$pValues);

  keepRows = rowSums(p1<rowPThreshold) > 0;
  keepCols = colSums(p1 < colPThreshold) > 0;

  if (!any(keepRows)  || !any(keepCols))
  {
    printFlush("No rows or columns remained after restriction. Returning without generating a plot.");
    return(NULL);
  }

  p1 = p1[keepRows, keepCols, drop = FALSE];
  counts1 = counts1[keepRows, keepCols, drop = FALSE];

  if (!is.null(dropRowPattern))
  {
    keepRows = multiGrep(dropRowPattern, rownames(counts1), invert = TRUE);
    if (length(keepRows)==0)
    {
      printFlush("No rows or columns remained after restriction. Returning without generating a plot.");
      return(NULL);
    }
    p1 = p1[keepRows, , drop = FALSE];
    counts1 = counts1[keepRows, , drop = FALSE];
  }

  mat = -log10(p1);
  mat[mat > 20] = 20;


  setNames0 = dataSetNames(WGCNACollection)[keepCols];
  setNames = spaste(sub(").+", "", setNames0), ")")
  moduleLabels = as.numeric(multiSub(c("M\\.*", " .+"), c("", ""), setNames));
  if (any(is.na(moduleLabels)))
  {
    setNames = sub(":.+", "", setNames0);
    moduleLabels = as.numeric(multiSub(c("M\\.*", " .+"), c("", ""), setNames));
  }

  if (is.null(useModules)) useModules = moduleLabels;

  if (length(useModules)==1 && is.numeric(useModules) && useModules < 1)
      useModules = moduleLabels[colSums(p1 < useModules) > 0];

  keepCols2 = moduleLabels %in% useModules;
  if (!is.null(plotDevice))
  {
    if (is.null(height)) height = baseHeight + heightPerRow*nrow(p1);
    if (is.null(width)) width = baseWidth + widthPerCol*min(maxColsPerPage, sum(keepCols2));
    device = match.fun(plotDevice);
    device(file = file.path(plotDir,
                     spaste("enrichmentInOlderWGCNAModules-Matrix", prependDash(analysisName), 
                            prependDash(plotFileID), ".pdf")),
           width = width, height = height);
    on.exit(dev.off());
  }

  rowOrder = c(grep("Downregulated|down.for", rownames(counts1)), grep("Upregulated|up.for", rownames(counts1)));
  yLabels = formatLabels(
              prettifyStrings.multi(rownames(counts1), list(analysis$prettifyList.plots, prettifyList)),
                              maxCharPerLine = maxCharPerLine);

  par(mar = mar)

  capitalize = function(s) { substr(s,1,1) = toupper(substr(s, 1,1)); s }
  labeledHeatmap.multiPage(mat[rowOrder, keepCols2, drop = FALSE],
                 yLabels = sapply(yLabels[rowOrder], capitalize),
                 xLabels = spaste("ME", labels2colors(moduleLabels[keepCols2])),
                 xSymbols = setNames[keepCols2],
                 textMatrix = counts1[rowOrder, keepCols2, drop = FALSE],
                 colors = blueWhiteRed(100)[50:100],
                 zlim = c(0, 30),
                 setStdMargins = FALSE, maxColsPerPage = maxColsPerPage, ...)
}
 
#=========================================================================================================
#
# plotGeneExpressionByGroup 
#
#=========================================================================================================

GVN.plotGeneExpressionByGroup = function(expr,
                                index = match(entrez, colnames(expr)),
                                entrez = geneAnnotation$Entrez[ match(symbol, geneAnnotation$Symbol)],
                                geneAnnotation = geneAnnotationFromEntrez(colnames(expr), organism = "mouse"),
                                symbol,
                                genotype,
                                wt = "WT",
                                background,
                                perturbation,
                                prettifyList,
                                genotypeOrder = c(wt, background, perturbation,
                                               spaste(background, "/", perturbation)),
                                col = "white", border = 1, lines.col = "grey", cex.text = 1,
                                ylab = spaste(symbol, " normalized expression"),
                                ...)
{
  data = expr[, index];
  groupBarplot(data= data, group = prettifyStrings(genotype, prettifyList), 
    groupOrder = prettifyStrings(genotypeOrder, prettifyList),
    col = col, border = border,
    lines.col = lines.col, cex.text = cex.text,
    ylab = ylab, 
    ...)
}


GVN.plotGeneExpressionBoxplot = function(expr,
                                index = match(entrez, colnames(expr)),
                                entrez = geneAnnotation$Entrez[ match(symbol, geneAnnotation$Symbol)],
                                geneAnnotation = geneAnnotationFromEntrez(colnames(expr), organism = "mouse"),
                                symbol,
                                genotype,
                                wt = "WT",
                                background,
                                perturbation,
                                prettifyList,
                                genotypeOrder = c(wt, background, perturbation,
                                               spaste(background, "/", perturbation)),
                                col = "white", border = 1,
                                ylab = spaste(symbol, " normalized expression"),
                                ...)
{
  data = expr[, index];
  data.lst = lapply(genotypeOrder, function(g) data[genotype==g]);

  labeledBoxplot(x = data.lst, 
    col = col, border = border,
    ylab = ylab, names = prettifyStrings(genotypeOrder, prettifyList),
    ...)
}


#=========================================================================================================
#
# Rescue statistics from a table of results
#
#=========================================================================================================

GNV.rescueStatistics = function(results, background, perturbation, 
    pertBack = spaste(background, "/", perturbation),
    wt = "WT",
    genotypeSep = " vs. ",
    ZPrefix = "Z for ", pPrefix = "p for ", FDRprefix = "FDR for ",
    Z.bg.name = spaste(ZPrefix, background, genotypeSep, wt),
    pert.bg.interName = spaste("interaction of ", perturbation, " with ", background),
    signifThresholdName = " (FDR < 0.1)",
    signifThreshold = list(p = 0.05, pAdj = 0.1),
    PPrefix = "p for ",
    FDRPrefix = "FDR for ",
    rescueColName = "significant rescue",
    exacerbationColName = "significant exacerbation"
)
{
  names(results) = sub("0\\.10", "0.1", names(results));
  Z.pertBack.name = spaste(ZPrefix, pertBack, genotypeSep, background);

  Z.bg = results[[Z.bg.name]];
  Z.pertBack = results[[Z.pertBack.name]];

  cor.rescue = -cor(Z.bg, Z.pertBack, use = 'p');

  rescueColName.x = spaste(rescueColName, signifThresholdName);
  exacerbationColName.x = spaste(exacerbationColName, signifThresholdName);

  rescueInd = results[[rescueColName.x]];
  exacerbationInd = results[[exacerbationColName.x]];

  PCol1 = results[[spaste(PPrefix, perturbation, genotypeSep, wt)]];
  FDRCol1 = results[[spaste(FDRPrefix, perturbation, genotypeSep, wt)]];

  PCol2 = results[[spaste(PPrefix, pertBack, genotypeSep, background)]];
  FDRCol2 = results[[spaste(FDRPrefix, pertBack, genotypeSep, background)]];

  PCol3 = results[[spaste(PPrefix, pert.bg.interName)]];
  FDRCol3 = results[[spaste(FDRPrefix, pert.bg.interName)]];
  browser()

  n1 = sum(replaceMissing( PCol1 < signifThreshold$p & FDRCol1 < signifThreshold$pAdj));
  n2 = sum(replaceMissing( PCol2 < signifThreshold$p & FDRCol2 < signifThreshold$pAdj));
  n3 = sum(replaceMissing( PCol3 < signifThreshold$p & FDRCol3 < signifThreshold$pAdj));
  if (is.null(rescueInd)) browser()
  list(genomewideCor = cor.rescue, 
       rescued.minus.exacerbated = sum(replaceMissing(rescueInd)) - sum(replaceMissing(exacerbationInd)),
       nRescued = sum(replaceMissing(rescueInd)),
       nExacerbated = sum(replaceMissing(exacerbationInd)),
       nDEGenes.for.perturbation.vs.WT = n1,
       nDEGenes.for.perturbation.Q.vs.Q = n2,
       nDEGenes.for.interaction.background.with.Q = n3);
}

GNV.findResults = function(rootPath, genes)
{
  tmpFile = tempfile();
  command = spaste("find ", rootPath, " -type d > ", tmpFile)

  lines = readLines(tmpFile);
  unlink(tmpFile);

  baseNames = sub(".*/", "", lines);
  candidates = lapply(spaste("[^A-Za-z0-9]", genes, "[^A-Za-z0-9]"), grep, baseNames);

  stop("Not finished.")
}
  

#============================================================================================================
#
# For automation of entire pre-processing: create a function that loads and formats the data for
# preprocessing. This function assumes a single data set. If several data sets are contained within one file,
# use separate code to split the data up.
#
#============================================================================================================

GNV.newPreprocessingInfo = function(
  perturbation,
  year, 
  projectID,
  organism = "mouse",
  background = "Q140",
  tissue = "Striatum",
  age = 6, 
  #sampleNamePattern, 
  inputBaseDir = spaste(getBaseDir(), "/Data/HuntingtonsDisease/Mouse/GNVProject/", year, "-", projectID, "-",
                                     perturbation, ".x.", background),
  outputBaseDir.all = spaste(getBaseDir(), "/Data/HuntingtonsDisease/Mouse/GNVProject/CommonReanalysis/"),
  outputBaseDir = spaste(outputBaseDir.all, year, "-", projectID, "-",
                                     perturbation, ".x.", background, ".", age, "m.", tissue),
  exprDir.in = file.path(inputBaseDir, "Expression"),
  exprDir.in.raw = file.path(exprDir.in, "010-RawData"), 

  exprFile = "metaReadCount_ensembl.txt.bz2",
  inFileSep = "\t",
  geneIDType = "Ensembl", ## Can be one of Ensembl, Symbol, Entrez
  geneIDCol = 1, # Can be either index or column name

  sampleDir.in = file.path(inputBaseDir, "SampleAnnotation"),
  sampleDir.in.raw = file.path(sampleDir.in, "010-AsSupplied"),
  sampleFile = NULL,

  manuallyExcludedSamples = NULL,

  annotIndex = c(1:5),
  genotypeTranslation = NULL,

  setName = spaste(perturbation, ".", tissue, ".", age, "m"),
  setName.pretty = spaste(perturbation, " HET KO x ", background, ", ", age, " month ", tolower(tissue)),
  amNameModifier = "",
  amNameSource = spaste("Nan Wang/Yang lab "),
  amAnalysisName = spaste(perturbation, " HET KO in ", background, 
                     " and wild type mice at ", age, " months, ", tolower(tissue), ",", 
                     prependPrefix(" ", appendSuffix(amNameModifier, ",")),
                     prependPrefix(" ", amNameSource), year),
  amAnalysisDescription = spaste("Bulk RNA-seq from ", age, " month old mouse ", tolower(tissue),
             ", ", perturbation, " x ", background, " cross and controls."),
  amDataSource = "Data from Nan Wang/Yang lab and UCLA sequencing core.",

  sampleAnnot.sampleCol = "Sample.Name",
  sampleAnnot.genotypeCol = "Genotype",
  sampleAnnot.genderCol = "Gender",

  phenoCols = c(sampleAnnot.genotypeCol, sampleAnnot.genderCol),

  genderCheckGene = "Xist",
  plotDir = "Plots",

  outlierRemovalZ = 6,
  sva.useNFactors = 0,
  covariates = character(0))
{

  if (is.null(sampleFile))
  {
    candidates = c(
      spaste("UNGC_", year, "-", projectID, "_sample_list_for_", perturbation, "-", background, ".csv"),
      spaste(year, "-", projectID, "_sample_list_for_", perturbation, ".csv"),
      spaste("UNGC_", year, "-", projectID, "_sample_list_", perturbation, "-", background, ".csv"),
      spaste("UNGC_", year, "-", projectID, "_sample_sheet_", perturbation, "-", background, ".csv"),
      spaste("UNGC_", year, "-", projectID, "_", perturbation, "-", background, "_sample_list_.csv"),
      spaste("UNGC_", year, "-", projectID, "_", perturbation, "-", background, "_mice_list_.csv"));
    foundFile = NULL;
    for (file in candidates)
      if (file.exists(file.path(sampleDir.in.raw, file))) foundFile = c(foundFile, file);

    if (length(foundFile)==0)
    {
      candidates = list.files(sampleDir.in.raw, pattern = "csv$");
      if (length(candidates)==1) 
        foundFile = candidates;
    }

    if (length(foundFile)==0) 
        stop("Could not find an existing file among these candidates: \n", 
             paste(candidates, collapse = "\n"));
    if (length(foundFile) > 1)
      warning(immediate = TRUE, spaste("Found multiple candidate files:\n", paste(foundFile, collapse= "\n")));
    sampleFile = foundFile[1];
    printFlush(spaste("Will use sample file ", sampleFile));
  }
  out = list(perturbation = perturbation,
       year = year,
       projectID = projectID,
       organism = organism,
       background = background,
       tissue = tissue,
       age = age,
      # sampleNamePattern = sampleNamePattern,
       inputBaseDir = inputBaseDir,
       outputBaseDir = outputBaseDir,
       exprDir.in = exprDir.in,
       exprDir.in.raw = exprDir.in.raw,
       exprFile = exprFile,
       exprDir.out = file.path(outputBaseDir, "Expression"),
       inFileSep = inFileSep,
       geneIDType = geneIDType,
       geneIDCol = geneIDCol,
       sampleDir.in = sampleDir.in,
       sampleDir.in.raw = sampleDir.in.raw,
       sampleDir.out = file.path(outputBaseDir, "SampleAnnotation"),
       sampleFile = sampleFile,
       manuallyExcludedSamples = manuallyExcludedSamples,
       annotIndex =  annotIndex,
       genotypeTranslation = genotypeTranslation,
       setName = setName,
       setName.pretty = setName.pretty,
       amAnalysisName = amAnalysisName,
       amAnalysisDescription = amAnalysisDescription,
       amDataSource = amDataSource,
       sampleAnnot.sampleCol = sampleAnnot.sampleCol,
       sampleAnnot.genotypeCol = sampleAnnot.genotypeCol,
       sampleAnnot.genderCol = sampleAnnot.genderCol,
       phenoCols = phenoCols,
       genderCheckGene = genderCheckGene,
       plotDir = plotDir,
       outlierRemovalZ = outlierRemovalZ,
       sva.useNFactors = sva.useNFactors,
       covariates = covariates);
  class(out) = c("preprocessingInfo", class(out));
  out;
}
       
       
preprocessingData = function(preprocessingInfo, checkBackground = TRUE, verbose = 1, indent = 0)
{      
 spaces = indentSpaces(indent);
 dataList = with(preprocessingInfo,
 {
  if (verbose > 0) printFlush(spaste(spaces, "Working on ", setName));
  data = read.table(bzfile(file.path(exprDir.in.raw, exprFile)),
            quote = '"', sep = inFileSep, header = TRUE);
  counts0 = t(data[, -annotIndex]);
  colnames(counts0) = data[, geneIDCol];
  geneAnnot0 = data[, annotIndex];

  if (geneIDType=="Ensembl")
  {
    geneAnnot.db = retrieveEnsembleData(type = "Gene", organism = "mouse", verbose = 0);
    geneAnnot = cbind(geneAnnot0, 
                   Start = translateUsingTable(geneAnnot0[, geneIDCol], geneAnnot.db[, c("EnsID", "EnsStart")]),
                   End = translateUsingTable(geneAnnot0[, geneIDCol], geneAnnot.db[, c("EnsID", "EnsEnd")]),
                   Entrez = translateUsingTable(geneAnnot0[, geneIDCol], geneAnnot.db[, c("EnsID", "Entrez")]));
    geneAnnot$Symbol = convert2symbol(entrez = geneAnnot$Entrez, organism = "mouse");
  } else if (geneIDType=="Symbol") {
    entrez = convert2entrez(geneAnnot0[, geneIDCol], organism = organism);
    annot1 = geneAnnotationFromEntrez(entrez, organism = organism, includePosition = TRUE);
    geneAnnot = cbind(geneAnnot0, annot1);
  } else if (geneIDType=="Entrez") {
    geneAnnot = geneAnnot0;
  }
  names(geneAnnot) = multiSub(c("Chromosum", "Chromosom.*"), c("Chr", "Chr"), names(geneAnnot));
  sampleAnnot0 = read.csv(file.path(sampleDir.in.raw, sampleFile));

  if (!is.null(manuallyExcludedSamples))
  {
    printFlush("Excluding marked samples ", paste(manuallyExcludedSamples, collapse = ", "));
    if (any(!manuallyExcludedSamples %in% sampleAnnot0[[sampleAnnot.sampleCol]]))
      stop("Not all 'manuallyExcludedSamples' are present in sample annotation: something is wrong.");
  }

  keepSamples = !sampleAnnot0[[sampleAnnot.sampleCol]] %in% manuallyExcludedSamples
  sampleAnnot02 = sampleAnnot0[keepSamples, ];

  exprSamples = multiGSub(c("^X", "[. _-]"), c("", ""), rownames(counts0));
  saSamples = multiGSub(c("^X", "[. _-]"), c("", ""), sampleAnnot02[[sampleAnnot.sampleCol]]);
  expr2sa = match(exprSamples, saSamples);
  if (verbose > 0) printFlush(spaste(spaces, "..number of unmatched expression samples: ", sum(is.na(expr2sa))));
  unmatchedSamples = exprSamples[is.na(expr2sa)];
  if (any(is.na(expr2sa)))
    printFlush(spaste(spaces, " Unmatched expression samples: ", paste(unmatchedSamples, collapse = ", ")));
  fin = is.finite(expr2sa);

  sampleAnnot1 = data.frame(sampleAnnot02[expr2sa[fin], ], exprSampleName = rownames(counts0)[fin],
       set.short = rep(setName, sum(fin)));
  counts = counts0[sampleAnnot1$exprSampleName, ];

  unmatchedAnnotationSamples = sampleAnnot02[[sampleAnnot.sampleCol]] [ !saSamples %in% exprSamples ];

  sampleAnnot = sampleAnnot1;
  if (!is.null(preprocessingInfo$genotypeTranslation))
    sampleAnnot[[sampleAnnot.genotypeCol]] = 
      translateUsingTable( sampleAnnot1[[sampleAnnot.genotypeCol]], genotypeTranslation);
  sampleAnnot$wellLetter = gsub("[0-9]", "", sampleAnnot1$Well.Address);
  sampleAnnot$wellNumber = as.numeric(sub("[[:alpha:]]*", "", sampleAnnot1$Well.Address));
  
  # Load sex-associated genes

  dir = file.path(getBaseDir(), "HuntingtonsDisease/IndividualAnalyses/070-CHDIAllelicSeries",
                  "CommonAnalysis-030-mRNA-all/060-StandardAnalysisForGender/Results");
  file = spaste("associationWithSex-", capitalize(tissue), ".csv.gz");
  sexAssoc = read.csv(gzfile(file.path(dir, file)));

  Zs = sexAssoc[, grep("stat", colnames(sexAssoc))];
  FDRs = sexAssoc[, grep("FDR", colnames(sexAssoc))];
  sameDir = abs(rowSums(sign(Zs), na.rm = TRUE))==ncol(Zs);
  #table(sameDir, rowSums(FDRs < 0.05, na.rm = TRUE) > 1)
  select = rowSums(FDRs < 0.05, na.rm = TRUE) > 1 & sameDir;
  sexAssocGenes = sexAssoc$Entrez[select];

  # Check gender via a quick nearest neighbor predictor of sex from Xist

  gene = "Xist"
  heCols = heColumns(counts, minValue = 1, minProportion = 0.25, useCountsPerMillion = TRUE)
  ds1 = DESeqDataSetFromMatrix.PL(t(counts[, heCols]),
                                    colData = sampleAnnot,
                                    design = ~1);
  ds1 = estimateSizeFactors(ds1);
  expr1 = log2(t(counts(ds1, normalized = TRUE)) + 1);
  index = match(gene, geneAnnot$Symbol[heCols]);
  cn1 = expr1[, index]

  # Check that gender column was interpreted correctly... all female (F) will be interpreted  as logical with
  # FALSE 
  if (is.logical(sampleAnnot[[sampleAnnot.genderCol]])) 
    sampleAnnot[[sampleAnnot.genderCol]] = rep("F", nrow(sampleAnnot));

  group = apply(sampleAnnot[, c(sampleAnnot.genderCol, sampleAnnot.genotypeCol)], 
                1, base::paste, collapse = ".");
  groupOrder = sort(unique(group));

  suppressWarnings(dir.create(plotDir, recursive = TRUE, showWarnings = FALSE));

  #groupOrder = sort(unique(group));
  pdf(file = spaste(plotDir, "/", setName, "-XistExpressionVsGroup.pdf"), wi = 10, he = 4)
  #sizeGrWindow(10, 4);
  boxplot(cn1~group)

  dev.off();

  #if (length(unique(sampleAnnot[[sampleAnnot.genderCol]])) > 1)
  normals = tapply(cn1, sampleAnnot[[sampleAnnot.genderCol]], median)
  predictedSex = names(normals)[ apply(   do.call(cbind, lapply(normals, function(x) abs(cn1-x))),
                                        1, which.min)]
  if (verbose > 0)
  {
    printFlush(spaste(spaces, "..sex check:"));
    print(table(predictedSex, sampleAnnot[[sampleAnnot.genderCol]]));
  }
  sampleAnnot$predictedGender = predictedSex;
  # Also predict background status from the data. 
  library(randomGLM);
  bg = as.numeric(grepl(background, sampleAnnot[[sampleAnnot.genotypeCol]]));

  if (checkBackground)
  {
    prediction = randomGLM(x = expr1, y = bg, nBags = 50, verbose = 0, nThreads = 2)$predictedOOB;
    if (verbose > 0)
    {
      printFlush(spaste(spaces, "..background genotype check:"));
      print(table(prediction, bg));
    }
    #if (any(prediction!=Q140)) 
    samplesWithSuspiciousQGenotype = sampleAnnot[[sampleAnnot.sampleCol]] [prediction!=bg];
    sampleAnnot$suspiciousQGenotype = as.numeric(prediction!=bg);
  } else {
    samplesWithSuspiciousQGenotype = character(0);
    sampleAnnot$suspiciousQGenotype = rep(0, nrow(sampleAnnot));
  }

  list(counts = counts, sampleAnnot = sampleAnnot, geneAnnot = geneAnnot,
       sexAssocGenes = sexAssocGenes, samplesWithSuspiciousQGenotype = samplesWithSuspiciousQGenotype)

  });
  out = c(preprocessingInfo, dataList);
  class(out) = c("preprocessingData", class(preprocessingInfo));

  out;
}

