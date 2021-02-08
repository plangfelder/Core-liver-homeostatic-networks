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

humanData = loadAsList(file = file.path("../060-NetworkAnalysis-Suppli2019-HumanNAFLD-GSE126848/RData",
                                        "dataForModulePreservation-humanLiver.RData"));

mouseData = loadAsList(file = file.path("../030-NetworkAnalysis/RData/dataForModulePreservation.RData"));

# For plotting: also load the maps between module labels and color labels

humanInfo = loadAsList(file = file.path("../060-NetworkAnalysis-Suppli2019-HumanNAFLD-GSE126848/RData",
                                        "eLabels.RData"));
humanLabelMap = mtd.apply(humanInfo$eLabels2, function(x) apply(x[, c("module", "labelForColor")], 2, as.numeric));

mouseInfo = loadAsList(file = file.path("../030-NetworkAnalysis/RData/eLabels.RData"));
mouseLabelMap = mtd.apply(mouseInfo$eLabels2, function(x) apply(x[, c("module", "labelForColor")], 2, as.numeric));

labelMap = c(mouseLabelMap, humanLabelMap);


# Map human data to mouse Entrez and keep only mapped genes

humanEntrez = lapply(humanData$geneAnnotation.ana, getElement, "Entrez");
humanEntrez.mouse = lapply(humanEntrez, mapEntrez, orgFrom = "human", orgTo = "mouse");
mapped = lapply(humanEntrez.mouse, function(x) !is.na(x));
lapply(mapped, table)
#FALSE  TRUE 
# 1512 13979 

humanExpr = unlist(removeListNames(mymapply(function(x, .entrez) 
{
  mtd.setColnames(mtd.subset(x, , !is.na(.entrez)), .entrez[!is.na(.entrez)]);
}, humanData$expr.ana, humanEntrez.mouse)), recursive = FALSE);

humanWeights = unlist(removeListNames(mymapply(function(x, .entrez)  
{
  mtd.setColnames(mtd.subset(x, , !is.na(.entrez)), .entrez[!is.na(.entrez)]);
}, humanData$weights.ana, humanEntrez.mouse)), recursive = FALSE);


humanLabels = mtd.mapply(function(.labels, .mapped)
  c(.labels)[.mapped],
humanData$labels.x, c(mapped, mapped)) 

# Restrict mouse data to genes with valid Entrez ID

mouseValidEntrez = lapply(mouseData$geneAnnotation.ana, function(.ga) !is.na(.ga$Entrez));

mouseExpr = unlist(removeListNames(mymapply(function(x, ga)
{
  mtd.setColnames(mtd.subset(x, , !is.na(ga$Entrez)), ga$Entrez[!is.na(ga$Entrez)]);
}, mouseData$expr.ana, mouseData$geneAnnotation.ana)), recursive = FALSE);

mouseWeights = unlist(removeListNames(mymapply(function(x, ga)
{
  mtd.setColnames(mtd.subset(x, , !is.na(ga$Entrez)), ga$Entrez[!is.na(ga$Entrez)]);
}, mouseData$weights.ana, mouseData$geneAnnotation.ana)), recursive = FALSE);

mouseLabels = mtd.mapply(function(.labels, .valid)
  c(.labels)[.valid],
mouseData$labels.x, c(mouseValidEntrez, mouseValidEntrez))

# Put together data for preservation

presExpr = c(mouseExpr, mouseExpr, humanExpr, humanExpr);
names(presExpr) = c("mouse.main", "mouse.sub", "human.main", "human.sub");

presWeights = c(mouseWeights, mouseWeights, humanWeights, humanWeights);
names(presWeights) = names(presExpr);

presLabels = c(multiData2list(mouseLabels), multiData2list(humanLabels));

names(presLabels) = names(labelMap) = names(presExpr);

setNames = names(presExpr);
setNames.pretty = spaste(sub("\\.", " ", setNames), " modules");
setOrganism = sub("\\..*", "", setNames);

#===================================================================================================================
#
# Run calculation if data cannot be loaded
#
#===================================================================================================================

refNetworks = c(1:4);
testNetworks = list(3,4,1,2);
if (!checkRDataFile("RData/mp.RData"))
{
   print(system.time( {
   mp = modulePreservation(
      multiData = presExpr,
      multiColor = presLabels,
      multiWeights = presWeights,
      networkType = "signed hybrid",
      referenceNetworks = refNetworks,
      testNetworks = testNetworks,
      nPermutations = 500,
      savePermutedStatistics = TRUE,
      permutedStatisticsFile = "RData/permutedStats.RData",
      parallelCalculation = FALSE,
      verbose = 4);
   }));
   save(mp, file = "RData/mp.RData");
}


#===================================================================================================================
#
# Get results
#
#===================================================================================================================

getTable = function(object, comp1, comp2, ref, test, dropInfo, infoCols = 1, name, dropRowNames = character(0))
{
  tab = object[[comp1]] [[comp2]] [[ref]] [[test]];
  tab2 = tab[, -infoCols];
  keepRows = !(rownames(tab) %in% dropRowNames);
  colnames(tab2) = spaste(colnames(tab2), ".in.", name);
  if (!dropInfo) tab2 = data.frame(module = rownames(tab), tab[, infoCols, drop = FALSE], 
                                 tab2);
  tab2 = tab2[keepRows, ]
  rn = rownames(tab2)
  num = suppressWarnings(!any(is.na(as.numeric(rn))));
  if (num) tab2 = tab2[order(as.numeric(rn)), ];
  tab2;
}

# Write out preservation statistic tables

for (ref in refNetworks)
{
  pres.obs = pres.Z = acc.obs = acc.Z = list();
  allStats = list();
  for (ts in testNetworks[[ref]])
  {
    refName = spaste("ref.", setNames[ref]);
    testName = spaste("inColumnsAlsoPresentIn.", setNames[ts]);
    
    pres.obs1 = getTable(mp, "preservation", "observed", refName, testName, 
                         dropInfo = FALSE, name = setNames[ts], dropRowNames = c("0", "0.1"));
    pres.Z1 = getTable(mp, "preservation", "Z", refName, testName, 
                         dropInfo = FALSE, name = setNames[ts], dropRowNames = c("0", "0.1"));
    acc.obs1 = getTable(mp, "accuracy", "observed", refName, testName, 
                         dropInfo = FALSE, name = setNames[ts], dropRowNames = c("0", "0.1"));
    acc.Z1 = getTable(mp, "accuracy", "Z", refName, testName, 
                         dropInfo = FALSE, name = setNames[ts], dropRowNames = c("0", "0.1"));
    
    pres.obs = c(pres.obs, pres.obs1)
    pres.Z = c(pres.Z, pres.Z1)
    acc.obs = c(acc.obs, acc.obs1);
    acc.Z = c(acc.Z, acc.Z1);

    allStats = c(allStats, list(dropDuplicatedColumns(cbind(pres.obs1, pres.Z1, acc.obs1, acc.Z1), matchColnames = TRUE)));
  }

  allStats = dropDuplicatedColumns(interleave(allStats, nameBase = rep("", length(allStats)), sep = "", check.names = FALSE),
                   matchColnames = TRUE)
  allStats = allStats[, grep("moduleSize\\.[1-9]", names(allStats), invert = TRUE)];

  colOrder = order(1-multiGrepl(c("module$", "moduleSize", "medianRank.pres", "Zsummary.pres"), names(allStats)));

  allStats2 = allStats[, colOrder]

  write.csv(signifNumeric(allStats2, 3), 
            file = spaste("Results/preservationStatistics-", setNames[ref], ".in.", 
                          paste(setNames[testNetworks[[ref]] ], collapse = "."), ".csv"),
            row.names = FALSE);
}

# Plot 1: standard plots of Zsummary and preservationRank.summary vs. module size 
plotSize = 6;
limStretch = 0.06;
thresholds = c(4, 8);
thresholdColors = c("blue", "red");

for (ref in refNetworks)
{
  pdf(file = spaste("Plots/modulePreservation-scatterplots-", setNames[ref], ".pdf"), wi = plotSize, he = plotSize);
  scpp(3);
  for (ts in testNetworks[[ref]])
  {
    obs = mp$preservation$observed[[spaste("ref.", setNames[ref])]] [[spaste("inColumnsAlsoPresentIn.", setNames[ts])]];

    modules1 = as.numeric(rownames(obs));
    keep = modules1 >=1
    obs = obs[keep, ];
    sizes = obs$moduleSize;
    modules1 = modules1[keep];
    #sizeGrWindow(6,6)
    plot(sizes, obs$medianRank.pres, 
         xlim = stretchLim(sizes, limStretch), ylim = stretchLim(obs$medianRank.pres, limStretch),
         ylab = "Median preservation rank", xlab = "Module size",
         pch = 21, bg = labels2colors(translateUsingTable(modules1, labelMap[[ref]]$data)),
         cex = 2, 
         main = spaste("Median preseration rank\n", setNames.pretty[ref], " in ", setOrganism[ts]));
    labelPoints2(x = sizes, y = obs$medianRank.pres, labels = spaste("M", modules1),
                 ratio.pointToChar = 0.6);

    plot(log10(sizes), obs$medianRank.pres, 
         xlim = stretchLim(log10(sizes), limStretch), ylim = stretchLim(obs$medianRank.pres, limStretch),
         ylab = "Median preservation rank", xlab = "log10(module size)",
         pch = 21, bg = labels2colors(translateUsingTable(modules1, labelMap[[ref]]$data)),
         cex = 2,
         main = spaste("Median preseration rank\n", setNames.pretty[ref], " in ", setOrganism[ts]));
    labelPoints2(x = log10(sizes), y = obs$medianRank.pres, labels = spaste("M", modules1),
                 ratio.pointToChar = 0.6);
    Z = mp$preservation$Z[[spaste("ref.", setNames[ref])]] [[spaste("inColumnsAlsoPresentIn.", setNames[ts])]];
    Z = Z[keep, ];
    plot(sizes, Z$Zsummary.pres, 
         xlim = stretchLim(sizes, limStretch), ylim = stretchLim(Z$Zsummary.pres, limStretch),
         ylab = "Preservation Zsummary", xlab = "Module size",
         pch = 21, bg = labels2colors(translateUsingTable(modules1, labelMap[[ref]]$data)),
         cex = 2,
         main = spaste("Preservation Zsummary\n", setNames.pretty[ref], " in ", setOrganism[ts]));
    for (th in 1:length(thresholds)) abline(h = thresholds[th], col = thresholdColors[th], lty = 2);
    labelPoints2(x = sizes, y = Z$Zsummary.pres, labels = spaste("M", modules1),
                 ratio.pointToChar = 0.6);
    plot(log10(sizes), Z$Zsummary.pres, 
         xlim = stretchLim(log10(sizes), limStretch), ylim = stretchLim(Z$Zsummary.pres, limStretch),
         ylab = "Preservation Zsummary", xlab = "log10(module size)",
         pch = 21, bg = labels2colors(translateUsingTable(modules1, labelMap[[ref]]$data)),
         cex = 2,
         main = spaste("Preservation Zsummary\n", setNames.pretty[ref], " in ", setOrganism[ts]));
    for (th in 1:length(thresholds)) abline(h = thresholds[th], col = thresholdColors[th], lty = 2);
    labelPoints2(x = log10(sizes), y = Z$Zsummary.pres, labels = spaste("M", modules1),
                 ratio.pointToChar = 0.6);
  }
  dev.off();
}


#===============================================================================================================
#
# Plot overlap heatmaps of module labels
#
#===============================================================================================================

plotPairs = list(main = c("mouse.main", "human.main"), sub = c("mouse.sub", "human.sub"))

for (pp in 1:length(plotPairs))
{
  ppp = plotPairs[[pp]];
  mexpr1 = presExpr[ppp];
  common = multiIntersect(mtd.apply(mexpr1, colnames, returnList = TRUE));
  common2ind = mtd.apply(mexpr1, function(x) match(common, colnames(x)));
  labels1 = mtd.mapply(function(l, i) l[i], presLabels[ppp], common2ind)

  calculateAndPlotModuleOverlaps(
    labels1[[1]]$data,
    labels1[[2]]$data,
    colLabelPrefix = spaste(multiSubr(c(".main", ".sub"), ppp[2]), " M"),
    rowLabelPrefix = spaste(multiSubr(c(".main", ".sub"), ppp[1]), " M"),
    colLabelTranslationForColor = labelMap[[ ppp[2] ]]$data,
    rowLabelTranslationForColor = labelMap[[ ppp[1] ]]$data,
    plotFile = spaste("Plots/moduleOverlaps-", names(plotPairs)[pp], ".pdf"),
    threshold = 1e-2,
    mar = c(8, 10, 2, 1),
    main = spaste("Overlaps of mouse and human ", names(plotPairs)[pp], " modules"),
    separatorInterval = 3);
}
