source("../Functions/networkFunctions-extras-19.R")
source("../Functions/labelPoints2-01.R");
source("../Functions/outlierRemovalFunctions.R")
source("../Functions/preprocessing-General-013.R")
source("../Functions/GNVFunctions-016-02.R")
source("../Functions/individualAnalysis-General-008-02.R");

library(anRichment)

dir.create("RData", recursive = TRUE);
dir.create("Results", recursive = TRUE);
dir.create("Plots", recursive = TRUE);


# Load data

humanData = loadAsList(file = file.path("../060-NetworkAnalysis-Suppli2019-HumanNAFLD-GSE126848/RData",
                                        "dataForModulePreservation-humanLiver.RData"));

mouseData = loadAsList(file = file.path("../030-NetworkAnalysis/RData/dataForModulePreservation.RData"));

TCGAData = loadAsList(file = file.path("../130-NetworkAnalysis-TCGAPhase/RData/dataForModulePreservation-TCGA.RData"));

# For plotting: also load the maps between module labels and color labels

humanInfo = loadAsList(file = file.path("../060-NetworkAnalysis-Suppli2019-HumanNAFLD-GSE126848/RData",
                                        "eLabels.RData"));
humanLabelMap = mtd.apply(humanInfo$eLabels2, function(x) apply(x[, c("module", "labelForColor")], 2, as.numeric));

mouseInfo = loadAsList(file = file.path("../030-NetworkAnalysis/RData/eLabels.RData"));
mouseLabelMap = mtd.apply(mouseInfo$eLabels2, function(x) apply(x[, c("module", "labelForColor")], 2, as.numeric));

TCGAInfo = loadAsList(file = file.path("../130-NetworkAnalysis-TCGAPhase/RData/eLabels.RData"));
TCGALabelMap = mtd.apply(TCGAInfo$eLabels2, function(x) apply(x[, c("module", "labelForColor")], 2, as.numeric));

labelMap0 = c(mouseLabelMap, humanLabelMap, TCGALabelMap);
names(labelMap0) = sub("humanLiver", "HumanNAFLDLiver", names(labelMap0));


# Map human data to mouse Entrez and keep only mapped genes

humanEntrez = c(
                 lapply(humanData$geneAnnotation.ana, getElement, "Entrez"),
                 lapply(TCGAData$geneAnnotation.ana, getElement, "Entrez"))

humanEntrez.mouse = lapply(humanEntrez, mapEntrez, orgFrom = "human", orgTo = "mouse");
mapped = lapply(humanEntrez.mouse, function(x) !is.na(x));
lapply(mapped, table)
#FALSE  TRUE 
# 1512 13979 
#[[2]]
#
#FALSE  TRUE 
#13065 14696 
#
#[[3]]
#
#FALSE  TRUE 
# 2818 14749 
#
#[[4]]
#
#FALSE  TRUE 
# 4585 15285 


humanExpr = unlist(removeListNames(mymapply(function(x, .entrez) 
{
  mtd.setColnames(mtd.subset(x, , !is.na(.entrez)), .entrez[!is.na(.entrez)]);
}, c(humanData$expr.ana, TCGAData$expr.ana), humanEntrez.mouse)), recursive = FALSE);

humanWeights = unlist(removeListNames(mymapply(function(x, .entrez)  
{
  mtd.setColnames(mtd.subset(x, , !is.na(.entrez)), .entrez[!is.na(.entrez)]);
}, c(humanData$weights.ana, TCGAData$weights.ana), humanEntrez.mouse)), recursive = FALSE);


humanLabels0 = mtd.mapply(function(.labels, .mapped)
  c(.labels)[.mapped],
c(humanData$labels.x, TCGAData$labels.x), mapped[c(1,1,2,3,4,2,3,4)]); 

humanLabels = humanLabels0[c(1,3,4,5,2,6,7,8)];

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
names(presExpr) = c(spaste(names(mouseExpr), ".main"), spaste(names(mouseExpr), ".sub"),
                    spaste(names(humanExpr), ".main"), spaste(names(humanExpr), ".sub"));
 
presWeights = c(mouseWeights, mouseWeights, humanWeights, humanWeights);
names(presWeights) = names(presExpr);

presLabels = c(multiData2list(mouseLabels), multiData2list(humanLabels));

names(presLabels) = names(labelMap) = names(presExpr);

setNames = names(presExpr);
setNames.pretty = spaste(multiSub(c("\\.", "NAFLDLiver", "^Liver"), c(" ", " NAFLD liver", "Mouse liver"), setNames), 
                         " modules");
setOrganism = c("human", "mouse")[ 1 + grepl("^Liver", setNames)];

labelMap = labelMap0[setNames]

#===================================================================================================================
#
# Run calculation if data cannot be loaded
#
#===================================================================================================================

nSets = length(presExpr);
refNetworks = c(1:nSets);
testNetworks = c(list( c(4,5,6), c(8,9,10), c(4,5,6)),
                    listRep(c(1,3), 3), 
                    list(c(8,9,10)),
                    listRep(c(2,7), 3));

testNetworks.hr = lapply(testNetworks, function(i) setNames[i]);
names(testNetworks.hr) = setNames;

#mymapply(function(x, y) isTRUE(all.equal(x, y)), testNetworks.hr, testNetworks.hr.old)

if (!checkRDataFile("RData/mp.RData"))
{
  mp = list();
  for (ref in refNetworks)
  {
    print(system.time( {
    mp[[ref]] = modulePreservation(
      multiData = presExpr,
      multiColor = presLabels,
      multiWeights = presWeights,
      networkType = "signed hybrid",
      referenceNetworks = refNetworks[ref],
      testNetworks = testNetworks[ref],
      nPermutations = 200,
      randomSeed = ref*2 + 1,
      savePermutedStatistics = TRUE,
      permutedStatisticsFile = spaste("RData/permutedStats-", ref, "-", setNames[ref], ".RData"),
      parallelCalculation = FALSE,
      verbose = 4);
    }));
    save(mp, file = "RData/mpBySet.RData");
  }
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
    
    pres.obs1 = getTable(mp[[ref]], "preservation", "observed", refName, testName, 
                         dropInfo = FALSE, name = setNames[ts], dropRowNames = c("0", "0.1"));
    pres.Z1 = getTable(mp[[ref]], "preservation", "Z", refName, testName, 
                         dropInfo = FALSE, name = setNames[ts], dropRowNames = c("0", "0.1"));
    acc.obs1 = getTable(mp[[ref]], "accuracy", "observed", refName, testName, 
                         dropInfo = FALSE, name = setNames[ts], dropRowNames = c("0", "0.1"));
    acc.Z1 = getTable(mp[[ref]], "accuracy", "Z", refName, testName, 
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
    obs = mp[[ref]]$preservation$observed[[spaste("ref.", setNames[ref])]] [[spaste("inColumnsAlsoPresentIn.", setNames[ts])]];

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
         cex = 2,  cex.main = 1,
         main = spaste("Median preseration rank\n", setNames.pretty[ref], " in ", setNames.pretty[ts]));
    labelPoints2(x = sizes, y = obs$medianRank.pres, labels = spaste("M", modules1),
                 ratio.pointToChar = 0.6);

    plot(log10(sizes), obs$medianRank.pres, 
         xlim = stretchLim(log10(sizes), limStretch), ylim = stretchLim(obs$medianRank.pres, limStretch),
         ylab = "Median preservation rank", xlab = "log10(module size)",
         pch = 21, bg = labels2colors(translateUsingTable(modules1, labelMap[[ref]]$data)),
         cex = 2, cex.main = 1,
         main = spaste("Median preseration rank\n", setNames.pretty[ref], " in ", setNames.pretty[ts]));
    labelPoints2(x = log10(sizes), y = obs$medianRank.pres, labels = spaste("M", modules1),
                 ratio.pointToChar = 0.6);
    Z = mp[[ref]]$preservation$Z[[spaste("ref.", setNames[ref])]] [[spaste("inColumnsAlsoPresentIn.", setNames[ts])]];
    Z = Z[keep, ];
    plot(sizes, Z$Zsummary.pres, 
         xlim = stretchLim(sizes, limStretch), ylim = stretchLim(Z$Zsummary.pres, limStretch),
         ylab = "Preservation Zsummary", xlab = "Module size",
         pch = 21, bg = labels2colors(translateUsingTable(modules1, labelMap[[ref]]$data)),
         cex = 2, cex.main = 1,
         main = spaste("Preservation Zsummary\n", setNames.pretty[ref], " in ", setNames.pretty[ts]));
    for (th in 1:length(thresholds)) abline(h = thresholds[th], col = thresholdColors[th], lty = 2);
    labelPoints2(x = sizes, y = Z$Zsummary.pres, labels = spaste("M", modules1),
                 ratio.pointToChar = 0.6);
    plot(log10(sizes), Z$Zsummary.pres, 
         xlim = stretchLim(log10(sizes), limStretch), ylim = stretchLim(Z$Zsummary.pres, limStretch),
         ylab = "Preservation Zsummary", xlab = "log10(module size)",
         pch = 21, bg = labels2colors(translateUsingTable(modules1, labelMap[[ref]]$data)),
         cex = 2, cex.main = 1,
         main = spaste("Preservation Zsummary\n", setNames.pretty[ref], " in ", setNames.pretty[ts]));
    for (th in 1:length(thresholds)) abline(h = thresholds[th], col = thresholdColors[th], lty = 2);
    labelPoints2(x = log10(sizes), y = Z$Zsummary.pres, labels = spaste("M", modules1),
                 ratio.pointToChar = 0.6);
  }
  dev.off();
}

#===============================================================================================================
#
# Plot heatmaps of preservation Z summary
#
#===============================================================================================================

useRef = 1:nSets

for (ref in useRef)
{
  ts1 = testNetworks[[ref]] [1];
  obs = mp[[ref]]$preservation$observed[[spaste("ref.", setNames[ref])]] [[spaste("inColumnsAlsoPresentIn.", setNames[ts1])]];
  modules1 = as.numeric(rownames(obs));
  keep = modules1 >=1
  obs = obs[keep, ];
  sizes = obs$moduleSize;
  modules1 = modules1[keep];
  morder = order(as.numeric(modules1));
  nm = length(modules1);
  Zs = do.call(cbind, lapply(testNetworks[[ref]], function(ts)
  {
    Z = mp[[ref]]$preservation$Z[[spaste("ref.", setNames[ref])]] [[spaste("inColumnsAlsoPresentIn.", setNames[ts])]];
    Z$Zsummary.pres[keep];
  }));
  colnames(Zs) = setNames.pretty[ testNetworks[[ref]] ];
  Zs = Zs[morder, ]
  modules1 = modules1[morder]
  sizes2 = sizes[morder];
  lim = min(20, max(Zs, na.rm = TRUE));
  text = round(Zs, 1);
  text[Zs < 2] = "";
  xlabs = multiSubr(c(" (main|sub) modules"), setNames.pretty[testNetworks[[ref]]]);
  mw = marginWidth(xlabs);
  pdf(file = spaste("Plots/modulePreservation-heatmaps-", setNames[ref], ".pdf"), width = 4.5, 
      height = 0.5 + 0.8*mw$inch + 0.2 * nm);
  par(mar = c(0.8*mw$lines, 7, 2, 1));
  labeledHeatmap(Zs,
    xLabels = xlabs,
    yLabels = spaste("ME", labels2colors(translateUsingTable(modules1, labelMap[[ref]]$data))),
    ySymbols = spaste("M", modules1, " (", sizes2, ")"),
    colors = blueWhiteRed(100)[50:100],
    zlim = c(0, lim),
    textMatrix = text,
    setStdMargins = FALSE,
    main = spaste("Zsummary for ", setNames.pretty[ref], "                "),
    cex.main = 1);

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
