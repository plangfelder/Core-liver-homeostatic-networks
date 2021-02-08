source("../../Functions/networkFunctions-extras-20.R")
source("../../Functions/labelPoints2-01.R");
source("../../Functions/heatmap.wg.R");

source("../../Functions/outlierRemovalFunctions.R")
source("../../Functions/preprocessing-General-013.R")
source("../../Functions/GNVFunctions-015.R")
source("../../Functions/individualAnalysis-General-007-02.R");

library(anRichment)

dir.create("RData", recursive = TRUE);
dir.create("Results", recursive = TRUE);
dir.create("Plots", recursive = TRUE);


# Load data

x = load("../../010-Preprocessing/RData/preprocessedData.RData");

prettifyList = loadAsList("../../020-IndividualAnalysis/RData/prettifyList2.RData")[[1]]

multiCounts = multiExpr;

tissues = setNames
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


library("MODA");

# Run the three module identification schemes in MODA.

power = 8

dir.create("Results/MODA-hc", recursive = TRUE);
system2("rm", args = c("Results/MODA-hc/*"));
system.time( { moda.hc = WeightedModulePartitionHierarchical(expr.flat[[1]]$data, foldername = "Results/MODA-hc",
                                indicatename = "X", cutmethod = "Density", power = power) });

#   user  system elapsed 
#122.461   8.628 130.770 


dir.create("Results/MODA-louvain", recursive = TRUE)
system2("rm", args = c("Results/MODA-louvain/*"));
system.time( {moda.louvain = WeightedModulePartitionLouvain(expr.flat[[1]]$data, foldername = "Results/MODA-louvain",
                               indicatename = "X", GeneNames = multiGeneAnnot[[1]]$data$Symbol, maxsize = 2000, power = power) })

#  user  system elapsed 
# 32.053   4.135  34.981 

WeightedModulePartitionSpectral = function (datExpr, foldername, indicatename, GeneNames, power = 6, 
    nn = 10, k = 2) 
{
    dir.create(file.path("./", foldername), showWarnings = FALSE)
    ADJ1 = abs(cor(datExpr, use = "p"))^power
    W = TOMsimilarity(ADJ1)
    A = MODA:::make.affinity(W, nn)
    d <- apply(A, 1, sum)
    L <- diag(d) - A
    L <- diag(d^-0.5) %*% L %*% diag(d^-0.5)
    evL <- eigen(L, symmetric = TRUE)
    Z <- evL$vectors[, (ncol(evL$vectors) - k + 1):ncol(evL$vectors)]
    spc <- pam(Z, k)
    colorSpectralTom <- labels2colors(spc$cluster)
    intModules = table(colorSpectralTom)
    for (J in 1:length(intModules)) {
        idx <- which(colorSpectralTom == names(intModules)[J])
        DenseGenes = GeneNames[idx]
        densegenefile <- paste(foldername, "/DenseModuleGene_", 
            indicatename, "_", J, ".txt", sep = "")
        write.table(idx, file = paste(foldername, "/DenseModuleGeneID_", 
            indicatename, "_", J, ".txt", sep = ""), quote = FALSE, 
            row.names = FALSE, col.names = FALSE)
        write.table(DenseGenes, densegenefile, sep = "\n", col.names = FALSE, 
            row.names = FALSE, quote = FALSE)
    }
}


dir.create("Results/MODA-spectral", recursive = TRUE)
system2("rm", args = c("Results/MODA-spectral/*"));
system.time( {moda.spectral = WeightedModulePartitionSpectral(expr.flat[[1]]$data, foldername = "Results/MODA-spectral",
                    indicatename = "X", GeneNames = multiGeneAnnot[[1]]$data$Symbol, nn = 30, k = 20, power = power) })

#    user   system  elapsed 
#8656.928  151.905 5537.983 


library(minet)
system.time( {mim = build.mim(expr.flat[[1]]$data);});

#   user  system elapsed 
# 32.541   4.123  37.199 


system.time( {aracneNet = aracne(mim = mim);});

#     user   system  elapsed 
#17774.57    17.38 17939.87 

save(aracneNet, file = "RData/aracneNet.RData");

aracneNet2 = aracneNet;
aracneNet2[aracneNet>2] = 2;
tree = hclust(1-as.dist(aracneNet2), method = "a")

mods0 = cutreeDynamic(tree, distM = 1-aracneNet2, deepSplit = 4, minClusterSize = 20)

merge = mergeCloseModules(exprData = expr.flat[[1]]$data, colors = mods0, cutHeight = 0.2);
aracneMods = merge$colors

plotDendroAndColors(tree, labels2colors(aracneMods), dendroLabels = FALSE, hang = 0.01)


system.time( {mrNet = mrnet(mim = mim); })
#      user     system    elapsed 
# 77511.593    192.465 209409.804 
### elapsed is loooong because of interruptions
# but even the user time is over 20 hours...


save(mrNet, file = "RData/mrNet.RData");

#load(file = "RData/mrNet.RData");
gc()

tree = hclust(1-as.dist(mrNet), method = "a")

mods0 = cutreeDynamic(tree, distM = 1-mrNet, deepSplit = 4, minClusterSize = 20)

merge = mergeCloseModules(exprData = expr.flat[[1]]$data, colors = mods0, cutHeight = 0.2, relabel = TRUE);
mrMods = merge$colors

plotDendroAndColors(tree, labels2colors(mrMods), dendroLabels = FALSE, hang = 0.01)



# Put together a list of module labels

# Retrieve MODA modules

retrieveMODAmodules = function(expr, dir)
{
  labels = rep(0, ncol(expr));
  files = list.files(dir, pattern = "DenseModuleGeneID_X");
  for (f in files)
  {
    lab = as.numeric(multiSubr(c("DenseModuleGeneID_X_", ".txt"), f));
    index = read.table(file.path(dir, f), header = FALSE)[[1]];
    labels[index] = lab;
  }
  labels;
}

testLabels = list(MODA.hc = retrieveMODAmodules(expr.flat[[1]]$data, dir = "Results/MODA-hc"),
                  MODA.louvain = retrieveMODAmodules(expr.flat[[1]]$data, dir = "Results/MODA-louvain"),
                  MODA.spectral = retrieveMODAmodules(expr.flat[[1]]$data, dir = "Results/MODA-spectral"),
                  aracne = aracneMods,
                  mrnet = mrMods);

#====================================================================================================
#
# Evaluate concordance with our own WGCNA modules.
#
#====================================================================================================

refNet = read.csv("../../030-NetworkAnalysis/Results/networkAnalysisResults-Liver.sub.csv.gz", check.names = FALSE) 

eLabels = loadAsList("../../030-NetworkAnalysis/RData/eLabels.RData")[[1]][["Liver.sub"]]$data

modData = loadAsList("../../030-NetworkAnalysis/RData/labels.x.etc.RData");

labelTranslationForColor = data.frame(label = sort(unique(refNet$module)),
          colorLabel = modData$subLabelsForColors[[1]]$data[ match(sort(unique(refNet$module)), refNet$module)]);

adjRand = lapply(testLabels, mclust::adjustedRandIndex, refNet$module)

# The Rand indices are low.

# Test module overlaps.

t = table(refNet$module);
modSizes.ref = t[names(t)!="0"];
modSizes.test = lapply(testLabels, function(l) table(l[l!=0]));

ots = lapply(testLabels, overlapTable, labels1 = refNet$module, ignore = 0);
testNetNames = names(testLabels);
testNetNames.pretty = c("MODA-hc", "MODA-Louvain", "MODA-spectral", "ARACNE", "MRNET");
for (ts in 1:length(testLabels))
  plotOverlapHeatmap(ots[[ts]]$countTable, ots[[ts]]$pTable,
      colLabels = spaste("ME", labels2colors(as.numeric(colnames(ots[[ts]]$countTable)))),
      colSymbols = spaste(testNetNames.pretty[ts], " M", colnames(ots[[ts]]$countTable), " (",
                          modSizes.test[[ts]] [ match(colnames(ots[[ts]]$countTable), names(modSizes.test[[ts]]))], ")"),
      rowLabels = spaste("ME", labels2colors(translateUsingTable(rownames(ots[[ts]]$countTable), labelTranslationForColor))),
      rowSymbols = spaste("WGCNA M", rownames(ots[[ts]]$countTable), " (",
                    modSizes.ref [ match(rownames(ots[[ts]]$countTable), names(modSizes.ref))], ")"),
      main = spaste("Overlap of mouse liver WGCNA modules with ", testNetNames.pretty[ts], " modules             "),
      cex.main = 1,
      baseWidth = 3, baseHeight = 2.7,
      plotFile = spaste("Plots/overlapHeatmap-mouseWGCNA.vs.mouse", testNetNames[ts], ".pdf"),
      mar = c(10, 10, 2, 1),
      logpLimit = 100,
      threshold = 1e-20, separatorInterval = 3);


# Get module eigengenes and calculate correlations.

MEs.ref = moduleEigengenes(expr.flat[[1]]$data, refNet$module, excludeGrey = TRUE)$eigengenes;

MEs.test = lapply(testLabels, function(lab) moduleEigengenes(expr.flat[[1]]$data, lab, excludeGrey = TRUE)$eigengenes)



threshold = 0.7;
for (ts in 1:length(testLabels))
{
  cor1 = cor(MEs.ref, MEs.test[[ts]]);
  keepCols = colSums(abs(cor1)>threshold) > 0;
  keepRows = rowSums(abs(cor1)>threshold) > 0; 
  if (sum(keepRows)==0) next;
  cor1 = cor1[keepRows, keepCols];
  txt = round(cor1, 2);
  txt[abs(cor1) < 0.5] = "";
  colMod = as.numeric(sub("ME", "", colnames(cor1)));
  rowMod = as.numeric(sub("ME", "", rownames(cor1)));

  width = 3 + ncol(cor1) * 0.4;
  height = 2.8 + 0.2 * nrow(cor1);
  pdf(file = spaste("Plots/MECorrelations-mouseWGCNA.vs.mouse", testNetNames[ts], ".pdf"), width = width,
        height = height);
  par(mar = c(10,10,2,1));
  wgHeatmap(cor1, 
     colLabels.bottom = spaste(testNetNames.pretty[ts], " M", colMod, " (",
                    modSizes.test[[ts]] [ match(colMod, names(modSizes.test[[ts]]))], ")"),
     colColors.bottom = labels2colors(colMod),
     rowLabels.left = spaste("WGCNA M", rowMod, " (",
                    modSizes.ref [ match(rowMod, names(modSizes.ref))], ")"),
     rowColors.left = labels2colors(translateUsingTable(rowMod, labelTranslationForColor)),
     textMatrix = txt,
     colors = blueWhiteRed(100),
     zlim = c(-1, 1),
     cex.text = 0.9,
     main =  spaste("Correlations of eigengenes of WGCNA modules with ", testNetNames.pretty[ts], " modules"),
     verticalSeparator.interval = 3,
     horizontalSeparator.interval = 3);

  dev.off();
}


#=================================================================================================================
#
# Run enrichment analyses on the test modules.
#
#=================================================================================================================

x = load("../../020-IndividualAnalysis/RData/stdCollections.RData");

stdCollections.combined = do.call(mergeCollections, stdCollections);

identifiers = refNet$Entrez;

enr = mtd.mapply(function(lab1, ids)
     enrichmentAnalysis(classLabels= lab1, identifiers = ids,
                     refCollection = stdCollections.combined, ignoreLabels = 0, threshold = 5e-2,
                     nBestDataSets = 5,
                     getOverlapSymbols = TRUE, getOverlapEntrez = FALSE,
                     getDataSetDetails = FALSE, useBackground = "given",
                     maxReportedOverlapGenes = 10),
    testLabels, MoreArgs = list(ids = identifiers), mdmaVerbose = TRUE)

save(enr, file = "RData/enrichment-enr.RData");

#=================================================================================================================
#
# For each reference module, select one corresponding module from the test data
#
#=================================================================================================================


focusOn = c(1.1, 2, 3.1, 7);
nFocus = length(focusOn);
nTests = length(testLabels);
bestOverlapM.p = matrix(NA, nFocus, nTests);
bestOverlapM.n = matrix(NA, nFocus, nTests);
bestOverlapN = matrix(NA, nFocus, nTests);
bestOverlapP = matrix(NA, nFocus, nTests);
bestCorM = matrix(NA, nFocus, nTests);
bestCor = matrix(NA, nFocus, nTests);
bestModuleSize = matrix(NA, nFocus, nTests);

# Try seleting the module with strongest overlap p-value.
for (f in 1:nFocus) for (ts in 1:nTests)
{
  ff = focusOn[f]
  p1 = ots[[ts]]$pTable[ match(ff, rownames(ots[[ts]]$pTable)), ];
  bestOverlapP[f, ts] = min(p1);
  m1 = colnames(ots[[ts]]$pTable)[which.min(p1)];
  bestOverlapM.p[f, ts] = m1;

  bestModuleSize[f, ts] = modSizes.test[[ts]] [ match(m1, names(modSizes.test[[ts]] ))];
  
  n1 = ots[[ts]]$countTable[ match(ff, rownames(ots[[ts]]$pTable)), ];
  bestOverlapN[f, ts] = ots[[ts]]$countTable[ match(ff, rownames(ots[[ts]]$pTable)), match(m1, colnames(ots[[ts]]$pTable))];
  #bestOverlapN[f, ts] = max(n1);
  bestOverlapM.n[f, ts] = colnames(ots[[ts]]$countTable)[which.max(n1)];

  cor1 = cor(MEs.ref, MEs.test[[ts]]);
  dimnames(cor1) = lapply(dimnames(cor1), function(s) sub("ME", "", s));
  c1 = cor1[ match(ff, rownames(cor1)), ]
  bestCor[f, ts] = c1[ match(m1, colnames(ots[[ts]]$pTable))];
  #bestCor[f, ts] = max(c1);
  bestCorM[f, ts] = colnames(cor1)[which.max(c1)];
};

colnames(bestOverlapM.p) = colnames(bestOverlapM.n) = colnames(bestOverlapN) = 
    colnames(bestOverlapP) = colnames(bestCorM) = colnames(bestCor) = colnames(bestModuleSize) = testNetNames;
  

bestCorrespondenceInfo = data.frame.ncn(`WGCNA module` = focusOn, 
    `WGCNA module size` = as.numeric(modSizes.ref[match(focusOn, names(modSizes.ref))]),
    interleave(list(bestOverlapP, bestOverlapM.p, bestModuleSize, bestOverlapN, bestCor),
               nameBase = c("Top overlap P in ", "Top module in ", "Size of top module in ", 
                            "Overlap size, ", "Eigengene correlation, "),
               sep = "", baseFirst = TRUE, check.names = FALSE));

dir.create("Results", recursive = TRUE);

write.csv.nr(signifNumeric(bestCorrespondenceInfo, 3), file = "Results/bestModuleInfo.csv");

enr.ref = loadAsList("../../030-NetworkAnalysis/RData/enrichment-enr.RData")[[1]]
pMat.ref = enr.ref[[2]]$data$pValues

bestModMat = bestOverlapM.p;

for (ts in 1:length(testLabels))
{
  rows = match(focusOn, rownames(ots[[ts]]$countTable));
  plotOverlapHeatmap(ots[[ts]]$countTable[rows, ], ots[[ts]]$pTable[rows, ],
      colLabels = spaste("ME", labels2colors(as.numeric(colnames(ots[[ts]]$countTable)))),
      colSymbols = spaste(testNetNames.pretty[ts], " M", colnames(ots[[ts]]$countTable), " (",
                          modSizes.test[[ts]] [ match(colnames(ots[[ts]]$countTable), names(modSizes.test[[ts]]))], ")"),
      rowLabels = spaste("ME", labels2colors(translateUsingTable(rownames(ots[[ts]]$countTable)[rows], 
                                                                 labelTranslationForColor))),
      rowSymbols = spaste("WGCNA M", rownames(ots[[ts]]$countTable), " (",
                    modSizes.ref [ match(rownames(ots[[ts]]$countTable), names(modSizes.ref))], ")")[rows],
      main = spaste("Overlap of selected WGCNA modules with ", testNetNames.pretty[ts], " modules             "),
      cex.main = 1,
      baseWidth = 3, baseHeight = 2.7,
      plotFile = spaste("Plots/overlapHeatmap-mouseWGCNA.vs.mouse", testNetNames[ts], "-selectedWGCNA.pdf"),
      mar = c(10, 10, 2, 1),
      logpLimit = 100,
      threshold = 1e-20, separatorInterval = 3);
}



for (ts in 1:nTests)
{
  pdf(file = spaste("Plots/enrichmentScatterplots-", testNetNames[ts], ".pdf"), wi = 6, he = 6)
  scpp(2.5);
  for (m in 1:nFocus)
  {
    lp.ref = -log10(pMat.ref[, match(focusOn[m], colnames(pMat.ref))]);
    lp.test = -log10(enr[[ts]]$data$pValues[, match(bestModMat[m, ts], colnames(enr[[ts]]$data$pValues))]);

    thinnedScatterplot(as.numeric(lp.ref), as.numeric(lp.test), 
         xlab = spaste("-log enrichment p-value for WGCNA M", focusOn[m]),
         ylab = spaste("-log enrichment p-value for ", testNetNames.pretty[ts], " M", bestModMat[m, ts]),
         cex.lab = 1, cex.axis = 1, cex.main = 1,
         main = spaste("Enrichment significance for ", testNetNames.pretty[ts], " M", bestModMat[m, ts], 
                "\nvs. WGCNA M", focusOn[m], ": "),
         pch = 21, col = "grey25", bg = "grey75", verbose = TRUE, showProgress = TRUE);
  }

  dev.off();
}

# Plot a correlation matrix heatmap of the enrichment p-values for the best modules

for (ts in 1:nTests)
{
  pdf(file = spaste("Plots/enrichmentCorrelationHeatmaps-", testNetNames[ts], ".pdf"), wi = 6, he = 4)
  par(mar = c(8, 8, 2.5, 1));
  lp.ref = -log10(pMat.ref[, match(focusOn, colnames(pMat.ref))]);
  lp.test = -log10(enr[[ts]]$data$pValues[, match(bestModMat[, ts], colnames(enr[[ts]]$data$pValues))]);

  cr1 = cor(lp.ref, lp.test);

  wgHeatmap(cr1, 
     rowLabels.left = spaste("WGCNA M", focusOn),
     rowColors.left = labels2colors(translateUsingTable(focusOn, labelTranslationForColor)),
     colLabels.bottom = spaste(testNetNames.pretty[ts], " M", bestModMat[, ts]),
     colColors.bottom = labels2colors(bestModMat[, ts]),
     textMatrix = round(cr1, 2),
     colors = blueWhiteRed(100),
     zlim = c(-1, 1),
     main = spaste("Correlations of enrichment p-values\n", testNetNames.pretty[ts]),
     cex.main = 1, cex.lab = 0.8);

  dev.off();
}

#=====================================================================================================================
#
# Quick permutation study for the enrichment p-values
#
#=====================================================================================================================

nRuns = 10000;
permCor = numeric(nRuns)
set.seed(1)
pind = initProgInd();
for (r in 1:nRuns) 
{
  permCor[r] = cor(sample(as.numeric(lp.ref)), as.numeric(lp.test));
  pind = updateProgInd(r/nRuns, pind);
  if (r%%200 == 0) gc();
}


