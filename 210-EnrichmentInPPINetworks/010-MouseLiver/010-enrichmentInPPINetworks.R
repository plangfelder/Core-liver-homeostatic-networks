source("../../Functions/networkFunctions-extras-20.R")
source("../../Functions/labelPoints2-01.R");
source("../../Functions/heatmap.wg.R");

source("../../Functions/outlierRemovalFunctions.R")
source("../../Functions/preprocessing-General-013.R")
source("../../Functions/GNVFunctions-015.R")
source("../../Functions/individualAnalysis-General-007-02.R");

#dir.create("RData", recursive = TRUE);
dir.create("Results", recursive = TRUE);
#dir.create("Plots", recursive = TRUE);


library(anRichment)

# Load data

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


nGenes = nrow(refNet)
inSameModule = outer(refNet$module, refNet$module, `==`) + 0;
inModule =  as.numeric(inSameModule * matrix(refNet$module, nGenes, nGenes));

# Read BioGRID PPI network

bgNodes = read.csv(gzfile("../../Data/Annotation/PPI/BioGRID-physical-mouse-4.0.189-nodes.csv.gz"))
bgEdges = read.csv(gzfile("../../Data/Annotation/PPI/BioGRID-physical-mouse-4.0.189-edges.csv.gz"))

dup = duplicated(spaste(bgEdges[[1]], "-", bgEdges[[2]]));
bgEdges = bgEdges[!dup, ];

table(bgNodes$entrez %in% refNet$Entrez)

bgNodes2 = bgNodes[bgNodes$entrez %in% refNet$Entrez, ];
bgEdges2 = bgEdges[ bgEdges$index.1 %in% bgNodes2$index & bgEdges$index.2 %in% bgNodes2$index, ];

bgNodes2$index.expr = match(bgNodes2$entrez, refNet$Entrez)
bgEdges2$index.expr.1 = translateUsingTable(bgEdges2$index.1, bgNodes2[c("index", "index.expr")]);
bgEdges2$index.expr.2 = translateUsingTable(bgEdges2$index.2, bgNodes2[c("index", "index.expr")]);

bgMat = matrix(0, nGenes, nGenes);
bgMat[as.matrix(bgEdges2[c("index.expr.1", "index.expr.2")])] = 1;

sum(bgMat)/length(bgMat)
#[1] 0.0001435566

max(bgMat - t(bgMat)) ## The matrix is not symmetrix.
# [1] 1

sum(bgMat==1 & t(bgMat)==1)/length(bgMat)
# [1] 1.1419e-05
## It's also not approximately symmetric: only less than 1/10 of the non-zero entries have a corresponding transpose
## non-zero.

bgMat.flat = as.numeric(bgMat)

ids = 1:length(inModule)
lst = list(`BioGRID interactions` = ids[as.logical(bgMat)])
bgCollection = collectionFromGeneLists(entrez = lst, organism = "mouse");

gc();


# Run enrichment analysis

moduleLevels = setdiff(sort(unique(refNet$module)), 0);
nModLevels = length(moduleLevels)
overlapSizes = numeric(nModLevels);
modSizes = numeric(nModLevels);
modNGenes = sapply(moduleLevels, function(m) sum(refNet$module==m))
for (m in 1:nModLevels)
{
  printFlush(moduleLevels[m]);
  overlapSizes[m] = sum(inModule== moduleLevels[m] & bgMat.flat==1);
  modSizes[m] = modNGenes[m] * (modNGenes[m] - 1)
  #if (m%%5==0) print(gc());
}

pValues = lpValues = numeric(nModLevels)
n.bg = sum(bgMat.flat==1);
nAll = nGenes * (nGenes-1)

for (m in 1:nModLevels)
{
  pValues[m] = phyper(overlapSizes[m]-1, m = n.bg, n = nAll - n.bg, k = modSizes[m], lower.tail = FALSE, log.p = FALSE) 
  lpValues[m] = phyper(overlapSizes[m]-1, m = n.bg, n = nAll - n.bg, k = modSizes[m], lower.tail = FALSE, log.p = TRUE) 
}

enrTab = data.frame.ncn(Module = moduleLevels, `Module size` = modNGenes,
         `BioGRID ln P-value` = lpValues, `BioGRID p-value` = pValues, 
           #nPotentialInteractions = modSizes, 
           `BioGRID overlap size` = overlapSizes, 
           `BioGRID expected overlap size` = modSizes * n.bg/nAll, 
           `BioGRID enrichment ratio` = overlapSizes/(modSizes * n.bg/nAll))

enrTab.BG = enrTab;

write.csv.nr(signifNumeric(enrTab, 3), file = "Results/enrichmentOfMouseLiverModulesInBioGRIDInteractions.csv");
 
#=======================================================================================================================
#
# Run a similar analysis for STRING
#
#=======================================================================================================================

# Read STRING networkSTRING network

stringNodes = read.csv(gzfile("../../Data/Annotation/PPI/STRING-PPI-mouse-v11.0-nodes.csv.gz"))
stringEdges = read.csv(gzfile("../../Data/Annotation/PPI/STRING-PPI-mouse-v11.0-edges.csv.gz"));

dup = duplicated(spaste(stringEdges[[1]], "-", stringEdges[[2]]));
if (any(dup)) stringEdges = stringEdges[!dup, ];

sum(stringNodes$Entrez %in% refNet$Entrez)

stringNodes2 = stringNodes[stringNodes$Entrez %in% refNet$Entrez, ];
stringEdges2 = stringEdges[ stringEdges$Index.1 %in% stringNodes2$Index & stringEdges$Index.2 %in% stringNodes2$Index, ];

stringNodes2$Index.expr = match(stringNodes2$Entrez, refNet$Entrez)
stringEdges2$Index.expr.1 = translateUsingTable(stringEdges2$Index.1, stringNodes2[c("Index", "Index.expr")]);
stringEdges2$Index.expr.2 = translateUsingTable(stringEdges2$Index.2, stringNodes2[c("Index", "Index.expr")]);

stringMat = matrix(0, nGenes, nGenes);
stringMat[as.matrix(stringEdges2[c("Index.expr.1", "Index.expr.2")])] = stringEdges2$combined_score;

sum(stringMat>0)/length(stringMat)
# [1] 0.0251522

max(stringMat - t(stringMat)) ## The matrix is not symmetrix.
# [1] 756

stringMat.flat = as.numeric(stringMat)

thresholds = c(1, 400, 600, 800);
ids = 1:length(inModule)
lst = lapply(thresholds, function(th) ids[stringMat.flat >= th])
names(lst) = spaste("STRING (confidence >= ", thresholds, ")");
 
stringCollection = collectionFromGeneLists(entrez = lst, organism = "mouse");

gc();

nThresholds = length(thresholds);

# Run enrichment analysis

moduleLevels = setdiff(sort(unique(refNet$module)), 0);
nModLevels = length(moduleLevels)
overlapSizes = listRep(numeric(nModLevels), nThresholds);
modSizes = numeric(nModLevels);
modNGenes = sapply(moduleLevels, function(m) sum(refNet$module==m))
n.bg = numeric(nThresholds);
for (th in 1:nThresholds) 
{
  n.bg[th] = sum(stringMat.flat>=thresholds[th]);
  for (m in 1:nModLevels)
  {
    printFlush(moduleLevels[m]);
    overlapSizes[[th]][m] = sum(inModule== moduleLevels[m] & stringMat.flat >= thresholds[th]);
    modSizes[m] = modNGenes[m] * (modNGenes[m] - 1)
    #if (m%%5==0) print(gc());
  }
  gc();
}

pValues = lpValues = matrix(NA, nModLevels, nThresholds)
nAll = nGenes * (nGenes-1)

for (th in 1:nThresholds) for (m in 1:nModLevels)
{
  pValues[m, th] = phyper(overlapSizes[[th]][m]-1, m = n.bg[th], n = nAll - n.bg[th], k = modSizes[m], 
                          lower.tail = FALSE, log.p = FALSE)
  lpValues[m, th] = phyper(overlapSizes[[th]][m]-1, m = n.bg[th], n = nAll - n.bg[th], k = modSizes[m], 
                          lower.tail = FALSE, log.p = TRUE)
}

colnames(pValues) = spaste("STRING p-value at threshold ", thresholds);
colnames(lpValues) = spaste("STRING log p-value at threshold ", thresholds);
overlapSizes.mat = do.call(cbind, overlapSizes);
colnames(overlapSizes.mat) = spaste("STRING overlap at threshold ", thresholds);
expected.mat = outer(modSizes, n.bg, `*`)/nAll
colnames(expected.mat) = spaste("STRING expected overlap at threshold ", thresholds);
enrichmentRatio.mat = overlapSizes.mat/expected.mat;
colnames(enrichmentRatio.mat) = spaste("STRING enrichment ratio at threshold ", thresholds);


enrTab = data.frame.ncn(Module = moduleLevels, 
           `Module size` = modNGenes, 
   interleave(list(lpValues, pValues, overlapSizes.mat, expected.mat, enrichmentRatio.mat), nameBase = rep("", 5), sep = ""));

dir.create("Results", recursive = TRUE);
write.csv.nr(signifNumeric(enrTab, 3), file = "Results/enrichmentOfMouseLiverModulesInSTRINGInteractions.csv");


