# Standard analysis; only difference is that here we use DESeq2 to run the analysis.

source("../Functions/networkFunctions-extras-20.R")
source("../Functions/labelPoints2-01.R");
source("../Functions/outlierRemovalFunctions.R")
source("../Functions/preprocessing-General-013-02.R")
source("../Functions/GNVFunctions-016-02.R")
source("../Functions/individualAnalysis-General-008-03.R");

library(anRichment)
library(DESeq2)

dir.create("RData", recursive = TRUE);
dir.create("Plots", recursive = TRUE);


# Load data

x = load("../110-Preprocessing-TCGA/RData/preprocessedData.RData");

rm(multiExpr)

nSets = length(multiCounts);
tissues = rep("Liver", nSets)
ltissues = tolower(tissues);

# Raise all weights below a threshold to the threshold. 
multiWeights0 = multiWeights;
minWeight = 0.2;
multiWeights = mtd.apply(multiWeights0, function(w) 
{
  w[replaceMissing(w<minWeight)] = minWeight;
  w;
});

topThresholds = list(._FDR.lessThan.0.05_. = list(pAdj = 0.05, p = 1),
                     ._p.lessThan.0.01_. = list(pAdj = 2, p = 0.01));

nTopThresholds = length(topThresholds);


prettifyList = list(
          c("significantMove", "rescueFraction", "statisticallyFullRescue",
            "significantRescue", "significantExacerbation", "differentially.expressed",
            "significantReversal", "significantEnhancement",
            "AllelicSeries", "BACHD.dN17", ".and.", "X5xFAD.TREM2.Plcg2", "X5xFAD.Plcg2",
            "X5xFAD.TREM2", "X5xFAD", "6and10M",
            "._FDR.lessThan.", "._p.lessThan.", "_.", "stat.Q.1", "ContinuousAge", "Age.",
            "inter.", ".with.", "allGNV", ".vs.", ".for.", ".in.", ".by.", "_", "MetabEtiology"),
          c("significant rescue or exacerbation", "rescue fraction", "statistically full rescue",
            "significant rescue", "significant exacerbation", "Differentially expressed",
            "significant reversal", "significant enhancement",
            "Allelic series", "BACHD-deltaN17", " and ", "5xFAD/TREM2/Plcg2", "5xFAD/Plcg2",
            "5xFAD/TREM2", "5xFAD", "6 and 10M",
            " (FDR < ", " (p < ", ")", "AS Z for Q1", "Continuous age", "Age ",
            "interaction of ", " with ", "all GNV", " vs. ", " for ", " in ", " by ", 
            " ", "metabolic etiology"));

prettifyList.plots = list(
          c("FDR < ", "p < "),
          c("FDR\U202F<\U202F", "p\U202F<\U202F"));


topThresholdNames.pretty = prettifyStrings(names(topThresholds), prettifyList)
topThresholdNames.plot = prettifyStrings(topThresholdNames.pretty, prettifyList.plots)

#=============================================================================================================
#
# Create designs and contrasts
#
#=============================================================================================================

# Retain columns from sample annotation that contain no more than 3/4 missing data.

minFractionPresent = 0.4;
multiSampleAnnot.present = mtd.apply(multiSampleAnnot, function(sa)
{
  nr = nrow(sa);
  num = sapply(sa, is.numeric);
  sa = as.data.frame(lapply(sa, function(x) if (is.numeric(x)) x else {x[x==""] = NA; x}), check.names = FALSE);
  nUnique = sapply(sa, function(x) length(unique(x[!is.na(x)])));
  keep = sapply(sa, function(x) sum(!is.na(x)) >= nr*minFractionPresent) & (nUnique > 1) & ( (nUnique < nrow(sa)) | num);
  sa = sa[, keep];
  # Special preprocessing of the trait Etiology, if present
  if ("Etiology" %in% names(sa))
  {
    unclear = is.na(sa$Etiology) | replaceMissing(sa$Etiology=="W/O Etiology") | grepl(".+ METAB", replaceMissing(sa$Etiology));
    metab = replaceMissing(sa$Etiology=="METAB");
    metab[unclear] = NA;
    sa$MetabEtiology = 0 + metab;
  }
  sa;
});

multiSampleAnnot.present.machine = mtd.apply(multiSampleAnnot.present, function(sa)
{
  sa = lapply(sa, function(x) if (is.numeric(x)) x else make.names(x));
  as.data.frame(setNames(sa, make.names(names(sa))))
})

#=============================================================================================================
#
# Decide on what to use as covariates. 
#
#=============================================================================================================

multiSampleAnnot.present.machine2 = mtd.apply(multiSampleAnnot.present.machine, function(sa)
{
  if (!"Center.of.sequencing" %in% names(sa)) return(sa);
  extras = lapply(sa[c("Center.of.sequencing", "Race.Category")], function(x)
  {
    x2 = x;
    t = table(x);
    n1 = t[x];
    x2[n1 < 5] = NA;
    x2[x=="NA."] = NA;
    x2;
  });
  names(extras) = spaste(c("Center.of.sequencing", "Race.Category"), "2");
  data.frame(sa, do.call(data.frame, extras));
});


candidateCovars = c("donor_sex", "donor_age_at_diagnosis", "Center.of.sequencing2", "Race.Category2");

prettifyList2 = mymapply(c, prettifyList, list(c("sequencing2", "Category2"), c("sequencing", "category")));

candidateCovars.machine = make.names(candidateCovars); 
ncc = length(candidateCovars);
covarDesign = mtd.apply(multiSampleAnnot.present.machine2, function(sa)
{
   covars1 = intersect(candidateCovars, names(sa));
   num = covars1[sapply(sa[covars1], is.numeric)];
   cat = setdiff(covars1, num);
   levels = lapply(sa[cat], function(x) levels(factor(x)));
   lastLevel = sapply(levels, function(x) x[length(x)]);
   firstLevel = sapply(levels, function(x) x[1]);
   covars2 = c(num, cat);
   data.frame(
     design = spaste("~", covars2),
     reduced = "~1",
     coefficientName = c(num, spaste(cat, "_", lastLevel, "_vs_", firstLevel)),
     contrast.name = NA, contrast.numerator = NA, contrast.denominator = NA,
     testName = covars2,
     testName.pretty = covars2,
     combineSignificantHits = FALSE,
     splitSignificantHitsBy = NA);
});

combineResults = mtd.apply(covarDesign, function(di) list(all = di$testName))

analyses.candidateCovars = list();

nSets = length(multiCounts);

for (set in 1:nSets)
 analyses.candidateCovars[[set]] = list(data = individualAssociation.general(
    data = multiCounts[[set]]$data,
    dataWeights = multiWeights[[set]]$data,
    pheno = multiSampleAnnot.present.machine2[[set]]$data,
    analysisName = spaste(setNames, "-candidateCovars")[set],
    designInfo = covarDesign[[set]]$data,
    combineResults = combineResults[[set]]$data,
    geneAnnotation.start = multiGeneAnnot[[set]]$data,
      analysisPipeline = "DESeq2",
      addToAnalysisManager = FALSE,
      combineColNameExtensions = "",
      resultDir = file.path("Results"),
      stopAt = "40-Gene lists",
      createCollectionOfSignificantGenes = FALSE,
      topThresholds = topThresholds,
      separateExtendedResults = FALSE,
      topList.minN = 0,
      prettifyList = prettifyList,
      extendPrettifyList = TRUE,
      keepFullDESeqObjects = FALSE,
      topList.machineColName = "GeneID",
      topList.humanColName = "Symbol",
      topList.backupColName = "GeneID",
      dropRepeated.exceptions = c("changePAdjTo", "F.df"),
      forceExport = TRUE,
      organism = "human",
      internalClassification = setNames,
      indent = 2, verbose = 2));

nSignif = mtd.apply(analyses.candidateCovars, getElement, "nSignif.all")

nSignif.FDR = mtd.apply(nSignif, getElement, "._FDR.lessThan.0.05_.");

save(nSignif, nSignif.FDR, file = "RData/nSignif-candidateCovars.RData");

table(multiSampleAnnot.present.machine[[2]]$data$Center.of.sequencing)

# Plot the numbers of DE genes and move on.

nSignif.all = mtd.rbindSelf(nSignif.FDR);
col.right = "#FF9070"; col.left = "skyblue";
border.right = "brown"; border.left = "blue";

separators = cumsum(unlist(mtd.apply(nSignif.FDR, nrow)));
separators2 = c(0, separators[-nSets]);

tt = 1;

library(Cairo)

CairoPDF(file = "Plots/candidateCovars-nSignificant.all.pdf", wi = 7, he = 3);
par(mar = c(0.2, 15, 2, 0.2));
    plot = twoSideBarplot(nSignif.all[, 1], nSignif.all[, 2],
    col.right = col.right, col.left = col.left,
    border.right = border.right, border.left = border.left,
    yLabels = prettifyStrings(rownames(nSignif.all), prettifyList2),
    barGap = 0.2,
    main = spaste("Numbers of significant genes", topThresholdNames.plot[tt], "           "),
    cex.main = 1.2,
    separatorPositions = separators, sep.col = "grey", sep.ext = TRUE)

    text(rep(plot$box[1]-0.98*plot$leftMargin, length(separators)), plot$ytop[separators2+1],
           setNames, adj = c(0,1), cex = 0.9, xpd = TRUE, font = 2)

dev.off();

# Adjust for age and sex. In principle could also adjust for race except that many of the samples have NA for race.
# It is not quite clear how come none of the genes are significant for center of sequencing, even after removing the
# low-frequency ones...

#=============================================================================================================
#
# Create designs and contrasts
#
#=============================================================================================================

selectedTraits = c("specimen_type", "specimen_interval", "MetabEtiology", "Disease Stage", "Liver fibrosis ishak score category", 
       "Tumor Size", "Vascular Invasion",
       "Cancer Type", "Cancer Type Detailed", "Numeric liver fibrosis ishak score category");

traitsAsCovars = c("donor_sex", "donor_age_at_diagnosis");

xx = mtd.apply(multiSampleAnnot.present, function(sa)
{
  printFlush("================================================================================");
  lapply(selectedTraits, function(tr) {printFlush(tr, ":"); print(table(sa[[tr]]))})
})

minCountAtLevel = 10;
# Binarize/preprocess selected traits
multiSampleAnnot.selected= mtd.apply(multiSampleAnnot.present, function(sa)
{
  if ("Liver fibrosis ishak score category" %in% names(sa))
  {
    sa$`Numeric liver fibrosis ishak score category` = as.numeric(factor(sa$`Liver fibrosis ishak score category`));
  }
  present = selectedTraits[selectedTraits %in% names(sa)];
  prettifyList = list(make.names(present), present);
  categorical = present[!sapply(sa[present], is.numeric)]; 
  num = setdiff(present, categorical);
  out = sa[num];
  if (length(categorical) > 0) 
  {
    sa.bin = binarizeCategoricalColumns(sa[categorical], 
               convertColumns = categorical, includePrefix = TRUE, checkNames = TRUE, minCount = minCountAtLevel,
               includePairwise = TRUE, includeLevelVsAll = FALSE);
    for (tr in categorical) 
    {
      levels1 = setdiff(unique(sa[[tr]]), NA)
      prettifyList = mymapply(c, prettifyList, list(make.names(levels1), levels1));
    }
    out = data.frame(out, sa.bin);
  }
  attr(out, "prettifyList") = prettifyList;
  out;
})

prettifyList3 = lapply(1:2, function(i) unlist(lapply(  
        mtd.apply(multiSampleAnnot.selected, attr, "prettifyList", returnList = TRUE), `[[`, i)));
  

multiSampleAnnot.covars = mtd.apply(multiSampleAnnot, function(sa)
{
  out = sa["donor_age_at_diagnosis"];
  out$M.vs.F = as.numeric(factor(sa$donor_sex))-1;
  out;
});

multiSampleAnnot.otherTraits = mtd.apply(multiSampleAnnot.present, function(sa)
{
  n0 = setdiff(names(sa), c(selectedTraits, traitsAsCovars));
  n1 = multiGrep(c("_sex$", "_id$", "donor_vital_status", 
                   "disease_status_last_followup", "_age_", "PrincipalComponent", "Sex", "donor_interval_of_last_followup",
                   "donor_age_at_enrollment", "Form.completion.date", 
                    "Number of Samples Per Patient", "Birth from Initial Pathologic Diagnosis Date", 
                 "specimen_", "ID$", "Age$", "Code$", "Sample Type", "Race.Category", "Center.of.sequencing",
                 "donor_tumour_stage_at_diagnosis", "donor_survival_time", 
                 "American.Joint.Committee.on.Cancer.Publication.Version.Type",
                 "New.Neoplasm.Event.Post.Initial.Therapy.Indicator", "Informed.consent.verified",
                 "Ethnicity.Category", "Person.Gender", "Tissue.Source.Site_9", 
                 "Tissue Prospective Collection Indicator_7", "Tissue Retrospective Collection Indicator_8",
                 "In PanCan Pathway Analysis"), n0, value = TRUE, invert = TRUE);
  printFlush("Dropped:", formatLabels(paste(setdiff(n0, n1), collapse = "; "), 80), "\n")
  sa = sa[n1];
  num = n1[sapply(sa, is.numeric)];
  categorical = setdiff(n1, num);
  prettifyList = list(make.names(n1), n1);
  out = sa[num];
  if (length(categorical) > 0)
  { 
    sa.bin = binarizeCategoricalColumns(sa[categorical],           
               convertColumns = categorical, includePrefix = TRUE, checkNames = TRUE, minCount = minCountAtLevel,
               includePairwise = TRUE, includeLevelVsAll = FALSE);
    for (tr in categorical)
    { 
      levels1 = setdiff(unique(sa[[tr]]), NA)
      prettifyList = mymapply(c, prettifyList, list(make.names(levels1), levels1));
    }
    sa.bin = sa.bin[, grep("Status.Not.Called", names(sa.bin), invert = TRUE)];
    out = data.frame(out, sa.bin);
  }
  attr(out, "prettifyList") = prettifyList;
  out;
})

prettifyList4 = lapply(1:2, function(i) unlist(lapply(
        mtd.apply(multiSampleAnnot.otherTraits, attr, "prettifyList", returnList = TRUE), `[[`, i)));

prettifyList5 = as.data.frame(mymapply(c, prettifyList2, prettifyList3, prettifyList4))
prettifyList5 = prettifyList5[!duplicated(prettifyList5[[1]]), ];
prettifyList5 = prettifyList5[order(-nchar(prettifyList5[[1]])), ];



covarNames = mtd.colnames(multiSampleAnnot.covars);
covarString = prependPrefix(" + ", paste(covarNames, collapse = " + "))

designInfo.selectedTraits = mtd.apply(multiSampleAnnot.selected, function(sa)
{
  traits = colnames(sa);
  if (length(traits) > 0)
  {
    data.frame(
       design = spaste("~", traits, covarString),
       coefficientName = traits,
       contrast.name = NA, contrast.numerator = NA, contrast.denominator = NA,
       testName = traits,
       testName.pretty = prettifyStrings(traits, prettifyList5),
       combineSignificantHits = FALSE,
       splitSignificantHitsBy = NA);
  } else NULL;
   
});

designInfo.otherTraits = mtd.apply(multiSampleAnnot.otherTraits, function(sa)
{
  traits = colnames(sa);
  if (length(traits) > 0)
  {
    data.frame(
       design = spaste("~", traits, covarString),
       coefficientName = traits,
       contrast.name = NA, contrast.numerator = NA, contrast.denominator = NA,
       testName = traits,
       testName.pretty = prettifyStrings(traits, prettifyList5),
       combineSignificantHits = FALSE,
       splitSignificantHitsBy = NA);
  } else NULL;

});


# Leave this analysis as is for now. Just save all sample annotation and design info versions for use in network analysis.
designInfo.covars = covarDesign;

save(multiSampleAnnot.present, multiSampleAnnot.present.machine, multiSampleAnnot.present.machine2,
     multiSampleAnnot.selected, multiSampleAnnot.otherTraits, multiSampleAnnot.covars,
     designInfo.covars, designInfo.selectedTraits, designInfo.otherTraits,  
     covarNames, covarString, prettifyList5, file = "RData/sampleAnnotation-designInfo-prettifyList.RData");

#===================================================================================================================
#
# Run the DE analysis on selected and other traits
#
#===================================================================================================================

x = load("../050-IndividualAnalysis-Suppli2019-HumanNAFLD-GSE126848/RData/stdCollections.RData");

mouseColl = loadAsList("../020-IndividualAnalysis/RData/collectionOfDEGenes-Liver.RData")[[1]];

allCollections = c(stdCollections, list(mouseLiver = convertCollectionToOrganism(mouseColl, organism = "human")))



traitSetNames = c("selected", "other");
designInfo.run = list(designInfo.selectedTraits, designInfo.otherTraits);
names(designInfo.run) = traitSetNames;
nTraitSets = length(traitSetNames);

multiSampleAnnot.run = mtd.mapply(cbind, multiSampleAnnot.covars, multiSampleAnnot.otherTraits, multiSampleAnnot.selected);
combineResults.run = lapply(designInfo.run, mtd.apply, function(di) list(all = di$testName))

analyses = listRep(list(), nTraitSets);
for (traitSet in 1:nTraitSets) for (set in 1:nSets) if (length(designInfo.run[[traitSet]][[set]]$data) > 0)
 analyses[[traitSet]] [[set]] = list(data = individualAssociation.general(
    data = multiCounts[[set]]$data,
    dataWeights = multiWeights[[set]]$data,
    intermediateResults = if (length(analyses[[traitSet]]) < set) NULL else analyses[[traitSet]][[set]]$data,
    pheno = multiSampleAnnot.run[[set]]$data,
    analysisName = spaste(setNames[set], "-", traitSetNames[traitSet]),
    designInfo = designInfo.run[[traitSet]] [[set]]$data,
    combineResults = combineResults.run[[traitSet]] [[set]]$data,
    geneAnnotation.start = multiGeneAnnot[[set]]$data,
      analysisPipeline = "DESeq2",
      addToAnalysisManager = FALSE,
      combineColNameExtensions = "",
      resultDir = file.path("Results"),
      dropTopTableColumns = c("t", "B", "genelist"),
      #stopAt = "40-Gene lists",
      createCollectionOfSignificantGenes = FALSE,
      topThresholds = topThresholds,
      separateExtendedResults = FALSE,
      topList.minN = 0,
      enrichment.minN = 0,
      enrichment.collections = allCollections,
      prettifyList = prettifyList,
      extendPrettifyList = TRUE,
      keepFullDESeqObjects = FALSE,
      topList.machineColName = "GeneID",
      topList.humanColName = "Symbol",
      topList.backupColName = "GeneID",
      forceExport = TRUE,
      organism = "human",
      internalClassification = setNames,
      indent = 2, verbose = 2));

save(analyses, file = "RData/analyses.RData");

nSignif = lapply(analyses, mtd.apply, getElement, "nSignif.all")

nSignif.FDR = lapply(nSignif, mtd.apply, getElement, "._FDR.lessThan.0.05_.");

save(nSignif, nSignif.FDR, file = "RData/nSignif.RData");

# Plot the numbers of DE genes 

nSignif.all = mtd.rbindSelf(nSignif.FDR);
col.right = "#FF9070"; col.left = "skyblue";
border.right = "brown"; border.left = "blue";

tt = 1;

library(Cairo)

for (traitSet in 1:nTraitSets) for (set in 1:length(nSignif.FDR[[traitSet]])) 
{
  ns1 = nSignif.FDR[[traitSet]] [[set]]$data;
  nr = nrow(ns1);
  labs = prettifyStrings(rownames(ns1), prettifyList2)
  mw = marginWidth(labs, device = "CairoPDF");
  CairoPDF(file = spaste("Plots/nSignificant-", setNames[set], "-", traitSetNames[traitSet], ".pdf"), 
           wi = 4 + mw$inch, he = 0.5 + 0.3 * nr);
  par(mar = c(0.2, mw$lines+0.2, 2, 0.2));
  plot = twoSideBarplot(ns1[, 1], ns1[, 2],
      col.right = col.right, col.left = col.left,
      border.right = border.right, border.left = border.left,
      yLabels = labs,
      barGap = 0.2,
      main = spaste("Numbers of signif. genes in ", setNames[set], topThresholdNames.plot[tt], "                           "),
      cex.main = 1);

  dev.off();
}

#===================================================================================================================
#
# prediction analysis for metabolic etiology: this is carried out in the network analysis code since that code has the
# appropriate (vst) data.
#
#===================================================================================================================


