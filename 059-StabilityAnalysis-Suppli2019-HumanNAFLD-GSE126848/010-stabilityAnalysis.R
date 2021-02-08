source("../Functions/networkFunctions-extras-18.R");

x = load("RData/dataForStabilityAnalysis.RData");

#==================================================================================================
#
# Consensus calculation
#
#==================================================================================================

#runID = as.numeric(Sys.getenv("SGE_TASK_ID"))-1;
#if (is.na(runID)) runID = 0;

nRuns.ana = 50;

for (runID in 0:49)
{

printFlush("===============================================================");
printFlush("Working on run", runID + 1);
printFlush("===============================================================");

ana = floor( runID/50) + 1;
run.ana = runID%%50 + 1;
seed = 9 + 2*runID 
#cacheDir = spaste("/u/flashscratch/p/plangfel/temp/", runID);
#cacheDir= spaste("RData/temp/", runID)
#runDir = cacheDir;

#dir.create(runDir, recursive = TRUE);
#dir.create(cacheDir, recursive = TRUE);

#if (FALSE)
#{
#  expr.save = expr.ana;
#  weights.save = weights.ana
#  sample = sample(13000, 1000);
#  expr.ana = lapply(expr.save, mtd.subset, NULL, sample);
#  weights.ana = lapply(weights.save, mtd.subset, NULL, sample);
#}
nAnalyses = length(expr.ana);

nRuns = 1
gc();
print(system.time( {
  shcm = sampledHierarchicalConsensusModules(
    expr.ana[[ana]], 
    multiWeights = weights.ana[[ana]],
    networkOptions = networkOptions.ana[[ana]],
    consensusTree = analysisTrees[[ana]],
    nRuns = nRuns,
    skipUnsampledCalculation = TRUE,
    saveRunningResults = FALSE,
    randomSeed = seed,
    checkSoftPower = TRUE,
  
    detectCutHeight = detectCutHeight,
    maxBlockSize = 30000,
    deepSplit = deepSplit[ana],
    minModuleSize = 20,
    useBranchEigennodeDissim = TRUE,
    mergeCutHeight = 0.2,
    minCoreKME = 0.4, minKMEtoStay = 0.2,
    iteratePruningAndMerging = TRUE,
    cacheDir = cacheDir,
    saveConsensusTOM = FALSE,
    saveIndividualTOMs = FALSE,
    #consensusTOMFilePattern = file.path(runDir, "consensusTOM-Run.%r-%a-block.%b.RData"),
    #individualTOMFilePattern = file.path(runDir, "individualTOM-Run.%r-Set%s-Block%b.RData"),
    useDiskCache = FALSE, verbose = 0);
}));


#save(shcm, file = spaste("RData/sampledHCM-run", prependZeros(runID, 3), ".RData"));

sampledLabels = shcm[[1]]$mods$labels;
save(sampledLabels, file = spaste("RData/sampledLabels-run", prependZeros(runID, 3), ".RData"));
}

stop("All done.");


