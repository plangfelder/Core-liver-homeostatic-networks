# Preprocessing of general expression data. Inspired by the RNA-seq preprocessing code, the DESeq2 variance
# stabilizing transformation is replaced by vsn from package vsn.

.preprocessing.General.version = "14.00";


### FIXME: when collapsing data, add the collapsed gene annotation and its file name to the output

# general input: 
#    . expression data in columns-are-genes or columns-are-samples format
#    . phenotypes, split into main groups of covariates and phenotypes of interest
#    . Possibly a version of phenotypes for plotting sample trees

# ideal pre-processing pipeline:
# . Adjust for confounders using weights that correctly take into account the design of the experiment
# . Remove outliers. To avoid biasing association studies, outlier removal should NOT be done in groups of 
#   phenotypes/covariates that were left over after adjustment. 
# . Optionally quantile normalize

# Version 002: 
#   .main change is that phenoIfIneterst and covariates are now merged into a single input data
#        frame/matrix. 

# Version 003:
#   . Added many additions
#   . Removal of PCs: at present the PCs are computed only once, on vst data before outlier removal.

# Version 003.01
#   . Fix sample tree plot file names.
#   . Added optional saving of individual outlier removal weights and weight factors, for now only for xdata.

# Version 004
#   . Changed the default directory for saving weights to exprDir.OSR
#   . Added weights for collapsed data

# Version 005
#   . For sample clustering and outlier removal, add the option to scale data first and make default value
#     TRUE, with restriction to top 8000 genes.

# Version 006
#   . sample tree 011 is drawn by clustering adjusted data (the same data that is used for outlier
#     identification)
#   . Change within 006: fixed bug where OSR details were not recorded.

# Version 007
#   . The restriction to HE columns, VST, adjustment and outlier removal are iterated until no further
#      outliers are removed. Most often this is not necessary but there have been cases where one removes 
#     a few outliers and the set of good genes changes a lot.
#   . The code should work without any reference to analysis manager when addToAnalysisManager is FALSE.
#   . PCs from filtered data are optionally added to sampleAnnotation and pheno in the final step where there
#     are no outliers.

# Version 008
#   . Adding optional RUV-generated factors as covariates

# Version 009
#   . Changing heColumns to by default scale all samples to common mean (which is not the same as sum when
#     missing data are present); for RNA-seq (count) pipeline also scale to counts per million for more
#     consistent interpretation. Also changing defaults for heColumns to be far more restrictive than before.

# Version 010
#   . Distance calculation for outlier removal now scales data as well, controlled by or.scaleData.
#   . Addign optional calculation of SVA factors and the option to use them for adjustment of data. SVA
#     factors and RUV factors are only calculated from the top probes.

# Version 011
#   . RUV and SVA factors are now by default calculated on scaled data; the behaviour is controlled by
#     argument scaleForSVD (that also controls scaling for PC calculation)

# Version 012
#   . Uses networkFunctions-extras-14.R which changed how covariates are handled in calculations of bicov
#     weights. Corresponding changes include removal of arguments bw.combineGroups and bw.minSamplesPerGroup.

# Version 013
#   . Adding the ability to do oulier removal without SVA/RUV factor adjustment. The reason is that presence
#     of outliers will cause the first PC/SVA/RUV factor to correct for the outlier, which makes no sense. In
#     fact, adjustment for PCs/SVA/RUV should be done before outlier removal only if a whole cluster of
#     samples gets called "outliers", preferably with several genotypes present in the cluster.
#   . Removing the distinction between pheno, plotPheno and sampleAnnotation; only sampleAnnotation will be used. 
#   . Adding optional PCA plots every time the sample tree is plotted.  
#   . PCA and sample tree plots are produced irrespective of whether the calculations were done or contained in
#     intermediate results
#   . PCA plots use border to indicate outlier status
#   . Argument amComponent can be used to specify other than "preprocessing", e.g., "preprocessing with 1 SVA" etc.
#   . Adding optional argument adj.fitToSamples that gets passed to empiricalBayesLM as fitToSamples. 
#   . New argument iorFromHEGenes, default FALSE for backward compatibility, can be used to restrict the IOR calculation
#     to just the HE genes. The reason this was added is that the variance stabilization used in calculation of IOR
#     weights depends on what all genes do, and sometimes the difference could be relatively large. 
#     Note that this is not tested, and may not work properly, for IOR calculations on collapsed data. 
#   . nTopProbes is adjusted if it is more than the number of variables left after HE filtering
#   . pca.colorColumn, pca.shapeColumn, pca.colorPrefix, pca.shapePrefix can now be vectors; one PCA plot is generated for
#     each entry.
#   . Input sample annotation and gene annotation filtered to "valid Entrez" is always saved, even if saveXData.entrez is
#     FALSE.

# Version 013-02: 
#   . bugfix fixes saveing objects into object manager with incorrect analysis component part of the name
#     vector.
#   . Adding output of more file and folder information, especially for IOR

# Version 014:
#   . Changing structure of output by combining all file names into a single component Files
#   . Changing defaults for directories for saving collapsed data, now they are separate from 
#     their OSR, IORep and IORem counterparts
#   . EBLM adjustment can be done in several steps.

histogramReferencePoints = function(expr, breaks = 200, log = TRUE)
{
  h = hist(as.matrix(expr), breaks = breaks, plot = FALSE);
  if (log) h$counts = log10(h$counts + 1);
  loc.min = min(h$mids[h$counts > 0]);
  loc.max = h$mids[which.max(h$counts)];
  loc.maxDiff = h$mids[which.max( h$counts[-1] - h$counts[-length(h$counts)])];
  list(loc.min = loc.min, loc.max = loc.max, loc.maxDiff = loc.maxDiff,
       hist = h);
}

heColumns = function(x, minValue, minProportion, useCountsPerMillion = FALSE, scaleSamples = TRUE, 
                     verbose = 1, ...)
{
  ns = nrow(x);
  if (scaleSamples)
  {
    if (any(x<0, na.rm = TRUE))
      warning(immediate. = TRUE, 
           "heColumns: scaling samples by sums may not make sense in the presence \n",
           "of negative data. Proceed with caution or stop and rerun with scaleSamples = FALSE.");
    sums = rowSums(x, na.rm = TRUE);
    meanSum = mean(sums, na.rm = TRUE);
    factors = (if (useCountsPerMillion) 1e6 else meanSum) /sums;
    factors[!is.finite(factors)] = 0;
    if (any(factors==0)) 
      stop("heColumns: have non-finite or zero scaling factors. Please check data and re-run.");

    x = x * factors;
  }
  out = goodGenes(x, minNSamples = min(WGCNA:::..minNSamples, nrow(x)), 
                  minNGenes = min(WGCNA:::..minNGenes, ncol(x)), verbose = 0) & 
        colSums(x >= minValue, na.rm = TRUE) >= minProportion * ns;
  if (verbose) printFlush("heColumns: keeping", sum(out), "variables.");
  out;
}


# Convert column names from gene symbols to Entrez.

convertColumnNames2entrez = function(x, organism, transpose = FALSE,
                  sourceSymbols = if (transpose) rownames(x) else colnames(x))
{
  entrez0 = convert2entrez(symbol = sourceSymbols, organism = "Mouse");
  keepGenes = !is.na(entrez0);
  if (transpose) x = t(x);
  x1 = x[, keepGenes];
  entrez1 = entrez0[keepGenes];
  t = table(entrez1);
  # Some genes may be duplicated
  keep2 = names(t)[t==1];
  x2 = x1[,match(keep2, entrez1)];
  entrez = keep2;
  colnames(x2) = entrez;
  attr(x2, "indexInOriginalData") = match(keep2, entrez0);
  x2;
}


writeSessionInfo = function(file)
{
  sink(file)
  sessionInfo();
  printFlush(spaste("RNA-seq preprocessing functions version", .preprocessing.General.version));
  sink(NULL)
}

imputeAndScaleData = function(expr, scale= TRUE)
{
  if (length(expr)==0) return(NULL);
  if (sum(is.na(expr)) > 0) {
     expr.imp = t(impute.knn(t(expr))$data);
  } else
      expr.imp = expr;
  if (scale) expr.imp = scale(expr.imp);
  expr.imp;
}


calculatePCs = function(expr, nPCs, scale = TRUE, PCColnamePrefix = "principalComponent.")
{
   if (length(expr)==0 || nPCs==0) return(NULL);
   expr.simp = imputeAndScaleData(expr, scale = scale);
   printFlush("   ..calculating SVD..")
   svd = svd(expr.simp, nu = nPCs, nv = 0);
   printFlush("     ..done.");
   PCs = setColnames(svd$u, spaste(PCColnamePrefix, 1:nPCs));
   rownames(PCs) = rownames(expr);
   PCs;
}

conditionalWriteAndAddToAnalysis = function(
  doSave, 
  object, 
  dir, fileName, quote = TRUE, row.names = FALSE,
  fileReadArgs = defaultReadArguments("csv"),
  amAnalysis, amName, amComponent, allowNonstandardObjectType = TRUE,
  verbose = 0, indent = 0)
{
   if (doSave && length(object) > 0)
   {
      dir.create(dir, recursive = TRUE, showWarnings = FALSE);
      write.csv(object, file = file.path(dir, fileName), quote = quote, row.names = row.names)
      if (!is.null(amAnalysis)) amAnalysis = am.addNewObjectToAnalysis(
        name = amName,
        fileName = file.path(dir, fileName),
        anaType = amComponent,
        fileReadArgs = fileReadArgs,
        analysis = amAnalysis,
        allowNonstandardObjectType = allowNonstandardObjectType, 
        allowNonstandardAnalysisType = TRUE,
        returnBoth = FALSE,
        verbose = verbose - 1, indent = indent);
   }
   amAnalysis;
}


adjustCovariates.EBLM = function(
   x,
   covariates,
   PCs = NULL,
   removedCovariates = NULL,
   retainedCovariates = NULL,
   fitToSamples = NULL,
   weights = NULL,
   useBicovWeights = TRUE,
   bw.groupBy = NULL,
   bw.groupsForMinWeightRestriction = NULL,
   bw.otherArgs = list(maxPOutliers = 0.10,
                       outlierReferenceWeight = 0.1,
                       minWeightInGroups = 0.9,
                       maxPropUnderMinWeight = 0.4,
                       defaultWeight = 1),
   method = c("EB", "OLS"),
   verbose = 0, indent = 0
)
{
  spaces = indentSpaces(indent);
    
  if (length(removedCovariates)==0) return(x);
  if (is.null(x)) return(x);

  method = match.arg(method);
  useEB = method=="EB";
  multiStep = is.list(removedCovariates);
  #if (multiStep && !is.list(removedCovariates)) 
  #  stop("When adjustment is multi-step, 'removedCovariates' must be a list.");
  if (multiStep && length(retainedCovariates) > 0 && !is.list(retainedCovariates))
    stop("When adjustment is multi-step, 'retainedCovariates' must be length-zero or a list.");
  if (multiStep && length(fitToSamples) > 0 && !is.list(fitToSamples))
    stop("When adjustment is multi-step, 'fitToSamples' must be length-zero or a list.");

  if (multiStep)
  { 
    nSteps = length(removedCovariates);
    if (length(retainedCovariates)>0 && length(retainedCovariates)!=nSteps)
     stop("When adjustment is multi-step, length of 'retainedCovariates' must be zero or ",
          "equal length of 'removedCovariates'.");
    if (length(fitToSamples)>0 && length(fitToSamples)!=nSteps)
     stop("When adjustment is multi-step, length of 'fitToSamples' must be zero or ",
          "equal length of 'removedCovariates'.");
  } else nSteps = 1;

  if (!is.list(removedCovariates)) removedCovariates = list(removedCovariates);
  if (length(retainedCovariates) > 0 && !is.list(retainedCovariates)) retainedCovariates = list(retainedCovariates);
  if (length(fitToSamples) > 0 && !is.list(fitToSamples)) fitToSamples = list(fitToSamples);

  for (step in 1:nSteps)
  {
    if (length(removedCovariates[[step]])==0) next;
    if (verbose > 0 && multiStep) printFlush(spaste(spaces, " adjustment step ", step));
    if (step==nSteps && !is.null(PCs))
    {
      removedCovariates[[step]] = c(removedCovariates[[step]], colnames(PCs));
      if (verbose > 0 && multiStep) printFlush(spaste(spaces, " ..with latent factors"));
      if (length(covariates) > 0) {
         covariates = cbind(covariates, PCs);
      } else 
         covariates = PCs;
    }
    removedIndex = match(removedCovariates[[step]], colnames(covariates));
    if (any(is.na(removedIndex)))
      stop("The following entries in 'removedCovariates' were\n ",
           "not found among colnames(covariates):\n",
           paste(removedCovariates[[step]][is.na(removedIndex)], collapse = ", "))

    if (length(retainedCovariates)>0 && length(retainedCovariates[[step]]) > 0)
    {
      retainedIndex = match(retainedCovariates[[step]], colnames(covariates));
      if (any(is.na(retainedIndex)))
        stop("The following entries in 'retainedCovariates' were\n ",
             "not found among colnames(covariates):\n",
             paste(retainedCovariates[is.na(retainedIndex)], collapse = ", "))
      retainedCovariates1 = covariates[, retainedIndex];
    } else 
      retainedCovariates1 = NULL;

    if (is.null(weights))
    {
      if (useBicovWeights)
      {
        weights = do.call(matrixBicovWeights,
                    c(list(data = x, covars = covariates,
                           groupBy = bw.groupBy,
                           groupsForMinWeightRestriction = bw.groupsForMinWeightRestriction),
                      bw.otherArgs));
      }
    }
    ebData = empiricalBayesLM(
      x,
      removedCovariates = covariates[, removedIndex],
      retainedCovariates = retainedCovariates1,
      fitToSamples = fitToSamples[[step]],
      weights = weights,
      weightType = "empirical",
      getOLSAdjustedData = !useEB,
      getResiduals = FALSE,
      getFittedValues = FALSE,
      getWeights = FALSE,
      getEBadjustedData = useEB,
      verbose = verbose, indent = indent);

    if (useEB)
    {
      x = ebData$adjustedData
    } else
      x = ebData$adjustedData.OLS
  }
  x;
}              
                  

allCollapseSteps = function()
{
  c("OSR", "IOReplaced", "IORemoved", "IORep.adj", "IORem.adj")
}


#======================================================================================================
#
# The first steps of preprocessing, separated out into a separate function because they are called
# from a couple different points in the workflow.
#
#======================================================================================================

### Important note: if principal component removal is carried out and the PCs are saved and recorded in the
### analysis, each call to HE.VST.Adjust overwrites the previous PCs; at the end, only the last set of PCs
### remain.

addNamedList = function(master, add, name)
{
  addList = list(add);
  names(addList) = name;
  c(master, addList);
}


HE.VST.Adjust = function(
  xdata,
  pheno,
  geneAnnotation,
  intermediateResults,

  normalizationPipeline = NULL,
  forceCalculation,
  stopAt = "End",
  amAnalysis,
  amComponent,

  # HE filter
  ir.HEFilter = "HEFilter",
  heColumnFnc,
  minValue,
  minProportion,
  he.useCPM,
  heColumnOtherArgs,

  # Data normalization and transformation
  ir.VST = "VST",
  # Variance stabilization
  useVarianceStabilization = TRUE,   # only useful for general data - VST is mandatory in Counts pipeline

  # Arguments for count data
  vst.design,
  vst.blind = FALSE,
  vst.calibrateLowestValue = TRUE,

  # Arguments for general data
  vst.calib = "affine",
  vst.keepModel = TRUE,

  #Arguments for approximate VST. Can be used within Counts as well as general pipeline.
  useApproximateVST = FALSE,
  approximateVST.minShift = 1,
  approximateVST.maxShift = 500,
  approximateVST.xWeightScale = 0,
  # Caution: if input data are not log transformed, use approximateVST.existingLogBase = NULL
  approximateVST.existingLogBase = 2,
  approximateVST.existingLogShift = 1,

  # Before outlier removal: should VST data be adjusted for covariates?
  # The covariates below must be names of the columns in pheno above so they can be
  # used for adjustment of IOMR data.
  # If removedCovariates is of length 0, adjustment is skipped.
  ir.VST.adj = "VST.adj",
  adj.removedCovariates = NULL,
  adj.retainedCovariates = NULL,
  adj.fitToSamples = NULL,
  adj.useMethod = c("EB", "OLS"),
  adj.useBicovWeights = TRUE,
  adj.nProbesForFactorAnalysis = -1,
  adj.calculateNPCs = 0,
  adj.removeNPCs = 0,

  PCColnamePrefix = "PrincipalComponent.",
  savePCs,
  PCDir,
  PCFile,
  amNameForPCs,
  scaleForSVD = TRUE,

  adj.RUV.useNFactors = 0,
  adj.RUV.calculateNFactors = 0,
  adj.RUV.groupVariables = NULL,
  RUVColnamePrefix = "RUVFactor.",
  RUVFactorFile,
  amNameForRUVFactors,

  adj.calculateSVAFactors = FALSE,
  adj.sva.useNFactors = 0,
  adj.svaModel,
  adj.svaReducedModel = NULL,
  svaColnamePrefix = "SVAFactor.",
  svaFactorFile,
  amNameForSVAFactors,

  bw.groupBy = NULL,
  bw.groupsForMinWeightRestriction = NULL,
  bw.otherArgs = list(maxPOutliers = 0.10,
                      outlierReferenceWeight = 0.1,
                      minWeightInGroups = 0.9,
                      maxPropUnderMinWeight = 0.4,
                      defaultWeight = 1),
  verbose = 1, indent = 0)

{
   spaces = indentSpaces(indent);
   recognizedPipelines = c("Counts", "General");
   normPipe = charmatch(tolower(normalizationPipeline), tolower(recognizedPipelines));

   if (is.na(normPipe))
       stop("If given, 'normalizationPipeline' must be one of\n   ",
            paste(recognizedPipelines, collapse = ", "));

   #----------------------------------------------------------------------------------------------
   # Check for non-varying and low-expressed genes
   #----------------------------------------------------------------------------------------------
   HEFilter = getElement(intermediateResults, ir.HEFilter);
   if (length(HEFilter)==0 || forceCalculation)
   {
     HEFilter = list();
     if (verbose > 0) printFlush(spaste(spaces, "Filtering out low-expressed genes..."))
     heCols = do.call(heColumnFnc, 
                      c(list(x = xdata, minValue = minValue, minProportion = minProportion,
                             useCountsPerMillion = he.useCPM, verbose = verbose -1),
                        heColumnOtherArgs));
     xdata.he = xdata[, heCols];
     geneAnnotation.he = geneAnnotation[heCols, ];
     HEFilter$xdata.he = xdata.he;
     HEFilter$geneAnnotation.he = geneAnnotation.he;
     HEFilter$heColumns = heCols; ## This is a logical vector.
     intermediateResults = addNamedList(intermediateResults, HEFilter, ir.HEFilter);
     forceCalculation = TRUE;
   } else {
     xdata.he = HEFilter$xdata.he;
     geneAnnotation.he = HEFilter$geneAnnotation.he;
   }

   if (stopAt=="20-Filtering")
   {
     #if (addToAnalysisManager) am.insertAnalysis(amAnalysis);
     return(intermediateResults);
   }

   #----------------------------------------------------------------------------------------------
   # Data normalization and transformation
   #----------------------------------------------------------------------------------------------
   VST = getElement(intermediateResults, ir.VST);
   if (length(VST)==0 || forceCalculation)
   {
     VST = list();
     if (verbose > 0) printFlush(spaste(spaces, "Normalization and variance stabilization..."))

     if (recognizedPipelines[normPipe]=="Counts")
     {
        expr.vst = varianceStabilizedData.PL(counts = xdata.he, pheno = pheno,
                     design = as.formula(vst.design), blind = vst.blind,
                     calibrateLowestValue = vst.calibrateLowestValue,
                     useApproximateVST = useApproximateVST,
                     approximateVST.minShift = approximateVST.minShift,
                     approximateVST.maxShift = approximateVST.maxShift,
                     approximateVST.xWeightScale = approximateVST.xWeightScale);
        VST$logShift = attr(expr.vst, "logShift");
     } else {
       if (useVarianceStabilization)
       {
         if (useApproximateVST)
         {
           expr.vst = approximateVST(xdata.he, min = approximateVST.minShift,
                    max = approximateVST.maxShift,
                    xIsLogTransformed = !is.null(approximateVST.existingLogBase),
                    xWeightScale = approximateVST.xWeightScale,
                    xLogBase = approximateVST.existingLogBase,
                    xLogShift = approximateVST.existingLogShift);
           VST$logShift = attr(expr.vst, "logShift");
         } else {
           vstFit = vsn2(t(xdata.he), calib = vst.calib);
           expr.vst = predict(vstFit, newdata = xdata);
           if (vst.keepModel) intermediateResults$VST$vsnModel = vstFit;
         }
       } else
         expr.vst = xdata.he;
     }
     VST$expr.vst = expr.vst;
     intermediateResults = addNamedList(intermediateResults, VST, ir.VST);
     forceCalculation = TRUE;
   } else {
     expr.vst = VST$expr.vst;
   }

   if (stopAt=="30-VST") 
   {
     #if (addToAnalysisManager) am.insertAnalysis(amAnalysis);
     return(intermediateResults);
   }

   #----------------------------------------------------------------------------------------------
   # Adjust for covariates, if requested
   #----------------------------------------------------------------------------------------------

   
   if (adj.nProbesForFactorAnalysis > 0)
   {
     if (adj.nProbesForFactorAnalysis > ncol(expr.vst)) adj.nProbesForFactorAnalysis = ncol(expr.vst);
     selectionStat = colMedians(expr.vst, na.rm = TRUE);
     adj.keepProbes = order(-selectionStat)[1:adj.nProbesForFactorAnalysis];
   } else adj.keepProbes = 1:ncol(expr.vst);
     
   expr.forFactorAnalysis = imputeAndScaleData(expr.vst[, adj.keepProbes], scale = scaleForSVD);
   
   if (adj.calculateNPCs < adj.removeNPCs)
     stop(" 'adj.calculateNPCs' must be at least 'adj.removeNPCs'.");
   if (adj.calculateNPCs > 0)
   {
     PCs.filteredInputData = calculatePCs(expr.forFactorAnalysis, nPCs = adj.calculateNPCs, scale = FALSE,
                                          PCColnamePrefix = PCColnamePrefix);
      amAnalysis = conditionalWriteAndAddToAnalysis(
          doSave = savePCs,
          object = dataForTable(PCs.filteredInputData, transpose = FALSE, IDcolName = "Sample"),
          dir = PCDir,
          file = PCFile, 
          fileReadArgs = list(row.names = "Sample"),
          amComponent = amComponent,
          amAnalysis = amAnalysis, amName = amNameForPCs, verbose = verbose-1, indent = indent + 1)
   }
   else
     PCs.filteredInputData = NULL;

   if (adj.removeNPCs > 0) {
      PCCovariates = PCs.filteredInputData[, 1:adj.removeNPCs, drop = FALSE] 
   } else PCCovariates = NULL;

   if (adj.RUV.calculateNFactors > 0)
   {
     if (is.null(adj.RUV.groupVariables)) 
       stop("Group variables for RUV ('adj.RUV.groupVariables') must be given.");
     if (any(!adj.RUV.groupVariables %in% colnames(pheno)))
       stop("All entries in 'adj.RUV.groupVariables' must correspond to colnames in 'pheno'.");
     
     group = apply(pheno[, match(adj.RUV.groupVariables, colnames(pheno)), drop = FALSE], 1, 
                   base::paste, collapse = ".");
     RUVdata = RUVs(t(expr.forFactorAnalysis), cIdx = 1:ncol(expr.forFactorAnalysis), 
                    k = adj.RUV.calculateNFactors, 
                    scIdx = makeGroups(group), round = FALSE, isLog = TRUE);
     RUVFactors = RUVdata$W;
     colnames(RUVFactors) = spaste(RUVColnamePrefix, 1:adj.RUV.calculateNFactors);
     amAnalysis = conditionalWriteAndAddToAnalysis(
          doSave = savePCs,
          object = dataForTable(RUVFactors, transpose = FALSE, IDcolName = "Sample"),
          dir = PCDir,
          file = RUVFactorFile,
          fileReadArgs = list(row.names = "Sample"),
          amComponent = amComponent,
          amAnalysis = amAnalysis, amName = amNameForRUVFactors, verbose = verbose-1, indent = indent + 1)
   } else
     RUVFactors = NULL;
     
   RUVCovariates = if (adj.RUV.useNFactors > 0) RUVFactors[, 1:adj.RUV.useNFactors] else NULL;

   if (adj.calculateSVAFactors)
   {
     model = model.matrix(as.formula(adj.svaModel), data = as.data.frame(pheno));
     if (!is.null(adj.svaReducedModel)) 
        model0 = model.matrix(as.formula(adj.svaReducedModel), data = as.data.frame(pheno)) else model0 = NULL;
     svaData = try(sva(t(expr.forFactorAnalysis), mod = model, mod0 = model0));
     if (inherits(svaData, "try-error")) {
        printFlush("call to sva returned an error. Dropping to browser.");
        browser()
     }
     printFlush("");
     svaFactors = svaData$sv;
     colnames(svaFactors) = spaste(svaColnamePrefix, 1:ncol(svaFactors));
     rownames(svaFactors) = rownames(expr.vst);
     amAnalysis = conditionalWriteAndAddToAnalysis(
          doSave = savePCs,
          object = dataForTable(svaFactors, transpose = FALSE, IDcolName = "Sample"),
          dir = PCDir,
          file = svaFactorFile,
          fileReadArgs = list(row.names = "Sample"),
          amComponent = amComponent,
          amAnalysis = amAnalysis, amName = amNameForSVAFactors, verbose = verbose-1, indent = indent + 1)
   } else svaFactors = NULL;

   svaCovariates = if (adj.sva.useNFactors > 0) svaFactors[, 1:adj.sva.useNFactors, drop = FALSE] else NULL
   # If we're not ading SVA factors, need to remove them from bw.groupBy as well.
   svaNames = colnames(svaFactors);
   removeSVANamesFromGroupBy = setdiff(svaNames, svaCovariates);
   bw.groupBy1 = setdiff(bw.groupBy, removeSVANamesFromGroupBy);

   VST.adj = getElement(intermediateResults, ir.VST.adj);
   if (length(VST.adj)==0 || forceCalculation)
   {
     VST.adj = list();
     if (verbose > 0) printFlush(spaste(spaces, "Adjusting data for unwanted covariates"))
     expr.adj = adjustCovariates.EBLM(
       expr.vst,
       pheno,
       PCs = cbind(PCCovariates, RUVCovariates, svaCovariates),
       removedCovariates = adj.removedCovariates,
       retainedCovariates = adj.retainedCovariates,
       fitToSamples = adj.fitToSamples,
       useBicovWeights = adj.useBicovWeights,
       bw.groupBy = bw.groupBy1,
       bw.groupsForMinWeightRestriction = bw.groupsForMinWeightRestriction,
       bw.otherArgs = bw.otherArgs,
       method = adj.useMethod,
       verbose = verbose - 1, indent = indent + 1);
     gc();
     VST.adj$expr.adj = expr.adj;
     VST.adj$PCs.filteredInputData = PCs.filteredInputData;
     VST.adj$RUVFactors = RUVFactors;
     VST.adj$svaFactors = svaFactors;
     VST.adj$keepProbes = adj.keepProbes;
     intermediateResults = addNamedList(intermediateResults, VST.adj, ir.VST.adj);
     forceCalculation = TRUE;
   } else {
     expr.adj = intermediateResults$VST.adj$expr.adj;
   }

   #if (addToAnalysisManager) am.insertAnalysis(amAnalysis);
   intermediateResults$forceCalculation = forceCalculation;
   list(res =intermediateResults, amAnalysis = amAnalysis); 
}


#======================================================================================================
#
# Convenience function
#
#======================================================================================================

lastHE.VST.OR = function(prepr, comp, name)
{
  n = length(prepr$HE.VST.OR);
  out = prepr$HE.VST.OR[[n]][[comp]][[name]];
  if (is.null(out)) warning("Component ", comp, "$", name, " is empty.");
  out;
}

#======================================================================================================
#
# Weights for distance calculations
#
#======================================================================================================

newWeightFromMediansArgs = function(power = 1)
{
  list(power = power);
}

weightsFromMedians = function(data, power)
{
  meds = replaceMissing(colMedians(data, na.rm = TRUE), 0);
  meds - min(meds)
}
  
   
#======================================================================================================
#
# Main preprocessing function.
#
#======================================================================================================

preprocessGeneral = function(
   # Main input: xdata
   xdata,
   organism,

   # choice of normalization pipeline
   # Choices: Counts, General.
   # Default: Counts for data that are integer.

   normalizationPipeline = NULL,

   # Flow control, for some measure of interactivity
   intermediateResults = list(),
   stopAt = c("End", "10-Entrez", "20-Filtering", "30-VST", "40-Adjustment", "50-Initial sample clustering",
              "60-Outlier sample removal", "80-Individual outlier removal", "90-Adjustment of IOR",
              "100-Collapse"),

   # analysis manager options
   addToAnalysisManager = TRUE,
   analysisName,
   analysisDescription,
   tissue,
   metaData = list(),
   dataSource,
   analysisName.short = analysisName,
   analysisName.pretty = analysisName,
   amComponent = "preprocessing",

   # Optional analysis ID. This will serve as a default for base of file names.
   analysisID = "",
   fileNamePrefix = if (analysisID=="") "" else spaste(analysisID, "-"),

   # Are xdata in genes-in-rows format (need to be transposed)?
   transpose = (ncol(xdata) < nrow(xdata)),
   # geneIDs
   geneIDs = if (transpose) rownames(xdata) else colnames(xdata),
   # Are identifiers symbols?
   idsAreSymbols = FALSE,

   # Phenotypes and covariates together in one data frame
   # A "passenger" sample annotation that will be restricted to non-outlier samples
   sampleAnnotation = NULL,
   # Phenotypes for plotting sample trees. 
   phenoColors,

   # "Passenger" gene annotation
   geneAnnotation = NULL,

   # User-specified filtering, before any analysis starts but after conversion to Entrez.
   geneFilterFnc = NULL,
   geneFilterArgs = list(),
   # preprocessing options
   # filtering: scale samples, use CPM and filter variables with less than 
   heColumnFnc = match.fun("heColumns"),
   minValue = 1, 
   minProportion = 0.25,
   he.useCPM = NULL,  # TRUE if using Counts pipeline, FALSE otherwise
   heColumnOtherArgs = list(scaleSamples = TRUE),

   # Variance stabilization
   useVarianceStabilization = TRUE,   # only useful for general data - VST is mandatory in Counts pipeline

   # Arguments for count data
   vst.design,
   vst.blind = FALSE,
   vst.calibrateLowestValue = TRUE,

   # Arguments for general data
   vst.calib = "affine",
   vst.keepModel = TRUE,

   #Arguments for approximate VST. Can be used within Counts as well as general pipeline.
   useApproximateVST = FALSE,
   approximateVST.minShift = 1,
   approximateVST.maxShift = 500,
   # Caution: if input data are not log transformed, use approximateVST.existingLogBase = NULL
   approximateVST.existingLogBase = 2,
   approximateVST.existingLogShift = 1,
   approximateVST.xWeightScale = 0,

   # Before outlier removal: should VST data be adjusted for covariates?
   # The covariates below must be names of the columns in sampleAnnotation above so they can be
   # used for adjustment of IOMR data.
   # If removedCovariates is of length 0, adjustment is skipped.
   adj.removedCovariates = NULL,
   adj.retainedCovariates = NULL,
   adj.fitToSamples = NULL,
   adj.useMethod = c("EB", "OLS"),
   adj.useBicovWeights = TRUE,

   # Important: if adj.adjustForFactorsBeforeOR is FALSE, only the output expression data will be adjusted
   # for factors.
   adj.adjustForFactorsBeforeOR = FALSE,
   adj.removeNPCs = 0,
   adj.calculateNPCs = 0,
   PCColnamePrefix = "PrincipalComponent.",
   scaleForSVD = TRUE,
   addPCsToSampleAnnotation = TRUE,

   adj.RUV.useNFactors = 0,
   adj.RUV.calculateNFactors = adj.RUV.useNFactors,
   adj.RUV.groupVariables = NULL,
   RUVColnamePrefix = "RUVFactor.",

   adj.calculateSVAFactors = FALSE,
   adj.sva.useNFactors = 0,
   adj.svaModel,
   adj.svaReducedModel = NULL,
   svaColnamePrefix = "SVAFactor.",

   bw.groupBy = NULL,
   bw.groupsForMinWeightRestriction = NULL,
   bw.otherArgs = list(maxPOutliers = 0.10,
                       outlierReferenceWeight = 0.1,
                       minWeightInGroups = 0.9,
                       maxPropUnderMinWeight = 0.4,
                       defaultWeight = 1),

   # Outlier removal options
   outlierRemovalZ = 8,
   or.robustScale = TRUE,
   or.nTopProbes = 8000,
   or.scaleData = TRUE,
   or.distWeightArgs = list(),

   # Individual outlier removal and replacement

   ior.replace = TRUE,
   ior.remove = FALSE,

   ior.replaceThreshold = 0.5, 
   ior.weightedReplace = TRUE,
   ior.removeThreshold = 0.1,
   restrictIORExprToHE = TRUE,
   restrictIORXDataToHE = FALSE,

   iorFromHEGenes = FALSE,  # If TRUE, entire IOR calculation will start from HE genes
   # Note that final IOR xdata or expr could conceivably contain different genes since HE filter on IOR data may result in
   # a subset of original HE genes  

   # Collapsing of final data. The idea is to select the representative features based on number of missing
   # values first, then on a user-specified criterion, but give priority to features that appear in the final
   # filtered data.

   collapseSteps = allCollapseSteps(),
   cl.probe2gene = NULL,
   cl.method = "MaxMean",

   # Plot file arguments
   plotDir = "Plots",
   plotFnc = "pdf",
   plotExt = ".pdf",
   stp.filePrefix = fileNamePrefix,

   # Sample tree plot (stp) arguments
   stp.width = 12,
   stp.height = 7,
   stp.mainBase,
   stp.mainSep = "\n",
   stp.maxCharPerLine = 30,
   stp.marAll = c(1, 11, 3, 1),
   stp.includeSampleLabels = FALSE,
   stp.mainPart1 = "var.-stabilized data",
   stp.mainPart2 = "adjusted data",
   stp.mainPart3 = "adjusted data, outlier samples removed",

   # PCA plot arguments
   pca.mainBase = stp.mainBase,
   pca.mainSep = "\n",
   pca.main2 = spaste("Corrected for ", pca.correctForColumn),
   pca.correctForColumn = NULL,  # This needs to be a column name in sample annotation
   pca.plotLegend = TRUE,
   pca.separateLegend = TRUE,
   pca.legendWidth = 1.6,
   pca.height = 4,
   pca.colorColumn = NULL,
   pca.shapeColumn = NULL,
   pca.colorLevels = list(),
   pca.shapeLevels = list(),
   pca.colorPrefix = rep("", length(pca.colorColumn)),
   pca.shapePrefix = rep("", length(pca.shapeColumn)),
   pca.filePrefix = fileNamePrefix,
   pca.pt.cex = 2.5, 
   pca.pt.bg = "grey70",   # Default value when color factor is not given
   pca.cex.main = 1,
   pca.includePointLabels = FALSE,
   pca.cex.pointLabels = 0.7,
   pca.includePVE = TRUE,
   pca.thinned = FALSE,

   prettifyList = NULL,

   # Saving of intermediate data
   saveDir.base,
   saveDir.expr = "Expression",
   saveDir.sampleAnnot = "SampleAnnotation",
   saveXData.entrez = TRUE,  ## Note: input sample annotation and (possibly filtered) input gene annotation are always saved. 
                             ## The sample annotation is small and gene annotation is typically needed in downstream analyses.
   IDcolName.out = if (is.null(cl.probe2gene)) "Entrez" else "ProbeID",
   IDcolName.out.collapsed = "Entrez",
 
   exprDir.entrez = file.path(saveDir.base, saveDir.expr, "015-XDataWithValidEntrez"),
   xdataFile.entrez = spaste(fileNamePrefix, "xdataWithValidEntrez.csv.gz"),
   dataFile.entrez = spaste(fileNamePrefix, "xdataWithValidEntrez-sampleAnnot.RData"),

   sampleDir.input = file.path(saveDir.base, saveDir.sampleAnnot, "020-AsSupplied"),
   PCFile.HEData = spaste(fileNamePrefix, "principalComponents-HEData.csv"),
   RUVFactorFile.HEData = spaste(fileNamePrefix, "RUVFactors-HEData.csv"),
   SVAFactorFile.HEData = spaste(fileNamePrefix, "SVAFactors-HEData.csv"),

   saveOSR = TRUE,
   exprDir.OSR = file.path(saveDir.base, saveDir.expr, "030-OutlierSamplesRemoved"),
   xdataFile.OSR= spaste(fileNamePrefix, "xdata.OSR.csv.gz"),
   exprFile.OSR = spaste(fileNamePrefix, "expr.OSR.csv.gz"),
   dataFile.OSR = spaste(fileNamePrefix, "expr-xdata-sampleAnnotation.OSR.RData"),
   sampleDir.OSR = file.path(saveDir.base, saveDir.sampleAnnot, "030-OutlierSamplesRemoved"),
   sampleFile.OSR = spaste(fileNamePrefix, "sampleAnnotation.OSR.csv"),
   PCFile.HEData.OSR = spaste(fileNamePrefix, "principalComponents-HEData.OSR.csv"),
   RUVFactorFile.HEData.OSR = spaste(fileNamePrefix, "RUVFactors-HEData.OSR.csv"),
   SVAFactorFile.HEData.OSR = spaste(fileNamePrefix, "SVAFactors-HEData.OSR.csv"),

   saveIOR = TRUE,
   exprDir.IORep = file.path(saveDir.base, saveDir.expr, "050-IndividualOutliersReplaced"),
   xdataFile.IORep = spaste(fileNamePrefix, "xdata.IORep.csv.gz"),
   exprFile.IORep = spaste(fileNamePrefix, "expr.IORep.csv.gz"),
   dataFile.IORep = spaste(fileNamePrefix, "expr-xdata-sampleAnnot.IORep.RData"),
   PCFile.IORep = spaste(fileNamePrefix, "principalComponents-IORepData.csv"), 

   exprDir.IORem = file.path(saveDir.base, saveDir.expr, "051-IndividualOutliersRemoved"),
   xdataFile.IORem = spaste(fileNamePrefix, "xdata.IORem.csv.gz"),
   exprFile.IORem = spaste(fileNamePrefix, "expr.IORem.csv.gz"),
   dataFile.IORem = spaste(fileNamePrefix, "expr-xdata-sampleAnnot.IORem.RData"),
   PCFile.IORem = spaste(fileNamePrefix, "principalComponents-IORemData.csv"), 

   # Saving of weights
   saveIORWeights = TRUE, ## This is effective even if neither of ior.replace and ior.remove are TRUE
   IORWeightDir = exprDir.OSR,
   xdataWeightFile.IORep = spaste(fileNamePrefix, "xdataWeightsForIOR.csv.gz"),
   xdataWeightFactorFile.IORep = spaste(fileNamePrefix, "xdataWeightFactorsForIOR.csv.gz"),
   weightFile.IORep = spaste(fileNamePrefix, "xdataWeightsAndFactors.RData"),

   # Saving of adjusted data
   saveIOR.adj = TRUE,
   exprDir.IORep.adj = file.path(saveDir.base, saveDir.expr, "060-IOReplaced-adjusted"),
   exprFile.IORep.adj = spaste(fileNamePrefix, "expr.IORep.adj.csv.gz"),
   dataFile.IORep.adj = spaste(fileNamePrefix, "expr-xdata-sampleAnnot.IORep.adj.RData"),

   exprDir.IORem.adj = file.path(saveDir.base, saveDir.expr, "061-IORemoved-adjusted"),
   exprFile.IORem.adj = spaste(fileNamePrefix, "expr.IORem.adj.csv.gz"),
   dataFile.IORem.adj = spaste(fileNamePrefix, "expr-xdata-sampleAnnot.IORem.adj.RData"),

   # Saving of collapsed data
   saveCollapsed = TRUE,
   exprDir.OSR.collapsed = file.path(saveDir.base, saveDir.expr, "040-OSR-Collapsed"),
   xdataFile.OSR.collapsed = spaste(fileNamePrefix, "xdata.OSR.collapsed.csv.gz"),
   exprFile.OSR.collapsed = spaste(fileNamePrefix, "expr.OSR.collapsed.csv.gz"),
   geneAnnotFile.xdata.OSR.collapsed = spaste(fileNamePrefix, "geneAnnot.xdata.OSR.collapsed.csv.gz"),
   geneAnnotFile.expr.OSR.collapsed = spaste(fileNamePrefix, "geneAnnot.expr.OSR.collapsed.csv.gz"),
   dataFile.OSR.collapsed = spaste(fileNamePrefix, "expr-xdata-sampleAnnotation.collapsed.OSR.RData"),

   exprDir.IORep.collapsed = file.path(saveDir.base, saveDir.expr, "070-IOReplaced-collapsed"),
   xdataFile.IORep.collapsed = spaste(fileNamePrefix, "xdata.IORep.collapsed.csv.gz"),
   exprFile.IORep.collapsed = spaste(fileNamePrefix, "expr.IORep.collapsed.csv.gz"),
   ### FIXME: saving into these files is not implemented yet.
   geneAnnotFile.xdata.IORep.collapsed = spaste(fileNamePrefix, "geneAnnot.xdata.IORep.collapsed.csv.gz"),
   geneAnnotFile.expr.IORep.collapsed = spaste(fileNamePrefix, "geneAnnot.expr.IORep.collapsed.csv.gz"),
   dataFile.IORep.collapsed = spaste(fileNamePrefix, "expr-xdata-sampleAnnot.IORep.collapsed.RData"),

   exprDir.IORem.collapsed = file.path(saveDir.base, saveDir.expr, "071-IORemoved-collapsed"),
   xdataFile.IORem.collapsed = spaste(fileNamePrefix, "xdata.IORem.collapsed.csv.gz"),
   exprFile.IORem.collapsed = spaste(fileNamePrefix, "expr.IORem.collapsed.csv.gz"),
   ### FIXME: saving into these files is not implemented yet.
   geneAnnotFile.xdata.IORem.collapsed = spaste(fileNamePrefix, "geneAnnot.xdata.IORem.collapsed.csv.gz"),
   geneAnnotFile.expr.IORem.collapsed = spaste(fileNamePrefix, "geneAnnot.expr.IORem.collapsed.csv.gz"),
   dataFile.IORem.collapsed = spaste(fileNamePrefix, "expr-xdata-sampleAnnot.collapsed.IORem.RData"),

   exprDir.IORep.adj.collapsed = file.path(saveDir.base, saveDir.expr, "080-IOReplaced-adjusted-collapsed"),
   exprFile.IORep.adj.collapsed = spaste(fileNamePrefix, "expr.IORep.collapsed.adj.csv.gz"),
   dataFile.IORep.adj.collapsed = spaste(fileNamePrefix, "expr-xdata-sampleAnnot.IORep.adj.collapsed.RData"),

   exprDir.IORem.adj.collapsed = file.path(saveDir.base, saveDir.expr, "080-IOReplaced-adjusted-collapsed"),
   exprFile.IORem.adj.collapsed = spaste(fileNamePrefix, "expr.IORem.adj.collapsed.csv.gz"),
   dataFile.IORem.adj.collapsed = spaste(fileNamePrefix, "expr-xdata-sampleAnnot.IORem.adj.collapsed.RData"),

   saveCollapsedIORWeights = TRUE, ## This is effective even if neither of ior.replace and ior.remove are TRUE
   IORWeightDir.collapsed = exprDir.OSR.collapsed,
   recalculateCollapsedWeights = FALSE,  ## Set to TRUE if collapseRows uses a summaries rather than representatives
   xdataWeightFile.IORep.collapsed = spaste(fileNamePrefix, "xdataWeightsForIOR.collapsed.csv.gz"),
   xdataWeightFactorFile.IORep.collapsed = spaste(fileNamePrefix, "xdataWeightFactorsForIOR.collapsed.csv.gz"),
   weightFile.IORep.collapsed = spaste(fileNamePrefix, "xdataWeightsAndFactors.collapsed.RData"),


   # Behavior
   forceCalculation = FALSE,
   verbose = 1,
   indent = 0

)
{
   # A few convenience functions, really macros
   plotSampleTree = function(tree, colors, file, sub, dendroLabels = FALSE)
   {
     suppressWarnings(dir.create(plotDir, recursive = TRUE, showWarnings = FALSE));
     plotFnc(file= spaste(plotDir, "/", stp.filePrefix, "-", file, plotExt), wi=stp.width, he = stp.height);
     on.exit(dev.off());
     plotDendroAndColors(tree,
                 colors$colors,
                 formatLabels(colors$legend, maxCharPerLine = stp.maxCharPerLine),
                 addGuide = TRUE, guideAll = TRUE, dendroLabels = dendroLabels,
                 guideHang = 0.10,
                 main = spaste(stp.mainBase, stp.mainSep, sub, "               "),
                 marAll = stp.marAll);
     #dev.off();
     #on.exit(NULL);
   }

   plotSampleTreeFromData = function(data, scale, colors, file, sub, ...)
   {
     tree = hclust(dist(if (scale) scale(data) else data), method = "a");
     plotSampleTree(tree, colors, file, sub, ...);
   }

   plotPCA.local = function(data, sampleAnnot, scale, file, sub, outlierIndicator = NULL,
                     pointLabels = NULL)
   {
     if (any(is.na(data))) 
     {
       gg = goodGenes(data);
       if (any(!gg)) data = data[, gg];
       if (any(is.na(data))) data = t(impute.knn(t(data))$data);
     }
     if (length(pca.correctForColumn) > 0)
     {
       correctCol = match(pca.correctForColumn, colnames(sampleAnnot));
       if (is.na(correctCol))
          stop("Invalid 'pca.correctForColumn': `", pca.correctForColumn, "' not found\n", 
               "among the column names of 'sampleAnnotation: \n",
               formatLabels(paste(colnames(sampleAnnot), collapse = ", "),
                  maxCharPerLine = options("width")$width-3, split = ", "));
       colData = sampleAnnot[, correctCol];
       data.corr = try(residuals(lm(data~colData, na.action = "na.exclude")));
       if (inherits(data.corr, "try-error")) 
       {
         printFlush("Adjusting in plotPCA.local returned an error. Dropping into browser.");
         browser();
       }
     } else { correctCol = NULL; data.corr = NULL; }
 
    main1 = spaste(pca.mainBase, pca.mainSep, sub, "               ");
    colorCol = pca.colorColumn;
    shapeCol = pca.shapeColumn;
    lc = length(pca.colorColumn);
    ls = length(pca.shapeColumn);
    if (lc*ls > 0)
    {
      if (lc!=ls)
        stop("When both 'pca.colorColumn' and 'pca.shapeColumn' are given, their lengths must be the same.");
    }
    if (lc > 0 && lc!=length(pca.colorPrefix))
      stop("When 'pca.colorColumn' is given, 'pca.colorPrefix' must have the same length.");
    
    if (ls > 0 &&  ls!=length(pca.shapePrefix))
      stop("When 'pca.shapeColumn' si given, 'pca.shapePrefix' must have the same length.");
    plotLegend = pca.plotLegend && lc > 0
    nPlots = max(1, ls, lc);
   
    for (cc in 1:nPlots)
    {
       if (length(pca.colorColumn)>0 && pca.colorColumn[cc]!="") 
       {
         if (length(pca.colorLevels)==0 || length(pca.colorLevels[[cc]])==0) 
         {
            pca.colorLevels1 = sort(unique(sampleAnnot[, pca.colorColumn[cc]]));
         } else pca.colorLevels1 = pca.colorLevels[[cc]];
         colorFactor = factor(sampleAnnot[, pca.colorColumn[cc]], levels = pca.colorLevels1);
       } else colorFactor = NULL;

       if (length(pca.shapeColumn) > 0 && pca.shapeColumn[cc]!="") 
       {
         if (length(pca.shapeLevels)==0 || length(pca.shapeLevels[[cc]])==0) 
         {
            pca.shapeLevels1 = sort(unique(sampleAnnot[, pca.shapeColumn[cc]]));
         } else pca.shapeLevels1 = pca.shapeLevels[[cc]];
         shapeFactor = factor(sampleAnnot[, pca.shapeColumn[cc]], levels = pca.shapeLevels1);
       } else shapeFactor = NULL;

       plotPCA.2panels(
          data1 = data,
          data2 = data.corr,
          colorFactor = colorFactor,
          shapeFactor = shapeFactor,
          borderFactor = outlierIndicator,
          pointLabels = if (pca.includePointLabels) pointLabels else NULL,
          colorPrefix = pca.colorPrefix[cc],
          shapePrefix = pca.shapePrefix[cc],
          scale = scale, device = plotFnc, dir = plotDir,
          file = spaste(pca.filePrefix, "-", file, if (nPlots > 1) spaste("-", cc) else "", plotExt),
          main1 = main1,
          main2 = pca.main2,
          includePVE = pca.includePVE,
          plotLegend = plotLegend,
          legendWidth = pca.legendWidth/pca.height,
          setLayout = TRUE,
          height = pca.height,
          width = pca.height*(1 + !is.null(data.corr)) + pca.legendWidth * pca.separateLegend * plotLegend,
          pt.cex = pca.pt.cex,
          pt.bg = pca.pt.bg,
          thinned = pca.thinned,
          cex.pointLabels = pca.cex.pointLabels,
          fig.cex = 1,
          mar = c(3.2, 3.2, pca.cex.main * (nLines(main1) + 0.2), 1),
          mgp = c(2, 0.7, 0),
          legendOnly = FALSE,
          prettifyList = prettifyList,
          cex.main = pca.cex.main);
     }
   }
   #----------------------------------------------------------------------------------------------
   # Start of actual code
   #----------------------------------------------------------------------------------------------
  
   stopAt = match.arg(stopAt);
   spaces = indentSpaces(indent);
   plotFnc = match.fun(plotFnc);

   if (adj.RUV.useNFactors > 0 && adj.sva.useNFactors > 0)
     stop("RUV and sva factors cannot be used at the same time to adjust data (but\n",
          "  they can be computed at the same time for future use.)");

   if (transpose) xdata = t(xdata);

   if (is.null(normalizationPipeline))
   {
     if (all(round(xdata)==xdata, na.rm = TRUE))
     {
       normalizationPipeline = "Counts"
     } else 
       normalizationPipeline = "General";
   }
   recognizedPipelines = c("Counts", "General");
   normPipe = charmatch(tolower(normalizationPipeline), tolower(recognizedPipelines));

   if (is.na(normPipe)) 
       stop("If given, 'normalizationPipeline' must be one of\n   ", 
            paste(recognizedPipelines, collapse = ", "));

   if (verbose > 1)
   {
      printFlush(spaste(spaces, "Basic configuration:"));
      printFlush(spaste(spaces, "  analysisName: ", if (addToAnalysisManager) analysisName else "(not applicable)"));
      printFlush(spaste(spaces, "  analysisID: ", analysisID));
      printFlush(spaste(spaces, "   fileNamePrefix: ", fileNamePrefix));
      printFlush(spaste(spaces, "   normalization pipeline: ", normalizationPipeline));
   }
   
   if (is.null(he.useCPM)) he.useCPM = normPipe=="Counts";

   if (addToAnalysisManager)
   {
      am.analysisID = am.analysisID(name = analysisName, throw = FALSE);
      if (is.na(am.analysisID))
      {
        amAnalysis = am.newAnalysis(name = analysisName,
           description = analysisDescription,
           analysisType = amComponent,
           organism = organism,
           tissue = tissue,
           shortName = analysisName.short,
           prettyName = analysisName.pretty,
           metaData = metaData,
           dataSource = dataSource,
           allowNonstandardType = TRUE);
      } else
        amAnalysis = am.getAnalysis(am.analysisID);
   } else amAnalysis = NULL;

   if (saveCollapsedIORWeights) collapseSteps = unique(c(collapseSteps, "OSR"));

   if (length(collapseSteps)==0) cl.probe2gene = NULL;

   if (length(cl.probe2gene) > 0)
   {
     if (any(! (collapseSteps %in% allCollapseSteps())))
        stop("Entries of 'collapseSteps' must be one of the following options:\n",
             formatStrings(paste(allCollapseSteps(), collapse = ", "), 80));

   }

   nSamples = nrow(xdata);
   if (length(adj.fitToSamples) > 0)
   {
     if (!is.list(adj.fitToSamples)) adj.fitToSamples = list(adj.fitToSamples);
     adj.fitToSamples = lapply(adj.fitToSamples, function(fit1)
     {
       if (is.numeric(fit1)) 
       {
         tmp = fit1;
         fit1 = rep(FALSE, nSamples);
         fit1[tmp] = TRUE;
       }

       if (length(fit1)!=nSamples)
         stop("When given, 'adj.fitToSamples' or each component should be logical\n",
              "with length equal number of samples in 'xdata'.");
       fit1;
     })
   }

   nGenes = ncol(xdata);
   if (is.null(or.nTopProbes) || or.nTopProbes < 1) or.nTopProbes = nGenes;

   if (!"Files" %in% names(intermediateResults)) intermediateResults$Files = list();

   if (length(intermediateResults$entrez)==0 || forceCalculation)
   {
     intermediateResults$entrez = list();
     if (idsAreSymbols) 
     {
       if (verbose > 0) printFlush(spaste(spaces, "Converting column IDs to Entrez..."));
       xdata.entrez = convertColumnNames2entrez(xdata, organism = organism,
                      transpose = FALSE, sourceSymbols = geneIDs);
       geneAnnotation.entrez = if (is.null(geneAnnotation)) NULL else
            geneAnnotation[attr(xdata.entrez, "indexInOriginalData"), ];

     } else {
       xdata.entrez = xdata;
       colnames(xdata.entrez) = geneIDs;
       geneAnnotation.entrez = geneAnnotation;
     }

     # If additional filtering is specified, filter at this point to exclude the genes from the entire
     # analysis.
     if (!is.null(geneFilterFnc))
     {
       if (verbose > 0) printFlush(spaste(spaces, "User-defined filter..."));
       fnc = match.fun(geneFilterFnc);
       userCols = do.call(fnc, c(list(x = xdata.entrez), geneFilterArgs));
       if (!is.logical(userCols) | length(userCols)!=ncol(xdata.entrez))
         stop("If given, 'geneFilterFnc' must return a logical vector with one element per input gene.");
       xdata.entrez = xdata.entrez[, userCols];
       attr(xdata.entrez, "indexInOriginalData") = attr(xdata.entrez, "indexInOriginalData")[userCols];
       geneAnnotation.entrez = geneAnnotation.entrez[userCols, ]
     }
     intermediateResults$entrez$xdata.entrez = xdata.entrez;
     intermediateResults$entrez$geneAnnotation.entrez = geneAnnotation;
     dir.create(exprDir.entrez, recursive = TRUE, showWarnings = FALSE);
     if (saveXData.entrez)
     {
       if (verbose > 1) printFlush(spaste(spaces, "  ...saving filtered data..."));
       write.csv(dataForTable(xdata.entrez, IDcolName = "Name"),
                 file = gzfile(file.path(exprDir.entrez, xdataFile.entrez)),
                 quote = TRUE, row.names = FALSE);
     }
     save(list = c(if (saveXData.entrez) "xdata.entrez" else character(0), "sampleAnnotation", "geneAnnotation.entrez"),
          file = file.path(exprDir.entrez, dataFile.entrez));
     intermediateResults$Files$exprDir.entrez = exprDir.entrez
     intermediateResults$Files$xdataFile.entrez = xdataFile.entrez;
     intermediateResults$Files$dataFile.entrez = dataFile.entrez;
     if (verbose > 1) printFlush(spaste(spaces, "  ...saving input sample annotation and gene annotation..."));
     if (addToAnalysisManager) amAnalysis = am.addMultipleObjectsToAnalysis(
        names = c(if (saveXData.entrez) "input data mapped to Entrez" else character(0),
                  "feature data for individual analysis", "input sample information"),
        fileName = file.path(exprDir.entrez, dataFile.entrez),
        fileObject = c(if (saveXData.entrez) "xdata.entrez" else character(0),
                     "geneAnnotation.entrez", "sampleAnnotation"),
        anaType = amComponent,
        allowNonstandardAnalysisType = TRUE,
        analysis = amAnalysis, verbose = verbose -1);
     forceCalculation = TRUE;
   } else {
     xdata.entrez = intermediateResults$entrez$xdata.entrez;
     geneAnnotation.entrez = intermediateResults$entrez$geneAnnotation.entrez;
   }

   if (stopAt=="10-Entrez") 
   {
     if (addToAnalysisManager) am.insertAnalysis(amAnalysis);
     return(intermediateResults);
   }

   
   
   #----------------------------------------------------------------------------------------------
   # Iterate HE, VST, adjustment and outlier removal until no further outliers are removed.
   #----------------------------------------------------------------------------------------------

   ORStep = 1;
   HE.VST.OR = intermediateResults$HE.VST.OR;
   if (length(HE.VST.OR)==0) HE.VST.OR = list();
   changed = TRUE;
   nSamples.OSR = nrow(xdata.entrez);
   xdata.OSR = xdata.entrez
   geneAnnotation.OSR = geneAnnotation.entrez;
   phenoColors.OSR = phenoColors;
   sampleAnnotation.OSR = sampleAnnotation;
   bw.groupsForMinWeightRestriction.OSR = bw.groupsForMinWeightRestriction

   while (changed)
   {
      if (verbose > 0) printFlush(spaste(spaces, "HE-VST-OSR step ", ORStep))
      if (length(HE.VST.OR) < ORStep) HE.VST.OR[[ORStep]] = list();
      phenoColors = phenoColors.OSR;
      sampleAnnotation = sampleAnnotation.OSR;
      tmp1 = HE.VST.Adjust(
        xdata = xdata.OSR, pheno = sampleAnnotation.OSR,
        geneAnnotation = geneAnnotation.OSR,
        intermediateResults = HE.VST.OR[[ORStep]],

        normalizationPipeline = normalizationPipeline,
        forceCalculation = forceCalculation,
        stopAt = stopAt,

        amAnalysis = amAnalysis,
        amComponent = amComponent,

        ir.HEFilter = "HEFilter",
        heColumnFnc = heColumnFnc, minValue = minValue, minProportion = minProportion,
        he.useCPM = he.useCPM,
        heColumnOtherArgs = heColumnOtherArgs,

        ir.VST = "VST",
        useVarianceStabilization = useVarianceStabilization,

        vst.design = vst.design, vst.blind = vst.blind,
        vst.calibrateLowestValue = vst.calibrateLowestValue,

        vst.calib = vst.calib, vst.keepModel = vst.keepModel,

        useApproximateVST = useApproximateVST,
        approximateVST.minShift = approximateVST.minShift,
        approximateVST.maxShift = approximateVST.maxShift,
        approximateVST.existingLogBase = approximateVST.existingLogBase,
        approximateVST.existingLogShift = approximateVST.existingLogShift,
        approximateVST.xWeightScale = approximateVST.xWeightScale,

        ir.VST.adj = "VST.adj",
        adj.removedCovariates = adj.removedCovariates, adj.retainedCovariates = adj.retainedCovariates,
        adj.useMethod = adj.useMethod, adj.fitToSamples = adj.fitToSamples, adj.useBicovWeights = adj.useBicovWeights,
        adj.nProbesForFactorAnalysis = or.nTopProbes,
        adj.removeNPCs = if (adj.adjustForFactorsBeforeOR) adj.removeNPCs else 0,
        adj.calculateNPCs = adj.calculateNPCs, 
        PCColnamePrefix = PCColnamePrefix, 
        savePCs = saveOSR || (!is.null(cl.probe2gene) && ("OSR" %in% collapseSteps)),
        amNameForPCs = spaste("principal components of filtered ", if (ORStep==1) "" else "and OR ", "data"),
        PCFile = if (ORStep==1) PCFile.HEData else PCFile.HEData.OSR,
        PCDir = if (ORStep==1) sampleDir.input else sampleDir.OSR,
        
        scaleForSVD = scaleForSVD,
        adj.RUV.useNFactors = if (adj.adjustForFactorsBeforeOR) adj.RUV.useNFactors else 0,
        adj.RUV.calculateNFactors = adj.RUV.calculateNFactors,
        adj.RUV.groupVariables = adj.RUV.groupVariables,
        RUVColnamePrefix = RUVColnamePrefix,
        RUVFactorFile = if (ORStep==1) RUVFactorFile.HEData else RUVFactorFile.HEData.OSR,
        amNameForRUVFactors = spaste("RUV factors of filtered ", if (ORStep==1) "" else "and OR ", "data"),

        adj.calculateSVAFactors = adj.calculateSVAFactors,
        adj.sva.useNFactors = if (adj.adjustForFactorsBeforeOR) adj.sva.useNFactors else 0,
        adj.svaModel = adj.svaModel,
        adj.svaReducedModel = adj.svaReducedModel,
        svaColnamePrefix = svaColnamePrefix,
        svaFactorFile = if (ORStep==1) SVAFactorFile.HEData else SVAFactorFile.HEData.OSR,
        amNameForSVAFactors = spaste("SVA factors of filtered ", if (ORStep==1) "" else "and OR ", "data"),


        bw.groupBy = bw.groupBy,
        bw.groupsForMinWeightRestriction = bw.groupsForMinWeightRestriction.OSR,
        bw.otherArgs = bw.otherArgs,
        verbose = verbose, indent = indent + 1);
      HE.VST.OR[[ORStep]] = tmp1$res;
      amAnalysis = tmp1$amAnalysis;

      forceCalculation = HE.VST.OR[[ORStep]]$forceCalculation;
      expr.vst = HE.VST.OR[[ORStep]]$VST$expr.vst;
      expr.adj = HE.VST.OR[[ORStep]]$VST.adj$expr.adj;
      geneAnnotation.he = HE.VST.OR[[ORStep]]$HEFilter$geneAnnotation.he;
      
      if (stopAt%in% c("20-Filtering", "30-VST", "40-Adjustment"))
      {
        if (addToAnalysisManager) am.insertAnalysis(amAnalysis);
        intermediateResults$HE.VST.OR = HE.VST.OR;
        return(intermediateResults);
      }
   
      if (stopAt=="50-Initial sample clustering")
      {
        if (addToAnalysisManager) am.insertAnalysis(amAnalysis);
        intermediateResults$HE.VST.OR = HE.VST.OR;
        return(intermediateResults);
      }

      #----------------------------------------------------------------------------------------------
      # Iterative outlier removal
      #----------------------------------------------------------------------------------------------

      if (length(HE.VST.OR[[ORStep]]$OSR)==0 || forceCalculation)
      {
         HE.VST.OR[[ORStep]]$OSR = list();
         if (verbose > 0) printFlush(spaste(spaces, "  Iterative outlier removal..."))
         data.OSR = iterativeOutlierRemoval.withDetails(expr.adj, 
            Zcut = outlierRemovalZ, useProbes = HE.VST.OR[[ORStep]]$VST.adj$keepProbes, 
            robustScale = or.robustScale, 
            scaleData = or.scaleData);
         expr.OSR = data.OSR$ORData;

         outliers = data.OSR$removedOutliers;

         if (verbose > 1)
         {
            printFlush(paste(spaces, "  outliers: "));
            print(table(outliers))
         }
         xdata.OSR = xdata.OSR[!outliers, , drop = FALSE];
         sampleAnnotation.OSR = sampleAnnotation.OSR[!outliers, , drop = FALSE];
         phenoColors.OSR$colors = phenoColors.OSR$colors[!outliers, , drop = FALSE];
         bw.groupsForMinWeightRestriction.OSR = bw.groupsForMinWeightRestriction.OSR[!outliers];
         if (length(adj.fitToSamples)>0) adj.fitToSamples = lapply(adj.fitToSamples, `[`, !outliers);
        
         # If there are no more outliers, add the PCs (which were computed before outlier removal) to the
         # sample data (sampleAnnotation) that will be in the output, to facilitate downstream analysis. 

         if (!any(outliers))
         {
           if (adj.calculateNPCs > 0 && addPCsToSampleAnnotation)
             sampleAnnotation.OSR = cbind(sampleAnnotation.OSR, HE.VST.OR[[ORStep]]$VST.adj$PCs.filteredInputData);
           if (adj.RUV.calculateNFactors > 0 && addPCsToSampleAnnotation)
             sampleAnnotation.OSR = cbind(sampleAnnotation.OSR, HE.VST.OR[[ORStep]]$VST.adj$RUVFactors);
           if (adj.calculateSVAFactors && addPCsToSampleAnnotation)
             sampleAnnotation.OSR = cbind(sampleAnnotation.OSR, HE.VST.OR[[ORStep]]$VST.adj$svaFactors);

           if (adj.removeNPCs > 0) {
             PCCovariates = HE.VST.OR[[ORStep]]$VST.adj$PCs.filteredInputData[, 1:adj.removeNPCs, drop = FALSE]
           } else PCCovariates = NULL;

           if (adj.RUV.useNFactors > 0) {
             RUVCovariates = HE.VST.OR[[ORStep]]$VST.adj$RUVFactors[, 1:adj.RUV.useNFactors, drop = FALSE]
           } else RUVCovariates = NULL;

           if (adj.sva.useNFactors > 0) {
             svaCovariates = HE.VST.OR[[ORStep]]$VST.adj$svaFactors[, 1:adj.sva.useNFactors, drop = FALSE]
           } else svaCovariates = NULL;

           if (!adj.adjustForFactorsBeforeOR)
              expr.OSR = adjustCovariates.EBLM(
                 expr.OSR,
                 sampleAnnotation.OSR,
                 PCs = cbind(PCCovariates, RUVCovariates, svaCovariates),
                 removedCovariates = adj.removedCovariates,
                 retainedCovariates = adj.retainedCovariates,
                 fitToSamples = adj.fitToSamples,
                 useBicovWeights = adj.useBicovWeights,
                 bw.groupBy = bw.groupBy,
                 bw.groupsForMinWeightRestriction = bw.groupsForMinWeightRestriction.OSR,
                 bw.otherArgs = bw.otherArgs,
                 method = adj.useMethod,
                 verbose = verbose - 1, indent = indent + 1); 
         }

         forceCalculation = TRUE;
         HE.VST.OR[[ORStep]]$OSR$expr.OSR = expr.OSR;
         HE.VST.OR[[ORStep]]$OSR$xdata.OSR = xdata.OSR;
         HE.VST.OR[[ORStep]]$OSR$outliers = outliers;
         HE.VST.OR[[ORStep]]$OSR$adj.fitToSamples = adj.fitToSamples;
         HE.VST.OR[[ORStep]]$OSR$phenoColors.OSR = phenoColors.OSR;
         HE.VST.OR[[ORStep]]$OSR$sampleAnnotation.OSR = sampleAnnotation.OSR;
         HE.VST.OR[[ORStep]]$OSR$bw.groupsForMinWeightRestriction.OSR = bw.groupsForMinWeightRestriction.OSR;
         HE.VST.OR[[ORStep]]$OSR$OSRDetails = data.OSR$details;
         HE.VST.OR[[ORStep]]$OSR$or.topProbes = data.OSR$topProbes;  # This is a numeric index
      } else {
         expr.OSR = HE.VST.OR[[ORStep]]$OSR$expr.OSR;
         xdata.OSR = HE.VST.OR[[ORStep]]$OSR$xdata.OSR;
         sampleAnnotation.OSR = HE.VST.OR[[ORStep]]$OSR$sampleAnnotation.OSR;
         outliers = HE.VST.OR[[ORStep]]$OSR$outliers;
         phenoColors.OSR = HE.VST.OR[[ORStep]]$OSR$phenoColors.OSR;
         bw.groupsForMinWeightRestriction.OSR = HE.VST.OR[[ORStep]]$OSR$bw.groupsForMinWeightRestriction.OSR;
      }
      # Produce plots
      colors.x = phenoColors;
      colors.x$colors = cbind(phenoColors$colors, Outlier = outliers+1);
      colors.x$legend = c(phenoColors$legend, "Outlier");

      #plotSampleTreeFromData(expr.vst[, HE.VST.OR[[ORStep]]$OSR$or.topProbes], scale = or.scaleData, phenoColors,
      #            file = spaste("sampleTree-Step.", ORStep, "-010-NormalizedData"),
      #            sub = "var.-stabilized data", dendroLabels = if (stp.includeSampleLabels) NULL else FALSE);

      #zz1 = HE.VST.OR[[ORStep]]$OSR$OSRDetails$Step.1$Z;
      #keep = rank(-zz1) < 7;
      plotPCA.local(expr.vst[, HE.VST.OR[[ORStep]]$OSR$or.topProbes], sampleAnnotation, scale = or.scaleData, 
                  outlierIndicator = outliers + 1, 
                  pointLabels = ifelse(outliers, rownames(expr.vst), ""),
                  file = spaste("PCA-Step.", ORStep, "-010-NormalizedData"),
                  sub = stp.mainPart1)


      plotSampleTreeFromData(expr.vst[, HE.VST.OR[[ORStep]]$OSR$or.topProbes], scale = or.scaleData, colors.x,
                  file = spaste("sampleTree-Step.", ORStep, "-011-NormalizedDataWithOutlierIndicator"),
                  sub = stp.mainPart1, dendroLabels = if (stp.includeSampleLabels) NULL else FALSE);

      if (length(adj.removedCovariates) > 0)
      {
        plotSampleTreeFromData(expr.adj[, HE.VST.OR[[ORStep]]$OSR$or.topProbes], scale = or.scaleData, colors.x,
                    file = spaste("sampleTree-Step.", ORStep,  
                                  "-015-NormalizedAndAdjustedDataWithOutlierIndicator"),
                    sub = stp.mainPart2, dendroLabels = if (stp.includeSampleLabels) NULL else FALSE);

        plotPCA.local(expr.adj[, HE.VST.OR[[ORStep]]$OSR$or.topProbes], sampleAnnotation, scale = or.scaleData, 
                    outlierIndicator = outliers + 1,
                    file = spaste("PCA-Step.", ORStep, "-015-NormalizedAndAdjustedData"),
                    sub = stp.mainPart2)
      }

      plotSampleTreeFromData(expr.OSR[, HE.VST.OR[[ORStep]]$OSR$or.topProbes], scale = or.scaleData, phenoColors.OSR,
                  file = spaste("sampleTree-Step.", ORStep, "-020-OutlierSamplesRemoved"),
                  sub = stp.mainPart3,
                  dendroLabels = if (stp.includeSampleLabels) NULL else FALSE);
      plotPCA.local(expr.OSR[, HE.VST.OR[[ORStep]]$OSR$or.topProbes], sampleAnnotation.OSR, scale = or.scaleData, 
                  file = spaste("PCA-Step.", ORStep, "-020-OutlierSamplesRemoved"),
                  pointLabels = rownames(expr.OSR), 
                  sub = stp.mainPart3)
      heCols = HE.VST.OR[[ORStep]]$HEFilter$heColumns;
      or.useProbes = HE.VST.OR[[ORStep]]$OSR$or.topProbes;
      or.useProbes.all = c(1:ncol(xdata.OSR))[ heCols] [or.useProbes];
      or.useProbes.all.logical = c(1:ncol(xdata.OSR))%in% or.useProbes.all;
      changed = nrow(xdata.OSR) < nSamples.OSR;
      nSamples.OSR = nrow(xdata.OSR);
      ORStep = ORStep + 1;
   }
   intermediateResults$HE.VST.OR = HE.VST.OR;
   if (saveOSR)
   {
     dir.create(exprDir.OSR, recursive = TRUE, showWarnings = FALSE);
     write.csv(dataForTable(xdata.OSR, transpose = TRUE, IDcolName = IDcolName.out),
               file = gzfile(file.path(exprDir.OSR, xdataFile.OSR)),
               row.names = FALSE, quote = FALSE);
     write.csv(dataForTable(expr.OSR, transpose = TRUE, IDcolName = IDcolName.out),
               file = gzfile(file.path(exprDir.OSR, exprFile.OSR)),
               row.names = FALSE, quote = FALSE);
     dir.create(sampleDir.OSR, recursive = TRUE, showWarnings = FALSE);
     write.csv(sampleAnnotation.OSR, 
               file = file.path(sampleDir.OSR, sampleFile.OSR),
               row.names = FALSE, quote = TRUE);
     save(expr.OSR, xdata.OSR, 
          sampleAnnotation.OSR, geneAnnotation.he,
          geneAnnotation.entrez,
          file = file.path(exprDir.OSR, dataFile.OSR));
     intermediateResults$Files$exprDir.OSR = exprDir.OSR;
     intermediateResults$Files$xdataFile.OSR = xdataFile.OSR;
     intermediateResults$Files$exprFile.OSR = exprFile.OSR;
     intermediateResults$Files$sampleDir.OSR = sampleDir.OSR;
     intermediateResults$Files$sampleFile.OSR = sampleFile.OSR;
     intermediateResults$Files$dataFile.OSR = dataFile.OSR;
     if (addToAnalysisManager) amAnalysis = am.addMultipleObjectsToAnalysis(
       name = c("OR data for individual analysis", "OR data for WGCNA",
                "OR sample information", "OR feature data for WGCNA"),
       fileName = file.path(exprDir.OSR, dataFile.OSR),
       fileObject = c("xdata.OSR", "expr.OSR", "sampleAnnotation.OSR", "geneAnnotation.he"),
       anaType = amComponent,
       allowNonstandardAnalysisType = TRUE,
       analysis = amAnalysis, verbose = verbose - 1);

   }
   if ("OSR" %in% collapseSteps)
   {
     collapseXData = list(OSR = xdata.OSR);
     collapseExpr = list(OSR = expr.OSR);
   } else {
     collapseXData = list();
     collapseExpr = list();
   }

   if (stopAt=="60-Outlier sample removal")
   {
     if (addToAnalysisManager) am.insertAnalysis(amAnalysis);
     return(intermediateResults);
   }

   #----------------------------------------------------------------------------------------------
   # Removal or replacement of individual outlier measurements
   #----------------------------------------------------------------------------------------------

   haveIORWeights = FALSE;
   if (ior.replace || ior.remove || saveIORWeights)
   {
      if (length(intermediateResults$IOR)==0 || forceCalculation)
      {
        intermediateResults$IOR = list();
        if (verbose > 0) printFlush(spaste(spaces, "Removal/replacement of individual outlier measurements..."))
        if (iorFromHEGenes) 
        {
          inputForIOR = xdata.OSR[, heCols];
        } else 
          inputForIOR = xdata.OSR;
        if (recognizedPipelines[normPipe]=="Counts") {
           iorData = do.call(outlierReplacementData, c(list(
                       counts = inputForIOR,
                       pheno = sampleAnnotation.OSR,
                       design = vst.design,
                       blind = vst.blind,
                       calibrateLowestValue = vst.calibrateLowestValue,
                       groupBy = bw.groupBy,
                       groupsForMinWeightRestriction = bw.groupsForMinWeightRestriction.OSR,
                       DESeqArgs = list(minReplicatesForReplace = Inf, quiet = TRUE),
                       replaceThreshold = ior.replaceThreshold,
                       weightedReplace = ior.weightedReplace,
                       returnWeights = TRUE,
                       returnVSTdata = TRUE,
                       useApproximateVST = useApproximateVST,
                       approximateVST.minShift = approximateVST.minShift,
                       approximateVST.maxShift = approximateVST.maxShift,
                       approximateVST.xWeightScale = approximateVST.xWeightScale),
                       bw.otherArgs));
            names(iorData) = sub("Counts", "Data", names(iorData) );
        } else {
            iorData = do.call(outlierReplacementData.general, c(list(
                       data = inputForIOR,
                       pheno = sampleAnnotation.OSR,
                       useVarianceStabilization = useVarianceStabilization,
                       groupBy = bw.groupBy,
                       groupsForMinWeightRestriction = bw.groupsForMinWeightRestriction.OSR,
                       replaceThreshold = ior.replaceThreshold,
                       weightedReplace = ior.weightedReplace,
                       returnWeights = TRUE,
                       returnVSTdata = TRUE,
                       approximateVST.minShift = approximateVST.minShift,
                       approximateVST.maxShift = approximateVST.maxShift,
                       approximateVST.existingLogBase = approximateVST.existingLogBase,
                       approximateVST.existingLogShift = approximateVST.existingLogShift,
                       approximateVST.xWeightScale = approximateVST.xWeightScale),
                       bw.otherArgs));
        }

        weightsForIOR = iorData$replacementWeights;
        weightFactorsForIOR = iorData$replacementFactors;
        haveIORWeights = TRUE;
        if (saveIORWeights)
        {
           if (verbose > 1) printFlush("    Saving outlier suppression weights...");
           dir.create(IORWeightDir, recursive = TRUE, showWarnings = FALSE);
           write.csv(dataForTable(weightsForIOR, transpose = TRUE, IDcolName = IDcolName.out),
                  file = gzfile(file.path(IORWeightDir, xdataWeightFile.IORep)),
                  row.names = FALSE, quote = FALSE);
           write.csv(dataForTable(weightFactorsForIOR, transpose = TRUE, IDcolName = IDcolName.out),
                  file = gzfile(file.path(IORWeightDir, xdataWeightFactorFile.IORep)),
                  row.names = FALSE, quote = FALSE);
           f = file.path(IORWeightDir, weightFile.IORep);
           save(weightsForIOR, weightFactorsForIOR, file = f);
           if (addToAnalysisManager) amAnalysis = am.addMultipleObjectsToAnalysis(
                name = c("IO replacement weights for data for individual analysis", 
                         "IO replacement weight factors for data for individual analysis"),
                fileName = f,
                fileObject = c("weightsForIOR", "weightFactorsForIOR"),
                anaType = amComponent,
                allowNonstandardAnalysisType = TRUE,
                analysis = amAnalysis, verbose = verbose - 1);

           intermediateResults$Files$IORWeightDir = IORWeightDir;
           intermediateResults$Files$xdataWeightFile.IORep = xdataWeightFile.IORep;
           intermediateResults$Files$weightFile.IORep = weightFile.IORep;
           intermediateResults$Files$xdataWeightFactorFile.IORep = xdataWeightFactorFile.IORep
        }
        
        xdata.IORemoved = iorData$replacedData;
        expr.IORemoved = iorData$replacedVST;
        remove = weightsForIOR <= ior.removeThreshold;
        xdata.IORemoved[remove] = NA;
        expr.IORemoved[remove] = NA;
        heGenes.IOR = do.call(heColumnFnc, 
                           c(list(x = xdata.IORemoved, minValue = minValue, minProportion = minProportion),
                             heColumnOtherArgs));
        if (restrictIORExprToHE) expr.IORemoved = expr.IORemoved[, heGenes.IOR];
        intermediateResults$IOR$heGenes.IOR = heGenes.IOR;
        xdata.IOReplaced = if (restrictIORXDataToHE) 
                               iorData$replacedData[, heGenes.IOR] else iorData$replacedData ;

        expr.IOReplaced = if (restrictIORExprToHE) iorData$replacedVST[, heGenes.IOR] else iorData$replacedVST;
        or.useProbes.IOR0 = if (iorFromHEGenes) or.useProbes else or.useProbes.all.logical;
        or.useProbes.IOR = if (restrictIORExprToHE) or.useProbes.IOR0[heGenes.IOR] else or.useProbes.IOR0;
        if (is.logical(or.useProbes.IOR)) or.useProbes.IOR = which(or.useProbes.IOR);

        intermediateResults$IOR$weightsForIOR = weightsForIOR;
        intermediateResults$IOR$weightFactorsForIOR = weightFactorsForIOR;
        geneAnnotation.IOR.expr0 = if (iorFromHEGenes) geneAnnotation.entrez[heCols, ] else
                                          geneAnnotation.entrez;
        geneAnnotation.IOR.expr = if (restrictIORExprToHE) geneAnnotation.IOR.expr0[heGenes.IOR, ] else 
                                  geneAnnotation.IOR.expr0;
        geneAnnotation.IOR.xdata = if (restrictIORExprToHE) geneAnnotation.IOR.expr0[heGenes.IOR, ] else 
                                  geneAnnotation.IOR.expr0;
        intermediateResults$IOR$geneAnnotation.IOR.expr = geneAnnotation.IOR.expr
        intermediateResults$IOR$geneAnnotation.IOR.xdata = geneAnnotation.IOR.xdata

        if (ior.replace)
        {
          intermediateResults$IOR$xdata.IOReplaced = xdata.IOReplaced;
          intermediateResults$IOR$expr.IOReplaced = expr.IOReplaced;
          plotSampleTreeFromData(expr.IOReplaced[, or.useProbes.IOR], scale = or.scaleData, phenoColors.OSR, 
                  file = "sampleTree-030-IndividualOutliersReplaced",
                  sub = "VST data, individual outliers removed");
          if ("IOReplaced" %in% collapseSteps) 
          {
            collapseXData = c(collapseXData, list(IOReplaced = xdata.IOReplaced));
            collapseExpr = c(collapseExpr, list(IOReplaced = expr.IOReplaced));
          }
          if (saveIOR)
          {
             dir.create(exprDir.IORep, recursive = TRUE, showWarnings = FALSE);
             write.csv(dataForTable(xdata.IOReplaced, transpose = TRUE, IDcolName = IDcolName.out),
                  file = gzfile(file.path(exprDir.IORep, xdataFile.IORep)),
                  row.names = FALSE, quote = FALSE);
             write.csv(dataForTable(expr.IOReplaced, transpose = TRUE, IDcolName = IDcolName.out),
                  file = gzfile(file.path(exprDir.IORep, exprFile.IORep)),
                  row.names = FALSE, quote = FALSE);
             save(expr.IOReplaced, xdata.IOReplaced, 
                  sampleAnnotation.OSR, geneAnnotation.IOR.xdata, 
                  geneAnnotation.IOR.expr,
                  file = file.path(exprDir.IORep, dataFile.IORep));
             intermediateResults$IOR$exprDir.IORep = exprDir.IORep;
             intermediateResults$IOR$xdataFile.IORep = xdataFile.IORep;
             intermediateResults$IOR$exprFile.IORep = exprFile.IORep;
             intermediateResults$IOR$dataFile.IORep = dataFile.IORep;
             if (addToAnalysisManager) amAnalysis = am.addMultipleObjectsToAnalysis(
                name = c("IO replaced data for individual analysis", "IO replaced data for WGCNA",
                         "OR sample information",
                         "IO removed feature data for individual analysis",
                         "IO removed feature data for WGCNA"),
                fileName = file.path(exprDir.IORep, dataFile.IORep),
                fileObject = c("xdata.IOReplaced", "expr.IOReplaced", "sampleAnnotation.OSR",
                               "geneAnnotation.IOR.xdata", "geneAnnotation.IOR.expr"),
                anaType = amComponent,
                allowNonstandardAnalysisType = TRUE,
                analysis = amAnalysis, verbose = verbose - 1);

             intermediateResults$Files$exprDir.IORep = exprDir.IORep;
             intermediateResults$Files$xdataFile.IORep = xdataFile.IORep;
             intermediateResults$Files$exprFile.IORep = exprFile.IORep;
             intermediateResults$Files$dataFile.IORep = dataFile.IORep;
          }
        }

        if (ior.remove)
        {
          if (restrictIORXDataToHE) xdata.IORemoved = xdata.IORemoved[, heGenes.IOR];
          intermediateResults$IOR$xdata.IORemoved = xdata.IORemoved;
          intermediateResults$IOR$expr.IORemoved = expr.IORemoved;
          plotSampleTreeFromData(expr.IORemoved[, or.useProbes.IOR], scale = or.scaleData, phenoColors.OSR, 
                  file = "sampleTree-031-IndividualOutliersRemoved",
                  sub = "VST data, individual outliers removed");
          if ("IORemoved" %in% collapseSteps) 
          {
            collapseXData = c(collapseXData, list(IORemoved = xdata.IORemoved));
            collapseExpr = c(collapseExpr, list(IORemoved = expr.IORemoved));
          }

          if (saveIOR)
          {
             dir.create(exprDir.IORem, recursive = TRUE, showWarnings = FALSE);
             printFlush(file.path("Saving IORem data into file ", dataFile.IORem));
             write.csv(dataForTable(xdata.IORemoved, transpose = TRUE, IDcolName = IDcolName.out),
                  file = gzfile(file.path(exprDir.IORem, xdataFile.IORem)),
                  row.names = FALSE, quote = FALSE);
             write.csv(dataForTable(expr.IORemoved, transpose = TRUE, IDcolName = IDcolName.out),
                  file = gzfile(file.path(exprDir.IORem, exprFile.IORem)),
                  row.names = FALSE, quote = FALSE);
             save(expr.IORemoved, xdata.IORemoved, 
                  sampleAnnotation.OSR, geneAnnotation.IOR.xdata, geneAnnotation.IOR.expr,
                  file = file.path(exprDir.IORem, dataFile.IORem));
             intermediateResults$IOR$exprDir.IORem = exprDir.IORem;
             intermediateResults$IOR$xdataFile.IORem = xdataFile.IORem;
             intermediateResults$IOR$exprFile.IORem = exprFile.IORem;
             intermediateResults$IOR$dataFile.IORem = dataFile.IORem;
             if (addToAnalysisManager) amAnalysis = am.addMultipleObjectsToAnalysis(
                name = c("IO removed data for individual analysis", "IO removed data for WGCNA",
                         "OR sample information",
                         "IO removed feature data for individual analysis",
                         "IO removed feature data for WGCNA"),
                fileName = file.path(exprDir.IORem, dataFile.IORem),
                fileObject = c("xdata.IORemoved", "expr.IORemoved", "sampleAnnotation.OSR",
                               "geneAnnotation.IOR.xdata", "geneAnnotation.IOR.expr"),
                anaType = amComponent,
                allowNonstandardAnalysisType = TRUE,
                analysis = amAnalysis, verbose = verbose - 1);
             intermediateResults$Files$exprDir.IORem = exprDir.IORem;
             intermediateResults$Files$xdataFile.IORem = xdataFile.IORem
             intermediateResults$Files$exprFile.IORem = exprFile.IORem;
             intermediateResults$Files$dataFile.IORep = dataFile.IORep;
          }
        }
        forceCalculation = TRUE;
      } else {
        weightsForIOR = intermediateResults$IOR$weightsForIOR;
        weightFactorsForIOR = intermediateResults$IOR$weightFactorsForIOR;
        
        haveIORWeights = !is.null(weightsForIOR) && !is.null(weightFactorsForIOR)
        heGenes.IOR = intermediateResults$IOR$heGenes.IOR;
        geneAnnotation.IOR.expr = intermediateResults$IOR$geneAnnotation.IOR.expr;
        geneAnnotation.IOR.xdata = intermediateResults$IOR$geneAnnotation.IOR.xdata;
        if (ior.replace)
        {
          xdata.IOReplaced = intermediateResults$IOR$xdata.IOReplaced
          expr.IOReplaced = intermediateResults$IOR$expr.IOReplaced;
          if ("IOReplaced" %in% collapseSteps) 
          {
            collapseXData = c(collapseXData, list(IOReplaced = xdata.IOReplaced));
            collapseExpr = c(collapseExpr, list(IOReplaced = expr.IOReplaced));
          }
        }
        if (ior.remove)
        {
          xdata.IORemoved = intermediateResults$IOR$xdata.IORemoved
          expr.IORemoved = intermediateResults$IOR$expr.IORemoved
          if ("IORemoved" %in% collapseSteps)
          {
            collapseXData = c(collapseXData, list(IORemoved = xdata.IORemoved));
            collapseExpr = c(collapseExpr, list(IORemoved = expr.IORemoved));
          }
        }
      }

      if (stopAt=="80-Individual outlier removal")
      {
        if (addToAnalysisManager) am.insertAnalysis(amAnalysis);
        return(intermediateResults);
      }

      #----------------------------------------------------------------------------------------------
      # Adjustment of IOR data
      #----------------------------------------------------------------------------------------------

      if (forceCalculation || length(intermediateResults$IOR.adj)==0)
      {
        intermediateResults$IOR.adj = list();
        if (verbose > 0) printFlush(spaste(spaces, "Adjusting IOR data..."));
        data1 = list(replaced = if (ior.replace) expr.IOReplaced else NULL,
                     removed = if (ior.remove) expr.IORemoved else NULL);

        if (adj.removeNPCs > 0)
        {
          PCs.IOR = lapply(data1, calculatePCs, nPCs = adj.removeNPCs, scale = scaleForSVD,
                                          PCColnamePrefix = PCColnamePrefix)
          for (type in c("replaced", "removed")[c(ior.replace, ior.remove)]) 
          {
             amAnalysis = conditionalWriteAndAddToAnalysis(
                 doSave = saveIOR.adj || (!is.null(cl.probe2gene) && any(grepl("IOR", collapseSteps))),
                 object = dataForTable(PCs.IOR[[type]], transpose =FALSE, IDcolName = "Sample"),
                 dir = sampleDir.OSR, 
                 fileName = if (type=="removed") PCFile.IORem else PCFile.IORep, 
                 amAnalysis = amAnalysis, 
                 fileReadArgs = list(row.names = "Sample"),
                 amName = spaste("principal components of IO ", type, " data for WGCNA"), 
                 amComponent = amComponent,
                 verbose = verbose-1, indent = indent + 1);
             if (saveIOR.adj) 
             {
               intermediateResults$Files$PCDir.IOR = sampleDir.OSR;
               if (type=="removed") intermediateResults$Files$PCFile.IORem = PCFile.IORem else
                     intermediateResults$Files$PCFile.IORep = PCFile.IORep;
             }
          }
        } else
          PCs.IOR = vector(mode = "list", length = 2);

        data.adj = mymapply(adjustCovariates.EBLM,
                     data1, 
                     PCs = PCs.IOR,
                  MoreArgs = list(
                     sampleAnnotation.OSR,
                     removedCovariates = adj.removedCovariates,
                     retainedCovariates = adj.retainedCovariates,
                     fitToSamples = adj.fitToSamples,
                     useBicovWeights = adj.useBicovWeights,
                     bw.groupBy = bw.groupBy,
                     bw.groupsForMinWeightRestriction = bw.groupsForMinWeightRestriction.OSR,
                     bw.otherArgs = bw.otherArgs,
                     method = adj.useMethod,
                     verbose = verbose - 1, indent = indent + 1));
        gc();
        expr.IORep.adj = data.adj$replaced;
        expr.IORem.adj = data.adj$removed;
        if (!is.null(expr.IORem.adj)) 
          plotSampleTreeFromData(expr.IORem.adj[, or.useProbes.IOR], scale = or.scaleData, phenoColors.OSR, 
                  file = "sampleTree-041-IORemoved-adjusted",
                  sub = "VST data, IO removed, adjusted");
        if (!is.null(expr.IORep.adj))
          plotSampleTreeFromData(expr.IORep.adj[, or.useProbes.IOR], scale = or.scaleData, phenoColors.OSR, 
                  file = "sampleTree-040-IOReplaced-adjusted",
                  sub = "VST data, IO replaced, adjusted");
        intermediateResults$IOR.adj$expr.IORep.adj = expr.IORep.adj;
        intermediateResults$IOR.adj$expr.IORem.adj = expr.IORem.adj;
        if (ior.remove && "IORem.adj" %in% collapseSteps) 
           collapseExpr = c(collapseExpr, list(IORem.adj = expr.IORem.adj));
        if (ior.replace && "IORep.adj" %in% collapseSteps) 
           collapseExpr = c(collapseExpr, list(IORep.adj = expr.IORep.adj));
        if (saveIOR.adj)
        {
           if (ior.remove) 
           {
             printFlush(file.path("Saving IORem.adj data into file ", dataFile.IORem.adj));
	     dir.create(exprDir.IORem.adj, recursive = TRUE, showWarnings = FALSE);
	     write.csv(dataForTable(expr.IORem.adj, transpose = TRUE, IDcolName = IDcolName.out),
		file = gzfile(file.path(exprDir.IORem.adj, exprFile.IORem.adj)),
		row.names = FALSE, quote = FALSE);
	     save(expr.IORem.adj, xdata.IORemoved, 
		sampleAnnotation.OSR,
		file = file.path(exprDir.IORem.adj, dataFile.IORem.adj));
	     intermediateResults$IOR.adj$exprDir.IORem.adj = exprDir.IORem.adj;
	     intermediateResults$IOR.adj$exprFile.IORem.adj = exprFile.IORem.adj;
	     intermediateResults$IOR.adj$dataFile.IORem.adj = dataFile.IORem.adj;
             amAnalysis = am.addMultipleObjectsToAnalysis(
                name = c("adjusted IO removed data for individual analysis", "adjusted IO removed data for WGCNA",
                         "OR sample information",
                         "IO removed feature data for WGCNA",
                         "IO removed feature data for individual analysis"),
                fileName = file.path(exprDir.IORem, dataFile.IORem),
                fileObject = c("xdata.IORem.adj", "expr.IORem.adj", "sampleAnnotation.OSR",
                               "geneAnnotation.IOR.expr", "geneAnnotation.IOR.xdata"),
                anaType = amComponent,
                allowNonstandardAnalysisType = TRUE,
                analysis = amAnalysis, verbose = verbose - 1);
              intermediateResults$Files$exprDir.IORem.adj = exprDir.IORem.adj;
              intermediateResults$Files$exprFile.IORem.adj = exprFile.IORem.adj;
              intermediateResults$Files$dataFile.IORem = dataFile.IORem;
           }
           if (ior.replace)
           {
             dir.create(exprDir.IORep.adj, recursive = TRUE, showWarnings = FALSE);
             write.csv(dataForTable(expr.IORep.adj, transpose = TRUE, IDcolName = IDcolName.out),
		file = gzfile(file.path(exprDir.IORep.adj, exprFile.IORep.adj)),
		row.names = FALSE, quote = FALSE);
	     save(expr.IORep.adj, xdata.IOReplaced, 
		sampleAnnotation.OSR, geneAnnotation.IOR.expr, geneAnnotation.IOR.xdata,
		file = file.path(exprDir.IORep.adj, dataFile.IORep.adj));
	     intermediateResults$IOR.adj$exprDir.IORep.adj = exprDir.IORep.adj;
	     intermediateResults$IOR.adj$exprFile.IORep.adj = exprFile.IORep.adj;
	     intermediateResults$IOR.adj$dataFile.IORep.adj = dataFile.IORep.adj;
             amAnalysis = am.addMultipleObjectsToAnalysis(
                name = c("adjusted IO replaced data for individual analysis", 
                         "adjusted IO replaced data for WGCNA",
                         "OR sample information",
                         "IO removed feature data for WGCNA",
                         "IO removed feature data for individual analysis"),
                fileName = file.path(exprDir.IORep, dataFile.IORep),
                fileObject = c("xdata.IORep.adj", "expr.IORep.adj", "sampleAnnotation.OSR",
                               "geneAnnotation.IOR.expr", "geneAnnotation.IOR.xdata"),
                anaType = amComponent,
                allowNonstandardAnalysisType = TRUE,
                analysis = amAnalysis, verbose = verbose - 1);
              intermediateResults$Files$exprDir.IORep.adj = exprDir.IORep.adj;
              intermediateResults$Files$exprFile.IORep.adj = exprFile.IORep.adj;
              intermediateResults$Files$dataFile.IORep = dataFile.IORep;
           }
	}
        forceCalculation = TRUE;
      } else {
        expr.IORep.adj = intermediateResults$IOR.adj$expr.IORep.adj;
        expr.IORem.adj = intermediateResults$IOR.adj$expr.IORem.adj;
        if (ior.remove && "IORem.adj" %in% collapseSteps) 
           collapseExpr = c(collapseExpr, list(IORem.adj = expr.IORem.adj));
        if (ior.replace && "IORep.adj" %in% collapseSteps) 
           collapseExpr = c(collapseExpr, list(IORep.adj = expr.IORep.adj));
      }

      if (stopAt=="90-Adjustment of IOR")
      {
        if (addToAnalysisManager) am.insertAnalysis(amAnalysis);
        return(intermediateResults);
      }
   } 

   #----------------------------------------------------------------------------------------------
   # Collapsing of data
   #----------------------------------------------------------------------------------------------

   if ((length(cl.probe2gene) > 0) && (forceCalculation || length(intermediateResults$collapse)==0))
   {
     # Figure out the representative for each gene.
     # Run collapseRows on the last entry in collapseXData and collapseExpr.
     intermediateResults$collapse = list();
     if (length(collapseExpr)==0) stop("Internal error: length of 'collapseExpr' is 0.");

     dataForBaseCollapse = list(expr = t(collapseExpr[[ length(collapseExpr)]]),
                                xdata = if (length(collapseXData)==0) NULL else
                                           t(collapseXData[[length(collapseXData)]]));
     dataForBaseCollapse = dataForBaseCollapse[ sapply(dataForBaseCollapse, length) > 0];
     baseProbeIDs = lapply(dataForBaseCollapse, rownames);
     baseGenes = lapply(baseProbeIDs, translateUsingTable, cl.probe2gene);
     #keepRows = lapply(baseGenes, function(x) !is.na(x));
     #dataForBaseCollapse = mymapply(function(x, i) x[i, ], dataForBaseCollapse, keepRows);
     if (verbose > 0) printFlush(spaste(spaces, "Running collapseRows"));
     baseCollapse = mymapply(collapseRows, dataForBaseCollapse, 
                      rowGroup = baseGenes,
                      rowID = baseProbeIDs,
                      method = cl.method);

     if (length(dataForBaseCollapse)==1)
     { 
       representatives = baseCollapse[[1]]$group2row
     } else {
       allGenes = multiUnion(baseGenes);
       allGenes = allGenes[!is.na(allGenes)];
       #allGenes = allGenes[allGenes!=""]
       if (any( !allGenes %in% baseCollapse[[2]]$group2row[, 1]))
         stop("Internal error. Some of 'allGenes' are not among genes in 'baseCollapse[[2]]$group2row'.");
       representatives = baseCollapse[[2]]$group2row;
       exprGenes2xdataGenes = match(baseCollapse[[1]]$group2row[, 1], representatives[, 1]);
       representatives[ exprGenes2xdataGenes, 2] = baseCollapse[[1]]$group2row[, 2];
     }

     collapsedExpr = lapply(collapseExpr, function(x)
     {
       repre2x = match(representatives[, 2], colnames(x));
       keep = is.finite(repre2x);
       out = x[ , repre2x[keep]];
       setColnames(out, representatives[keep, 1]);
     });

     collapsedXData = lapply(collapseXData, function(x)
     {
       repre2x = match(representatives[, 2], colnames(x));
       keep = is.finite(repre2x);
       out = x[ , repre2x[keep]];
       setColnames(out, representatives[keep, 1]);
     });
     intermediateResults$Files = list();
     if (saveCollapsedIORWeights)
     {
       if (recalculateCollapsedWeights || !haveIORWeights)
       {
         if (verbose > 0) printFlush(spaste(spaces, "Collapsed individual outlier weights..."));
         if (recognizedPipelines[normPipe]=="Counts") {
            iorData.collapsed = do.call(outlierReplacementData, c(list(
                        counts = collapsedXData[["OSR"]],
                        pheno = sampleAnnotation.OSR,
                        design = vst.design,
                        blind = vst.blind,
                        calibrateLowestValue = vst.calibrateLowestValue,
                        groupBy = bw.groupBy,
                        groupsForMinWeightRestriction = bw.groupsForMinWeightRestriction.OSR,
                        DESeqArgs = list(minReplicatesForReplace = Inf, quiet = TRUE),
                        replaceThreshold = ior.replaceThreshold,
                        weightedReplace = ior.weightedReplace,
                        useApproximateVST = useApproximateVST,
                        approximateVST.minShift = approximateVST.minShift,
                        approximateVST.maxShift = approximateVST.maxShift,
                        approximateVST.xWeightScale = approximateVST.xWeightScale,
                        returnWeights = TRUE,
                        returnVSTdata = TRUE),
                        bw.otherArgs));
         } else {
             iorData.collapsed = do.call(outlierReplacementData.general, c(list(
                        data = collapsedXData[["OSR"]],
                        pheno = sampleAnnotation.OSR,
                        useVarianceStabilization = useVarianceStabilization,
                        groupBy = bw.groupBy,
                        groupsForMinWeightRestriction = bw.groupsForMinWeightRestriction.OSR,
                        vsnArguments = list(calib = vst.calib),
                        replaceThreshold = ior.replaceThreshold,
                        weightedReplace = ior.weightedReplace,
                        useApproximateVST = useApproximateVST,
                        approximateVST.minShift = approximateVST.minShift,
                        approximateVST.maxShift = approximateVST.maxShift,
                        approximateVST.xWeightScale = approximateVST.xWeightScale,
                        returnWeights = TRUE,
                        returnVSTdata = TRUE),
                        bw.otherArgs));
         }
         weightsForIOR.collapsed = iorData.collapsed$replacementWeights;
         weightFactorsForIOR.collapsed = iorData.collapsed$replacementFactors;
       } else {
         repre2x = match(representatives[, 2], colnames(weightsForIOR));
         keep = is.finite(repre2x);
         weightsForIOR.collapsed = setColnames(weightsForIOR[, repre2x[keep]], representatives[keep, 1]);
         weightFactorsForIOR.collapsed = 
             setColnames(weightFactorsForIOR[, repre2x[keep]], representatives[keep, 1]);
       }
       # Save the weights and weight factors
       dir.create(IORWeightDir.collapsed, recursive = TRUE, showWarnings = FALSE);
       write.csv(dataForTable(weightsForIOR.collapsed, transpose = TRUE, IDcolName = IDcolName.out),
              file = gzfile(file.path(IORWeightDir.collapsed, xdataWeightFile.IORep.collapsed)),
              row.names = FALSE, quote = FALSE);
       write.csv(dataForTable(weightFactorsForIOR.collapsed, transpose = TRUE, IDcolName = IDcolName.out),
              file = gzfile(file.path(IORWeightDir.collapsed, xdataWeightFactorFile.IORep.collapsed)),
              row.names = FALSE, quote = FALSE);
       f = file.path(IORWeightDir.collapsed, weightFile.IORep.collapsed);
       save(weightsForIOR.collapsed, weightFactorsForIOR.collapsed, file = f);
       amAnalysis = am.addMultipleObjectsToAnalysis(
            name = c("collapsed IO replacement weights for data for individual analysis",
                     "collapsed IO replacement weight factors for data for individual analysis"),
            fileName = f,
            fileObject = c("weightsForIOR.collapsed", "weightFactorsForIOR.collapsed"),
            anaType = amComponent,
            allowNonstandardAnalysisType = TRUE,
            analysis = amAnalysis, verbose = verbose - 1);
       intermediateResults$Files$IORWeightDir.collapsed = IORWeightDir.collapsed;
       intermediateResults$Files$xdataWeightFile.IORep.collapsed = xdataWeightFile.IORep.collapsed;
       intermediateResults$Files$xdataWeightFactorFile.IORep.collapsed = xdataWeightFactorFile.IORep.collapsed;
       intermediateResults$Files$weightFile.IORep.collapsed = weightFile.IORep.collapsed;
     } else {
       weightsForIOR.collapsed = NULL;
       weightFactorsForIOR.collapsed = NULL;
     }
     
     geneAnnotation.entrez.collapsed = 
                   geneAnnotation.entrez[ geneAnnotation.entrez[, 1]%in% representatives[, 2], ];
     geneAnnotation.he.collapsed = 
                   geneAnnotation.he[ geneAnnotation.he[, 1]%in% representatives[, 2], ];
     intermediateResults$collapse$exprList = collapsedExpr;
     intermediateResults$collapse$xdataList = collapsedXData;
     intermediateResults$collapse$weightsForIOR.collapsed = weightsForIOR.collapsed;
     intermediateResults$collapse$weightFactorsForIOR.collapsed = weightFactorsForIOR.collapsed;

     collapsedGeneAnnot = function(geneAnnotation, collData, representatives)
     {
       reprProbe = translateUsingTable(colnames(collData), representatives);
       geneAnnotation[ match(reprProbe, geneAnnotation[, 1]), ];
     }

     if (saveCollapsed)
     {
       saveList = character(0);
       amNames = character(0);
       if ("OSR" %in% names(collapsedXData))
       {
         dir.create(exprDir.OSR.collapsed, recursive = TRUE, showWarnings = FALSE);
         xdata.OSR.collapsed = collapsedXData$OSR;
         geneAnnot.xdata.OSR.collapsed = 
             collapsedGeneAnnot( geneAnnotation.entrez, xdata.OSR.collapsed, representatives)
         write.csv(dataForTable(xdata.OSR.collapsed, transpose = TRUE, IDcolName = IDcolName.out.collapsed),
                   file = gzfile(file.path(exprDir.OSR.collapsed, xdataFile.OSR.collapsed)),
                   row.names = FALSE, quote = FALSE);
         write.csv.nr(geneAnnot.xdata.OSR.collapsed,
                   file = gzfile(file.path(exprDir.OSR.collapsed, geneAnnotFile.xdata.OSR.collapsed)))
         saveList = c(saveList, "xdata.OSR.collapsed", "geneAnnot.xdata.OSR.collapsed")
         amNames = c(amNames, "OR collapsed data for individual analysis", 
                              "OR collapsed feature data for individual analysis");
         intermediateResults$Files$exprDir.OSR.collapsed = exprDir.OSR.collapsed;
         intermediateResults$Files$xdataFile.OSR.collapsed = xdataFile.OSR.collapsed;
         intermediateResults$Files$geneAnnotFile.xdata.OSR.collapsed = geneAnnotFile.xdata.OSR.collapsed;
       }
       if ("OSR" %in% names(collapsedExpr))
       {
         dir.create(exprDir.OSR.collapsed, recursive = TRUE, showWarnings = FALSE);
         expr.OSR.collapsed = collapsedExpr$OSR;
         geneAnnot.expr.OSR.collapsed =
             collapsedGeneAnnot( geneAnnotation.entrez, expr.OSR.collapsed, representatives)
         write.csv(dataForTable(expr.OSR.collapsed, transpose = TRUE, IDcolName = IDcolName.out.collapsed),
                   file = gzfile(file.path(exprDir.OSR.collapsed, exprFile.OSR.collapsed)),
                   row.names = FALSE, quote = FALSE);
         write.csv.nr(geneAnnot.expr.OSR.collapsed,
                   file = gzfile(file.path(exprDir.OSR.collapsed, geneAnnotFile.expr.OSR.collapsed)))
         saveList = c(saveList, "expr.OSR.collapsed", "geneAnnot.expr.OSR.collapsed")
         amNames = c(amNames, "OR collapsed data for WGCNA", "OR collapsed feature data for WGCNA");
         intermediateResults$Files$exprDir.OSR.collapsed = exprDir.OSR.collapsed;
         intermediateResults$Files$exprFile.OSR.collapsed = exprFile.OSR.collapsed;
         intermediateResults$Files$geneAnnotFile.expr.OSR.collapsed = geneAnnotFile.expr.OSR.collapsed;
       }
       if (length(saveList) > 0)
       { 
         # Add the sample data to OM once more since it may not have been saved before (e.g. if saveOSR is
         # FALSE) 
         saveList = c(saveList, "sampleAnnotation.OSR");
         amNames = c(amNames, "OR sample information");
         save(list = saveList,
            file = file.path(exprDir.OSR.collapsed, dataFile.OSR.collapsed));
         intermediateResults$collapse$exprDir.OSR.collapsed = exprDir.OSR.collapsed;
         intermediateResults$collapse$dataFile.OSR.collapsed = dataFile.OSR.collapsed;
         if (addToAnalysisManager) amAnalysis = am.addMultipleObjectsToAnalysis(
           name = amNames,
           fileName = file.path(exprDir.OSR.collapsed, dataFile.OSR.collapsed),
           fileObject = saveList,
           anaType = amComponent,
           allowNonstandardAnalysisType = TRUE,
           analysis = amAnalysis, verbose = verbose - 1);
         intermediateResults$Files$dataFile.OSR.collapsed = dataFile.OSR.collapsed;
       }
     }

     if (saveCollapsed && ior.remove)
     {
        saveList = character(0);
        amNames = character(0); 
        if ("IORemoved" %in% names(collapsedXData))
        {
          dir.create(exprDir.IORem.collapsed, recursive = TRUE, showWarnings = FALSE);
          xdata.IORemoved.collapsed = collapsedXData$IORemoved;
          geneAnnot.xdata.IORem.collapsed =
             collapsedGeneAnnot( geneAnnotation.entrez, xdata.IORemoved.collapsed, representatives)
          write.csv(dataForTable(xdata.IORemoved.collapsed, transpose = TRUE, IDcolName = IDcolName.out.collapsed),
             file = gzfile(file.path(exprDir.IORem.collapsed, xdataFile.IORem.collapsed)),
             row.names = FALSE, quote = FALSE);
          saveList = c(saveList, "xdata.IORemoved.collapsed", "geneAnnot.xdata.IORem.collapsed");
          amNames = c(amNames, "IO removed collapsed data for individual analysis", 
                         "IO removed collapsed feature data for individual analysis");
          intermediateResults$Files$exprDir.IORem.collapsed = exprDir.IORem.collapsed;
          intermediateResults$Files$xdataFile.IORem.collapsed = xdataFile.IORem.collapsed;
        }
        if ("IORemoved" %in% names(collapsedExpr))
        {
          dir.create(exprDir.IORem.collapsed, recursive = TRUE, showWarnings = FALSE);
          expr.IORemoved.collapsed = collapsedExpr$IORemoved;
          geneAnnot.expr.IORem.collapsed =
             collapsedGeneAnnot( geneAnnotation.entrez, expr.IORemoved.collapsed, representatives)
          write.csv(dataForTable(expr.IORemoved.collapsed, transpose = TRUE, IDcolName = IDcolName.out.collapsed),
             file = gzfile(file.path(exprDir.IORem.collapsed, exprFile.IORem.collapsed)),
             row.names = FALSE, quote = FALSE);
          saveList = c(saveList, "expr.IORemoved.collapsed", "geneAnnot.expr.IORem.collapsed");
          amNames = c(amNames, "IO removed collapsed data for WGCNA", 
                        "IO removed collapsed feature data for WGCNA");
          intermediateResults$Files$exprDir.IORem.collapsed = exprDir.IORem.collapsed;
          intermediateResults$Files$exprFile.IORem.collapsed = exprFile.IORem.collapsed;
        }
        if (length(saveList) > 0)
        {
          saveList = c(saveList, "sampleAnnotation.OSR");
          amNames = c(amNames, "OR sample information");
          save(list = saveList, 
               file = file.path(exprDir.IORem.collapsed, dataFile.IORem.collapsed));

          intermediateResults$collapse$exprDir.IORem.collapsed = exprDir.IORem.collapsed;
          intermediateResults$collapse$xdataFile.IORem.collapsed = xdataFile.IORem.collapsed;
          intermediateResults$collapse$exprFile.IORem.collapsed = exprFile.IORem.collapsed;
          intermediateResults$collapse$dataFile.IORem.collapsed = dataFile.IORem.collapsed;
          if (addToAnalysisManager) amAnalysis = am.addMultipleObjectsToAnalysis(
             name = amNames,
             fileName = file.path(exprDir.IORem.collapsed, dataFile.IORem.collapsed),
             fileObject = saveList,
             anaType = amComponent,
             allowNonstandardAnalysisType = TRUE,
             analysis = amAnalysis, verbose = verbose - 1);
          intermediateResults$Files$dataFile.IORem.collapsed = dataFile.IORem.collapsed;
        }
     }

     if (saveCollapsed && ior.replace)
     {
        saveList = character(0);
        amNames = character(0); 
        if ("IOReplaced" %in% names(collapsedXData))
        {
          dir.create(exprDir.IORep.collapsed, recursive = TRUE, showWarnings = FALSE);
          xdata.IOReplaced.collapsed = collapsedXData$IOReplaced;
          geneAnnot.xdata.IORem.collapsed =
             collapsedGeneAnnot( geneAnnotation.entrez, xdata.IOReplaced.collapsed, representatives)
          write.csv(dataForTable(xdata.IOReplaced.collapsed, transpose = TRUE, IDcolName = IDcolName.out.collapsed),
             file = gzfile(file.path(exprDir.IORep.collapsed, xdataFile.IORep.collapsed)),
             row.names = FALSE, quote = FALSE);
          saveList = c(saveList, "geneAnnot.xdata.IORem.collapsed", "xdata.IOReplaced.collapsed");
          amNames = c(amNames, "IO removed collapsed feature data for individual analysis",
                               "IO replaced collapsed data for individual analysis");
          intermediateResults$Files$exprDir.IORep.collapsed = exprDir.IORep.collapsed;
          intermediateResults$Files$xdataFile.IORep.collapsed = xdataFile.IORep.collapsed;
        }
        if ("IOReplaced" %in% names(collapsedExpr))
        {
          dir.create(exprDir.IORep.collapsed, recursive = TRUE, showWarnings = FALSE);
          expr.IOReplaced.collapsed = collapsedExpr$IOReplaced;
          geneAnnot.expr.IORem.collapsed =
             collapsedGeneAnnot( geneAnnotation.entrez, expr.IOReplaced.collapsed, representatives)
          write.csv(dataForTable(expr.IOReplaced.collapsed, transpose = TRUE, IDcolName = IDcolName.out.collapsed),
             file = gzfile(file.path(exprDir.IORep.collapsed, exprFile.IORep.collapsed)),
             row.names = FALSE, quote = FALSE);
          saveList = c(saveList, "expr.IOReplaced.collapsed", "geneAnnot.expr.IORem.collapsed");
          amNames = c(amNames, "IO replaced collapsed data for WGCNA", 
                               "IO replaced collapsed feature data for WGCNA");
          intermediateResults$Files$exprDir.IORep.collapsed = exprDir.IORep.collapsed;
          intermediateResults$Files$exprFile.IORep.collapsed = exprFile.IORep.collapsed;
        }
        if (length(saveList) > 0)
        {
          saveList = c(saveList, "sampleAnnotation.OSR");
          amNames = c(amNames, "OR sample information");
          save(list = saveList,
               file = file.path(exprDir.IORep.collapsed, dataFile.IORep.collapsed));
          intermediateResults$collapse$exprDir.IORep.collapsed = exprDir.IORep.collapsed;
          intermediateResults$collapse$xdataFile.IORep.collapsed = xdataFile.IORep.collapsed;
          intermediateResults$collapse$exprFile.IORep.collapsed = exprFile.IORep.collapsed;
          intermediateResults$collapse$dataFile.IORep.collapsed = dataFile.IORep.collapsed;
          if (addToAnalysisManager) amAnalysis = am.addMultipleObjectsToAnalysis(
             name = amNames,
             fileName = file.path(exprDir.IORep.collapsed, dataFile.IORep.collapsed),
             fileObject = saveList,
             anaType = amComponent,
             allowNonstandardAnalysisType = TRUE,
             analysis = amAnalysis, verbose = verbose - 1);
          intermediateResults$Files$dataFile.IORep.collapsed = dataFile.IORep.collapsed;
        }
     }

     if (saveCollapsed && "IORep.adj" %in% names(collapsedExpr))
     {
       dir.create(exprDir.IORep.adj.collapsed, recursive = TRUE, showWarnings = FALSE);
       write.csv(dataForTable(collapsedExpr$IORep.adj, transpose = TRUE, IDcolName = IDcolName.out.collapsed),
          file = gzfile(file.path(exprDir.IORep.adj.collapsed, exprFile.IORep.adj.collapsed)),
          row.names = FALSE, quote = FALSE);
       expr.IORep.adj.collapsed = collapsedExpr$IORep.adj;
       geneAnnot.expr.IORem.collapsed =
             collapsedGeneAnnot( geneAnnotation.entrez, expr.IORep.adj.collapsed, representatives)
       save(expr.IORep.adj.collapsed, sampleAnnotation.OSR, geneAnnot.expr.IORem.collapsed,
          file = file.path(exprDir.IORep.adj.collapsed, dataFile.IORep.adj.collapsed));
       intermediateResults$collapse$exprDir.IORep.adj.collapsed = exprDir.IORep.adj.collapsed;
       intermediateResults$collapse$exprFile.IORep.adj.collapsed = exprFile.IORep.adj.collapsed;
       intermediateResults$collapse$dataFile.IORep.adj.collapsed = dataFile.IORep.adj.collapsed;
       if (addToAnalysisManager) amAnalysis = am.addMultipleObjectsToAnalysis(
          name = c("adjusted IO replaced collapsed data for WGCNA",
                   "IO replaced collapsed feature data for WGCNA",
                   "OR sample information"),  
          fileName = file.path(exprDir.IORep.adj.collapsed, dataFile.IORep.adj.collapsed),
          fileObject = c("expr.IORep.adj.collapsed", "geneAnnot.expr.IORep.collapsed",
                   "sampleAnnotation.OSR"),
          anaType = amComponent,
          allowNonstandardAnalysisType = TRUE,
          analysis = amAnalysis, verbose = verbose - 1);
       intermediateResults$Files$exprDir.IORep.adj.collapsed = exprDir.IORep.adj.collapsed;
       intermediateResults$Files$exprFile.IORep.adj.collapsed = exprFile.IORep.adj.collapsed;
       intermediateResults$Files$dataFile.IORep.adj.collapsed = dataFile.IORep.adj.collapsed;
     }

     if (saveCollapsed && "IORem.adj" %in% names(collapsedExpr))
     {
       dir.create(exprDir.IORem.adj.collapsed, recursive = TRUE, showWarnings = FALSE);
       write.csv(dataForTable(collapsedExpr$IORem.adj, transpose = TRUE, IDcolName = IDcolName.out.collapsed),
          file = gzfile(file.path(exprDir.IORem.adj.collapsed, exprFile.IORem.adj.collapsed)),
          row.names = FALSE, quote = FALSE);
       expr.IORem.adj.collapsed = collapsedExpr$IORem.adj;
       geneAnnot.expr.IORem.collapsed =
             collapsedGeneAnnot( geneAnnotation.entrez, expr.IORem.adj.collapsed, representatives)
       save(expr.IORem.adj.collapsed, sampleAnnotation.OSR, geneAnnot.expr.IORem.collapsed,
          file = file.path(exprDir.IORem.adj.collapsed, dataFile.IORem.adj.collapsed));
       intermediateResults$IOR.adj$exprDir.IORem.adj.collapsed = exprDir.IORem.adj.collapsed;
       intermediateResults$IOR.adj$exprFile.IORem.adj.collapsed = exprFile.IORem.adj.collapsed;
       intermediateResults$IOR.adj$dataFile.IORem.adj.collapsed = dataFile.IORem.adj.collapsed;
       if (addToAnalysisManager) amAnalysis = am.addMultipleObjectsToAnalysis(
          name = c("adjusted IO removed collapsed data for WGCNA",
                   "IO removed collapsed feature data for WGCNA",
                   "OR sample information"),
          fileName = file.path(exprDir.IORem.adj.collapsed, dataFile.IORem.adj.collapsed),
          fileObject = c("expr.IORem.adj.collapsed", "geneAnnot.expr.IORem.collapsed",
                   "sampleAnnotation.OSR"),
          anaType = amComponent,
          allowNonstandardAnalysisType = TRUE,
          analysis = amAnalysis, verbose = verbose - 1);
       intermediateResults$Files$exprDir.IORem.adj.collapsed = exprDir.IORem.adj.collapsed;
       intermediateResults$Files$exprFile.IORem.adj.collapsed = exprFile.IORem.adj.collapsed;
       intermediateResults$Files$dataFile.IORem.adj.collapsed = dataFile.IORem.adj.collapsed;
     }
     forceCalculation = TRUE;
     if (stopAt=="100-Collapse")
     {
       if (addToAnalysisManager) am.insertAnalysis(amAnalysis);
       return(intermediateResults);
     }

   } else {
     collapsedExpr = intermediateResults$collapse$exprList;
     collapsedXData = intermediateResults$collapse$xdataList;
     weightFactorsForIOR.collapsed = intermediateResults$collapse$weightFactorsForIOR.collapsed;
     weightsForIOR.collapsed = intermediateResults$collapse$weightsForIOR.collapsed;
   }
   if (addToAnalysisManager) am.insertAnalysis(amAnalysis);
   intermediateResults;
}


