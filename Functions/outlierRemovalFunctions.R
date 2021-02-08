#===================================================================================================
#
# Attempt a more involved outlier identification and removal
#
#===================================================================================================

robustScale = function(x, var = TRUE, newMedians = NULL, newMADs = NULL)
{
  x = as.matrix(x);
  medians = colQuantileC(x, 0.5);
  if (is.null(newMedians))
  {
    newMedians = rep(0, ncol(x));
  } else if (ncol(x)!=length(newMedians)) {
     stop("Length of 'newMedians' must equal the number of columns in 'x'.");
  } else {
     if (any(!is.finite(newMedians))) stop("All entries in 'newMedians' must be finite.");
  }

  if (is.null(newMADs))
  {
    newMADs = rep(1, ncol(x));
  } else if (ncol(x)!=length(newMADs)) {
     stop("Length of 'newMADs' must equal the number of columns in 'x'.");
  } else {
     if (any(!is.finite(newMADs))) stop("All entries in 'newMADs' must be finite.");
  }

  if (var)
  {
    mad = apply(x, 2, mad, na.rm = TRUE)
    if (any(mad==0) > 0)
    {
       printFlush("Warning in robustScale: some mad = 0. Will substitute standard deviation.");
       sd = apply(x[, mad==0, drop = FALSE], 2, sd);
       mad[mad ==0 ] = sd;
       mad[mad==0] = 1;
    }
    rsx = (x - matrix(medians, nrow(x), ncol(x), byrow = TRUE)) / 
            matrix(mad, nrow(x), ncol(x), byrow = TRUE) * 
            matrix(newMADs, nrow(x), ncol(x), byrow = TRUE)
  } else 
    rsx = (x - matrix(medians, nrow(x), ncol(x), byrow = TRUE))

  rsx = rsx + matrix(newMedians,  nrow(x), ncol(x), byrow = TRUE)
  rsx;
}

equalizeGroupMeans = function(x, group, trim = 0.2)
{
  nCols = ncol(x);
  if (any(is.na(group))) stop("Some entries of 'group' are missing");
  groupLevels = sort(unique(group));
  nGroups = length(groupLevels)
  means = apply(x, 2, tapply, group, mean, trim = trim, na.rm = TRUE);
  if (any (is.na(means))) stop("Some means are missing.");
  averageMeans = colMeans(means, na.rm = TRUE);
  meanShift =  means-matrix(averageMeans, nrow(means), ncol(means), byrow = TRUE);
  for (g in 1:nGroups)
  {
    inGroup = which(group==groupLevels[g]);
    nInGroup = length(inGroup)
    x[inGroup, ] = x[inGroup, ] - matrix(meanShift[g, ], nrow = nInGroup, ncol = nCols, byrow = TRUE);
  }
  x;
}


outlierSamples = function( expr, Z, simFnc = "dist", simOptions = "use = 'p'")
{
  if (simFnc =="dist")
  {
    distM = as.matrix(dist(expr))
  } else {
    sim = eval(parse(text = paste(simFnc, "(t(expr) ", prepComma(simOptions), ")")));
    distM = 1-sim;
  }
  diag(distM) = NA;
  meanDist = robustScale(colSums(distM, na.rm = TRUE))

  list(indicator = meanDist > Z, actualZ = meanDist);
}

dynamicZ = function(n, dynamic = TRUE, baseZ = 2.5, baseN = 100, minZ = 1.5)
{
  if (dynamic)
  {
    Z = max(minZ, baseZ * log10(n)/log10(baseN));
  } else
    Z = baseZ;

  Z;
}

  
removeOutliersAndBatches = function (
   expr, batchLabels, simFnc = "dist", simOptions = "use = 'p'", 
   nDependentZ = TRUE, baseZ = 2.5, baseN = 100, minZ = 1.5,
   useCombat = TRUE, covariates = NULL, numCovariates = NULL,
   removeBatchOutliers = TRUE,
   removeCombinedOutliers = TRUE,
   prior.plots = TRUE, par.prior = TRUE,
   verbose = 1, indent = 0)
{

  if (any(is.na(expr)))
  {
     printFlush("Warning in removeOutliersAndBatches: imputing missing data.");
     expr = t(impute.knn(t(expr))$data);
  }
  nGenes = ncol(expr);
  nSamples = nrow(expr);
  if (nSamples != length(batchLabels)) 
     stop("Number of samples in expr must equal length of batchLabels.");

  spaces = indentSpaces(indent);

  # Step 1: flag outliers in each sample cluster separately
  outliers = rep(FALSE, nSamples)

  batchZ = rep(0, nSamples);

  if (removeBatchOutliers)
  {
    batchNames = sort(unique(batchLabels));
    nBatches = length(batchNames);
    for (b in 1:nBatches)
    {
      batchSamples = c(1:nSamples) [ batchLabels==batchNames[b] ]
      batchSize = length(batchSamples);
      Z = dynamicZ(batchSize, nDependentZ, baseZ, baseN, minZ)
      o = outlierSamples(expr[batchSamples, ], Z, simFnc = simFnc,
                                                simOptions = simOptions);
      outliers [ batchSamples] = o$indicator;
      batchZ [batchSamples] = o$actualZ;
                                   
      if (verbose > 0 && sum(outliers[batchSamples]) > 0)
       printFlush(paste(spaces, "Step 1: in batch", batchNames[b], 
                        "used Z:", signif(Z, 2), " to remove samples", 
                        paste(rownames(expr)[batchSamples] [outliers[batchSamples] ], collapse = ", ")));
    }
  }

  expr1 = expr[!outliers, ];
  labels1 = batchLabels[!outliers];
  nSamples1 = nrow(expr1);

  printFlush(paste(spaces, "..removing batch effects..."));

  # Step 2: use Combat or robust standardization to standardize all batches to the same mean. 

  if (!is.null(covariates))  
  {
    covariates = as.matrix(covariates);
    mod = model.matrix(~., data = as.data.frame(covariates[!outliers, , drop = FALSE ]))
  } else {
    mod = model.matrix(~1, data = as.data.frame(labels1));
  }

  if (useCombat)
  {
    if (length(unique(labels1)) > 1)
    {
      tmp = sva::ComBat(t(expr1), batch = labels1, mod = mod,
                        prior.plots = prior.plots, par.prior = par.prior);
      stdExpr = t(tmp);
    } else 
      stdExpr = expr1;
  } else {
    # Use robust standardization to standardize all batches to the same mean. 
    # Leave variance  unchanged.
    batchSizes = table(labels1);
    biggestBatch = names(batchSizes)[ which.max(batchSizes) ];
    bbMedians = colQuantileC(expr1[labels1==biggestBatch, ], 0.5);
    stdExpr = expr1;
    for (b in 1:nBatches)
    {
      batchSamples = c(1:nSamples1) [ labels1==batchNames[b] ]
      stdExpr[batchSamples, ] = robustScale( expr1[batchSamples, ], var = FALSE);
    }
  }
  
  # Step 3: remove outliers from the combined data:

  if (removeCombinedOutliers)
  {
    o2 = outlierSamples(stdExpr, dynamicZ(nSamples1, nDependentZ, baseZ, baseN, minZ), 
                                                 simFnc = simFnc, simOptions = simOptions);
    outliers2 = o2$indicator;
    Z.all = o2$actualZ;
    if (verbose > 0 && sum(outliers2) > 0)
       printFlush(paste(spaces, "Step 2: removing samples", 
                        paste(rownames(expr1)[outliers2], collapse = ", ")));
  } else {
    outliers2 = rep(FALSE, nSamples1);
    Z.all = rep(NA, nSamples1);
  }
  
  expr2 = stdExpr[!outliers2, ];

  if (!useCombat)
  {
     # Step 4: restore the medians of the largest batch into all samples
     expr2 = expr2 + matrix(bbMedians, nrow(expr2), nGenes, byrow = TRUE);
  }
  outliers2.inAll = rep(TRUE, nSamples);
  outliers2.inAll[!outliers] = outliers2;
  Z.inAll = rep(0, nSamples);
  Z.inAll[!outliers] = Z.all;
  # That's the final result.
  list(expr = expr2, batchOutliers = outliers, allOutliers = outliers2.inAll, 
                     batchZ = batchZ, combinedZ = Z.inAll, stdExpr = stdExpr);
}




#===============================================================================================
#
# Iterative outlier removal
#
#===============================================================================================

robustScale.1 = function(x)
{
  med = median(x, na.rm = TRUE);
  mad = mad(x);
  (x-med)/mad;
}

iterativeOutlierRemoval.withDetails = function(data, Zcut, 
   nTopProbes = NULL, 
   useProbes = NULL,
   robustScale = FALSE,
   scaleData = FALSE)
{
  n = nrow(data);
  outliers = rep(FALSE, n);
  keepIndex = c(1:n);
  if (is.null(useProbes))
  {
    if (is.null(nTopProbes))
    {
      useProbes = 1:ncol(data);
    } else {
      mean = colMedians(data, na.rm = TRUE);
      useProbes = order(-mean, na.last = TRUE)[1:nTopProbes];
    } 
  }
  done = FALSE;
  details = list();
  while (!done)
  {
    x = if (scaleData) scale(data[, useProbes]) else data[, useProbes];
    dist = as.matrix(dist(x));
    Z.1 = if (robustScale) robustScale.1(colSums(dist)) else scale(colSums(dist));
    outliers.1 = Z.1 > Zcut;
    details = c(details, list(data.frame(sample = rownames(data), Z = Z.1, outliers = outliers.1)));
    if (any(outliers.1))
    {
      outliers[keepIndex] = outliers.1;
      keepIndex = keepIndex[!outliers.1];
      data = data[!outliers.1, ];
    } else done = TRUE
  }
  names(details) = spaste("Step.", 1:length(details));
  list(ORData = data,
       topProbes = useProbes,
       removedOutliers = outliers,
       details = details)
}


iterativeOutlierRemoval = function(data, Zcut, nTopProbes = NULL, robustScale = FALSE,
                           getDetails = FALSE)
{
  lst = iterativeOutlierRemoval.withDetails(data = data, Zcut = Zcut, nTopProbes = nTopProbes, 
                 robustScale = robustScale);
  data = lst$ORData;
  attr(data, "removedOutliers") = lst$removedOutliers;
  if (getDetails) attr(data, "ORDetails") = lst$details;
  data;
}
