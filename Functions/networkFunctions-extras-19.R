# Extras not contained in the R packages.

# Changes in 19:
#   . associationScreening.singleDesign now returns model comparison stats only if there is an empty string ("") entry in
#     termsOfInterest; the statistics (F in coefficient column, p in p column, t calculated from p and DOF) go in that row. The
#     commonStatistics part of the output is dropped.

# Changes to 18:
#   . ablineThresholdFromFDR is now less misleading when none of the FDRs are below the threshold.

# Changes to 17:
#   . Change in dependentColumns makes it more resistant to numeric overflows and it now works with constant non-zero
#     columns.

# Changes to 16:
#   . removing multiGrep, multiGrepl, multiSub, multiGSub which are now part of WGCNA
#   . Changing names in output of marginWidth
#   . outlierReplacementData.general switched from vsn to approximateVST and fix bugs in using VST

# Change to 15:
#   . removing binarizeCategoricaCovariates, binarizeCategoricalColumns and variations that have been moved to
#     the WGCNA package

# Change to 14: 
#   . matrixBicovWeights changed in how it handles grouping variables. Instead of assuming factors and
#   equalizing expression among all combinations, the function will fit a linear model and 
#   the bicov weights will be calculated from residuals of a
#   linear model where data is regressed on the 'pheno' covariates according to 

#   . In-version change: adapting function rearrangeInternalCollections to new tagging HDSigDB in anRichment.

# Change to 13: 
#   . Adapating standardCollections to new anRichment; 
#   . Bugfix in dependentColumns: instead of cov(data) we now use t(data)%*%data. This could change results
#   of some calculations although I don't think this is an issue with any calculation using verion 13.

# Change from 11 to 12: changed arguments to associationScreening to make the input more similar to that of
# test specification for DESeq2.

# Bugg fix: replacementData.general was rounding replacement values for general data which does not make
# sense.


# Change from 10 to 11: fixed a bug in and substantially improved matrixBicovWeights and outlier
# identification and relacement functions that rely on it.
# Also: in removing individual outliers, the VST on the final IOMR data now includes calibration of lowest
# value.

# changing 09 to 10: standardScreeningBinaryTrait switched default for consistentSigns to TRUE.

# reason for changing from 08 to 09: change in kIM and kME that is not backward-compatible.

#--------------------------------------------------------------------------------------------

# Load the requisite libraries


libraryLocs = .libPaths();

loadPackages = c("MASS", "class", "cluster", "impute", "Hmisc",
"matrixStats") #, "mclust");

#loadPackages = c("MASS", "class", "cluster", "impute", "Hmisc", 
#"foreach", "doParallel", "matrixStats",
#"sva", "preprocessCore", "anRichment") #, "mclust");

for (pack in loadPackages)
{
  if (!require(pack, character.only = TRUE, quietly = TRUE))
  {
    ok = FALSE;
    for (loc in libraryLocs) ok = ok | require(pack, lib.loc = loc, 
                                               character.only = TRUE, quietly = TRUE);
    if (!ok) warning(paste("Could not find package", pack, 
                        "either in the standard library nor in the personal libraries\n",
                        paste("     ", libraryLocs, collapse = "\n")));
  }
}


library(WGCNA)

getBaseDir = function()
{ 
 file.path(sub("/Work/.*", "", getwd()), "Work");
}

#source(file.path(getBaseDir(), "CommonFunctions/labelPoints2-01.R"));
options(stringsAsFactors = FALSE, deparse.max.lines = 30, max.print = 1000);

#-----------------------------------------------------------------------------------------------
# Trait selection based on independence and significance
# Assumes that, just like with the PCs, the Traits have the same columns in each dataset (though the
# sample sets need not be the same).

SelectTraits = function(Traits, BranchCut = 0.25, SelectOnSignificance = FALSE, PCs = NULL, 
                        SignifThres = 0.03, Impute = FALSE, verbose = 1, indent = 0)
{
  spaces = indentSpaces(indent);
  if (verbose>0) printFlush(paste(spaces, "SelectTraits: Selecting from ", dim(Traits[[1]]$data)[2],
                                   "traits."));
  No.Sets = length(Traits);
  TDiss = 1-cor(Traits[[1]]$data, use = "pairwise.complete.obs");
  if (No.Sets>1) for (set in 2:No.Sets)
  {
     TDiss = pmax(TDiss, 1-cor(Traits[[set]]$data, use = "pairwise.complete.obs"));
  }
  h = hclust(as.dist(TDiss), method = "average");
  TMods = ModuleNumber(h, CutHeight = BranchCut, MinSize = 1);
  No.TMods = nlevels(as.factor(TMods));
  SelTraits = vector(mode="list", length = No.Sets);
  for (set in 1:No.Sets)
  {
    TData = Traits[[set]]$data;
    TData[is.na(TData)] = 0;
    TPCs = ModulePrinComps1(TData, as.factor(TMods), Impute = Impute, verbose = 0, 
                            GetConformity = TRUE);
    SelTraits[[set]] = list(data = TPCs$PrinComps);
    for (tmod in 1:No.TMods)
    {
      if (sum(TMods==tmod)>1)
      {
        rnk = order(-TPCs$ModuleConformity[TMods==tmod]);
        SelTraits[[set]]$data[, tmod] = (TData[, TMods==tmod])[, rnk[1]];
        names(SelTraits[[set]]$data)[tmod] = (names(TData)[TMods==tmod])[rnk[1]];
      } else {
        SelTraits[[set]]$data[, tmod] = TData[, TMods==tmod];
        names(SelTraits[[set]]$data)[tmod] = names(TData)[TMods==tmod];
      }
    }
  }

  if (verbose>0) printFlush(paste(spaces, "SelectTraits: Clustering led to ", dim(SelTraits[[1]]$data)[2],
                                   "traits."));
 
  if (SelectOnSignificance)
  {
    # Reduce further: calculate cor.tests for each ME with each trait in each set; 
    # keep only traits that have at least one cor.test$p.value below a threshold

    if (is.null(PCs)) stop("PCs must be given when SelectOnSignificance is requested.");

    No.Mods = dim(PCs[[1]]$data)[2];
    if (is.null(No.Mods)) 
       stop("Given PCs do not appear to have the correct structure (vector of list",
            "with \'data\' component being a matrix whose columns are PC vectors");

    No.Traits = dim(SelTraits[[1]]$data)[2];
    
    SelectTrait = rep(FALSE, times = No.Traits);
    
    for (trait in 1:No.Traits)
      for (mod in 1:No.Mods)
      {
        Significant = TRUE;
        for (set in (1:No.Sets))
        {
          ct = cor.test(PCs[[set]]$data[, mod], SelTraits[[set]]$data[, trait]); 
          if (ct$p.value>SignifThres) Significant = FALSE;
        }
        if (Significant) SelectTrait[trait] = TRUE;
      }

    for (set in 1:No.Sets)
    {
      SelTraits[[set]]$data = SelTraits[[set]]$data[, SelectTrait];
    }
    
    # print(paste("No. of selected traits expected by chance:", No.Mods * No.Traits * TraitThres));
  }
  
  # Re-cluster the significant traits for diagnostic purposes
  if (sum(SelectTrait)>1)
  {
    TDiss = 1-cor(SelTraits[[1]]$data, use = "pairwise.complete.obs");
    if (No.Sets>1) for (set in 2:No.Sets)
    {
       TDiss = pmax(TDiss, 1-cor(SelTraits[[set]]$data, use = "pairwise.complete.obs"));
    }
    newh = hclust(as.dist(TDiss), method = "average");
  } else {
    newh = NULL;
  }

  if (verbose>0) 
  {
    printFlush(paste(spaces, "SelectTraits: Selected", sum(SelectTrait), "traits: "));
    printFlush(paste(spaces, paste(names(SelTraits[[1]]$data), collapse = ", ")));
  }
                  

  list(No.SelectedTraits = sum(SelectTrait), Traits = SelTraits, ClusterTree = h, NewClusterTree = newh);
}

#----------------------------------------------------------------------------------------------
# 
# EvaluateClusters
#
#----------------------------------------------------------------------------------------------
# Input: distance matrix and cluster labels. 0 in Labels means the particular point is a
# singleton (not assigned to a cluster).
# Output: a list of cluster quality indicators.

EvaluateClusters = function(DistM, Labels, RefLabels = NULL, Sample = FALSE, SampleProp = 1,
                            verbose = 2, indent = 0)
{

  spaces = indentSpaces(indent);

  d = dim(DistM);
  No.Points = length(Labels);

  if (d[1]!=d[2]) stop("Distance matrix DistM must be square.");

  if (d[1]!=No.Points) stop("Dimension of DistM incompatible with number of cluster labels");

  if (verbose>0) print(paste(spaces, "Evaluating clusters...")); 

  # Calculate indices relating to the internal and external links

  if (Sample)
  {
    if ( (SampleProp<=0) | (SampleProp>=1) )
      stop(paste("Incorrect SampleProp parameter given:", SampleProp));
    No.Sampled = as.integer(No.Points * SampleProp);
    Sampled = sample(x = No.Points, size = No.Sampled);
    DistM = DistM[Sampled, Sampled];
    Labels = Labels[Sampled];
    RefLabels = RefLabels[Sampled];
    No.Points = No.Sampled;
  }

  InternalLinks = matrix(0, nrow = No.Points, ncol = No.Points);
  for (point in 1:No.Points) if (Labels[point]!=0)
    InternalLinks[, point] = ifelse(Labels==Labels[point], Labels, 0);

  IntLinks = as.vector(InternalLinks[upper.tri(InternalLinks)]);
  Dist = as.vector(DistM[upper.tri(DistM)]);
  No.Internal = sum(IntLinks!=0);
  No.External = length(IntLinks) - No.Internal;

  if ( (No.Internal > 0) & (No.External>0) )
  {
     ord = order(Dist, c(1:length(Dist)));

     # ...number of "misplaced" internal and external links

     s1 = Dist[ord[No.Internal]];
     s2 = Dist[ord[No.Internal+1]];
     if (s1==s2)
     {
       MaxInt = s1 + 1e-10;
       MinExt = s1 - 1e-10;
     } else {
       MaxInt = s1;
       MinExt = s2;
     }
     No.Misplaced = sum( ( (IntLinks>0) & (Dist>MaxInt) ) | ( (IntLinks==0) & (Dist<MinExt) ) );

     # ...weight index of the internal and external links

     DistInt = sum(Dist[IntLinks>0]);
     DistExt = sum(Dist[IntLinks==0]);
     DistSm = sum(Dist[ord <= No.Internal]);
     DistLg = sum(Dist[ord > No.Internal]);
     DistIndex = DistInt/DistSm * DistLg/DistExt;
  } else {
     No.Misplaced = 0;
     DistIndex = 0;
  }

  if (!is.null(RefLabels))
  {
    if (length(RefLabels)!=No.Points)
       stop("Length of given RefLabels incompatible with number of points.");
    
    RefInternalLinks = matrix(0, nrow = No.Points, ncol = No.Points);
    for (point in 1:No.Points) if (RefLabels[point]!=0)
      RefInternalLinks[, point] = ifelse(RefLabels==RefLabels[point], RefLabels, 0);

    RefIntLinks = as.vector(RefInternalLinks[upper.tri(RefInternalLinks)]);
 
    No.Agreed = sum( xor(IntLinks>0, RefIntLinks==0));
    No.Disagreed = sum( xor(IntLinks==0, RefIntLinks==0));

    RandIndex = No.Agreed/(No.Agreed+No.Disagreed);

    # Make each unassigned label unique
    XLabels = Labels;
    NUnassigned = sum(XLabels==0);
    start = max(XLabels);
    XLabels[XLabels==0] = start + c(1:NUnassigned);

    # Same for the reference labels 
    XRefLabels = RefLabels;
    NRefUnassigned = sum(XRefLabels==0);
    start = max(XRefLabels);
    XRefLabels[XRefLabels==0] = start + c(1:NRefUnassigned);
 
    AdjRandIndex = adjustedRandIndex(XLabels, XRefLabels);  # This requires package mclust
  } else {
    RandIndex = 0;
    AdjRandIndex = 0;
  }

  list(No.Misplaced = No.Misplaced, DistIndex = DistIndex, RandIndex = RandIndex, 
       AdjRandIndex = AdjRandIndex);
}


# Biweight correlation. The original functions are from
# http://www.unt.edu/benchmarks/archives/2001/december01/rss.htm and that in
# turn seems to be based on Wilcox (1997, page 197).
# The actually useful versions are PL's rewrites that use block code so the
# functions can be used with matrices as inputs. PL's function could be
# improved since it does some calculations unnecessarily twice when y==NULL.

bicov.original<-function(x, y, na.rm = FALSE){
mx <- median(x, na.rm = na.rm)
my <-median(y, na.rm = na.rm)
ux <- abs((x - mx)/(9 * qnorm(0.75) * mad(x, na.rm = na.rm)))
uy <- abs((y - my)/(9 * qnorm(0.75) * mad(y, na.rm = na.rm)))
aval <- ifelse(ux <= 1, 1, 0)
bval <- ifelse(uy <= 1, 1, 0)
top <- sum(aval * (x - mx) * (1 - ux^2)^2 * bval * (y - my) * (1 - uy^2)^2, na.rm = na.rm)
top <- sum(!is.na(x) & !is.na(y)) * top
botx <- sum(aval * (1 - ux^2) * (1 - 5 * ux^2), na.rm = na.rm)
boty <- sum(bval * (1 - uy^2) * (1 - 5 * uy^2), na.rm = na.rm)
bi <- top/(botx * boty)
bi
}

bicor.original<-function(x, y, na.rm = FALSE){
x <-as.numeric(x)
y<-as.numeric(y)
bicov(x,y, na.rm = na.rm)/(sqrt(bicov(x,x, na.rm = na.rm)*bicov(y,y, na.rm = na.rm)))
}

# This function should now work for both vector and matrix x and y (all combinations);
# Big WARNING: The use = pairwise.complete.obs simply works by removing NAs independently in x and
# y!!!!!

# Weighted version of bicov and bicor

bicovw = function(x, xWeight, y = NULL, yWeight = NULL, robustX = TRUE, robustY = TRUE, 
                  use = "all.obs", diagOnly = FALSE)
{
  x = as.matrix(x);
  xWeight = as.matrix(xWeight);
  if (!all.equal(dim(x), dim(xWeight)))
     stop("x and xWeight must have the same dimensions.");

  na.method <- pmatch(use, c("all.obs", "pairwise.complete.obs"))
  if (is.na(na.method)) 
    stop(paste("Unrecognized parameter 'use'. Recognized values are \n",
         "'all.obs', 'pairwise.complete.obs'"))
  na.rm = (na.method==2);
  if (is.null(y)) {
    if (robustX)
    {
      #mx <- apply(x, 2, median, na.rm = na.rm)
      mx = medianw(x, xWeight);
      mxMat = matrix(mx, nrow = nrow(x), ncol = ncol(x), byrow = TRUE);
      madxMat = matrix(apply(x, 2, mad, na.rm = na.rm), nrow = nrow(x), ncol = ncol(x), byrow = TRUE);
      ux <- as.matrix(abs((x - mxMat)/(9 * qnorm(0.75) * madxMat)))
    } else {
      mx = meanw(x, xWeight)
      mxMat = matrix(mx, nrow = nrow(x), ncol = ncol(x), byrow = TRUE);
      ux = matrix(0, nrow = nrow(mxMat), ncol = ncol(mxMat));
    }
    aval  = ifelse(ux <= 1, 1, 0) * xWeight;
    botx <- apply(as.matrix(aval * (1 - ux^2) * (1 - 5 * ux^2)), 2, sum, na.rm = na.rm);
    ux[is.na(ux)] = 1;
    aval[is.na(aval)] = 0;
    if (diagOnly)	
    {
      # diagOnly makes sense particularly for y=NULL because it is used in the
      # normalization in bicor.
      topFact = apply(!is.na(x), 2, sum)
      x[is.na(x)] = mxMat[is.na(x)];
      top <- apply(as.matrix(aval * (x - mxMat) * (1 - ux^2)^2)^2, 2, sum); 
      top <- top * topFact;
      bi <- top/((botx^2) * (1-1/nrow(x)));
    } else { 
      topFact = t(!is.na(x)) %*% !is.na(x) 
      x[is.na(x)] = mxMat[is.na(x)];
      top <- t(as.matrix(aval * (x - mxMat) * (1 - ux^2)^2)) %*% 
               as.matrix((aval * (x - mxMat) * (1 - ux^2)^2))
      top <- top * topFact;
      bi <- top/(botx %o% (botx *(1-1/nrow(x)) ))
    }
    bi
  } else {
    y = as.matrix(y);
    yWeight = as.matrix(yWeight);
    if (robustX)
    {
      #mx <- apply(x, 2, median, na.rm = na.rm)
      mx = medianw(x, xWeight)
      mxMat = matrix(mx, nrow = nrow(x), ncol = ncol(x), byrow = TRUE);
      madxMat = matrix(apply(x, 2, mad, na.rm = na.rm), nrow = nrow(x), ncol = ncol(x), byrow = TRUE);
      ux <- as.matrix(abs((x - mxMat)/(9 * qnorm(0.75) * madxMat)))
    } else {
      mx = meanw(x, xWeight)
      mxMat = matrix(mx, nrow = nrow(x), ncol = ncol(x), byrow = TRUE);
      ux = matrix(0, nrow = nrow(mxMat), ncol = ncol(mxMat));
    }
    if (robustY)
    {
      #my <- apply(y, 2, median, na.rm = na.rm)
      my = medianw(y, yWeight)
      myMat = matrix(my, nrow = nrow(y), ncol = ncol(y), byrow = TRUE);
      madyMat = matrix(apply(y, 2, mad, na.rm = na.rm), nrow = nrow(y), ncol = ncol(y), byrow = TRUE);
      uy <- as.matrix(abs((y - myMat)/(9 * qnorm(0.75) * madyMat)))
    } else {
      my = meanw(y, yWeight);
      myMat = matrix(my, nrow = nrow(y), ncol = ncol(y), byrow = TRUE);
      uy = matrix(0, nrow = nrow(myMat), ncol = ncol(myMat));
    }
    aval <- ifelse(ux <= 1, 1, 0) * xWeight;
    bval <- ifelse(uy <= 1, 1, 0) * yWeight;
    botx <- apply(as.matrix(aval * (1 - ux^2) * (1 - 5 * ux^2)), 2, sum, na.rm = na.rm)
    boty <- apply(as.matrix(bval * (1 - uy^2) * (1 - 5 * uy^2)), 2, sum, na.rm = na.rm)
    ux[is.na(ux)] = 1;
    uy[is.na(uy)] = 1;
    aval[is.na(aval)] = 0;
    bval[is.na(bval)] = 0;
    topFact = t(!is.na(x)) %*% !is.na(y) 
    x[is.na(x)] = mxMat[is.na(x)];
    y[is.na(y)] = myMat[is.na(y)];
    top <- t(as.matrix(aval * (x - mxMat) * (1 - ux^2)^2)) %*% as.matrix((bval * (y - myMat) * (1 - uy^2)^2))
    top <- top * topFact;
    bi <- top/(botx %o% (boty * (1-1/nrow(x))))
    bi
  }
}

# weighted version of biweight mid-corellation
#
bicorw<-function(x, xWeight, y = NULL, yWeight = NULL, robustX = TRUE, robustY = TRUE, use = "all.obs")
{
  na.method <- pmatch(use, c("all.obs", "pairwise.complete.obs"))
  if (is.na(na.method)) 
    stop(paste("Unrecognized parameter 'use'. Recognized values are \n",
         "'all.obs', 'pairwise.complete.obs'"))
  if (!all.equal(dim(x), dim(xWeight)))
    stop("x and xWeight must have the same dimensions");
  if (is.null(y)) 
  {
    bcx = sqrt(bicov(x, robustX = robustX, use = use, diagOnly = TRUE));
    bicov(x, robustX = robustX, use = use)/(bcx %o% bcx);
  } else
  {
    if (!all.equal(dim(y), dim(yWeight)))
       stop("y and yWeight must have the same dimensions");
    bicov(x, y, robustX = robustX, robustY = robustY, use = use)/
                (sqrt(bicov(x, robustX = robustX, use = use, diagOnly = TRUE)) %o%
                 sqrt(bicov(y, robustX = robustY, use = use, diagOnly = TRUE)))
  }
}

meanw = function(x, w, na.rm = TRUE)
{
  x = as.matrix(x);
  w = as.matrix(w);
  if (!all.equal(dim(x), dim(w)))
    stop("Dimensions of x and w must be the same.");

  if (na.rm)
  {
    w[is.na(w)] = 0;
    w[is.na(x)] = 0;
  }
  s = colSums(x*w, na.rm = na.rm);
  sw = colSums(w, na.rm = na.rm);

  if (sum(sw==0)>0) return(NA)

  s/sw;
}

sdw = function(x, w)
{
  x = as.matrix(x);
  w = as.matrix(w);
  if (!all.equal(dim(x), dim(w)))
    stop("Dimensions of x and w must be the same.");

  if (na.rm)
  {
    w[is.na(w)] = 0;
    w[is.na(x)] = 0;
  } 
  s = colSums(x*w, na.rm = TRUE)
  ss = colSums(x*x*w, na.rm = TRUE)
  sw = colSums(w, na.rm = TRUE);
  stop("not finished")
}



medianw = function(x, w, na.rm = FALSE)
{
  x = as.matrix(x);
  w = as.matrix(w);
  if (!all.equal(dim(x), dim(w)))
    stop("Dimensions of x and w must be the same.");

  if (na.rm)
  {
    w[is.na(w)] = 0;
    w[is.na(x)] = 0;
  }

  if (sum(w<0)>0) stop("Weights w must be non-negative.");

  order = apply(x, 2, order, na.last = TRUE);
  for (i in 1:ncol(x)) w[, i] = w[order[, i], i];

  nSamples = nrow(w);
  csw = apply(w, 2, cumsum)
  sw = apply(w, 2, sum);
  if (sum(sw==0)>0)
    stop("Some columns of w have no non-zero entries.")
  # csw are cumulative sums of weights
  csw = csw * matrix(nSamples/csw[nSamples, ], nrow = nSamples, ncol = ncol(w), byrow = TRUE);

  midw = nSamples/2;
  # midh will be the index in each column where csw rises above midw. Can be anywehere between 1 and
  # nSamples.
  midh = apply(csw-midw>0, 2, match, x = TRUE)

  # Couldn't figure out a way to do this in block notation.
  medians = rep(0, ncol(x));
  for (i in 1:ncol(x))
    if (midh[i]==1) {
      medians[i] = x[order[1, i], i];
    } else {
      w0 = abs(csw[midh[i]-1]); x0 = x[order[midh[i]-1, i], i];
      w1 = csw[midh[i]]; x1 = x[order[midh[i], i], i];
      medians[i] = (w0*x1 + w1* x0)/(w1+w0);
    }
  medians;
}

corw = function(x, y = NULL, wx, wy = NULL, verbose = 0)
{
  if (is.null(y))
  {
    y = x;
    wy = wx;
    symm = TRUE;
  } else 
    symm = FALSE;

  nx = ncol(x);
  ny = ncol(y);

  if (is.null(wy)) stop("If 'y' is given, 'wy' must also be explicitly given.");

  if (nrow(x)!=nrow(y)) stop("'x' and 'y' do not have the same number of rows.");

  out = matrix(NA, nx, ny);
  dimnames(out) = list(colnames(x), colnames(y));
  if (verbose > 0) pind = initProgInd();
  for (i in 1:nx) for (j in (if (symm) i else 1):ny)
  { 
    fin = is.finite(x[, i]) & is.finite(y[, j]) & is.finite(wx[, i]) & is.finite(wy[, j])
    x1 = x[fin, i];
    y1 = y[fin, j];
    w1 = wx[fin, i]*wy[fin, j];
    mx = weighted.mean(x1, w1);
    my = weighted.mean(y1, w1);
    x1 = x1-mx;
    y1 = y1-my;
    num = weighted.mean(x1 * y1, w1);
    denom = sqrt( weighted.mean(x1^2, w1) * weighted.mean(y1^2, w1));
    out[i,j] = num/denom;
    if (symm) out[j,i] = num/denom;
    if (verbose > 0) pind = updateProgInd(i/nx, pind);
  }
  if (verbose > 0) printFlush("");
  out;
}
    

propVarExplained = function(datExpr, colors, MEs, corFnc = "cor", corOptions = "use = 'p'")
{
  fc = as.factor(colors);
  mods = levels(fc);
  nMods = nlevels(fc);
  if (nMods!=ncol(MEs))
    stop(paste("Input error: number of distinct 'colors' differs from\n", 
               "the number of module eigengenes given in ME."));

  if (ncol(datExpr)!=length(colors))
    stop("Input error: number of probes (columns) in 'datExpr' differs from the length of goven 'colors'.");

  if (nrow(datExpr)!=nrow(MEs))
    stop("Input error: number of observations (rows) in 'datExpr' and 'MEs' differ.");

  PVE = rep(0, nMods);

  col2MEs = match(mods, substring(names(MEs), 3));

  if (sum(is.na(col2MEs))>0)
    stop("Input error: not all given colors could be matched to names of module eigengenes.");

  for (mod in 1:nMods)
  {
    modGenes = c(1:nGenes)[as.character(colors)==mods[mod]];
    corExpr = parse(text = paste(corFnc, "(datExpr[, modGenes], MEs[, col2MEs[mod]],",
                                 corOptions, ")"));
    PVE[mod] = mean(as.vector(eval(corExpr)^2));
  }

  names(PVE) = paste("PVE", mods, sep = "");
  PVE
}
 
cutreeDynamic.1 = function(hierclust, maxTreeHeight=1, deepSplit=TRUE, minModuleSize=50, 
                           minAttachModuleSize=100, nameallmodules=FALSE, useblackwhite=FALSE)
{
  if ( (minAttachModuleSize!=100) | (nameallmodules!=FALSE) | (useblackwhite!=FALSE) )
    waring("cutreeDynamic.1: Parameters minAttachModuleSize, nameallmodules, useblackwhite are ignored.");
  labels2colors(cutreeDynamicTree(hierclust, maxTreeHeight, deepSplit, minModuleSize));
}

#=================================================================================================
#
# Dendrogram orderings
#
#=================================================================================================


dendrogramOrderings = function(dendro, maxMerge = 18, verbose = 1, indent = 0)
{

  spaces = indentSpaces(indent)

  if (verbose>0) printFlush(paste(spaces, "Calculating dendrogramOrderings.."));
  nMerge = length(dendro$height);

  nSingl = nMerge + 1;
  if (nMerge < 2) return(list(orders = matrix(1:nSingl, nSingl, 1),
                              rankings = matrix(1:nSingl, nSingl, 1)));
  mleft = dendro$merge[,1];
  mright = dendro$merge[,2];

  if ((verbose>0) & (nMerge>maxMerge)) 
     printFlush(paste("..FYI: limiting number of flipped merges to", maxMerge)); 

  # First step: get the list of L and R singletons below each merge, i.e. for each merge and each
  # singleton the array sides[, merge, side] gives singletons to the side side of merge merge.
  # The number of those sigletons is recorded in nSides[merge, side]

  sides = array(0, dim = c(nSingl, ncol = nMerge, 2));

  nSides = matrix(0, nrow = nMerge, ncol = 2);

  traceBack = rep(0, times = nMerge);
  curSide = rep(1, times = nMerge);	# 1=left, 2=right, 3=all done.

  nTrace = 1;		# Number of traces includes current merge.
  traceBack[1] = nMerge;
  curMerge = nMerge;

  if (verbose>1) printFlush(paste(spaces, " ..finding `sides' of each merge.."));
  while ((curSide[nMerge]<3) | (nTrace > 1))
  {
    # print(paste("nTrace:", nTrace, "traceBack:", paste(traceBack, collapse = ", ")));
    # print(paste("curMerge:", curMerge, "curSide", paste(curSide, collapse = ", ")));
    if (curSide[curMerge] < 3)
    { 
      # print(paste("Processing side", curSide[curMerge]));
      # Both sides not yet processed.
      if (dendro$merge[curMerge, curSide[curMerge]] < 0)
      {
        # Have a sigleton on the current side. Record the singleton as being on the current side of the
        # current merge and all merges in the traceback.
        # print("Singleton.");
        singleton = -dendro$merge[curMerge, curSide[curMerge]];
        curSide[curMerge] = curSide[curMerge] + 1;
        iTrace = nTrace;
        while (iTrace > 0)
        {
          merge = traceBack[iTrace];
          mergeSide = curSide[merge] - 1;
          i = nSides[merge, mergeSide] + 1;
          sides[i, merge, mergeSide] = singleton;
          nSides[merge, mergeSide] = i;
          iTrace = iTrace - 1;
        }
        # print("nSides after:"); print(nSides);
        # print("sides[,,1]:"); print(sides[,,1]);
        # print("sides[,,2]:"); print(sides[,,2]);
      } else {
        # Record traceback and do the next merge.
        # print("Branch.")
        x = curSide[curMerge];
        curSide[curMerge] = x+1;
        curMerge = dendro$merge[curMerge, x];
        nTrace = nTrace + 1;
        traceBack[nTrace] = curMerge;
      }
    } else {
      # print("Taking one step back.");
      # Both sides were already processed. Return one step back.
      if (nTrace > 1)
      {
        nTrace = nTrace - 1;
        curMerge = traceBack[nTrace];
      } else {
        printFlush("All done, apparently!");
      }
    }
    # print("Press enter"); scan();
  }

  # Generate all orderings by recursively flipping left and right sides of each merge.

  # It is much easier to flip rankings than orderings.

  nRankings = 2**(min(nMerge, maxMerge)-1);
  ranks = matrix(0, nrow = nSingl, ncol = nRankings);

  cuRanking = 1;  # Ranking being currently worked on

  # Convert the dendorgram order into ranks
  dendrRanking = rep(0, nSingl);
  dendrRanking[dendro$order] = c(1:nSingl);

  # Rankings used as a base before flipping for every merge

  baseRankings = matrix(rep(dendrRanking, nMerge), nrow = nSingl, ncol = nMerge);

  flipState = rep(0, nMerge);

  # The very last merge can be left out.
  merge = nMerge-1;

  if (verbose>1) printFlush(paste(spaces, " ..generating rankings.."));
  counter = 0; countedMerge = nMerge - 5;

  while (merge<nMerge) 
  {
    # print(paste("merge:", merge, "flipState:", paste(flipState, collapse = ", ")));
    if (flipState[merge]==0)
    {
      # First pass through this merge.
      flipState[merge] = 1;
      # Copy base order from above
      if (merge<nMerge) baseRankings[,merge] = baseRankings[,merge+1];
      if ((merge>1) & (merge>nMerge-maxMerge+1))
      {
        # Reset everything below
        flipState[1:(merge-1)] = 0;
        # ...and go one step below
        merge = merge-1;
      } else {
        # No merges below: save current ranking
        ranks[,cuRanking] = baseRankings[, merge];
        cuRanking = cuRanking + 1;
        # Perform flip at the merge
        nLeft = nSides[merge, 1];
        nRight = nSides[merge, 2];
        left = sides[1:nLeft, merge, 1];
        right = sides[1:nRight, merge, 2];
        flippedRanking = baseRankings[, merge];
        #print(paste("Flip: left=", paste(left, collapse = ","), "; right=", 
                     #paste(right, collapse = ", ")));
        #print(paste("Ranks before:", paste( flippedRanking, collapse = ", ")));
        flippedRanking[left] = baseRankings[left, merge] + nRight;
        flippedRanking[right] = baseRankings[right, merge] - nLeft;
        #print(paste("Ranks after", paste( flippedRanking, collapse = ", ")));
        # Save the flipped order
        ranks[, cuRanking] = flippedRanking;
        cuRanking = cuRanking + 1;
        flipState[merge] = 2;
      }
    } else if (flipState[merge]==1)
    {
      if ((verbose > 2) & (merge == countedMerge))
      { 
        printFlush(paste(spaces, "  ..counter on counted merge:", counter, "of",
                         2**(nMerge-countedMerge-1)));
        counter = counter + 1;
      }
      flipState[merge] = 2;
      # perform the flip of this merge
      nLeft = nSides[merge, 1];
      nRight = nSides[merge, 2];
      left = sides[1:nLeft, merge, 1];
      right = sides[1:nRight, merge, 2];
      flippedRanking = baseRankings[, merge];
      #print(paste("Flip: left=", paste(left, collapse = ","), "; right=", 
                   #paste(right, collapse = ", ")));
      #print(paste("Ranks before:", paste( flippedRanking, collapse = ", ")));
      flippedRanking[left] = baseRankings[left, merge] + nRight;
      flippedRanking[right] = baseRankings[right, merge] - nLeft;
      #print(paste("Ranks after", paste( flippedRanking, collapse = ", ")));
      # Reset all flips below
      flipState[1:(merge-1)] = 0;
      # Set new base order
      baseRankings[, merge] = flippedRanking;
      # and go to next lower merge
      merge = merge -1;
    } else {
      # Both flips on this merge were performed and saved... go one step back.
      merge = merge + 1;
    }
  }

  # collectGarbage();

  if (verbose>1) printFlush(paste(spaces, " ..converting rankings into oderings.."));

  
# The following lines do not work for large dendrograms, presumably because subsetting creates a copy of
# the whole object.

  orders = ranks;
  for (i in 1:nRankings) orders[ranks[, i], i] = c(1:nSingl);

# Replace it with slightly more complicated but block code.
# Want to convert everything into a flat vector.

  #BatchSize = as.integer(1000000/nSingl);

  #nBatch = as.integer( (nRankings-1)/BatchSize ) + 1;
  #for (batch in 1:nBatch)

  #nData = nRankings * nSingl;
  #flatRanks = as.integer( c(0:(nData-1)) / nSingl ) * nSingl + as.vector(ranks);
  #orders = flatRanks;
  #orders[flatRanks] = rep(c(1:nSingl), times = nRankings);
  #dim(orders) = c(nSingl, nRankings);

  list(orders = orders, ranks = ranks);
}

evaluateOrderings = function(orders, dstMat)
{
  n = ncol(dstMat);
  weightMat = matrix(1, n, n);
  for (orderDist in 0:(n-1))
  {
    index = cbind(1:(n-orderDist), (1+orderDist):n);
    weightMat[index] = 1/(orderDist+2);
    weightMat[ cbind(index[, 2], index[, 1])] = 1/(orderDist+2);
  }

  weightedDist = apply(orders, 2, function(x) { sum(dstMat[x, x] * weightMat) })
  weightedDist;
}


bestTree = function(dstMat, method = "a", ...)
{
  tree = hclust(as.dist(dstMat), method = method);
  orderings = dendrogramOrderings(tree, ...);
  printFlush("Evaluating orderings...");
  weightedDists = evaluateOrderings(orderings$orders, dstMat);
  best = which.min(weightedDists);
  tree$order = orderings$orders[, best];
  tree;
}



#-----------------------------------------------------------------------------------------------
#
# outlierMeasure
#
#-----------------------------------------------------------------------------------------------

# This function calculates the measure of being an outlier for each sample in each column and averages
# them. It returns a vector of the length = nrow(x) giving the mean measure of outlierishness for each
# sample.   
        
outlierMeasure = function(x)
{
    x = as.matrix(x);
    mx = apply(x, 2, median, na.rm = TRUE)
    mxMat = matrix(mx, nrow = nrow(x), ncol = ncol(x), byrow = TRUE);
    madxMat = matrix(apply(x, 2, mad, na.rm = TRUE), nrow = nrow(x), ncol = ncol(x), byrow = TRUE);
    ux = as.matrix(abs((x - mxMat)/(9 * qnorm(0.75) * madxMat)))
    aval = ifelse(ux <= 1, 1, 0)
    outMeas = ux^2;   # Could multiply by 1-aval as well 
    apply(outMeas, 1, mean, na.rm = TRUE);
}


#=======================================================================================================
#
# Legend with automaitc placement
#
#=======================================================================================================

# preferCorner is not implemented yet.

placeLegend = function(points.x, points.y, nPositions = 20,
                       nPositions.x = nPositions, nPositions.y = nPositions,
                       tryNCol = 1, pointSize = 0.02, 
                       maxOverlap = length(points.x), preferCorner = TRUE,
                       limitBox = par("usr"),
                       preferenceBox = limitBox,
                       ...)
{
  if (is.null(tryNCol)) tryNCol = ncol;
  nTryNCol = length(tryNCol)
  stat = overlap = array(NA, dim = c(nPositions, nPositions, nTryNCol))
  tryPositions.x = tryPositions.y = list();
  legendBox = list();
  points.x = as.numeric(unlist(points.x));
  points.y = as.numeric(unlist(points.y));
  nPoints = length(points.y);
  repX = ceil(nPoints/length(points.x));
  if (repX != nPoints/length(points.x)) warning(immediate. = TRUE,
      "Length of vectorized 'points.y' is not an integer multiple of the length of 'points.x'.");
  points.x = rep(points.x, repX)[1:nPoints];
  xRange = limitBox[2] - limitBox[1];
  yRange = limitBox[4] - limitBox[3];
  for (nc in 1:nTryNCol)
  {
    legendBox[[nc]] = legend(x = "center", ..., ncol = tryNCol[nc], plot = FALSE, xjust = 0, yjust = 0)$rect
    if (xRange < legendBox[[nc]]$w) next;
    if (yRange < legendBox[[nc]]$h) next;
    tryPositions.x[[nc]] = seq(from = limitBox[1], to = limitBox[2] - legendBox[[nc]]$w, length.out = nPositions);
    tryPositions.y[[nc]] = seq(from = limitBox[3], to = limitBox[4] - legendBox[[nc]]$h, length.out = nPositions);
    for (tpx in 1:nPositions) for (tpy in 1:nPositions)
    {
      x.min = tryPositions.x[[nc]] [tpx] - pointSize * xRange;
      x.max = x.min + legendBox[[nc]]$w + 2*pointSize * xRange;
      y.min = tryPositions.y[[nc]] [tpy] - pointSize * yRange;
      y.max = y.min + legendBox[[nc]]$h + 2*pointSize * yRange;
      xdist = (points.x < x.min) * (x.min-points.x) + (points.x > x.max) * (points.x-x.max);
      ydist = (points.y < y.min) * (y.min-points.y) + (points.y > y.max) * (points.y-y.max);
      stat[tpy, tpx, nc] = 
       nPoints * 2*(xRange + yRange)^2 *
         sum(points.x <= x.max & points.x >= x.min & points.y >= y.min & points.y <= y.max, na.rm = TRUE) +
       sum((1/((xRange + yRange)^2) + xdist^2 + ydist^2)^(-1), na.rm = TRUE) + 
       nPoints^2 * 2*(xRange + yRange) * ( 
          (tryPositions.x[[nc]] [tpx] < preferenceBox[1]) * (preferenceBox[1] - tryPositions.x[[nc]] [tpx]) + 
          (tryPositions.x[[nc]] [tpx] + legendBox[[nc]]$w > preferenceBox[2]) * 
                  (tryPositions.x[[nc]] [tpx] + legendBox[[nc]]$w - preferenceBox[2]) + 
          (tryPositions.y[[nc]] [tpy] < preferenceBox[3]) * (preferenceBox[3] - tryPositions.y[[nc]] [tpy]) + 
          (tryPositions.y[[nc]] [tpy] + legendBox[[nc]]$h > preferenceBox[4]) * 
                  (tryPositions.y[[nc]] [tpy] + legendBox[[nc]]$h - preferenceBox[4]));
      overlap[tpy, tpx, nc] = sum(points.x <= x.max & points.x >= x.min & points.y >= y.min &
                                points.y <= y.max, na.rm = TRUE)
    }
  }
  if (all(is.na(stat))) return(list(doPlot = FALSE, x = preferenceBox[1], y = preferenceBox[3], ncol = 1));
  best = which.min(stat);
  best1 = (best-1) %% (nPositions * nPositions);
  nc.best = floor( (best-1) / (nPositions * nPositions))+1;
  ncol = tryNCol[ nc.best];
  x = tryPositions.x[[nc.best]][ floor((best1)/nPositions) + 1];
  y = tryPositions.y[[nc.best]][ (best1)%% nPositions + 1] + legendBox[[nc.best]]$h;
  doPlot = c(overlap)[best]<=maxOverlap

  list(doPlot = doPlot, x = x, y = y, ncol = ncol);
}

pointsFillingLegend = function(legendBox, nx = 40, ny = 20)
{
  x0 = seq(from = legendBox$left, to = legendBox$left + legendBox$w, length.out = nx);
  y0 = seq(from = legendBox$top-legendBox$h, to = legendBox$top, length.out = ny);
  list(x = rep(x0, each = ny), y = rep(y0, nx));
}
  
legendClean = function(x, y=NULL, ..., ncol = 1, points.x = NULL, points.y = NULL, nPositions = 20,
                       tryNCol = 1, pointSize = 0.02, opacity = 0.7, 
                       maxOverlap = length(points.x), preferCorner = TRUE)
{
  doPlot = TRUE;
  if (x=="auto")
  {
    if (is.null(points.x) || is.null(points.y))
      stop("When 'x' is 'auto', arguments 'points.x' and 'points.y' must be given.");
    place = placeLegend(points.x, points.y, nPositions = nPositions, 
              tryNCol = tryNCol, pointSize = pointSize, maxOverlap = maxOverlap, preferCorner = preferCorner,
              ...);
  } else place = list(x = x, y = y, ncol = ncol, doPlot = TRUE);
  box = legend(x = place$x, y = place$y, nc = place$ncol, ..., plot = FALSE)$rect
  if (place$doPlot)
  {
    #text(box$left, box$top-box$h, "Test", cex = 2, adj = c(0,0));
    #rect(box$left, box$top-box$h, box$left +box$w, box$top, col = rgb(1,1,1,opacity), border = "white");
    legend(x = place$x, y = place$y, ..., ncol = place$ncol, bg = rgb(1,1,1,opacity));
  }
  invisible(box);
}

initLib = function(libname, force = FALSE)
{
  os = R.Version()$os;
  recOS = c("linux-gnu", "mingw32");
  ios = match(os, recOS);
  if (is.na(ios)) stop("Unrecognized OS.");
  exts = c(".so", ".dll");
  filename = paste(libname, exts[ios], sep="");
  if (force && is.loaded(filename)) dyn.unload(filename);
  if (!is.loaded(filename)) dyn.load(filename);
}


# Plot enrichment barplot of modules in a yes/no indicator variable

plotEnrichmentBarplot = function(indicator, moduleLabels, moduleLevels, moduleColors,
                xlab = "Module number",
                ylab = "-log(p-value)",
                cex.lab = 1,
                ...)
{
   nGenes = length(moduleLabels);
   geneColors = moduleLabels[indicator];
   presentMods = sort(unique(geneColors))
   nPresMods = length(presentMods)
   pValues = rep(NA, nPresMods)
   nIntGenes = rep(0, nPresMods)
   xGeneColors = rep(-1, nGenes)
   xGeneColors[indicator] = geneColors;
   for (imod in 1:nPresMods)
   {
     mod = presentMods[imod]
     testMatrix = matrix( c(
          sum(xGeneColors>-1 & moduleLabels==mod),
          sum(xGeneColors==-1 & moduleLabels==mod),
          sum(xGeneColors>-1 & moduleLabels!=mod),
          sum(xGeneColors==-1 & moduleLabels!=mod)), 2, 2)
     pValues[imod] = fisher.test(testMatrix, alternative = "greater")$p.value;
     nIntGenes[imod] = testMatrix[1,1];
   }

   pValues[pValues==0] = 1e-300;

   max = max(-log10(pValues), -log10(0.01 / nPresMods) + 1);
   min = 0;

   mp = barplot(-log10(pValues), col = moduleColors[match(presentMods, moduleLevels)], names.arg = FALSE,
                xlab = xlab, ylab = ylab, ylim = c(min, max), ...)
   text(x = mp, y = pmax(-log10(pValues) + 1, rep(-log10(0.01 / nPresMods) + 1, length(pValues))),
        labels =  nIntGenes, adj = c(0.5, 0), xpd = TRUE, cex = 0.8 )
   text(x = mp, y = min - (max-min)*0.05, labels =presentMods, adj = c(0.5, 1), xpd = TRUE, cex = cex.lab )
   abline(h = -log10(0.01 / nPresMods), col = "red")
}
 

# Row-wise consensus: must have equal signs
# i.e. consensus across columns, returning one value for each row (for higher-dimensional arrays, returns
# an array of one dimension lower of consensi across the last dimension)
# attempts to ignore NA's 

rowConsensus = function(data)
{
  nDim = length(dim(data))
  if (nDim < 2) stop("'data' must be a matrix or array with at least 2 dimensions.");
  if (nDim==2) data = as.matrix(data);
  dimd = dim(data);
  nc = dimd[length(dimd)];
  sign = apply(sign(data), c(1:(nDim-1)), sum, na.rm = TRUE)
  sign[is.na(sign)] = 0;
  sign[abs(sign) < nc] = 0;
  consensus = apply(abs(data), c(1:(nDim-1)), min, na.rm = TRUE) * sign/nc;
  consensus;
}

plotHistogramAndCF = function(data, nTicks = NULL, ...)
{
   h = hist(data, ...);
   ymax = par("yaxp")[2];
   n = par("yaxp")[3];
   #if (!(n %in% c(2,4,5,10)))
   if (is.null(nTicks)) nTicks = n;
   sum = sum(h$counts);
   lines(h$mids, cumsum(h$counts) / sum * ymax, col="red");
   axis(side = 4, at = c(0:nTicks)/nTicks * ymax, labels = signif(c(0:nTicks)/nTicks, 2))
   addGrid();
   invisible(h);
}

# The [bi]corAndPvalue functions have been moved to WGCNA package

#===============================================================================================
#
# Plot of a set of variables vs. genomic location.
#
#===============================================================================================

.plotQTLtypes = c("lines", "points");

plotQTL = function(data, chr, bp, colors = 1, type = "lines", bpBasedSpacing = FALSE,
                   lty = 1, lwd = 1, ylim = NULL, xLabOffset = 0.04, 
                   gridColor = "grey", grid.lty = 2, ... )
{

  itype = charmatch(type, .plotQTLtypes);
  if (is.na(itype))
    stop(paste("Unrecognized 'type'. Recognized values are", paste(.plotQTLtypes, collapse = ", ")));

  nchr = as.numeric(chr);
  order = order(nchr, bp);
  data = as.matrix(data);
  nCols = ncol(data);
  nSNPs = length(order)
  x = data[order, , drop = FALSE];
  if (bpBasedSpacing)
  {
    chromosomeNumbers = sort(unique(nchr));
    nChromo = length(chromosomeNumbers);
    chromoLength = tapply(bp, nchr, max);
    dividers = cumsum(chromoLength);
    totalLength = sum(chromoLength);
    snpGlobalPositions = bp;
    addLength = 0; index = 1;
    for (ch in chromosomeNumbers)
    {
      snpGlobalPositions[ch==nchr] = bp[ch==nchr] + addLength;
      addLength = addLength + chromoLength[index];
      index = index + 1;
    }
    plotX = snpGlobalPositions;
    breaksX = c(0.5, dividers + 0.5);
    breaks = dividers[-length(dividers)] + 0.5;
  } else {
    plotX = c(1:nSNPs);
    cx = chr[order];
    cx2 = cx[-1];
    breaks = which(cx[-nSNPs]!=cx2) + 0.5;
    breaksX = c(0.5, breaks, nSNPs + 0.5);
  }
  mids = (breaksX[-1] + breaksX[-length(breaksX)])/2;
   
  if (is.null(colors)) colors = matrix(1, nSNPs, nCols);
  if ((length(colors)==nCols) || (length(colors)==1)) 
     colors = matrix(colors, nSNPs, nCols, byrow = TRUE);
  if (length(lty==1)) lty = rep(lty, nCols);
  if (length(lwd==1)) lwd = rep(lwd, nCols);


  if (is.null(ylim))
  {
     min = min(data, na.rm = TRUE);
     max = max(data, na.rm = TRUE);
     ylim = c(min, max);
  }

  if (itype==1)
  {
    plot(plotX, x[, 1], col = colors[, 1], type = "l", lty = lty[1], lwd = lwd[1], 
         ylim = ylim, xaxt = "none", xlab = "Chromosome", xaxs = "i", ...)

    if (nCols > 1) for (col in 2:nCols)
      lines(plotX, x[, col], col = colors[, col], lty = lty[col], lwd = lwd[col]);
  } else {
    plot(plotX, x[, 1], col = colors[, 1], pch = ".",
         ylim = ylim, xaxt = "none", xlab = "Chromosome", xaxs = "i", ...)

    if (nCols > 1) for (col in 2:nCols)
      points(plotX, x[, col], col = colors[, col], pch = ".");
  }


  box = par("usr");
  for (b in breaks)
    lines(x = rep(b, 2), y = box[3:4], col = gridColor, lty = grid.lty );

  text(mids, rep(box[3] - xLabOffset * (box[4]-box[3]), length(mids)), sort(unique(chr)), xpd = TRUE);

  invisible(plotX);
}

distributionOnSNPs2 = function(dataList, snpChr, snpBp, window, filterType = "gaussian", fun = sum, 
                     normalize = FALSE)
{
  nData = length(dataList);
  nSNPs= length(snpChr);
  distribution = matrix(0, nSNPs, nData);
  snpChrLevs = sort(unique(snpchr));

  ft = charmatch(filterType, .filterTypes);
  if (is.na(ft))
    stop(paste("Unrecognized filter type: ", filterType,
               "\nRecognized values are ", paste(.filterTypes, collapse = ", ")))

  for (m in 1:nData)
  {
    sourceChr = dataList[[m]]$chr;
    sourceBp = dataList[[m]]$bp;
    for (c in snpChrLevs)
    {
      chrSnps = snpChr==c;

      targetBp = snpBp[chrSnps];
      sourceBpC = sourceBp[sourceChr==c];
      finiteBp = is.finite(sourceBpC);
      sourceBpC = sourceBpC[finiteBp];

      # outer(x,y,`-`)[i,j] = x[i]-y[j]
      distMat = outer(sourceBpC, targetBp, `-`);
      dataMat = matrix(dataList[[m]]$data[sourceChr==c][finiteBp], nrow(distMat), ncol(distMat));

      if (ft==1)
      {
         weight = exp(-distMat^2/window^2);
      } else {
         weight = matrix(as.numeric(abs(distMat) <= window), nrow(distMat), ncol(distMat));
      }
      d = apply(weight * dataMat, 2, fun);
      if (normalize)
      {
         sw = apply(weight, 2, fun);
         d = d/sw;
         d[sw==0] = NA;
      }
      distribution[chrSnps, m] = d;
    }
  }
  distribution;
}

#================================================================================================
#
# qtlPeaks
#
#================================================================================================

# Find peaks in data (qtl) that have no near higher QTL (within tol). 
# Returns a logical vector that indicates whether
# each entry is a peak or not.

qtlPeaks = function(qtl, chr, bp, minQTL, window, tol = 0.001)
{
  order = order(chr, bp)
  nSNPs = length(qtl);
  qtlL = c(qtl[-1], 0);
  qtlR = c(0, qtl[-nSNPs]);

  chrL = c(chr[-1], chr[nSNPs]);
  chrR = c(chr[1], chr[-nSNPs]);

  peaks = (qtl >= qtlL | chrL != chr) & (qtl >= qtlR | chr != chrR)

  peaks[is.na(peaks)] = FALSE;
  peaks[peaks][qtl[peaks] < minQTL] = FALSE;

  # remove neighboring peaks that have the exact same LOD
  nPeaks2 = sum(peaks);
  delPeaks2 = rep(FALSE, nPeaks2);
  qtl2 = qtl[peaks];
  delPeaks2[qtl2[-1]==qtl2[-nPeaks2]] = TRUE;
  # This looks a bit confusing :)
  peaks[peaks][delPeaks2] = FALSE;

  peakP = qtl[peaks];
  nPeaks = sum(peaks);
  peakchr = chr[peaks];
  peakBp = bp[peaks];

  # Remove peaks that are not highest in their neighborhood
  delPeaks = rep(FALSE, nPeaks);
  for (p in 1:nPeaks)
  {
    near = chr==peakchr[p] & abs(peakBp[p] - bp) <= window
    maxNearQTL = max(qtl[near], na.rm = TRUE);
   if (maxNearQTL > peakP[p] + tol) delPeaks[p] = TRUE;
  }

  peaks[peaks][delPeaks] = FALSE;

  peaks;
}

#===============================================================================================
#
# labeled boxplot
#
#===============================================================================================

labeledBoxplot = function(x, g = NULL, names = NULL, srt = 45, adj = c(1, 0.5), namesColor = 1, cex.names = 1, 
                          verbose = FALSE, yOffsetForXLabels = 0.9,
                          addScatterplot = FALSE,
                          pt.cex = 0.8, pch = 21, pt.col = "blue", pt.bg = "skyblue",
                          randomSeed = 31425, jitter = NULL, jitterAmount = 0.4,
                          pointLabels = NULL, nLabel = 2,
                          ratio.pointToChar = 0.2, pointLabel.cex = 1, pointLabel.col = 1,
                           ...)
{
  if (!is.null(g))
  {
    x = tapply(x, g, identity);
    if (is.null(names)) names = names(x);
    if (!is.null(pointLabels)) pointLabels = unlist(tapply(pointLabels, g, identity));
  }
  x = as.list(x);
  # Drop all missing
  if (verbose)
  {
    lengths = lapply(x, length);
    n = length(x);
    g = rep(1:n, lengths);
    xx = as.vector(unlist(x));
    out = verboseBoxplot(xx, g, names = FALSE, ...);
  } else 
    out = boxplot(x, names = FALSE, ...);
  n = length(x);
  box = par("usr");
  xText = c(1:n);
  yText = rep(box[3] - yOffsetForXLabels * strheight("M"),  n);
  pointLabels = unlist(pointLabels);
  text(xText, yText, names, srt = srt, adj = adj, col = namesColor, cex = cex.names, xpd = TRUE);
  if (addScatterplot)
  {
    if (exists(".Random.seed"))
    {
      savedSeed = .Random.seed;
      set.seed(randomSeed);
      on.exit(.Random.seed <<- savedSeed);
    }

    nPoints = length(unlist(x))
    pch = WGCNA:::.extend(pch, nPoints);
    pt.col = WGCNA:::.extend(pt.col, nPoints);
    pt.bg = WGCNA:::.extend(pt.bg, nPoints);
    pt.cex = WGCNA:::.extend(pt.cex, nPoints);

    present = lapply(x, function(x1) is.finite(x1))
    x = mymapply(`[`, x, present);
    groupLengths = sapply(x, length);
    present.flat = unlist(present);
    nPoints = length(unlist(x))
    if (!is.null(pointLabels)) pointLabels = pointLabels[present.flat];
    pch = pch[present.flat];
    pt.col = pt.col[present.flat];
    pt.bg = pt.bg[present.flat];
    pt.cex = pt.cex[present.flat];

    set.seed(randomSeed)  # so we can make the identical plot again
    jitterFnc = function(n, amount)
    {
      seq(from=-amount, to=amount, length.out = n)[sample(n)]
    }
    if (is.null(jitter)) jitter = jitterAmount*tanh(out$n/10) else jitter = WGCNA:::.extend(pch,jitter);
    points.x = rep(1:ncol(out$stats),out$n) + unlist(mymapply(jitterFnc, out$n, jitter));
    #points.x = jitter(rep(1:ncol(out$stats),out$n), jitter);

    points.y = unlist(x);
    points(points.x, points.y,
           pch=pch, col=pt.col, bg = pt.bg, cex = pt.cex)

    out = c(out, list(points.x = points.x, points.y = points.y));

    if (!is.null(pointLabels))
    {
      indices = mymapply(function(n, b) b-n + c(1:n), n = groupLengths, b = cumsum(groupLengths));
      x = lapply(x, function(xx) xx[!is.na(xx)]);
      ranks.up = lapply(x, rank);
      ranks.dn = lapply(x, function(x) rank(-x));
      keep.up.0 = mymapply(function(r, i) i[r<=nLabel], ranks.up, indices);
      keep.dn.0 = mymapply(function(r, i) i[r<=nLabel], ranks.dn, indices)
      keep = unlist(c(keep.up.0, keep.dn.0));
      labels.plot = rep("", nPoints);
      labels.plot[keep] = pointLabels[keep];
      labelPoints2(points.x, points.y, labels.plot, ratio.pointToChar = ratio.pointToChar, cex = pointLabel.cex,
                  col = pointLabel.col, pt.cex = pt.cex);
    } 
  }
  invisible(out);
}


#===============================================================================================
#
# labeled plot
#
#===============================================================================================

labeledPlot = function(x, y, names, srt = 45, adj = c(1, 0.5), namesColor = 1, cex.names = 1, ...)
{
  plot(x, y, xaxt = "n", ...);
  n = length(x);
  box = par("usr");
  xText = c(1:n);
  yText = rep(box[3] - 0.02 * (box[4] - box[3]), n);

  text(xText, yText, names, srt = srt, adj = adj, col = namesColor, cex = cex.names, xpd = TRUE);
}


#===============================================================================================
#
# labeled barplot, version 2 (somewhat different aim than the WGCNA function)
#
#===============================================================================================

labeledBarplot2 = function(x, names, g = NULL, horiz = FALSE, 
                           srt = if (horiz) 0 else 45, 
                           adj = if (horiz) c(1, 0.5) else c(1, 0.5), 
                           namesColor = 1, cex.names = 1, 
                           addGuide = TRUE, guideColor = "grey30", guideLty = 2, 
                           labelOffset = 1.5, ...)
{

  if (!is.null(g)) {
    mp = verboseBarplot(x, g, names.arg = NA, horiz = horiz, ...);
    heights = attr(mp, "height");
  } else {
    mp = barplot(x, names.arg = NA, horiz= horiz, ...);
    heights = x;
    attr(mp, "height") = x;
    attr(mp, "stdErr") = rep(0, length(x));
  }

  box = par("usr");

  if (horiz)
  {
    yText = mp;
    xText = rep(box[1] - labelOffset * strwidth("W"), length(mp));
  } else {
    xText = mp;
    yText = rep(box[3] - labelOffset * strheight("M"), length(mp));
  }

  text(xText, yText, names, srt = srt, adj = adj, col = namesColor, cex = cex.names, xpd = TRUE);
  if (addGuide) for (i in 1:length(mp))
    if (horiz)
    {
      lines(c(min(0, heights[i]), box[1] - 0.02 * (box[2] - box[1])),  rep(mp[i], 2), 
            col = guideColor, lty = guideLty);
    } else {
      lines(rep(mp[i], 2), c(min(0, heights[i]), box[3] - 0.02 * (box[4] - box[3])), 
            col = guideColor, lty = guideLty);
    }

  invisible(mp);
  
}


#===============================================================================================
#
# Gaussian smoothing filter
#
#===============================================================================================

smoothGauss = function(x,y, xtest, width, normalize = TRUE)
{
  y[!is.finite(y)] = NA;
  dist = outer(x, xtest, "-");
  weight = exp(-dist^2/(2*width^2));
  yMat = matrix(y, length(x), length(xtest));

  weight[is.na(yMat)] = 0;
  if (normalize)
  {
     val = apply(yMat * weight, 2, sum, na.rm = TRUE)/apply(weight, 2, sum)
  } else
     val = apply(yMat * weight, 2, sum, na.rm = TRUE)
  val;
}

#===============================================================================================
#
# Convert date in the form month/day/year into number of days
# Accepts the separator (usually /, -, or .) and the position of day, month, and year in the date
#
#===============================================================================================

linDate = function(dateStr, separator= "/", day = 2, month = 1, year = 3)
{
  # Here we assume date specified as month/day/year

  daysInMonth = c(31,28,31,30,31,30,31,31,30,31,30,31);
  sumDIM = c(0, cumsum(daysInMonth));
  daysInYear = sum(daysInMonth);

  split = strsplit(dateStr, split = separator, fixed = TRUE);
  nSubj = length(split);
  date = rep(NA, nSubj);
  for (s in 1:nSubj)
  {
    date[s] = as.numeric(split[[s]][day]) + sumDIM[as.numeric(split[[s]][month])] +
              as.numeric(split[[s]][year]) * daysInYear;
  }
  date;
}

#================================================================================================
#
# linearize date
#
#================================================================================================

linearizeDate = function(dates, sep = "-", subtractMin = TRUE, yearPosition = 3, monthPosition = 1)
{
  daysInMonth = c(31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31);
  add = c(0, cumsum(daysInMonth));
  dayPosition = c(1,2,3)[ -c(yearPosition, monthPosition)];

  splitDOB = strsplit(dates, split = sep, fixed = TRUE);
  for (d in 1:length(dates))
     if (length(splitDOB[[d]])!=3) splitDOB[[d]] = c(NA, NA, NA);

  tab = as.matrix(as.data.frame(splitDOB));
  month = suppressWarnings(as.numeric(as.character(tab[monthPosition, ])));

  if (!any(!is.na(month)) && any(is.na(month)))
  {
    monthcode = substring(as.character(tab[monthPosition, ]), 1, 3);
    month = match(tolower(monthcode), 
             c("jan", "feb", "mar", "apr", "may", "jun", "jul", "aug", "sep", "oct", "nov", "dec"));
  }

  lin = as.numeric(as.character(tab[yearPosition, ])) * 365 +
        add[month] + as.numeric(as.character(tab[dayPosition, ]));
  lin - as.numeric(subtractMin) * min(lin, na.rm = TRUE);
}


#===============================================================================================
#
# Convert age in the form xxyxxm to number of years. Assumes all months have the same length.
#
#===============================================================================================

linAge = function(ageStr)
{
   # Assume ageStr contains age in the form xxyxxm where xx are numbers and y and m are year and month
   # delimiters

   split = strsplit(ageStr, split = c("y"), fixed = TRUE)

   nSubj = length(split);
   age = rep(NA, nSubj);
   for (s in 1:nSubj)
   {
     age[s] = as.numeric(split[[s]][1]);
     if (length(split[[s]])==2)
     {
        months = strsplit(split[[s]][[2]], split = "m", fixed = TRUE)[[1]];
        age[s] = age[s] + as.numeric(months)/12;
     }
   }
   age;
}

#=========================================================================================================
#
# Simulation of a causal model. From Old/GeneExpressionSimulation/.
# Note: the causal model simulation has been modified to produce a consistent interpretation of the path
# coefficients.
#
#=========================================================================================================

#===================================================================================================== 
#
# partial correlation of columns of x conditioned on A
#
#=====================================================================================================

matrixPCor = function(x, A, corFnc = "cor", corOptions = "use = 'p'")
{
  cx = eval(parse(text = spaste(corFnc, "(x, ", corOptions, ")")));
  ca = as.numeric(eval(parse(text = spaste(corFnc, "(x, A, ", corOptions, ")"))));

  pc = ( cx - ca %*% t(ca)) /  sqrt( (1-ca^2) %*% t(1-ca^2) )
  pc;
};

zeo = function(x, A, corFnc = "cor", corOptions = "use = 'p'")
{
  n = ncol(x);
  cx = eval(parse(text = spaste(corFnc, "(x, ", corOptions, ")")));
  ca = as.numeric(eval(parse(text = spaste(corFnc, "(x, A, ", corOptions, ")"))));
  caMatCol = matrix(ca, n, n);
  caMatRow = matrix(ca, n, n, byrow = TRUE);
  zeo = (caMatCol - cx * caMatRow) / sqrt( (1-caMatRow^2) * (1-cx^2) );
  diag(zeo) = NA;
  zeo;
}

#===================================================================================================
#
# plotSampled
#
#===================================================================================================

plotSampled = function(x, y = NULL, sample = NULL, col = 1, bg = 0, pch = 1, cex = 1, ...)
{
  n = length(x);

  if (length(col) < length(x)) 
        col = rep(col, ceiling(length(x)/length(col)))
  if (length(pch) < length(x)) 
        pch = rep(pch, ceiling(length(x)/length(pch)))
  if (length(cex) < length(x)) 
        cex = rep(cex, ceiling(length(x)/length(cex)))
  if (length(bg) < length(x)) 
        bg = rep(bg, ceiling(length(x)/length(bg)))

  if (is.null(sample))
  {
     sample = 1:n;
  } else if (length(sample)==1) {
     sample = sample(n, sample, replace = FALSE);
  }

  if (length(sample) < 1) stop("If 'sample' is given, it cannot be empty.");

  plot(x[sample], y = if (is.null(y)) y else y[sample], col= col[sample], bg = bg[sample], pch = pch[sample], 
       cex = cex[sample], ...);
}


#===================================================================================================
#
# Thinned scatterplot
#
#===================================================================================================

thinnedScatterplot = function(x, y, bins = 100, thinTo = 5, log = "", verbose = FALSE, 
                              showProgress = FALSE, showPoints = NULL, debug = FALSE, 
                              stretchLims = 0, stretchLim.x = stretchLims, stretchLim.y = stretchLims,
                              xlim = NULL, ylim = NULL,
                              excludePoints = NULL,
                              ...)
{
  if (showProgress)
    printFlush("  ..flattening data and removing missing values...");

  x = as.numeric(x);
  y = as.numeric(y);

  n = length(x);
  if (n!=length(y)) stop("'x' and 'y' must have the same length.");

  xOriginal = x;
  yOriginal = y;

  if (is.null(showPoints)) showPoints = rep(FALSE, n);

  if (length(showPoints)!=n) 
    stop("If given, 'showPoints' must be a logical vector with same length as 'x'.");

  if (length(excludePoints) > 0 && is.logical(excludePoints)) excludePoints = which(excludePoints);

  if (length(excludePoints) > 0)
  {
    xOriginal.include = xOriginal[-excludePoints];
    yOriginal.include = yOriginal[-excludePoints];
  } else {
    xOriginal.include = xOriginal;
    yOriginal.include = yOriginal;
  }

  keep.base = is.finite(x) & is.finite(y);
  keep = keep.base;
  x = x[keep];
  y = y[keep];
  showPoints = showPoints[keep];

  n = length(x);

  if (length(grep("x", log))>0)
  {
    cutX = log(x)
  } else
    cutX = x;

  if (length(grep("y", log))>0)
  {
    cutY = log(y)
  } else
    cutY = y;

  if (showProgress)
    printFlush("  ..cutting data into bins...");

  xInt = cut(cutX, breaks = bins);
  yInt = cut(cutY, breaks = bins);

  xInt.list = tapply(1:n, xInt, identity);
  yInt.list = tapply(1:n, yInt, identity);

  if (showProgress)
    printFlush("  ..generating counts in bins...");

  t = table(xInt, yInt);
  thin = which(t > thinTo);

  thinX = ((thin - 1) %% bins) + 1;
  thinY = (floor( thin-1)/bins) + 1;

  keep.ind = rep(TRUE, n);
  nThin = length(thin);
  if (showProgress) 
    pind = initProgInd("  ..thinning bins..");
   
  if (nThin > 0) for (th in 1:nThin)
  {
    pointsInBin = intersect(xInt.list[[ thinX[th] ]], yInt.list[[ thinY[th] ]]);
    nBin = length(pointsInBin);
    # printFlush(spaste("nBin: ", nBin, ", thinTo: ", thinTo));
    candidates1 = !(showPoints[pointsInBin]);
    sc1 = sum(candidates1);
    nKeep1 = nBin-sc1;
    if (sc1 > nBin-thinTo+nKeep1)
    {
       keep.ind[sample(pointsInBin[candidates1], nBin-thinTo+nKeep1)] = FALSE;
    } else
       keep.ind[pointsInBin[candidates1]] = FALSE;
    if (showProgress)
       pind = updateProgInd(th/nThin, pind);
  }
  keep[keep] = keep.ind;

  if (length(excludePoints)>0) keep[excludePoints] = FALSE;

  plotPriority = as.numeric(keep);
  plotPriority[keep.base][showPoints] = 2;

  if (showProgress) 
    printFlush("\n ..plotting..");

  if (debug) browser()

  if (is.null(xlim) && stretchLim.x > 0)
    xlim = stretchLim(xOriginal.include, stretchLim.x)
    
  if (is.null(ylim) && stretchLim.y > 0)
    ylim = stretchLim(yOriginal.include, stretchLim.y)

  if (verbose) {
     verboseScatterplot(xOriginal, yOriginal, log = log, sample = which(keep), 
                        plotPriority = plotPriority, xlim = xlim, ylim = ylim, ...) 
  } else 
     plotSampled(xOriginal, yOriginal, log = log, sample = which(keep), 
                 plotPriority = plotPriority, xlim = xlim, ylim = ylim, ...);

  invisible(list(x = xOriginal[keep], y = yOriginal[keep], keptPoints = which(keep)));
}

addScatterplotCounts = function(
  x, y,
  avoid.x = x, avoid.y = y,
  cuts.x, cuts.y,
  addFisherTestPValues = TRUE,
  alternative = "greater",
  edgeOnly = TRUE,
  cornerOnly = addFisherTestPValues,
  printZeroCounts = FALSE,
  prefix = "n=",
  ...,
  nPositions = 20, tryNCol = c(1, 1 + addFisherTestPValues), pointSize = 0.02, maxOverlap = 10,
  y.intersp = 0.2 + 0.6 * addFisherTestPValues,
  opacity = 0.7)
{
  if (any(!is.finite(cuts.x) || !is.finite(cuts.y)))
  {
    warning("addScatterplotCounts: some 'cuts.x' or 'cuts.y' are not finite; skipping adding counts.");
    return(NULL);
  }
  fin = is.finite(x) & is.finite(y);
  x = x[fin];
  y = y[fin];
  #cuts.x = cuts.x[is.finite(cuts.x)];
  #cuts.y = cuts.y[is.finite(cuts.y)];
  bin.x = cut(x, c(-Inf, cuts.x, Inf), labels = FALSE)
  bin.y = cut(y, c(-Inf, cuts.y, Inf), labels = FALSE)

  nBins.x = length(cuts.x) + 1;
  nBins.y = length(cuts.y) + 1;
  tab = WGCNA:::.table2.allLevels(bin.x, bin.y, levels.x = 1:nBins.x, levels.y = 1:nBins.y);

  box = par("usr");# x1, x2, y1, y2

  binLims.x1 = c(box[1], pmax(cuts.x, rep(box[1], length(cuts.x))));  
  binLins.x2 = c(pmin(cuts.x, rep(box[2], length(cuts.x))), box[2]);
  binLims.y1 = c(box[3], pmax(cuts.y, rep(box[3], length(cuts.y))));  
  binLins.y2 = c(pmin(cuts.y, rep(box[4], length(cuts.y))), box[4]);

  for (bx in 1:nBins.x) for (by in 1:nBins.y) 
    if ((tab[bx, by] > 0 || printZeroCounts) && 
        (!edgeOnly || (bx %in% c(1, nBins.x) || (by %in% c(1, nBins.y)))) &&
        (!cornerOnly || (bx %in% c(1, nBins.x) && (by %in% c(1, nBins.y))))) #&&
        #(!cornerOnly || (bx %in% c(1, nBins.x) && (by %in% c(1, nBins.y)))))
  {
    legend = spaste(prefix, tab[bx, by])
    if (addFisherTestPValues && (bx %in% c(1, nBins.x)) && (by %in% c(1, nBins.y)))  # Only add Fisher p-value in corners
    {
      ftp = fisher.test(WGCNA:::.table2.allLevels(bin.x==bx, bin.y==by, levels.x = c(FALSE, TRUE),
                                                 levels.y = c(FALSE, TRUE)),
                        alternative = alternative)$p.value;
      legend = c(legend, if (ftp<1e-200) "p<1e-200" else spaste("p=", signif(ftp, 1)));
    }
    place = try(placeLegend(avoid.x, avoid.y, nPositions = nPositions,
                       tryNCol = tryNCol, pointSize = pointSize,
                       maxOverlap = maxOverlap, preferCorner = TRUE,
                       preferenceBox = c(binLims.x1[bx], binLins.x2[bx],
                                    binLims.y1[by], binLins.y2[by]),
                       limitBox = par("usr"),
                       legend = legend,
                       x.intersp = 0, y.intersp = y.intersp,
                       ...))
     legendClean(place$x, place$y, legend = legend, ncol = place$ncol, 
                       x.intersp = 0, y.intersp = y.intersp, opacity = opacity,
                       ...);
  }
} 

annotatedScatterplot = function(
  x, y, bins = 100, thinTo = 5, log = "", 
  verbose = TRUE,
  showProgress = FALSE, showPoints = NULL, debug = FALSE,

  pointLabels = NULL,
  labelPointsVersion = 2,
  nLabel = if (labelPointsVersion==2) 7 else 5,
  nConsider = if (labelPointsVersion==2) 1000 else 20,
  scaleForSelection = TRUE,
  cex.labels = 0.9,
  forceLabel = NULL,

  cex.lab = 1,
  cex.axis = 1,
  cex.main = 1.2,
  col = 1,
  bg = 0,
  pch = 1,

  pt.cex = 0.5,

  addLegend = TRUE,
  legendPos = "auto",
  leg.pt.cex = 2*pt.cex,
  leg.cex = 0.7*cex.lab,
  pointAttr = NULL, ## Needed for the legend for now... see GNV functions.

  addCounts = TRUE,
  addOverlapPValues = FALSE,
  counts.edgeOnly = TRUE,
  counts.cornerOnly = addOverlapPValues,

  stretchLims = 0.05 * (addLegend | addCounts),
  stretchLim.x = stretchLims, stretchLim.y = stretchLims,

  addAblines = TRUE,
  fdr.x = NULL,
  fdr.y = NULL,
  ab.thresholdType = c("FDR", "p"),
  ab.threshold = if (ab.thresholdType=="p") 0.025 else 0.1,
  addAxisLines = FALSE,
  axisLines.col = "grey",
  excludeIndex = NULL,

  counts.prefix = "n=",
  counts.opacity = 0.7,

  legend = NULL,
  legend.x = "auto", legend.y = NULL,
  legend.args = list(),

  ...

)
{
  x0 = x; y0 = y;
  ab.thresholdType = match.arg(ab.thresholdType);
  scp = thinnedScatterplot(x, y,
                     bins = bins, thinTo = thinTo, log = log, showProgress = showProgress,
                     verbose = verbose,
                     cex.lab = cex.lab, cex.axis = cex.axis, cex.main = cex.main,
                     col = col, bg = bg, pch = pch, cex = pt.cex,
                     showPoints = showPoints,
                     stretchLim.x = stretchLim.x, stretchLim.y = stretchLim.y,
                     debug = debug, excludePoints = excludeIndex,
                     ...);
  if (ab.thresholdType=="p") {
    at.x = at.y = -qnorm(ab.threshold/2);
  } else {
    at.x = ablinePositionFromFDR(Z = x0, fdr = fdr.x, threshold = ab.threshold);
    at.y = ablinePositionFromFDR(Z = y0, fdr = fdr.y, threshold = ab.threshold);
  }

  if (addAblines) addAblines(at.x = at.x, at.y = at.y);
  if (addAxisLines) addAblines(0, col = axisLines.col);
  if (nLabel > 0)
  {
    labInfo = if (labelPointsVersion==2) 
                  labelExtremePoints2(scp$x, scp$y, pointLabels[scp$keptPoints], cex = cex.labels,
                       ratio.pointToChar = 0.33*pt.cex, nLabel = nLabel, nConsider = nConsider, forceLabel = forceLabel,
                       scaleForSelection = scaleForSelection) else
                  labelExtremePoints(scp$x, scp$y, pointLabels[scp$keptPoints], cex = cex.labels,
                       offs = 0.04, nLabel = nLabel, nConsider = nConsider, forceLabel = forceLabel,
                       scaleForSelection = scaleForSelection);
    avoid2 = fillLabelSpaceWithPoints(labInfo$x, labInfo$y, labInfo$label, cex = cex.labels);
  } else avoid2 = list(x = numeric(0), y = numeric(0));
  pfl = list(x = numeric(0), y = numeric(0));
  if (addLegend)
  {
    if (!is.null(pointAttr))
    {
       leg = GNV.addRescueLegend(pointAttr, legendPos, points.x = c(scp$x, avoid2$x), points.y = c(scp$y, avoid2$y),
                           pt.cex = leg.pt.cex, cex = leg.cex);
       pfl = pointsFillingLegend(leg);
    };
    if (!is.null(legend))
    {
       leg = do.call(legendClean, c(list(x = legend.x, y = legend.y, legend = legend,
                         points.x = x, points.y = y), legend.args));
       pfl = pointsFillingLegend(leg);
    }
  } 

  if (addCounts && addAblines)
  {
    addScatterplotCounts(x = x, y = y,
                         prefix = counts.prefix,
                         avoid.x = c(scp$x, pfl$x, avoid2$x), avoid.y = c(scp$y, pfl$y, avoid2$y),
                         cuts.x = c(-at.x, at.x), cuts.y = c(-at.y, at.y),
                         cex = leg.cex, addFisherTestPValues = addOverlapPValues,
                         edgeOnly = counts.edgeOnly, 
                         cornerOnly = counts.cornerOnly,
                         opacity = counts.opacity);
  }
}


#===================================================================================================
#
# verbosePairs
#
#===================================================================================================

verbosePairs = function(x, mar1 = c(2,2,1,1), names = colnames(x), cex.names = 3, sample = NULL,
                        cex.lab = 1, cex.main = 1, cex.axis = 1, breaks = 100, corFnc = "cor",
                        histLog = FALSE, verbose = TRUE, thinned = FALSE, cex = 1, 
                        addGrid.scatterplots = FALSE, 
                        afterScatterplot = NULL, 
                        pointLabels = NULL, label.offs = 0.07, label.cex = 1, ...)
{
  x = as.matrix(x);
  nSets = ncol(x);
  par(mfrow = c(nSets, nSets));
  par(mar = mar1)
  if (is.null(names)) names = c(1:nSets);

  for (s1 in 1:nSets) for (s2 in 1:nSets)
  {
    if (s1==s2)
    {
      #plot(c(0,1), c(0,1), type = "n", xaxt = "n", yaxt = "n");
      #text(0.5, 0.5, adj = c(0.5, 0.5), cex = cex.names, labels = names[s1]);
      if (corFnc=="bicor")
      {
        mean = median(x[, s1], na.rm = TRUE);
        sd = mad(x[, s1], na.rm = TRUE);
        mainLine = spaste("\nmedian: ", signif(mean, 2), ", mad: ", signif(sd, 2));
      } else {
        mean = mean(x[, s1], na.rm = TRUE);
        sd = sd(x[, s1], na.rm = TRUE);
        mainLine = spaste("\nmean: ", signif(mean, 2), ", sd: ", signif(sd, 2));
      }
      par(mar = mar1 + c(0,0,2,0))
      h = hist(x[, s1], breaks = breaks, plot = FALSE)
      if (histLog) h$counts = log10(h$counts + 1);
      plot(h, main = spaste(names[s1], mainLine), xlab = "", cex.lab = cex.lab,
           cex.axis = cex.axis, cex.main = cex.main);
    } else {
      par(mar = mar1)
      if (thinned)
      {
        thinnedScatterplot(x[, s2], x[, s1], cex.lab = cex.lab,
                           cex.axis = cex.axis, cex.main = cex.main, xlab = "", ylab = "",
                           corFnc = corFnc, cex = cex, verbose = verbose, ...);
      } else if (verbose) 
      {
        verboseScatterplot(x[, s2], x[, s1], sample = sample, cex.lab = cex.lab,
                           cex.axis = cex.axis, cex.main = cex.main, xlab = "", ylab = "",
                           corFnc = corFnc, cex = cex, ...);
      } else {
        plotSampled(x[, s2], x[, s1], sample = sample, cex.lab = cex.lab,
                           cex.axis = cex.axis, cex.main = cex.main, xlab = "", ylab = "", 
                           cex = cex, ...);
      }
      if (addGrid.scatterplots)
          addGrid(vert = TRUE, horiz = TRUE, linesPerTick = 2);
      if (!is.null(pointLabels)) labelPoints(x[, s2], x[, s1], labels = pointLabels, cex = label.cex, offs = label.offs);
      eval(afterScatterplot)
    }
  }
}

verbosePairs.matched = function(multiData, varNames = NULL,
                        mar1 = c(2,2,1,1), setNames = names(multiData), cex.names = 3, sample = NULL,
                        cex.lab = 1, cex.main = 1, cex.axis = 1, breaks = 100, corFnc = "cor",
                        histLog = FALSE, verbose = TRUE, thinned = FALSE, cex = 1,
                        addGrid.scatterplots = FALSE,
                        afterScatterplot = NULL, ...)
{
  nSets = length(multiData);
  if (is.null(varNames)) varNames = mtd.apply(multiData, names);
  if (is.null(setNames)) setNames = c(1:nSets);
  par(mfrow = c(nSets, nSets));
  par(mar = mar1)

  for (s1 in 1:nSets) for (s2 in 1:nSets)
  {
    if (s1==s2)
    {
      if (corFnc=="bicor")
      {
        mean = median(multiData[[s1]]$data, na.rm = TRUE);
        sd = mad(multiData[[s1]]$data, na.rm = TRUE);
        mainLine = spaste("\nmedian: ", signif(mean, 2), ", mad: ", signif(sd, 2));
      } else {
        mean = mean(multiData[[s1]]$data, na.rm = TRUE);
        sd = sd(multiData[[s1]]$data, na.rm = TRUE);
        mainLine = spaste("\nmean: ", signif(mean, 2), ", sd: ", signif(sd, 2));
      }
      par(mar = mar1 + c(0,0,2,0))
      h = hist(multiData[[s1]]$data, breaks = breaks, plot = FALSE)
      if (histLog) h$counts = log10(h$counts + 1);
      plot(h, main = spaste(setNames[s1], mainLine), xlab = "", cex.lab = cex.lab,
           cex.axis = cex.axis, cex.main = cex.main);
    } else {
      par(mar = mar1)
      common = intersect(varNames[[s1]]$data, varNames[[s2]]$data);
      nc = length(common)
      x1 = multiData[[s1]]$data[match(common, varNames[[s1]]$data)];
      x2 = multiData[[s2]]$data[match(common, varNames[[s2]]$data)];
     
      if (thinned & nc > 10)
      {
         thinnedScatterplot(x2, x1, cex.lab = cex.lab,
                            cex.axis = cex.axis, cex.main = cex.main, xlab = "", ylab = "",
                            corFnc = corFnc, cex = cex, verbose = verbose, ...);
      } else if (verbose & nc > 2)
      {
        verboseScatterplot(x2, x1, sample = sample, cex.lab = cex.lab,
                           cex.axis = cex.axis, cex.main = cex.main, xlab = "", ylab = "",
                           corFnc = corFnc, cex = cex, ...);
      } else {
        plotSampled(x2, x1, sample = sample, cex.lab = cex.lab,
                           cex.axis = cex.axis, cex.main = cex.main, xlab = "", ylab = "",
                           cex = cex, ...);
      }
      if (addGrid.scatterplots)
          addGrid(vert = TRUE, horiz = TRUE, linesPerTick = 2);
      eval(afterScatterplot)
    }
  }
}


#===================================================================================================
#
# generalPairs
#
#===================================================================================================

histForGeneralPairs = function(x, name, ...)
{
  mean = mean(x, na.rm = TRUE);
  sd = sd(x, na.rm = TRUE);
  mainLine = spaste("\nmedian: ", signif(mean, 2), ", mad: ", signif(sd, 2));

  main = spaste(name, mainLine);
  hist(x, ..., main = main);
}


generalPairs = function(x, mar1 = c(2,2,1,1), names = colnames(x), cex.names = 3, 
                        cex.lab = 1, cex.main = 1, cex.axis = 1, 
                        diagonalFnc = "histForGeneralPairs",  diagonalOptions = list(breaks = 100), 
                        upperFnc, upperOptions, 
                        lowerFnc = upperFnc, lowerOptions = upperOptions)
{
  x = as.matrix(x);
  nSets = ncol(x);
  par(mfrow = c(nSets, nSets));
  par(mar = mar1)
  if (is.null(names)) names = c(1:nSets);

  diagonalFnc = match.fun(diagonalFnc);
  upperFnc = match.fun(upperFnc);
  lowerFnc = match.fun(lowerFnc);
  

  diagonalOptions = c(diagonalOptions, list(cex.lab = cex.lab, cex.main = cex.main, cex.axis = cex.axis));
  upperOptions = c(upperOptions, list(cex.lab = cex.lab, cex.main = cex.main, cex.axis = cex.axis));
  lowerOptions = c(lowerOptions, list(cex.lab = cex.lab, cex.main = cex.main, cex.axis = cex.axis));
  for (s1 in 1:nSets) for (s2 in 1:nSets)
  {
    if (s1==s2)
    {
      par(mar = mar1 + c(0,0,2,0))
      diagonalOptions$x = x[, s1];
      diagonalOptions$name = names[s1];
      do.call(diagonalFnc, diagonalOptions);
      #h = hist(x[, s1], breaks = breaks, plot = FALSE)
      #if (histLog) h$counts = log10(h$counts + 1);
      #plot(h, main = spaste(names[s1], mainLine), xlab = "", cex.lab = cex.lab,
      #     cex.axis = cex.axis, cex.main = cex.main);
    } else {
      par(mar = mar1)
      if (s2>s1) 
      {
        upperOptions1 = c(list(x = x[, s2], y = x[, s1]), upperOptions);
        do.call(upperFnc, upperOptions1);
      } else {
        lowerOptions1 = c(list(x = x[, s2], y = x[, s1]), lowerOptions);
        do.call(lowerFnc, lowerOptions1);
      }
    }
  }
}


#=======================================================================================================
#
# invertEmpiricalDistribution
#
#=======================================================================================================

# Inverts an empirical distribution

invertEmpiricalDistribution = function(data, nBreaks = 100, smoothWindow = 3, plot = FALSE)
{
  if (!is.null(dim(data))) stop("'data' must be a vector.");
  n = length(data);
  h = hist(data, breaks = nBreaks, plot = plot);


  step = h$mids[2] - h$mids[1];
  smooth = smoothGauss(h$mids, h$counts, h$mids, smoothWindow * step);
  smooth[smooth < 0] = 0;
  sums = c(0, cumsum(smooth))/max(cumsum(smooth));

  if (plot) lines(h$mids, smooth, col = "red")
  if (plot) lines(h$breaks, sums * max(h$counts), col = "blue")
  
  y = h$breaks
  x = sums;
  # fit = lm(y~ns(x, df = nBreaks));

  if (FALSE)
  {
    nInvBreaks = nBreaks;  
    invX = seq(from = 0, to=1, length.out = nInvBreaks + 1)
    inverse = predict(fit, newdata = data.frame(x=invX))

    inverse[1] = h$breaks[1];
    inverse[nInvBreaks + 1] = h$breaks[length(h$breaks)];

    if (plot) lines(inverse, invX*max(h$counts), col = "green");

    plot(sums, h$breaks, type = "l", col = 1);
    lines(invX, inverse+1, col = 2);
  }

  #fit;

  list(x = x, y = y);
}

# Make a slightly modified version of the function approx to supress a warning

approx = function (x, y = NULL, xout, method = "linear", n = 50, yleft, 
    yright, rule = 1, f = 0, ties = mean) 
{
    x <- xy.coords(x, y)
    y <- x$y
    x <- x$x
    nx <- length(x)
    method <- pmatch(method, c("linear", "constant"))
    if (is.na(method)) 
        stop("invalid interpolation method")
    stopifnot(is.numeric(rule), (lenR <- length(rule)) >= 1, 
        lenR <= 2)
    if (lenR == 1) 
        rule <- rule[c(1, 1)]
    if (any(na <- is.na(x) | is.na(y))) {
        ok <- !na
        x <- x[ok]
        y <- y[ok]
        nx <- length(x)
    }
    if (!identical(ties, "ordered")) {
        if (length(ux <- unique(x)) < nx) {
#            if (missing(ties)) 
#                warning("collapsing to unique 'x' values")
            y <- as.vector(tapply(y, x, ties))
            x <- sort(ux)
            nx <- length(x)
        }
        else {
            o <- order(x)
            x <- x[o]
            y <- y[o]
        }
    }
    if (nx <= 1) {
        if (method == 1) 
            stop("need at least two non-NA values to interpolate")
        if (nx == 0) 
            stop("zero non-NA points")
    }
    if (missing(yleft)) 
        yleft <- if (rule[1] == 1) 
            NA
        else y[1L]
    if (missing(yright)) 
        yright <- if (rule[2] == 1) 
            NA
        else y[length(y)]
    stopifnot(length(yleft) == 1, length(yright) == 1, length(f) == 
        1)
    if (missing(xout)) {
        if (n <= 0) 
            stop("'approx' requires n >= 1")
        xout <- seq.int(x[1L], x[nx], length.out = n)
    }
    y <- .C("R_approx", as.double(x), as.double(y), as.integer(nx), 
        xout = as.double(xout), as.integer(length(xout)), as.integer(method), 
        as.double(yleft), as.double(yright), as.double(f), NAOK = TRUE, 
        PACKAGE = "stats")$xout
    list(x = xout, y = y)
}

# convenience function

evalInvDist = function(invDist, xnew)
{
  approx(invDist$x, invDist$y, xout = xnew)$y;
}


#===============================================================================================
#
# grey2red
#
#===============================================================================================

grey2red = function(n, base, gamma)
{
  red = seq(from=base^gamma, to=1, length.out = n)^(1/gamma)
  green = blue = seq(from = base^gamma, to=0, length.out = n)^(1/gamma);
  col = rgb(red, green, blue, maxColorValue = 1); 
}


# Example of grey2red:

if (FALSE)
{
  par(mfrow = c(5,1))
  par(mar = c(1,3,1,1))
  n= 100
  barplot(rep(1, n), col = grey2red(n, 0, 1))
  barplot(rep(1, n), col = grey2red(n, 1, 1))
  barplot(rep(1, n), col = grey2red(n, 0.5, 1))
  barplot(rep(1, n), col = grey2red(n, 0.5, 0.2))
  barplot(rep(1, n), col = grey2red(n, 0.5, 5.0))
}



#===============================================================================================
#
# Network plot for generating VisANT-like plots
#
#===============================================================================================

#transform to radians

toRadians = function(angle, degrees = TRUE)
{
  if (degrees) angle = angle * pi/180;
  angle;
}

toDegrees = function(angle, radians = TRUE)
{
  if (radians) angle = angle / pi * 180;
  angle;
}

normalizeAngle = function(angle, degrees = TRUE)
{
  if (degrees)
  {
     step = 360; min = 0;
  } else {
     step = 2*pi; min = 0;
  }

  for (i in 1:length(angle))
  {
    while (angle[i] < min) angle[i] = angle[i] + step;
    while (angle[i] > step + min) angle[i] = angle[i] - step;
  }
  angle;
}

scaleToRange = function(x, min, max, scale = TRUE)
{
  if (scale)
  {
    y = (x-min(x, na.rm = TRUE))/(max(x, na.rm = TRUE) - min(x, na.rm = TRUE)) * (max - min) + min;
  } else {
    y = x;
    y[!is.na(x)] = min;
  }
  y;
}

# Angles here are assumed in degrees
# This should return the distance from center of a symbol to the center of a label
pointToLabelDist = function(labelWidth, labelHeight, pointRadius, alignAngle, labelAngle, degrees = TRUE)
{
  alignAngle = toRadians(alignAngle, degrees);
  labelAngle = toRadians(labelAngle, degrees);

  # Align and label angles cannot go against each other since that would put the label into the point
  # symbol.
  if (cos(alignAngle - labelAngle) < 0) labelAngle = labelAngle + pi;

  al.m.lab = alignAngle - labelAngle
  ratio = labelHeight/labelWidth;
  epsilon = atan(ratio) * sign( sin(al.m.lab) );
  gamma = al.m.lab - epsilon; 

  if (abs(tan(al.m.lab)) > ratio + 2*pointRadius/labelWidth)
  {
     # "Scenario A": long (width) side touches the circle
     return ( (pointRadius + 0.5 * labelHeight)/abs(sin(al.m.lab)) );
  } else if (abs(tan(al.m.lab)) > ratio) {
     # "Scenario B": corner touches the circle
     a = 0.5 * sqrt(labelWidth^2 + labelHeight^2);
     delta.m.al = asin( a*sin(gamma) / pointRadius );
     kappa = pi - gamma - delta.m.al;
     return ( pointRadius * sin(kappa) / sin(gamma) );
  } else {
     # Scenario C: short (height) side touches
     return ( (pointRadius + 0.5 * labelWidth) / abs(cos(gamma)) );
  }
}

labelPosition = function(labelWidth, labelHeight, pointRadius, alignAngle, labelAngle, degrees = TRUE)
{
  dst = pointToLabelDist(labelWidth, labelHeight, pointRadius, alignAngle, labelAngle, degrees);
  alignAngle = toRadians(alignAngle, degrees);
  c( dst * cos(alignAngle), dst * sin(alignAngle) );
}

drawRectangle = function(center.x, center.y, width, height, angle, degrees = TRUE, ...)
{
  corners = matrix( c(-width/2, -height/2, 
                      width/2, -height/2,
                      width/2, height/2,
                      -width/2, height/2,
                      -width/2, -height/2), 2, 5)
  angle = toRadians(angle, degrees);
  rotMat = matrix(c( cos(angle), sin(angle), -sin(angle), cos(angle)), 2, 2);

  points = rotMat %*% corners + matrix(c(center.x, center.y), 2, 5);
  lines(points[1, ], points[2, ], ...);
}

drawEllipse = function(center.x, center.y, width, height = width, angle = 0, startAngle = 0, 
                       stopAngle = 360, degrees = TRUE, nSteps = 100,...)
{
  pixelAngles = toRadians(seq(from=startAngle, to = stopAngle, length.out = nSteps + 1), degrees)
  pixels = matrix( c(width * cos(pixelAngles), 
                     height * sin(pixelAngles)), 2, nSteps + 1, byrow = TRUE);

  angle = toRadians(angle, degrees);
  rotMat = matrix(c( cos(angle), sin(angle), -sin(angle), cos(angle)), 2, 2)
  points = rotMat %*% pixels +
            matrix(c(center.x, center.y), 2, nSteps + 1);
  lines(points[1, ], points[2, ], ...);
}



  

# This function needs x and y for the points.

networkPlot = function(
  adjacency,
  labels,
  pos.x, pos.y, 
  alignAngles,
  labelAngles,
  cex.labels,
  cex.points,
  lineAngles = NULL, # Not used for now
  maxAdj = NULL,
  colors = NULL,
  startNewPlot = TRUE,
  plotBox = c(-1, 1, -1, 1),
  variable.line.width = TRUE,
  min.line.width = 1,
  max.line.width = 5,
  nPlotLines = NULL,
  pch = 21,
  labelColors = "black",
  pointColors = "black",
  pointBg = "black",
  pointLwd = 1,
  xLabelOffset = 0.01,
  yLabelOffset = 0.01,
  degrees = TRUE,
  ...)
{

  if (startNewPlot)
    plot(plotBox[1:2], plotBox[3:4], axes = FALSE, type = "n", xlab = "", ylab = "", ...);

  # plot(c(-1-xMargin,1+xMargin), c(-1-yMargin,1+yMargin), axes = FALSE, type = "n", xlab = "", ylab = "", ...) 
  checkAdjMat(adjacency, min = -1)
  n = length(labels);

  adjx = adjacency
  adjx[is.na(adjx)] = 0;

  if (length(pch)==1) pch = rep(pch, n);
  if (length(labelColors)==1) labelColors = rep(labelColors, n);
  if (length(pointColors)==1) pointColors = rep(pointColors, n);
  if (length(pointBg)==1) pointBg = rep(pointBg, n);
  if (length(xLabelOffset)==1) xLabelOffset = rep(xLabelOffset, n);
  if (length(yLabelOffset)==1) yLabelOffset = rep(yLabelOffset, n);
  if (length(cex.labels)==1) cex.labels = rep(cex.labels, n);
  if (length(cex.points)==1) cex.points = rep(cex.points, n);
  if (length(alignAngles)==1) alignAngles = rep(alignAngles, n);
  if (length(labelAngles)==1) labelAngles = rep(labelAngles, n);

  diag(adjx) = 0;
  diag(adjacency) = 0;
  maxA = max(abs(adjx));
  if (!is.null(maxAdj)) if (maxA<maxAdj) maxA = maxAdj;
  if (sum(adjx < 0) > 0)
  {
     if (is.null(colors)) colors = blueWhiteRed(100);
     adjCol = numbers2colors(adjacency, signed = TRUE, colors = colors, lim = c(-maxA, maxA));
  } else {
     if (is.null(colors)) colors = blueWhiteRed(100)[50:100];
     adjCol = numbers2colors(adjacency, signed = FALSE, colors = colors, lim = c(0, maxA));
  }


  ltA = adjacency;
  diag(ltA) = NA;
  ltA[upper.tri(ltA)] = NA;

  adjOrder = order(c(abs(ltA)))
  rows = row(adjacency)[adjOrder];
  cols = col(adjacency)[adjOrder];

  nLines = n*(n-1)/2;
  if (is.null(nPlotLines)) 
  {
    startLine = 1
  } else { 
    startLine = nLines - nPlotLines + 1;
  }
  for (line in startLine:nLines)
  {
    n1 = rows[line];
    n2 = cols[line];
    a = adjacency[n1, n2];
    normA = abs(a)/maxA;

    w = min.line.width;
    if (variable.line.width)
      w = min.line.width + (max.line.width - min.line.width) * normA;

    #pRadius1 = par("cxy") * cex.points[n1]/35;  # Emprical fudge factor..
    #pRadius2 = par("cxy") * cex.points[n2]/35;
    lineLen = sqrt( (pos.x[n1] - pos.x[n2])^2 + (pos.y[n1] - pos.y[n2])^2);
    x1 = pos.x[n1] #+ pRadius1[1] * (x[n2] - x[n1]) / lineLen
    y1 = pos.y[n1] #+ pRadius1[1] * (y[n2] - y[n1]) / lineLen
    x2 = pos.x[n2] #+ pRadius2[1] * (x[n1] - x[n2]) / lineLen
    y2 = pos.y[n2] #+ pRadius2[1] * (y[n1] - y[n2]) / lineLen

    lines(c(x1,x2),c(y1, y2), lwd = w, col = adjCol[n1, n2]);
  }

  x = pos.x;
  y = pos.y;

  for (node in 1:n)
    points(x[node], y[node], pch = pch[node], cex = cex.points[node], 
           bg = pointBg[node], col = pointColors[node], lwd = pointLwd);

  for (node in 1:n)
  {
    cex = cex.labels[node];
    textWidth = strwidth(labels[node], cex = cex);
    textHeight = strheight(labels[node], cex = cex);
    pRadius = par("cxy") * cex.points[node]/5  ;  # Emprical fudge factor..
    effPointRadius = sqrt(mean(pRadius^2));
    labelShift = labelPosition(textWidth+2*xLabelOffset[node], 
                               textHeight + 2*yLabelOffset[node], effPointRadius, 
                               alignAngles[node], labelAngles[node], degrees);
    #printFlush(paste("Node:", node, " labelShift: ", paste(signif(labelShift, 2), collapse = ", "),
    #               ", width, height: ", paste(signif(c(textWidth, textHeight), 2), collapse = ",")));
    #  Make sure the label will read from left to right
    labAng = toDegrees(normalizeAngle(labelAngles[node], degrees), !degrees);

    if (labAng > 90 & labAng < 270) labAng = labAng - 180;

    text(x[node] + labelShift[1], y[node] + labelShift[2],
         labels = labels[node], adj = c(0.5, 0.5), 
         cex = cex, col = labelColors[node], srt = labAng, xpd = TRUE);
  }

}


# Debug
if (FALSE)
{
labelWidth = textWidth+2*xLabelOffset[node];
labelHeight = textHeight + 2*yLabelOffset[node];
pointRadius = effPointRadius;
alignAngle= alignAngles[node];
labelAngle = labelAngles[node];
}

#===============================================================================================
#
# Circle plot for generating VisANT-like plots
#
#===============================================================================================

circlePlot = function(
  adjacency,
  labels,
  order = NULL,
  maxAdj = NULL,
  colors = NULL,
  startNewPlot = TRUE,
  plotBox = c(-1, 1, -1, 1),
  center = c(0,0), 
  radii = c(0.8, 0.8),
  startAngle = 0,
  variable.cex.labels = TRUE,
  min.cex.labels = 1,
  max.cex.labels = 1.5,
  variable.cex.points = TRUE,
  min.cex.points = 1,
  max.cex.points = 3,
  variable.line.width = TRUE,
  min.line.width = 1,
  max.line.width = 5,
  nPlotLines = NULL,
  pch = 21,
  labelColors = "black",
  pointColors = "black",
  pointBg = "black",
# [x,y]Margin arguments have been taken out. Use radii instead.
  xLabelOffset = 0.01,
  yLabelOffset = 0.01,
  variableLabelAngle = TRUE,
  labelAngleMultiplier = 0.5,
  
  ...)
{

  checkAdjMat(adjacency, min = -1)
  n = length(labels);
  angles = seq(from = startAngle, to = startAngle + 2*pi * (1-1/n), length.out = n);
  x = center[1] + radii[1] * sin(angles);  # This is intentional; top should correspond to angle=0
  y = center[2] + radii[2] * cos(angles);

  adjx = adjacency
  adjx[is.na(adjx)] = 0;
  connectivity = colSums(abs(adjx))-diag(adjx)
  minConn = min(connectivity, na.rm = TRUE);
  maxConn = max(connectivity, na.rm = TRUE);

  if (is.null(order)) order = order(-connectivity, na.last = TRUE)

  if (length(pch)==1) pch = rep(pch, n);
  if (length(labelColors)==1) labelColors = rep(labelColors, n);
  if (length(pointColors)==1) pointColors = rep(pointColors, n);
  if (length(pointBg)==1) pointBg = rep(pointBg, n);
  if (length(xLabelOffset)==1) xLabelOffset = rep(xLabelOffset, n);
  if (length(yLabelOffset)==1) yLabelOffset = rep(yLabelOffset, n);

  oLabs = labels[order]
  oLColors = labelColors[order];
  oPColors = pointColors[order];
  oPBg = pointBg[order];
  oConn = connectivity[order];
  oAdj = adjx[order, order];
  oPch = pch[order];

  alignAngles = normalizeAngle(toDegrees(-angles) + 90); 
  alignAngles[alignAngles > 270] = alignAngles[alignAngles > 270] - 360;
  if (variableLabelAngle)
  {
     labelAngles = ifelse(alignAngles<=90, alignAngles*labelAngleMultiplier, 
                                           180 + (alignAngles - 180)*labelAngleMultiplier)
     alignAngles = labelAngles;
  } else
     labelAngles = 0;

  actualCexPts = scaleToRange(oConn, min.cex.points, max.cex.points, variable.cex.points)
  #actualCexPts = rep(min.cex.points, n);
  #if (variable.cex.points)
  #     actualCexPts = min.cex.points + (max.cex.points - min.cex.points) * 
  #                                (oConn - minConn)/(maxConn - minConn)

  cex.labels = scaleToRange(oConn, min.cex.labels, max.cex.labels, variable.cex.labels);
  #cex.labels = rep(min.cex.labels, n);
  #if (variable.cex.labels)
  #     cex.labels[node] = min.cex.labels + (max.cex.labels - min.cex.labels) *
  #                               (oConn - minConn)/(maxConn - minConn)

  networkPlot(adjacency[order, order],
              oLabs, 
              x, y,
              alignAngles,
              labelAngles,
              cex.labels = cex.labels,
              cex.points = actualCexPts,
              maxAdj = maxAdj,
              colors = colors,
              startNewPlot = startNewPlot,
              plotBox = plotBox,
              variable.line.width = variable.line.width,
              min.line.width = min.line.width,
              max.line.width = max.line.width,
              nPlotLines = nPlotLines,
              pch = oPch,
              labelColors = oLColors,
              pointColors = oPColors,
              pointBg = oPBg,
              xLabelOffset = xLabelOffset[order],
              yLabelOffset = yLabelOffset[order],
              degrees = TRUE, ...);
            
}

# Example of circle plot:

if (FALSE)
{

   sizeGrWindow(8,8)
   par(mfrow = c(1,1));
   nS = 100;
   nn = 30;
   mod = simulateModule(rnorm(nS), nn);
   adjacency = cor(mod)^3;
   
   order = NULL
   
   labels = paste("Gene", c(1:nn));
   circlePlot(adjacency, labels, order, variable.cex.labels = FALSE, radii = c(0.6, 0.6));

   # Plot two circles in one plot

   circlePlot(adjacency, labels, order, variable.cex.labels = FALSE, center = c(-0.5, -0.5), 
              radii = c(0.35, 0.35));

   circlePlot(adjacency, labels, order, startNewPlot = FALSE, 
              variable.cex.labels = FALSE, center = c(0.5, 0.5), 
              radii = c(0.35, 0.35));

   
}


#=======================================================================================================
#
# Hierarchical circle network plot
#
#=======================================================================================================

alternatingAngles = function(n, startAngle, degrees = TRUE)
{
  nm = n-1;
  steps = c(0, rep(c(1:ceil(nm/2)), rep(2, ceil(nm/2))) * rep(c(1, -1), ceil(nm/2)))
  steps = steps[1:n];

  startAngle = toDegrees(startAngle, !degrees);
  angles = normalizeAngle(startAngle + 360/n * steps);
  if (!degrees) angles = toRadians(angles);
  angles 
}


networkPlot.hierarchical = function(
  adjacency,
  labels,
  clusteringAdjacency = adjacency,
  maxAdj = NULL,
  colors = NULL,
  startNewPlot = TRUE,
  plotBox = c(-1, 1, -1, 1),
  center = c(0,0),
  outerRadii = c(0.5, 0.5),

  variableInnerRadius = TRUE,
  minInnerRadii = c(0.1, 0.1),
  maxInnerRadii = c(0.25, 0.25),

  minClusterSize = 3,
  deepSplit = 1,

  angleOfLargestCluster = -90,
  angularSizeOffset = 3,
  variable.cex.labels = TRUE,
  min.cex.labels = 1,
  max.cex.labels = 1.5,
  variable.cex.points = TRUE,
  min.cex.points = 1,
  max.cex.points = 3,
  variable.line.width = TRUE,
  min.line.width = 1,
  max.line.width = 5,
  nPlotLines = NULL,
  pch = 21,
  labelColors = "black",
  pointColors = "black",
  pointBg = "black",
# [x,y]Margin arguments have been taken out. Use radii instead.
  xLabelOffset = 0.01,
  yLabelOffset = 0.01,
  variableLabelAngle = TRUE,
  variableLabelAngleMinSize = 4,
  showSkeleton = FALSE,
  ...)
{
   
  checkAdjMat(adjacency, min = -1)
  n = length(labels);

  adjx = clusteringAdjacency
  adjx[is.na(adjx)] = 0;
  tree = hclust(as.dist(1-adjx), method = "a");

  # Get cluster order around the big circle. 
  # Room for improvement: instead of this could get eigennodes, cluster them and order their dendrogram
  # using dendrogramOrderings 
  clusters.noPam = cutreeDynamic(tree, distM = 1-adjx,
                                 minClusterSize = minClusterSize, 
                                 deepSplit = deepSplit,
                                 pamStage = FALSE);
  clusterOrder = unique(clusters.noPam[tree$order]);
  clusterOrder = clusterOrder[clusterOrder!=0];


  # Here I want all points in a cluster
  clusters = cutreeDynamic(tree, distM = 1-adjx,
                           minClusterSize = minClusterSize,
                           deepSplit = deepSplit,
                           pamStage = TRUE, cutHeight = 2);

  #plotDendroAndColors(tree, cbind(clusters.noPam, clusters))
  clusterSizes = table(clusters)
  nClusters = length(clusterSizes);

  if (nClusters==1)
  {
    circlePlot(
        adjacency,
        labels,
        order = NULL,
        maxAdj = maxAdj,
        colors = colors,
        startNewPlot = startNewPlot,
        plotBox = plotBox,
        center = center,
        radii = radii,
        startAngle = startAngle,
        variable.cex.labels = variable.cex.labels,
        min.cex.labels = min.cex.labels,
        max.cex.labels = max.cex.labels,
        variable.cex.points = variable.cex.points,
        min.cex.points = min.cex.points,
        max.cex.points = max.cex.points,
        variable.line.width = variable.line.width,
        min.line.width = min.line.width,
        max.line.width = max.line.width,
        nPlotLines = nPlotLines,
        pch = pch,
        labelColors = labelColors,
        pointColors = pointColors,
        pointBg = pointBg,
        xLabelOffset = xLabelOffset,
        yLabelOffset = yLabelOffset,
        variableLabelAngle = variableLabelAngle,
        showSkeleton = FALSE,
        ...);
    return ;
  }

  effectiveSizes = clusterSizes + angularSizeOffset 
  if (variableInnerRadius & (max(clusterSizes) > min(clusterSizes))) 
  {
    ratios = effectiveSizes / max(effectiveSizes);
    if (min(ratios) < max(minInnerRadii/maxInnerRadii))
       ratios = scaleToRange(ratios, max(minInnerRadii/maxInnerRadii), 1);
    innerRadii = matrix( ratios, nClusters, 2)  * 
                 matrix( maxInnerRadii, nClusters, 2, byrow = TRUE);
  } else
    innerRadii = matrix( minInnerRadii, nClusters, 2, byrow = TRUE);

  # This may not work well in certain cases but for now will do...
  meanApproxInnerRadius = sqrt(rowMeans(innerRadii^2));

  clusterSlot = match(names(clusterSizes), clusterOrder);

  angleRanges = meanApproxInnerRadius/sum(meanApproxInnerRadius) * 360;
  angleRanges.ordered = angleRanges[ clusterOrder ]
  circleAngles1 = c(0, cumsum(angleRanges.ordered)) + c(angleRanges.ordered/2, 0);
  circleAngles = circleAngles1[1:nClusters] - circleAngles1[2]/2;
  
  largestCluster = which.max(clusterSizes);
  circleAngles = circleAngles - circleAngles[ clusterSlot[largestCluster] ] + angleOfLargestCluster;

  if (showSkeleton) 
  {
     addGrid(v = TRUE)
     abline(h=center[2]);
     abline(v=center[1])
  }

  n = ncol(adjacency);
  x = rep(NA, n);
  y = rep(NA, n);
  kIn = rep(0, n);
  alignAngles = rep(0, n);

  if (showSkeleton) drawEllipse(center[1], center[2], outerRadii[1], outerRadii[2]);
  clusterLevels = names(clusterSizes);
  for (c in 1:nClusters)
  {
    inCluster = clusters==clusterLevels[c];
    kIn[inCluster] = colSums(abs(adjacency[inCluster, inCluster]));
    clusterAngle = circleAngles[ clusterSlot[c] ];
    incOrder = rank(-kIn[inCluster]);
    incAngles = alternatingAngles( sum(inCluster), 180 +clusterAngle)[incOrder];
    clustCenter = center + outerRadii* c(cos(toRadians(clusterAngle)), sin(toRadians(clusterAngle)));
    x[inCluster] = clustCenter[1] + innerRadii[c, 1] * cos(toRadians(incAngles));  
    y[inCluster] = clustCenter[2] + innerRadii[c, 2] * sin(toRadians(incAngles));
    alignAngles[inCluster] = normalizeAngle(incAngles);
    if (showSkeleton) 
      {
        drawEllipse(clustCenter[1], clustCenter[2], innerRadii[c, 1], innerRadii[c, 2]);
        points(clustCenter[1], clustCenter[2], pch = 5)
      }
  }

  cex.points = scaleToRange(kIn, min.cex.points, max.cex.points, variable.cex.points)
  cex.labels = scaleToRange(kIn, min.cex.labels, max.cex.labels, variable.cex.labels);

  alignAngles[alignAngles > 270] = alignAngles[alignAngles > 270] - 360;
  labelAngles = rep(0, n);
  for (c in 1:nClusters)
  {
    inCluster = clusters==clusterLevels[c];
    if (variableLabelAngle & sum(inCluster) >= variableLabelAngleMinSize)
    {
       alignAngles[inCluster] = ifelse(alignAngles[inCluster]<=90, alignAngles[inCluster]/2, 
                                          180 + (alignAngles[inCluster] - 180)/2)
       labelAngles[inCluster] = alignAngles[inCluster];
    } 
  }

  networkPlot(
        adjacency,
        labels,
        pos.x = x, pos.y = y,
        alignAngles = alignAngles,
        labelAngles = labelAngles,
        cex.labels = cex.labels,
        cex.points = cex.points,
        lineAngles = NULL, # Not used for now
        maxAdj = maxAdj,
        colors = colors,
        startNewPlot = startNewPlot & !showSkeleton,
        plotBox = plotBox,
        variable.line.width = variable.line.width,
        min.line.width = min.line.width,
        max.line.width = max.line.width,
        nPlotLines = nPlotLines,
        pch = pch,
        labelColors = labelColors,
        pointColors = pointColors,
        pointBg = pointBg,
        xLabelOffset = xLabelOffset,
        yLabelOffset = yLabelOffset,
        degrees = TRUE,
        ...);
  invisible(list(x = x, y = y, clusters = clusters, tree =tree));
}

#=======================================================================================================
#
# Force-directed network layout, my own version
#
#=======================================================================================================

# strong adjacencies: strong atractive force until a short distance away
# moderate and weak adjacencies: progressively weakening repulsion until middle distance, then no force.

forceFunction = function(
  x,y,
  similarity, 
  areaRadius,
  minDistance,

  radialForceFactor = 10,
  contactForceFactor = 10,

  potential.a = areaRadius/2, 
  potential.min = 0.2,
  potential.B = 1,
  potential.beta = 2
)
{
  nPoints = length(x);
  force.x = force.y = rep(0, nPoints)
  radii = sqrt(x^2 + y^2);

  # Forces that keep points from running away
  radialFactor = pmax( rep(0, nPoints), radii - areaRadius);
  radialForce = radialForceFactor * radialFactor;
  radialForceVec = -radialForce * cbind(x,y)/radii;
  radialForceVec[!is.finite(radialForceVec)] = 0;

  # inter-point distances
  relativeX = outer(x, x, `-`);
  relativeY = outer(y, y, `-`);
  relativeDist = sqrt(relativeX^2 + relativeY^2);

  zeroMat = matrix(0, nPoints, nPoints)

  # Interpoint forces

  pairForce.contact = contactForceFactor * pmax(zeroMat, sign(minDistance - relativeDist));

  if.a = minDistance + (areaRadius - minDistance) * (1-similarity);
  V0 = potential.min + potential.B * similarity^potential.beta
  pairForce.pairs = -V0 * (relativeDist-if.a)/areaRadius;

  pairForce.norm = pairForce.contact + pairForce.pairs;
  diag(pairForce.norm) = 0;
  pairForce = sapply(list(x = relativeX, y = relativeY), function(delta) 
  {
     out = pairForce.norm * delta/relativeDist;
     diag(out) = 0;
     rowSums(out);
  })
  forceVec = pairForce + radialForceVec;
  forceNorm = sqrt(forceVec[, 1]^2 + forceVec[, 2]^2)
  list(forceVec = forceVec, forceNorm = forceNorm, dist = relativeDist);
}

networkLayout.forceDirected = function(
  similarity, radius = 1, 
  minDistance = radius/30,
  x = NULL, y = NULL,
  randomSeed = 123445,
  radialForceFactor = 10,
  contactForceFactor = 10,

  potential.a = radius/2,
  potential.min = 0.2,
  potential.B = 1,
  potential.beta = 2,

  timeStepFactor = 2,
  plotProgress = TRUE,
  maxIterations = 100,
  sleepInterval = 0.03,

  point.col = 1,
  point.bg = NULL,
  point.cex = 1,
  point.lwd = 1
  )
{
  checkAdjMat(similarity, min = 0, max = Inf);
  nPoints = ncol(similarity);

  if (plotProgress)
  {
    x11();
    par(mar = c(0,0,0,0));
    plot(0,0, type = "n", xlim = c(-radius, radius), ylim = c(-radius, radius), axes = FALSE, frame = FALSE);
  }

  if (is.null(x))
  {
    if (exists(".Random.seed"))
    {
      savedSeed = .Random.seed;
      set.seed(randomSeed)
      on.exit(.Random.seed <<- savedSeed);
    }
    x = runif(nPoints, min = -radius, max = radius);
    y = runif(nPoints, min = -radius, max = radius);
  }

  if (plotProgress)
  {
    box = par("usr")
    if (is.null(point.bg)) point.bg = standardColors()[1:nPoints]
    points(x, y, pch = 21, bg = point.bg);
  }

  it = 1;
  while (it <= maxIterations)
  {
    forceData = forceFunction(x, y, similarity,
       areaRadius = radius,
       minDistance = minDistance,
       radialForceFactor = radialForceFactor,
       contactForceFactor = contactForceFactor,
       potential.a = potential.a,
       potential.min = potential.min,
       potential.B = potential.B,
       potential.beta = potential.beta);

    timeStep = minDistance/forceData$forceNorm/timeStepFactor;
    timeStep[!is.finite(timeStep)] = mean(timeStep[is.finite(timeStep)]);
    shift = forceData$forceVec * timeStep;
    if (plotProgress) rect(box[1], box[3], box[2], box[4], col = "white", border = 0);
    x = x + shift[, 1];
    y = y + shift[, 2];
    if (plotProgress) points(x, y, pch = 21, bg = point.bg, cex = point.cex, lwd = point.lwd, col = point.col);
    it = it + 1;
    if (plotProgress) Sys.sleep(sleepInterval);
  }

  list(x = x, y = y);
}

if (FALSE)
{
  # test the force-directed network layout

  n = 3;
  data = matrix(rnorm(10*n), 10, n);
  similarity = abs(cor(data));
  similarity[similarity<0] = 0;

  similarity[1, 2] = similarity[2, 1] = 0.9

  radius = 1
  minDistance = radius/30;
  x = c(0, -1, 1.01); y = c(0, -1, 0.99);
  randomSeed = 123445;
  radialForceFactor = 10;
  contactForceFactor = 10;

  potential.a = radius/2;
  potential.min = 0.2;
  potential.B = 1;
  potential.beta = 2;

  timeStepFactor = 2;
  plotProgress = TRUE
  maxIterations = 100


  nModules = 5;
  nSamples = 50;
  nGenes = 200;
  mat = matrix(0, nModules, nModules)
  mat[1, 2] = 0.4;
  mat[1, 3] = 0.6;
  mat[4, 5] = 0.7;
  eigenNet = simulateEigengeneNetwork(mat, anchorIndex = c(2,3,5), 
          anchorVectors = matrix(rnorm(nSamples * 3), nSamples, 3));
  data = simulateDatExpr(eigenNet$eigengenes, nGenes, modProportions = c(0.3, 0.2, 0.12, 0.1, 0.08, 0.1), 
                          signed = TRUE);

  sim = adjacency(data$datExpr, type = "signed hybrid", power = 2);

  layout = networkLayout.forceDirected(
     similarity = tanh(5*sim), radius = 1,
     minDistance = radius/30,
     x = NULL, y = NULL,
     randomSeed = 123445,
     radialForceFactor = 10,
     contactForceFactor = 10,

     potential.a = radius/2,
     potential.min = 0.2,
     potential.B = 1,
     potential.beta = 2,

     timeStepFactor = 5,
     plotProgress = TRUE,
     maxIterations = 500,
     sleepInterval = 0.010,
     point.bg = labels2colors(data$allLabels)
     )

  plot(layout$x, layout$y, pch = 21, bg = labels2colors(data$allLabels));

  
}

    
#=======================================================================================================
#
# Force-directed network layout, second version
#
#=======================================================================================================

# Adjacencies are converted to "optimal" distances; forces point to the optimal distance for each pair of points.
# In addition, there is a contact term  that prevents points from overlapping, and a boundary repulsion term.
# These may actually be unnecessary but put them in anyway.
# The area is assumed to be rectangular.

# Linear "negative" scale: scale the negative of x to distances
linearNegativeScale = function(x, minx = min(x, na.rm = TRUE), maxx = max(x, na.rm = TRUE), minDist = 0, maxDist, ...)
{
  (maxx-x)/(maxx-minx) * (maxDist-minDist) + minDist;
}


forceFunction2.defaultArgs = function()
{
  list(minDist.pairs = 0.02,
  minDist.boundary = 0.02,

  minDist.coeff = 10,
  boundaryDist.coeff = 10,

  similarityToDesiredDistance.fnc = linearNegativeScale,
  similarityToDesiredDistance.args = list(),

  pairwise.coeff = 1,
  pairwise.power = 1);
}

forceFunction2 = function(
  x,y,
  similarity, 
  width = 1,
  height = 1,
  minDist.pairs = 0.02,
  minDist.boundary = 0.02,

  minDist.coeff = 10,
  boundaryDist.coeff = 10,

  similarityToDesiredDistance.fnc = linearNegativeScale,
  similarityToDesiredDistance.args = list(),

  pairwise.coeff = 1,
  pairwise.power = 1
  
)
{
  nPoints = length(x);
  force.x = force.y = rep(0, nPoints)
  zeroVec = rep(0, nPoints)

  Rx = width/2;
  Ry = height/2;

  # For x and y outside of the rectangle, set the boundary distance to something small so they get pushed in
  # Here I assume that the points are not far outside of the boundary.
  boundaryDist.x = abs(Rx-x);
  boundaryDist.x[ abs(x) > Rx ] = minDist.boundary/5;
  boundaryDist.y = abs(Ry-y);
  boundaryDist.y[ abs(y) > Ry ] = minDist.boundary/5;

  # Absolute values of the boundary contact force
  boundaryFx.0 = boundaryDist.coeff * pmax(zeroVec, 1/boundaryDist.x - 1/minDist.boundary);
  boundaryFy.0 = boundaryDist.coeff * pmax(zeroVec, 1/boundaryDist.y - 1/minDist.boundary);
  # Add the appropriate sign to always point towards the center
  boundaryFx = -sign(x) * boundaryFx.0;
  boundaryFy = -sign(y) * boundaryFy.0;

  # inter-point distances
  relativeX = outer(x, x, `-`);
  relativeY = outer(y, y, `-`);
  relativeDist = sqrt(relativeX^2 + relativeY^2);
  diag(relativeDist) = 1;  ## This ensures that self-forces will be zero

  # Regularized inter-point distances. If points are too close, set the distance to a non-zero constant.
  relativeDist.trunc = relativeDist;
  relativeDist.trunc[relativeDist < minDist.pairs/5] = minDist.pairs/5;

  zeroMat = matrix(0, nPoints, nPoints)

  # Interpoint contact force magnitude
  pairF.contact.0 = minDist.coeff * pmax(zeroMat, 1/relativeDist.trunc - 1/minDist.pairs);
  # interpoint contact force components
  pairFx.contact = rowSums(pairF.contact.0 * relativeX/relativeDist);
  pairFy.contact = rowSums(pairF.contact.0 * relativeY/relativeDist);
  
  # interpoint similarity forces
  # optimal pairwise distances from similarities 
  optimalIPDist = do.call(similarityToDesiredDistance.fnc,
     c(list(x = similarity, minDist = 0, maxDist = sqrt(width^2 + height^2)),
       similarityToDesiredDistance.args));

  diag(optimalIPDist) = diag(relativeDist);

  pairFx.sim = rowSums(sign(optimalIPDist-relativeDist) * pairwise.coeff * abs(optimalIPDist-relativeDist)^pairwise.power  *
                     relativeX/relativeDist)
  
  pairFy.sim = rowSums(sign(optimalIPDist-relativeDist) * pairwise.coeff * abs(optimalIPDist-relativeDist)^pairwise.power  *
                     relativeY/relativeDist);

  forceVec = cbind(x = boundaryFx + pairFx.contact + pairFx.sim, y = boundaryFy + pairFy.contact + pairFy.sim);
  forceNorm = sqrt(forceVec[, 1]^2 + forceVec[, 2]^2)
  list(forceVec = forceVec, forceNorm = forceNorm, dist = relativeDist);
}



networkLayout.functionDirected = function(
  similarity, width = 1, height = 1,
  x = NULL, y = NULL,
  randomSeed = 123445,

  forceFnc = forceFunction2,
  forceArgs = list(),

  minDistance = width/50,
  timeStepFactor = 2,
  noise = minDistance/4/timeStepFactor,
  plotProgress = TRUE,
  maxIterations = 100,
  sleepInterval = 0.03,

  point.col = 1,
  point.bg = NULL,
  point.cex = 1,
  point.lwd = 1
  )
{
  checkAdjMat(similarity, min = -Inf, max = Inf);
  nPoints = ncol(similarity);

  if (plotProgress)
  {
    x11();
    par(mar = c(0,0,0,0));
    plot(0,0, type = "n", xlim = c(-width/2, width/2), ylim = c(-height/2, height/2), axes = FALSE, frame = FALSE);
  }

  if (exists(".Random.seed"))
  {
    savedSeed = .Random.seed;
    set.seed(randomSeed)
    on.exit(.Random.seed <<- savedSeed);
  }

  if (is.null(x))
  {
    x = runif(nPoints, min = -width/2, max = width/2);
    y = runif(nPoints, min = -height/2, max = height/2);
  }

  if (plotProgress)
  {
    box = par("usr")
    if (is.null(point.bg)) point.bg = standardColors()[1:nPoints]
    points(x, y, pch = 21, bg = point.bg, lwd = point.lwd, cex = point.cex, col = point.col);
  }
  it = 1;
  while (it <= maxIterations)
  {
    forceData = do.call(forceFnc, c(list(x = x, y = y, similarity = similarity, width = width, height = height),
        forceArgs));

    timeStep = minDistance/forceData$forceNorm/timeStepFactor;
    timeStep[!is.finite(timeStep)] = mean(timeStep[is.finite(timeStep)]);
    shift = forceData$forceVec * timeStep + noise * matrix(runif(2*nPoints, min = -1, max = 1), nPoints, 2);
    if (plotProgress) rect(box[1], box[3], box[2], box[4], col = "white", border = 0);
    x = x + shift[, 1];
    y = y + shift[, 2];
    if (plotProgress) points(x, y, pch = 21, bg = point.bg, cex = point.cex, lwd = point.lwd, col = point.col);
    it = it + 1;
    if (plotProgress) Sys.sleep(sleepInterval);
  }

  list(x = x, y = y);
}

if (FALSE)
{
  # Test the layout
 n = 3;
  data = matrix(rnorm(10*n), 10, n);
  similarity = abs(cor(data));
  similarity[similarity<0] = 0;

  similarity[1, 2] = similarity[2, 1] = 0.9

  radius = 1
  minDistance = radius/30;
  x = c(0, -1, 1.01); y = c(0, -1, 0.99);
  randomSeed = 123445;

  diag(similarity) = 0;
  layout = networkLayout.functionDirected(
     similarity = similarity,
     x = NULL, y = NULL,
     randomSeed = 123445,

     timeStepFactor = 5,
     plotProgress = TRUE,
     maxIterations = 500,
     sleepInterval = 0.010,
     #point.bg = labels2colors(data$allLabels)
     )

  set.seed(2)
  nModules = 5;
  nSamples = 50;
  nGenes = 200;
  mat = matrix(0, nModules, nModules)
  mat[1, 2] = 0.4;
  mat[1, 3] = 0.6;
  mat[4, 5] = 0;
  eigenNet = simulateEigengeneNetwork(mat, anchorIndex = c(2,3,5),
          anchorVectors = matrix(rnorm(nSamples * 3), nSamples, 3));
  data = simulateDatExpr(eigenNet$eigengenes, nGenes, modProportions = c(0.3, 0.2, 0.12, 0.1, 0.08, 0.1),
                          signed = TRUE, minCor = 0.5);

  sim = cor(data$datExpr)

  diag(sim) = 0;
  args = forceFunction2.defaultArgs();
  args$similarityToDesiredDistance.args$minx = -1.5;

  layout = networkLayout.functionDirected(
     similarity = sim,
     x = NULL, y = NULL,
     randomSeed = 12,

     forceArgs = args,
     timeStepFactor = 3,
     plotProgress = TRUE,
     noise = 0.005,
     maxIterations = 500,
     sleepInterval = 0.0050,
     point.bg = labels2colors(data$allLabels)
     )

  plot(layout$x, layout$y, pch = 21, bg = labels2colors(data$allLabels));


}


   

#=======================================================================================================
#
# Interval overlap
#
#=======================================================================================================

# Interval overlap: 
#    return -1 if [xl, xh] does not intersect with [yl, yh]
#    return 1 if there's partial overlap between the intervals
#    return 2 if there's partial overlap and the ymid is in the [xl, xh] interval
#    return 3 if [yl, yh] is a subset of [xl, xh]
#    return 4 if [xl, xh] is a subset of [yl, yh]
#    return 5 if [xl, xh] is a subset of [yl, yh] and ymid is in [xl, xh]

intervalOverlap = function(xl, xh, yl, yh, ymid = NULL)
{
  n = length(xl);
  m = length(yl);
  result = matrix(-1, m, n);
  xl = as.double(xl);
  xh = as.double(xh);
  yl = as.double(yl);
  yh = as.double(yh);

  a = outer(yl, xl, `-`) * outer(yl, xh, `-`) <= 0; 
  if (!is.null(ymid))
  {
    ymid = as.double(ymid);
    b = outer(ymid, xl, `-`) * outer(ymid, xh, `-`) <= 0;
  } else 
    b = array(FALSE, dim = dim(a));

  c = outer(yh, xl, `-`) * outer(yh, xh, `-`) <= 0;
  d = outer(yh, xl, `-`) * outer(yl, xl, `-`) <= 0;
  e = outer(yh, xh, `-`) * outer(yl, xh, `-`) <= 0;

  result[xor(a,c)] = 1 + as.numeric(b[xor(a,c)]);
  result[a&c] = 3;
  result[d&e] = 4 + as.numeric(b[d&e]);
  t(result);
}

#=======================================================================================================
#
# eps - produce eps plots
#
#=======================================================================================================

eps = function(...) { postscript(..., horizontal = FALSE, onefile = FALSE, paper = "special"); }

#========================================================================================================
#
# Array of vertical barplots
#
#=======================================================================================================
    
barplotArray = function(mat, 
                        lim = NULL,
                        commonScale = FALSE, 
                        colors = "grey40",
                        border = "black",
                        barGap = 0,
                        ablines = NULL,
                        abColors = NULL,
                        ablines.lty = NULL,
                        textMat = NULL,
                        cex.text = 0.8,
                        leftThreshold = 0.8,
                        alignRight = TRUE,
                        col.text = 1,
                        textMargin = 0.04/ncol(mat),
                        margin = 0.1, 
                        textMat.minXPosition = 0,
                        main = "",
                        xlab = "", ylab = "",
                        cex.main = 1.4
                        
                        )
{
  nCol = ncol(mat);
  nRow = nrow(mat);

  if (length(colors)==1 | length(colors)==nRow)
  {
    colors = matrix(colors, nRow, nCol)
  } else if (length(colors)==nCol) {
    colors = matrix(colors, nRow, nCol, byrow = TRUE)
  } else if (!isTRUE(all.equal(dim(colors), c(nRow, nCol))))
    stop(spaste("'colors' must be either a single color, a vector of length nrow or ncol(mat),", 
                "or a matrix of the same dimensions as mat."));

  if (length(border)==1 | length(border)==nRow)
  {
    border = matrix(border, nRow, nCol)
  } else if (length(border)==nCol) {
    border = matrix(border, nRow, nCol, byrow = TRUE)
  } else if (!isTRUE(all.equal(dim(border), c(nRow, nCol))))
    stop(spaste("'border' must be either a single color, a vector of length nrow or ncol(mat),", 
                "or a matrix of the same dimensions as mat."));

  plot(c(0,1), c(0,1), xlim = c(0,1), ylim = c(0,1), axes = FALSE, main = main, cex.main = cex.main,
       type = "n", xlab = xlab, ylab = ylab);

  box = par("usr");
  xMin = box[1];
  xMax = box[2];
  yMin = box[3];
  yMax = box[4];

  if (is.null(lim))
  {
    if (commonScale)
    {
      matMax = rep(max(mat, na.rm = TRUE), nCol);
      matMin = rep(min(mat, na.rm = TRUE), nCol);
    } else {
      matMax = apply(mat, 2, max, na.rm = TRUE);
      matMin = apply(mat, 2, min, na.rm = TRUE);
    }
  } else { 
    if (length(lim)==2) lim = matrix(lim, 2, nCol);
    if (!all.equal(dim(lim), c(2, nCol)))
      stop(paste("if 'lim' is given, it must be either a vector of length 2 or a matrix of 2 rows and",
                 "same number of columns as 'mat'"));
    matMax = lim[2, ];
    matMin = lim[1, ];
  }

  matMin[matMin > 0] = 0; 
  mat[is.na(mat)] = 0;

  x0 = matrix(c(0:(nCol-1))/nCol * (xMax - xMin) + xMin, nRow,nCol, byrow = TRUE);
  x1 = matrix(c(1:nCol)/nCol * (xMax - xMin) + xMin, nRow,nCol, byrow = TRUE);

  y0 = matrix(c((nRow-1):0)/nRow * (yMax - yMin)+ yMin, nRow,nCol);
  y1 = matrix(c(nRow:1)/nRow * (yMax - yMin) + yMin, nRow,nCol);

  y1 = y1 - (y1-y0) * barGap;

  if (!is.null(ablines))
  {
    if (is.null(abColors)) abColors = rep(1, length(ablines))
    if (is.null(ablines.lty)) ablines.lty = rep(1, length(ablines))
  }
  
  out = list(xleft = x0, xright = x1, ytop = y0, ybottom = y1,
             x0 = rep(0, nCol),
             box = par("usr"),
             xMid = (x0[1, ] + x1[1, ])/2,
             yMid = (y0[, 1] + y1[, 1])/2);

  for (c in 1:nCol)
  {
    scale = (x1[1, c] - x0[1, c])/(matMax[c] - matMin[c])/(1+margin);
    cx0 = x0[, c] - matMin[c] * scale; # Position of the zero line. Note that matMin by def. cannot positive.
    cx1 = x0[, c] + mat[, c] * scale; # Position of the end of the bar.
    rect(cx0, y1[, c], cx1, y0[, c], col = colors[, c], border = border[, c]);

    if (!is.null(ablines))
      for (al in 1:length(ablines))
        if (ablines[al] >= matMin[c] && ablines[al] <= matMax[c])
          lines(rep(x0[1, c] + ablines[al] * scale, 2), c(y0[nRow, c], y1[1, c]), col = abColors[al],
                lty = ablines.lty[al])
    if (c > 1) lines(x = c(x0[1, c], x0[1, c]), y = c(y0[1, 1], y1[nRow, 1]), col = 1);
    out$x0[c] = cx0[1];
    if (!is.null(textMat))
    {
      onLeft = (cx1-x0[, c])/(x1[, c]-x0[, c]) > leftThreshold
      ytext = (y1[, c]  + y0[, c])/2
      if (alignRight)
      {
        xtext = x1[, c] - textMargin;
        xtext[onLeft] = pmax(cx0, cx1)[onLeft]  - textMargin
        adj.x = rep(1, nRow)
      } else {
        adj.x = ifelse(onLeft, 1, 0);
        xtext = pmax(cx0, cx1) + ifelse(onLeft, -1, 1) * textMargin;
        scaledMinXPosition = x0[1, c] + textMat.minXPosition * scale
        xtext[xtext < scaledMinXPosition + textMargin] = scaledMinXPosition + textMargin
      }
      for (r in 1:nRow)
        text(xtext[r], ytext[r], textMat[r, c], col = col.text, cex = cex.text, adj = c(adj.x[r], 0.5));
    }
  }

  invisible(out);
}


#---------------------------------------------------------------------------------------------------------
# labeledBarplotArray.R
#---------------------------------------------------------------------------------------------------------

#--------------------------------------------------------------------------
#
# .reverseRows = function(Matrix)
#
#--------------------------------------------------------------------------
#


.reverseRows = function(Matrix)
{
  ind = seq(from=dim(Matrix)[1], to=1, by=-1);
  Matrix[ind,];
  #Matrix
}

.reverseVector = function(Vector)
{
  ind = seq(from=length(Vector), to=1, by=-1);
  Vector[ind];
  #Vector
}
  
#--------------------------------------------------------------------------
#
# labeledBarplotArray
#
#--------------------------------------------------------------------------
# This function plots a heatmap of the specified matrix 
# and labels the x and y axes wit the given labels.
# It is assumed that the number of entries in xLabels and yLabels is consistent 
# with the dimensions in.
# If colorLabels==TRUE, the labels are not printed and instead interpreted as colors --
#  -- a simple symbol with the appropriate color is printed instead of the label.
# The x,yLabels are expected to have the form "..color" as in "MEgrey" or "PCturquoise".
# xSymbol, ySymbols are additional markers that can be placed next to color labels

#barplotArray = function(mat,
#                        lim = NULL,
#                        commonScale = FALSE,
#                        colors = "grey40",
#                        border = "black",
#                        ablines = NULL,
#                        abColors = NULL,
#                        ablines.lty = NULL,
#                        textMat = NULL,
#                        cex.text = 0.8,
#                        leftThreshold = 0.8,
#                        col.text = 1,
#                        textMargin = 0.04/ncol(mat),
#                        margin = 0.1,
#                        main = "",
#                        xlab = "", ylab = ""
#
#                        )


labeledBarplotArray = function (mat, 
                           signed = FALSE,
                           xLabels, yLabels = NULL, 
                           xSymbols = NULL, ySymbols = NULL, 
                           colorLabels = NULL, 
                           xColorLabels = FALSE, yColorLabels = FALSE,
                           checkColorsValid = TRUE,
                           xLabelsPosition = "bottom",
                           xLabelsAngle = 45,
                           xLabelsAdj = 1,
                           xColorWidth = 0.05,
                           yColorWidth = 0.05,
                           fixedColors = NULL,
                           colors = if (signed) greenWhiteRed(100) else greenWhiteRed(100)[50:100],
                           invertColors = FALSE, 
                           cex.lab = NULL, 
                           cex.lab.x = cex.lab,
                           cex.lab.y = cex.lab,
                           colors.lab.x = 1,
                           colors.lab.y = 1,
                           ... ) 
{
  if (!is.null(colorLabels)) {xColorLabels = colorLabels; yColorLabels = colorLabels; }
  
  if (is.null(yLabels) & (!is.null(xLabels)) & (dim(mat)[1]==dim(mat)[2])) 
    yLabels = xLabels; 

  if (checkColorsValid)
  {
    xValidColors = !is.na(match(substring(xLabels, 3), colors()));
    yValidColors = !is.na(match(substring(yLabels, 3), colors()));
  } else {
    xValidColors = rep(TRUE, length(xLabels));
    yValidColors = rep(TRUE, length(yLabels));
  }

  if (sum(xValidColors)>0) xColorLabInd = c(1:length(xLabels))[xValidColors]
  if (sum(!xValidColors)>0) xTextLabInd = c(1:length(xLabels))[!xValidColors]

  if (sum(yValidColors)>0) yColorLabInd = c(1:length(yLabels))[yValidColors]
  if (sum(!yValidColors)>0) yTextLabInd = c(1:length(yLabels))[!yValidColors]

  xLabPos = charmatch(xLabelsPosition, c("bottom", "top"));
  if (is.na(xLabPos))
    stop("Argument 'xLabelsPosition' must be (a unique abbreviation of) 'bottom', 'top'");

  if (is.null(colors)) colors = heat.colors(30);
  if (invertColors) colors = .reverseVector(colors);

  if (is.null(fixedColors))
  {
    fixedColors = numbers2colors(mat, colors = colors, signed = signed, lim = c(0, max(mat)));
  }
   

  # Call the barplotArray function

  labPos = barplotArray(mat, colors = fixedColors,
                        ...)

  # Draw the labels

  nxlabels = length(xLabels)
  plotbox = labPos$box;
  xmin = plotbox[1]; xmax = plotbox[2]; ymin = plotbox[3]; yrange = plotbox[4]-ymin;
  ymax = plotbox[4]; xrange = xmax - xmin;

  xspacing = labPos$xMid[2] - labPos$xMid[1];
  yspacing = abs(labPos$yMid[2] - labPos$yMid[1]);

  nylabels = length(yLabels)
  offsetx = par("cxy")[1] / 3;
  offsety = par("cxy")[2] / 3;
  # Transform fractional widths into coordinate widths
  xColW = min(xmax - xmin, ymax - ymin) * xColorWidth;
  yColW = min(xmax - xmin, ymax - ymin) * yColorWidth;
  if (sum(!xValidColors)>0)
  {
    xLabYPos = ifelse(xLabPos==1, ymin - offsety, ymax + offsety)
    if (is.null(cex.lab)) cex.lab = 1;
    text(labPos$xMid[xTextLabInd] , xLabYPos, srt = xLabelsAngle, 
          adj = xLabelsAdj, labels = xLabels[xTextLabInd], xpd = TRUE, cex = cex.lab.x, col = colors.lab.x)
  }
  if (sum(xValidColors)>0)
  {
    baseY = ifelse(xLabPos==1, ymin-offsety-xColW, ymax + offsety + xColW);
    deltaY = ifelse(xLabPos==1, xColW, -xColW);
    rect(xleft = labPos$xMid[xColorLabInd] - xspacing/2, ybottom = baseY,
         xright = labPos$xMid[xColorLabInd] + xspacing/2, ytop = baseY + deltaY,
         density = -1,  col = substring(xLabels[xColorLabInd], 3), 
         border = substring(xLabels[xColorLabInd], 3), xpd = TRUE)
    if (!is.null(xSymbols))
      text ( labPos$xMid[xColorLabInd], baseY - sign(deltaY)* offsety, xSymbols[xColorLabInd], 
             adj = xLabelsAdj, 
             xpd = TRUE, srt = xLabelsAngle, cex = cex.lab.x, col = colors.lab.x);
  }
  if (sum(!yValidColors)>0)
  {
    if (is.null(cex.lab)) cex.lab = 1;
    text(xmin - offsetx, labPos$yMid[yTextLabInd], srt = 0, 
         adj = c(1, 0.5), labels = yLabels[yTextLabInd], xpd = TRUE, cex = cex.lab.y, col = colors.lab.y )
  } 
  if (sum(yValidColors)>0)
  {
    rect(xleft = xmin- yColW - offsetx, ybottom = labPos$yMid[yColorLabInd] - yspacing/2,
         xright = xmin- offsetx, ytop = labPos$yMid[yColorLabInd] + yspacing/2, 
         density = -1,  col = substring(yLabels[yColorLabInd], 3), 
         border = substring(yLabels[yColorLabInd], 3), xpd = TRUE)
    if (!is.null(ySymbols))
      text (xmin- yColW - 2*offsetx, 
            labPos$yMid[yColorLabInd], ySymbols[yColorLabInd], 
            adj = c(1, 0.5), xpd = TRUE, cex = cex.lab.y, col = colors.lab.y);
  }

  axis(1, labels = FALSE, tick = FALSE)
  axis(2, labels = FALSE, tick = FALSE)
  axis(3, labels = FALSE, tick = FALSE)
  axis(4, labels = FALSE, tick = FALSE)
}

#==================================================================================================
#
# fixLabels
#
#=================================================================================================

# Replace (some) spaces in the given labels by newlines so that the length of each line is no more than
# maxCharPerLine


fixLabels = function(...)
{
  formatLabels(...)
}

nLines = function(s, eol = "\\n")
{
  sapply( strsplit( s, split = "\\n"),  length)
}


#==================================================================================================
#
# shortenLabels
#
#=================================================================================================

# Truncate labels at the last space before given maximum length, add ... if the label is shortened.


shortenLabels = function(labels, maxLength = 25, minLength = 10, split = " ", fixed = TRUE)
{
  n = length(labels);
  if (length(split) > 0)
  {
    newLabels= rep("", n);
    split = strsplit(labels, split = split, fixed = fixed);
    for (l in 1:n)
    {
      nl = "";
      s = 1; len = 0;
      while (s <= length(split[[l]]) && len <= maxLength)
      {
        newLen = len + nchar(split [[l]] [s]);
        if (len < minLength | newLen < maxLength)
        {
          nl = paste(nl, split[[l]] [s]);
          len = nchar(nl);
          s = s+1;
        } else {
          nl = spaste(nl, "...");
          len = maxLength + 1;
        }
      }
      newLabels[l] = nl;
    }
    substring(newLabels, 2);
  } else {
    shorten = nchar(labels) > maxLength;
    add = ifelse(shorten, "...", "");
    spaste(substring(labels, 1, maxLength), add);
  }
}


#===============================================================================================
#
# stretchLimits
#
#===============================================================================================

# Calculate minimum and maximum and stretch the limits by a given factor

stretchLim = function(x, rx, rx.low = rx, rx.high = rx)
{
  minx = min(x, na.rm = TRUE);
  maxx = max(x, na.rm = TRUE);
  c(minx - rx.low * (maxx-minx), maxx + rx.high * (maxx-minx));
}

stretchLimits = function(x, y, rx = 0.03, ry = 0.03)
{
  list(xlim = stretchLim(x, rx), ylim = stretchLim(y, ry));
}


#===============================================================================================
#
# nearestSNP
#
#===============================================================================================

# Find the nearest SNP and the corresponding distance to given gene positions.

nearestSNP = function(snpChr, snpBp, geneChr, geneBp)
{
  snpHasAnno = !is.na(snpChr) & !is.na(snpBp);
  geneHasAnno = !is.na(geneChr) & !is.na(geneBp);

  snpChr[!snpHasAnno] = -131;  # Hopefully no one will ever have chromosome -131 in their data. 
  geneChr[!geneHasAnno] = -132;  # Same for -132
  nGenes = length(geneChr);
  validChr = snpChr[snpChr!=-131]
  chrNames = sort(unique(validChr));
  nChr = length(chrNames);
  nSNPs = length(snpChr);
  nearestSNP = rep(NA, nGenes);
  for (chr in 1:nChr)
  {
    print(paste("chromosome", chrNames[chr]));
    snps = snpChr==chrNames[chr]
    pos = snpBp[snps];
    pos1 = pos[-1];
    breaks = (pos[-length(pos)] + pos1)/2;

    genes = c(1:nGenes)[geneChr==chrNames[chr] ];
    genePos = geneBp[genes];
    nearest = findInterval(genePos, breaks);
    nearestSNP[genes] = c(1:nSNPs)[snps][nearest+1];
  }
  nearestSNPdist = geneBp - snpBp[nearestSNP];
  list(nearestSNP = nearestSNP, nearestSNPdist = nearestSNPdist);
}

#===============================================================================================
#
# peakRegions
#
#===============================================================================================

# This function finds peaks and peak regions (defined as areas around the peak where the function drops less
# than a certain amount). Meant for QTL analysis.


peakRegions = function(qtl, chromo, basePair, minLOD, lodDrop, minPeakDistance)
{
  nSNPs = length(qtl);
  peaks = qtlPeaks(qtl, chromo, basePair, minQTL = minLOD, window = minPeakDistance);
  peakLocs = c(1:nSNPs)[peaks];
  nPeaks = sum(peaks)

  # For each peak, find the nearest base pairs up- and down-stream where the LOD score goes down by 1.

  peakChromo = chromo[peakLocs];
  peakBp = basePair[peakLocs];

  fromBp = rep(0, nPeaks);
  toBp = rep(0, nPeaks);
  fromLoc = rep(0, nPeaks);
  toLoc = rep(0, nPeaks);
  zzzzz = 10123;
  for (p in 1:nPeaks)
  {
    loc = peakLocs[p];
    baseLod = qtl[loc];
    while (loc > 0 && chromo[loc]==peakChromo[p] && qtl[loc] > baseLod - lodDrop)
      loc = loc -1;

    if (loc < 1) loc = 1;
  
    if (chromo[loc]!=peakChromo[p]) fromBp[p] = 0 else fromBp[p] = basePair[loc];
    fromLoc[p] = loc;
  
    loc = peakLocs[p];
    while (loc <= nSNPs && chromo[loc]==peakChromo[p] && qtl[loc] > baseLod - lodDrop)
      loc = loc +1;
  
    if (loc > nSNPs) loc = nSNPs;
    toLoc[p] = loc;
    if (chromo[loc]!=peakChromo[p]) toBp[p] = 1e10 else toBp[p] = basePair[loc];
  }

  list(peakIndicator = peaks, peakIndex = peakLocs, 
       peakChromo = peakChromo, peakBp = peakBp,
       peakStartInd = fromLoc,
       peakEndInd = toLoc,
       peakStartBp = fromBp,
       peakEndBp = toBp);
}


#===============================================================================================
#
# makeCross
#
#===============================================================================================

# prepare cross data
# First column in both genotypes and phenotypes must be sample id
# Second column in phenotypes must be sex

makeCross = function(genotypes, traits, markerNames = NULL, 
                       chromo = NULL, basePair = NULL, fileNameBase = NULL,
                       outlierZ = 3, genoCodes = c(1,2,3), alleles = c("B", "C"))
{

  commonSamples = intersect(genotypes[, 1], traits[, 1]);
  if (length(commonSamples) < 10)
    stop("Something's wrong: less than 10 samples common between genotypes and traits.");
  genotypes = genotypes[match(commonSamples, genotypes[, 1]), ];
  traits = traits[match(commonSamples, traits[, 1]), ];
  
  if (is.null(fileNameBase))
  {
    delete = TRUE
    genoFile = tempfile("genotypes", getwd());
    phenFile = tempfile("phenotypes", getwd());
  } else {
    delete = FALSE
    genoFile = spaste(fileNameBase, "-genotypes.csv");
    phenFile = spaste(fileNameBase, "-phenotypes.csv");
  }

  if (is.null(markerNames))
  {
    # Determine marker names from colnames of genotypes
    markerNames = colnames(genotypes);
    if (is.null(chromo)) # Assume marker names have the form marker.chrN.BpNNNNNN
    {
      split = sapply(strsplit(markerNames[-1], split = ".", fixed = TRUE), I);
      markerNames = split[1, ];
      chromo = as.numeric(substring( split[2, ], 4));
      basePair = as.numeric(substring( split[3, ], 3));
    }
  }

  SNPinfo = rbind( c("", chromo), c("", basePair));
  colnames(SNPinfo) =c(colnames(genotypes)[1], markerNames);
  colnames(genotypes) = colnames(SNPinfo);

  genoForQTL = rbind(SNPinfo, genotypes);

  write.csv(genoForQTL, genoFile, row.names = FALSE, quote = FALSE);

  if (outlierZ > 0)
  {
     stdTraits = scale(apply(traits[, -c(1:2)], 2, as.numeric)); # Do not scale mouse ID and sex:)
     outliers = abs(stdTraits) > outlierZ
     traits2 = as.data.frame(traits);
     traits2[, -c(1:2)][outliers] = NA;
  } else 
     traits2 = traits;
  
  write.csv(traits2, phenFile, row.names = FALSE, quote = FALSE);

  cross = read.cross(format = "csvs",
                     genfile = genoFile,
                     phefile = phenFile,
                     genotypes = genoCodes, alleles = alleles);
   
  if (delete)
  {
    unlink(genoFile);
    unlink(phenFile);
  }
  list(cross = cross, markerNames = markerNames, chromo = chromo, basePair = basePair);
}


#================================================================================================
#
# add x- and y- error bars to a plot
#
#================================================================================================

addErrorBars.xy = function(x, y, sdx, sdy, barWidth = 0.01, lwd = 1, col = 1)
{
  box = par("usr");
  barWidth.x = (box[2] - box[1])*barWidth;
  barWidth.y = (box[4] - box[3])*barWidth;

  n = length(x);

  if (length(y)!=n) stop("Lengths of x and y must be the same.");
  if (length(sdx)!=n) stop("Lengths of x and sdx must be the same.");
  if (length(sdy)!=n) stop("Lengths of x and sdy must be the same.");

  for (i in 1:n) if ( is.finite(x[i]) && is.finite(y[i]))
  {
    if (is.finite(sdx[i]))
    {
      segments(x[i]-sdx[i], y[i], x[i] + sdx[i], y[i], lwd = lwd, col = col);
      segments(x[i]+sdx[i], y[i] - barWidth.y/2, x[i]+sdx[i], y[i] + barWidth.y/2, lwd = lwd, col = col);
      segments(x[i]-sdx[i], y[i] - barWidth.y/2, x[i]-sdx[i], y[i] + barWidth.y/2, lwd = lwd, col = col);
    }
    if (is.finite(sdy[i]))
    {
      segments(x[i], y[i]-sdy[i], x[i], y[i]+sdy[i], lwd = lwd, col = col);
      segments(x[i]-barWidth.x/2, y[i] - sdy[i], x[i]+barWidth.x/2, y[i] - sdy[i], lwd = lwd, col = col);
      segments(x[i]-barWidth.x/2, y[i] + sdy[i], x[i]+barWidth.x/2, y[i] + sdy[i], lwd = lwd, col = col);
    }
  }
}

#================================================================================================
#
# verboseMultiplot: overlay several verboseScatterplots
#
#================================================================================================

verboseMultiplot = function(x, y, sample = NULL, corFnc = "cor", corOptions = "use = 'p'", 
    colors = c(1:length(x)),
    main = "", xlab = NA, ylab = NA, cex = 1, cex.axis = 1.5, 
    cex.lab = 1.5, cex.main = 1.5, abline = FALSE, abline.color = colors, 
    abline.lty = 1, corLabel = corFnc, xlim = NULL, ylim = NULL, 
    pch = rep(1, length(x)), ...) 
{
    if (is.na(xlab)) 
        xlab = as.character(match.call(expand.dots = FALSE)$x)
    if (is.na(ylab)) 
        ylab = as.character(match.call(expand.dots = FALSE)$y)
    if (mode(x)!="list") { x = list(x); y = list(y); }
    nSets = length(x);
    if (nSets!=length(y))
      stop("Number of components in 'x' and 'y' must be the same.");
    cor = corp = rep(NA, nSets);
    for (set in 1:nSets)
    {
      corExpr = parse(text = paste(corFnc, "(x[[set]], y[[set]] ", prepComma(corOptions), ")"))
      cor[set] = signif(eval(corExpr), 2)
      corp1 = signif(corPvalueStudent(cor[set], sum(is.finite(x[[set]]) & is.finite(y[[set]]))), 2)
      if (corp1 < 10^(-200)) 
          corp[set] = "<1e-200"
      else corp[set] = spaste("=", corp1);
    }
    if (!is.na(corLabel)) {
        mainX = paste(main, " ", corLabel, "=", cor, ", p", corp, 
            sep = "")
      } else 
         mainX = main

    x1 = unlist(x);
    y1 = unlist(y);
    fin = is.finite(x1) & is.finite(y1);
    if (is.null(xlim))
    {
      minx = min(x1[fin]);
      maxx = max(x1[fin]);
      xlim = c(minx, maxx);
    }
    if (is.null(ylim))
    {
      miny = min(y1[fin]);
      maxy = max(y1[fin]);
      ylim = c(miny, maxy);
    }

    for (set in 1:nSets)
    {
      if (!is.null(sample)) 
      {
        if (length(sample) == 1) {
            sample = sample(length(x[[set]]), sample)
        }
      } else
        sample = c(1:length(x[[set]]));
      if (set==1)
      {
        plot(x[[set]] [sample], y[[set]] [sample], main = mainX, xlab = xlab, 
            ylab = ylab, cex = cex, cex.axis = cex.axis, cex.lab = cex.lab, 
            cex.main = cex.main, xlim = xlim, ylim = ylim, col = colors[set], 
            pch = pch[set], ...)
      } else 
        points(x[[set]] [sample], y[[set]] [sample], col = colors[set], pch = pch[set])
 
      if (abline) 
      {
        fit = lm(y[[set]] ~ x[[set]])
        abline(reg = fit, col = abline.color[set], lty = abline.lty)
      }
    }
}


# Test:

if (FALSE)
{
  x = list(rnorm(20), rnorm(20)+2, rnorm(10)-1);
  y = list(rnorm(20), rnorm(20)+2, rnorm(10)-1);
  y[[1]] = y[[1]] +x[[1]];
  y[[2]] = 0.5 * y[[2]] +x[[2]];
  y[[3]] = y[[3]] - x[[3]];

  verboseMultiplot(x, y, abline = TRUE)
  
}

#=========================================================================================================
#
# Plot module heatmap and eigengene barplot
#
#=========================================================================================================

# This needs finishing

moduleHeatmapAndEigengeneBarplot = function(
   expr, colors = greenWhiteRed(50, gamma = 0.3), geneOrder = NULL, sampleOrder = NULL, 
   setLayout = TRUE, layoutHeights = c(0.6, 0.4), marHeatmap = c(0,3.8,3,1), marBarplot = c(2,3.8,0,1),
   maxNShowGenes = 500, scale = TRUE, main = "" )
{
   if (setLayout) layout(matrix(c(1,2), 2, 1), heights= layoutHeights);
   nGenes = ncol(expr);
   if (nGenes > maxNShowGenes )
   {
     sample = sample(c(1:nGenes), size = maxNShowGenes, replace = FALSE);
     modExpr = expr[, sample];
     nGenes = maxNShowGenes;
   } else
     modExpr = expr;
   if (scale) modExpr = scale(modExpr);
   if (is.null(geneOrder)) geneOrder = hclust(as.dist(1-cor(modExpr, use="p")), method = "average")$order
   if (is.null(sampleOrder)) sampleOrder = hclust(dist(modExpr), method = "average")$order;
   maxZ = max(abs(modExpr), na.rm = TRUE);
   exprColors = numbers2colors(modExpr, signed = TRUE, lim = c(-maxZ, maxZ), colors = colors);
   par(mar = marHeatmap);
   eigengene = moduleEigengenes(modExpr, rep(1, nGenes))$eigengenes[, 1];
   mp = barplot(as.vector(eigengene[sampleOrder]), col = "white", border = "white", axisnames = FALSE,
                main = main, 
                axes = FALSE);
   plotbox = par("usr");
   xstep = mp[2]-mp[1]; xLeft = mp - xstep/2; xRight = mp + xstep/2;

   nrows = ncol(modExpr);
   yRange = plotbox[4]-plotbox[3];
   yBot = plotbox[3] + c(0:(nrows-1)) * yRange/nrows;
   yTop = yBot + yRange/nrows;

   for (sample in 1:nrow(modExpr))
   {
     rect(xleft = rep(xLeft[sample], nrows), xright = rep(xRight[sample], nrows),
          ybottom = yBot, ytop = yTop, col = exprColors[ sampleOrder[sample], geneOrder], 
          border = exprColors[ sampleOrder[sample], geneOrder]);
   }
   par(mar = marBarplot);
   eigengene[is.na(eigengene)] = 0;
   colors = numbers2colors(eigengene, signed = TRUE, colors = colors)
   barplot(as.vector(eigengene[sampleOrder]), col = colors[sampleOrder],
           xlab = "", ylab = "Eigengene expression");
}


#========================================================================================================
#
# Convenient prediction accuracy functions
#
#========================================================================================================

naiveAccuracyFnc = function(response)
{
  max(table(response))/sum(table(response))
}

accuraciesFnc = function(confusion)
{
  c (diag(confusion) / rowSums(confusion), sum(diag(confusion))/sum(confusion) );
}


table2.allLevels = function(x, y, levels.x = sort(unique(x)), levels.y = sort(unique(y)), setNames = FALSE)
{
  nx = length(levels.x);
  ny = length(levels.y);
  t = table(x, y);

  out = matrix(0, nx, ny);
  if (setNames)
  {
    rownames(out) = levels.x;
    colnames(out) = levels.y;
  }
  out[ match(rownames(t), levels.x), match(colnames(t), levels.y) ] = t;
  out;
}

#========================================================================================================
#
# Log-transformation of data
#
#========================================================================================================

# this function replaces negative values

logTrafo = function(data, base = 2, zeroOrNegativeAction = c("NA", "min"), minCoeff = 0.5)
{
  origDim = dim(data);
  origNames = dimnames(data);
  data = as.matrix(data)
  action = match.arg(zeroOrNegativeAction)
  missing = is.na(data);
  data[data<=0] = NA;
  if (action=="min")
  {
    mins = apply(data,2,min, na.rm = TRUE);
    minMat = matrix(mins*minCoeff, nrow(data), ncol(data), byrow = TRUE);
    data[is.na(data)] = minMat[is.na(data)];
  }
  out = logb(data, base);
  out[missing] = NA;
  dim(out) = origDim;
  dimnames(out) = origNames;
  out;
}

# These two functions shift individual columns of data where necessary to make a log of negative value
# meaningful.

logTrafoProcessedData.1 = function(x, base = 2)
{
  if (any(x < 0, na.rm = TRUE))
    x = x-min(x, na.rm = TRUE);

  small = x<base;
  small[is.na(small)] = TRUE;
  x[ small ] = x[small]/base;
  x[!small] = logb( x[!small], base);
  x;
}

logTrafoProcessedData = function(data, base=2)
{
  out = apply(data, 2, logTrafoProcessedData.1, base=base);
  dimnames(out) = dimnames(data);
  out;
}


#========================================================================================================
#
# LinerarOrLog-transformation of data
#
#========================================================================================================

linearOrLogTrafo = function(data)
{
  fin = is.finite(data);
  data[!fin] = 0;
  above2 = data > 2;
  below2 = !above2;
  data[below2] = data[below2]/2
  data[above2] = log2(data[above2]);
  data[!fin] = NA;
  data;
}



#=====================================================================================================
#
# removeOutliers (using inter-sample connectivity)
#
#=====================================================================================================

outliers = function(expr, adjacencyFnc = adjacency, adjacencyOptions = list(),
                          Z.min = -3, returnZ = FALSE)
{
  adjacencyFnc = match.fun(adjacency);
  adjacencyOptions$datExpr = t(expr);
  adj = do.call(adjacencyFnc, adjacencyOptions)
  Z = scale(colSums(adj)-1);
  if (returnZ)
  {
     list(outliers = Z<Z.min, Z = Z)
  } else 
     Z<Z.min
}

#=====================================================================================================
#
# goGenes
#
#=====================================================================================================

# Returns a list of gene identifiers in the given GO category.

GOgenes = function(
  termNames = NULL,
  organism = "human",
  ontologies = c("BP", "CC", "MF"),
  evidence = c("IMP", "IGI", "IPI", "ISS", "IDA", "IEA", "TAS", "NAS", "ND", "IC"),
  includeOffspring = TRUE, 
  verbose = 2, indent = 0)
{
   organisms = c("human", "mouse", "rat", "malaria", "yeast", "fly", "bovine", "worm", "canine",
                 "zebrafish", "chicken");
   allEvidence =  c("IMP", "IGI", "IPI", "ISS", "IDA", "IEA", "TAS", "NAS", "ND", "IC");
   allOntologies = c("BP", "CC", "MF");

   spaces = indentSpaces(indent);
   orgInd = pmatch(organism, organisms);
   if (is.na(orgInd))
     stop(paste("Unrecognized 'organism' given. Recognized values are ",
                paste(organisms, collapse = ", ")));

   if (length(evidence)==0)
     stop("At least one valid evidence code must be given in 'evidence'.");
   if (length(ontologies)==0)
     stop("At least one valid ontology code must be given in 'ontology'.");

   evidInd = pmatch(evidence, allEvidence);
   if (sum(is.na(evidInd))!=0)
     stop(paste("Unrecognized 'evidence' given. Recognized values are ",
                paste(allEvidence, collapse = ", ")));

   ontoInd = pmatch(ontologies, allOntologies);
   if (sum(is.na(ontoInd))!=0)
     stop(paste("Unrecognized 'ontologies' given. Recognized values are ",
                paste(allEvidence, collapse = ", ")));

   orgCodes = c("Hs", "Mm", "Rn", "Pf", "Sc", "Dm", "Bt", "Ce", "Cf", "Dr", "Gg");
   orgExtensions = c(rep(".eg", 4), ".sgd", rep(".eg", 6));
   reverseMap = c(rep(".egGO2EG", 4), ".sgdGO2ORF", rep(".egGO2EG", 6))

   missingPacks = NULL;
   packageName = paste("org.", orgCodes[orgInd], orgExtensions[orgInd], ".db", sep="");
   if (!require(packageName, character.only = TRUE))
     missingPacks = c(missingPacks, packageName);

   if (!require(GO.db))
     missingPacks = c(missingPacks, "GO.db");

   if (!is.null(missingPacks))
     stop(paste("Could not load the requisite package(s)",
           paste(missingPacks, collapse = ", "), ". Please install the package(s)."))

   if (verbose > 0)
   {
     printFlush(paste(spaces, "goGenesInCategory: loading annotation data..."));
   }

   egGO = eval(parse(text = paste(packageName, ":::org.", orgCodes[orgInd], orgExtensions[orgInd],
                                  "GO", sep = "")));

   Go2eg = eval(parse(text = paste("AnnotationDbi::as.list(", packageName, ":::org.", orgCodes[orgInd],
                                           reverseMap[orgInd],")", sep = "")));
   nTerms = length(Go2eg);

   goInfo = eval(parse(text = "AnnotationDbi::as.list(GO.db:::GOTERM)"));
   if (length(goInfo) > 0)
   {
      orgGoNames = names(Go2eg);
      dbGoIDs = as.character(sapply(goInfo, GOID));
      dbGoOntologies = as.character(sapply(goInfo, Ontology));
      dbGoTerm = as.character(sapply(goInfo, Term));
   } else {
      dbGoIDs = "";
   }

   goOffSpr = list();
   if (includeOffspring)
   {
     goOffSpr[[1]] = eval(parse(text = "AnnotationDbi::as.list(GO.db:::GOBPOFFSPRING)"));
     goOffSpr[[2]] = eval(parse(text = "AnnotationDbi::as.list(GO.db:::GOCCOFFSPRING)"));
     goOffSpr[[3]] = eval(parse(text = "AnnotationDbi::as.list(GO.db:::GOMFOFFSPRING)"));
   }
   term2info = match(names(Go2eg), names(goInfo));
   orgOntologies = dbGoOntologies[term2info];
   orgIDs = dbGoIDs[term2info];
   orgTerms = dbGoTerm[term2info];

   if (is.null(termNames)) termNames = orgTerms;

   term2orgTerm = match(tolower(termNames), tolower(orgTerms));
   if (any(is.na(term2orgTerm)))
   {
      missingTermNames = term2orgTerm[is.na(term2orgTerm)];
      term2orgID = match(missingTermNames, dbGoIDs);
      if (any(is.na(term2orgID)))
         printFlush(spaste("Warning in GOGenesInCategory: the following terms were found ",
                           "in neither GO names nor GO term IDs:\n",
                           paste(missingTermNames[is.na(term2orgID)], collapse = ", ")));
      term2orgTerm[is.na(term2orgTerm)] = term2orgID;
   }

   keepInputTerms = is.finite(term2orgTerm);
   nInputTerms = length(termNames);

   termCodes = vector(mode="list", length = nInputTerms);
   for (t in 1:nInputTerms) termCodes[[t]] = NA;
   names(termCodes) = termNames;

   termIDs = rep(NA, nInputTerms);
   termIDs[keepInputTerms] = orgIDs[term2orgTerm[keepInputTerms]];

   collectGarbage();
   nExpandLength = 0;
   blockSize = 1000; # For a more efficient concatenating of offspring genes
   nAllInTerm = rep(0, nTerms);

   if (verbose > 0)
   {
      cat(paste(spaces, " ..preparing term lists..."));
      if (nInputTerms > 10) pind = initProgInd();
   }

   for (t in 1:nInputTerms) if (keepInputTerms[t] && !is.na(Go2eg[[ term2orgTerm[t] ]] [[1]]))
   {
      c = term2orgTerm[t];
      te = as.character(names(Go2eg[[c]])); # Term evidence codes
      tc = Go2eg[[c]]; # Gene codes
      if (includeOffspring)
      {
        termOffspring = NULL;
        for (ont in 1:length(goOffSpr))
        {
          term2off = match(names(Go2eg)[c], names(goOffSpr[[ont]]))
          if (!is.na(term2off))
            termOffspring = c(termOffspring, goOffSpr[[ont]][[term2off]]);
        }
        if (length(termOffspring)>0)
        {
           maxLen = blockSize;
           tex = rep("", maxLen);
           tcx = rep("", maxLen);
           ind = 1;
           len = length(te);
           tex[ ind:len ] = te;
           tcx[ ind:len ] = tc;
           ind = len + 1;
           o2go = match(termOffspring, as.character(names(Go2eg)));
           o2go = o2go[is.finite(o2go)]
           if (length(o2go)>0) for (o in 1:length(o2go)) if (!is.na(Go2eg[[o2go[o]]][[1]]))
           {
             #printFlush(paste("Have offspring for term", c, ": ", names(Go2eg)[c], 
             #           Term(goInfo[[term2info[c]]])));
             newc = Go2eg[[o2go[o]]];
             newe = names(newc);
             newl = length(newe);
             if ((len + newl) > maxLen)
             {
               nl = len + newl;
               maxLen = blockSize * ceiling( nl/blockSize);
               tex = c(tex, rep("", maxLen - length(tex)));
               tcx = c(tcx, rep("", maxLen - length(tex)));
               nExpandLength = nExpandLength + 1;
             }
             tex[ind:(len + newl)] = newe;
             tcx[ind:(len + newl)] = newc;
             ind = ind + newl;
             len = len + newl;
           }
           te = tex[1:len];
           tc = tcx[1:len];
        }
      }
      use = is.finite(match(te, evidence));
      termCodes[[t]] = unique(as.character(tc[use]));
      if (nInputTerms > 10 && verbose > 0) pind = updateProgInd(t/nInputTerms, pind);
   }

   if (nInputTerms > 10 && verbose > 0)
   {
      pind = updateProgInd(1, pind);
      printFlush("");
   }

   attr(termCodes, "GOtermID") = termIDs

   termCodes;
}


#====================================================================================================
#
# topEnrichment
#
#====================================================================================================

# Enrichment of top genes in a given list.

topEnrichment = function(ranking, referenceList, 
                         nTopGenes = NULL, geneIDs = NULL, direction = c("increasing", "decreasing"),
                         chiSquare = FALSE,
                         verbose = 1)
{
  direction = match.arg(direction);
  ranking = as.matrix(ranking);
  if (is.null(geneIDs)) geneIDs = rownames(ranking);
  if (is.null(geneIDs)) stop("geneIDs must be specified if 'ranking' does not have row names.");
  if (length(intersect(referenceList, geneIDs))==0)
     stop("There are no common genes between referenceList and geneIDs.");

  if (is.null(nTopGenes)) nTopGenes = c(10, seq(from=20, to=1000, by=20))
  sign = ifelse(direction=="increasing", 1, -1);
  ranking = apply(sign * ranking, 2, rank, na.last = TRUE);

  nRankings = ncol(ranking);
  nNumbers = length(nTopGenes);

  inReference = geneIDs %in% referenceList;
  if (length(unique(inReference))==1) 
    stop("referenceList contains either no or all genes in 'ranking'.");

  enrichmentStat = matrix(0, nNumbers, nRankings);

  if (verbose > 0) pind = initProgInd();
  for (r in 1:nRankings)
    for (n in 1:nNumbers)
    {
      inTop = ranking[, r] <= nTopGenes[n];
      if (length(unique(inTop)) > 1)
      {
        if (chiSquare) {
          enrichmentStat[n, r] = chiSquare.table(table(inReference, inTop));
        } else {
          enrichmentStat[n, r] = -log10(fisher.test(table(inReference, inTop), 
                                      alternative = "greater")$p.value);
        }
      }
      if (verbose > 0) pind = updateProgInd(((r-1)*nNumbers + n)/(nRankings * nNumbers), pind);

    }

  if (verbose > 0) printFlush("");

  colnames(enrichmentStat) = colnames(ranking);
  rownames(enrichmentStat) = spaste("nTopGenes.", nTopGenes);
  attr(enrichmentStat, "nTopGenes") = nTopGenes;
  enrichmentStat;
}


#====================================================================================================
#
# chiSquare.table
#
#====================================================================================================

chiSquare.table = function(tab)
{
  if (sum(tab<5) > 0)
    return(qchisq(fisher.test(tab, alt = "greater")$p.value), df = 1, lower.tail = FALSE);

  n = sum(tab);
  expected = outer(rowSums(tab), colSums(tab))/n;
  chsq = sum( (expected - tab) ^2 / expected );
  chsq;
}

#====================================================================================================
#
# topEnrichment
#
#====================================================================================================

# Enrichment of top genes in a given list.




#====================================================================================================
#
# plotEnrichments
#
#====================================================================================================

plotEnrichments = function(enrichments, counts = attr(enrichments, "nTopGenes"),
                           colors = c(1:ncol(enrichments)),
                           ylim = NULL, pch = c(1:ncol(enrichments)), 
                           lty = rep(1, ncol(enrichments)), 
                           plotLegend = TRUE,
                           leg.position = "auto",
                           legend = colnames(enrichments),
                           cex.legend = 1,
                           legend.ncol = 1,
                           ...)
{
  max = max(enrichments, na.rm = TRUE);
  if (is.null(ylim)) ylim = c(0, max)
  plot(counts, enrichments[, 1], col = colors[1], type = 'l', lty = lty[1], ylim = ylim, ...);
  points(counts, enrichments[, 1], col = colors[1], pch = pch[1], bg = colors[1]);
  nPlots = ncol(enrichments);
  if (nPlots > 1) for (c in 2:nPlots)
  {
    lines(counts, enrichments[, c], col = colors[c], lty = lty[c]);
    points(counts, enrichments[, c], col = colors[c], pch = pch[c], bg = colors[c]);
  }

  if (plotLegend)
    legendClean(leg.position, legend = legend, lty = lty, pch = pch, col = colors, pt.bg = colors, 
                cex = cex.legend, ncol = legend.ncol,
                points.x = rep(counts, nPlots), points.y = c(enrichments));
}

plotEnrichments.barplot = function(enrichments, counts = attr(enrichments, "nTopGenes"),
                           nTop = as.integer(length(counts)/5),
                           colors = c(1:ncol(enrichments)),
                           ylim = NULL, 
                           legend = colnames(enrichments),
                           barplot = TRUE,
                           top = c("highest", "lowest"),
                           reverseYAxis = NULL,
                           horiz = FALSE,
                           # Ornaments: braced and delimeted groups
                           # Note: ornaments won't work for horizontal barplots
                           braceEdges = NULL,
                           braceText = NULL,
                           braceText.col = 1,
                           braceText.cex = 1,
                           braceSepLine = TRUE,
                           braceSepLine.lty = 2,
                           braceSepLine.col = "darkgrey",
                           braceHeight = 0.05,
                           braceSpace = 0.25,
                           bracePower = 13,
                           braceCol = braceText.col,
                           # Simple groups
                           groupEdges = NULL,
                           groupText = NULL, 
                           groupText.col = 1,
                           groupText.cex = 1,
                           groupSpace = 0.1,
                           ...)
{

  top = match.arg(top);
  if (is.null(reverseYAxis)) reverseYAxis = top=="lowest";

  ordEnr = apply(enrichments, 2, sort, decreasing = (top=="highest"));
  topEnr = ordEnr[1:nTop, ];
  means = colMeans(topEnr);
  stdErr = apply(topEnr, 2, stdErr)
  nBars = ncol(enrichments);

  max = max(means + stdErr, na.rm = TRUE);
  if (is.null(ylim))
  {
    if (barplot) 
    {
       ylim = c(0, max)
    } else {
       ylim = range(topEnr);
       if (reverseYAxis) ylim = rev(ylim);
    }
  }

  yLimExt = 0;

  if (!is.null(braceEdges))
  {
    yLimExt = yLimExt + braceSpace;
    braceStart = 1/(1+braceSpace);
  } else 
    braceSpace = 0;

  if (!is.null(groupEdges))
  {
    yLimExt = yLimExt + groupSpace;
    groupMiddle = (1 + 0.5 * groupSpace)/(1+groupSpace+braceSpace);
    braceStart = (1+groupSpace)/(1+groupSpace + braceSpace);
  } else
    groupSpace = 0;

  ylim[2] = ylim[2] + (ylim[2]-ylim[1]) * yLimExt;
  

  x = c(topEnr);
  g = rep(c(1:nBars), rep(nTop, nBars));

  if (barplot) {
    if (horiz)
    {
      mp = labeledBarplot2(x, g, col = colors,  xlim = ylim, names = legend, horiz = horiz, ...);
    } else {
      mp = labeledBarplot2(x, g, col = colors,  ylim = ylim, names = legend, horiz = horiz, ...);
    }
    yTop = means+stdErr;
  } else {
    bxp = labeledBoxplot(as.list(as.data.frame(topEnr)), names = legend, namesColor = colors, verbose = TRUE, 
                   border = colors, ylim = ylim, ...);
    yTop = apply(bxp$stats, 1, max, na.rm = TRUE);
    mp = c(1:ncol(topEnr));
  }

  addOrnaments(mp, yTop,
    braceEdges = braceEdges,
    braceText = braceText,      
    braceText.col = braceText.col,     
    braceText.cex = braceText.cex,     
    braceSepLine = braceSepLine,   
    braceSepLine.lty = braceSepLine.lty,  
    braceSepLine.col = braceSepLine.col,
    braceHeight = braceHeight,    
    braceSpace = braceSpace,     
    bracePower = bracePower,       
    braceCol = braceCol,
    # Simple groups        
    groupEdges = groupEdges,     
    groupText = groupText,      
    groupText.col = groupText.col,     
    groupText.cex = groupText.cex,     
    groupSpace = groupSpace);

  invisible(list(topEnrichment = topEnr, midpoints = mp, ylim = ylim));
}

# Add ornaments to a barplot/enrichment barplot

addOrnaments = function(
    midPoints,
    heights,
    braceEdges = NULL,
    braceText = NULL,
    braceText.col = 1,
    braceText.cex = 1,
    braceSepLine = TRUE,
    braceSepLine.lty = 2,
    braceSepLine.col = "darkgrey",
    braceHeight = 0.05,
    braceSpace = 0.25,
    bracePower = 13,
    braceCol = braceText.col,
    # Simple groups
    groupEdges = NULL,
    groupText = NULL,
    groupText.col = 1,
    groupText.cex = 1,
    groupSpace = 0.1)
{

  if (!is.null(braceEdges))
  {
    braceStart = 1/(1+braceSpace);
  } else
    braceSpace = 0;

  if (!is.null(groupEdges))
  {
    groupMiddle = (1 + 0.5 * groupSpace)/(1+groupSpace+braceSpace);
    braceStart = (1+groupSpace)/(1+groupSpace + braceSpace);
  } else
    groupSpace = 0;

  nBars = length(midPoints);

  box = par("usr");
  ymin = box[3];
  ymax = box[4]
  yrange = ymax - ymin;
  barSpace = midPoints[2]-midPoints[1];
  if (!is.null(braceEdges))
  {
    braceText.cex = .extend(braceText.cex, nBars);
    braceText.col = .extend(braceText.col, nBars);
    nBraces = ncol(braceEdges);
    for (b in 1:nBraces)
    {
      st = braceEdges[1, b]; en = braceEdges[2, b]; 
      if (!is.null(braceCol))
        brace.horizontal(midPoints[ st ]-0.3*barSpace, midPoints[ en ] + 0.3*barSpace, 
                         ymin + yrange * braceStart, ymin + yrange * (braceStart + braceHeight), 
                         power = bracePower, col = braceCol);
      if (braceSepLine && braceEdges[2, b] < nBars)
         abline(v = (midPoints[en] + midPoints[en +1])/2, col = braceSepLine.col, lty = braceSepLine.lty);
      text((midPoints[st] + midPoints[en])/2, ymax, braceText[b], adj = c(0.5, 1), cex = braceText.cex[b],
           col = braceText.col[b]);
    }
  }

  if (!is.null(groupEdges))
  {
    groupText.col = .extend(groupText.col, nBars)
    groupText.cex = .extend(groupText.cex, nBars)
    nGroups = ncol(groupEdges);
    for (gr in 1:nGroups)
    {
       st = groupEdges[1, gr]; en = groupEdges[2, gr];
       text((midPoints[st] + midPoints[en])/2, 
            max(0, heights[st:en]) + 0.5 * yrange * groupSpace/(1+groupSpace + braceSpace), 
            groupText[gr], adj = c(0.5, 0.5),
            cex = groupText.cex[gr], col = groupText.col[gr]);
    }
  }
    
}


.extend = function(x, n)
{
  nRep = ceil(n/length(x));
  rep(x, nRep)[1:n];
}

#======================================================================================================
#
# Validation success
#
#======================================================================================================

# This form of the function returns the mean observed value and mean rank of observed values for a given
# vector of top number of genes in predicted.
# In other words: for each column of predicted: the predicted values are ranked and nTopGenes is selected;
# the mean of observed and rank(observed) for these values are calculated. 

validationSuccess = function(predicted, observed, nTopGenes, rankPredicted = TRUE,
                             direction = c("increasing", "decreasing"))
{
  direction = match.arg(direction);
  rankSign = ifelse(direction=="decreasing", -1, 1);
  nNGenes = length(nTopGenes)
  nMethods = 2
  methodNames = c("AverageObserved", "AverageObservedRank");

  predicted = as.matrix(predicted);
  nRankings = ncol(predicted);
  if (is.null(colnames(predicted)))
  {  
    rankingNames= spaste("Ranking.", c(1:nRankings));
  } else
    rankingNames = colnames(predicted);

  success = array(NA, dim = c(nNGenes, nRankings, nMethods));
  observedRank = rank(observed * rankSign, na.last = "keep");

  if (rankPredicted) {
    predictedRank = apply(predicted * rankSign, 2, rank);
  } else 
    predictedRank = predicted;
 
  for (ing in 1:nNGenes) for (r in 1:nRankings) 
  {
    ng = nTopGenes[ing];
    keep = predictedRank[, r] <=ng;
    success[ing, r, 1] = mean(observed[keep], na.rm = TRUE);
    success[ing, r, 2] = mean(observedRank[keep], na.rm = TRUE);
  }
  
  dimnames(success) = list(spaste("nTopGenes.", nTopGenes), rankingNames, methodNames);
  success = success[ , , , drop = TRUE];
  attr(success, "nTopGenes") = nTopGenes
  success;
}

# The second form of the function returns only the success measured on correlations but also returns the
# actual values that were averaged. 

validationSuccess.ext = function(predicted, observed, nTopGenes, rankPredicted = TRUE,
                                 direction = c("increasing", "decreasing"))
{
  direction = match.arg(direction);
  rankSign = ifelse(direction=="decreasing", -1, 1);
  nNGenes = length(nTopGenes)

  predicted = as.matrix(predicted);
  nRankings = ncol(predicted);
  if (is.null(colnames(predicted)))
  {
    rankingNames= spaste("Ranking.", c(1:nRankings));
  } else
    rankingNames = colnames(predicted);

  success = array(NA, dim = c(nNGenes, nRankings));

  if (rankPredicted) {
    predictedRank = apply(predicted * rankSign, 2, rank);
  } else
    predictedRank = predicted;

  topObservedValues = matrix(0,max(nTopGenes), nRankings);

  for (r in 1:nRankings)
  {
    for (ing in 1:nNGenes) 
    {
      ng = nTopGenes[ing];
      keep = predictedRank[, r] <=ng;
      success[ing, r] = mean(observed[keep], na.rm = TRUE);
    }
    ng = max(nTopGenes);
    order = order(predicted[, r] * rankSign);
    topObservedValues[, r] = observed[order[1:ng]];
  }

  dimnames(success) = list(spaste("nTopGenes.", nTopGenes), rankingNames);
  colnames(topObservedValues) = rankingNames;
  #success = success[ , , drop = TRUE];
  attr(success, "nTopGenes") = nTopGenes
  list(success = success, topObservedValues = topObservedValues, nTopGenes = nTopGenes);
}

#====================================================================================================#
#
# replace NAs by a finite value
#
#====================================================================================================

replaceNA = function(x, replaceBy=-9)
{
  x[is.na(x)] = replaceBy;
  x;
}

#====================================================================================================
#
# brace
#
#====================================================================================================

# draw a brace in a plot area
# For now only a horizontal brace

brace.horizontal = function(xLeft, xRight, yBottom, yTop, color=1, lwd=3, fine = 100, power = 9)
{
  xMid = (xLeft + xRight)/2
  xMM1 = (3*xLeft + xRight)/4;
  xMM2 = (xLeft + 3*xRight)/4;
  yMid = (yBottom + yTop)/2;

  xx = seq(from = -1, to=1, length.out = fine);
  yy = xx^power;

  x = c( xMM1 + (xMM1-xLeft) * xx, xMM2 + (xMM2 - xMid) * xx);
  y = yMid + (yMid-yBottom) * c(yy, rev(yy));

  lines(x, y, col = color, lwd = lwd);
} 

#=====================================================================================================
#
# Pretty names for meta-analysis output
#
#=====================================================================================================

prettyNamesForMetaAnalysis = function(names)
{
  translation = matrix( c(
        "Z.equalWeights",  "Stouffer (equal wts)",
        "Z.RootDoFWeights",  "Stouffer (sqrt wts)",
        "Z.DoFWeights", "Stouffer (dof wts)",
        "pValueHighScale.equalWeights", "Scale (equal wts)",
        "pValueHighScale.RootDoFWeights", "Scale (sqrt wts)",
        "pValueHighScale.DoFWeights", "Scale (dof wts)",
        "consensus", "Consensus",
        "weightedAverage.equalWeights", "Mean (equal wts)",
        "weightedAverage.RootDoFWeights", "Mean (sqrt wts)",
        "weightedAverage.DoFWeights", "Mean (dof wts)",
        "meta.Z.equalWeights", "Stouffer (equal wts)",
        "meta.Z.RootDoFWeights", "Stouffer (sqrt wts)",
        "meta.Z.DoFWeights", "Stouffer (dof wts)",
        "pValueHighScale.equalWeights", "Scale (equal wts)",
        "pValueHighScale.RootDoFWeights", "Scale (sqrt wts)",
        "pValueHighScale.DoFWeights", "Scale (dof wts)",
        "pValueHigh.equalWeights", "Rank (equal wts)",
        "pValueHigh.RootDoFWeights", "Rank (sqrt wts)",
        "pValueHigh.DoFWeights", "Rank (dof wts)"
        ), ncol = 2, byrow = TRUE)

  names2trans = match(names, translation[, 1]);
  out = names;
  out[is.finite(names2trans)] = translation[ names2trans[is.finite(names2trans)], 2];
  out;
}

#=====================================================================================================
#
# prependZeros
#
#=====================================================================================================
# prepend as many zeros as necessary to fill number to a certain width. Assumes an integer input.

########### This function is now in WGCNA package.


#=====================================================================================================
#
# translateUsingTable. Now capable of handling nested lists.
#
#=====================================================================================================

translateUsingTable = function(x, translationTable, keepUntranslated = FALSE)
{
  if (is.atomic(x))
  {
    out = translationTable[ match(x, translationTable[, 1]), 2]
    if (keepUntranslated) out[is.na(out)] = x[is.na(out)];
    out;
  } else lapply(x, translateUsingTable, translationTable, keepUntranslated = keepUntranslated);
}

translationTable = function(from, to)
{
  if (length(from)!= length(to)) stop("Length of 'from' and 'to' must be the same.");

  tab = as.matrix(table(from, to));
  if (ncol(tab)!=nrow(tab)) {
     printFlush("Error in translationTable: table is not 1-to-1.")
     if (ncol(tab)<10 && nrow(tab) < 10) print(tab);
     stop();
  }

  nonZeros.row = rowSums( tab!=0 ) 
  nonZeros.col = colSums( tab!=0 )

  if (any(c(nonZeros.row, nonZeros.col) != 1)) {  
     printFlush("Error in translationTable: table is not 1-to-1.")
     if (ncol(tab)<10 && nrow(tab) < 10) print(tab);
     stop();
  }

  unique.from = sort(unique(from));
  tt = cbind(unique.from, to [ match(unique.from, from)] );
  colnames(tt) = c("from", "to");
  tt;
}
  
#======================================================================================================
#
# Convert a multi-level factor variable to a set of binary variables
#
#======================================================================================================

multiLev2numeric = function(x, nameBase, nameSep = ".", minFraction = 0.05, exclude = c("", "NA"))
{
  t = table(x);
  minn = ceil(length(x) * minFraction);
  keepLevels = names(t)[t >= minn & ! (names(t) %in% exclude)];
  if (sum(t<minn) > 0) addOthers = TRUE else addOthers = FALSE;
  nLevels = length(keepLevels)
  nVars = nLevels;

  numer = matrix(NA, length(x), nVars + addOthers)

  for (l in 1:nLevels)
    numer[, l] = as.numeric(x==keepLevels[l]);

  if (addOthers) numer[, nVars + 1] = as.numeric(!(x%in%keepLevels));

  numer[x %in% exclude, ] = NA;

  if (addOthers) keepLevelsX = c(keepLevels, "other") else keepLevelsX = keepLevels;

  colnames(numer) = gsub("[ /]", ".", spaste(nameBase, nameSep, keepLevelsX));

  numer;
}


# This function attempts to determine whether a vector is numeric in the sense that coercing it to numeric
# will not lead to information loss.

isNumericVector = function(x, naStrings = c("NA", "NULL", "NO DATA"))
{
  if (is.numeric(x)) return(TRUE)

  x[x%in% naStrings] = NA
  x.num = suppressWarnings(as.numeric(x));
  missing = is.na(x.num);
  t = table(x[missing])
  if (length(t) ==0 ) return (TRUE)
  if (length(t)>1) return(FALSE)
  #if (all(missing)) return(TRUE) else return(FALSE);
  return(FALSE);
}

# Split a variable into several variables by levels of a covariate. The levels can be combined.

restrictVariableByCovariateLevels = function(var, covar, varName, covarName,
                                      covarLevels = sort(unique(covar)),
                                      covarLevelNames = covarLevels,
                                      nameSep = spaste(".in.", covarName, if (covarName=="") "" else "."),
                                      levelSep = ".", returnMatrix = FALSE)
{
  n = length(covarLevels);
  if (n<1) stop("No 'covarLevels' given.");
  allCovarLevels = sort(unique(covar));
  lapply(covarLevels, function(cl) if (!all(cl %in% allCovarLevels)) 
    warning(immediate. = TRUE, "Some entries in 'covarLevels' are not present in 'covar'."));
  out = list();
  for (cl in covarLevels)
  {
    vec = var;
    vec[!covar %in% cl] = NA;
    out = c(out, list(vec));
  }
  names(out) = spaste(varName, nameSep, sapply(covarLevels, paste, collapse = levelSep));
  if (returnMatrix) as.matrix(as.data.frame(out)) else as.data.frame(out);
}

restrictVariablesByCovariateLevels = function(x, covar, varNames = colnames(x), ...)
{
  if (length(dim(x))!=2) stop("'x' must be a two-dimensional object.");
  if (length(varNames)!=ncol(x)) stop("Length of 'varNames' must equal the number of columns of 'x'.");
  do.call(cbind, removeListNames(mymapply(restrictVariableByCovariateLevels,
            var = x, varName = varNames, MoreArgs = list(covar = covar, ...))));
}
         
# If center is FALSE, vectors consisting of 0 and 1s will be turned to -1 and 1
centerBinaryColumns = function(x, center = TRUE)
{
  if (mode(x)=="numeric")
  {
    out = apply(as.matrix(x), 2, function(x1)
    {
      x2 = x1[is.finite(x1)];
      if (length(x2)==0) return(x1);
      if (all(x2 %in% c(0,1))) 
      {
        out = x1-mean(x2) 
        if (!center) out = sign(out)
      } else out = x1;
      out;
    });
    dim(out) = dim(x);
    dimnames(out) = dimnames(x);
    out;
  } else if (mode(x)=="list") {
    out = lapply(x, function(x1)
    {
      if (!is.numeric(x1)) return(x1);
      x2 = x1[is.finite(x1)];
      if (length(x2)==0) return(x1);
      if (all(x2 %in% c(0,1))) 
      {
        out = x1-mean(x2) 
        if (!center) out = sign(out)
      } else out = x1;
      out;
    });
    if (is.data.frame(x)) out = data.frame(out, check.names = FALSE);
    out;
  }
}




#=======================================================================================================
#
# normalize
#
#=======================================================================================================

# Normalize a continuous variable to normal distribution

normalize = function(x, ties.method = "average" )
{
  rank = rank(x, ties.method = ties.method, na.last = "keep");
  qnorm( (rank-0.5) / max(rank, na.rm = TRUE));
}

#=======================================================================================================
#
# Convenience functions for Fisher test and kruskal.test
#
#=======================================================================================================

fisherTestPvalue = function(x, y, alternative = "two.sided", default = NA)
{
  keep = !is.na(x) & !is.na(y);
  tab = table(x[keep], y[keep]);
  tryCatch( fisher.test(x[keep], y[keep], alternative = alternative)$p.value, error = function(e) default);
}

matrixFisherTestPvalue = function(x, y, alternative = "two.sided")
{
  x = as.matrix(x);
  y = as.matrix(y);

  nx = ncol(x);
  ny = ncol(y);
  mat = matrix(NA, nx, ny);
  for (ix in 1:nx) for (iy in 1:ny)
    mat[ix, iy] = fisherTestPvalue(x[, ix], y[, iy]);
  colnames(mat) = colnames(y);
  rownames(mat) = colnames(x);
  finite = is.finite(mat);
  mat[finite][mat[finite] > 1] = 1;
  mat;
}

#=======================================================================================================
#
# multi-column labeled heatmap
#
#=======================================================================================================

# Multi-column display for tall, narrow labeled heatmaps.

labeledHeatmap.multiCol = function(Matrix, 
                           nColumns,
                           margins = NULL,
                           columnContent = NULL,
                           setLayout = TRUE,
                           widths = NULL,
                           xLabels = NULL, yLabels = NULL, 
                           xSymbols = NULL, ySymbols = NULL, 
                           textMatrix = NULL, 
                           colors.lab.x = 1,
                           colors.lab.y = 1,
                           plotLegend = TRUE, 
                           main = "",
                           cex.main = 1.4,
                           cex.lab = 1,
                           cex.lab.x = cex.lab,
                           cex.lab.y = cex.lab,
                           marginAdd = c(0.5, 0.5, 0.5, 0.5),
                           zlim = NULL,
                           colors = NULL,
                           mainHeight = 0.1,
                           titleMargin = 2,
                           parSettings = NULL,
                           ... ) 
{

  Matrix = as.matrix(Matrix);
  nRows = nrow(Matrix);
  if (nColumns > nRows) stop("Number of columns cannot be larger than number of rows of 'Matrix'.");
  if (is.null(columnContent))
    columnContent = allocateJobs(nRows, nColumns);

  if (is.null(widths)) widths = rep(1, nColumns);

  #if (setLayout)
    #layout( matrix( c( rep(1, nColumns), 2:(nColumns+1)), 2, nColumns, byrow = TRUE), 
    #        heights = c(mainHeight, 1-mainHeight), width = widths);

  par(oma = c(0, 0, titleMargin, 0));

  if (setLayout)
    layout( matrix( 1:nColumns, 1, nColumns, byrow = TRUE), width = widths);

  if (!is.null(parSettings)) do.call(par, parSettings);

  ind  = 0;
  colIndex = unlist(lapply(columnContent, function(x) { ind <<- ind + 1; return(rep(ind, length(x))) } ))

  if (!is.null(yLabels))
  {
    yLabelsX = tapply(yLabels, colIndex, I)
  } else 
    yLabelsX = rep(NULL, nColumns);

  if (!is.null(ySymbols))
  {
    ySymbolsX = tapply(ySymbols, colIndex, I)
  } else 
    ySymbolsX = rep(NULL, nColumns);

 
  if (is.null(margins))
  {
    stop("Margins must be specified.")
    # try to determine margins automatically. For now base the widths on yLabels but this is not correct in
    # general.
    # For now also assume that xLabels are positioned at the bottom and are rotated 45 degrees.
    # This code needs re-working since it assumes knowledge of user units which are not really available
    # yet.
    bottom = rep (1/sqrt(2) * (strheight(xLabels) + strwidth(xLabels)) + marginAdd[1] , nColumns);
    top = rep( strheight(main, cex = cex.main)/par("cxy") + marginAdd[3], nColumns);
    left = sapply( lapply(yLabelsX, strwidth, cex = cex.lab.y)/par("cxy")[1], max) + marginAdd[2];
    right = rep(marginAdd[4], nColumns);
    margins = as.list( as.data.frame( matrix ( c(bottom, left, top, right), 4, nColumns, byrow = TRUE)));
  }

  if (is.null(zlim))
  {
    max = max(abs(Matrix), na.rm = TRUE);
    min = min(Matrix, na.rm = TRUE);
    if (any(Matrix < 0)) zlim = c(-max, max) else zlim = c(min, max);
  }

  mainCol = as.integer((nColumns + 1)/2);

  #par(mar = c(0,0,0,0));
  #plot(c(0,1), c(0,1), type = "n", axes = FALSE, frame = FALSE, xlab = "", ylab = "")
  #text(0.5, 0.5, main, cex = cex.main);

  for (col in 1:nColumns)
  {
    par(mar = margins[[col]]);
    WGCNA::labeledHeatmap( Matrix[columnContent[[col]], , drop = FALSE],
                    xSymbols = xSymbols, 
                    xLabels = xLabels,
                    ySymbols = ySymbolsX[[col]],
                    yLabels = yLabelsX[[col]],
                    setStdMargins = FALSE,
                    main = "",
                    cex.lab = cex.lab,
                    cex.lab.x = cex.lab.x,
                    cex.lab.y = cex.lab.y,
                    zlim = zlim,
                    colors = colors,
                    textMatrix = if (!is.null(textMatrix)) 
                         textMatrix[ columnContent[[col]], , drop = FALSE] else NULL,
                    plotLegend = (col==nColumns),
                    ...);
    if (col==mainCol)
      title(main, cex.main = cex.main, outer = TRUE);
  }
}


corHeatmapWithDendro = function(corMat, 
                             maxMerge = 18,
                             setLayout = TRUE,
                             dendroWidth = 0.25,
                             mar.main = c(5, 7, 2, 1),
                             leftMar.dendro = 1,
                             roundTextToDigits = 2,
                             showTextAbove = 0.2,
                             xLabels = colnames(corMat),
                             yLabels = rownames(corMat),
                             xSymbols = NULL,
                             ySymbols = NULL,
                             colors = blueWhiteRed(100),
                             zlim = c(-max(abs(corMat), na.rm = TRUE), max(abs(corMat), na.rm = TRUE)),
                             tree = NULL,
                             colors.lab.x = 1, colors.lab.y = 1, 
                             replaceMissingBy = 0,
                             pValues = NULL,
                             showPValues = TRUE,
                             pValuesBelow = TRUE,
                             pThreshold = 0.05,
                             ...)
                             
{
  corMat.rep = corMat;
  corMat.rep[is.na(corMat)] = replaceMissingBy;
  if (is.null(tree))
    tree = bestTree(1-corMat.rep, maxMerge = maxMerge )
  order1 = tree$order;
  tree.rev = tree;
  tree.rev$order = rev(tree$order);

  n = length(order1) 
  if (length(colors.lab.x)==n) colors.lab.x = colors.lab.x[order1];
  if (length(colors.lab.y)==n) colors.lab.y = colors.lab.y[order1];
  if (setLayout) layout(matrix(c(1,2), 1, 2), widths= c(dendroWidth, 1-dendroWidth));
  marDendro = mar.main; marDendro[2] = leftMar.dendro; marDendro[4] = 0;
  par(mar = marDendro)
  plotDendrogram(tree.rev, labels = FALSE, hang = -1, horiz = TRUE, reverseDirection = FALSE, 
                 adjustRange = TRUE);
  par(mar = mar.main);
  textMat0 =  round(corMat[order1, order1], roundTextToDigits);
  textMat = textMat0;
  if (!is.null(pValues) && showPValues)
  {
    if (pValuesBelow) textMat = spaste(textMat, "\n", signif(pValues[order1, order1], 1)) else
        textMat = spaste(textMat, " (", signif(pValues[order1, order1], 1), ")");
    par(lheight = 0.85);
  }
  if (is.null(pValues))
  {
     textMat[abs(textMat0) < showTextAbove] = "";
  } else
     textMat[pValues[order1, order1] > pThreshold] = "";
  labeledHeatmap(corMat[order1, order1],
               xLabels = xLabels[order1],
               yLabels = yLabels[order1],
               xSymbols = if (is.null(xSymbols)) NULL else xSymbols[order1],
               ySymbols = if (is.null(ySymbols)) NULL else ySymbols[order1],
               textMatrix = textMat,
               colors = colors,
               zlim = zlim,
               setStdMargins = FALSE, 
               colors.lab.x = colors.lab.x, colors.lab.y = colors.lab.y, ...);
  invisible(list(tree=tree, tree.rev = tree.rev))
}


heatmapWithDendro = function(
  mat,
  plotRowTree = TRUE,
  plotColTree = TRUE,
  rowTree = NULL,
  colTree = NULL,
  maxMerge = 18,
  setLayout = TRUE,
  dendroWidth = 0.25,
  dendroHeight = 0.25,
  mar.main = c(5, 7, 2, 1),
  leftMar.dendro = 1,
  topMar.dendro = 2,
  textMat = NULL,
  showText =TRUE,
  roundTextToDigits = 2,
  showTextAbove = 0.2,
  xLabels = colnames(mat),
  yLabels = rownames(mat),
  xSymbols = NULL,
  ySymbols = NULL,
  colors = blueWhiteRed(100),
  zlim = c(-max(abs(mat), 1, na.rm = TRUE), max(abs(mat), 1, na.rm = TRUE)),
  colors.lab.x = 1, colors.lab.y = 1, 
  colTree.shift = 0.5,
  main = "", ...)

{
  nr = nrow(mat);
  nc = ncol(mat);
  if (plotRowTree)
  {
    if (is.null(rowTree))
      rowTree = bestTree(as.matrix(dist(mat)), maxMerge = maxMerge )
  }
  if (plotColTree)
  {  
     if (is.null(colTree))
       colTree = bestTree(as.matrix(dist(t(mat))), maxMerge = maxMerge )
  }
  if (!is.null(rowTree)) rowOrder1 = rowTree$order else rowOrder1 = 1:nr
  if (!is.null(colTree)) colOrder1 = colTree$order else colOrder1 = 1:nc
  rowTree.rev = rowTree;
  rowTree.rev$order = rev(rowTree$order);

  if (length(colors.lab.x)==nc) colors.lab.x = colors.lab.x[colOrder1];
  if (length(colors.lab.y)==nr) colors.lab.y = colors.lab.y[rowOrder1];
  if (setLayout)
    switch( plotRowTree  + 2*plotColTree + 1, 
       layout(matrix(1, 1, 1)),
       layout(matrix(c(1, 2), 1, 2), widths= c(dendroWidth, 1-dendroWidth)),
       layout(matrix(c(1, 2), 2, 1), heights = c(dendroHeight, 1-dendroHeight)),
       layout(matrix(c(0,1,2,3), 2, 2), widths= c(dendroWidth, 1-dendroWidth), heights = c(dendroHeight, 1-dendroHeight)));

  #if (plotColTree) mar.main[3] = 0;
  marRowDendro = mar.main; marRowDendro[2] = leftMar.dendro; marRowDendro[4] = 0;
  marColDendro = mar.main; marColDendro[3] = topMar.dendro; marColDendro[1] = 0;
  if (plotRowTree)
  {
    par(mar = marRowDendro)
    plotDendrogram(rowTree.rev, labels = FALSE, hang = -1, horiz = TRUE, reverseDirection = FALSE,
                   adjustRange = TRUE);
  }
  if (plotColTree)
  {
    par(mar = marColDendro);
    plotDendrogram(colTree, labels = FALSE, hang = -1, horiz = FALSE, reverseDirection = FALSE,
                   adjustRange = TRUE, shift.last = colTree.shift, main = main);
  }
  par(mar = mar.main);
  if (showText)
  {
    if (is.null(textMat))
    {
      textMat =  round(mat, roundTextToDigits);
      textMat[abs(textMat) < showTextAbove] = "";
    } 
  }

  if (!is.null(textMat))
  {
    if (is.null(dim(textMat))) dim(textMat) = dim(mat);
    if (!isTRUE(all.equal(dim(textMat), dim(mat)))) stop("'textMat' and 'mat' must have the same dimensions."); 
  }

  labeledHeatmap(mat[rowOrder1, colOrder1, drop = FALSE],
               xLabels = xLabels[colOrder1],
               yLabels = yLabels[rowOrder1],
               xSymbols = xSymbols[colOrder1],
               ySymbols = ySymbols[rowOrder1], 
               textMatrix = textMat[rowOrder1, colOrder1, drop = FALSE],
               colors = colors,
               #showRows = rowOrder1,
               #showCols = colOrder1,
               zlim = zlim,
               setStdMargins = FALSE,
               colors.lab.x = colors.lab.x, colors.lab.y = colors.lab.y, 
               main = if (plotColTree) "" else main, ...);
  invisible(list(rowTree=rowTree, rowTree.rev = rowTree.rev, colTree = colTree))
}



#========================================================================================================
#
# Fast matrix Student and Kruskal tests and basic statistics for standard screening for binary/categorical
# traits
#
#========================================================================================================

# This function computes basic statistics (mean, sum, sum of squares, std. dev., number of present
# observations) for each column in x stratified by factor g. It is roughly equivalent to
# apply(x, tapply, g, {sum, sum^2, mean, nPresent, stdDev}, na.rm = TRUE) but should be faster since it
# explicitly loops only over the levels of g, not over the columns of x.

# The group argument 'g' can be basically anything (character, numeric, logical, etc)

# The function also works for a single level g (the default), in which case it can simplify the matrix
# outputs to vectors. 

colStatsByGroup = function(x, g = NULL, assignNames = TRUE, getSquares = TRUE, simplify = TRUE)
{
  x = as.matrix(x)
  if (is.null(g)) g = rep(1, nrow(x));
  gLevels = sort(unique(g));
  nLevels = length(gLevels)
  nCols = ncol(x);
  groupSum = groupNPresent = matrix(0, nLevels, nCols);
  if (assignNames) {
    dimnames(groupSum) = dimnames(groupNPresent) = list ( gLevels, colnames(x));
  }
  if (getSquares) groupSumSq = groupSum;
  x.nna = x;
  nax = is.na(x);
  finx = !nax
  x.nna[nax] = 0;
  for (gl in 1:nLevels)
  {
    indicator = as.numeric(!is.na(g) & g==gLevels[gl]);
    groupSum[gl, ] = t(indicator) %*% x.nna;
    if (getSquares) groupSumSq[gl, ] = t(indicator) %*% x.nna^2;
    groupNPresent[gl, ] = t(indicator) %*% finx;
  }
  groupMean = groupSum/groupNPresent;
  if (getSquares) {
    num = groupSumSq - groupMean^2 * groupNPresent;
    num[num < 0] = 0;
    groupSD = sqrt( num / (groupNPresent - 1) );
    #if (any(!is.finite(groupSD))) browser();
    groupSE = groupSD/sqrt(groupNPresent);
  } else {
    groupSumSq = groupSD = groupSE = NULL;
  }
  list(sum = groupSum[, , drop = simplify], sumSq = groupSumSq[, , drop = simplify], 
       mean = groupMean[, , drop = simplify], stdDev = groupSD[, , drop = simplify], 
       stdErr = groupSE[, , drop = simplify],
       nPresent = groupNPresent[, , drop = simplify]);
}

colMeansByGroup = function(x, g = NULL, w = NULL, assignNames = TRUE, simplify = TRUE)
{
  x = as.matrix(x);
  if (is.null(w)) w = matrix(1, nrow(x), ncol(x)) else w = as.matrix(w);
  if (is.null(g)) g = rep(1, nrow(x));
  g = as.factor(g);
  gLevels = levels(g);
  nLevels = length(gLevels)
  nCols = ncol(x);
  groupSum = groupSumW = matrix(0, nLevels, nCols);
  if (assignNames) {
    dimnames(groupSum) = list ( gLevels, colnames(x));
  }
  w.nna = w;
  x.nna = x;
  nax = is.na(x) | is.na(w);
  finx = !nax
  w.nna[nax] = 0;
  x.nna[nax] = 0;
  wx = w.nna * x.nna;
  for (gl in 1:nLevels)
  {
    indicator = as.numeric(!is.na(g) & g==gLevels[gl]);
    groupSum[gl, ] = t(indicator) %*% wx;
    groupSumW[gl, ] = t(indicator) %*% w.nna;
  }
  groupMean = groupSum/groupSumW;
  groupMean;
}

if (FALSE)
{
  mat = matrix(1:100, 10, 10);
  w = matrix(1, 10, 10);
  w[1, ] = 0.5;

  grp = rep(c(1, 2), each = 5)

  colMeansByGroup(mat, grp, w)
}


# Fast Student single-sample t-test. Assumes that the null-mean is zero.

matrixSingleSampleTTest = function(x, alternative, colStats = NULL)
{
  if (is.null(colStats))
    colStats = colStatsByGroup(x, assignNames = FALSE, getSquares = TRUE, simplify = TRUE);

  t = colStats$mean/colStats$stdDev * sqrt(colStats$nPresent);
  df = colStats$nPresent - 1;

  if (alternative=="two.sided")
  {
    p = 2*pt(abs(t), df = df, lower.tail = FALSE);
    Z = -sign(t)*qnorm(p/2);
  } else {
    p = pt(t, df = df, lower.tail = alternative=="less");
    Z = qnorm(p, lower.tail = alternative=="less")
  }

  cbind(statistic = t, parameter = df, p.value = p, Z = Z)
}

# Matrix binomial test. Tests whether entries ineach column in x are equally likely to be positive and
# negative. Can be viewed as a single-sample analog to Kruskal-Wallis or a non-parametric version of
# single-sample t-test.

matrixSingleSampleBinomialTest = function(x, alternative)
{
  nLarger = colSums(x > 0, na.rm = TRUE)
  nSmaller = colSums(x < 0, na.rm = TRUE)

  nAll = nLarger + nSmaller 

  if (alternative=="two.sided")
  {
    x = pmin(nLarger, nSmaller);
  } else if (alternative=="less") {
    x = nLarger
  } else if (alternative=="greater") {
    x = nSmaller
  } else stop(spaste("Error in matrixSingleSampleBinomialTest: Unrecognized 'alternative' ", alternative));

  p = pbinom(x, nAll, prob = 0.5, lower.tail = TRUE);
  Z = -qnorm(p);
  if (alternative=="two.sided") {
    p = 2*p;
    p[p>1] = 1;
  }
  cbind(p.value = p, nLess = nSmaller, nGreater = nLarger, Z = Z);
}


# Fast Student t test for a single indicator of group membership. Assumes g has two levels. If column stats
# have already been calculated, they can be supplied and it will speed up the calculation further.
# This function only handles the unpaired test.

# here g can be any two-level variable.

# Note: the sign of the t statistic follows R convention and is positive if x[level 1] is greater  than
# x[level 2]

matrixTTest.singleG = function(x, g, 
                               alternative,
                               reverseGroups = FALSE,
                               var.equal = FALSE,
                               colStats = NULL)
{
  gFinite = is.finite(g);
  nGLevels = length(unique(g[gFinite]));
  nCols = ncol(x);
  if (!any(gFinite) || nGLevels < 2)
  {
    out = matrix(NA, nCols, 4);
    colnames(out) = TTestStatNames();
    return(out);
  }

  if (is.null(colStats))
    colStats = colStatsByGroup(x, g, assignNames = FALSE, getSquares = TRUE, simplify = FALSE);

  groupMean = colStats$mean;
  nPresent = colStats$nPresent;
  groupSD = colStats$stdDev;
  if (var.equal)
  {
     df = colSums(nPresent) - 2;
     denom = sqrt( colSums(1/nPresent) * colSums( (nPresent -1) * groupSD^2 ) / df );
  } else {
     meanVar.1 = groupSD^2 / nPresent
     meanVarSums = colSums(meanVar.1);
     denom = sqrt(meanVarSums);
     df = meanVarSums^2 / colSums( meanVar.1^2/(nPresent -1 ))
  }
  denom[!is.finite(denom) | denom < 0] = NA;
  if (reverseGroups)
  {
    g.num = 2;
    g.denom = 1;
  } else {
    g.num = 1;
    g.denom = 2;
  }

  t = (groupMean[g.num, ] - groupMean[g.denom, ])/denom;
  t[!is.finite(t)] = NA;

  if (alternative=="two.sided")
  {
    p = 2*pt(abs(t), df = df, lower.tail = FALSE);
    Z = -sign(t)*qnorm(p/2);
  } else {
    p = pt(t, df = df, lower.tail = alternative=="less");
    Z = qnorm(p, lower.tail = alternative=="less")
  } 

  cbind(statistic = t, parameter = df, p.value = p, Z = Z);
}

TTestStatNames = function()
{
  c("statistic", "parameter", "p.value", "Z");
}

# Fast matrix Kruskal test. The test is "matricized" in x. This workhorse function assumes that g is a
# vector.

# Convention: the sign is positive when the spearman correlation of x and g is positive.
# Note that this is the opposite of the sign convention for the t-test which follows the R convention and
# is positive if group 1 is higher than group 2 (negative correlation).

# Note: because of the use of is.finite, this function requires that the group argument 'g' be numeric. A
# factor derived from numbers may work as well.

# uses function rowRanks from matrixStats

matrixKruskalTest.singleG = function(x, g, getSign, getAreaUnderROC)
{
  gFinite = is.finite(g);
  nGLevels = length(unique(g[gFinite]));
  nCols = ncol(x);
  if (!any(gFinite) || nGLevels < 2)
  {
    out = matrix(NA, nCols, 4 + getSign + getAreaUnderROC);
    colnames(out) = KTstatNames(getSign, getAreaUnderROC)
    return(out);
  }

  if (any(!gFinite))
  {
    x = x[gFinite, , drop = FALSE];
    g = g[gFinite];
  }

  nSamples = nrow(x);

  r = colRanks(as.matrix(x), ties.method = "average", preserveShape = TRUE);
  # r = as.matrix(apply(x, 2, rank, na.last = "keep", ties.method = "average"));

  meanR = colMeans(r, na.rm = TRUE);

  gLevels = sort(unique(g));
  nLevels = length(gLevels)
  if (nLevels > 2 && (getSign | getAreaUnderROC))
    stop("Sign can only be determined for 2-group comparisons.");

  nCols = ncol(x);

  groupStats = colStatsByGroup(r, g, assignNames = FALSE, getSquares = FALSE, simplify = FALSE);
  groupNPresent = groupStats$nPresent;
  groupMeans = groupStats$mean;
  groupMeans[!is.finite(groupMeans)] = 0;
  nPresent = colSums(groupNPresent)
  nLevels.byCol = colSums(groupNPresent > 0);

  num = colSums((groupMeans - matrix(meanR, nrow = nLevels, ncol = nCols, byrow = TRUE))^2 * 
                 groupNPresent, na.rm = TRUE);
  denom = colSums( (r - matrix(meanR, nSamples, nCols, byrow = TRUE))^2, na.rm = TRUE)

  K = (nPresent -1) * num/denom;
  K[!is.finite(K)] = NA;
  df = nLevels.byCol - 1;
  p = pchisq(K, df = df, lower.tail = FALSE);
  sign = 2 * ( as.numeric(groupMeans[1, ] < groupMeans[2, ]) - 0.5 )
  Z = -qnorm(p/2) * sign;
  out = cbind(statistic = K, parameter = df, p.value = p, Z = Z);

  if (!is.null(colnames(x))) rownames(out) = colnames(x);

  if (getSign)
     out = cbind(out, sign = sign);

  if (getAreaUnderROC)
  {
    S = groupStats$sum[1, ] - groupNPresent[1, ] * (groupNPresent[1, ] + 1) / 2; 
    auc = S/(groupNPresent[1, ] * groupNPresent[2, ]); 
    # Note: the sign of auc-0.5 should be consistent with "sign" above in this convention.
    areaUnderROC = auc # ifelse(auc > 0.5, auc, 1-auc);
    areaUnderROC[!is.finite(areaUnderROC)] = NA;
    out = cbind(out, areaUnderROC = areaUnderROC);
  }
  out;
}
  
# Can also take a matrix for g but the test is not matricized in g because 
# missing data cannot be handled correctly for all columns of g at the same time (rankings would be
# incorrect). 

matrixTest = function(x, g, testFnc, statNames, ...)
{
  testFnc = match.fun(testFnc);
  d = dim(g);
  x = as.matrix(x);
  nCols = ncol(x);
  nSamples.g = if (is.null(d)) length(g) else d[1];
  if (nSamples.g !=nrow(x)) 
    stop ("Number of obervations in 'g' must equal the number of observations in 'x'.");
  nStats = length(statNames);
  # Check if there are multiple columns in g
  if (length(d) > 1)
  {
     #For multiple columns, call self on each column and assemble output
     out = array(NA, dim = c(nCols, nStats, d[2]));
     dimnames(out) = list(colnames(x), statNames, colnames(g));
     for (gc in 1:d[2])
       out[, , gc] = testFnc(x, as.numeric(factor(g[, gc])), ...);
  } else {
     out = as.data.frame(testFnc(x, g, ...));
     rownames(out) = colnames(x);
     colnames(out) = statNames;
  }
  out;
}

KTstatNames = function(getSign, getAreaUnderROC)
{
  statNames = c("statistic", "parameter", "p.value", "Z");
  if (getSign) statNames = c(statNames, "sign");
  if (getAreaUnderROC) statNames = c(statNames, "areaUnderROC");
  statNames;
}

matrixTTest = function(x, g,      
                       reverseGroups = FALSE,
                       alternative = c("two.sided", "less", "greater"),
                       var.equal = FALSE)
{
  alternative = match.arg(alternative)
  matrixTest(x, g, testFnc = matrixTTest.singleG, statNames = TTestStatNames(),
             reverseGroups = reverseGroups, alternative = alternative, var.equal = var.equal);
}

matrixKruskalTest = function(x, g, getSign = FALSE, getAreaUnderROC = FALSE)
{
  matrixTest(x, g, testFnc = matrixKruskalTest.singleG, 
             statNames = KTstatNames(getSign, getAreaUnderROC), 
             getSign = getSign, getAreaUnderROC = getAreaUnderROC);
}

# Area under ROC. 
areaUnderROC = function(x, g)
{
  levels = sort(unique(g));
  ng = table(g);
  r = rank(x, na.last = "keep");
  S = sum(r[g==levels[1] ]) - ng[1] * (ng[1] + 1)/2;
  max(1 - S/prod(ng));
}



#==================================================================================================
#
# standardScreening family of functions re-written to be faster, more flexible and require less maintenance.
#
#==================================================================================================

# A re-write of Steve's function that runs much faster, including AUC calculations.

nSSBTstats = function(getQvalues, getFDR, getZstatistics, getAreaUnderROC, ...)
{
  14 + 2*(as.numeric(getQvalues) + as.numeric(getFDR) + as.numeric(getZstatistics)) + getAreaUnderROC;
}

SSBTstatNames = function(getQvalues, getFDR, getZstatistics, getAreaUnderROC, corName, levelNames = NULL, 
                         dataIsLogTransformed, reportLogFoldChange, ...)
{
  if (dataIsLogTransformed)
  {
    logMod = "Log"
    SElogMod = ".Log"
  } else {
    logMod = ""
    SElogMod = ""
  }
  if (is.null(levelNames)) 
    levelNames = spaste("Group.", 1:2);
  out = c(corName, "t.Student", "p.Student");
  if (getZstatistics) out = c(out, "Z.Student");
  if (getQvalues) out = c(out, "q.Student");
  if (getFDR) out = c(out, "FDR.Student");
  out = c(out, if (reportLogFoldChange) "log2FoldChange" else "foldChange", 
          spaste("mean", logMod, ".group.", levelNames), 
          spaste("stdErr", SElogMod, ".group.", levelNames),
          spaste("nPresent.group.", levelNames), "nPresent.bothGroups",
          "stat.Kruskal", "stat.Kruskal.signed", "p.Kruskal");
  if (getZstatistics) out = c(out, "Z.Kruskal");
  if (getQvalues) out = c(out, "q.Kruskal");
  if (getFDR) out = c(out, "FDR.Kruskal");
  if (getAreaUnderROC) out = c(out, "areaUnderROC");
  out;
}

SSSstatNames = function(getQvalues, getFDR, getZstatistics, dataIsLogTransformed, ...)
{
  if (dataIsLogTransformed)
  {
    logMod = "Log"
    SElogMod = ".Log"
  } else {
    logMod = ""
    SElogMod = ""
  }

  out = c("t.Student", "p.Student");
  if (getZstatistics) out = c(out, "Z.Student");
  if (getQvalues) out = c(out, "q.Student");
  if (getFDR) out = c(out, "FDR.Student");
  out = c(out, spaste("mean", logMod), spaste("stdErr", SElogMod), "nPresent", "p.Binomial");
  if (getZstatistics) out = c(out, "Z.Binomial");
  if (getQvalues) out = c(out, "q.Binomial");
  if (getFDR) out = c(out, "FDR.Binomial");
  c(out, "nNegative", "nPositive");
}


nSSSstats = function(...)
{
  length(SSSstatNames(...));
}

nSSNTstats = function(getQvalues, getFDR, getAreaUnderROC, getDetails, ...)
{
  4 + getQvalues + getFDR + getAreaUnderROC + 2*getDetails;
}

SSNTstatNames = function(getQvalues, getFDR, corName, getAreaUnderROC, getDetails, getLMCoeffs, ...)
{
  out = c(corName, "Z", "p.Student");
  if (getQvalues) out = c(out, "q.Student");
  if (getFDR) out = c(out, "FDR.Student");
  if (getAreaUnderROC) out = c(out, "areaUnderROC");
  out = c(out, "nPresent");

  if (getDetails) out = c(out, "mean", "stdDev");
  if (getLMCoeffs) out = c(out, "coefficient")

  out;
}

# Workhorse function for binary traits
SSBT.singleTrait = function(data, trait, 
           corFnc = cor, corOptions = list(use = 'p'), corName = "cor",
           tTestAlternative,
           var.equal, 
           getQvalues,
           getFDR,
           getZstatistics,
           getAreaUnderROC,
           dataIsLogTransformed,
           logBase,
           logShift,
           levelNames = NULL,
           traitName = "",
           consistentSigns = TRUE,
           reportLogFoldChange = FALSE )
{
  f = factor(trait);
  trait = as.numeric(f);
  levels = levels(f);
  if (is.null(levelNames)) levelNames = levels;

  # Fudge levelNames for the case where trait is constant. 
  # I want to return missing data for this case and need levelNames of length 2.
  if (length(levelNames)<2) levelNames = c(0,1);

  nLevels = length(levels);
  if (nLevels > 2)
    stop("The sample trait contains more than 2 levels. Please input a binary variable.")
        
  if (length(trait)!=nrow(data))
    stop("Length of the sample trait does not equal the number of observations (rows) in 'data'.")

  nVars = ncol(data);

  out = matrix(NA, nVars, nSSBTstats(getQvalues, getFDR, getZstatistics, getAreaUnderROC = getAreaUnderROC));
  colnames(out) = SSBTstatNames(getQvalues, getFDR, getZstatistics, getAreaUnderROC = getAreaUnderROC, 
                                corName, levelNames, dataIsLogTransformed, reportLogFoldChange);

  if (nLevels<2 )
  {
    printFlush(spaste("Warning in standardScreeningBinaryTrait.singleTrait:\n   the trait ", 
                      traitName, if (nLevels==0) " contains only missing values." else " is constant."));
    return(out);
  }

  corFnc = match.fun(corFnc);
  corOptions$y = trait;
  corOptions$x = data;
  corPearson = as.numeric(do.call(corFnc, corOptions));

  sign = if (consistentSigns) -1 else 1;
  #print(paste("Consistent signs: ", consistentSigns));
  colStats = colStatsByGroup(data, trait, assignNames = FALSE, getSquares = TRUE);
  tTest = matrixTTest.singleG(data, trait, 
                              reverseGroups = consistentSigns,
                              var.equal = var.equal, alternative = tTestAlternative,
                              colStats = colStats)

  out[, c(1:3)] = cbind(cor = corPearson,
                   t.Student = tTest[, "statistic"],
                   #parameter.Student = tTest[, "parameter"],
                   p.Student = tTest[, "p.value"]);

  ind = 4;
  if (getZstatistics)
  {
    out[, ind] = tTest[, "Z"];
    ind = ind + 1;
  }

  if (getQvalues)
  {
    out[, ind] = qvalue.restricted(tTest[, "p.value"]);
    ind = ind + 1;
  }

  if (getFDR)
  {
    out[, ind] = p.adjust(tTest[, "p.value"], method = "fdr");
    ind = ind + 1;
  }

  # Basic statistics
  mean1 = colStats$mean[1, ];
  mean2 = colStats$mean[2, ];

  if (dataIsLogTransformed)
  {
    mean1.t = logBase^mean1-logShift;
    mean2.t = logBase^mean2-logShift;
  } else {
    mean1.t = mean1;
    mean2.t = mean2;
  }
  if (reportLogFoldChange)
  {
    if (dataIsLogTransformed)
    {
      fc =  mean1-mean2;
    } else
      fc = log2(mean1/mean2);
  } else
    fc = ifelse(mean1.t > mean2.t, mean1.t/mean2.t, -mean2.t/mean1.t);

  if (consistentSigns) fc = -fc;
     
  out[, ind:(ind + 7)] = cbind( 
               fc,
               cbind(mean1, mean2), 
               t(colStats$stdDev/sqrt(colStats$nPresent)),
               t(colStats$nPresent), colSums(colStats$nPresent));
  ind = ind + 8;

  # Kruskal test (including area under ROC). Note that matrixKruskalTest.singleG uses a convention opposite
  # to that of t-test for the sign, so here trait needs to be multiplied by -sign.

  kruskalTest = matrixKruskalTest.singleG(data, -sign * trait, getSign = TRUE, 
                                          getAreaUnderROC = getAreaUnderROC);
  out[, ind:(ind + 2)] = cbind(kruskalTest[, "statistic"],
                               kruskalTest[, "statistic"] * kruskalTest[, "sign"],
                               kruskalTest[, "p.value"]);

  ind = ind + 3;
  if (getZstatistics)
  {
    out[, ind] = kruskalTest[, "Z"]
    ind = ind + 1;
  }

  if (getQvalues)
  {
      out[, ind] = qvalue.restricted(kruskalTest[, "p.value"]);
      ind = ind + 1;
  }

  if (getFDR)
  {
      out[, ind] = p.adjust(kruskalTest[, "p.value"], method = "fdr");
      ind = ind + 1;
  }

  if (getAreaUnderROC) out[, ind] = 1-kruskalTest[, "areaUnderROC"];

  out;
}

# Workhorse function for single-sample screening (i.e., single-sample t-test and "Kruskal" test, which is
# really a binomial test). 

singleSampleScreening = function(data,
           alternative,
           getQvalues,
           getFDR,
           getZstatistics,
           dataIsLogTransformed,
           ...)
{
  nVars = ncol(data);

  out = matrix(NA, nVars, nSSSstats(getQvalues, dataIsLogTransformed));
  colnames(out) = SSSstatNames(getQvalues, dataIsLogTransformed);

  nSamples = nrow(data);
  colStats = colStatsByGroup(data, assignNames = FALSE, getSquares = TRUE);
  tTest = matrixSingleSampleTTest(data, alternative = alternative,
                              colStats = colStats)

  out[, c(1:2)] = cbind(
                   t.Student = tTest[, "statistic"],
                   p.Student = tTest[, "p.value"]);

  ind = 3;

  if (getZstatistics)
  {
    out[, ind] = tTest[, "Z"];
    ind = ind + 1;
  }

  if (getQvalues)
  {
    out[, ind] = qvalue.restricted(tTest[, "p.value"]);
    ind = ind + 1;
  }

  out[, ind:(ind + 2)] = cbind(colStats$mean, 
                               colStats$stdDev/sqrt(colStats$nPresent),
                               colStats$nPresent);
  ind = ind + 3;

  # Binomial test for whether positive and negative entries are equally likely
  bt = matrixSingleSampleBinomialTest(data, alternative)

  out[, ind] = bt[, "p.value"];
  ind = ind + 1;
  if (getZstatistics)
  {
    out[, ind] = bt[, "Z"];
    ind = ind + 1;
  }
  if (getQvalues)
  {
    out[, ind] = qvalue.restricted(bt[, "p.value"]);
    ind = ind + 1;
  }
  if (getFDR)
  {
    out[, ind] = p.adjust(bt[, "p.value"], mthod = "fdr");
    ind = ind + 1;
  }
  out[, ind:(ind+1)] = bt[, c("nLess", "nGreater")]
  out;
}

# Workhorse screening function for numeric traits
SSNT.singleTrait = function(data, trait,
                            corFnc,
                            corOptions,
                            corName,
                            alternative,
                            getQvalues,
                            getFDR,
                            getAreaUnderROC,
                            getDetails,
                            getLMCoeffs,
                            traitName = "",
                            exactZ = TRUE)
{
  
  corFnc = match.fun(corFnc);
  corOptions$y = trait;
  corOptions$x = data;
  cp = do.call(corFnc, corOptions);
  corPearson = as.numeric(cp);

  nVars = ncol(data);

  finMat = !is.na(data)
  nPresent = as.numeric( t(finMat) %*% (!is.na(as.matrix(trait))) );

  T = sqrt(nPresent - 2) * corPearson/sqrt(1 - corPearson^2)
  T[!is.finite(T)] = NA;
  if (alternative == "two.sided") {
       p = 2 * pt(abs(T), nPresent - 2, lower.tail = FALSE)
       Zp = sign(T) * qnorm(p/2, lower.tail = FALSE);
  }
  else if (alternative == "less") {
      p = pt(T, nPresent - 2, lower.tail = TRUE)
      Zp = sign(T) * qnorm(p, lower.tail = FALSE);
  }
  else if (alternative == "greater") {
      p = pt(T, nPresent - 2, lower.tail = FALSE)
      Zp = sign(T) * qnorm(p, lower.tail = FALSE);
  }
  pvalueStudent = as.numeric(p);

  if (exactZ)
  {
    Z = Zp
  } else 
    Z = 0.5 * log( (1+corPearson)/(1-corPearson) ) * sqrt(nPresent -2 );
  Z[!is.finite(Z)] = NA;
  if (getAreaUnderROC) 
  {
    AreaUnderROC = rep(NA, nVars);
    for (i in 1:ncol(data))
      AreaUnderROC[i] = rcorr.cens(data[, i], trait, outx = T)[[1]]
  }

  if (getQvalues)
    q.Student = qvalue.restricted(pvalueStudent);

  if (getFDR)
    fdr.Student = p.adjust(pvalueStudent, method = "fdr");

  if (getDetails)
  {
    means = colMeans(data, na.rm = TRUE);
    sd = colSds(data, na.rm = TRUE);
    details = cbind(mean = means, stdDev = sd);
  }

  output = cbind(corPearson,
                 Z,
                 if (getQvalues) cbind(pvalueStudent, q.Student) else pvalueStudent,
                 if (getFDR) fdr.Student else NULL,
                 if (getAreaUnderROC) cbind(AreaUnderROC, nPresent) else nPresent);
  if (getDetails) 
    output = cbind(output, details);

  if (getLMCoeffs)
  {
    keepSamples = is.finite(trait);
    if (sum(keepSamples) >= 2)
    {
      coeffList = linearModelCoefficients(data[keepSamples, ,drop = FALSE], trait[keepSamples]);
      coeffs = as.numeric(coeffList$coefficients[2, ]);
    } else coeffs = rep(NA, nVars);
    output = cbind(output, coeffs);
  }
    
  colnames(output) = SSNTstatNames(getQvalues = getQvalues, getFDR = getFDR, corName = corName, 
                                   getAreaUnderROC = getAreaUnderROC, getDetails = getDetails,
                                   getLMCoeffs);
  output
}

# Generic screening function: the workhorse for multi-trait screening

standardScreening.generic = function(data, trait,
           screeningFnc.singleTrait,
           nStatsFnc,
           statNameFnc,
           ...,
           verbose = 1,
           indent = 0,
           callFncOneByOne = TRUE,
           doCollectGarbage = FALSE,
           flatOutput = FALSE,
           addIDcolumn = TRUE,
           traitNameSep = ".for.",
           traitNames = NULL,
           varNames = NULL,
           additionalVarInfo = NULL)
{

  spaces = indentSpaces(indent);
 
  screeningFnc.singleTrait = match.fun(screeningFnc.singleTrait);
  nStatsFnc = match.fun(nStatsFnc);
  statNameFnc = match.fun(statNameFnc);

  nVars = ncol(data);
  if (is.null(varNames)) varNames = colnames(data);
  if (is.null(varNames)) varNames = spaste("Variable.", prependZeros(1:nVars))

  ID = data.frame(ID = varNames);
  if (!is.null(additionalVarInfo))
  {
    d = dim(additionalVarInfo)
    if (length(d) < 2) stop("If given, 'additionalVarInfo' must be a 2-dimensional array or data frame.");
    if (d[1]!=nVars) 
      stop("If givem, the number of rows in 'additionalVarInfo' must equal the number\n",
           "  of variables (columns) in components of 'multiExpr'.");
    if (addIDcolumn)
    {
      ID = cbind(ID, additionalVarInfo);
    } else
      ID = additionalVarInfo;
  }

  if (length(varNames)!=nVars)
    stop("Length of 'varNames' does not match the number of variables (columns) in 'data'.");

  dots = list(...)

  d = dim(trait);
  data = as.matrix(data);
  nVars = ncol(data);

  nSamples = if (is.null(d)) length(trait) else d[1];
  if (nSamples !=nrow(data))
    stop ("Number of observations in 'trait' (", nSamples, 
          ") must equal the number of observations in 'data'.");
  nStats = nStatsFnc(...);

  # Check if there are multiple columns in g
  if (length(d) > 1)
  {
     nTraits = d[2];
     if (nTraits==0) stop("Number of columns in 'trait' is zero ('trait' is empty).");
     nTraits.out = nTraits;
     if (is.null(traitNames))
     {
       if (is.null(colnames(trait)))
       {
         traitNames = spaste("Trait.", 1:nTraits);
         colnames(trait) = traitNames;
       } else 
         traitNames = colnames(trait);
     } else {
       if (callFncOneByOne)
       {
         if (length(traitNames)!=nTraits)
           stop("Length of 'traitNames' must equal the number of traits (columns of 'trait').");
       } else
         nTraits.out = length(traitNames); 
     }
     # For multiple columns, call the single-trait function on each column and assemble output
     # Put output into a 3-dimensional array; if flat output is requested, 
     # turn it into a 2-dimensional output at the end.
     if (callFncOneByOne)
     {
       out = array(NA, dim = c(nVars, nStats, nTraits));
       statNames = rep("", nStats * nTraits);
       for (t in 1:nTraits)
       {
         if (verbose > 0)
           printFlush(paste(spaces, "Working on ", traitNames[t]));
         if (doCollectGarbage) collectGarbage();
         tmp = screeningFnc.singleTrait(data = data, trait = trait[, t], ..., traitName = traitNames[t]);
         out[, , t] = tmp;
         statNames[ ((t-1) * nStats + 1):(t*nStats)] = colnames(tmp);
       }
     } else {
        out = screeningFnc.singleTrait(data = data, trait = trait, ...);
        statNames = rep(statNameFnc(...), nTraits.out);
     }
     if (flatOutput)
     {
       dim(out) = c(nVars, nStats *nTraits.out);
       colnames(out) = spaste(statNames, traitNameSep, rep( traitNames, rep(nStats, nTraits.out)));
       rownames(out) = varNames;
       out = as.data.frame(out);
       if (addIDcolumn) 
         out = data.frame(ID, out);
     } else
       dimnames(out) = list(varNames, statNameFnc(...), traitNames);
  } else {
     out = as.data.frame( screeningFnc.singleTrait(data = data, trait = trait, ...));
     if (addIDcolumn)
         out = data.frame(ID, out);
  }
  out;
}
  

# User-level standard screening functions. Basically wrappers for standardScreening.generic.

standardScreeningBinaryTrait = function(data, trait,
           corFnc = cor, corOptions = list(use = 'p'), corName = "cor",
           tTestAlternative = c("two.sided", "less", "greater"),
           var.equal = FALSE,
           getQvalues = FALSE,
           getFDR = FALSE,
           getZstatistics = FALSE,
           getAreaUnderROC = TRUE,
           dataIsLogTransformed = FALSE,
           logBase = 2,
           logShift = 1,
           consistentSigns = TRUE,
           reportLogFoldChange = FALSE,

           verbose = 1,
           indent = 0,
           doCollectGarbage = FALSE,

           flatOutput = FALSE,
           addIDcolumn = TRUE,
           traitNameSep = ".for.",
           traitNames = NULL,
           varNames = NULL,
           additionalVarInfo = NULL)
{
  tTestAlternative = match.arg(tTestAlternative);
  standardScreening.generic(data, trait, 
               screeningFnc.singleTrait = SSBT.singleTrait,
               statNameFnc = SSBTstatNames,
               nStatsFnc = nSSBTstats,

               corFnc = corFnc,
               corOptions = corOptions,
               corName = corName,
               tTestAlternative = tTestAlternative,
               var.equal = var.equal,
               getQvalues = getQvalues,
               getFDR = getFDR,
               getZstatistics = getZstatistics,
               getAreaUnderROC = getAreaUnderROC,
               dataIsLogTransformed = dataIsLogTransformed,
               logBase = logBase,
               logShift = logShift,
               consistentSigns = consistentSigns,
               reportLogFoldChange = reportLogFoldChange,

               verbose = verbose,
               indent = indent,
               callFncOneByOne = TRUE,
               doCollectGarbage = doCollectGarbage,

               flatOutput = flatOutput,
               addIDcolumn = addIDcolumn,
               traitNameSep = traitNameSep,
               traitNames = traitNames,
               varNames = varNames,
               additionalVarInfo = additionalVarInfo);
}

# Quick test of consistency of signs and group names in standardScreeningBinaryTrait
if (FALSE)
{
  data = matrix(rep(c(0,1), each = 10), nrow = 20, ncol = 5);
  trait = rep(c(0,1), c(8, 12));

  standardScreeningBinaryTrait(data, trait, flatOutput = TRUE, dataIsLogTransformed = TRUE, 
                                reportLogFoldChange = TRUE)
}

standardScreeningNumericTrait = function(data, trait,
         corFnc = cor,
         corOptions = list(use = 'p'),
         alternative = c("two.sided", "less", "greater"),
         getQvalues = TRUE,
         getFDR = FALSE,
         getAreaUnderROC = TRUE,
         getDetails = FALSE,
         getLMCoeffs = FALSE,
         exactZ = TRUE,

         verbose = 1,
         indent = 0,
         doCollectGarbage = FALSE,

         corName = "cor",
         flatOutput = FALSE,
         addIDcolumn = TRUE,
         traitNameSep = ".for.",
         traitNames = NULL,
         varNames = NULL,
         additionalVarInfo = NULL)
{
  alternative = match.arg(alternative);

  standardScreening.generic(data, trait,
               screeningFnc.singleTrait = SSNT.singleTrait,
               statNameFnc = SSNTstatNames,
               nStatsFnc = nSSNTstats,

               corFnc = corFnc,
               corOptions = corOptions,
               corName = corName,
               alternative = alternative,
               getQvalues = getQvalues,
               getFDR = getFDR,
               getAreaUnderROC = getAreaUnderROC,
               getDetails = getDetails,
               getLMCoeffs = getLMCoeffs,
               exactZ = exactZ,

               verbose = verbose, 
               indent = indent,
               callFncOneByOne = TRUE,
               doCollectGarbage = doCollectGarbage,

               flatOutput = flatOutput,
               addIDcolumn = addIDcolumn,
               traitNameSep = traitNameSep,
               traitNames = traitNames,
               varNames = varNames,
               additionalVarInfo = additionalVarInfo);

}

standardScreeningSingleSample = function(data, 
               alternative = c("two.sided", "less", "greater"), 
               getQvalues = TRUE, 
               getFDR = FALSE,
               getZstatistics = FALSE, 
               dataIsLogTransformed, 
               addIDcolumn = TRUE, 
               varNames = NULL,
               additionalVarInfo = NULL)
{
  alternative = match.arg(alternative);

  # Need a dummy trait for standardScreening.generic
  trait = rep(1, nrow(data));

  standardScreening.generic(data, trait,
               screeningFnc.singleTrait = singleSampleScreening,
               statNameFnc = SSSstatNames,
               nStatsFnc = nSSSstats,

               alternative = alternative,
               getQvalues = getQvalues,
               getFDR = getFDR,
               getZstatistics = getZstatistics,
               dataIsLogTransformed = dataIsLogTransformed,

               callFncOneByOne = TRUE,
               addIDcolumn = addIDcolumn,
               varNames = varNames,
               additionalVarInfo = additionalVarInfo);
}


#=================================================================================================
#
# Standard screening with covariates
#
#=================================================================================================

dropNames = function(x)
{
  dimnames(x) = NULL;
  x;
}

.last = function(x) x[length(x)];

correlationLikeStatistic = function(Z, t = NULL, n)
{
  if (is.null(t))
  {
    p = pnorm(Z, lower.tail = TRUE);
    t = qt(p, df = n, lower.tail = TRUE);
  }
  t/sqrt(n+ t^2)
}

termStatNames.associationScreening = function(getQvalues,
   getFDR, getZstatistics, getDoF, getCLS, getCor = FALSE)
{
  c("Coefficient", "t", if (getZstatistics) "Z" else character(0),
    "p", if (getQvalues) "q" else character(0), if (getFDR) "FDR" else character(0),
    if (getCor) "cor" else character(0),
    if (getDoF) "DoF" else character(0),
    if (getCLS) "CLS" else character(0));
}

nTermStatistics.associationScreening = function(getQvalues,
   getFDR, getZstatistics, getDoF, getCLS, getCor = FALSE)
{
  length(termStatNames.associationScreening(getQvalues =getQvalues, getFDR =getFDR,
                    getZstatistics = getZstatistics, getDoF = getDoF, getCLS = getCLS, getCor = getCor));
}

commonStatNames.associationScreening = function(getQvalues,
   getFDR, getZstatistics, getDoF, getCLS)
{
  c(if (getZstatistics) "Z" else character(0), "p", if (getQvalues) "q" else character(0), 
    if (getFDR) "FDR" else character(0),
    if (getDoF) "DoF" else character(0),
    if (getCLS) "CLS" else character(0))
}

nCommonStatistics.associationScreening = function(getQvalues,
   getFDR, getZstatistics, getDoF, getCLS)
{
  length(commonStatNames.associationScreening(getQvalues =getQvalues, getFDR =getFDR,
                    getZstatistics = getZstatistics, getDoF = getDoF, getCLS = getCLS));
}



#### This function is not finished.
DependentCovariates = function(modelMat, termsOfInterest)
{
  nPred = ncol(modelMat);
  termIndex = match(termsOfInterest, colnames(modelMat));
  depCols = dependentColumns(modelMat, testOrder = c(termIndex, setdiff(1:nPred, termIndex)),
                             returnColumnNames = TRUE);
  if (any(depCols %in% termsOfInterest))
    stop("Some of the 'termsOfInterest' are among the auto-removed predictors.\n",
         " Please check and correct the design or data.");

  if (any(! (depCols %in% colnames(predictors))))
      stop("Some of the dependent columns are interaction terms. These cannot\n",
           "be auto-removed at this time; please remove manually:\n",
           depCols[! (depCols %in% colnames(predictors))]);
  predictors = predictors[, ! (names(predictors) %in% depCols)];
}

# This function now assumes to have just one termOfInterest and testName.

associationScreening.singleDesign = function(
   response, 
   predictors, 
   full, 
   # reduced can be NULL to skip the LRT-type output
   reduced = NULL,

   weights = NULL,
   fitFnc = "lm",
   fitOptions = list(),

   autoRemoveDependentCovariates = FALSE,

   termsOfInterest = NULL, # Caution: needs to be specified explicitly when reduced is NULL
   testNames = termsOfInterest,

   pEssentiallyZero = 1e-200,

   getQvalues = FALSE,
   getFDR = TRUE,
   getZstatistics = TRUE,
   getDoF = FALSE,
   getCLS = FALSE,
   getCor = FALSE,

   responseNames = colnames(response),
   addTestNameToColumnNames = TRUE,
   testNameSep = ".for.",
   modelComparisonSuffix = spaste(".for.", paste(testNames, collapse = "|")),

   flatOutput = TRUE,
   interleaveOutput = FALSE,
   addID = TRUE,
   additionalResponseInfo = NULL,

   gcInterval = 5000,
   verbose = 2, indent = 0)
{
  spaces = indentSpaces(indent);
  response = as.matrix(response);
  nResponses = ncol(response);

  if (is.null(responseNames)) responseNames = spaste("Response.", 1:nResponses);

  nSamples = nrow(response);
  predictors = as.data.frame(predictors)

  if (is.null(rownames(predictors)))
    rownames(predictors) = spaste("Sample.", c(1:nrow(predictors)));

  if (nrow(predictors)!=nSamples) stop("'response' and 'predictors' must have the same number of rows.");

  if (length(weights)==1 && is.atomic(weights) && is.na(weights)) weights = NULL;
  if (!is.null(weights))
  {
    if (length(weights)==nSamples) weights = matrix(weights, nrow = nSamples, ncol = nResponses);
    weights = as.matrix(weights);
    if (!isTRUE(all.equal(dim(response), dim(weights)))) 
      stop("If given, 'weights' must have the same dimensions as 'response' or as one column in 'response'");
  }

  fitFnc = match.fun(fitFnc);

  # Get the terms for which we want to record test statistics
  terms.full = attr( terms(as.formula(full)), "term.labels");
  nTerms.full = length(terms.full);
  haveReduced = !is.null(reduced) && reduced!="";

  termsOfInterest.valid = termsOfInterest!="";
  validTermsOfInterest = termsOfInterest[termsOfInterest.valid];
  if (haveReduced)
  {
    terms.reduced = attr( terms(as.formula(reduced)), "term.labels");
    if (length(setdiff(terms.full, terms.reduced))==0) 
       stop("Reduced model contains all explanatory variables of the full model.");
    #if (length(termsOfInterest)==0) termsOfInterest = setdiff(terms.full, terms.reduced);
    p.common = rep(NA, nResponses) 
  } else {
    if (length(termsOfInterest)==0) 
       stop("'termsOfInterest' must be explicitly given when 'reduced' is NULL.");
  }
  nTerms = length(termsOfInterest);
  if (is.null(testNames)) testNames = termsOfInterest;

  predictorModelMat = model.matrix(formula(full), data = predictors);
  if (!all(validTermsOfInterest %in% colnames(predictorModelMat)))
    stop("Some 'termsOfInterest' are not among column names of model matrix:\n  ",
         paste(validTermsOfInterest[!validTermsOfInterest %in% colnames(predictorModelMat)], collapse = ", "));

  if (autoRemoveDependentCovariates)
  {
    nPred = ncol(predictorModelMat);
    termIndex = match(validTermsOfInterest, colnames(predictorModelMat));
    depCols = dependentColumns(predictorModelMat, testOrder = c(termIndex, setdiff(1:nPred, termIndex)),
                               returnColumnNames = TRUE);
    if (any(depCols %in% validTermsOfInterest))
      stop("Some of the 'termsOfInterest' are among the auto-removed predictors.\n",
           " Please check and correct the design or data.");

    if (any(! (depCols %in% colnames(predictors))))
      stop("Some of the dependent columns are interaction terms. These cannot\n",
           "be auto-removed at this time; please remove manually:\n",
           depCols[! (depCols %in% colnames(predictors))]);
    predictors = predictors[, ! (names(predictors) %in% depCols)];
  }

  if (length(testNames)!=length(termsOfInterest))
    stop("If given, 'testNames' must have the same length as 'termsOfInterest'.");

  termStatNames = termStatNames.associationScreening(getQvalues =getQvalues,
                                      getFDR =getFDR, getZstatistics = getZstatistics,
                                      getDoF = getDoF, getCLS = getCLS, getCor = getCor);

  outColIndex = match(termStatNames.associationScreening(getQvalues = FALSE, getFDR = FALSE, 
                     getZstatistics = FALSE, getDoF = FALSE, getCLS = FALSE, getCor = getCor), termStatNames);

  stats.fromFits = array(NA, dim = c(nResponses, 3+getCor, nTerms));
    # The stats are coefficient, t and p-value.
  nDoF = rep(NA, nResponses);
  if (verbose > 1) pind = initProgInd(spaste(spaces, "associationScreening: fitting.."));
  counter = 0;
  for (r in 1:nResponses)
  {
    y1 = response[, r];
    sd1 = var(y1, na.rm = TRUE);
    
    if (is.finite(sd1) && sd1 > 0)
    {
      modelData = data.frame(predictors, response[, r, drop= FALSE]);
      full.1 = formula(paste(.last(names(modelData)), full));
      mm = model.matrix(full.1, data = modelData);
      nDoF[r] = sum(rowSums(is.na(mm))==0) - nTerms.full -1 ;
      keepSamples = match(rownames(mm), rownames(modelData));
      modelData.keep = modelData[keepSamples, ];
      if (haveReduced) reduced.1 = formula(paste(.last(names(modelData)), reduced));
      if (!is.null(weights)) w1 = weights[keepSamples, r] else w1 = NULL;

      xx = try({
        fit.full = do.call(fitFnc, c(list(formula = full.1, data = modelData.keep, weights = w1), 
                                          fitOptions));
        if (haveReduced)
        {
          fit.reduced = do.call(fitFnc, c(list(formula = reduced.1, data = modelData.keep, weights = w1), 
                                            fitOptions));
          aov = anova(fit.full, fit.reduced)
          p.common[r] = aov[2, ncol(aov)]; 
          commonRow = which(termsOfInterest=="");
          if (length(commonRow) > 0)
          {
            stats.fromFits[r, 3, commonRow] = aov$`Pr(>F)`[2];
            stats.fromFits[r, 1, commonRow] = aov$F[2];
            stats.fromFits[r, 2, commonRow] = qt(stats.fromFits[r, 3, commonRow], nDoF[r]);
          }
        }
        if (any(termsOfInterest.valid))
        {
          coeffs = summary(fit.full)$coefficients;
          rows = match(validTermsOfInterest, rownames(coeffs));
          rows.out = match(validTermsOfInterest, termsOfInterest);
          stats.fromFits[r,1:3 , rows.out] = t(coeffs[rows, c(1,3,4), drop = FALSE]);
          if (getCor)
          {
            cr1 = cor(response[keepSamples, r, drop= FALSE], 
                      mm[ , setdiff(colnames(mm), "(Intercept)"), drop = FALSE], use = 'p');
            stats.fromFits[r, 4, rows.out] = c(cr1[, match(validTermsOfInterest, colnames(cr1))]);
          }
        }
      })
      if (inherits(xx, "try-error")) 
      {
        printFlush("Error in associationScreening.singleDesign, droppping into browser. "); browser();
      }
    }
    if (verbose > 1 && as.integer(r*100/nResponses)>counter) 
    {
      pind = updateProgInd(r/nResponses, pind);
      counter = as.integer(r*100/nResponses)
    }
    if (r%%gcInterval==0) gc();
  }
  if (verbose > 1) printFlush("");

  # Put together output for each coefficient

  stats.out = lapply(testNames, function(name, mat, colNames)
  {
    if (addTestNameToColumnNames)
    {
      colnames(mat) = spaste(colNames, testNameSep, name);
    } else colnames(mat) = colNames; 
    rownames(mat) = responseNames;
    mat
  }, matrix(NA, nResponses, nTermStatistics.associationScreening(getQvalues =getQvalues, 
                                      getFDR =getFDR, getZstatistics = getZstatistics,
                                      getDoF = getDoF, getCLS = getCLS, getCor = getCor)),
     termStatNames);

  pCol = match("p", termStatNames);
  qCol = match("q", termStatNames);
  zCol = match("Z", termStatNames);
  fdrCol = match("FDR", termStatNames);
  coeffCol = match("Coefficient", termStatNames);
  DoFCol = match("DoF", termStatNames);
  clsCol = match("CLS", termStatNames);
  for (t in 1:nTerms)
  {
    stats.out[[t]] [, outColIndex] = stats.fromFits[, , t];
    # Check p column for essentially zero values
    p1 = stats.out[[t]][, pCol];
    p1[p1 < pEssentiallyZero] = pEssentiallyZero;
    stats.out[[t]][, pCol] = p1;
    Z1 = sign(replaceMissing(stats.out[[t]][, coeffCol], 1)) * qnorm(stats.out[[t]] [, pCol]/2, lower.tail = FALSE);
    if (getZstatistics) stats.out[[t]] [, zCol] = Z1;

    if (getQvalues)
      stats.out[[t]] [, qCol] = qvalue.restricted(stats.out[[t]] [, pCol]);

    if (getFDR)
      stats.out[[t]] [, fdrCol] = p.adjust(stats.out[[t]] [, pCol], method = "fdr");

    if (getDoF)
      stats.out[[t]] [, DoFCol] = nDoF;

    if (getCLS)
      stats.out[[t]] [, clsCol] = correlationLikeStatistic(Z = Z1, n = nDoF);

    if (!termsOfInterest.valid[t])
      colnames(stats.out[[t]]) = sub("^Coefficient", "F", colnames(stats.out[[t]]))
  }

  # Put together output for model significance statistics

  if (flatOutput)
  {
    if (interleaveOutput) {
      out = interleave(stats.out, nameBase = rep("", nTerms), sep = "")
    } else out = do.call(cbind, stats.out)
    if (!is.null(additionalResponseInfo)) 
      out = data.frame(additionalResponseInfo, out);
    if (addID) out = data.frame(ID = responseNames, out);
  } else {
    out = stats.out;
  }

  out;
}

associationScreening = function(
   response, predictors,
   full,
   # reduced can be NULL to skip the LRT-type output
   reduced = NULL,
   termsOfInterest = NULL, # Caution: needs to be specified explicitly when reduced is NULL
   testNames = termsOfInterest,
   ...,
   getFDR = TRUE,
   getQvalues = FALSE,
   FDRFromCombinedPValues = TRUE,
   flatOutput = TRUE,
   interleaveOutput = FALSE,
   addID = TRUE,
   additionalResponseInfo = NULL,
   responseNames = colnames(response),

   gcInterval = 5000,
   verbose = 2, indent = 0)
{
  nDesigns = length(full);
  haveReduced = length(reduced)>0;

  response = as.matrix(response);
  nResponses = ncol(response);

  if (is.null(responseNames)) responseNames = spaste("Response.", 1:nResponses);

  if (!(getFDR || getQvalues)) FDRFromCombinedPValues = FALSE;
  if (haveReduced)
  {
     if (length(reduced)!=length(full))
        stop("Lengths of 'full' and 'reduced' designs must be the same.");
     if (!is.null(termsOfInterest) && length(termsOfInterest)!=nDesigns)
        stop("When given, length of 'termsOfInterest' must equal length of 'full'.");
  } else {
     if (length(termsOfInterest)!=nDesigns)
        stop("When given, length of 'termsOfInterest' must equal length of 'full'.");
  }

  as.sd = list();
  for (ds in 1:nDesigns)
  {
    if (verbose > 0) printFlush("Working on full design ", full[ds]);
    as.sd[[ds]] = associationScreening.singleDesign(
       response = response, 
       predictors = predictors,
       full = full[ds],
       reduced = if (haveReduced) reduced[ds] else NULL,
       termsOfInterest = if (is.null(termsOfInterest)) NULL else termsOfInterest[[ds]],
       testNames = if (is.null(testNames)) NULL else testNames[[ds]],
       responseNames = responseNames,
       getFDR = getFDR, getQvalues = getQvalues,
       ...,
       flatOutput = FALSE,
       gcInterval = gcInterval,
       verbose = verbose -1, indent = indent + 1)
  }

  selectASComp = rep(1, nDesigns);
  #as.sd.keep = mymapply(function(x, i) x[[i]], as.sd, selectASComp);
  as.sd.keep = lapply(as.sd, `[[`, 1);
  if (FDRFromCombinedPValues)
  {
    if (FALSE) # Debugging
    {
    allPValues = try(unlist(lapply(as.sd.keep, function(x) 
      x[, grepSingleOrError("^p", colnames(x))])));
    if (inherits(allPValues, "try-error")) browser()
    }
    allPValues = unlist(lapply(as.sd.keep, function(x)
      x[, grepSingleOrError("^p", colnames(x))]));
    outerIndex = unlist(mymapply(function(x, i) rep(i, nrow(x)), as.sd.keep, 1:nDesigns));
    if (getFDR)
    {
      pAdj = p.adjust(allPValues, method = "fdr");
      columns = lapply(as.sd.keep, function(x) grepSingleOrError("^FDR", colnames(x)));
      for (ds in 1:nDesigns) 
        as.sd[[ds]] [[selectASComp[ds] ]] [, columns[[ds]]] = pAdj[outerIndex==ds];
    }
    if (getQvalues)
    {
      pAdj = qvalue.restricted(allPValues, trapErrors = FALSE)
      columns = lapply(as.sd.keep, function(x) grepSingleOrError("^q", colnames(x)));
      for (ds in 1:nDesigns) 
        as.sd[[ds]] [[ts]] [, columns[[ds]] ] = pAdj[outerIndex==ds];
    }
  }
    
  if (flatOutput)
  {
    as.sd = lapply(as.sd, function(x) do.call(cbind, x));
    if (interleaveOutput) {
      out = interleave(as.sd, nameBase = rep("", nTerms), sep = "")
    } else out = do.call(cbind, as.sd)
    out = as.data.frame(out);
    if (!is.null(additionalResponseInfo))
      out = data.frame(additionalResponseInfo, out);
    if (addID) out = data.frame(ID = responseNames, out);
  } else {
    out = as.sd;
    if (addID) out = lapply(out, function(x) data.frame(ID = responseNames, x));
    if (!is.null(additionalResponseInfo))
      out = lapply(out, function(x) data.frame(additionalResponseInfo, x));
  }
  out;
}

#===================================================================================================
#
# Meta-analysis based on associationScreening.singleDesign
#
#===================================================================================================

metaAnalysisWithCovariates.oneVariable = function(
   multiResponse, 
   multiPredictors,

   full, 
   # reduced can be NULL to skip the LRT-type output
   reduced = NULL,

   termOfInterest,
   testName = termOfInterest,

   useLimma = FALSE,

   observationWeights = NULL,

   fitFnc = "lm",
   fitOptions = list(),

   getQvalues = FALSE,
   getFDR = TRUE,
   getZstatistics = TRUE,
   getDoF = FALSE,
   #getCoeffSignificance = is.null(reduced),
   getCLS = FALSE,

   addTestNameToColumnNames = TRUE,
   testNameSep = ".for.",
   metaAnalysisWeights = c("RootDoF", "Equal", "DoF"),
   setNameSep = ".in.",

   setNames = names(response),
   additionalVarInfo = NULL,

   autoRemoveDependentCovars = FALSE,

   gcInterval = 5000,
   verbose = 2, indent = 0)
{

   # Remove rows in which any of the predictors contain missing values
   #keepSamples = mtd.apply(multiPredictors, function(x) rowSums(is.na(x))==0, returnList = TRUE);

   #multiResponse = mtd.subset(multiResponse, keepSamples);
   #multiPredictors = mtd.subset(multiPredictors, keepSamples);
   #if (!is.null(observationWeights)) observationWeights = mtd.subset(observationWeights, keepSamples);

   # prepare for running screening
   nSets = nSets(multiResponse);

   if (is.null(observationWeights)) observationWeights = listRep( list(data = NA), nSets);

   if (useLimma)
   {
      screening = mtd.mapply(limmaTest.fromDesigns,
              response = multiResponse,
              predictors = multiPredictors,
          MoreArgs = list(
              designs = full,
              termsOfInterest = termOfInterest,
              testNames = testName,
              autoRemoveDependentCovariates = autoRemoveDependentCovars,
              getDoF = TRUE,
              getCLS = getCLS,
              collectGarbage = gcInterval < checkSets(multiResponse)$nGenes,
              sort.by = "none",
              addTestName = addTestNameToColumnNames, sep = testNameSep,
              additionalVarInfo = NULL, includeColumnIDs = FALSE,
              flatOutput = TRUE, interleave = FALSE,
              verbose = verbose - 1));
      modelComparisonSuffix = spaste(testNameSep, termOfInterest);
   } else {
      modelComparisonSuffix = spaste(testNameSep, termOfInterest);
      screening = mtd.mapply(associationScreening,
                 response = multiResponse,
                 predictors = multiPredictors,
                 weights = observationWeights,
             MoreArgs = list(
                 full = full, 
                 reduced = reduced,
                 fitFnc = fitFnc,
                 fitOptions = fitOptions,
                 autoRemoveDependentCovariates = autoRemoveDependentCovars,
                 termsOfInterest = termOfInterest, 
                 testNames = testName,
                 getQvalues = getQvalues,
                 getFDR = getFDR,
                 getZstatistics = TRUE,
                 getDoF = TRUE,
                 getCLS = getCLS,
                 addTestNameToColumnNames = addTestNameToColumnNames,
                 testNameSep = testNameSep,
                 modelComparisonSuffix = modelComparisonSuffix,
                 flatOutput = TRUE,
                 addID = FALSE,
                 interleaveOutput = FALSE,
                 additionalResponseInfo = NULL,
                 gcInterval = gcInterval,
                 verbose = verbose , indent = indent + 1));
    }
    Zstats = mtd.apply(screening, getElement, spaste("Z", modelComparisonSuffix), mdaSimplify = TRUE);
    DoF = mtd.apply(screening, getElement, spaste("DoF", modelComparisonSuffix), mdaSimplify = TRUE);

    if (is.character(metaAnalysisWeights)) {
        metaAnalysisWeights = match.arg(metaAnalysisWeights)
        if (metaAnalysisWeights == "Equal") {
          metaAnalysisWeights = rep(1, nSets)
        } else if (metaAnalysisWeights == "RootDoF") {
          metaAnalysisWeights = sqrt(DoF)
        } else if (metaAnalysisWeights == "DoF")
          metaAnalysisWeights = DoF
    } else {
      if (length(metaAnalysisWeights)!=nSets)
        stop("If given as numbers, length of 'metaAnalysisWeights' must equal\n",
             "   the number of sets in 'multiResponse'.");
    }
    ma = metaAnalysis.simple(Zstats, setWeights = metaAnalysisWeights);

    screening = mtd.mapply(function(x, name) setColnames(x, spaste(colnames(x), setNameSep, name)),
                       screening, setNames);
    colnames(ma) = spaste(colnames(ma), testNameSep, termOfInterest);

    names(screening) = NULL;
    if (!getDoF) 
      screening = mtd.apply(screening, function(x) x[, grep("^DoF", colnames(x), invert = TRUE)]);

    out = cbind(ID = mtd.colnames(multiResponse),
                ma,
                do.call(cbind, multiData2list(screening)));

    out;
}


metaAnalysisWithCovariates = function(
   multiResponse,
   multiPredictors,

   full,
   # reduced can be NULL to skip the LRT-type output
   reduced,

   termOfInterest,
   testNames = termOfInterest,

   useLimma = FALSE,

   observationWeights = NULL,
   fitFnc = "lm",
   fitOptions = list(),

   getQvalues = FALSE,
   getFDR = TRUE,
   getZstatistics = TRUE,
   getDoF = FALSE,
   getCLS = FALSE,

   addTestNameToColumnNames = TRUE,
   testNameSep = ".for.",

   metaAnalysisWeights = c("RootDoF", "Equal", "DoF"),
   setNameSep = ".in.",
   setNames = names(multiResponse),
   additionalVarInfo = NULL,

   flatOutput = TRUE,
   addIDColumn = TRUE,
   interleave = FALSE,

   autoRemoveDependentCovars = FALSE,

   gcInterval = 5000,
   verbose = 2, indent = 0)
{
  if (!is.character(full)) stop("'full' must be a character vector, not formula(s).");

  nDesigns = length(full)
  if (!is.null(reduced))
  {
    if (!is.character(reduced)) stop("'reduced' must be a character vector, not formula(s).");
    if (length(reduced)!=nDesigns)
      stop("If given, length of 'reduced' must be the same as length of 'full'.");
  }

  if (length(termOfInterest)!=nDesigns)
    stop("Length of 'termOfInterest' must equal length of 'full'.");

  if (length(testNames)!=nDesigns)
    stop("Length of 'testNames' must equal length of 'full'.");


  nResponses = checkSets(multiResponse)$nGenes;
  if (!is.null(additionalVarInfo))
  {
    if (nrow(additionalVarInfo)!=nResponses) 
      stop("If 'additionalVarInfo' is given, its number of rows must\n",
           "equal the number of variables (columns) in 'multiResponse'");
  }

  out = lapply(1:nDesigns, function(d)
  {
    if (verbose > 0) printFlush(spaste("Working on design", full[d]));

    metaAnalysisWithCovariates.oneVariable(
         multiResponse = multiResponse,
         multiPredictors = multiPredictors,

         full = full[d],
         reduced = if (!is.null(reduced)) reduced[d] else NULL,

         useLimma = useLimma,
         observationWeights = observationWeights,
         fitFnc = fitFnc,
         fitOptions = fitOptions,

         termOfInterest = termOfInterest[[d]],
         testName = testNames[d],

         getQvalues = getQvalues,
         getFDR = getFDR,
         getZstatistics = getZstatistics,
         getDoF = getDoF,
         getCLS = getCLS,

         addTestNameToColumnNames = TRUE,
         testNameSep = testNameSep,

         metaAnalysisWeights = metaAnalysisWeights,
         setNameSep = setNameSep,

         setNames = setNames,
         additionalVarInfo = NULL,

         autoRemoveDependentCovars = autoRemoveDependentCovars,
         gcInterval = gcInterval,
         verbose = verbose-1, indent = indent+1);
  });

  if (flatOutput)
  {
    names(out) = NULL;
    IDs = out[[1]]$ID;
    out = lapply(out, function(x) x[, -1, drop = FALSE]);
    out2 = if (interleave) interleave(out, nameBase = rep("", nDesigns), sep = "") else
                  do.call(cbind, out);

    if (!is.null(additionalVarInfo)) out2 = cbind(as.data.frame(additionalVarInfo), out2);
    if (addIDColumn) out2 = data.frame(ID = IDs, out2);
  } else {
    out2 = out;
    if (!is.null(additionalVarInfo)) 
        out2 = lapply(out2, function(x) data.frame(ID = x[, 1], additionalVarInfo, x[, -1]))
  }

  out2;
}

#===================================================================================================
#
# Association with covariates
#
#===================================================================================================

significanceWithCovariates.1trait = function(expr, trait, covariates, traitName = "Trait")
{
  nExpr = ncol(expr);
  if (is.null(covariates))
  {
    null = lm(expr~1, data = data.frame(v1 = rep(1, nrow(expr))));
  } else {
    null = lm(expr~., data = as.data.frame(covariates));
  }
  trait.df = data.frame(trait = as.numeric(trait));
  names(trait.df) = traitName;
  if (is.null(covariates))
  {
    df.full = trait.df;
  } else {
   df.full = cbind(trait.df, covariates);
  }
  full = lm(expr~., data = df.full);
  null.sum = summary(null);
  full.sum = summary(full);

  rsq.null = sapply(null.sum, getElement, "r.squared");
  rsq.full = sapply(full.sum, getElement, "r.squared");

  coeff.trait = sapply(full.sum, function(x) x$coefficients[2, 1]);

  # Significance is the multi-variate generalization of correlation of trait with variables in expr
  significance = sign(coeff.trait) * sqrt(rsq.full - rsq.null);

  # P-value will be determined from F statistics
  RSS.null = sapply(null.sum, function(s) sum(s$residuals^2));
  RSS.full = sapply(full.sum, function(s) sum(s$residuals^2));

  dof.full = sapply(full.sum, function(s) s$df[2]);

  Fs = (RSS.null - RSS.full) * dof.full/RSS.full;

  p = mapply(pf, Fs, 1, dof.full, lower.tail = FALSE);
  Z = qnorm(p, lower.tail = FALSE) * sign(coeff.trait);

  list(GS = significance, FStatistic = Fs, p = p, Z = Z);
}

significanceWithCovariates = function(expr, traits, covariates)
{
  nTraits = ncol(traits);

  lst = list();
  for (t in 1:nTraits)
    lst[[t]] = significanceWithCovariates(expr, traits[, t], covariates, traitName = colnames(traits)[t]);

  GS = sapply(lst, getElement, "GS");
  FStatistic = sapply(lst, getElement, "FStatistic");
  p = sapply(lst, getElement, "p");
  Z = sapply(lst, getElement, "Z");

  colnames(GS) = colnames(FStatistic) = colnames(p) = colnames(Z) = colnames(traits);
  rownames(GS) = rownames(FStatistic) = rownames(p) = rownames(Z) = colnames(expr);
  
  list(GS = GS, FStatistic = FStatistic, p = p, Z = Z);
}



#====================================================================================================
#
# Multivariate standard screening
#
#====================================================================================================

multivariateSignificance = function(data, traits, designs, suppressOutputFor = "(Intercept)",
                                         simplifyOutput = FALSE, outputListNames = designs)
{
  if (!is.character(designs)) stop("'designs' must be a character vector (not a formula).");
  nDesigns = length(designs);
  nSamples= nrow(data);
  t.out = p.out = z.out = list();
  data = as.matrix(data);
  for (d in 1:nDesigns)
  {
    dm = model.matrix(formula(designs[d]), data = traits);
    keepSamples = rowSums(is.na(dm))==0;
    if (sum(keepSamples)!=nSamples)
    {
       data1 = data[keepSamples, , drop = FALSE];
       traits1 = traits[keepSamples, , drop = FALSE];
    } else {
       data1 = data;
       traits1 = traits;
    }
    nVariables = ncol(dm);
    fit = lm(formula(paste("data1", designs[d])), data = traits1);
    sum = summary(fit);
    coeffs = lapply(sum, getElement, "coefficients");

    # t-statistics
    t0 = as.numeric(sapply(coeffs, function(c) c[, 3]))
    dim(t0) = c(nVariables, ncol(data));
    t = t(t0);
    colnames(t) = colnames(dm);
    rownames(t) = colnames(data);

    keepCols = !colnames(t) %in% suppressOutputFor;
    t = t[, keepCols, drop = FALSE];
    

    # p-values 
    p0 = as.numeric(sapply(coeffs, function(c) c[, 4]))
    dim(p0) = c(nVariables, ncol(data));
    p = t(p0)[, keepCols, drop = FALSE];
 
    
    # Z statistics
    z = -qnorm(p/2) * sign(t);
    dimnames(z) = dimnames(p) = dimnames(t);

    t.out[[d]] = t;
    p.out[[d]] = p;
    z.out[[d]] = z;
  }

  if (!is.null(outputListNames)) names(t.out) = names(p.out) = names(z.out) = outputListNames;

  if (nDesigns==1 && simplifyOutput) {
     t.out = t.out[[1]];
     p.out = p.out[[1]];
     z.out = z.out[[1]];
  }

  list(t = t.out, p = p.out, Z = z.out);
}

#=================================================================================================
#
# Return standard statistics from limma's eBayes
#
#=================================================================================================

ordinaryLimmaTestStatistics = function(fit, coefName = NULL, interleave = TRUE, orderBy = NULL)
{
  nCoeffs0 =ncol(fit$coefficients);
  if (!is.null(coefName)) {
    coefIndex = match(coefName, colnames(fit$coefficients));
    if (any(is.na(coefIndex))) stop("The name ", coefName, " is not among the coefficient names in 'fit'.");
  } else 
    coefIndex = 1:nCoeffs0;

  nCoeffs = length(coefIndex);
  ordinary.t <- (fit$coefficients / fit$stdev.unscaled / fit$sigma)[, coefIndex, drop = FALSE];
  df.residual <- fit$df.residual;
  df.mat = matrix(df.residual, length(df.residual), nCoeffs);
  p.value <- 2 * pt(-abs(ordinary.t), df = df.mat)
  fdr = as.matrix(apply(p.value, 2, p.adjust, method = "fdr"));
  colnames(fdr) = coefName;
  means = fit@.Data[[8]];
  if (!is.null(orderBy))
  {
    order = switch(orderBy, p = order(p.value), t = order(ordinary.t));
    if (is.null(order)) stop("'orderBy' must be one of 'p' or 't'.");
    out = list(coefficients = fit$coefficients[order, coefIndex, drop = FALSE],
                           t = ordinary.t[order, ,drop = FALSE],
                           p = p.value[order, ,drop = FALSE],
                           FDR = fdr[order, ,drop = FALSE],
                           mean = means[order]);
  } else {
    out = list(coefficients = fit$coefficients[, coefIndex, drop = FALSE],
                            t = ordinary.t,
                            p = p.value,
                            FDR = fdr,
                            mean = means);
  }
  if (interleave) 
    out = data.frame(mean = out$mean, 
                     interleave(out[1:4], nameBase = c("Coef", "t", "p", "FDR")));

  out;
} 

#==================================================================================================
#
# Conveniece wraper for limma linear model fit
#
#==================================================================================================

limmaTest = function(response, pheno, covar, phenoName = "trait", collectGarbage = FALSE,
                     sort.by = "none",
                     addPhenoName = TRUE, sep = ".for.",
                     additionalVarInfo = NULL, includeColumnIDs = TRUE,
                     flatOutput = TRUE, interleave = FALSE,
                     verbose = 0)
{
  if (!is.null(dim(pheno)))
  {
     out = mymapply(function(ph, name) limmaTest(response, ph, covar, phenoName = name,
                                           sort.by = sort.by,
                                           includeColumnIDs = includeColumnIDs, 
                                           additionalVarInfo = NULL, 
                                           addPhenoName = addPhenoName | flatOutput, 
                                           sep = sep, verbose = verbose),
                      pheno, colnames(pheno));
     if (flatOutput)
     {
       names(out) = NULL;
       index = c((includeColumnIDs+1): ncol(out[[1]]))
       col1 = out[[1]][, 1];
       out = lapply(out, `[`, , index);
       if (interleave)
       {
          interleave(out, nameBase = rep("", length(out)), sep = "") 
       } else
          do.call(cbind, out);
       if (includeColumnIDs) out = data.frame(ID = col1, out);
       if (!is.null(additionalVarInfo)) out = data.frame(additionalVarInfo, out);
     } else
       names(out) = colnames(pheno);
     return(out);
  }

  if (verbose > 0) printFlush("Working on ", phenoName);
  keepSamples = !is.na(pheno);
  predictors = data.frame(pheno = pheno, covar)[keepSamples, ];
  predictors = dropConstantColumns(predictors)
  designMat = model.matrix(~., data = predictors);
  
  if (collectGarbage) gc();
  fits = lmFit(t(response[keepSamples, ]), design = designMat)
  if (collectGarbage) gc();

  efits = eBayes(fits)
  tt = topTable(efits, coef = "pheno", sort.by = sort.by, number = Inf,
                genelist = if(includeColumnIDs) efits$genes else NULL);
  Z = ZfromP(tt$t, tt$P.Value);
  tt = cbind(tt, Z= Z);
  colnames(tt)[ (includeColumnIDs+1): ncol(tt)] = spaste(colnames(tt)[ (includeColumnIDs+1): ncol(tt)],
                                                           sep, phenoName);
  if (!is.null(additionalVarInfo)) tt = cbind(additionalVarInfo, tt);
  tt;

}


limmaTest.fromDesigns = function(response, predictors, 
                     designs,
                     termsOfInterest,
                     testNames = termsOfInterest,
                     autoRemoveDependentCovariates = FALSE,
                     getDoF = FALSE,
                     getCLS = FALSE,
                     collectGarbage = FALSE,
                     sort.by = "none",
                     addTestName = TRUE, sep = ".for.",
                     additionalVarInfo = NULL, includeColumnIDs = TRUE,
                     flatOutput = TRUE, interleave = FALSE,
                     verbose = 0)
{
  require("limma");
  if (!is.character(designs)) stop("'designs' must be a character vector.");

  if (length(designs)!=length(termsOfInterest))
     stop("Length of 'design' and 'termOfInterest' must be the same.");

  if (length(designs) > 1)
  {
     out = mymapply(function(ds, traitName, testName) 
               limmaTest.fromDesigns(response, predictors, designs = ds,
                                     termsOfInterest = traitName,
                                     testNames = testName,
                                     sort.by = sort.by,
                                     includeColumnIDs = includeColumnIDs,
                                     additionalVarInfo = NULL,
                                     addTestName = addTestName | flatOutput,
                                     sep = sep, verbose = verbose),
                      designs, termsOfInteres, testNames);
     if (flatOutput)
     {
       names(out) = NULL;
       index = c((includeColumnIDs+1): ncol(out[[1]]))
       col1 = out[[1]][, 1];
       out = lapply(out, `[`, , index);
       if (interleave)
       {
          interleave(out, nameBase = rep("", length(out)), sep = "")
       } else
          do.call(cbind, out);
       if (includeColumnIDs) out = data.frame(ID = col1, out);
       if (!is.null(additionalVarInfo)) out = data.frame(additionalVarInfo, out);
     } else
       names(out) = colnames(predictors);
     return(out);
  }

  if (is.null(rownames(predictors))) 
     rownames(predictors) = spaste("Sample.", c(1:nrow(predictors)));

  if (verbose > 0) printFlush("Working on ", testNames[1]);
  designMat0 = model.matrix(as.formula(designs), data = predictors);
  keepSamples = rownames(predictors) %in% rownames(designMat0);

  if (autoRemoveDependentCovariates)
  {
    nPred = ncol(designMat0);
    termIndex = match(termsOfInterest, colnames(designMat0));
    depCols = dependentColumns(designMat0, testOrder = c(termIndex, setdiff(1:nPred, termIndex)),
                               returnColumnNames = TRUE);
    if (any(depCols %in% termsOfInterest))
      stop("Some of the 'termsOfInterest' are among the auto-removed predictors.\n",
           " Please check and correct the design or data.");
    designMat = designMat0[, !colnames(designMat0 %in% depCols)]
  } else
    designMat = designMat0;

  if (collectGarbage) gc();
  fits = lmFit(t(response[keepSamples, ]), design = designMat)
  if (collectGarbage) gc();

  efits = eBayes(fits)
  tt = topTable(efits, coef = termsOfInterest, sort.by = sort.by, number = Inf,
                genelist = if(includeColumnIDs) efits$genes else NULL);
  Z = ZfromP(tt$t, tt$P.Value);
  tt = cbind(tt, Z= Z);
  if (getDoF | getCLS)
  {
    nDoF = colSums(!is.na(response[keepSamples, ])) - ncol(designMat);
    if (getDoF) tt = cbind(tt, DoF = nDoF);
  }
  if (getCLS)
    tt = cbind(tt, CLS = correlationLikeStatistic(Z, n = nDoF));
    
  if (addTestName)
    colnames(tt)[ (includeColumnIDs+1): ncol(tt)] = spaste(colnames(tt)[ (includeColumnIDs+1): ncol(tt)],
                                                           sep, testNames);
  if (!is.null(additionalVarInfo)) tt = cbind(additionalVarInfo, tt);
  tt;
}

  

#=================================================================================================
#
# simplifySSBT
#
#=================================================================================================

# Simplify the results of SSBT: keep either Student or Kruskal statistics, depending on whether nPresent is
# above a threshold.

keepStudentOrKruskal = function(ssResults, studentPattern = "Student", kruskalPattern = "Kruskal",
                                nPresentPattern = "nPresent", traitSep = "\\.for\\.", nPresentThreshold = 5,
                                thresholdQuantile = 0.1)
{

  # This is not finished and doesn't work.
  nRes = ncol(ssResults);
  traits = rep("", nRes);
  statNames = rep("", nRes);

  colsToKeep = rep(TRUE, nRes);
  colsToKeep [grep(studentPattern, colnames(ssResults))] = FALSE;
  colsToKeep [grep(kruskalPattern, colnames(ssResults))] = FALSE;

  split = strsplit(colnames(ssResults), split = traitSep, fixed = FALSE);
  lengths = sapply(split, length);
  keep = lengths == 2;

  nameMat = sapply(split[keep], I);
  traits[keep] = nameMat[2, ];
  statNames[keep] = nameMat[1, ];

  traitIndex = rep(NA, nRes);
  traitIndex[keep] = as.numeric(factor(nameMat[2, ], levels = unique(nameMat[2, ])));

  numberCols = grep(nPresentPattern, statNames);
  numberQuantiles = colQuantileC( ssResults[, numberCols], p = thresholdQuantile)
  keepKruskal = tapply(numberQuantiles, traitIndex[numberCols], function(x) { min(x) >= nPresentThreshold });

} 


#=================================================================================================
#
# metaAnalysis
#
#=================================================================================================

# This is an adaptation of the WGCNA function metaAnalysis to the new standardScreening functions.

metaAnalysis = function(multiExpr, multiTrait, 
                        binary = NULL,
                        #consensusQuantile = 0,
                        reportWeights = c("Equal", "RootDoF", "DoF"),
                        metaAnalysisWeights = NULL,

                        alternative = c("two.sided", "less", "greater"),
                        corFnc = cor, corOptions = list(use = 'p'), corName = "cor",
                        getQvalues = FALSE,
                        getFDR = FALSE,
                        getAreaUnderROC = FALSE,
                        exactZ = TRUE,
                        getFoldChange = FALSE,
                        getMeanAndSD = FALSE,
                        useRankPvalue = TRUE,
                        rankPvalueOptions = list(),
                        dataIsLogTransformed = FALSE,
                        logBase = 2,
                        logShift = 1,
                        consistentSigns = TRUE,
                        setNames = NULL, 
                        setNameSep = ".in.",
                        var.equal = FALSE, 
                        metaKruskal = FALSE,
                        # At present this controls whether missing Zs are considered to be 0 (na.exclude) or
                        # whether they are considered missing at random (na.omit).
                        # I should rename this argument to something else though.
                        na.action = "na.exclude",

                        varNames = NULL,
                        additionalVarInfo = NULL,
                        addIDcolumn = TRUE,

                        traitNames = NULL,
                        traitNameSep = ".for.",

                        doCollectGarbage = FALSE)
{

  weightTypes = c("Equal", "RootDoF", "DoF", "User");

  size = checkSets(multiExpr);
  nSets = size$nSets;

  nVars = size$nGenes;

  reportWeights = unique(reportWeights)
  if (!all(tolower(reportWeights) %in% tolower(weightTypes))) 
    stop("Some of the supplied 'reportWeights' are not recognized: ", 
         paste( reportWeights[!tolower(reportWeights) %in% tolower(weightTypes)], collapse = ", "));

  if (!is.null(additionalVarInfo))
  {
    d = dim(additionalVarInfo)
    if (length(d) < 2) stop("If given, 'additionalVarInfo' must be a 2-dimensional array or data frame.");
    if (d[1]!=nVars) 
      stop("If givem, the number of rows in 'additionalVarInfo' must equal the number\n",
           "  of variables (columns) in components of 'multiExpr'.");
  }

  if (is.null(varNames)) varNames = mtd.colnames(multiExpr);
  if (is.null(varNames)) varNames = spaste("Variable.", prependZeros(1:nVars))

  if (length(varNames)!=nVars)
    stop("Length of 'varNames' does not match the number of variables (columns) in 'data'.");

  for (set in 1:nSets)
    multiTrait[[set]] $ data = as.matrix(multiTrait[[set]] $ data);

  alternative = match.arg(alternative);

  tSize = checkSets(multiTrait);
  if (tSize$nGenes > 1)
  {
    # Call self recursively for each individual trait, put results together and return
    nTraits = tSize$nGenes;
    if (is.null(traitNames)) traitNames = mtd.colnames(multiTrait);
    if (is.null(traitNames)) traitNames = spaste("Trait.", c(1:nTraits));
    traitNames= make.unique(traitNames);
    if (length(traitNames)!=nTraits)
      stop("Length of 'traitNames' must equal the number of traits in 'multiTrait'.");
    for (t in 1:nTraits)
    {
       printFlush(spaste("Working on trait ", traitNames[t], " (", t, " of ", nTraits, ")"));
       mTrait1 = mtd.subset(multiTrait, colIndex = t)
       if (doCollectGarbage) collectGarbage();
       ma1 = metaAnalysis(multiExpr, mTrait1, 
                          binary = binary,
                          metaAnalysisWeights = metaAnalysisWeights,
                          alternative = alternative,
                          reportWeights = reportWeights,
                          corFnc = corFnc, corOptions = corOptions, corName = corName,
                          getQvalues = getQvalues,
                          getFDR = getFDR,
                          getAreaUnderROC = getAreaUnderROC,
                          exactZ = exactZ,
                          getFoldChange = getFoldChange,
                          getMeanAndSD = getMeanAndSD,
                          useRankPvalue = useRankPvalue,
                          rankPvalueOptions = rankPvalueOptions,
                          dataIsLogTransformed = dataIsLogTransformed,
                          logBase = logBase,
                          logShift = logShift,
                          consistentSigns = consistentSigns,
                          setNames= setNames,
                          setNameSep = setNameSep,
                          var.equal = var.equal,
                          metaKruskal = metaKruskal,
                          na.action = na.action,
                          varNames = varNames,
                          traitNames = traitNames[t],
                          traitNameSep = traitNameSep,
                          addIDcolumn = addIDcolumn && t==1,
                          additionalVarInfo = if (t==1) additionalVarInfo else NULL);
       if (t==1) {
         out = ma1;
       } else 
         out = cbind(out, ma1);
     }
     return(out);
  }
                           
  if (size$nSets!=tSize$nSets)
     stop("The number of sets in 'multiExpr' and 'multiTrait' must be the same.");

  if (!all.equal(size$nSamples, tSize$nSamples))
     stop("Numbers of samples in each set of 'multiExpr' and 'multiTrait' must be the same.");

  #if (!is.finite(consensusQuantile) || consensusQuantile < 0 || consensusQuantile > 1)
  #   stop("'consensusQuantile' must be between 0 and 1.");

  if (is.null(setNames))
     setNames = names(multiExpr);

  if (is.null(setNames))
     setNames = spaste("Set_", c(1:nSets));

  if (is.null(binary)) binary = WGCNA:::.isBinary(multiTrait);

  if (!is.null(metaAnalysisWeights))
  {
    if (length(metaAnalysisWeights)!=nSets)
      stop("Length of 'metaAnalysisWeights' must equal the number of sets in 'multiExpr'.")
    if (any (!is.finite(metaAnalysisWeights)) || any(metaAnalysisWeights < 0))
      stop("All weights in 'metaAnalysisWeights' must be positive.");
  } else 
    metaAnalysisWeights = rep(1, size$nSets);

  if (length(logBase)==1)
    logBase = rep(logBase, size$nSets);

  if (length(logBase)!=size$nSets) 
    stop("Incorrect 'logBase' length: must be 1 or equal to the number of sets in 'multiExpr'.");

  if (length(logShift)==1)
    logShift = rep(logShift, size$nSets);

  if (length(logShift)!=size$nSets) 
    stop("Incorrect 'logShift' length: must be 1 or equal to the number of sets in 'multiExpr'.");

  setResults = list();

  for (set in 1:size$nSets)
  {
    if (binary)
    {
      setResults[[set]] = SSBT.singleTrait(multiExpr[[set]]$data,
                            as.vector(multiTrait[[set]]$data),  
                            tTestAlternative = alternative,
                            getQvalues = getQvalues, var.equal = var.equal,
                            getFDR = getFDR,
                            corFnc = corFnc, corOptions = corOptions, corName = corName,
                            dataIsLogTransformed = dataIsLogTransformed,
                            logBase = logBase[set],
                            logShift = logShift[set],
                            traitName = mtd.colnames(multiTrait),
                            getZstatistics = TRUE,
                            getAreaUnderROC = getAreaUnderROC,
                            consistentSigns = consistentSigns,
                            reportLogFoldChange = TRUE);
      if (metaKruskal) 
      {
        metaStat = "Z.Kruskal";
      } else {
        metaStat = "Z.Student";
      }
      metaEffects = c(corName, "log2FoldChange");
      ## FIXME: should also include the means, or at least one of them.
    } else {
      setResults[[set]] = SSNT.singleTrait(multiExpr[[set]]$data,
                            as.vector(multiTrait[[set]]$data), getQvalues = getQvalues,
                            getFDR = getFDR,
                            alternative = alternative,
                            corFnc = corFnc, corOptions = corOptions, corName = corName, 
                            getAreaUnderROC = getAreaUnderROC, getDetails = getMeanAndSD, exactZ = exactZ,
                            getLMCoeffs = getFoldChange);
      metaStat = "Z";
      metaEffects = c(if (getMeanAndSD) "mean" else character(0), corName, 
                      if (getFoldChange)"coefficient" else character(0));
    }
  }



  comb = as.matrix(interleave(setResults, 
              nameBase = spaste(if (!is.null(traitNames)) spaste(traitNameSep, traitNames[1]) else character(0),
                                setNameSep, setNames),
              baseFirst = FALSE, sep = ""));

  #comb = NULL;
  #for (set in 1:nSets)
  #{
  #  if (set==1) 
  #  {
  #    comb = setResults[[set]];
  #    colNames= colnames(comb);
  #    nColumns = ncol(comb);
  #    colnames(comb) = spaste("X", c(1:nColumns));
  #  } else {
  #    xx = setResults[[set]];
  #    colnames(xx) = spaste("X", c(1:nColumns));
  #    comb = rbind(comb, xx);
  #  }
  #}
#
#  # Re-arrange comb:
#
#  comb = matrix(as.matrix(as.data.frame(comb)), size$nGenes, nColumns * nSets);
#
#  colnames(comb) = spaste( rep( colNames, rep(nSets, nColumns)), 
#                           if (!is.null(traitNames)) spaste(traitNameSep, traitNames[1]) else "",
#                           setNameSep, rep(setNames, nColumns));
#
  # Get the "effects" (statistics) which should be averaged to find a meta-analysis effect
  effectsForMA = lapply(metaEffects, function(name) 
  {
    out = do.call(cbind, lapply(setResults, function(x) x[, name]));
    colnames(out) = setNames;
    out;
  });

  # Find the columns from which to do meta-analysis
  statCols = grep(spaste("^", metaStat), colnames(comb));
  if (length(statCols)==0) stop("Internal error: no columns for meta-analysis found. Sorry!");
  if (length(statCols)!=nSets) 
    stop(spaste("Columns for meta-analysis are not unique. This could be due to an internal error\n",
                "  or because a trait or set name matches one of the names of the meta-analysis \n",
                "  statistics. To facilitate solving this problem, here is the relevant information:\n",
                "  meta statistic search pattern: ^", metaStat, "\n",
                "  Column names of the combined results in which meta statistic is searched for: \n",
                fixLabels(paste( colnames(comb), collapse = ", "), maxCharPerLine = 60, split = " ")));
 
  setZ = comb[, statCols, drop = FALSE];
  setZ[!is.finite(setZ)] = NA;
  colnames(setZ) = spaste("Z", if (!is.null(traitNames)) spaste(traitNameSep, traitNames[1]) else "",
                          setNameSep, setNames);
  nObsCols = grep(if (binary) "^nPresent.bothGroups" else "^nPresent", colnames(comb));
  nObs = comb[, nObsCols, drop = FALSE];

  useWeightIndex = match(tolower(reportWeights), tolower(weightTypes));
  powers = c(0, 0.5, 1, NA)[useWeightIndex];
  nMeta = length(useWeightIndex);

  metaNames = c("equalWeights", "RootDoFWeights", "DoFWeights", "userWeights")[useWeightIndex];
  if (nMeta==1) metaNames = "metaAnalysis";

  metaResults = NULL;
  for (m in 1:nMeta)
  {
    if (useWeightIndex[m] <=3) 
    {
      weights = nObs^powers[m]
    } else
      weights = matrix( metaAnalysisWeights, size$nGenes, nSets, byrow = TRUE);

    if (!all.equal(dim(setZ), dim(weights))) 
       stop("Internal error: inconsistent dimensions of setZ and weights.");

    if (na.action=="na.omit")
      weights[is.na(setZ)] = NA;

    metaZ = rowSums( setZ * weights, na.rm = TRUE) / sqrt(rowSums(weights^2, na.rm = TRUE))
    p.meta = 2*pnorm(abs(metaZ), lower.tail = FALSE);
    if (getQvalues)
    {
      q.meta = qvalue.restricted(p.meta);
      meta1 = cbind(metaZ, p.meta, q.meta)
      colnames.1 = c("Z.", "p.", "q.");
    } else {
      q.meta = NULL;
      meta1 = cbind(metaZ, p.meta);
      colnames.1 = c("Z.", "p.");
    }
    if (getFDR)
    {
      fdr.meta = p.adjust(p.meta, method = "fdr");
      meta1 = cbind(meta1, fdr.meta);
      colnames.1 = c(colnames.1, "FDR.");
    }
    effectMA = do.call(cbind, lapply(effectsForMA, function(x)
         rowSums( x * weights, na.rm = TRUE) / rowSums(weights, na.rm = TRUE)));

    meta1 = cbind(meta1, effectMA);
    colnames.1 = c(colnames.1, spaste(metaEffects, "."));

    colnames(meta1) = spaste(colnames.1, metaNames[m], 
                             if (!is.null(traitNames)) spaste(traitNameSep, traitNames[1]) else "")
    
    metaResults = cbind(metaResults, meta1);
  }

  # Use rankPvalue to produce yet another meta-analysis

  rankMetaResults = NULL;
  if (useRankPvalue)
  {
    rankPvalueOptions$datS = as.data.frame(setZ);
    if (is.na(match("calculateQvalue", names(rankPvalueOptions))))
      rankPvalueOptions$calculateQvalue = getQvalues;
    for (m in 1:nMeta)
    {
      if (useWeightIndex[m] <=3) {
        weights = nObs^powers[m]
      } else
        weights = matrix( metaAnalysisWeights, size$nGenes, nSets, byrow = TRUE);

      # rankPvalue requires a vector of weights... so compress the weights to a vector.
      # Output a warning if the compression loses information.
      nDifferent = apply(weights, 2, function(x) {length(unique(x)) });
      if (any(nDifferent)>1)
        printFlush(paste("Warning in metaAnalysis: rankPvalue requires compressed weights.\n", 
                         "Some weights may not be entirely accurate."));
      rankPvalueOptions$columnweights = colMeans(weights, na.rm = TRUE);
      rankPvalueOptions$columnweights = rankPvalueOptions$columnweights / sum(rankPvalueOptions$columnweights)
      rp = do.call(rankPvalue, rankPvalueOptions);
      colnames(rp) = spaste(colnames(rp), ".", metaNames[m],
                            if (!is.null(traitNames)) spaste(traitNameSep, traitNames[1]) else "");
      rankMetaResults = cbind(rankMetaResults, as.matrix(rp));
    }
  }

  # Put together the output

  out = list(ID = if (addIDcolumn) varNames else NULL,
             additionalVarInfo,
             metaResults,
             rankMetaResults,
             comb,
             NULL);   # The last NULL is necessary so the line below works even if nothing else is NULL

  out = do.call(data.frame, out[ -(which(sapply(out,is.null),arr.ind=TRUE))])

  out;
}

#==================================================================================================
#
# simple meta-analysis function
#
#==================================================================================================

metaAnalysis.simple = function(statistics, pValues = NULL, 
                               Zs = NULL,
                               method = c("normal", "chi square"),
                               setWeights = NULL, includeMeanStat = FALSE,
                               meanName = "mean", 
                               getStdEffectSize = FALSE,
                               nSamples = NULL)
{
  statistics = as.matrix(statistics);
  method = match.arg(method);
  nGenes = nrow(statistics);
  nSets = ncol(statistics);
  if (is.null(setWeights)) setWeights = rep(1, nSets);
  if (is.null(dim(setWeights)) & nSets!=length(setWeights))
     stop("Numbers of columns in 'statistics' and length of 'setWeights' differ.")

  if (!is.null(dim(setWeights)) & !isTRUE(all.equal(dim(statistics), dim(setWeights))))
     stop("Dimensions of 'statistics' and 'setWeights' must agree.");

  if (is.null(Zs) && is.null(pValues)) {
     Zs = statistics 
  } else  if (is.null(Zs)) {
     pValues = as.matrix(pValues);
     if (nSets!=ncol(statistics)) stop("Numbers of columns in 'statistics' and 'pValues' differ.");
     if (nGenes!=nrow(pValues)) stop("Numbers of rows in 'statistics' and 'pValues' differ.");
     Zs = sign(statistics) * qnorm(pValues/2, lower.tail = FALSE);
  }
  if (is.null(dim(setWeights)))
  {
    weightMat = matrix(setWeights, nGenes, nSets, byrow = TRUE);
  } else
    weightMat = setWeights;

  weightMat[is.na(Zs)] = NA;
  allMissing = rowSums(is.finite(weightMat))==0;
  if (method=="normal")
  {
    Z.meta = rowSums(Zs*weightMat, na.rm = TRUE)/ sqrt(rowSums(weightMat^2, na.rm = TRUE));
    Z.meta[allMissing] = NA;
    p.meta = pnorm(abs(Z.meta), lower.tail = FALSE)*2
  } else {
    ch = rowSums(Zs^2, na.rm = TRUE);
    ch[allMissing] = NA;
    n = rowSums(is.finite(Zs));
    p.meta = pchisq(ch, df = n, lower.tail = FALSE);
    Z.meta = qnorm(p.meta/2, lower.tail = FALSE)
  }
  FDR.meta = p.adjust(p.meta, method = "fdr");
  out = data.frame(Z.meta = Z.meta, p.meta = p.meta, FDR.meta = FDR.meta);
  if (includeMeanStat) 
  {
    out = data.frame(XXX = rowSums(statistics*weightMat, na.rm = TRUE)/rowSums(weightMat, na.rm = TRUE), 
                     out);
    colnames(out) = sub("XXX", meanName, colnames(out));
  }
  if (getStdEffectSize)
  {
    if (length(nSamples)!=nSets)
      stop("Length of 'nSamples' must equal number of columns in 'statistics'. Actual values: ", 
           length(nSamples), ", ", nSets);
    ns = (!is.na(weightMat) + 0) %*% as.matrix(nSamples);
    es = tanh(Z.meta/sqrt(ns-3))
    out = data.frame(out, stdEffectSize = es, effectiveNSamples = ns);
  }
  out;
}


  

#=================================================================================================
#
# rbind to a matrix
#
#=================================================================================================

rbindToMatrix = function(...)
{
  as.matrix(rbind(...));
}

#=================================================================================================
#
# signif on numeric columns of a data frame
#
#=================================================================================================

signifNumeric = function(x, digits, fnc = "signif")
{
  x = as.data.frame(x);
  isNumeric = sapply(x, is.numeric);
  isDecimal = isNumeric;
  if (any(isNumeric)) {
    isDecimal[isNumeric] = sapply(x[, isNumeric, drop = FALSE], function(xx) { any(round(xx)!=xx, na.rm = TRUE)});
  } else browser()
  fnc = match.fun(fnc);
  x[, isDecimal] = do.call(fnc, list(x = x[, isDecimal], digits = digits));
  x;
}


#=================================================================================================
#
# remove PCs from data
#
#=================================================================================================

removePCs = function(data, nRemove, multiDataOutput = FALSE, includeOriginal = FALSE)
{
  nNRemove = length(nRemove);
  nRemove.sort = sort(nRemove);
  maxNRemove = max(nRemove);
  out = vector(mode = "list", length = nNRemove + includeOriginal);
  if (includeOriginal)
  {
    out[[1]] = list(data = data);
    shift = 1;
  } else 
    shift = 0;
  svd = svd(data, nu = maxNRemove, nv = 0);
  PCs = svd$u;
  for (nr in 1:nNRemove)
  {
    x = as.data.frame(svd$u[, 1:nRemove[nr]]);
    colnames(x) = spaste("PC", 1:nRemove[nr]);
    out[[nr + shift]] = list(data = residuals(lm(data~., data = x)));
    dimnames(out[[nr]]$data) = dimnames(data);
  } 
  names = rep("Removed.0.PCs", nNRemove + includeOriginal);
  names[c(1:nNRemove) + shift] = spaste("Removed.", prependZeros(nRemove), ".PCs");
  names(out) = names;
  if (!multiDataOutput)
  {
    if (nNRemove==1) {
       out = out[[1]]$data
    } else {
       out = multiData2list(out);
    }
  }

  list(data = out, PCs = PCs);
}

#=================================================================================================
#
# shorten the results of a GO enrichment analysis
#
#=================================================================================================

shortenGOenrichmentTable = function(tab, 
          keep = c("module", "modSize", "rank", "BonferoniP", "fracOfBkgrModSize", "termOntology",
                   "termName"),
          signifCols = c(4,5),
          signifDigits = 2,
          shortenCol = match("termName", keep),
          maxStrLen = 60,
          newColNames = c("Mod", "Size", "Rank", "p.Bonf", "Frac", "Ont", "Term"))
{
  short = tab[, match(keep, colnames(tab))];
  short[ , signifCols] = signif(short[, signifCols], signifDigits)
  if (is.finite(shortenCol)) short[, shortenCol] = substring(short[, shortenCol], 1, maxStrLen);
  colnames(short) = newColNames;
  rownames(short) = NULL;
  short;
}

printGOenrichmentTable = function(tab,  
            keep = c("module", "modSize", "rank", "BonferoniP", "fracOfBkgrModSize", "termOntology",
                   "termName"),
          signifCols = c(4,5),
          signifDigits = 2,
          shortenCol = match("termName", keep),
          maxStrLen = 60,
          newColNames = c("Mod", "Size", "Rank", "p.Bonf", "Frac", "Ont", "Term"))
{
  tab.print = shortenGOenrichmentTable(tab, keep = keep, signifCols = signifCols, 
                          signifDigits = signifDigits, shortenCol = shortenCol, 
                          maxStrLen = maxStrLen, newColNames = newColNames);

  isNumeric = sapply(tab.print, is.numeric);
  sepDF = as.data.frame(matrix(NA, 1, ncol(tab.print)));
  sepDF[, !isNumeric] = "---------"
  names(sepDF) = names(tab.print);

  split = list();
  levels = unique(tab.print[, 1]);
  nLevels = length(levels);
  for (l in 1:nLevels)
  {
    split[[2*l-1]] = tab.print[ tab.print[, 1]==levels[l], ];
    if (l<nLevels) split[[2*l]] = sepDF;
  }

  final = do.call(rbind, split);
  rownames(final) = NULL;
  print(final)

  invisible(final);
}

#=================================================================================================
#
# Format a data frame for saving as a plain-text table
#
#=================================================================================================

dataForTable = function(x, transpose = TRUE, IDcolName = "ID", check.names = FALSE)
{
  if (transpose) out = t(x) else out = x;
  out = data.frame(ID = rownames(out), out, check.names = check.names);
  colnames(out)[1] = IDcolName;
  out;
}

# "Inverse": load a table from read.table

loadTable = function(transpose = FALSE, convertToMatrix = FALSE,
                     sql = FALSE, removeQuotes = FALSE, ...)
{
  if (sql) {
    data = read.csv.sql(...)
    if (removeQuotes)
    {   
      data = as.data.frame(lapply(data, function(x) { gsub('"', '', x, fixed = TRUE) }));
      data = as.data.frame(lapply(data, function(x) { gsub("'", '', x, fixed = TRUE) }));
    }
  } else {
    data = read.table(...);
  }

  names = data[, 1];
  if (transpose)
  {
    out = t(data[, -1, drop = FALSE]);
    colnames(out) = names;
  } else {
    out = data[, -1, drop = FALSE];
    rownames(out) = names;
  }

  if (convertToMatrix) out = as.matrix(out);
  out;
}

#=================================================================================================
#
# Turn given columns into signs and replace names
#
#=================================================================================================

turnColumnsIntoSign = function(x, turnTarget, nameReplacement, fixed = FALSE)
{
  cols = grep(turnTarget, colnames(x), fixed = fixed);
  if (length(cols)==0)
  {
    printFlush("Warning in turnColumnsIntoSign: The search pattern 'turnTarget'");
    printFlush("  was not found among the column names of x.");
    printFlush("  The input x is returned unchanged.");
    x;
  }
  x[, cols] = apply(x[, cols, drop = FALSE], 2, sign);
  oldNames = colnames(x)[cols];
  colnames(x)[ cols ] = gsub(turnTarget, nameReplacement, oldNames, fixed = fixed);
  x;
}

#=================================================================================================
#
# multiPlot
#
#=================================================================================================

addErrorBars.2sided = function (x, means, upper, lower=NULL, width = strwidth("II"), ...) 
{
    if (!is.numeric(means) | !is.numeric(x) || !is.numeric(upper) ) {
        stop("All arguments must be numeric")
    }
    ERR1 <- upper
    ERR2 <- lower
    for (i in 1:length(means)) {
        segments(x[i], means[i], x[i], ERR1[i], ...)
        segments(x[i] - width/2, ERR1[i], x[i] + width/2, ERR1[i], ...)
        if (!is.null(ERR2))
        {
          segments(x[i], means[i], x[i], ERR2[i], ...)
          segments(x[i] - width/2, ERR2[i], x[i] + width/2, ERR2[i], ...)
        }
    }
}


multiPlot = function( x = NULL, y = NULL, data = NULL,
                      columnX = NULL, columnY = NULL,
                      stdErr = NULL,
                      barHigh = NULL, barLow = NULL,
                      type = "p",
                      xlim = NULL, ylim = NULL, 
                      pch = 1, col = 1, bg = 0, lwd = 1, lty = 1,
                      lineColor = col,
                      cex = 1, barColor = 1,
                      addGrid = FALSE, linesPerTick = NULL, 
                      horiz = TRUE, vert = FALSE, gridColor = "grey30", gridLty = 3,
                      errBar.lwd = 1,
                      plotBg = NULL,
                      newPlot = TRUE,
                      dropMissing = TRUE,
                      ...)
{

  getColumn = function(data, column)
  {
     if (!is.numeric(column)) column = match(column, colnames(data));
     data[, column];
  }

  expand = function(x, n)
  {
    if (length(x) < n) x = rep(x, ceiling(n/length(x)));
    x[1:n];
  }

  if (!is.null(data))
  {
    if (is.null(columnX)) stop("'columnX' must be given.");
    if (is.null(columnY)) stop("'columnY' must be given.");
    
    x = lapply(data, getColumn, columnX);
    y = lapply(data, getColumn, columnY);
  }
  

  if (is.null(x) | is.null(y)) stop("'x' and 'y' or 'data' must be given.");

  if (mode(x)=="numeric") x = as.list(as.data.frame(as.matrix(x)));
  if (mode(y)=="numeric") y = as.list(as.data.frame(as.matrix(y)));

  if (!is.null(stdErr) && mode(stdErr)=="numeric") stdErr = as.list(as.data.frame(as.matrix(stdErr)));
  if (!is.null(barHigh) && mode(barHigh)=="numeric") barHigh = as.list(as.data.frame(as.matrix(barHigh)));
  if (!is.null(barLow) && mode(barLow)=="numeric") barLow = as.list(as.data.frame(as.matrix(barLow)));

  nx = length(x);
  ny = length(y);

  if (nx==1 && ny>1) 
  {
    for (c in 2:ny) x[[c]] = x[[1]];
    nx = length(x);
  }

  if (nx!=ny) stop("Length of 'x' and 'y' must be the same.");

  if (length(stdErr)>0 && length(stdErr)!=ny) stop("If given, 'stdErr' must have the same length as 'y'.");

  if (!is.null(stdErr))
  {
     barHigh = mymapply(`+`, y, stdErr);
     barLow = mymapply(`-`, y, stdErr);
  }
     
  if (length(barHigh)>0 && length(barHigh)!=ny) stop("If given, 'barHigh' must have the same length as 'y'.");
  if (length(barLow)>0 && length(barLow)!=ny) stop("If given, 'barLow' must have the same length as 'y'.");

  if (!is.null(barHigh) && is.null(barLow)) 
    barLow = mapply(function(m, u) 2*m-u, y, barHigh, SIMPLIFY = FALSE);

  pch = expand(pch, nx);
  col = expand(col, nx);
  bg = expand(bg, nx);
  lineColor = expand(lineColor, nx);
  lwd = expand(lwd, nx);
  lty = expand(lty, nx);
  cex = expand(cex, nx);
  barColor = expand(barColor, nx);

  if (is.null(xlim)) xlim = range(x, na.rm = TRUE) 
  if (is.null(ylim)) ylim = range(c(y, barLow, barHigh), na.rm = TRUE)

  if (newPlot) 
     plot(x[[1]], y[[1]], xlim = xlim, ylim = ylim, pch = pch[1], col = col[1],
       bg = bg[[1]], lwd = lwd[1], lty = lty[1], cex = cex[1], ..., type = "n");

  if (!is.null(plotBg))
  {
    if (length(plotBg)==1) plotBg = expand(plotBg, nx);
    box = par("usr");
    for (i in 1:nx)
    {
      if (i==1) xl = box[1] else xl = (x[[1]] [i] + x[[1]] [i-1])/2;
      if (i==nx) xr = box[2] else xr = (x[[1]] [i+1]+x[[1]] [i])/2;
      rect(xl, box[3], xr, box[4], border = plotBg[i], col = plotBg[i]);
    }
  }

  if (addGrid)
    addGrid(linesPerTick = linesPerTick, horiz = horiz, vert = vert, col = gridColor, lty = gridLty);

  if (!is.null(barHigh))
    for (p in 1:nx)
      addErrorBars.2sided(x[[p]], y[[p]], barHigh[[p]], barLow[[p]], col = barColor[p],
                          lwd = errBar.lwd);

  if (type %in% c("l", "b")) for (p in 1:nx)
  {
    if (dropMissing) present = is.finite(x[[p]]) & is.finite(y[[p]]) else
       present = rep(TRUE, length(x[[p]]));
    lines(x[[p]][present], y[[p]][present], lwd = lwd[p], lty = lty[p], col = lineColor[p]);
  }
  if (type %in% c("p", "b")) for (p in 1:nx)
      points(x[[p]], y[[p]], pch = pch[p], col = col[p], bg = bg[p], cex = cex[p])
}

#===================================================================================================
#
# dataFrame2colors
#
#===================================================================================================

# An extension of lables2colors and numers2colors: decide type automatically

dataFrame2colors = function(x, ordinalColumns = NULL, maxOrdinalLevels = 8,
                            maxLegendLength = 60, useLabels2colors = TRUE, ...)
{
  x = as.data.frame(x);
  nCols = ncol(x);
  
  legend.long = legend.std = legend.short = names(x);
  colors = NULL;

  if (!is.null(ordinalColumns))
  {
    if (is.numeric(ordinalColumns))
    {  
      if (any(ordinalColumns < 1 | ordinalColumns > nCols)) stop("'ordinalColumns' is out of range.");
    } else {
      ordinalColumns = match(ordinalColumns, names(x));
      if (any(is.na(ordinalColumns)))
        stop(spaste("Some entries in 'ordinalColumns' could not be matched to columns names in x."));
    }
  }

  systemColors = c("black", "red", "green", "blue", "cyan", "magenta", "yellow", "grey");

  if (is.null(ordinalColumns))
  {
    nUnique = sapply(x, function(x1) { length(unique(x1)) });
    isNumeric = sapply(x, is.numeric);
    ordinalColumns = nUnique <= maxOrdinalLevels | !isNumeric;
  } else if (is.numeric(ordinalColumns))
    ordinalColumns = c(1:nCols) %in% ordinalColumns;

  for (c in 1:nCols)
  {
    x1 = x[[c]];
    fx1 = as.factor(x1);
    if (all(is.na(x1))) 
    {
       col1 = rep("grey", length(x1))
    } else {
       unique1 = levels(fx1);
       if (ordinalColumns[c])
       {
         fac = as.numeric(fx1);
         if (useLabels2colors) 
         {
           col1 = labels2colors(fac)
           colorLevels = standardColors(length(unique1));
         } else {
           if (length(unique) > length(systemColors))
             warning(spaste("Warning in dataFrame2colors: number of ordinal levels is larger\n",
                            "   than number of system colors (8). Set useLabels2colors = TRUE\n",
                            "   to obtain a larger pallette."));
           colorLevels = systemColors[ 1:min(length(systemColors), length(unique)) ];
           col1 = fac;
         }
         legend.long[c] = spaste(names(x)[c], ": ",
                          paste ( spaste( unique1, " (", colorLevels, ")"), 
                                      collapse = ", "));
         legend.std[c] = spaste(names(x)[c], ": ", paste(unique1, collapse = ", "));
       } else {
           col1 = numbers2colors(x1, ...);
       }
    }
    if (is.null(colors)) {
       colors = data.frame(col1);
    } else 
       colors = data.frame(colors, col1);
  }

  names(colors) = names(x);

  finalLegend = legend.long;
  tooLong = nchar(finalLegend) > maxLegendLength
  finalLegend[tooLong] = legend.std[tooLong];
  tooLong = nchar(finalLegend) > maxLegendLength
  finalLegend[tooLong] = legend.short[tooLong];

  list(colors = colors, legend = finalLegend, 
       allLegends = data.frame(long = legend.long, standard = legend.std, short = legend.short));
}

#=======================================================================================================
#
# findInterval
#
#=======================================================================================================

# this function basically finds the interval that each genomic position is contained in (or NA if none of
# the intervals match).

# Note: chromosomal positions are assumed to be counted equally in query and interval positions.
# Intervals include both start and end positions (this is different from the UCSC BED format
# conventions).

findInterval = function(queryChromosome, queryPosition, intervalChromosome, intervalStart, intervalEnd,
                        maxDist = 0, querySorted = FALSE, intervalSorted = FALSE, 
                        verbose = 1, indent = 0)
{

  spaces = indentSpaces(indent);

  if (!is.numeric(queryChromosome))
    stop("'queryChromosome' must be numeric. Use a number other than existing numbers for X and/or Y. \n",
         "Make sure the chromosome numbers in 'queryChromosome' and 'intervalChromosome' are the same.");
  if (!is.numeric(intervalChromosome))
    stop("'intervalChromosome' must be numeric. Use a number other than existing numbers for X and/or Y. \n",
         "Make sure the chromosome numbers in 'queryChromosome' and 'intervalChromosome' are the same.");

  if (length(queryChromosome)!= length(queryPosition))
     stop("Lengths of 'queryChromosome' and 'queryPosition' must be the same.");
  
  if ( (length(intervalChromosome)!= length(intervalStart)) | 
       (length(intervalStart)!= length(intervalEnd)) )
     stop("Lengths of 'intervalChromosome', 'intervalStart' and 'intervalEnd' must be the same.");

  nQuery = length(queryPosition);
  nIntervals = length(intervalChromosome);

  # The spacer is at least 2x bigger than the maximum size of a chromosome, so the closest interval on one
  # chromosome cannot be on another chromosome unless the chromosome does not have intervals on it.
  chrSpacer = 2 * 10 ^ ( floor(log10(max(queryPosition, intervalEnd))) + 1)

  queryPosX = queryPosition + queryChromosome * chrSpacer
  intervalStartX = c(-1, intervalStart + intervalChromosome * chrSpacer, 
                         chrSpacer * (max(intervalChromosome, queryChromosome) + 2));
  intervalEndX = c(-1, intervalEnd + intervalChromosome * chrSpacer, 
                       chrSpacer * (max(intervalChromosome, queryChromosome) + 2));

  intervalChromosomeX = c(-1, intervalChromosome, max(intervalChromosome, queryChromosome) + 2)

  if (querySorted)
  {
    queryOrder = 1:nQuery;
  } else {
    queryOrder = order(queryPosX);
    queryPosX = queryPosX[queryOrder];
    queryChromosome = queryChromosome[queryOrder];
  }

  if (intervalSorted)
  {
    intervalOrder = 1:(nIntervals+1);
  } else {
    intervalOrder = order(intervalStartX);
    # Add a fake first and last interval so I don't have to worry about being out of range
    intervalStartX = intervalStartX[intervalOrder];
    intervalEndX = intervalEndX[intervalOrder];
    intervalChromosomeX = intervalChromosomeX[intervalOrder];
  }

  refIndex = 2;
  interval = rep(NA, nQuery); 
  if (verbose > 0) {
    cat(spaste(spaces, "Determining intervals"));
    pind = initProgInd();
    step = max(floor(nQuery/100), 1)
  }

  #cbind(intervalStartX, intervalEndX)
  diffChrDist = maxDist + 1;
  for (q in 1:nQuery)
  {
    pos.1 = queryPosX[q];
    while (pos.1 > intervalEndX[refIndex]) refIndex = refIndex + 1;
    if (pos.1 >= intervalStartX[refIndex])
    {
      # The position is within the interval.
      interval[q] = refIndex-1;
    } else {
      # Position is below the start of the interval. Check this interval and the one before.
      lowIndex = refIndex-1;
      dists = c(if (queryChromosome[q] == intervalChromosomeX[refIndex]) 
                          intervalStartX[refIndex] - pos.1 else diffChrDist,
                if (queryChromosome[q] == intervalChromosomeX[lowIndex]) 
                          pos.1 - intervalEndX[lowIndex] else diffChrDist);
      min = min(dists);
      which = which.min(dists);
      if (min <= maxDist)
        interval[q] = c(refIndex, lowIndex)[which] - 1;
    }
    if (verbose > 0) 
    {
      if (q %% step==0) pind = updateProgInd(q/nQuery, pind);
    }
  }
  if (verbose > 0) {pind = updateProgInd(q/nQuery, pind); printFlush("")}

  out = rep(NA, nQuery);
  out[queryOrder] = intervalOrder[interval]

  #cbind(queryChromosome, queryPosition, queryPosX, out)
  #cbind(intervalChromosome, intervalStart, intervalEnd)

  out;
}

# merge overlapping intervals

mergeOverlappingIntervals = function(data, chrCol, startCol, endCol)
{
  n = nrow(data);

  chr = data[, chrCol];
  start = data[, chrCol];
  end = data[, chrCol];

  overlap = chr[-1] == chr[-n] & start[-1] < end[-n];

  if (sum(overlap)==0)
    return(data[, c(chrCol, startCol, endCol)]);

  index = which(overlap);
  remove = rep(FALSE, n);

  for (i in index)
  {
    iLow = i;
    while (remove[iLow]) iLow = iLow-1;
    iHigh = i+1;
    while (iHigh <= n && start[iHigh] < end[iLow]) iHigh = iHigh + 1;
    iHigh = iHigh - 1;
    end[iLow] = max(end[iLow], end[iHigh]);
    for (ii in (i+1):iHigh) remove[ii] = TRUE;
  }

  out = data.frame(chr = chr[!remove], start = start[!remove], end = end[!remove])
  names(out) = colnames(data)[ c( chrCol, startCol, endCol)];
  out;
}

# The query positions are given in chromosome and position; the reference intervals are given in
# referenceStates. 

chromatinStates = function(chromosome, position, referenceStates, 
                           names = spaste("chr.", chromosome, "..bp.", position),
                           shift = 0)
{
  nRefs = length(referenceStates);
  nQuery = length(position);

  states = matrix("", nQuery, nRefs);

  chromosome[chromosome=="X"] = 23;
  chromosome = as.numeric(chromosome);

  order = order(chromosome, position);

  chr = chromosome[order];
  pos = position[order];

  for (r in 1:nRefs)
  {
    cat(spaste("Working on cell line ", r, ": "));

    int = findInterval(chr, pos, referenceStates[[r]]$chromosomeNumber, 
                       referenceStates[[r]]$Start + shift,
                       referenceStates[[r]]$End, querySorted = TRUE, intervalSorted = FALSE);

    finite = is.finite(int);
    states[finite, r] = referenceStates[[r]]$State[int[finite]];
  }
  rownames(states) = names;
  colnames(states) = names(referenceStates);

  states[order, ] = states;
  states;
}


# general function to retrieve information for given positions from a BED-like data frame

findBEDinfo = function(chromosome, position, bed.df, bedChrCol = "Chromosome", bedStartCol = "Start",
                       bedEndCol = "End", chrStripPattern = "chr", fixed = TRUE,
                       maxDist = 0, getColumns = NULL,
                       querySorted = FALSE, bedSorted = FALSE,
                       names = spaste("chr.", chromosome, "..bp.", position),
                       verbose = 1, indent = 0)
{
  if (is.character(bedChrCol)) bedChrCol = match(bedChrCol, colnames(bed.df));
  if (!is.finite(bedChrCol)) stop(spaste("Chromosomal column not found in 'bed.df'."));

  if (is.character(bedStartCol)) bedStartCol = match(bedStartCol, colnames(bed.df));
  if (!is.finite(bedStartCol)) stop(spaste("Interval start column not found in 'bed.df'."));

  if (is.character(bedEndCol)) bedEndCol = match(bedEndCol, colnames(bed.df));
  if (!is.finite(bedEndCol)) stop(spaste("Interval end column not found in 'bed.df'."));

  bedChr = bed.df[, bedChrCol];

  if (!is.null(chrStripPattern) && chrStripPattern!="")
    bedChr = gsub(chrStripPattern, "", bedChr, fixed = fixed);

  chrLevels = sort(unique(c(chromosome, bedChr)));

  bedChrNum = match(bedChr, chrLevels);
  chrNum = match(chromosome, chrLevels);

  if (is.null(getColumns))
  {
     getColumns = c(1:ncol(bed.df))[-c(bedChrCol, bedStartCol, bedEndCol)];
  } else {
     if (!is.numeric(getColumns)) getColumns = match(getColumns, colnames(bed.df));
     if (any(is.na(getColumns))) stop("'getColumns' contains invalid entries.");
     if (any(getColumns > ncol(bed.df)) | any(getColumns < 1)) 
       stop("'getColumns' contains entries that are out of range.");
  }

  # Note: here assuming that the start locations in the BED data frame are 0-based, and end locations are
  # 0-based locations of the first base outside of the interval (or 1-based locations of the last base).

  intervals = findInterval(chrNum, position, bedChrNum, bed.df[, bedStartCol]+1, bed.df[, bedEndCol],
                           maxDist = maxDist,
                           querySorted = querySorted, intervalSorted = bedSorted,
                           verbose = verbose, indent = indent);

  info = bed.df[intervals, getColumns]

  rownames(info) = names;

  info;
}

     

#=====================================================================================================
#
# Read a wigFix data file
#
#=====================================================================================================

# Note: according to http://genome.ucsc.edu/goldenPath/help/wiggle.html, the wiggle formats are 1-based,
# unlike BED formats (which are 0-based) .

readWigFix = function(file)
{
  content = scan(file = file, what = character(), sep = "\b", quote = "");
  blockDividers = grep("^fixedStep", content, fixed = FALSE);

  nLines = length(content);
  blockStart = blockDividers + 1;
  blockEnd = c(blockDividers[-1] -1, nLines);

  nBlocks = length(blockDividers);

  out = vector(length = nBlocks, mode = "list");
  for (b in 1:nBlocks)
  {
    header = content[ blockDividers[b] ];
    data = as.numeric(content[ blockStart[b]: blockEnd[b] ]);

    headerSplit = strsplit(header, "  *", fixed = FALSE)[[1]]
    if (headerSplit[1] != "fixedStep") 
      stop("Supplied file does not appear to be fixed-step Wiggle format.\n",
           "First line of block starting on line ", blockDividers[b], " does not start with 'fixedStep'.\n", 
           "The offending header: ", header);

    nValues = length(headerSplit)-1;
    valueList.0 = lapply(headerSplit[-1], function(x) {strsplit(x, "=", fixed = TRUE)[[1]] });
    valueList = vector(mode = "list", length = nValues);
    for (v in 1:nValues)
    {
      x = valueList.0[[v]] [2];
      xNum = suppressWarnings( as.numeric(x) );
      valueList[[v]] = if (is.finite(xNum)) xNum else x;
      names(valueList)[v] = valueList.0[[v]] [1];
    }

    out[[b]] = list(header = valueList,
                    data = data);
  }

  out;
}

#======================================================================================================
#
# blackWhite
#
#======================================================================================================

blackWhite = function(n, gamma = 1)
{
  level = ( c(0:(n-1))/(n-1) )^(1/gamma)
  rgb(level, level, level);
}

#======================================================================================================
#
# blackWhite
#
#======================================================================================================

blueWhiteRed.Lab = function(n, bias = 1)
{
  half = floor(n/2);
  c(colorRampPalette(c("blue", "white"), bias = bias, space = "Lab")(half),
    colorRampPalette(c("white", "red"), bias = bias, space = "Lab")(n-half+1)[-1]);
}

blueGreyRed.Lab = function(n, bias = 1, grey = "grey80")
{
  half = floor(n/2);
  c(colorRampPalette(c("blue", grey), bias = bias, space = "Lab")(half),
    colorRampPalette(c(grey, "red"), bias = bias, space = "Lab")(n-half+1)[-1]);
}

#======================================================================================================
#
# fast Euclidean distance calculation
#
#======================================================================================================

dist = function(x, method = "euclidean", diag = FALSE, upper = FALSE, p = 2)
{
  x = as.matrix(x);
  if (method!="euclidean") 
     return(stats::dist(x, method, diag, upper, p));
  if (any(is.na(x)))
  {
     fin = is.finite(x)+0;
     xf = x;
     xf[fin==0] = 0;

     square = xf %*% t(xf);
     square1 = (xf*xf) %*% t(fin);
     fSquare = fin %*% t(fin);
     as.dist(sqrt( (square1 + t(square1) - 2*square)/ (fSquare/max(fSquare))), diag = diag, upper = upper);
  } else {
    #printFlush("using fast Euclidean distance calculation...");
    square = x %*% t(x);
    diagonal = diag(square);
    tmp = outer(diagonal, diagonal, `+`) -2*square;
    tmp[tmp<0] = 0;
    out = as.dist(sqrt(tmp), diag = diag, upper = upper);
    #if (any(!is.finite(out))) browser();
    out;
  }
}


# Distance of two sets of row vectors. Does not work with missing data.

dist2 = function(x, y)
{
  xsq = rowSums(x^2)
  ysq = rowSums(y^2)
  xy = 2 * x %*% t(y);
  out = sqrt(outer(xsq, ysq, `+`) - xy);
  colnames(out) = rownames(y);
  rownames(out) = colnames(x);
  out;
}



#=====================================================================================================
#
# mean/median/whatever impute
#
#=====================================================================================================

imputeFunctionValue = function(x, fnc, ...)
{
  x = as.data.frame(x);
  useColumns = which(colSums(is.na(x)) > 0);
  if (!any(useColumns)) return(x);

  imputedValues = sapply(x[, useColumns], fnc, ...);
  for (c in 1:length(useColumns))
    x[ is.na(x[, useColumns[c] ]), useColumns[c]] = imputedValues[c]

  x;
}

#=====================================================================================================
#
# Imputation of extreme values by prediction
#
#=====================================================================================================

imputeExtremeValues = function(predictors, response, maxZ, minZ = -maxZ, ...)
{
  printFlush("Initial prediction...");

  rg = randomGLM(x = predictors, y = response, ..., keepModels = FALSE)
  predicted1 = rg$predictedOOB.response;
  Z = scale(predicted1-response);
  good = Z >= minZ & Z <= maxZ;
  if (all(good)) return(list (imputed = response, initialPredicted = predicted1, good = good));

  printFlush("Final prediction...");
  rg = randomGLM(x = predictors[good, ], y = response[good],
                 xtest = predictors[!good, ], ..., keepModels = FALSE);
  imputed = response;
  imputed[!good] = rg$predictedTest.response;
  predicted = response;
  predicted[good] = rg$predictedOOB.response;
  predicted[!good] = rg$predictedTest.response;
  list(imputed = imputed, initialPredicted = predicted1, predicted = predicted, good = good);
}

#====================================================================================================
#
# Imputation of missing values by voting linear predictor
#
#====================================================================================================

imputeByLVP = function(data, ...)
{
  nCols = ncol(data);
  imputed = data;
  pind = initProgInd();
  for (c in 1:nCols)
  {
     missing = is.na(data[, c]);
     if (any(missing)) 
     { 
        if (sum(missing) > length(missing)/2)
           stop("Column ", c, " has more than 50% missing values.");

        x = data[!missing, -c];
        xTest = data[missing, -c, drop = FALSE];
        y = data[!missing, c];

        pr = votingLinearPredictor(x, y, xTest, ..., verbose = 0);

        imputed[missing, c] = pr$predictedTest;
     }
     pind = updateProgInd(c/nCols, pind);
  }

  printFlush("");
  imputed;
}

#=====================================================================================================
#
# Adjust data for technical factors in the presence of important variables
#
#=====================================================================================================
  
adjustData = function(data, adjustedCovariates, retainedCovariates = NULL, shrink = TRUE)
{
  collectGarbage();
  adjustedCovariates = as.data.frame(adjustedCovariates);

  if (is.null(retainedCovariates))
    return(residuals(lm(data~., data = adjustedCovariates)));

  mm.adj = model.matrix(data~., data = adjustedCovariates);

  retainedCovariates = as.data.frame(retainedCovariates);

  allCovars = cbind(adjustedCovariates, retainedCovariates);

  collectGarbage();
  fit = lm(data~., data = allCovars);

  coeffs = as.matrix(fit$coefficients);
  coeffs.keep = coeffs[ -c(2:(ncol(mm.adj))), ];

  summ = summary(fit);

  stdErrors = sapply(summ, function(x, keep) {x$coefficients[keep, 2]}, -c(2:(ncol(mm.adj))))

  if (shrink)
    coeffs.keep[-1, ] = coeffs.keep[-1, ]/ ( 1 + stdErrors[-1, ]^2/coeffs.keep[-1, ]^2);

  covars.keep = model.matrix(data~., data = retainedCovariates);


  # coeffs.keep are the coefficients of the kept terms (retainedCovariates)
  res = residuals(fit);

  predicted = array(0, dim = dim(res));
  nCoeffs = nrow(coeffs.keep)

  collectGarbage();
  for (c in 1:nCoeffs)
  {
    predicted = predicted + outer(covars.keep[, c], coeffs.keep[c, ], '*')
  }

  out = predicted + res;
  dimnames(out) = dimnames(data);
  out
}

adjustData.wrong.1 = function(data, adjustedCovariates, retainedCovariates, fitFnc,
                        fitOptions = list())
{
  fitFnc = match.fun(fitFnc);
  fitOptions.1 = fitOptions;
  fitOptions.1$formula = as.formula(data~.);
  fitOptions.1$data = retainedCovariates;
  fit1 = do.call(fitFnc, fitOptions.1);
  res1 = residuals(fit1);
  pred1 = predict(fit1);
 
  fitOptions.2 = fitOptions;
  fitOptions.2$formula = as.formula(res1~.);
  fitOptions.2$data = adjustedCovariates;

  res2 = residuals(do.call(fitFnc, fitOptions.2)); 

  out = pred1 + res2;

  out
}

adjustData.wrong = function(data, adjustedCovariates, retainedCovariates, fitFnc = "lm",
                        fitOptions = list(), splitCall = fitFnc != "lm")
{
  adjustedCovariates = as.data.frame(adjustedCovariates);
  retainedCovariates = as.data.frame(retainedCovariates);
  x = splitCall;
  collectGarbage();
  if (splitCall && length(dim(data))==2)
  {
    nCol = ncol(data)
    out = data;
    pind = initProgInd();
    for (c in 1:nCol)
    {
      out[, c] = adjustData.1(data[, c], adjustedCovariates=adjustedCovariates, 
             retainedCovariates = retainedCovariates, fitFnc = fitFnc, fitOptions = fitOptions)
      pind = updateProgInd(c/nCol, pind);
    }
  } else
    out = adjustData.1(data, adjustedCovariates, retainedCovariates, fitFnc, fitOptions);

  collectGarbage();
  dimnames(out) = dimnames(data);
  out;
}
 

#=========================================================================================================
#
# intramodular and eigengene-based connectivity across modules
#
#=========================================================================================================

kIMandKME = function(expr, labels, MEs = NULL,
                     corAndPvalueFnc = "corAndPvalue", cpvOptions = list(),
                     corFnc = "cor", corOptions = "use = 'p'",
                     networkType = "signed hybrid",
                     power = 6,
                     ignoreLabel = 0,
                     geneAnnot.start = NULL,
                     geneAnnot.end = NULL,
                     altLabels = NULL,
                     includeLabelColumn =TRUE,
                     includeAltLabelColumn = TRUE,
                     labelColumnName = "ModuleLabel",
                     altLabelColumnName = "ModuleColor",
                     includeAltLabelInKMEColumns = FALSE,
                     getKIM = TRUE,
                     getZ.kME = FALSE,

                     # How many top hubs to return and how to label them
                     nTopHubs = 20,
                     hubIDs = colnames(expr),
                     sortHubs= c("none", "kME"),

                     verbose = 2, indent = 0)
{
   
   spaces = indentSpaces(indent);
   nGenes = ncol(expr);
   kIM = kIM.scaled = kME.own = rep(NA, nGenes);
   if (is.null(MEs)) MEs = moduleEigengenes(expr, labels, excludeGrey = TRUE, grey = ignoreLabel)$eigengenes;

   cpvOptions$x = expr;
   cpvOptions$y = MEs;
   cpvFnc.fnc= match.fun(corAndPvalueFnc);

   sortHubs = match.arg(sortHubs);

   kMEdata = do.call(cpvFnc.fnc, cpvOptions);
   kME = kMEdata[[1]];
   Z.kME = kMEdata$Z;
   
   labLevels = sort(unique(labels));
   labLevels = labLevels[ !labLevels %in% ignoreLabel ];
   nLabLevels = length(labLevels)

   colnames(kME) = gsub("ME", "kME.", colnames(kME));
   colnames(Z.kME) = gsub("ME", "Z.kME.", colnames(Z.kME));
   if (includeLabelColumn)
   {
     if (is.null(geneAnnot.start)) 
     {
       geneAnnot.start = data.frame(ModuleLabel = labels);
     } else
       geneAnnot.start = cbind(as.data.frame(geneAnnot.start), ModuleLabel = labels);

     names(geneAnnot.start) [ ncol(geneAnnot.start) ] = labelColumnName;
   }
       
   if (!is.null(altLabels))
   {
     if (includeAltLabelInKMEColumns)
     {
       colME = sub("^kME.", "", colnames(kME));
       positions = match(colME, labels);
       altColNames = altLabels[positions];
       colnames(kME) = paste(colnames(kME), altColNames, sep =".");
       colnames(Z.kME) = paste(colnames(Z.kME), altColNames, sep =".");
     }

     geneAnnot.start = cbind(geneAnnot.start, ModuleColor = altLabels);
     names(geneAnnot.start) [ ncol(geneAnnot.start) ] = altLabelColumnName;
   }

   topHubs = matrix("", nTopHubs, nLabLevels)
   colnames(topHubs) = labLevels;
   rownames(topHubs) = c(1:nTopHubs);

   if (is.null(hubIDs)) stop("Must have valiud hub IDs (usually column names of 'expr').");

   for (ll in labLevels)
   {
     if (verbose > 1) printFlush(paste(spaces, "Working on module", ll));
     inModule = c(1:nGenes)[ll==labels];
     if (getKIM)
     {
       adj = adjacency(expr[, inModule], type = networkType, corFnc= corFnc, corOptions = corOptions);
       kIM.1 = colSums(adj, na.rm = TRUE)-1;
       kIM[inModule] = kIM.1;
       kIM.scaled[inModule] = kIM.1/max(kIM.1, na.rm = TRUE);
     }
     
     kME.1 = kME [inModule, match(ll, substring(colnames(MEs), 3))];
     kME.own[inModule] = kME.1;
  
     # Find the top hubs by kME
     if (length(grep("unsigned", networkType)) > 0) kME.1 = abs(kME.1);
     if (sortHubs == "none")
     {
       rank = rank(-kME.1, ties.method = "average", na.last = TRUE);
       topHubs[c(1:sum(rank<=nTopHubs)), match(ll, labLevels)] = hubIDs[inModule] [rank<=nTopHubs];
     } else if (tolower(sortHubs)=="kme") {
       order = order(-kME.1, na.last = TRUE);
       n1 = min(nTopHubs, length(order));
       topHubs[1:n1, match(ll, labLevels)] = hubIDs[inModule] [order[1:n1]];
     } else stop("Invalid 'sortHubs'.");
   }
   out = data.frame(kME.own = kME.own);
   if (getKIM) out = cbind(out, kIM = kIM, kIM.scaled = kIM.scaled);
   kME.out = if (getZ.kME) interleave(list(kME, Z.kME), nameBase = c("", ""), sep = "") else kME;
   out = cbind(out, kME.out);
   if (!is.null(geneAnnot.start)) out = cbind(as.data.frame(geneAnnot.start), out);
   if (!is.null(geneAnnot.end)) out = cbind(out, geneAnnot.end);
   list(kTable = out, topHubs = topHubs);
}


# Consensus hubs out of output of consensusKME

consensusHubs = function(consKME, 
                         nHubs = 20,
                         moduleCol = "module",
                         selectBy = "meta.Z.RootDoFWeights.kME", 
                         select = c("positive", "negative", "both"),
                         returnColumns = "ID",
                         simpleOutput = TRUE
                         )
{
  selectCols = grep(selectBy, colnames(consKME));
  moduleLevels = gsub(selectBy, "", colnames(consKME)[selectCols]);

  ret.mat = consKME[, match(returnColumns, colnames(consKME)), drop = FALSE];
  select = match.arg(select);

  nModules = length(moduleLevels);
  stats = switch(select, positive = consKME[, selectCols], negative = -consKME[, selectCols],
                         both = abs(consKME[, selectCols]));
  moduleLabel = consKME[, match(moduleCol, colnames(consKME))];

  moduleMembership = lapply(moduleLevels, `==`, moduleLabel);
  orders = mymapply(function(mm, st)
    {
      st[!mm] = NA;
      order(st, decreasing = TRUE, na.last = TRUE);
    }, moduleMembership, stats);

  if (ncol(ret.mat)==1 && simpleOutput)
  {
    hubs = sapply(orders, function(i, data, n) data[i[1:n]], ret.mat[, 1], n=nHubs);
  } else
    hubs = lapply(orders, function(i, data, n) data[, i[1:n]], ret.mat, n=nHubs);

  hubs;
}

  
#============================================================================================
#
# plotEigengenesAndTraits
#
#============================================================================================
# note: formatLabels is included in JamsPlenTFun

plotEigengenesAndTraits = function(MEs, traits, order = NULL,
                                   setLayout = TRUE,
                                   meColors = blueWhiteRed(100),
                                   meNames = colnames(MEs),
                                   marAll = c(1, 12, 2, 1), 
                                   traitColors = dataFrame2colors(traits)$colors,
                                   traitNames = formatLabels(dataFrame2colors(traits)$legend, 
                                                          maxCharPerLine = marAll[2] * 2.5, 
                                                          split = " ", keepSplitAtEOL = FALSE),
                                   cex.main = 1.4, cex.lab = 1, cex.axis = 1,
                                   verticalLinePositions = NULL,
                                   col.verticalLines = "darkgrey",
                                   lty.verticalLines = 2,
                                   ylim = NULL,
                                   signed = TRUE,
                                   traitHeight = 1,
                                   ...)
{
  nSamples = nrow(MEs);
  nModules = ncol(MEs)

  if (is.null(order)) order = c(1:nSamples);

  if (setLayout) 
    layout(matrix(c(1:(nModules+1)), nModules + 1, 1), height = c(rep(1, nModules), traitHeight));

  par(mar = c(0, marAll[2:4]));
  for (me in 1:nModules)
  {
    out = barplot(MEs[order, me], col = numbers2colors(MEs[order, me], signed = signed, colors = meColors),
            names.arg = rep("", nSamples),
            main = meNames[me], ylab = "Eigengene\nexpression", xlab = "",
            cex.main = cex.main, cex.lab = cex.lab, cex.axis = cex.axis, space = 0.2, offset = 0,
            ylim = ylim)
    if (!is.null(verticalLinePositions))
    {
      step = out[2] - out[1];
      edges = c(out[1] - step/2, out + step/2);
      box = par("usr");
      for (vp in verticalLinePositions)
        lines(x = rep(edges[vp+1], 2), y = box[c(3:4)], col = col.verticalLines, 
              lty = lty.verticalLines);
    }
  }
  par(mar = c(marAll[1:2], 0.5, marAll[4]))
  par(lheight = 0.85);
  plotOrderedColors(order, traitColors, rowLabels = traitNames, startAt = 0.5,
                     ...);
  out
}


plotEigengenesAndTraits.multiPage = function(MEs, traits,
   meNames = names(MEs),

   # Paging options
   rowsPerPage = NULL, maxRowsPerPage = 8,

   # Argument controlling common options

   commonYLim = FALSE,

   # Further arguments to plotEigengenesAndTraits
   ylim = NULL,
   signed = TRUE,
   ...)
{
  nr = ncol(MEs);

  if (is.null(rowsPerPage))
  {
    nPages.rows = ceiling(nr/maxRowsPerPage);
    rowsPerPage = allocateJobs(nr, nPages.rows);
  } else
    nPages.rows = length(rowsPerPage);

  if (is.null(ylim) && commonYLim)
  {
    ylim = range(MEs, na.rm = TRUE)
    if (signed) ylim = c(-max(abs(ylim)), max(abs(ylim)));
  }

  for (p in 1:nPages.rows)
  {
    plotEigengenesAndTraits(MEs[, rowsPerPage[[p]]], traits,
              meNames = meNames[ rowsPerPage[[p]] ],
              ylim = ylim, signed = signed, ...);
  }
}







#===================================================================================================
#
# Create a summary table for Mike Palazzolo's database
#
#===================================================================================================

# At this point Mike wants the following: gene ID, marginal association (fold change, Student p-value,
# q-value), module label, module color, position in module (kME, rank of kME), module annotation, module
# significance for status (Student p-value), top 20 hubs in each module.

# My solution for now is to report all these statistics in each single set. For each analysis and each data
# set create a separate table.

PalazzoloTable = function(expr.1set, labels, altLabels, status = NULL, 
                          corFnc = cor, corOptions = list(use = 'p'),
                          dataIsLogTransformed = FALSE,
                          analysisName, setName, GEO.ID, 
                          enrichmentLabels,
                          ...)
{
  nGenes = ncol(expr.1set);

  kInfo = kIMandKME(expr.1set, labels, corFnc = corFnc, corOptions = corOptions, altLabels = altLabels,
                    ...);

  topHubs = apply(kInfo$topHubs, 2, paste, collapse = ", ");

  MEs = moduleEigengenes(expr.1set, labels, excludeGrey = TRUE)$eigengenes

  label2ME = match(labels, substring(colnames(MEs), 3));

  if (nlevels(factor(status[, 1]))<3)
  {
    printFlush("Calling standardScreeningBinaryTrait");
    ss = standardScreeningBinaryTrait(expr.1set, status, corFnc = corFnc, corOptions = corOptions,
                                      getQvalues = TRUE, dataIsLogTransformed = dataIsLogTransformed,
                                      flatOutput = TRUE, addIDcolumn = TRUE,
                                      consistentSigns = TRUE);

    keepCols = c(1, grep("foldChange", names(ss)), grep("[pq].Student", names(ss)));
    out = cbind(ss[keepCols], kInfo$kTable[, c(1:3)]);
    ss.me = standardScreeningBinaryTrait(MEs, status, corFnc = corFnc, corOptions = corOptions,
                                      getQvalues = FALSE,  dataIsLogTransformed = FALSE,
                                      flatOutput = TRUE, addIDcolumn = TRUE, consistentSigns = TRUE)

    meInfo = data.frame( t.Student.moduleEigenene.vs.diagnosis = ss.me$t.Student[ label2ME],
                         p.Student.moduleEigengene.vs.diagnosis = ss.me$p.Student[ label2ME],
                         enrichmentLabels[label2ME, ],
                         TopModuleHubs = topHubs[label2ME]);

  } else {
    printFlush("Calling standardScreeningNumericTrait");
    ss = standardScreeningNumericTrait(expr.1set, status, corFnc = corFnc, corOptions = corOptions,
                                      getQvalues = TRUE, 
                                      flatOutput = TRUE, addIDcolumn = TRUE);
    keepCols = c(1, grep("[a-z]*cor\\.", names(ss)), grep("[pq].Student", names(ss)));

    out = cbind(ss[keepCols], kInfo$kTable[, c(1:3)]);
    ss.me = standardScreeningNumericTrait(MEs, status, corFnc = corFnc, corOptions = corOptions,
                                      getQvalues = FALSE, getDetails = TRUE,
                                      flatOutput = TRUE, addIDcolumn = TRUE)
    meInfo = data.frame( cor.Student.moduleEigenene.vs.diagnosis = ss.me$cor[ label2ME],
                         p.Student.moduleEigengene.vs.diagnosis = ss.me$p.Student[ label2ME],
                         enrichmentLabels[label2ME, ],
                         TopModuleHubs = topHubs[label2ME]);

  }


  names(out)[ names(out)=="kME.own" ] = "NumericModuleMembership";

  colnames(enrichmentLabels) = spaste("ModuleEnrichment.", colnames(enrichmentLabels));

  analysisInfo = data.frame(analysisName = rep(analysisName, nGenes),
                            dataSet = rep(setName, nGenes),
                            GEO.ID = rep(GEO.ID, nGenes));

  cbind(out, meInfo, analysisInfo);
}

#===================================================================================================
#
# Create a summary table for Mike Palazzolo's database - version 2, 
#
#===================================================================================================

# Submit gene ID, marginal association (fold change, FDR)
# module label, module color, position in module (kME), module annotation, module
# significance for status (correlation and Student p-value),
# preservation values

# My solution for now is to report all these statistics in each single set. For each analysis and each data
# set create a separate table.

PalazzoloTable.MM = function(multiExpr, labels, altLabels, multiStatus,
                          statusName,
                          corFnc.signif = cor, corOptions.signif = list(use = 'p'),
                          corFnc.kME = cor, corOptions.kME = list(use = 'p'),
                          dataIsLogTransformed = FALSE,
                          analysisName, setName, GEO.ID, 
                          enrichmentLabels,
                          modulePreservationTable,
                          includeStandardScreening,
                          includeAnalysisInfo,
                          ...)
{
  nGenes = checkSets(multiExpr)$nGenes;
  nSets = nSets(multiExpr);

  analysisName.x = spaste(analysisName, "Module");

  kInfo.x = mtd.apply(multiExpr, kIMandKME, labels, corFnc = corFnc.kME, corOptions = corOptions.kME, 
                    altLabels = altLabels, getKIM = FALSE, nTopHubs = 2, ...);

  topHubs = mtd.apply(kInfo.x, function(x) apply(x$topHubs, 2, paste, collapse = ", "));

  kInfo = mtd.apply(kInfo.x, getElement, "kTable")

  addNames = spaste(".in.", names(multiExpr));
  if (nSets > 1) addNames = spaste(".", analysisName.x, addNames);
  for (set in 1:nSets)
    colnames(kInfo[[set]]$data)[-c(1:2)] = spaste(colnames(kInfo[[set]]$data)[-c(1:2)], addNames[set]);

  MEs = multiSetMEs(multiExpr, universalColors = labels, excludeGrey = TRUE);
  nModules = checkSets(MEs)$nGenes;

  label2ME = match(labels, substring(mtd.colnames(MEs), 3));
  label2Pres = match(labels, modulePreservationTable$module);

  kInfoTable = cbind(kInfo[[1]]$data[, c(1:2)], 
          moduleEnrichmentSummary = enrichmentLabels$enrichmentLabel[label2ME]);
  kInfoTable = cbind(kInfoTable, do.call(cbind, multiData2list(mtd.subset(kInfo, , 3, drop = FALSE))));

  colnames(kInfoTable)[c(1:3)] = spaste(analysisName.x, c(".Label", ".Color", ".EnrichmentSummary")); 

  multiStatus = mtd.apply(multiStatus, as.matrix);
  nLevels = mtd.apply(multiStatus, function(x) nlevels(factor(x[, 1])), mdaSimplify = TRUE);
  #colnames(enrichmentLabels) = spaste("ModuleEnrichment.", colnames(enrichmentLabels));

  IDs = mtd.colnames(multiExpr);
  out = data.frame(ID = IDs);

  if (all(nLevels <3))
  {
    if (includeStandardScreening) 
    {
      printFlush("Calling standardScreeningBinaryTrait");
      ss = mtd.mapply(standardScreeningBinaryTrait, multiExpr, multiStatus,
                      traitNames = spaste(statusName, ".in.", names(multiExpr)),
                      MoreArgs = list(corFnc = corFnc.signif, corOptions = corOptions.signif,
                                      getQvalues = FALSE, getFDR = TRUE, 
                                      dataIsLogTransformed = dataIsLogTransformed,
                                      flatOutput = TRUE, addIDcolumn = FALSE,
                                      consistentSigns = TRUE));

      keepCols = c(5,4);
      out = cbind(out, do.call(cbind, multiData2list(mtd.subset(ss, , keepCols, drop = FALSE))));
    }

    out = cbind(out, kInfoTable);
    ss.me = mtd.mapply(standardScreeningBinaryTrait, MEs, multiStatus, 
                       traitNames = spaste(statusName, ".in.", names(multiExpr)),
                       MoreArgs = list(corFnc = corFnc.signif, corOptions = corOptions.signif,
                                      getQvalues = FALSE, getFDR = TRUE, dataIsLogTransformed = FALSE,
                                      flatOutput = TRUE, addIDcolumn = FALSE, consistentSigns = TRUE))

  } else {
    if (includeStandardScreening)
    {
      printFlush("Calling standardScreeningNumericTrait");
      ss = mtd.mapply(standardScreeningNumericTrait, multiExpr, multiStatus, 
                      traitNames = spaste(statusName, ".in.", names(multiExpr)),
                      MoreArgs = list(corFnc = corFnc.status, corOptions = corOptions.status,
                                         getQvalues = FALSE, getFDR = TRUE, getAreaUnderROC = FALSE, 
                                         flatOutput = TRUE, addIDcolumn = TRUE));
      keepCols = sapply(c("[a-z]*cor\\.", "FDR.Student"), grep, names(ss[[1]]$data))
      out = cbind(out, do.call(cbind, multiData2list(mtd.subset(ss, , keepCols, drop = FALSE))));
    }
    out = cbind(out, kInfoTable);

    ss.me = mtd.mapply(standardScreeningNumericTrait, MEs, multiStatus, 
                       traitNames = spaste(statusName, ".in.", names(multiExpr)),
                       MoreArgs = list(corFnc = corFnc.signif, corOptions = corOptions.signif,
                                      getQvalues = FALSE, getDetails = TRUE, getFDR = TRUE,
                                      flatOutput = TRUE, addIDcolumn = TRUE));

  }

  for (set in 1:nSets)
  {
     colnames(ss.me[[set]]$data) = gsub(".for.", 
         spaste(".of.", analysisName, ".ModuleEigengene.for."),
         colnames(ss.me[[set]]$data), fixed = TRUE);
     corCol = grep("^cor", colnames(ss.me[[set]]$data));
     colnames(ss.me[[set]]$data)[corCol] = gsub(".for.", ".with.", colnames(ss.me[[set]]$data)[corCol], 
                                                fixed = TRUE);
     colnames(ss.me[[set]]$data)[corCol] = gsub("^cor", "robustCor", colnames(ss.me[[set]]$data)[corCol]);
  }

  keepCols = sapply(c("^robustCor\\.", "^FDR\\.Student"), grep, colnames(ss.me[[1]]$data));
  names(ss.me) = NULL;
  meInfo = do.call(cbind, multiData2list(mtd.subset(ss.me, , keepCols)));

  names(out) = sub("kME.own", "NumericModuleMembership.in.", names(out));

  colnames(modulePreservationTable) = spaste(analysisName.x, ".preservation.in.",
                                             colnames(modulePreservationTable));


  out = cbind(out, meInfo[label2ME,], modulePreservationTable[label2Pres, -1]);

  if (includeAnalysisInfo)
  {
     analysisInfo = data.frame(analysisName = rep(analysisName, nGenes),
                               dataSet = rep(setName, nGenes),
                               GEO.ID = rep(GEO.ID, nGenes));
     out = cbind(out, analysisInfo);
  }
  out;
}

#=======================================================================================================
#
# restrictToReferenceRows
#
#=======================================================================================================

restrictToReferenceRows = function(target, targetNames = rownames(target), 
                                   reference, referenceNames = rownames(reference), warn = TRUE)
{
  ref2target = match(referenceNames, targetNames);
  if (warn && any(is.na(ref2target)))
  {
    show = which(is.na(ref2target));
    shorten = FALSE;
    if (length(show) > 50) {shorten = TRUE; show = show[1:50];}
    warning(immediate. = TRUE, 
            spaste("Some reference identifiers could not be matched in target:\n",
                   paste(referenceNames[show], collapse = ", "),
                   if (shorten) ", ... (truncated)" else "."));
  }

  target[ ref2target, ];
}
  

#=======================================================================================================
#
# loadOrNULL
#
#=======================================================================================================

# This is a special function to provide "existingResults" arguments for mtd.[m]apply... functions. Attempts
# to load a specified file that is assumed to contain a single object; if the file is successfully loaded,
# the object is returned; otherwise it returns NULL.

loadOrNULL = function(file)
{
  suppressWarnings(x <- try(load(file), silent = TRUE));
  if (inherits(x, "try-error")) return(NULL);
  get(x);
}

checkRDataFile = function(file, expectedObjects = NULL)
{
  suppressWarnings(x <- try(load(file, envir = parent.frame()), silent = TRUE));
  if (inherits(x, "try-error")) return(FALSE);

  if (!is.null(expectedObjects))
    if (!all(expectedObjects %in% x)) return(FALSE);

  TRUE;
}

 

#=======================================================================================================
#
# readEXPfiles
#
#=======================================================================================================

# Read Affy EXP files.

readEXPfiles = function(files, sep = "\t", row.names = files)
{
  nFiles = length(files);

  names = list();
  values = list();

  for (f in 1:nFiles)
  {
    conn = file(description = files[f], open = "rt");
    lines0 = readLines(conn);
    close(conn);

    keepLines = grep(sep, lines0);
    lines = lines0[keepLines];
    split= sapply(strsplit(lines, split = sep), identity);

    lengths = sapply(split, length);
    split = sapply( split[lengths==2], identity);

    names[[f]] = split[1, ]
    values[[f]] = split[2, ]
  }

  allNames = multiUnion(names);
  nNames = length(allNames);

  out = matrix("NA", nFiles, nNames);
  colnames(out) = allNames;
  rownames(out) = row.names;

  for (f in 1:nFiles)
    out[ f, match(names[[f]], allNames)] = values[[f]];

  out;
}

#=======================================================================================================
#
# readAgilentRawFiles
#
#=======================================================================================================

# Read Agilent raw files.

# These files contain several sections that are separated by lines containing a single "*". Some of these
# contain only one line of data (i.e., one of each entry per array, others have an entry for each probe.

# The problem is that with a large enough data set, returning all data may not be practical.

processLinesToTable = function(lines, header = TRUE, sep = "\t", keepColumns = NULL, excludeColumns = NULL)
{
  file = tempfile("tmp");
  con = file(file);
  writeLines(lines, con = con);
  close(con);

  data = read.csv.sql(file = file, header = header, sep = sep);
  if (length(keepColumns)==0)
  {
    keepColumns = c(1:ncol(data));
  } else {
    if (is.logical(keepColumns))
      keepColumns = which(keepColumns);
    if (!is.numeric(keepColumns))
    {
      keepColumns.num = match(keepColumns, names(data));
      if (any(is.na(keepColumns.num)))
        stop("Some entries in 'keepColumns' were not found in the read data:\n",
             paste(keepColumns[is.na(keepColumns.num)], collapse = ", "));
      keepColumns = keepColumns.num;
    }
    if (any(keepColumns <1 | keepColumns > ncol(data)))
        stop("Some entries in 'keepColumns' are out or range.");
  }

  if (length(excludeColumns)!=0)
  {
    if (is.logical(excludeColumns))
      excludeColumns = which(excludeColumns);
    if (!is.numeric(excludeColumns))
    {
      excludeColumns.num = match(excludeColumns, names(data));
      if (any(is.na(excludeColumns.num)))
        stop("Some entries in 'excludeColumns' were not found in the read data:\n",
             paste(excludeColumns[is.na(excludeColumns.num)], collapse = ", "));
      excludeColumns = excludeColumns.num;
    }
    if (any(excludeColumns <1 | excludeColumns > ncol(data)))
        stop("Some entries in 'excludeColumns' are out or range.");
  }

  keepColumns = keepColumns[ ! keepColumns %in% excludeColumns ];
     
  data[, keepColumns];
}

readAgilentRawFiles = function(files, sep = "\t", sepLine = "*", keepColumns = NULL,
                               excludeColumns = NULL)
{
  nFiles = length(files);
  decompFnc = rep("file", nFiles);
  decompFnc[grep(".gz$", files)] = "gzfile";
  decompFnc[grep(".bz2$", files)] = "bzfile";
  
  data = vector(mode = "list", length = nFiles);

  for (f in 1:nFiles)
  {
    printFlush(spaste("Working on file ", f, " of ", nFiles, "(", files[f], ")"));
    data[[f]] = list();
    fnc = match.fnc(decompFnc);
    x = try( { 
      conn = fnc(description = files[f], mode= "rt");
      lines = readLines(conn);
      close(conn);

      separators = (lines == sepLine);
      nSeparators = sum(separators);
      nLines = length(lines);

      lineRangeFrom = c(1, which(separators))+1;
      lineRangeTo = c(which(separators), nLines+1)-1;

      nSections = nSeparators + 1;
      for (s in 1:nSections)
        data[[f]] [[s]] = processLinesToTable(lines[lineRangeFrom[s]:lineRangeTo[s]]);
    } )
  }    
  data;
}

# Reorganize the set of agilent data into individual data frames that have one quantity per probe or one
# row per array.

reformatAgilentRawData = function(data)
{
  nSamples = length(data);
  nSections = length(data[[1]]);
}


#========================================================================================================
#
# listRep
#
#========================================================================================================

# Works like the rep function in its simplest case (a scalar n)  but returns a list instead of a vector

listRep = function(data, n)
{
  out = list();
  if (n> 0) for (i in 1:n) out[[i]] = data;
  out;
}


#========================================================================================================
#
# Functions for creating a set of reference labels
#
#========================================================================================================

weightedOverlapTable = function(labels1, labels2, weights1, weights2, na.rm = TRUE, ignore = NULL,
                        levels1 = NULL, levels2 = NULL)
{
  labels1 = as.vector(labels1);
  labels2 = as.vector(labels2);
  if (na.rm)
  {
    keep = !is.na(labels1) & !is.na(labels2) &!is.na(weights1) & !is.na(weights2);
    labels1 = labels1[keep];
    labels2 = labels2[keep];
    weights1 = weights1[keep];
    weights2 = weights2[keep];
  }
  if (is.null(levels1))
  {
    levels1 = sort(unique(labels1));
    levels1 = levels1[!levels1 %in%ignore];
  }
  if (is.null(levels2))
  {
    levels2 = sort(unique(labels2));
    levels2 = levels2[!levels2 %in%ignore];
  }
  n1 = length(levels1);
  n2 = length(levels2);
  countMat = pMat = weightMat = matrix(NA, n1, n2);

  for (m1 in 1:n1)
    for (m2 in 1:n2)
    {
      m1Members = (labels1 == levels1[m1]);
      m2Members = (labels2 == levels2[m2]);
      m12 = which(m1Members & m2Members);
      #print(paste("table for levels", levels1[m1], levels2[m2]));
      #print(table(m1Members, m2Members));
      pMat[m1, m2] = fisher.test(m1Members, m2Members, alternative = "greater")$p.value;
      countMat[m1, m2] = length(m12);
      if (length(m12) > 0) weightMat[m1, m2] = mean(weights1[m12] * weights2[m12]);
    }

 dimnames(pMat) = list(levels1, levels2);
 dimnames(countMat) = list(levels1, levels2);
 dimnames(weightMat) = list(levels1, levels2);

 pMat[is.na(pMat)] = 1;
 weightMat[is.na(weightMat)] = 0;

 list(countTable = countMat, pTable = pMat, weightTable = weightMat)
}


moduleSimilarity = function(labels, weights = NULL, ignoreLabels = 0, verbose = 1, indent = 0)
{
  labels = as.matrix(labels);
  spaces = indentSpaces(indent);
  nSets = ncol(labels);
  levels = apply(labels, 2, function(l1) {
                   levels.0 = sort(unique(l1));
                   levels.1 = levels.0[!is.na(levels.0)];
                   levels.1[!levels.1 %in% ignoreLabels];
                 });
  if (!is.null(dim(levels))) levels = as.list(as.data.frame(levels));

  if (is.null(weights)) weights = array(1, dim = dim(labels));
  weights = as.matrix(weights);

  lengths = sapply(levels, length);
  indexStart = c(1, cumsum(lengths)+1)[1:nSets];
  indexEnd = cumsum(lengths);
  nAll = sum(sapply(levels, length));
  overlaps = overlapWeights = matrix(0, nAll, nAll);
  overlapP = matrix(1, nAll, nAll);
  if (nSets <2) stop("Number of columns in 'labels' must be at least 2.");
  nComps = (nSets-1)*nSets/2;
  count = 0;
  if (verbose > 0) pind = initProgInd(spaces);
  for (s1 in 1:(nSets-1)) for (s2 in (s1+1):nSets)
  {
    ot = weightedOverlapTable(labels[, s1], labels[, s2], weights[, s1], weights[, s2], 
                              na.rm = TRUE, ignore = ignoreLabels, 
                              levels1 = levels[[s1]], levels2 = levels[[s2]]);
    overlaps[ indexStart[s1]:indexEnd[s1], indexStart[s2]:indexEnd[s2] ] = ot$countTable;
    overlaps[ indexStart[s2]:indexEnd[s2], indexStart[s1]:indexEnd[s1] ] = t(ot$countTable);
    overlapP[ indexStart[s1]:indexEnd[s1], indexStart[s2]:indexEnd[s2] ] = ot$pTable;
    overlapP[ indexStart[s2]:indexEnd[s2], indexStart[s1]:indexEnd[s1] ] = t(ot$pTable);
    overlapWeights[ indexStart[s1]:indexEnd[s1], indexStart[s2]:indexEnd[s2] ] = ot$weightTable;
    overlapWeights[ indexStart[s2]:indexEnd[s2], indexStart[s1]:indexEnd[s1] ] = t(ot$weightTable);
    count = count + 1;
    if (verbose > 0) pind = updateProgInd(count/nComps, pind);
  }
  if (verbose > 0) printFlush("");

  min = min(overlapP[overlapP > 0]);
  overlapP[overlapP < min] = min;

  similarity = -log10(overlapP) * overlapWeights;

  moduleNames = do.call(c, mapply(function(set, levels) spaste("Set", set, ".", levels),
                                  1:nSets, levels, SIMPLIFY = FALSE));

  dimnames(similarity) = list(moduleNames, moduleNames);

  similarity;
}


#=======================================================================================================
#
# Gene set similarity
#
#=======================================================================================================

# Gene set similarity for now simply defined as the average overlap fraction

geneSetSimilarity = function(identifiers1, identifiers2, assumeUnique = FALSE)
{
  if (!assumeUnique)
  {
    identifiers1 = unique(identifiers1);
    identifiers2 = unique(identifiers2);
  }
  
  nCommon = length(intersect(identifiers1, identifiers2));
  n1 = length(identifiers1);
  n2 = length(identifiers2);

  (nCommon/n1 + nCommon/n2)/2;
}

geneSetSimilarityMatrix = function(collection)
{
  genes = geneLists(collection);
  nSets = length(genes);
  similarity = matrix(1, nSets, nSets);
  for (s1 in 1:(nSets-1)) for (s2 in (s1+1):nSets)
    similarity[s1, s2] = similarity[s2, s1] = geneSetSimilarity(genes[[s1]], genes[[s2]], 
                                                      assumeUnique = TRUE);

  colnames(similarity) = rownames(similarity) = names(genes);

  similarity;
}

  

# Function to create reference labels. For each gene, record the clusters it belongs to plus the sum of
# weights of it belonging to the cluster (the weights can be, e.g., kME values but default to 1).
# Then select the cluster with highest weight and possibly only keep labels where the highest weight is in
# some sense "much more" than the second largest weight.

# it is assumed that 'clusters' is numeric and 0 means no cluster.

referenceLabels = function(labels, clusters, clusLabel2set, clusLabel2modLabel, 
                           weights = NULL,
                           useSets = 1:ncol(labels),
                           minWeightDifference = 1, minWeightRatio = 1,
                           minModuleSize = 20,
                           verbose = 1, indent = 0)
{

   spaces = indentSpaces(indent);

   nSets = ncol(labels);
   nGenes = nrow(labels);

   if (is.null(weights)) weights = array(1, dim = dim(labels));

   if (!isTRUE(all.equal(dim(labels), dim(weights)))) 
      stop("If 'weights' are supplied, they must have the same dimensions as 'labels'.");

   useSets.log = (1:ncol(labels)) %in% useSets;

   # Convert the cluster input into a list in which every component will give the map from the set labels
   # to clusters.

   label2cluster = list();
   for (set in 1:nSets)
   {
     index.set = which(clusLabel2set==set);
     labels.set = clusLabel2modLabel[index.set];
     label2cluster[[set]] = data.frame(label = labels.set, cluster = clusters[index.set]);
   }    

   clustersByGene = array(0, dim = dim(labels));
   for (set in 1:nSets)
   {
     lab2cl = match(labels[, set], label2cluster[[set]]$label)
     fin = is.finite(lab2cl);
     clustersByGene[fin, set] = label2cluster[[set]]$cluster[ lab2cl[fin] ];
   }

   # Calculate weights for gene in each cluster each cluster and pick the reference labels

   minWR.inv = 1/minWeightRatio;
   refLabels.0 = rep(0, nGenes);
   pind = initProgInd(spaces);
   for (g in 1:nGenes)
   {
     c1 = clustersByGene[g, ];
     nz = c1!=0 & useSets.log;
     if (sum(nz) > 0)
     {
       clWts = tapply(weights[g, nz], c1[nz], sum, na.rm = TRUE);
       order = order(-clWts, na.last = TRUE)
       if (length(clWts) == 1 || (
         clWts[order[1]] - clWts[order[2]] >= minWeightDifference &&
              clWts[order[2]]/clWts[order[1]] <= minWR.inv) )
            refLabels.0[g] = as.numeric(names(clWts)[order[1]]);
     }
     if (verbose > 0 && g%%100==0) pind = updateProgInd(g/nGenes, pind);
   }
   if (verbose > 0) { pind = updateProgInd(1, pind); printFlush("")}

   # Normalize the reference labels and remove reference modules that are too small.

   refLabels.1 = normalizeLabels(refLabels.0, keepZero = TRUE);

   t = table(refLabels.1);
   delete = as.numeric(names(t)[t < minModuleSize]);

   refLabels.1[ refLabels.1 %in% delete ] = 0;

   refLabels.1;
}

splitLabelsByOrganism = function(entrez, labels, org.entrez = NULL, organisms, stopOnMissingOrg = TRUE,
                                 useHomology = TRUE)
{
  nOrganisms = length(organisms);
  nGenes = length(entrez);
  if (is.null(org.entrez))
  {
    org.entrez = rep(NA, nGenes);
    for (o in 1:nOrganisms)
    {
       orgGenes = allOrgGenes(organisms[o]);
       inOrg = entrez %in% orgGenes;
       org.entrez[inOrg] = o;
    }
  } else org.entrez = match(org.entrez, organisms)

  if (any(is.na(org.entrez)) && stopOnMissingOrg) 
    stop("Some 'entrez' could not be found in any of the given organisms.");

  out = list()
  for (o in 1:nOrganisms)
  {
    entrez.o = rep(NA, nGenes)
    direct = org.entrez==o;
    entrez.o[direct] = entrez[direct];
    mapped = direct
    for (o2 in (1:nOrganisms)[-o])
    {
      try({
          entrez.21 = mapEntrez(entrez[!mapped], orgFrom = organisms[o2], orgTo = organisms[o], 
                                useHomology = useHomology);
          fin = !is.na(entrez.21)
          entrez.o[!mapped ]  = entrez.21;
          mapped[!mapped] = fin;
         }, silent = TRUE);
    }
    out[[o]] = data.frame(Entrez = entrez.o[mapped], Label = labels[mapped]);
  }
  
  names(out) = organisms;
  out;
}


# Note: sourceLabels may be a matrix.

matchLabelsToReference = function(sourceEntrez, sourceLabels,
                                  targetEntrez, targetLabels, pThreshold = 5e-2,
                                  ignoreLabels = if (is.numeric(targetLabels)) 0 else "grey",
                                  extraLabels = if (is.numeric(targetLabels)) c(1:1000) else standardColors())
{

  targetLabels.src = replaceMissing(targetLabels[ match(sourceEntrez, targetEntrez)]);

  matchLabels(sourceLabels, targetLabels.src, pThreshold = pThreshold, 
              na.rm = TRUE,
              ignoreLabels = ignoreLabels,
              extraLabels = extraLabels);
}


  
matchLabelsToReference.2 = function(sourceLabels, targetLabels, IDname = "Entrez", 
                                    labelName = "Label", ...)
{
  IDCol.src = match(IDname, names(sourceLabels));
  IDCol.tgt = match(IDname, names(targetLabels));
  
  labelCol.src = match(labelName, names(sourceLabels));
  labelCol.tgt = match(labelName, names(targetLabels));
  
  matchLabelsToReference(sourceLabels[, IDCol.src], sourceLabels[, labelCol.src],
                         targetLabels[, IDCol.tgt], targetLabels[, labelCol.tgt], ...)
}


#========================================================================================================
#
# splitData
#
#========================================================================================================

# Split a data frame or a matrix by rows into a multiData structure.

splitData = function(data, split, nameBase = "")
{
  groups = levels(factor(split));
  nGroups = length(groups);

  splitData = list();
  for (g in 1:nGroups)
  {
     splitData[[g]] = list(data = data[split==groups[g], ]);
  }

  names(splitData) = spaste(nameBase, groups);
  splitData;
}


#=======================================================================================================
#
# selectByListIndex
#
#=======================================================================================================

selectByListIndex = function(data, indexList)
{
  lapply(indexList, function(index, data) data[index], data)
}

#========================================================================================================
#
# Wrappers for normalize.quantiles...
#
#========================================================================================================

normalize.quantiles.PL = function(data)
{
  out = t(normalize.quantiles(t(data)));
  dimnames(out) = dimnames(data);
  out;
}

normalize.quantiles.robust.PL = function(data, ...)
{
  out = t(normalize.quantiles.robust(t(data), ...));
  dimnames(out) = dimnames(data);
  out;
}

#========================================================================================================
#
# moduleLabelsFromHubs
#
#========================================================================================================

moduleLabelsFromHubs = function(expr, labels, nHubs, sep = ", ", ...)
{
  if (is.data.frame(expr) || is.matrix(expr))
  {
    multiSet = FALSE;
    expr = multiData(expr);
  } else
    multiSet = TRUE

  args = list(...)
  lengths = sapply(args, length);
  nSets = length(expr);

  adj = mtd.apply(expr, adjacency, ...);

  kim = do.call(cbind, mtd.apply(adj, function(...) intramodularConnectivity(..., scaleByMax = TRUE)[, 2], labels,
                                 returnList = TRUE));

  if (is.numeric(labels))
  {
     modules.0 = sort(as.numeric(unique(labels)));
     modules = modules.0[modules.0!=0];
  } else {
     modules.0 = sort(unique(modules));
     modules = modules.0[modules.0!="grey"];
  }

  nModules = length(modules);
  moduleLabels = character(nModules);

  for (m in 1:nModules)
  {
    kim.mod = rowSums(kim[labels==modules[m], , drop = FALSE], na.rm = TRUE);
    genes.mod = mtd.colnames(expr)[ labels==modules[m]];
    keep = which(rank(-kim.mod) <= nHubs);
    moduleLabels[m] = paste(genes.mod[keep], collapse = sep);
  }
  t = table(labels);
  t = t[match(modules, names(t))];
  moduleLabels = spaste("M.", modules, " (", as.numeric(t), "): ", moduleLabels);
  data.frame(module = modules, label = moduleLabels);
}


#=======================================================================================================
#
# Rearrange multiData to matrices by columns
#
#=======================================================================================================

multiData2matrixByColumns = function(multiData, setNames = names(multiData))
{
  size = checkSets(multiData);
  nCols = size$nGenes;

  if (!all(size$nSamples==size$nSamples[1])) 
    stop("This function requires that number of samples in each set be the same.");

  out = lapply(1:nCols, function(index, data) mtd.apply(data, function(x1, i) x1[, i], i=index, 
                        mdaSimplify = TRUE), multiData);

  names(out) = mtd.colnames(multiData);

  if (!is.null(setNames)) for (c in 1:nCols) 
     colnames(out[[c]]) = spaste(mtd.colnames(multiData)[c], ".", setNames);
  
  out;

}

#=======================================================================================================
#
# interleave (or cbind with same colname modifications)
#
#=======================================================================================================

interleave = function(matrices, nameBase = names(matrices), sep = ".", baseFirst = TRUE,
                      doInterleave = TRUE, check.names = TRUE)
{
  nMats = length(matrices)
  nCols = ncol(matrices[[1]]);

  dims = lapply(matrices, dim);

  if (any(nameBase!=""))
  {
    if (baseFirst)
    {
       for (m in 1:nMats) colnames(matrices[[m]]) = spaste(nameBase[m], sep, colnames(matrices[[m]]));
    } else {
       for (m in 1:nMats) colnames(matrices[[m]]) = spaste(colnames(matrices[[m]]), sep, nameBase[m]);
    }
  }
  matrices = lapply(matrices, function(mat) if (is.list(mat)) as.data.frame(mat) else mat)

  if (doInterleave)
  {
    out = as.data.frame(lapply(1:nCols, 
                               function(index, matrices) 
                                  as.data.frame(lapply(matrices, 
                                            function(x, i) x[, i, drop = FALSE], index), check.names = FALSE),
                               matrices), check.names = check.names);
  } else
    out = as.data.frame(do.call(cbind, removeListNames(matrices)));

  if (!is.null(rownames(matrices[[1]]))) rownames(out) = make.unique(rownames(matrices[[1]]))
  out;
} 

interleave.withoutNames = function(matrices, doInterleave = TRUE, check.names = TRUE)
{
  interleave(matrices, nameBase = rep("", length(matrices)), sep = "", doInterleave = doInterleave,
            check.names = check.names);
}
  
interleaveLists = function(lists, nameBase = names(lists), sep = ".", baseFirst = TRUE,
                      doInterleave = TRUE, check.names = TRUE)
{
  nMats = length(lists)
  nCols = length(lists[[1]]);
  if (any(nameBase!=""))
  {
    if (baseFirst)
    {
       for (m in 1:nMats) names(lists[[m]]) = spaste(nameBase[m], sep, names(lists[[m]]));
    } else {
       for (m in 1:nMats) names(lists[[m]]) = spaste(names(lists[[m]]), sep, nameBase[m]);
    }
  }
  lists = removeListNames(lists);
  if (doInterleave)
  {
    out = do.call(c, lapply(1:nCols, 
                       function(index, lists) do.call(c, lapply(lists, function(x, i) x[i], index)),
                       lists));
  } else
    out = do.call(c, lists);

  out;
} 

interleaveLists.withoutNames = function(lists, doInterleave = TRUE, check.names = TRUE)
{
  interleaveLists(lists, nameBase = rep("", length(lists)), sep = "", doInterleave = doInterleave,
            check.names = check.names);
}




interleave.gapped = function(matrices, nameBase = names(matrices), sep = ".", baseFirst = TRUE,
                             removePatternBeforeMatch = "", fixed = FALSE)
{
  nMats = length(matrices)
  nCols = ncol(matrices[[1]]);

  dims = lapply(matrices, dim);
  nAllCols = sum(sapply(matrices, ncol));
  if (any(sapply(matrices, nrow)!=nrow(matrices[[1]]))) 
     stop("All 'matrices' must have the same number of rows.");

  colnames1 = lapply(matrices, colnames);
  if (removePatternBeforeMatch!="")
  {
     colnames.match = lapply(colnames1, sub, removePatternBeforeMatch, "", fixed = fixed);
  } else
     colnames.match = colnames1;

  if (baseFirst)
  {
     for (m in 1:nMats) colnames(matrices[[m]]) = spaste(nameBase[m], sep, colnames(matrices[[m]]));
  } else {
     for (m in 1:nMats) colnames(matrices[[m]]) = spaste(colnames(matrices[[m]]), sep, nameBase[m]);
  }


  cn.duplicated = lapply(colnames.match, duplicated);
  lapply(cn.duplicated, function(dup) stopifnot(all(!dup)));

  out  = list();
  taken = lapply(matrices, function(x) rep(FALSE, ncol(x)));

  for (m in 1:nMats)
  {
    for (i in c(1:ncol(matrices[[m]]))[!taken[[m]] ])
    {
       n1 = colnames.match[[m]] [i]
       tab1 = matrices[[m]] [, i, drop = FALSE];
       if (m < nMats) for (j in (m+1):nMats)
       {
          col = which(!taken[[j]]) [match(n1, colnames.match[[j]] [!taken[[j]] ])];
          if (!is.na(col)) 
          {
            tab1 = cbind(tab1, matrices[[j]] [, col, drop = FALSE]);
            taken[[j]] [col] = TRUE;
          }
       }
       taken[[m]] [i] = TRUE;
       out = c(out, list(tab1));
    }
  }
  do.call(cbind, out);
}

  
        

#=======================================================================================================
#
# Plot dendrogram
#
#=======================================================================================================

plotDendrogram = function(tree, labels = NULL, 
                          horiz = FALSE, reverseDirection = FALSE,
                          hang = 0.1, xlab = "", ylab = "", 
                          cex.labels = 1, ..., adjustRange = FALSE,
                          shift.first = 0, shift.last = 0,
                          labelsAngle = if (horiz) 0 else 90)
{
  hang.gr = hang;
  if (hang < 0) hang.gr = 0.1;
  n = length(tree$order);
  heights = tree$height;
  range = range(heights);
  hang.scaled = hang.gr * (max(heights) - min(heights));
  range[1] = range[1] - hang.scaled;

  indexLim = c(0.5, n+0.5);
  if (adjustRange)
  {
    ctr = mean(indexLim);
    indexLim = ctr + (indexLim - ctr)/1.08;
  }
  nMerge = n-1;
  if (is.null(labels)) labels = tree$labels;
  if (is.null(labels)) labels = rep("", n);
  if (is.na(labels[1])) labels = rep("", n);
  if (is.logical(labels) && labels[1]=="FALSE") labels = rep("", n);

  if (horiz) 
  {
    plot(NA, NA, xlim = if (reverseDirection) range else rev(range), ylim = indexLim,
         axes = TRUE, yaxt = "none", frame = FALSE, type ="n", xlab = xlab, ylab = ylab, ...);
  } else {
    plot(NA, NA, ylim = if (reverseDirection) rev(range) else range, xlim = indexLim,
         axes = TRUE, xaxt = "none", frame = FALSE, type ="n", xlab = xlab, ylab = ylab, ...);
  }

  pin = par("pin");
  box = par("usr");
  factor.usr = if (horiz) (box[4] - box[3])/pin[2] else (box[2] - box[1])/pin[1];
  shift.first.usr = shift.first * factor.usr;
  shift.last.usr = shift.last * factor.usr;

  space = n -shift.first.usr - shift.last.usr;
  singleton.x = rep(NA, n);
  singleton.x[tree$order] = seq(from = 1+shift.first.usr, to = n-shift.last.usr, length.out = n);
  cluster.x = rep(NA, n);

  for (m in 1:nMerge)
  {
     o1 = tree$merge[m, 1]
     o2 = tree$merge[m, 2]
     h = heights[m];
     hh = if (hang>0) h-hang.scaled else range[1];
     h1 = if (o1 < 0) hh else heights[o1];
     h2 = if (o2 < 0) hh else heights[o2];

     x1 = if (o1 < 0) singleton.x[-o1] else cluster.x[o1]
     x2 = if (o2 < 0) singleton.x[-o2] else cluster.x[o2]

     cluster.x[m] = mean(c(x1, x2));

     if (!is.null(labels))
     {
       if (horiz)
       {
          if (o1 < 0) text(h1, x1, spaste(labels[-o1], " "), adj = c(0, 0.5), srt = labelsAngle,
                        cex = cex.labels, xpd = TRUE)
          if (o2 < 0) text(h2, x2, spaste(labels[-o2], " "), adj = c(0, 0.5), srt = labelsAngle,
                        cex = cex.labels, xpd = TRUE)
       } else {
          if (o1 < 0) text(x1, h1, spaste(labels[-o1], " "), adj = c(1, 0.5), srt = labelsAngle,
                        cex = cex.labels, xpd = TRUE)
          if (o2 < 0) text(x2, h2, spaste(labels[-o2], " "), adj = c(1, 0.5), srt = labelsAngle,
                        cex = cex.labels, xpd = TRUE)
       }
     }
  
     if (horiz)
     {
       lines(c(h1, h, h, h2), c(x1, x1, x2, x2));
     } else {
       lines(c(x1, x1, x2, x2), c(h1, h, h, h2));
     }
  }
}

#=======================================================================================================
#
# pValueCodes: concise indication of p-values
#
#=======================================================================================================

pValueCodes = function(p, lead = "|", trail = "", digits = 1)
{
  exponent = round(-log10(p));
  exponent[exponent>10^digits-1] = 10^digits-1;

  spaste(lead, exponent, trail);
}

pValueStars = function(p, thresholds = c(0.05, 0.01, 0.001))
{
  ifelse(p>thresholds[1], "",
      ifelse(p>thresholds[2], "*",
      ifelse(p>thresholds[2], "**", "***")));
}

pValueStars2 = function(p, char = "*", start = 1e-1, max = 5)
{
  p = replaceMissing(p, 1);
  n = floor(-log10(p/start)) + 1;
  n[n < 0] = 0;
  n[ p>0.05] = 0;
  n[n>max] = max;
  out = sapply(n, function(n1) paste(rep(char, n1), collapse = ""));
  dim(out) = dim(p);
  dimnames(out) = dimnames(p);
  out;
}



#======================================================================================================
#
# multiMerge
#
#======================================================================================================

multiMerge = function(argList, ...)
{
    len = length(argList)
    if (len == 0) 
        return(NULL)
    if (len == 1) 
        return(argList[[1]])
    out = argList[[1]]
    for (elem in 2:len) out = merge(out, argList[[elem]], ...)
    out
}

#======================================================================================================
#
# Split strings by split, put into a matrix filling empty values by NAs
#
#======================================================================================================

strsplit2matrix = function(strings, ...)
{
  split = strsplit(strings, ...);
  lens = sapply(split, length);
  max = max(lens);
  t(sapply(split, function(x, l) {if (length(x)<l) x = c(x, rep(NA, l-length(x))); x}, max));
}

removeLeadingX = function(x) sub("^X", "", x);

#=====================================================================================================
#
# Combined mymapply(mtd.mapply, ...)
#
#=====================================================================================================

mymapply = function(..., SIMPLIFY = FALSE) base::mapply(..., SIMPLIFY = SIMPLIFY)

mmmapply = function(FUN, ..., MoreArgs = NULL, MoreArgs.out = list()) 
   mymapply(mtd.mapply, ..., MoreArgs = c(list(FUN=FUN, MoreArgs = MoreArgs), MoreArgs.out));

#====================================================================================================
#
#
#
#====================================================================================================

as.Date = function(x, ...)
{
  if (is.numeric(x)) base::as.Date(x, origin = "1970-1-1") else base::as.Date(x, ...)
}

#=============================================================================================================
#
# drop repeated columns
#
#=============================================================================================================

dropRepeatedColumns = function(x, exceptions = character(0))
{
  nCols= ncol(x);
  ref = 1;
  remove = rep(FALSE, nCols);
  keep = multiGrepl(exceptions, colnames(x));
  while (ref < nCols)
  {
    if (!remove[ref]) 
      for (col in (ref+1):nCols) if (!keep[col])
      {
         if (isTRUE(all.equal(x[, ref], x[, col]))) remove[col] = TRUE;
         if (!remove[col] && is.numeric(x[, ref]) && is.numeric(x[, col]) && all(is.na(x[,ref])==is.na(x[, col])) &&
            all(replaceMissing(scale(x[, col]))==replaceMissing(scale(x[, ref])))) remove[col] = TRUE;
      }
    ref = ref + 1;
  }
  if (any(remove) > 0) x = x[, !remove];
  attr(x, "droppedColumnIndicator") = remove;
  x;
}

#=============================================================================================================
#
# dropConstantColumns
#
#====================================================================================================


dropConstantColumns = function(data, checkDuplicates = TRUE)
{
  keep = apply(data, 2, function(x) {x = x[!is.na(x)]; if (length(x)==0) FALSE else any(x!=x[1])});
  out = data[, keep, drop = FALSE];
  if (checkDuplicates)
  {
    keep2 = !duplicated(as.data.frame(t(out)));
    out = out[, keep2, drop = FALSE];
  }
  out;
}

#====================================================================================================
#
# dropDependentColumns
#
#====================================================================================================

# Assumes that none of the columns are all zeros. Assumes no missing data.

dependentColumns = function(data, testOrder = c(1:ncol(data)), SDthreshold = 0,
                            reverseTestOrder = FALSE, returnColumnNames = FALSE,
                            removeConstant = TRUE)
{
  if (any(is.na(data))) stop("This function does not work with missing data yet.");

  if (reverseTestOrder) testOrder = rev(testOrder);

  nc.all = ncol(data);
  keep.all = rep(FALSE, nc.all);
  cNames.all = colnames(data)

  if (removeConstant)
  {
    sds = apply(data, 2, sd, na.rm = TRUE)
    sds[is.na(sds)] = 0;
    drop = sds<=SDthreshold;
    data = data[, !drop, drop = FALSE];
  } else 
    drop = rep(FALSE, nc.all)

  index = 1:sum(!drop);
  testOrder = testOrder[ !drop];
  order = order(testOrder);
  testOrder[order] = index;

  if (is.data.frame(data))
  {
    data = as.matrix(data);
    convertToDF = TRUE;
  } else convertToDF = FALSE;
  nc = ncol(data);
  keep = rep(FALSE, nc);
  keep[testOrder[1]] = TRUE;
  ref = sum(data[, keep]^2);
  cNames = colnames(data)
  if (is.null(cNames)) cNames = rep("", nc);
  for (cc in testOrder[-1])
  {
    candidates = keep;
    candidates[cc] = TRUE;
    data1 = data[, candidates, drop = FALSE];
    d1= as.matrix(data1);
    rnk = qr(d1)$rank;
    if (rnk==ncol(d1)) keep[cc] = TRUE else
      printFlush(spaste("Dropping column ", cc, " (", cNames[cc], ")"));
  }
  keep.all[!drop] = keep;
  if (returnColumnNames) cNames.all[!keep.all] else !keep.all;
}


dropDependentColumns = function(data, testOrder = c(1:ncol(data)), SDthreshold = 0,
                                reverseTestOrder = FALSE, removeConstant = TRUE)
{
  depCols = dependentColumns(
      data = data,
      testOrder = testOrder, 
      SDthreshold = SDthreshold,
      reverseTestOrder = reverseTestOrder,
      returnColumnNames = FALSE, removeConstant = removeConstant)
  out = data[, !depCols, drop = FALSE]
  attr(out, "columnsRetainedFromOriginalData") = !depCols;
  out;
}

# A simpler version that simply drops duplicated columns.

dropDuplicatedColumns = function(data, testOrder = c(1:ncol(data)), reverseTestOrder = FALSE,
                                 scale = FALSE, matchColnames = FALSE)
{
  if (reverseTestOrder) testOrder = rev(testOrder);
  if (scale) data.test = scale(data) else data.test = data;
  nc = ncol(data);
  if (is.null(colnames(data))) colnames(data) = spaste("Col.", 1:nc);
  keep = rep(FALSE, nc);
  keep[testOrder[1]] = TRUE;
  cNames = colnames(data);
  if (is.null(cNames)) cNames = rep("", nc);
  for (cc in testOrder[-1])
  {
    candidates = keep;
    if (matchColnames) candidates = candidates & colnames(data.test)==colnames(data.test)[cc];
    if (any(candidates)) {
        keep[cc] = !any(apply(data.test[, candidates, drop = FALSE], 2, 
            function(x, y) all(x==y, na.rm = TRUE), data.test[, cc]));
    } else keep[cc] = TRUE;
  }
  out = data[, keep];
  attr(out, "columnsRetainedFromOriginalData") = keep;
  out;
}




#===================================================================================================
#
# Linear-log transformation
#
#===================================================================================================

linLogTrafo = function(data)
{
  out = data/2;
  out[ data > 2] = log2(data[data > 2]);
  out;
}

#===================================================================================================
#
# Module membership in own module
#
#===================================================================================================

consensusModuleMembership.own = function(labels, conKME, pattern = "Z.ModMembership.in.M.")
{
    nGenes = length(labels)
    label2column = match(spaste(pattern, labels), colnames(conKME))
    index = cbind(1:nGenes, label2column);
    conKME[index];
}




#=================================================================================================
#
# setPVE
#
#=================================================================================================

setPVE = function(data, x, corFnc = "cor", corOptions = list(use = 'p'), p = 0.95,
                  nPermutations = 0, randomSeed = 12345)
{
  corFnc = match.fun(corFnc);

  cr.obs = do.call(corFnc, c(list(x = data, y = x), corOptions))^2;
  average = colMeans(cr.obs, na.rm = TRUE);

  q = colQuantileC(cr.obs, p = p);

  observed = data.frame(average = average, quantile = q);
  rownames(observed) = colnames(x);

  if (nPermutations > 0)
  {
    if (!is.null(randomSeed)) set.seed(randomSeed);
    permuted = array(NA, dim = c(nPermutations, dim(cr.obs)));
    pind = initProgInd();
    nSamples = nrow(data);
    for (perm in 1:nPermutations)
    {
      sampleIndex = sample(1:nSamples);
      permuted[perm, , ] = do.call(corFnc, c(list(x = data[sampleIndex, ], y = x), corOptions))^2;
      pind = updateProgInd(perm/nPermutations, pind);
    }
    printFlush("");

    mean.perm = colMeans(permuted, dims = 1, na.rm =TRUE);
    sd.perm = apply(permuted, c(2,3), sd, na.rm = TRUE);

    Z.perm = (cr.obs - mean.perm)/sd.perm;
    Z.ave = colMeans(Z.perm, na.rm = TRUE)
    Z.q = colQuantileC(Z.perm,p = p);
    

    Z.perm = data.frame(Z.average = Z.ave, Z.quantile = Z.q);
    rownames(Z.perm) = colnames(x);

    list(observed = observed, Z = Z.perm);
  } else {
    observed;
  }
}

#=================================================================================================
#
# exportCausalNetworkToCytoscape
#
#=================================================================================================

exportCausalNetworkToCytoscape = function(fromNode, toNode,
        edgeWeight,
        edgeFile = NULL,
        nodeFile = NULL,
        weighted = TRUE,
        threshold = 0.5,
        nodeAttr = NULL,
        nodeNameColumn = "node",
        includeColNames = TRUE)
{
    edgeData = data.frame(
        fromNode = fromNode,
        toNode = toNode,
        weight = edgeWeight,
        direction = 1);
    nodes = unique(c(fromNode, toNode));
    nNodes = length(nodes)
    if (is.null(nodeAttr))
    {
       nodeAttr = data.frame(node = nodes);
    } else
       nodeAttr = nodeAttr[ match(nodes, nodeAttr[, match(nodeNameColumn, colnames(nodeAttr))]), ];

    if (!is.null(edgeFile)) 
        write.table(edgeData, file = edgeFile, quote = FALSE, 
            row.names = FALSE, col.names = includeColNames, sep = "\t")
    if (!is.null(nodeFile)) 
        write.table(nodeAttr, file = nodeFile, quote = FALSE, 
            row.names = FALSE, col.names = includeColNames, sep = "\t")
    list(edgeData = edgeData, nodeAttr = nodeAttr)

}


#=================================================================================================
#
# quadratic allocation of jobs to workers
#
#=================================================================================================

# Assume that the length of task i is proportional to i.

allocateJobs.quadratic = function(nTasks, nWorkers, reverseOrder = TRUE)
{
  nLeft = nTasks;
  allocation = vector(mode = "list", length = nWorkers);
  for (w in 1:nWorkers)
  {
    leftPerWorker= nLeft*(nLeft+1)/2/(nWorkers - w + 1);
    a = min(round(sqrt( 0.25 + nLeft * (nLeft+1) - 2*leftPerWorker) -0.5), nLeft-1);
    allocation[[w]] = if (reverseOrder) (nTasks + 1 -c(nLeft:(a+1))) else c((a+1):nLeft);
    nLeft = a;
  }
  allocation;
}

#=================================================================================================
#
# Parametric and non-parametric Permutation test p-values
#
#=================================================================================================

# It is assumed that dim(permuted) = c(nPermutations, dim(observed))

parametricPermutationPValues = function(observed, permuted, alternative = c("two.sided", "less", "greater"))
{
  ia = match.arg(alternative)
  do = dim(observed)
  if (is.null(do)) do = length(observed);

  if (!isTRUE(all.equal(do, dim(permuted)[-1]))) 
     stop("Dimensions of 'observed' and 'permuted' are not consistent.");

  means = colMeans(permuted, na.rm = TRUE);
  sds = apply(permuted, c(2:length(dim(permuted))), sd, na.rm = TRUE);

  Z = (observed-means)/sds;

  if (ia=="less")
  {
    p = pnorm(Z, lower.tail = TRUE)
  } else if (ia=="greater") {
    p = pnorm(Z, lower.tail = FALSE)
  } else 
    p = 2*pnorm(abs(Z), lower.tail = FALSE);

  list(p = p, Z = Z);
}

permutationPValues = function(observed, permuted, alternative = c("two.sided", "less", "greater"))
                              #getFDR = TRUE, scaleBy = c("none", "absMean", "absMedian"))
{
  ia = match.arg(alternative)
  do = dim(observed)
  if (is.null(do)) do = length(observed);

  if (!isTRUE(all.equal(do, dim(permuted)[-1])))
     stop("Dimensions of 'observed' and 'permuted' are not consistent.");

  permuted.x = permuted;
  dim(permuted.x) = c(dim(permuted)[1], prod(do));
  observed.x = matrix(as.numeric(observed), dim(permuted)[1], prod(do), byrow = TRUE);

  valueCounts = colSums(is.finite(permuted.x));

  if (ia=="two.sided")
  {
    counts = colSums(abs(observed.x) < abs(permuted.x), na.rm = TRUE); 
  } else if (ia=="less") {
    counts = colSums(observed.x > permuted.x, na.rm = TRUE)
  } else
    counts = colSums(observed.x < permuted.x, na.rm = TRUE);

  p = (counts+1)/(valueCounts + 1);
  dim(p) = dim(observed);

  if (length(do)==1) names(p) = names(observed) else dimnames(p) = dimnames(observed);
  p;
}

permutationThresholds = function(permuted, alternative = c("two.sided", "less", "greater"),
    prob = c(0.05))
{
  ia = match.arg(alternative)

  dp = dim(permuted);
  if (is.null(dp)) dp = c(length(permuted), 1);
  dp1 = dp[-1];

  dn = dimnames(permuted);
  if (is.null(dn)) dn = lapply(dp, function(i) spaste("X", 1:i));

  nThresholds = length(prob);
  if (any(!is.finite(prob))) stop("All entries in 'prob' must be finite.");
  if (any(prob<0 | prob > 1)) stop("All 'prob' must be between 0 and 1.");
  permuted.x = permuted;
  dim(permuted.x) = c(dp[1], prod(dp1));

  valueCounts = colSums(is.finite(permuted.x));

  if (ia=="two.sided")
  {
    out = do.call(cbind, lapply(prob, function(p) colQuantileC(abs(permuted.x), 1-p)));
  } else if (ia=="less") {
    out = do.call(cbind, lapply(prob, function(p) colQuantileC(permuted.x, p)));
  } else
    out = do.call(cbind, lapply(prob, function(p) colQuantileC(permuted.x, 1-p)));

  dim(out) = c(dp1, nThresholds);
  dimnames(out) = c(dn[-1], list(prob));
  out;
}



#=================================================================================================
#
# Permutation test 
#
#=================================================================================================

permutationTest = function(data, pheno,
                      associationFnc,
                      observedStats = NULL,
                      statPatterns = NULL,
                      fixed = FALSE,
                      nPermutations, 
                      randomSeed = 12345,
                      ...,
                      ptVerbose = 1,
                      getPermutedStats = TRUE)
{
  assocFnc = match.fun(associationFnc);

  if (is.null(observedStats))
  {
    observedStats = assocFnc(data, pheno, ...);
    #observedStats = assocFnc(data, pheno, flatOutput= TRUE);
    if (length(observedStats)==0) stop("The return value of 'associationFnc' is empty.");
  }

  if (length(dim(observedStats)) > 2) 
     stop("This function only works with at most two-dimensional 'observedStats'.\n",
          "The restriction also applies to output of 'associationFnc'.");

  if (is.null(dim(observedStats))) observedStats  = as.matrix(observedStats);
  if (!is.null(statPatterns)) 
  {
    statCols= multiGrep(statPatterns, colnames(observedStats),
                        fixed = fixed, sort = TRUE)
    if (length(statCols)==0) 
      stop("No columns in 'observedStats' or output of 'associationFnc' match 'statPatterns'.");

    observedStats = observedStats[, statCols, drop = FALSE]
  } else {
    statCols = 1:ncol(observedStats);
  }

  observedStats  = as.matrix(observedStats);

  nStats = ncol(observedStats);

  nFeatures = nrow(observedStats);

  if (is.null(dim(pheno))) pheno = as.matrix(pheno);
  nSamples = nrow(pheno);

  permStats = array(NA, dim = c(nPermutations, nFeatures, nStats));
  samplePermutations = matrix(0, nSamples, nPermutations);

  if (!is.null(randomSeed)) { set.seed(randomSeed); seedSaved = TRUE; }
  if (ptVerbose==1) pind = initProgInd();
  for (p in 1:nPermutations)
  {
    if (ptVerbose>1) printFlush(paste("Working on permutation", p));
    sample = sample(nSamples);
    samplePermutations[, p] = sample;
    permStats[p, , ] = as.matrix(assocFnc(data, pheno[sample, ,drop = FALSE], ...)[, statCols]);
    if (ptVerbose==1) pind = updateProgInd(p/nPermutations, pind);
  }
  out = list(parametric = parametricPermutationPValues(observedStats, permStats),
             nonParametric = permutationPValues(observedStats, permStats),
             samples = samplePermutations,
             permutedStatistics = if (getPermutedStats) permStats else NULL);
  out;
}



# Utility function for permutation tests of correlations
    
corZvalue = function(cpFnc = corAndPvalue, ...)
{
  cpFnc = match.fun(cpFnc);
  cpFnc(...)$Z;
}
 

#==================================================================================================
#
# twoSideBarplot
#
#==================================================================================================

textBox = function(x, y, text, adj, cex = 1, adjustBottom = TRUE, adjAmount = 0.3, ...)
{
  width = strwidth(text, cex = cex, ...)
  height = strheight(text, cex = cex, ...)
  extend = grepl("q|g|j|y", text);
  x1 = x-adj[1] * width;
  y1 = y-adj[2] * height;
  if (extend && adjustBottom)
  {
    y1 = y1-height * adjAmount;
    height = height * (1+adjAmount);
  }
  c(xmin = x1, ymin = y1, xmax = x1 + width, ymax = y1 + height);
}

if (FALSE)
{
  x11()
  plot(1:10, type = "n")
  text(2,2, "test 1", adj = c(0,0))
  box = textBox(2,2,"test 1", adj = c(0,0), cex = 1)
  rect(box[1], box[2], box[3], box[4])

  text(4,1, "Test WW 718", adj = c(0.5, 0.5));
  box = textBox(4,1, "Test WW 718", adj = c(0.5, 0.5));
  rect(box[1], box[2], box[3], box[4])

  text(1, 7,"gyj", adj = c(0.1, 0.8));
  box = textBox(1, 7,"gyj", adj = c(0.1, 0.8));
  rect(box[1], box[2], box[3], box[4])
}




twoSideBarplot = function(nleft, nright,
                          textLeft = nleft,
                          textRight = nright,
                          textLeft.up = NULL,
                          textLeft.dn = NULL,
                          textRight.up = NULL,
                          textRight.dn = NULL,
                          lineSpacing = 1.1,
                          border.left = 1,
                          border.right = 1,
                          col.left = 1,
                          col.right = 1,
                          bg = 0,
                          yLabels = NULL,
                          cex.lab.y = 1,
                          col.lab.y = 1,
                          separatorPositions = NULL,
                          sep.lwd = 1,
                          sep.lty = 1,
                          sep.col = 1,
                          sep.ext = FALSE,

                        lim = NULL,
                        minLim = NULL,
                        barGap = 0,
                        insideThreshold = 0.8,
                        textMargin = 0.5,
                        labelMargin = 0.5,
                        text.minXPosition = 0,
                        main = "",
                        xaxt = "n",
                        xlab = "",
                        cex.main = 1.4,
                        cex.text = 1,
                        col.text = 1,
                        frame = FALSE,
                        symmetric = FALSE,
                          ...)
{
  n = length(nleft);
  if (n!=length(nright)) stop("lengths of nleft and 'nright' must be the same");

  if (any(nleft<0, na.rm = TRUE)) stop("'nleft' must all be non-negative");
  if (any(nright<0, na.rm = TRUE)) stop("'nright' must all be non-negative");

  checkOrExtend = function(x, n, name)
  {
    if (length(x)==1) x = rep(x, n);
    if (length(x)!=n) stop("'", name, "' must have length 1 or number of entries.");
    x;
  }

  col.left = checkOrExtend(col.left, n, "col.left"); 
  col.right = checkOrExtend(col.right, n, "col.right"); 
  bg = checkOrExtend(bg, n, "bg"); 
  border.left = checkOrExtend(border.left, n, "border.left"); 
  border.right = checkOrExtend(border.right, n, "border.right"); 
  cex.lab.y = checkOrExtend(cex.lab.y, n, "cex.lab.y"); 
  col.lab.y = checkOrExtend(col.lab.y, n, "col.lab.y"); 
  cex.text = checkOrExtend(cex.text, n, "cex.text");
  col.text = checkOrExtend(col.text, n, "col.text");

  #nright[is.na(nright)] = 0
  #nleft[is.na(nleft)] = 0
  max = max(nright, na.rm = TRUE);
  min = -max(nleft, na.rm = TRUE);

  if (symmetric) {max = max(max, -min); min = -max; }
  if (!is.null(minLim))
  {
    if (min > -minLim) min = -minLim;
    if (max < minLim) max = minLim;
  }

  if (is.null(lim)) lim = c(min, max);

  plot(lim, c(0,1),  xlim = lim, ylim = c(0,1), xaxt = xaxt, yaxt = "n", main = main, cex.main = cex.main,
       type = "n", xlab = xlab, ylab = "", frame = frame);

  textMargin.usr = textMargin * strwidth("M", cex = cex.text);
  labelMargin.usr = labelMargin * strwidth("M", cex = cex.text);

  box = par("usr");
  xMin = box[1];
  xMax = box[2];
  yMin = box[3];
  yMax = box[4];

  y0.all = c((n-1):0)/n * (yMax - yMin)+ yMin;
  y1.all = c(n:1)/n * (yMax - yMin) + yMin;

  y1 = y1.all - (y1.all-y0.all) * barGap/2
  y0 = y0.all + (y1.all-y0.all) * barGap/2

  nSeparators = length(separatorPositions);
  if (nSeparators > 0)
  {
    sep.col = checkOrExtend(sep.col, nSeparators, "sep.col");
    sep.lty = checkOrExtend(sep.lty, nSeparators, "sep.lty");
    sep.lwd = checkOrExtend(sep.lwd, nSeparators, "sep.lwd");
    sep.ext = checkOrExtend(sep.ext, nSeparators, "sep.lwd");
  }
  

  extension = par("mai")[2] * # left margin width in inches
              par("cxy")[1] / par("cin")[1]   # charcter size in user coordinates/character size in inches

  out = list(ytop = y1, ybottom = y0,
             box = par("usr"),
             yMid = (y0 + y1)/2,
             leftMargin = extension);

  # Background
  rect(xMin-extension, y1.all, xMax, y0.all, col = bg, border = bg, xpd = TRUE);

  # Separators
  if (nSeparators > 0)
  {
    for (al in 1:nSeparators)
    {
      lim1 = c(xMin, xMax);
      if (sep.ext[al]) lim1[1] = lim1[1] - extension;
      ysep1 = if (separatorPositions[al]==0) y1.all[1] else y0.all[separatorPositions[al]];
      lines(lim1, rep(ysep1, 2),
            col = sep.col[al],
            lty = sep.lty[al],
            lwd = sep.lwd[al], xpd = sep.ext[al]);
    }
  }

  
  # Bars
  rect(-nleft, y1, rep(0, n), y0, col = col.left, border = border.left);
  rect(nright, y1, rep(0, n), y0, col = col.right, border = border.right);

  # Text left and right
  nleft.text = replaceMissing(nleft);
  nleft.text[ nleft.text > -lim[1]] = -lim[1];
  nleft.text[ nleft.text < -text.minXPosition] = -text.minXPosition;
  nright.text = nright;
  nright.text[ nright > lim[2]] = lim[2];
  nright.text[ nright.text < text.minXPosition] = text.minXPosition;
  ytext = (y1  + y0)/2
  width = lim[2] - lim[1];
                   
  if (!is.null(textLeft))
  {
      inside = nleft.text/(-lim[1]) > insideThreshold
      adj.x = ifelse(inside, 0, 1);
      shiftDir = 2*(inside-0.5)
      textLeft.x = -nleft.text +  textMargin.usr *shiftDir;
      out = c(out, list(textLeft.x = nleft.text +  textMargin.usr *shiftDir, 
                        textLeft.adj.x = adj.x,
                        textLeft.cex = cex.text,
                        textLeft.box = do.call(cbind, mymapply(textBox,
                            x = textLeft.x, y = ytext,
                            text = textLeft, cex = cex.text, 
                            adj = lapply(adj.x, c, 0.5)))));
      if (!is.null(textLeft.up) | !is.null(textLeft.dn))
      {
        warning("Output text boxes are inaccurate for split text.");
        if (!is.null(textLeft.up)) 
          for (i in 1:n)
            text(textLeft.x[i], ytext[[i]] + lineSpacing/2 * strheight("W"),
               textLeft.up[[i]], col = col.text[i], cex = cex.text[i], adj = c(adj.x[i], 0.5), xpd = TRUE);

        if (!is.null(textLeft.dn)) 
          for (i in 1:n)
            text(textLeft.x[i], ytext[[i]] - lineSpacing/2 * strheight("W"),
               textLeft.dn[[i]], col = col.text[i], cex = cex.text[i], adj = c(adj.x[i], 0.5), xpd = TRUE);
      } else 
      for (i in 1:n)
        text(textLeft.x[i], ytext[[i]], 
           textLeft[[i]], col = col.text[i], cex = cex.text[i], adj = c(adj.x[i], 0.5), xpd = TRUE);
  }
  if (!is.null(textRight))
  {
    inside = nright.text/lim[2] > insideThreshold
    adj.x = ifelse(inside, 1, 0);
    shiftDir = -2*(inside-0.5)
    textRight.x = nright.text +  textMargin.usr *shiftDir;
    out = c(out, list(textRight.x = nright.text +  textMargin.usr *shiftDir, 
                        textRight.adj.x = adj.x,
                        textRight.cex = cex.text,
                        textRight.box = do.call(cbind, mymapply(textBox,
                            x = textRight.x, y = ytext,
                            text = textRight, cex = cex.text, 
                            adj = lapply(adj.x, c, 0.5)))));
    if (!is.null(textRight.up) | !is.null(textRight.dn))
    {
      if (!is.null(textRight.up))
        for (i in 1:n)
          text(textRight.x[i], ytext[[i]] + lineSpacing/2 * strheight("W"), 
             textRight.up[[i]], col = col.text[i], cex = cex.text[i], adj = c(adj.x[i], 0.5), xpd = TRUE);
      if (!is.null(textRight.dn))
        for (i in 1:n)
          text(textRight.x[i], ytext[[i]] - lineSpacing/2 * strheight("W"), 
             textRight.dn[[i]], col = col.text[i], cex = cex.text[i], adj = c(adj.x[i], 0.5), xpd = TRUE);
    } else
      for (i in 1:n)
        text(textRight.x[i], ytext[[i]], 
           textRight[[i]], col = col.text[i], cex = cex.text[i], adj = c(adj.x[i], 0.5), xpd = TRUE);
  }
  # Labels
  if (!is.null(yLabels))
  {
    mapply(text, x = rep(xMin-labelMargin.usr, n), y = ytext, labels = yLabels, cex = cex.lab.y, col = col.lab.y,
           MoreArgs = list(xpd = TRUE, adj = c(1, 0.5)));
  }
  invisible(out);
}


#==================================================================================================
#
# labeledBarplot3
#
#==================================================================================================

checkOrExtend = function(x, n, name)
{
  if (missing(name)) name = as.character(match.call(expand.dots = FALSE)$x);q
  if (length(x)==0) return(x);
  if (length(x)==1) x = rep(x, n);
  if (length(x)!=n) stop("'", name, "' must have length 1 or number of entries.");
  x;
}



labeledBarplot3 = function(x,
                          errors = NULL,
                          twoSidedErrors = TRUE,
                          text = NULL,
                          border = 1,
                          col = 1,
                          bg = 0,
                          yLabels = NULL,
                          cex.lab.y = 1,
                          col.lab.y = 1,
                          separatorPositions = NULL,
                          separatorSpace = 0,  ### In user units, which run from 0 to 1.
                          sep.lwd = 1,
                          sep.lty = 1,
                          sep.col = 1,
                          sep.ext = 0,
                          labelGap = 0.5,

                        lim = NULL,
                        barGap = 0,
                        insideThreshold = 0.8,
                        textMargin = 0.5,
                        margin = 0.1, 
                        textMat.minXPosition = 0,
                        main = "",
                        xlab = "",
                        cex.main = 1.4,
                        cex.text = 1,
                        col.text = 1,
                        symmetric = FALSE,
                        yaxt = "n",
                        frame = TRUE,

                        formatYLabels = FALSE,
                        yLabels.maxFracExtWidth = 0.90,
                        yLabels.maxLines = Inf,
                        ...)
{
  n = length(x);

  col = checkOrExtend(col, n, "col"); 
  bg = checkOrExtend(bg, n, "bg"); 
  border = checkOrExtend(border, n, "border"); 
  cex.lab.y = checkOrExtend(cex.lab.y, n, "cex.lab.y"); 
  col.lab.y = checkOrExtend(col.lab.y, n, "col.lab.y"); 
  cex.text = checkOrExtend(cex.text, n, "cex.text");
  col.text = checkOrExtend(col.text, n, "col.text");

  max = max(x, na.rm = TRUE);
  min = min(x, na.rm = TRUE);

  if (max < 0) max = 0;
  if (min > 0) min = 0;

  if (symmetric) {max = max(max, -min); min = -max; }

  if (is.null(lim)) lim = c(min, max);

  plot(lim, c(0,1),  xlim = lim, ylim = c(0,1), , main = main, cex.main = cex.main,
       type = "n", xlab = xlab, ylab = "", yaxt = yaxt, xaxs = "i", frame = frame,...);

  textMargin.usr = textMargin * strwidth("M", cex = cex.text);

  box = par("usr");
  xMin = box[1];
  xMax = box[2];
  yMin = box[3];
  yMax = box[4];

  nSeparators = length(separatorPositions);
  barSpace = (yMax - yMin) - nSeparators * separatorSpace;
  barStep = barSpace/n;
  nSepBelowBar = sapply(1:n, function(i) sum(separatorPositions >= i));

  y0.all = c((n-1):0) * barStep+ yMin + separatorSpace * nSepBelowBar;
  y1.all = c(n:1) * barStep + yMin + separatorSpace * nSepBelowBar;
 
  y1 = y1.all - barStep * barGap/2
  y0 = y0.all + barStep * barGap/2

  if (nSeparators > 0)
  {
    sep.col = checkOrExtend(sep.col, nSeparators, "sep.col");
    sep.lty = checkOrExtend(sep.lty, nSeparators, "sep.lty");
    sep.lwd = checkOrExtend(sep.lwd, nSeparators, "sep.lwd");
    sep.ext = checkOrExtend(sep.ext, nSeparators, "sep.ext");
  }
  

  extension = par("mai")[2] * # left margin width in inches
              par("cxy")[1] / par("cin")[1]   # charcter size in user corrdinates/character size in inches

  y.separators = c(y1.all[1]+separatorSpace, y0.all)[separatorPositions+1];

  if (!isTRUE(all.equal(y.separators, sort(y.separators, decreasing = TRUE))))
    browser("incorrect order.")

  out = list(ytop = y0, ybottom = y1,
             box = par("usr"),
             yMid = (y0 + y1)/2,
             leftMargin = extension,
             yTop.separators = y.separators,
             yMid.separators = y.separators - separatorSpace/2);

  # Background
  #cond = try(any(bg!=0));
  #if (is.na(cond) || inherits(cond, "try-error")) browser()
  if (any(bg!=0)) rect(xMin-extension, y1.all, xMax, y0.all, col = bg, border = bg, xpd = TRUE);

  # Separators
  if (nSeparators > 0)
  {
    for (al in 1:nSeparators)
    {
      lim1 = c(xMin, xMax);
      if (sep.ext[al]>0) lim1[1] = lim1[1] - extension*sep.ext[al];
      lines(lim1, rep(y.separators[al], 2),
            col = sep.col[al],
            lty = sep.lty[al],
            lwd = sep.lwd[al], xpd = TRUE);
    }
  }


  # Re-draw the plot frame
  if (frame) rect(xMin, yMin, xMax, yMax) else lines(c(0,0), c(yMin, yMax));
  # Bars
  rect(x, y1, rep(0, n), y0, col = col, border = border);

  # Text 
  ytext = (y1  + y0)/2
  width = lim[2] - lim[1];
  if (!is.null(text))
  {
      sign = sign(x);
      inside.right = ( x > 0 & (x-lim[1])/width > insideThreshold) 
      inside.left = ( x < 0 & (lim[2]-x)/width > insideThreshold) 
      inside = inside.left | inside.right;
      #inside[x<0] = !inside[x<0];
      adj.x = ifelse(inside, 1, 0)
      adj.x[x<0] = 1-adj.x[x<0];
      shiftDir = -2*(inside-0.5)*sign;
      for (i in 1:n)
        text(x[i] +  textMargin.usr *shiftDir[i], ytext[i], 
           text[i], col = col.text[i], cex = cex.text[i], adj = c(adj.x[i], 0.5));
  }
  # Labels
  if (!is.null(yLabels))
  {
    if (formatYLabels)
    {
      yLabels = formatLabels(yLabels, maxWidth = yLabels.maxFracExtWidth * extension,
                             maxLines = yLabels.maxLines);
    }
    labelGap.usr = labelGap * strwidth("M");
    text(rep(xMin, n)-labelGap.usr, ytext, yLabels, cex = cex.lab.y, xpd = TRUE, col = col.lab.y,
         adj = c(1, 0.5))
  }
  invisible(out);
}



labeledBarplot3.multiPage = function(x,
    text = NULL,
    border = 1,
    col = 1,
    bg = 0,
    yLabels = NULL,
    cex.lab.y = 1,
    col.lab.y = 1,
    cex.text = 1,
    col.text = 1,
    separatorPositions = NULL,
    sep.lwd = 1,
    sep.lty = 1,
    sep.col = 1,
    sep.ext = 0,

    lim = NULL,
    main = "",
    symmetric = FALSE,

     # Paging options
    rowsPerPage = NULL, maxRowsPerPage = 50,
    addPageNumberToMain = TRUE,
    pageNumberStart = " (page ",
    pageNumberEnd = ")",

    # Optional, user-supplied plot embellishment
    afterFnc = NULL,
    afterArgs = list(),

    ...)
{
  n = length(x);

  if (is.null(rowsPerPage))
  {
    nPages.rows = ceiling(nr/maxRowsPerPage);
    rowsPerPage = allocateJobs(nr, nPages.rows);
  } else
    nPages.rows = length(rowsPerPage);

  max = max(x, na.rm = TRUE);
  min = min(x, na.rm = TRUE);

  if (max < 0) max = 0;
  if (min > 0) min = 0;

  if (symmetric) {max = max(max, -min); min = -max; }

  if (is.null(lim)) lim = c(min, max);

  col = checkOrExtend(col, n, "col");
  bg = checkOrExtend(bg, n, "bg");
  border = checkOrExtend(border, n, "border");
  cex.lab.y = checkOrExtend(cex.lab.y, n, "cex.lab.y");
  col.lab.y = checkOrExtend(col.lab.y, n, "col.lab.y");
  cex.text = checkOrExtend(cex.text, n, "cex.text");
  col.text = checkOrExtend(col.text, n, "col.text");

  nSeparators = length(separatorPositions);
  if (nSeparators > 0)
  {
    sep.col = checkOrExtend(sep.col, nSeparators, "sep.col");
    sep.lty = checkOrExtend(sep.lty, nSeparators, "sep.lty");
    sep.lwd = checkOrExtend(sep.lwd, nSeparators, "sep.lwd");
    sep.ext = checkOrExtend(sep.ext, nSeparators, "sep.ext");
  }

  page = 1;

  multiPage = nPages.rows > 1

  for (page.row in 1:nPages.rows)
  {
    rows = rowsPerPage[[page.row]];
    if (nSeparators > 0)
    {
      keep.hs = (separatorPositions+1) %in% rows;
    } else
      keep.hs = numeric(0)

    main.1 = main;
    if (addPageNumberToMain & multiPage) main.1 = spaste(main, pageNumberStart, page, pageNumberEnd);

    plotCoords = labeledBarplot3(x[rows],
	  text = text[rows],
	  border = border[rows],
	  col = col[rows],
	  bg = bg[rows],
	  yLabels = yLabels[rows],
	  cex.lab.y = cex.lab.y[rows],
	  col.lab.y = col.lab.y[rows],
	  separatorPositions = separatorPositions[keep.hs]-min(rows) +1,
	  sep.lwd = sep.lwd[keep.hs],
	  sep.lty = sep.lty[keep.hs],
	  sep.col = sep.col[keep.hs],
	  sep.ext = sep.ext[keep.hs],

	  lim = lim,
	  main = main.1,
	  symmetric = symmetric,
	  ...)

     if (!is.null(afterFnc))
     {
       afterFnc = match.fun(afterFnc);
       do.call(afterFnc, c(list(
          plotCoords = plotCoords,
          rows = rows,
          separatorIndex = if (length(keep.hs)>0) which(keep.hs) else numeric(0),
          x = x[rows],
          text = text[rows],
          border = border[rows],
          col = col[rows],
          bg = bg[rows],
          yLabels = yLabels[rows],
          cex.lab.y = cex.lab.y[rows],
          col.lab.y = col.lab.y[rows],
          cex.text = cex.text[rows],
          col.text = col.text[rows],
          separatorPositions = separatorPositions[keep.hs]-min(rows) + 1,
          sep.lwd = sep.lwd[keep.hs],
          sep.lty = sep.lty[keep.hs],
          sep.col = sep.col[keep.hs],
          sep.ext = sep.ext[keep.hs]),
          afterArgs));
     }
     page = page + 1;

  }
  invisible(NULL);
}


if (FALSE)
{
  n = 30;
  x = 1:n;

  separators = c(0, which(x %% 5==0))
  sizeGrWindow(14, 6);
  par(mfrow = c(1,4));
  labeledBarplot3.multiPage(x,
            text = x,
            yLabels = spaste("r ", x),
            rowsPerPage = allocateJobs(n, 4),
            main = "test", separatorPositions = separators,
            border = 1, col = 1:n, barGap = 0.5, separatorSpace = 0.1);

}


    







#==================================================================================================
#
# Generic plot annotation function
#
#==================================================================================================

annotatePlot = function(
  at.x, at.y,

  plotbox = par("usr"),

  xLeft = NULL, xRight = NULL,
  yTop = NULL, yBottom = NULL,

  xLabels, yLabels = NULL, 
  xSymbols = NULL, ySymbols = NULL, 
  colorLabels = NULL, 
  xColorLabels = FALSE, yColorLabels = FALSE,
  checkColorsValid = TRUE,
  invertColors = FALSE, 
  xLabelsPosition = "bottom",
  xLabelsAngle = 45,
  xLabelsAdj = 1,
  xColorWidth = 0.05,
  yColorWidth = 0.05,
  cex.lab = NULL, 
  cex.lab.x = cex.lab,
  cex.lab.y = cex.lab,
  colors.lab.x = 1,
  colors.lab.y = 1,
  bg.lab.x = NULL,
  bg.lab.y = NULL,

  # Separator line specification                   
  verticalSeparator.x = NULL,
  verticalSeparator.col = 1,  
  verticalSeparator.lty = 1,
  verticalSeparator.lwd = 1,
  verticalSeparator.ext = 0,

  horizontalSeparator.y = NULL,
  horizontalSeparator.col = 1,  
  horizontalSeparator.lty = 1,
  horizontalSeparator.lwd = 1,
  horizontalSeparator.ext = 0,
  ... ) 
{
  if (!is.null(colorLabels)) {xColorLabels = colorLabels; yColorLabels = colorLabels; }

  nCols = length(at.x);
  nRows = length(at.y);
 
  if (!is.null(xLabels) && nCols!=length(xLabels))
    stop("If given, 'xLabels' must have the same length as 'at.x'."); 

  if (checkColorsValid)
  {
    xValidColors = !is.na(match(substring(xLabels, 3), colors()));
    yValidColors = !is.na(match(substring(yLabels, 3), colors()));
  } else {
    xValidColors = rep(TRUE, length(xLabels));
    yValidColors = rep(TRUE, length(yLabels));
  }

  if (sum(xValidColors)>0) xColorLabInd = c(1:length(xLabels))[xValidColors]
  if (sum(!xValidColors)>0) xTextLabInd = c(1:length(xLabels))[!xValidColors]

  if (sum(yValidColors)>0) yColorLabInd = c(1:length(yLabels))[yValidColors]
  if (sum(!yValidColors)>0) yTextLabInd = c(1:length(yLabels))[!yValidColors]

  xLabPos = charmatch(xLabelsPosition, c("bottom", "top"));
  if (is.na(xLabPos))
    stop("Argument 'xLabelsPosition' must be (a unique abbreviation of) 'bottom', 'top'");

  nxlabels = length(xLabels)
  xmin = plotbox[1]; xmax = plotbox[2]; ymin = plotbox[3]; yrange = plotbox[4]-ymin;
  ymax = plotbox[4]; xrange = xmax - xmin;
  if (is.null(xLeft)) xLeft = at.x - xrange/nCols/2;
  if (is.null(xRight)) xRight = at.x + xrange/nCols/2;
  if (is.null(yTop)) yTop = at.y + yrange/nRows/2;
  if (is.null(yBottom)) yBottom = at.y - yrange/nRows/2;

  xspacing = at.x[2] - at.x[1];
  yspacing = abs(at.y[2] - at.y[1]);

  nylabels = length(yLabels)
  offsetx = par("cxy")[1] / 3;
  offsety = par("cxy")[2] / 3;
  # Transform fractional widths into coordinate widths
  xColW = min(xmax - xmin, ymax - ymin) * xColorWidth;
  yColW = min(xmax - xmin, ymax - ymin) * yColorWidth;

  if (any(xValidColors)) offsety = offsety + xColW;
  if (any(yValidColors)) offsetx = offsetx + yColW;

  # Create the background for column and row labels.

  extension.left = par("mai")[2] * # left margin width in inches
                   par("cxy")[1] / par("cin")[1]   # charcter size in user corrdinates/character size in inches

  extension.bottom = par("mai")[1] * 
                   par("cxy")[2] / par("cin")[2]- # charcter size in user corrdinates/character size in inches
                      offsety   
                     
  extension.top = par("mai")[3] * 
                   par("cxy")[2] / par("cin")[2]-   # charcter size in user corrdinates/character size in inches
                     offsety

  figureBox = par("usr");
  figXrange = figureBox[2] - figureBox[1];
  figYrange = figureBox[4] - figureBox[3];
  if (!is.null(bg.lab.x))
  {
    bg.lab.x = .extend(bg.lab.x, nCols);
    if (xLabPos==1)
    {
      y0 = ymin;
      ext = extension.bottom;
      sign = 1;
    } else {
      y0 = ymax;
      ext = extension.top;
      sign = -1;
    }
    figureDims = par("pin");
    angle = xLabelsAngle/180*pi;
    ratio = figureDims[1]/figureDims[2] * figYrange/figXrange;
    ext.x = -sign * extension.bottom * 1/tan(angle)/ratio;
    ext.y = sign * extension.bottom * sign(sin(angle))
    for (c in 1:nCols)
       polygon(x = c(xLeft[c], xLeft[c], xLeft[c] + ext.x, xRight[c] + ext.x, xRight[c], xRight[c]),
               y = c(y0, y0-sign*offsety, y0-sign*offsety - ext.y, y0-sign*offsety - ext.y, 
                     y0-sign*offsety, y0), 
               border = bg.lab.x[c], col = bg.lab.x[c], xpd = TRUE);
  }

  if (!is.null(bg.lab.y))
  {
    bg.lab.y = .extend(bg.lab.y, nCols);
    reverseRows = TRUE;
    if (reverseRows)
    {
      bg.lab.y = rev(bg.lab.y);
    }
    for (r in 1:nRows)
      rect(xmin-extension.left, yBottom[r], xmin, yTop[r],
           col = bg.lab.y[r], border = bg.lab.y[r], xpd = TRUE);
  }

  # Write out labels
  if (sum(!xValidColors)>0)
  {
    xLabYPos = ifelse(xLabPos==1, ymin - offsety, ymax + offsety)
    if (is.null(cex.lab)) cex.lab = 1;
    text(at.x[xTextLabInd] , xLabYPos, srt = xLabelsAngle, 
          adj = xLabelsAdj, labels = xLabels[xTextLabInd], xpd = TRUE, cex = cex.lab.x, col = colors.lab.x)
  }
  if (sum(xValidColors)>0)
  {
    baseY = ifelse(xLabPos==1, ymin-offsety, ymax + offsety);
    deltaY = ifelse(xLabPos==1, xColW, -xColW);
    rect(xleft = at.x[xColorLabInd] - xspacing/2, ybottom = baseY,
         xright = at.x[xColorLabInd] + xspacing/2, ytop = baseY + deltaY,
         density = -1,  col = substring(xLabels[xColorLabInd], 3), 
         border = substring(xLabels[xColorLabInd], 3), xpd = TRUE)
    if (!is.null(xSymbols))
      text ( at.x[xColorLabInd], baseY - sign(deltaY)* offsety, xSymbols[xColorLabInd], 
             adj = xLabelsAdj, 
             xpd = TRUE, srt = xLabelsAngle, cex = cex.lab.x, col = colors.lab.x);
  }
  if (sum(!yValidColors)>0)
  {
    if (is.null(cex.lab)) cex.lab = 1;
    text(xmin - offsetx, at.y[yTextLabInd], srt = 0, 
         adj = c(1, 0.5), labels = yLabels[yTextLabInd], xpd = TRUE, cex = cex.lab.y, col = colors.lab.y )
  } 
  if (sum(yValidColors)>0)
  {
    rect(xleft = xmin- offsetx, ybottom = at.y[yColorLabInd] - yspacing/2,
         xright = xmin- offsetx+yColW, ytop = at.y[yColorLabInd] + yspacing/2, 
         density = -1,  col = substring(yLabels[yColorLabInd], 3), 
         border = substring(yLabels[yColorLabInd], 3), xpd = TRUE)
    if (!is.null(ySymbols))
      text (xmin+ yColW - 2*offsetx, 
            at.y[yColorLabInd], ySymbols[yColorLabInd], 
            adj = c(1, 0.5), xpd = TRUE, cex = cex.lab.y, col = colors.lab.y);
  }

  # Draw separator lines, if requested

  if (!is.null(verticalSeparator.x))
  {
    nLines = length(verticalSeparator.x);
    vs.col = .extend(verticalSeparator.col, nLines);
    vs.lty = .extend(verticalSeparator.lty, nLines);
    vs.lwd = .extend(verticalSeparator.lwd, nLines);
    vs.ext = .extend(verticalSeparator.ext, nLines);
    if (any(verticalSeparator.x < 0 | verticalSeparator.x > nCols))
      stop("If given. 'verticalSeparator.x' must all be between 0 and the number of columns.");
    x.lines = ifelse(verticalSeparator.x>0, xRight[verticalSeparator.x], xLeft[1]);
    for (l in 1:nLines)
      lines(rep(x.lines[l], 2), c(ymin, ymax), col = vs.col[l], lty = vs.lty[l], lwd = vs.lwd[l]);

    angle = xLabelsAngle/180*pi;
    if (xLabelsPosition =="bottom") 
    {
      sign = 1;
      y0 = ymin;
    } else {
      sign = -1;
      y0 = ymax;
    }
    figureDims = par("pin");
    ratio = figureDims[1]/figureDims[2] * figYrange/figXrange;
    ext.x = -sign * extension.bottom * 1/tan(angle)/ratio;
    ext.y = sign * extension.bottom * sign(sin(angle))
    for (l in 1:nLines)
         lines(c(x.lines[l], x.lines[l], x.lines[l] + vs.ext * ext.x), 
               c(y0, y0-sign*offsety, y0-sign*offsety - vs.ext * ext.y),  
                 col = vs.col[l], lty = vs.lty[l], lwd = vs.lwd[l], xpd = TRUE);
  }

  if (!is.null(horizontalSeparator.y))
  {
    if (any(horizontalSeparator.y < 0 | horizontalSeparator.y > nRows))
      stop("If given. 'horizontalSeparator.y' must all be between 0 and the number of rows.");
    reverseRows = TRUE;
    if (reverseRows) 
    {
      horizontalSeparator.y = nRows - horizontalSeparator.y+1;
      y.lines = ifelse( horizontalSeparator.y <=nRows, yBottom[horizontalSeparator.y], yTop[nRows]);
    } else {
      y.lines = ifelse( horizontalSeparator.y > 0, yBottom[horizontalSeparator.y], yTop[1]);
    }
    nLines = length(horizontalSeparator.y);
    vs.col = .extend(horizontalSeparator.col, nLines);
    vs.lty = .extend(horizontalSeparator.lty, nLines);
    vs.lwd = .extend(horizontalSeparator.lwd, nLines);
    vs.ext = .extend(horizontalSeparator.ext, nLines);
    for (l in 1:nLines)
      lines(c(xmin-vs.ext[l]*extension.left, xmax), rep(y.lines[l], 2), 
            col = vs.col[l], lty = vs.lty[l], lwd = vs.lwd[l], xpd = TRUE);
  }

}


#=====================================================================================================
#
# colorScale
#
#=====================================================================================================

colorScale = function(n, lowEnd, highEnd, gamma = 1 )
{
  if (length(lowEnd)==1) lowEnd = as.numeric(col2rgb(lowEnd));
  if (length(highEnd)==1) highEnd = as.numeric(col2rgb(highEnd));
  if (length(lowEnd)!=3 || !is.numeric(lowEnd) || any(!is.finite(lowEnd))) 
     stop("'lowEnd' must be a color specified either as a single value or an RGB triplet.")

  if (length(highEnd)!=3 || !is.numeric(highEnd) || any(!is.finite(highEnd))) 
     stop("'highEnd' must be a color specified either as a single value or an RGB triplet.")

  if (any(lowEnd> 1)) lowEnd = lowEnd/255;
  if (any(highEnd> 1)) highEnd = highEnd/255;

  middle = c(1, 1, 1)
  half = as.integer(n/2)
  if (n%%2 == 0) {
      index1 = c(1:half)
      index2 = c(1:half) + half
      frac1 = ((index1 - 1)/(half - 1))^(1/gamma)
      frac2 = rev(frac1)
  }
  else {
      index1 = c(1:(half + 1))
      index2 = c(1:half) + half + 1
      frac1 = (c(0:half)/half)^(1/gamma)
      frac2 = rev((c(1:half)/half)^(1/gamma))
  }
  cols = matrix(0, n, 3)
  for (c in 1:3) {
      cols[index1, c] = lowEnd[c] + (middle[c] - lowEnd[c]) * frac1
      cols[index2, c] = highEnd[c] + (middle[c] - highEnd[c]) * frac2
  }
  rgb(cols[, 1], cols[, 2], cols[, 3], maxColorValue = 1)
}

#=====================================================================================================
#
# plotPreservationLegend
#
#=====================================================================================================

plotPreservationLegend = function(palette, labels, cex.labels = 1, ...)
{
  n = length(palette);
  plot(c(0, 1), c(1,1), type = "n", xlim = c(-1, 1), ylim = c(0,1), axes = FALSE, frame = TRUE,
       xlab = "", ylab = "",...);
  xmin = par("usr")[1];
  xmax = par("usr")[2];
  ymin = par("usr")[3];
  ymax = par("usr")[4]; 
  xl = seq(from = xmin, to = xmax, length.out = n+1)[1:n];
  xr = seq(from = xmin, to = xmax, length.out = n+1)[2:(n+1)];
  rect(xl, ymin, xr, ymax, col = palette, border = palette);
  nLabels = length(labels); 
  xText = seq(from=xmin, to = xmax, length.out = nLabels);
  for (l in 1:nLabels)
    text(xText[l], -0.40, labels[l], adj = c((xText[l]-xmin)/(xmax-xmin), 1), xpd = TRUE,
         cex = cex.labels);
}

plotPreservationLegend2 = function(palette, labels, cex.labels = 1, xlim = c(-1, 1), ...)
{
  n = length(palette);
  plot(c(0, 1), c(1,1), type = "n", xlim = xlim, ylim = c(0,1), yaxt = "n", frame = TRUE,
       xlab = "", ylab = "",...);
  xmin = par("usr")[1];
  xmax = par("usr")[2];
  ymin = par("usr")[3];
  ymax = par("usr")[4]; 
  xl = seq(from = xmin, to = xmax, length.out = n+1)[1:n];
  xr = seq(from = xmin, to = xmax, length.out = n+1)[2:(n+1)];
  rect(xl, ymin, xr, ymax, col = palette, border = palette);
  nLabels = length(labels); 
  xText = seq(from=xmin, to = xmax, length.out = nLabels);
  for (l in 1:nLabels)
    text(xText[l], -2*strheight("W"), labels[l], adj = c((xText[l]-xmin)/(xmax-xmin), 1), xpd = TRUE,
         cex = cex.labels);
}



#=====================================================================================================
#
# whittle down connections in a dense module for plotting purposes.
#
#=====================================================================================================

segmentDist = function(network, labels)
{
  nSegments = length(unique(labels));
  segments = lapply(1:nSegments, function(x) which(labels==x));
  sizes = lapply(1:nSegments, function(x) sum(labels==x));

  segSim = matrix(0, nSegments, nSegments);
  nearestSeg = nearestSeg2 = matrix(0, nSegments, nSegments);
  for (s1 in 1:(nSegments-1)) for (s2 in (s1+1):nSegments)
  {
    net12 = network[ segments[[s1]], segments[[s2]] ];
    segSim[s1, s2] = segSim[s2, s1] = max(net12, na.rm = TRUE);
    which = which.max(net12, na.rm = TRUE)
    which1 = ((which-1) %% sizes[s1]) + 1;
    which2 = floor( (which-1) / sizes[s1]);
    nearestSeg[s1, s2] = which1;
    nearestSeg[s2, s1] = which2;
  }
}

nearestSegments =function(net, labels)
{
  nSegments = length(unique(labels));
  segs = lapply(1:nSegments, function(x) which(labels==x));
  sizes = lapply(1:nSegments, function(x) sum(labels==x));
  max = -1;
  for (s1 in 1:(nSegments-1)) for (s2 in (s1+1):nSegments)
  {
    net12 = net[ segs[[s1]], segs[[s2]] ];
    sim = max(net12, na.rm = TRUE);
    if (sim > max)
    {
      max = sim;
      segment1 = s1;
      segment2 = s2;
    }
  }
  c(segment1, segment2);
}

nearestNodes = function(network, labels, s1, s2)
{
    segment1 = which(labels==s1);
    segment2 = which(labels==s2);
    which = which.max(network[segment1, segment2])
    n1 = length(segment1);
    which1 = ((which-1) %% n1 ) + 1;
    which2 = floor( (which-1) / n1)+1;
    c(segment1[which1], segment2[which2]);
}

whittleConnections = function(network, minConnections = 1, maxConnections = 10, threshold = 0,
                              connectDisconnectedSegments = TRUE)
{
  diag(network) = -1;
  out = network;
  # Drop below-threshold connections
  out[network<threshold] = 0;

  ranks = apply(-network, 2, rank, na.last = TRUE, ties.method = "first");
  ranks2 = apply(-out, 2, rank, na.last = TRUE, ties.method = "first");
  # Drop excessive connections
  drop = ranks2 > maxConnections;
  out[drop ] = 0;
  out[t(drop)] = 0;
  # Restore at least one connection for each node
  restore = ranks <= minConnections;
  out[restore ] = network[restore];
  out[t(restore) ] = network[t(restore)];

  # Determine connected segments and connect them if needed

  repeat {
    dst = as.dist(1-out);
    tree = hclust(dst, method = "single");
    height = (max(dst[dst!=1])+1)/2
    segments = cutreeStatic(tree, cutHeight = height, minSize = 1);
    nSegments = length(unique(segments));
    if (nSegments==1) break;
    # Connect the two nearest segments.
    ns = nearestSegments(network, segments);
    nn = nearestNodes(network, segments, ns[1], ns[2]);
    out[ nn[1], nn[2] ] = out[nn[2], nn[1]] = network[ nn[1], nn[2] ];
    if ( network[ nn[1], nn[2] ] == 0) stop("Two disconnected segments in original networ... sorry!");
  }  

  out;
}

#============================================================================================================#
#
# addAlternatingNewlines
#
#=============================================================================================================

addAlternatingNewlines = function(text, cycle = 3)
{
  n = length(text);
  out = character(n)
  for (i in 1:n)
  {
    nBefore = (i-1)%%cycle;
    nAfter = cycle-nBefore-1;
    out[i] = paste(paste(rep("\n", nBefore), collapse = ""),
                   text[i],
                   paste(rep("\n", nAfter), collapse = ""), sep = "", collapse = "");
  }
  out;
}



#============================================================================================================
#
# Functions for relating miRNA data to gene data
#
#============================================================================================================
sdx = function(x)
{
  x = x[is.finite(x)];
  if (length(x) < 2) return(0);
  sd(x);
}

splitMIRNAs = function(name)
{
  prefix = substring(name, 1, 4);
  numbers = substring(name, 5);

  split = strsplit(numbers, split = "/", fixed = TRUE);
  spaste(prefix, split[[1]]);
}

# This currently assumes that the target organism is mouse.

loadMirnaTargets = function(pathToMirnaData)
{
  targets = list();
  # Targetscan from a Bioconductor library
  if (require(targetscan.Mm.eg.db))
  {
    miRNAmasters = processAnnotationObject(targetscan.Mm.egTARGETS, uniqueOnly = FALSE);
    allMIRNAs = unique(unlist(miRNAmasters))
    table(substring(allMIRNAs, 1, 4))
    mirnaTargets.targetscan = as.data.frame(do.call(rbind, mapply(function(vec, name)
                 {
                    vec2 = unique(unlist(lapply(vec, splitMIRNAs)));
                    cbind(miRNA = vec2, TargetEntrez = rep(name, length(vec2)))
                 }, miRNAmasters, names(miRNAmasters))));
    mirnaTargets.targetscan$TargetEntrez = as.numeric(as.character(mirnaTargets.targetscan$TargetEntrez));
    targets = c(targets, list(targetScan = mirnaTargets.targetscan));
  }

  # MicroRNA.org (target sites by miRanda, scores by mirSVR

  mirTargetData = read.table(
     gzfile(file.path(pathToMirnaData, "MicroRNA.org/mouse_predictions_S_C_aug2010.txt.bz2")),
     sep = "\t", header = TRUE)

  mirnaTargets.mr.org = mirTargetData[, c(2,3)];
  names(mirnaTargets.mr.org) = c("miRNA", "TargetEntrez");
  mirnaTargets.mr.org$miRNA = gsub("^mmu-", "", mirnaTargets.mr.org$miRNA);

  targets = c(targets, list(MicroRNA.org = mirnaTargets.mr.org));

  # microCosm/miRBase targets

  mcTargetData = read.table(
     bzfile(file.path(pathToMirnaData, "microCosm/V5/microCosm-Mmu-simplified.txt.bz2")),
     sep = "\t", header = TRUE)

  keep = grep("mmu-", mcTargetData$miRNA);

  mirnaTargets.mc.org = mcTargetData[keep, c(1,3)];
  mirnaTargets.mc.org[, 1] = gsub('"', '', mirnaTargets.mc.org[, 1])
  mirnaTargets.mc.org[, 1] = gsub('^mmu-', '', mirnaTargets.mc.org[, 1])

  targets = c(targets, list(microCosm = mirnaTargets.mc.org));


  # mirTarBase

  mtbTargetData = read.csv(
     bzfile(spaste(pathToMirnaData, "/MirTarBase/Version4.5-2013-11-01/mmu_MTI.csv.bz2")),
     sep = ",", header = TRUE)

  mtb2 = mtbTargetData[, c(2,5)];
  mtb2$miRNA = gsub('^mmu-', '', mtb2$miRNA);
  colnames(mtb2)[2] = "TargetEntrez";
  targets = c(targets, list(mirTarBase = mtb2));

  targets;
}

collectionOfMirnaTargets = function(targetInformation)
{
  if (!is.null(dim(targetInformation))) targetInformation = list(mirnaTargets = targetInformation);
  ndfs = length(targetInformation);
  geneSets = list();
  groups = list(newGroup("miRNA targets", "miRNA target genes compiled from various databases", 
                         "Downloaded and formatted by Peter Langfelder"));
  index = 1;
  for (ti in 1:ndfs)
  {
    uniqueMirnas = unique(targetInformation[[ti]]$miRNA);
    for (m in uniqueMirnas)
    {
      geneSets = c(geneSets, list(newGeneSet(
                  geneEntrez = targetInformation[[ti]]$TargetEntrez[targetInformation[[ti]]$miRNA==m ],
                  geneEvidence = "other",
                  geneSource = names(targetInformation)[ti],
                  ID = spaste("mirTar.", prependZeros(index, 6)),
                  name = spaste("Targets of ", m, " (", names(targetInformation)[ti], ") "),
                  description = spaste("Targets of ", m, " as suggested by ", names(targetInformation)[ti]),
                  source = names(targetInformation)[ti],
                  organism = "mouse", 
                  internalClassification = c("miRNA targets", names(targetInformation)[ti]),
                  groups = c("miRNA targets", spaste("miRNA targets from ", names(targetInformation)[ti]))
                  )));
      index = index + 1;
    }
    groups = c(groups, list(newGroup(spaste("miRNA targets from ", names(targetInformation)[ti]),
                                     spaste("miRNA targets from ", names(targetInformation)[ti]),
                                     names(targetInformation)[ti])));
  }
  
  newCollection(geneSets, groups);
}

stripMirnaModifierSuffix = function(mname, restrict = TRUE,
                             strip = c("3p", "5p", "*"))
{
  split = strsplit(mname, split = "-", fixed = TRUE);
  if (restrict)
  {
    out = sapply(split, function(x)
    {
      x = x[ !x%in% strip];
      spaste(x, collapse = "-");
    });
  } else
      out = sapply(split, function(x) if (length(x) > 1) spaste(x[1], "-", x[2]) else x[1]);
  if ("*" %in% strip) out= sub("*", "", out, fixed = TRUE);
  out;
}

splitMultiMirnaEntries = function(tab, mirCol, verbose = 1)
{
  if (!is.numeric(mirCol)) mirCol = match(mirCol, colnames(tab));
  mir0 = unique(tab[, mirCol])

  tabByMir = lapply(mir0, function(m) tab[  replaceMissing(tab[, mirCol]==m), ]);
  suffix = multiSub(c("hsa-", "\\*", "-[35]p$", "Mirlet[0-9]*", "Mir[0-9]*", "-", "[0-9]"), 
                    c("", "", "", "", "", "", ""), mir0);
  doSplit = nchar(suffix) > 1;

  if (verbose) print(cbind(mir0[doSplit], suffix[doSplit]))

  for (i in which(doSplit))
  {
    tab1 = tabByMir[[i]];
    m1 = multiSub(spaste(suffix[i], "$"), c(""), mir0[i]);
    suffixSplit = sapply(1:nchar(suffix[i]), function(char) substring(suffix[i], char, char));
    tab2 = do.call(rbind, lapply(suffixSplit, function(s1)
    {
      tab.loc = tab1;
      tab.loc[, mirCol] = spaste(m1, s1);
      tab.loc;
    }));
    tabByMir[[i]] = tab2;
  }
  do.call(rbind, tabByMir);
}
  
#======================================================================================================
#
# Collections from Enrichr genes
#
#======================================================================================================

agingCollectionFromEnrichr = function(enrichrFile, IDBase, nameBase, descriptionBase,
                                 enrichrClass, outputOrganism = NULL)
{
  if (grepl("\\.bz2$", enrichrFile)) {
    con = bzfile(enrichrFile);
  } else
    con = file(enrichrFile);
  fileLines = readLines(con);
  nSets = length(fileLines);

  split = strsplit(fileLines, split = "\t", fixed = TRUE);
  firstString = sapply(split, `[`, 1);
  firstStrSplit = strsplit(firstString, "_", fixed = TRUE);
  organism = sapply(firstStrSplit, `[`, 1);
  first2 = sapply(firstStrSplit, function(s) spaste(s[-1], collapse = "_"));
  IDs = spaste(IDBase, prependZeros(1:nSets, 5));
  setNames = spaste(nameBase, first2);
  setDescriptions = spaste(descriptionBase, first2);

  symbols = lapply(split, function(s) gsub(",[0-9.]*", "", s[-1]));
  printFlush("Converting symbols to Entrez...");
  humanEntrez = convert2entrez(organism = "human", symbol = symbols);
  printFlush("Converting human entrez to specified organisms...");
  orgEntrez = list();
  pind = initProgInd();
  for (set in 1:nSets)
  {
    orgEntrez[[set]] = mapEntrez(entrez = humanEntrez[[set]], 
                               orgTo = organism[set], orgFrom = "human");
    pind = updateProgInd(set/nSets, pind);
  }

  printFlush();

  geneSets = mymapply(newGeneSet, 
               geneEntrez = orgEntrez,
               ID = IDs,
               name = setNames,
               description = setDescriptions,
               organism = organism,
             MoreArgs = list(
               geneSource = "Enrichr",
               geneEvidence = "other",
               source = spaste("Enrichr, ", enrichrClass),
               internalClassification = c("Enrichr", enrichrClass),
               groups = c("Enrichr", enrichrClass)));

  groups = list(newGroup("Enrichr", "Gene sets provided by Enrichr", "Enrichr"),
                newGroup(enrichrClass, spaste(enrichrClass, " gene sets, provided by Enrichr"), "Enrichr"));

  out = newCollection(geneSets, groups);
  if (!is.null(outputOrganism)) 
  {
     printFlush(spaste("Converting collection to ", outputOrganism))
     out = convertCollectionToOrganism(out, organism = outputOrganism);
  }

  out;
}

EncodeCollectionFromEnrichr = function(enrichrFile, IDBase, nameBase, descriptionBase,
                                 enrichrClass, outputOrganism = NULL)
{
  if (grepl("\\.bz2$", enrichrFile)) {
    con = bzfile(enrichrFile);
  } else
    con = file(enrichrFile);
  fileLines = readLines(con);
  nSets = length(fileLines);

  split = strsplit(fileLines, split = "\t", fixed = TRUE);
  firstString = sapply(split, `[`, 1);
  firstStrSplit = strsplit(firstString, "_", fixed = TRUE);
  TFName = sapply(firstStrSplit, `[`, 1);
  CellLineName = sapply(firstStrSplit, `[`, 2);
  genome = sapply(firstStrSplit, `[`, 3);

  organism = gsub("[0-9]", "", genome);
  organism = multiSub(c("hg"), c("hs"), organism);

  IDs = spaste(IDBase, prependZeros(1:nSets, 5));
  setNames = spaste(nameBase, TFName, " in ", CellLineName);
  setDescriptions = spaste(descriptionBase, first2);

  symbols = lapply(split, function(s) gsub(",[0-9.]*", "", s[-1]));
  printFlush("Converting symbols to Entrez...");
  humanEntrez = convert2entrez(organism = "human", symbol = symbols);
  printFlush("Converting human entrez to specified organisms...");
  orgEntrez = list();
  pind = initProgInd();
  for (set in 1:nSets)
  {
    orgEntrez[[set]] = mapEntrez(entrez = humanEntrez[[set]],
                               orgTo = organism[set], orgFrom = "human");
    pind = updateProgInd(set/nSets, pind);
  }

  printFlush();

  geneSets = mymapply(newGeneSet,
               geneEntrez = orgEntrez,
               ID = IDs,
               name = setNames,
               description = setDescriptions,
               organism = organism,
             MoreArgs = list(
               geneSource = "Enrichr",
               geneEvidence = "other",
               source = spaste("Enrichr, ", enrichrClass),
               internalClassification = c("Enrichr", enrichrClass),
               groups = c("Enrichr", enrichrClass)));

  groups = list(newGroup("Enrichr", "Gene sets provided by Enrichr", "Enrichr"),
                newGroup(enrichrClass, spaste(enrichrClass, " gene sets, provided by Enrichr"), "Enrichr"));

  out = newCollection(geneSets, groups);
  if (!is.null(outputOrganism))
  {
     printFlush(spaste("Converting collection to ", outputOrganism))
     out = convertCollectionToOrganism(out, organism = outputOrganism);
  }

  out;
}


genericEnrichrCollection = function(enrichrFile, IDBase, nameBase, descriptionBase,
                                 enrichrClass, outputOrganism = NULL)
{
  if (grepl("\\.bz2$", enrichrFile)) {
    con = bzfile(enrichrFile);
  } else
    con = file(enrichrFile);
  fileLines = readLines(con);
  nSets = length(fileLines);

  split = strsplit(fileLines, split = "\t", fixed = TRUE);
  firstString = sapply(split, `[`, 1);

  IDs = spaste(IDBase, prependZeros(1:nSets, 5));
  setNames = spaste(nameBase, firstString);
  setDescriptions = spaste(descriptionBase, firstString);

  symbols = lapply(split, function(s) gsub(",[0-9.]*", "", s[-1]));
  printFlush("Converting symbols to Entrez...");
  humanEntrez = convert2entrez(organism = "human", symbol = symbols);
  geneSets = mymapply(newGeneSet,
               geneEntrez = humanEntrez,
               ID = IDs,
               name = setNames,
               description = setDescriptions,
             MoreArgs = list(
               organism = "human",
               geneSource = "Enrichr",
               geneEvidence = "other",
               source = spaste("Enrichr, ", enrichrClass),
               internalClassification = c("Enrichr", enrichrClass),
               groups = c("Enrichr", enrichrClass)));

  groups = list(newGroup("Enrichr", "Gene sets provided by Enrichr", "Enrichr"),
                newGroup(enrichrClass, spaste(enrichrClass, " gene sets, provided by Enrichr"), "Enrichr"));

  out = newCollection(geneSets, groups);
  if (!is.null(outputOrganism))
  {
     printFlush(spaste("Converting collection to ", outputOrganism))
     out = convertCollectionToOrganism(out, organism = outputOrganism);
  }
  out;
}


#======================================================================================================
#
# Collection from modules
#
#======================================================================================================

# description pattern can include %<colName> where <colName> is a column name in moduleInfo


collectionFromModules = function(labels, identifiers, moduleInfo,
   kME = NULL, kMEThreshold = 0.5,
   moduleCol = "module",
   descriptionCol = "enrichmentLabel",
   groupCol = NULL,
   addLabelToDescription = FALSE,
   IDBase,
   descriptionPattern, source, 
   organism, internalClassification, 
   groupList, 
   groups = sapply(groupList, getElement, "name"),
   lastModified = Sys.Date(),
   format = "%Y-%m-%d",
   alternateNames = character(0),
   externalDB = "",
   externalAccession = "",
   webLink = "",
   excludeLabels = 0,
   labelPrefix = "M")
{
   if (! moduleCol %in% names(moduleInfo)) stop("'moduleCol' must be one of names(moduleInfo)."); 
   if (! descriptionCol %in% names(moduleInfo)) stop("'descriptionCol' must be one of names(moduleInfo)."); 
   if (!is.null(groupCol) && ! groupCol %in% names(moduleInfo)) 
     stop("When given, 'groupCol' must be one of names(moduleInfo)."); 
   labelLevels = setdiff(sort(unique(labels)), excludeLabels);
   rows = match(labelLevels, moduleInfo[[moduleCol]]);
   if (any(is.na(rows))) stop("Not all module labels are in 'moduleInfo'.");

   setNames = moduleInfo[[descriptionCol]];
   if (addLabelToDescription) setNames = spaste(labelPrefix, moduleInfo[[moduleCol]], ": ", setNames);
   nModules = length(labelLevels)
   if (!is.null(groupCol)) {
     moduleGroups = strsplit(moduleInfo[[groupCol]], split = "|", fixed = TRUE)
   } else {
     moduleGroups = listRep(groups, nrow(moduleInfo));
   }
   sets = list();
   for (m in 1:nModules)
   {
     description1 = multiSub(spaste("%", colnames(moduleInfo)), as.character(unlist(moduleInfo[m, ])),
                             descriptionPattern);
     name1 = setNames[m];
     mm = moduleInfo[[moduleCol]] [ m ];
     if (!is.null(kME))
     {
       keep = labels==mm & kME >= kMEThreshold;
       description1 = spaste(description1, " Module members were restricted to kME >=", kMEThreshold);
       name1 = spaste(name1, " (kME >= ", kMEThreshold, ")");
     } else {
       keep = labels==mm;
     }
     keep = keep & !is.na(identifiers);
     sets[[m]] = newGeneSet(
         geneEntrez = identifiers[keep],
         geneEvidence = "HEP",
         geneSource = source,
         ID = spaste(IDBase, ".", prependZeros(m, nchar(nModules))),
         name = name1,
         description = description1,
         source = source,
         organism = organism,
         groups = unique(moduleGroups[[m]]),
         internalClassification = c(internalClassification, setNames[m]),
         lastModified =lastModified,
         format = fromat,
         alternateNames = alternateNames,
         externalDB = externalDB,
         externalAccession = externalAccession,
         webLink =webLink)
   }
   newCollection(sets, groupList);
}


#======================================================================================================
#
# unroll (flatten) a list and add index
#
#======================================================================================================

indexedFlattenedList = function(lst)
{
  n = length(lst);
  lengths = sapply(lst, length);
  index = do.call(c, mymapply(rep, 1:n, lengths));
  data.frame(index = index, data = unlist(lst));
}


#======================================================================================================
#
# complete multiData to same columns
#
#======================================================================================================

mtd.completeToCommonColumns = function(mtd, fillValue = NA, returnMatrix = NULL,
                               colOrder = NULL)
{
  if (is.null(returnMatrix)) 
  {
    returnMatrix = all(mtd.apply(mtd, is.matrix, mdaSimplify = TRUE));
    printFlush(paste("FYI: the returned object will contain", 
                     if (returnMatrix) "matrices." else "data frames."));
  }
  colUnion0 = multiUnion(mtd.apply(mtd, colnames, returnList = TRUE));
  if (any(!colOrder %in% colUnion0))
  {
    warning("Dropping entries of 'colOrder' that are not present in column names in 'mtd'.")
    colOrder = colOrder[colOrder %in% colUnion0];
  }
  colOrder = c(colOrder, setdiff(colUnion0, colOrder));
  
  order = match(colOrder, colUnion0);
  colUnion = colUnion0[order];

  nSamples = mtd.apply(mtd, nrow, mdaSimplify = TRUE);
  nCols = length(colUnion);
  # Copy all attributes and other information
  out = mtd;
  for (set in 1:length(mtd))
  {
    out1 = matrix(fillValue, nSamples[set], nCols);
    if (!returnMatrix) out1 = as.data.frame(out1);
    colnames(out1) = colUnion;
    rownames(out1) = rownames(mtd[[set]]$data);
    out1[ , match(colnames(mtd[[set]]$data), colUnion)] = mtd[[set]]$data;
    out[[set]]$data = out1;
  }
  out;
}

# Similar, but here we restrict to common columns.

mtd.restrictToCommonColumns = function(mtd)
{
  cols = multiIntersect(mtd.apply(mtd, colnames, returnList = TRUE));
  mtd.subset(mtd,, cols, permissive = TRUE);
}

#=========================================================================================================
#
# multiData.extendToCommon: an older version of mtd.completeToCommonColumns above.
#
#=========================================================================================================

multiData.extendToCommonColumns = function(multiData, verbose = 0)
{
  commonIDs = multiUnion(multiData2list(mtd.apply(multiData, colnames)));
  mtd.apply(multiData, function(data, IDs)
   {
     out = as.data.frame(matrix(NA, nrow(data), length(IDs)));
     rownames(out) = rownames(data);
     colnames(out) = IDs;
     out2data = match(IDs, colnames(data));
     inData = is.finite(out2data);
     out[, inData] = data[, out2data[inData]];
     out;
   }, commonIDs, mdaVerbose = verbose);
}


#======================================================================================================
#
# Restrict multiData expression data by minimum expression
#
#======================================================================================================

# Caution: this function assumes a multiData in the strict sense (same columns in all data sets).

mtd.highExpressedInAtLeastASubset = function(mtd, minValue, minNSamples, minNSets = 1)
{
  checkSets(mtd);
  he = mtd.apply(mtd, function(x) colSums(x>=minValue) >= minNSamples, mdaSimplify = TRUE);
  he.all = rowSums(he + 0) >= minNSets;
  if (all(he.all)) {
    mtd;
  } else 
    mtd.subset(mtd, , he.all);
}




#=====================================================================================================
#
# Combine model covariates from different data sets assuming that the confounders are
# distinct between data sets.
# These functions binarize categorical covariates to -1 and 1 to distinguish these values from 0 which really
# means the variable does not apply. This is necessary for proper centering before using in regression models.
#
#=====================================================================================================

combineCovariates.continuous = function(mtd, mtdNames = names(mtd), commonCovars = NULL,
                                binarizeCommonCovars = commonCovars, 
                                center = TRUE, ...)
{
  mtd.apply(mtd, function(x)
    if (!all(commonCovars %in% colnames(x))) 
       stop("Some 'commonCovars' are not present in all sets:\n",
            paste(commonCovars[!commonCovars %in% colnames(x)], collapse = ", ")));
  n = length(mtd);
  if (is.null(mtdNames)) mtdNames = spaste("Set.", 1:n);
  mtdNames = make.unique(mtdNames);
  mtd3 = mtd.mapply(function(x, name) 
  { 
      separate = !colnames(x) %in% commonCovars;
      if (any(separate))
          colnames(x)[separate] = spaste(colnames(x)[separate], ".", name); 
      if (center) 
        x[, separate] = apply(x[, separate, drop = FALSE], 2, function(x) as.numeric(scale(x, scale = FALSE)));
      x;
  }, mtd, mtdNames)
  out = mtd.rbindSelf(mtd.completeToCommonColumns(mtd3, fillValue = 0));
  if (length(binarizeCommonCovars)>0)
    out = binarizeCategoricalColumns(out, considerColumns = binarizeCommonCovars, minCount = 1,
                     includeLevelVsAll = TRUE, dropFirstLevelVsAll = TRUE, nameSep = "", 
                     val1 = -1, val2 = 1, ...);
  out;
}

combineCovariates.mixed = function(mtd, mtdNames = names(mtd), commonCovars = NULL, 
                                   binarizeCommonCovars = commonCovars, center = TRUE, ...)
{
  common = mtd.apply(mtd, function(x) colnames(x) %in% commonCovars);
  separate = mtd.apply(common, `!`);
  # Check that none of commonCovars have been specified in vain...
  allColNames = unique(unlist(mtd.apply(mtd, colnames, returnList = TRUE)));
  if (any(!commonCovars %in% allColNames))
    warning("combineCovariates.mixed: these 'commonCovars' are not colnames in mtd:\n",
            paste(commonCovars[!commonCovars %in% allColNames], collapse = ", "));

  mtd2 = mtd.mapply(binarizeCategoricalColumns, data = mtd, considerColumns = separate,
                   MoreArgs = list(..., includeLevelVsAll = TRUE, minCount = 1,
                          val1 = -1, val2 = 1,
                          dropUninformative = FALSE, dropFirstLevelVsAll = TRUE, nameSep = ""));

  combineCovariates.continuous(mtd2, mtdNames, commonCovars = commonCovars, 
                               binarizeCommonCovars = binarizeCommonCovars, 
                               center = center, ...);
}

centerNonZeroEntries = function(x)
{
  if (!is.null(dim(x))) return(apply(x, 2, centerNonZeroEntries));
  nonZero = replaceMissing(x!=0);
  x[nonZero] = c(scale(x[nonZero], scale = FALSE))
  x;
}
  

#=====================================================================================================
#
# Replace counts in a DESeq object based on the bicov weight of their VST
#
#=====================================================================================================

collapseGroups = function(groups)
{
  groups = as.matrix(groups);
  nc = ncol(groups);
  groups.fact = apply(groups, 2, function(x) as.numeric(factor(x)));
  nLevels = apply(groups.fact, 2, function(x) max(x));
  multipliers = c(1, cumprod(nLevels))[1:nc];
  group.mult = sweep((groups.fact-1), 2, multipliers, FUN = '*');
  rowSums(group.mult, na.rm = TRUE)+1;
}

deseqBicovWeights = function(object, blind = FALSE, 
  useApproximateVST = FALSE,
  approximateVST.minShift = 1,
  approximateVST.maxShift = 500,
  approximateVST.xWeightScale = 0,
  ...)
{
  if (useApproximateVST)
  {
    object = estimateSizeFactors(object);
    x = t(counts(object, normalized = TRUE));
    vst = approximateVST(x, min = approximateVST.minShift, max = approximateVST.maxShift, 
                         xWeightScale = approximateVST.xWeightScale);
  } else 
    vst = t(assay(varianceStabilizingTransformation(object, blind = blind)));
  matrixBicovWeights(vst, colData(object), ...)
}

deseqBicovWeightsFromMatrix = function( 
   counts,
   pheno,
   design,

   blind = FALSE,
   ...)
{
  object = DESeqDataSetFromMatrix.PL(t(counts), colData = pheno, design = design);
  deseqBicovWeights(object, blind = blind, ...);
}


matrixBicovWeights = function(data,
                             covars,
                             maxPOutliers = 0.05,
                             outlierReferenceWeight = 0.1,
                             groupBy = NULL,
                             #combineGroups = FALSE,
                             #minSamplesPerGroup = 10,
                             groupsForMinWeightRestriction = NULL,
                             minWeightInGroups = 0,
                             maxPropUnderMinWeight = 1,
                             defaultWeight = 1,
                             getFactors = FALSE)
{
  if (is.data.frame(data)) data = as.matrix(data);
  defaultFactor = sqrt(1-sqrt(defaultWeight));
  maxFactorInGroups = sqrt(1-sqrt(minWeightInGroups))
  if (length(groupBy)==0)
  {  
     factors = bicovWeightFactors(data, maxPOutliers = maxPOutliers, outlierReferenceWeight = outlierReferenceWeight,
                            defaultFactor = 1-defaultWeight);
  } else {
     col = match(groupBy, names(covars));
     if (any(is.na(col))) stop("Column(s)\n                 ", paste(groupBy[is.na(col)], collapse = ", "),
                               "\n not found in colData of 'object'.");
     res = empiricalBayesLM(data, removedCovariates = covars[, col], automaticWeights = "none",
              getOLSAdjustedData = FALSE,
              getFittedValues = FALSE,
              getEBadjustedData = FALSE)$residuals.OLS;
     dimnames(res) = dimnames(data);
     factors = bicovWeightFactors(res, maxPOutliers = maxPOutliers,
                              outlierReferenceWeight = outlierReferenceWeight, defaultFactor = defaultFactor)
  }  

  factors = abs(factors)

  # If groupsForMinWeightRestriction is a data frame or matrix (assumed to be 1-column), turn it into a vector
  if (!is.null(groupsForMinWeightRestriction)) 
     groupsForMinWeightRestriction = c(unlist(groupsForMinWeightRestriction));

  if (length(groupsForMinWeightRestriction)>0)
  {
    if (length(groupsForMinWeightRestriction)!=nrow(data)) 
       stop("'groupsForMinWeightRestriction' must be a vector of length equal the number of rows in 'data'.");
  }

  if (length(groupsForMinWeightRestriction)>0 && minWeightInGroups > 0)
  {
     fixCols = rep(FALSE, ncol(factors));
     refFactor = rep(maxFactorInGroups, ncol(factors))
     for (g in sort(unique(groupsForMinWeightRestriction)))
     {
       # Check that appropriate quantile of factor for each group is not below minWeightInGroups
       factorQuant = colQuantileC(factors[groupsForMinWeightRestriction==g, ], p = 1-maxPropUnderMinWeight)
       fix = factorQuant > maxFactorInGroups;
       fixCols[fix] = TRUE;
       refFactor[fix] = pmax(factorQuant[fix],  refFactor[fix]);
     }
     refFactorMat = matrix(refFactor, nrow(factors), ncol(factors), byrow = TRUE);
     factors = factors/refFactorMat * maxFactorInGroups;
     weights = bicovWeightsFromFactors(factors, defaultWeight = defaultWeight);
     attr(weights, "scaledColumnsToMeetMinWeight") = which(fixCols);
     attr(weights, "scaleFactorsToMeetMinWeights") = maxFactorInGroups/refFactorMat
  } else {
     weights = bicovWeightsFromFactors(factors, defaultWeight = defaultWeight)
  }
  
  if (getFactors) list(weights = weights, factors = factors) else weights;
}

  
replaceOutlierCounts = function(object, blind = FALSE, maxPOutliers = 0.05,
                                useApproximateVST = FALSE,
                                approximateVST.minShift = 1,
                                approximateVST.maxShift = 500,
                                approximateVST.xWeightScale = 0,
                                outlierReferenceWeight = 0.1, 
                                weightedReplace = TRUE,
                                replaceThreshold = if (weightedReplace) 0.5 else 0.1, 
                                returnDESeqDataSet = FALSE,
                                groupBy = NULL,
                                #combineGroups = FALSE,
                                #minSamplesPerGroup = 10,
                                groupsForMinWeightRestriction = NULL,
                                minWeightInGroups = 0,
                                maxPropUnderMinWeight = 1
                                )
{
  
  weights = deseqBicovWeights(object, blind = blind, 
                              useApproximateVST = useApproximateVST,
                              approximateVST.minShift = approximateVST.minShift,
                              approximateVST.maxShift = approximateVST.maxShift,
                              approximateVST.xWeightScale = approximateVST.xWeightScale,
                              maxPOutliers = maxPOutliers,
                              outlierReferenceWeight = outlierReferenceWeight,
                              groupBy = groupBy, #combineGroups = combineGroups,
                              #minSamplesPerGroup = minSamplesPerGroup,
                              groupsForMinWeightRestriction = groupsForMinWeightRestriction,
                              minWeightInGroups = minWeightInGroups,
                              maxPropUnderMinWeight = maxPropUnderMinWeight);

  groupBy = attr(weights, "group");

  replace = weights<replaceThreshold;
  normCounts = t(counts(object, normalized = TRUE));
  normFactors = object@colData$sizeFactor

  replacedCounts = t(counts(object))
  for (g in sort(unique(groupBy)))
  {
    inGroup = which(groupBy==g);
    replaceCols = which(colSums(replace[inGroup, , drop = FALSE]) > 0);

    if (length(replaceCols) > 0) for (c in replaceCols)
    {
      r1 = replace[inGroup, c];
      if (!all(r1))
      {
         rval = weighted.mean(normCounts[inGroup[!r1], c], weights[inGroup[!r1], c], na.rm = TRUE)*
                    normFactors[inGroup[r1]];
         if (weightedReplace)
         {
           w1  = weights[inGroup[r1], c]/replaceThreshold;
           if (any(w1>1)) stop("Internal error: invalid scaled weights. Sorry!");
           replacedCounts[inGroup[r1], c] = round(w1 * replacedCounts[inGroup[r1], c] + (1-w1) * rval);
         } else 
           replacedCounts[inGroup[r1], c] = round(rval);
      }
    }
  }
  if (returnDESeqDataSet) {
     # Force re-estimation of size factors (!)
     estimateSizeFactors(DESeqDataSetFromMatrix.PL(t(replacedCounts), design = design(object), 
                         colData = colData(object)));
  } else
     replacedCounts;
}

replicateIndex = function(x)
{
  n= length(x);
  levels = sort(unique(x));
  out = rep(1, n);
  for (l in levels)
  {
    out[x==l] = c(1:sum(x==l));
  } 
  out;
}



# All-in-one replacement function for most common use

outlierReplacementData = function(
   counts,
   pheno,
   design,

   DESeqArgs = list(minReplicatesForReplace = Inf),

   replaceThreshold = 0.5, 
   weightedReplace = TRUE, 
   maxPOutliers = 0.10,
   outlierReferenceWeight = 0.1,
   groupBy = NULL,
   #minSamplesPerGroup = 10,
   #combineGroups = TRUE,
   groupsForMinWeightRestriction = NULL,
   minWeightInGroups = 0.9,
   maxPropUnderMinWeight = 0.4,

   returnWeights = TRUE,
   returnVSTdata = TRUE,
   blind = FALSE,
   calibrateLowestValue = TRUE,
   useApproximateVST = FALSE,
   approximateVST.minShift = 1,
   approximateVST.maxShift = 500,
   approximateVST.xWeightScale = 0,
 
   keepGenes = NULL,
   ...)  # The dots are currently not used, only there for compatibility with other functions
{
  expr.ds = DESeqDataSetFromMatrix.PL(countData = t(counts),
                     colData = pheno,
                     design = as.formula(design));

  if (useApproximateVST) {
    deseq1 = estimateSizeFactors(expr.ds);
  } else 
    deseq1 = do.call(DESeq, c(list(object = expr.ds), DESeqArgs));
    
  counts.repo = replaceOutlierCounts(deseq1, 
           useApproximateVST =  useApproximateVST,
           approximateVST.minShift = approximateVST.minShift,
           approximateVST.maxShift = approximateVST.maxShift,
           approximateVST.xWeightScale = approximateVST.xWeightScale,
           replaceThreshold = replaceThreshold, 
           weightedReplace = weightedReplace, 
           maxPOutliers = maxPOutliers,
           outlierReferenceWeight = outlierReferenceWeight,
           groupBy = groupBy,
           #minSamplesPerGroup = minSamplesPerGroup,
           #combineGroups = combineGroups,
           returnDESeqDataSet = FALSE, 
           groupsForMinWeightRestriction = groupsForMinWeightRestriction,
           minWeightInGroups = minWeightInGroups,
           maxPropUnderMinWeight = maxPropUnderMinWeight);

   dimnames(counts.repo) = dimnames(counts);

   if (returnWeights)
   {
      weightData = deseqBicovWeights(deseq1,
           useApproximateVST =  useApproximateVST,
           approximateVST.minShift = approximateVST.minShift,
           approximateVST.maxShift = approximateVST.maxShift,
           approximateVST.xWeightScale = approximateVST.xWeightScale,
           maxPOutliers = maxPOutliers,
           outlierReferenceWeight = outlierReferenceWeight,
           groupBy = groupBy,
           #minSamplesPerGroup = minSamplesPerGroup,
           #combineGroups = combineGroups,
           groupsForMinWeightRestriction = groupsForMinWeightRestriction,
           minWeightInGroups = minWeightInGroups,
           maxPropUnderMinWeight = maxPropUnderMinWeight,
           getFactors = TRUE);
      bicovWeightsForReplacement = weightData$weights;
      bicovFactorsForReplacement = weightData$factors;
      dimnames(bicovFactorsForReplacement) = dimnames(bicovWeightsForReplacement) = dimnames(counts);
   } else {
     bicovWeightsForReplacement = NULL;
     bicovFactorsForReplacement = NULL;
   }

   # Variance-stabilize the replaced counts
   if (returnVSTdata)
   {
     #ds.repl= DESeqDataSetFromMatrix.PL(
     #                countData = t(counts.repo),
     #                colData = pheno,
     #                design = as.formula(design));
#
#     expr.vs.repl = varianceStabilizingTransformation(ds.repl, blind = blind);
#
#     expr.repo.all = t(assay(expr.vs.repl));
#     dimnames(expr.repo.all) = dimnames(counts.repo);

     expr.repo.all = varianceStabilizedData.PL(counts.repo, pheno = pheno, design = design,
                            blind = blind,
                            calibrateLowestValue = calibrateLowestValue,
                            useApproximateVST = useApproximateVST, 
                            approximateVST.minShift = approximateVST.minShift,
                            approximateVST.maxShift = approximateVST.maxShift,
                            approximateVST.xWeightScale = approximateVST.xWeightScale);

                                            
     if (!is.null(keepGenes))
     {
        if (is.logical(keepGenes) || is.numeric(keepGenes))
        {
          expr.repo = expr.repo.all[, keepGenes]
        } else 
          expr.repo = expr.repo.all[, match(keepGenes, colnames(expr.repo.all))];
     } else
        expr.repo = expr.repo.all;
   } else 
      expr.repo = NULL;

   list(replacedCounts = counts.repo,
        replacementWeights = bicovWeightsForReplacement,
        replacementFactors = bicovFactorsForReplacement,
        replacedVST = expr.repo);
}

# Outlier replacement function for general data, not necessarily counts.
# In the current version the replacement is performed 

outlierReplacementData.general = function(
   data,
   pheno,

   useVarianceStabilization = TRUE,
   approximateVST.minShift = 1,
   approximateVST.maxShift = 500,
   approximateVST.existingLogBase = 2,
   approximateVST.existingLogShift = 1,   
   approximateVST.xWeightScale = 0,

   replaceThreshold = 0.5, 
   weightedReplace = TRUE, 
   maxPOutliers = 0.10,
   outlierReferenceWeight = 0.1,
   groupBy = NULL,
   #minSamplesPerGroup = 10,
   #combineGroups = TRUE,
   groupsForMinWeightRestriction = NULL,
   minWeightInGroups = 0.9,
   maxPropUnderMinWeight = 0.4,

   returnWeights = TRUE,
   returnVSTdata = TRUE,
   keepGenes = NULL,
   ...)  # The dots are currently not used, only there for compatibility with other functions
{

  if (useVarianceStabilization)
  {
    expr.vst = approximateVST(x = data, min = approximateVST.minShift, max = approximateVST.minShift,
                 xIsLogTransformed = !is.null(approximateVST.existingLogBase), 
                 xLogBase = approximateVST.existingLogBase, xLogShift = approximateVST.existingLogShift, epsilon = 0.01,
                 xWeightScale = approximateVST.xWeightScale);
  } else 
    expr.vst = data;

  weightData = matrixBicovWeights(expr.vst,
           covars = pheno,
           maxPOutliers = maxPOutliers,
           outlierReferenceWeight = outlierReferenceWeight,
           groupBy = groupBy, #combineGroups = combineGroups,
           #minSamplesPerGroup = minSamplesPerGroup,
           groupsForMinWeightRestriction = groupsForMinWeightRestriction,
           minWeightInGroups = minWeightInGroups,
           maxPropUnderMinWeight = maxPropUnderMinWeight,
           getFactors = TRUE);

  weights = weightData$weights;
  factors = weightData$factors;
  dimnames(weights) = dimnames(factors) = dimnames(data);
  groupBy = attr(weights, "group");

  replace = weights<replaceThreshold;

  replacedData = data;
  for (g in sort(unique(groupBy)))
  {
    inGroup = which(groupBy==g);
    replaceCols = which(colSums(replace[inGroup, , drop = FALSE]) > 0);
    if (length(replaceCols) > 0) for (c in replaceCols)
    {
      r1 = replace[inGroup, c];
      if (!all(r1))
      {
         rval = weighted.mean(data[inGroup[!r1], c], weights[inGroup[!r1], c], na.rm = TRUE)
         if (weightedReplace)
         {
           w1  = weights[inGroup[r1], c]/replaceThreshold;
           if (any(w1>1)) stop("Internal error: invalid scaled weights. Sorry!");
           replacedData[inGroup[r1], c] = w1 * replacedData[inGroup[r1], c] + (1-w1) * rval;
         } else
           replacedData[inGroup[r1], c] = rval;
      }
    }
  }

  # Variance-stabilize the replaced data
  if (returnVSTdata)
  {
    if (useVarianceStabilization)
    {
      replacedVST.all = approximateVST(replacedData, min = approximateVST.minShift, max = approximateVST.minShift,
                 xIsLogTransformed = !is.null(approximateVST.existingLogBase), 
                 xLogBase = approximateVST.existingLogBase, xLogShift = approximateVST.existingLogShift, epsilon = 0.01,
                 xWeightScale = approximateVST.xWeightScale)
    } else
      replacedVST.all = replacedData;

    if (!is.null(keepGenes))
    {
        if (is.logical(keepGenes) || is.numeric(keepGenes))
        {
          replacedVST = replacedVST.all[, keepGenes]
        } else 
          replacedVST = replacedVST.all[, match(keepGenes, colnames(replacedVST.all))];
     } else
        replacedVST = replacedVST.all;
   } else 
      replacedVST = NULL;

   list(replacedData = replacedData,
        replacementWeights = weights,
        replacementFactors = factors,
        replacedVST = replacedVST);
}


calibrateLowestValue = function(data, reference, cut = NULL, zero = 0.01, subtractMin = TRUE)
{
  data.f = as.numeric(as.matrix(data));
  ref.f = as.numeric(as.matrix(reference));

  if (length(data.f)!=length(ref.f))
    stop("'data' and 'reference' must have the same dimensions.");

  if (!is.null(cut))
  {
    keep = ref.f < cut
  } else 
    keep = rep(TRUE, length(ref.f));

  minInd = ref.f <= zero;
  keep[minInd] = FALSE;

  if (sum(keep) < 5)
  {
    data.predicted = data.f
  } else {
    data.reg = data.f[keep];
    ref.reg = ref.f[keep];

    fit = lm(data.reg~ref.reg);

    predicted = predict(fit, newdata = data.frame(ref.reg = ref.f[minInd]));

    data.predicted = data.f;
    data.predicted[minInd] = predicted;
  }
  if (subtractMin) data.predicted = data.predicted - min(data.predicted, na.rm = TRUE);
  dim(data.predicted) = dim(data);
  dimnames(data.predicted) = dimnames(data);

  data.predicted;
}

varianceStabilizedData.PL = function(
  counts, pheno, design, blind = FALSE,
  calibrateLowestValue = TRUE,
  cut = 2.5, zero = 0.1, subtractMin = FALSE,
  useApproximateVST = FALSE,
  approximateVST.minShift = 1,
  approximateVST.maxShift = 500,
  approximateVST.xWeightScale = 0,
  xWeightScale = 0
)
{
  ds= DESeqDataSetFromMatrix.PL(
                     countData = t(counts),
                     colData = pheno,
                     design = as.formula(design));

  ds = estimateSizeFactors(ds);
  if (useApproximateVST)
  {
    x = t(counts(ds, normalized = TRUE));
    expr = approximateVST(x, approximateVST.minShift, max = approximateVST.maxShift,
      xIsLogTransformed = FALSE, xWeightScale = approximateVST.xWeightScale)
  } else {
    vst = varianceStabilizingTransformation(ds, blind = blind);
    expr = t(assay(vst));
    dimnames(expr) = dimnames(counts);
    if (calibrateLowestValue)
    {
      counts.log = log2(t(counts(ds, normalized = TRUE))+1);
      expr = calibrateLowestValue(expr, counts.log,
                       cut = cut, zero = zero, subtractMin = subtractMin);
    }
  }
  expr;
}

DESeqDataSetFromMatrix.PL = function(countData, colData, ...)
{
  if (ncol(countData)!=nrow(colData)) browser()
  rownames(colData) = colnames(countData);
  DESeqDataSetFromMatrix(countData, colData, ...);
}



#========================================================================================================
#
# replaceMissing
#
#========================================================================================================

replaceMissing = function(x, replaceWith)
{
  if (missing(replaceWith))
  {
    if (is.logical(x)) {
      replaceWith = FALSE
    } else if (is.numeric(x)) {
      replaceWith = 0;
    } else if (is.character(x)) {
      replaceWith = ""
    } else stop("Need 'replaceWith'.");
  }
  x[is.na(x)] = replaceWith;
  x;
}

removeMissing = function(x) x[!is.na(x)];


#========================================================================================================
#
# leftovers of multiGSub, multiSub
#
#========================================================================================================

multiSubr = function(patterns, x, ...) multiSub(patterns, rep("", length(patterns)), x, ...);
multiGSubr = function(patterns, x, ...) multiGSub(patterns, rep("", length(patterns)), x, ...);
    
grep.first = function(patterns, x, ...)
{
  sapply(patterns, function(p) min(grep(p, x, ...)));
}
  
#========================================================================================================
#
# Network screening (NTPS) functions
#
#========================================================================================================

getColFromDataFrame = function(df, col)
{
  if (length(col)!=1) stop("This function works with a scalar 'col' only.")
  if (is.character(col))
  {
    ncol = match(col, colnames(df))
  } else ncol = as.numeric(col);

  if (!is.finite(ncol)) stop("Column ", col, " could not be found among column names of 'df' which are:\n",
                             formatLabels(paste(colnames(df), collapse = "; "), 80));

  if (ncol < 0 || ncol > ncol(df)) stop("'col' is out of range.");

  if (ncol==0) 
  {
    if (is.null(rownames(df))) 
      stop("'df' does not have rownames.");
    out = rownames(df);
  } else
    out = df[, ncol];
  out;
}

# Convert an unsigned variable and a separate sign to a Z statistic 
  
signedScale = function(x, dir)
{
  sign = sign(dir);
  Z = scale(x*sign);
  zero = mean(Z[x==0], na.rm = TRUE);
  Z = Z-zero
  Z[is.na(Z)] = 0;
  Z;
}

 
signedScale.withMatching = function(x, Z,
                   x.IDcol = 1,
                   x.xCol = 2,
                   Z.IDcol = 1,
                   Z.Zcol = 2
                   )
{
  ID.x = getColFromDataFrame(x, x.IDcol);
  ID.Z = getColFromDataFrame(Z, Z.IDcol);

  x.matched = getColFromDataFrame(x, x.xCol)[match(ID.Z, ID.x)];
  data.frame(ID = ID.Z, x.matched = x.matched, Z = signedScale(x.matched, getColFromDataFrame(Z, Z.Zcol)));
}

# Combine Z statistics

combineZ = function(Zs, weights = NULL)
{
  args = as.matrix(as.data.frame(Zs));
  if (is.null(weights)) weights = rep(1, ncol(args));
  weightMat = matrix(weights, nrow(args), ncol(args), byrow = TRUE);
  finArgs = is.finite(args);
  weightMat[!finArgs] = 0;  # Missing observation do not contribute to the weighted average
  hasFinite = rowSums(finArgs) > 0;
  coeffs = 1/sqrt(apply(args, 2, var, na.rm = TRUE));
  out = rowSums(matrix(coeffs, nrow(args), ncol(args), byrow = TRUE) * weightMat * args, na.rm = TRUE)/
          sqrt(sum(coeffs^2 * weights^2, na.rm = TRUE));
  out[!hasFinite] = NA
  out;
}

# One NTPS score per module.

NTPSbyModule = function(kME, 
                        association, 
                        causalConnectivity = NULL, 
                        moduleSummary, 
                        weights = NULL,
                        zKMEPattern = "Z.ModMembership.in.M.",
                        associationCol = "Z.meta.for.continuousQ",
                        moduleAssociationCol = "Z.metaAnalysis.for.Q",
                        causalConnectivityCol = "Z",
                        # output options
                        rankKMEPattern = sub("Z", "rank", zKMEPattern))
{
  kMECols = grep(zKMEPattern, colnames(kME));
  kME = kME[, kMECols, drop = FALSE];
  nModules = ncol(kME);
  moduleLabels = as.numeric(sub(zKMEPattern, "", colnames(kME)));
  moduleLabels2 = moduleSummary$module;
  if (!isTRUE(all.equal(moduleLabels, moduleLabels2))) stop("Module labels disagree.");

  if (is.null(weights)) weights = rep(1, 3);

  haveCausal = !is.null(causalConnectivity);
  if (haveCausal)
    zCausal = getColFromDataFrame(causalConnectivity, causalConnectivityCol);

  if (nrow(association)!=nrow(kME)) 
    stop("Length of 'association' and number of rows in 'kME' must be the same.");

  if (haveCausal && nrow(association)!=nrow(causalConnectivity)) 
    stop("Length of 'association' and 'zCasual' must be the same.");

  moduleAssn = getColFromDataFrame(moduleSummary, moduleAssociationCol);
  moduleSign = sign(moduleAssn);

  associationZ = getColFromDataFrame(association, associationCol);

  ntps0 = kME;
  rank0.signKME = rank0.absKME = kME;
  colnames(ntps0) = spaste("Z.NTPS.for.M.", moduleLabels);
  for (m in 1:nModules)
  {
    d1= if (!haveCausal) data.frame(associationZ, moduleSign[m]* kME[, m]) else
                    data.frame(associationZ, moduleSign[m]* kME[, m], zCausal);
    ntps0[, m] = combineZ(d1, weights = weights[1:ncol(d1)]);
    rank0.signKME[, m] = rank(-moduleSign[m]*ntps0[, m], na.last = TRUE)
    rank0.absKME[, m] = rank(-abs(ntps0[, m]), na.last = TRUE)
  }
  colnames(rank0.signKME) = spaste("rank.NTPS.for.M.", moduleLabels);
  colnames(rank0.absKME) = spaste("rankAbs.NTPS.for.M.", moduleLabels);
  rank.kME = apply(-kME, 2, rank);
  colnames(rank.kME) = gsub(zKMEPattern, rankKMEPattern, colnames(kME));
  if (haveCausal) colnames(causalConnectivity) = c("ID", "causalConnectivity", "Z.causal");
  if (haveCausal)
  {
     causal.rank = apply(causalConnectivity[, -1], 2, function(x) rank(replaceMissing(x, 0), na.last = TRUE));
     causal.absRank = apply(causalConnectivity[, -1], 2, function(x) 
                               rank(-abs(replaceMissing(x, 0)), na.last = TRUE));
     causal.x = interleave(list(causalConnectivity[, -1], causal.rank, causal.absRank),
                           nameBase = c("", "rank.", "rankAbs."), sep = "")[, -2];
     data.frame(Z.Association = associationZ, rank.Z.Association = rank(associationZ, na.last = TRUE),
                rankAbs.Z.Association = rank(-abs(associationZ), na.last = TRUE), 
                causal.x,
        interleave(list(ntps0, rank0.signKME, rank0.absKME, kME, rank.kME), nameBase = rep("", 5), sep = ""))
  } else
     data.frame(Z.Association = associationZ, rank.Z.Association = rank(associationZ, na.last = TRUE),
                rankAbs.Z.Association = rank(-abs(associationZ), na.last = TRUE),
        interleave(list(ntps0, rank0.signKME, rank0.absKME, kME, rank.kME), nameBase = rep("", 5), sep = ""));
}

# One combined NTPS score for all modules.

NTPS.singleScoreForAllModules = function(kME,
                        association,
                        causalConnectivity = NULL,
                        moduleSummary,
                        weights = NULL,
                        kmePower = 6,
                        zKMEPattern = "Z.ModMembership.in.M.",
                        associationCol = "Z.meta.for.continuousQ",
                        moduleAssociationCol = "Z.metaAnalysis.for.Q",
                        causalConnectivityCol = "Z")
{
  kMECols = grep(zKMEPattern, colnames(kME));
  kME = kME[, kMECols, drop = FALSE];
  nModules = ncol(kME);
  moduleLabels = as.numeric(sub(zKMEPattern, "", colnames(kME)));
  moduleLabels2 = moduleSummary$module;
  if (!isTRUE(all.equal(moduleLabels, moduleLabels2))) stop("Module labels disagree.");

  if (is.null(weights)) weights = rep(1, 3);

  haveCausal = !is.null(causalConnectivity);
  if (haveCausal)
    zCausal = getColFromDataFrame(causalConnectivity, causalConnectivityCol);

  nGenes = nrow(kME);
  if (nrow(association)!=nGenes)
    stop("Length of 'association' and number of rows in 'kME' must be the same.");

  if (haveCausal && nrow(association)!=nrow(causalConnectivity))
    stop("Length of 'association' and 'zCasual' must be the same.");

  moduleAssn = getColFromDataFrame(moduleSummary, moduleAssociationCol);
  moduleSign = sign(moduleAssn);

  associationZ = getColFromDataFrame(association, associationCol);

  geneWeights = abs(kME)^kmePower
  moduleAssnMat = matrix(moduleAssn, nGenes, nModules, byrow = TRUE);
  
  moduleContribution = rowSums( kME * moduleAssnMat * geneWeights)/rowSums(geneWeights)

  if (haveCausal)
  {
     Zcomponents = data.frame(Z.associationWithQ = associationZ,
                           Z.kME.moduleSignificance = moduleContribution,
                           Z.causal = zCausal);
  } else {
     Zcomponents = data.frame(Z.associationWithQ = associationZ,
                              Z.kME.moduleSignificance = moduleContribution);
  }
  df.Z = data.frame(Z.NTPS = combineZ(Zcomponents, weights = weights[1:ncol(Zcomponents)]), Zcomponents);
  #df.Z.r = df.Z;
  #if (haveCausal)
  #  df.Z.r$Z.causal = replaceMissing(df.Z.r$Z.causal, 0);
  df.Ord = apply(df.Z, 2, rank, na.last = TRUE);
  colnames(df.Ord) = sub("^Z.", "rank.", colnames(df.Z));
  df.rankAbs = apply(df.Z, 2, function(x) rank(-abs(x), na.last = TRUE));
  colnames(df.rankAbs) = sub("^Z.", "rankAbs.", colnames(df.Z));
  interleave(list(df.Z, df.Ord, df.rankAbs), nameBase = rep("", 3), sep = "");
}







#======================================================================================================
#
# Insert columns into a data frame
#
#======================================================================================================

insertColumns = function(data, columnNames, insertAfter, fillValue = NA, returnDataFrame = TRUE)
{
  n = length(columnNames);
  if (length(insertAfter)!= n) stop("Lengths of 'columnNames' and 'insertAt' must be the same.");

  if (is.numeric(insertAfter)) 
  {
    after = insertAfter;
  } else {
    after = match(insertAfter, colnames(data));
    after[insertAfter=="^"] = 0;
    after[insertAfter=="$"] = ncol(data);
  }
  if (any(!is.finite(after))) stop("all entries in 'insertAfter' must be finite.");
  if (any(after < 0 | after > ncol(data))) 
     stop("all entries in 'insertAfter' must be between 0 and the number of columns in 'data'.");

  if (any(duplicated(after))) stop("Some entries in 'insertAfter' are duplicated; all must be unique.");
  order = order(after);
  columnNames = columnNames[order];
  after = after[order];
  index = 0;
  out.list = list();
  for (i in 1:n)
  {
    out1 = matrix(fillValue, nrow(data), length(columnNames[[i]]));
    colnames(out1) = columnNames[[i]];
    if (returnDataFrame) out1 = as.data.frame(out1);

    if (index < after[i]) out.list = c(out.list, list(data[, (index+1):after[i], drop = FALSE]));
    out.list = c(out.list, list(out1));
    index = after[i];
  }
  if (index < ncol(data)) out.list = c(out.list, list(data[, (index+1):ncol(data), drop = FALSE]));

  if (returnDataFrame) {
    out = do.call(data.frame, out.list)
  } else 
    out = do.call(cbind, out.list);

  out;
}
  
  
  
#==========================================================================================================
#
# Filter out "not good" genes
#
#========================================================================================================== 

filterSamplesGenes = function(data, ...)
{
  gsg = goodSamplesGenes(data, ...);
  if (gsg$allOK) data else data[gsg$goodSamples, gsg$goodGenes];
}

filterSamplesGenesMS = function(multiData, ...)
{
  gsg = goodSamplesGenesMS(multiData, ...);
  if (gsg$allOK) multiData else mtd.subset(multiData, gsg$goodSamples, gsg$goodGenes);
}

#=========================================================================================================
#
# compress methylation data that are outside the usual interval [0, 1] to the interval
#
#=========================================================================================================

compressRange = function(data, rangeMin, rangeMax, dontExpand = TRUE)
{
  data = as.matrix(data);
  dataMin = min(data, na.rm = TRUE);
  dataMax = max(data, na.rm = TRUE);

  if (dontExpand)
  {
    if (rangeMin < dataMin) rangeMin = dataMin;
    if (rangeMax > dataMax) rangeMax = dataMax;
  }

  groups = base::findInterval(data, vec = c(rangeMin, rangeMax), rightmost.closed = TRUE) + 1
  groupSizes = tabulate(groups)
  heights = c(0, cumsum(groupSizes)/sum(groupSizes)) * (rangeMax - rangeMin) + rangeMin;

  breaks = c(min(dataMin, rangeMin), rangeMin, rangeMax, max(dataMax, rangeMax));

  # Adjust data separately in each of the three groups
  for (g in 1:3)
  {
    keep = replaceMissing( groups==g);
    if (any(keep))
      data[keep] = (data[keep] - breaks[g])/(breaks[g+1] - breaks[g]) * ( heights[g+1] - heights[g]) + 
                                      heights[g];
  }

  data;
}


#=========================================================================================================
#
# "Anti-geometric" mean: log of the mean of exponentiated values
#
#=========================================================================================================

mean.exp = function(x, na.rm = TRUE)
{
  log2(mean(2^x, na.rm = TRUE));
}


#=========================================================================================================
#
# loadOrInstall.BioC
#
#=========================================================================================================

loadOrInstall.BioC = function(pkgs, lib.loc = .Library, lib = .Library)
{
  success = sapply(pkgs, require, lib.loc = lib.loc, character.only = TRUE);
  if (any(!success))
  {
    install = pkgs[!success]
    source("http://bioconductor.org/biocLite.R");
    biocLite(install, suppressUpdates=TRUE, suppressAutoUpdate=TRUE, lib.loc = lib.loc, lib = lib);
    success2 = sapply(install, require, lib.loc = lib, character.only = TRUE);
    if (any(!success2))
      stop("Could not install or load package(s) ", paste(install[!success2], collapse = ", "));
  }
}

#=========================================================================================================
#
# histLog
#
#=========================================================================================================

histLog = function(x, ..., plot = TRUE)
{
  h = hist(x, ..., plot = FALSE);
  h$counts = log10(h$counts + 1);
  if (plot) plot(h, ...);
  h;
}

#=========================================================================================================
#
# scientific format: note the application in text() below.
#
#=========================================================================================================

prettyExpFormat = function(x, lead = NULL, brackets = FALSE)
{
  mult = "\u00D7";
  x = as.character(x);
  split = strsplit(x, split = "e");
  base = sapply(split, `[`, 1);
  exp = sapply(split, `[`, 2);
  exp[is.na(exp)] = "0";
  if (length(lead)>0 && length(lead)!=length(x)) 
    stop("When 'lead' is given, it must have the same length as 'x'.");
  if (length(lead) > 0) brackets = TRUE;
  if (length(lead)==0) lead = rep("", length(x));
  lead = as.character(lead);
  
  if (brackets)
  {
    out = mymapply(function(b, e, l) 
       { 
         b = as.numeric(b); e = as.numeric(e);
         out2 = if (e==0) substitute(group("(", b, ")"),  list(b=b)) else 
                 substitute(group("(", b %*% 10^e, ")"), list(b=b, e=e)) 
         if (l!="") out2 = substitute(pref~out, list(pref = as.numeric(l), out = out2));
         as.expression(out2);
       }, base, exp, lead);
  } else
    out = mymapply(function(b, e, lead) {b = as.numeric(b); e = as.numeric(e);
               if (e==0) b else substitute(b %*% 10^e, list(b=b, e=e))}, base, exp, lead);
  out;
}
  

if (FALSE)
{
plot(1:10, type = "n")
mymapply(function(x, y, t) text(x, y, t), 1:10, 1:10, prettyExpFormat(1e-8*(1:10)))
}


#=========================================================================================================
#
# Pick one representative probe for each gene from methylation data.
#
#=========================================================================================================

getEntrezFromIllumina450Annotation = function(annot)
{
  warning(spaste("Unless reproducing old results, please use the corresponding function \n",
          "   entrezFromMultiSymbol in anRichment."), immediate. = TRUE);
  if (is.data.frame(annot)) symbol = annot$UCSCRefGeneName else symbol = annot;
  if (is.null(symbol)) stop("Could not find gene symbols in 'annot'.");

  split = lapply(strsplit(symbol, split = ";"), unique)
  symbols.indexed = indexedFlattenedList(split);
  
  entrez.indexed = convert2entrez(organism = "human", symbol = symbols.indexed$data);

  var = tapply(entrez.indexed, symbols.indexed$index, var, na.rm = TRUE)
  var[is.na(var)] = 0;
  keep = as.numeric(names(var)[var==0]);
  
  entrez.all = rep(NA, length(symbol));
  keep2index = match(keep, symbols.indexed$index)
  entrez.all[keep] = entrez.indexed[keep2index]; 
  entrez.all;
}


#========================================================================================================
#
# Summarize several enrichment analyses together.
#
#========================================================================================================

# Given several enrichment analyses for the same set of query classes, produce a unified table of highest
# enriched terms. Also cluster gene sets and report the clustering results.

summarizeEnrichmentTables = function(enrMtd, collection, minKeep = 10,
                                     deepSplit = 0, minClusterSize = 3, ...)
{
  if (!isMultiData(enrMtd, strict = FALSE)) enrMtd = list2multiData(enrMtd);
  pTables = mtd.apply(enrMtd, getElement, "pValues", returnList = TRUE);
  commonModules = multiUnion(lapply(pTables, colnames));
  pTables2 = lapply(pTables, function(x) x[, match(commonModules, colnames(x))]);
  pMin = do.call(pmin, pTables2);
  setNames = dataSetNames(collection);
  shortNames = sapply(collection$dataSets, getElement, "shortName");
  nSets = length(setNames);
  modules = colnames(pMin);
  nModules = length(modules);
  keepSets = apply(pMin, 2, function(x)
  {
     order = order(x);
     x.ordered = x[order];
     keep = sort(unique(c(1:minKeep, which(x[order] < 0.05/nSets))));
     order[keep];
  })

  if (!is.null(dim(keepSets))) keepSets = as.data.frame(keepSets);

  printFlush("Calculating gene set similarity matrices. This may take a few moments.");
  geneSetSimilarities = lapply(keepSets, function(ks)
  {
    col1 = collection;
    col1$dataSets = col1$dataSets[ks];
    dataSetSimilarity(col1,
                 type = "product",
                 namesFrom = "name",
                 productNormalization = "mean",
                 TOMType = "none", verbose = 0); 
  });
  trees = lapply(geneSetSimilarities, function(x)  hclust(as.dist(1-x), method = "average"));
  clusters = mymapply(function(t, s) cutreeDynamic(t, distM = 1-s, verbose = 0,
                       deepSplit = deepSplit, minClusterSize = minClusterSize, ... ),
                     trees, geneSetSimilarities);

  geneSetConnectivities = mymapply(function(sim, mod) 
             intramodularConnectivity(sim, mod, scaleByMax = TRUE)$kWithin,
                                 geneSetSimilarities, clusters);

  # For each module and each cluster of enriched gene sets, find the centroid (highest connectivity) gene set.

  centroids = mymapply(function(gsc, cl, names)
    {
      clL = sort(unique(cl));
      clL = clL[clL!=0];
      centroids = rep("", length(clL));
      for (c in clL)
      {
        inCl = which(cl==c);
        centroids[c] = names[ inCl[ which.max(gsc[inCl])]];
      }
      centroids;
    }, geneSetConnectivities, clusters, lapply(geneSetSimilarities, colnames)); 

  # Put together output
  
  if (FALSE)
  {
    n1 = shortenStrings(colnames(geneSetSimilarities[[1]]), maxLength = 50)
    plotDendroAndColors(trees[[1]], clusters[[1]], dendroLabels = n1, cex.dendroLabels = 0.5)
  }

  keepPvalues = mymapply(function(col, index)
   {
      out = sapply(pTables2, function(x) x[index, col]);
      colnames(out) = spaste("p.", colnames(out));
      which = colnames(out)[ apply(out, 1, which.min)];
      min =  apply(out, 1, min, na.rm = TRUE)
      data.frame(MinCol = which, minPvalue = min, out);
   }, 1:nModules, keepSets);

  moduleTables = mymapply(function(mod, index, pValues, 
                                   cl, ctrs)
  {
    n = length(index);
    cl[cl==0] = NA;
    out = data.frame(Module = rep(modules[mod], n),
          Rank = 1:n,
          geneSet = setNames[index],
          geneSetCluster = cl,
          clusterCentroid = ctrs[cl],
          pValues);
    rownames(out) = NULL;
    out;
  }, 1:nModules, keepSets, keepPvalues, clusters, centroids);

  do.call(rbind, moduleTables);
}


#========================================================================================================
#
# "intelligent" cbind
#
#========================================================================================================

# This cbind checks for common columns in the two data frames, and tries to merge them

icbind2 = function(x1, x2, checkValues)
{
  names1 = make.unique(make.names(colnames(x1)));
  names2 = make.unique(make.names(colnames(x2)));

  common = intersect(names1, names2);

  out = list();
  for (c in common)
  {
    c1 = match(c, names1);
    c2 = match(c, names2);
    d1 = x1[, c1];
    d2 = x2[, c2];
    d = d1;
    if (checkValues && !all(d1==d2, na.rm = TRUE))
      stop("Problem merging columns ", c, ": some non-missing values are different.\n",
           "      Column name in x1: ", colnames(x1)[c1], 
           "\n      Column name in x2: ", colnames(x2)[c2]);
    d[is.na(d1)] = d2[is.na(d1)];
    out = c(out, list(d));
  }
  names(out) = common;
  notCommon1 = !names1 %in% common;
  notCommon2 = !names2 %in% common;
  cbind(as.data.frame(out), x1[, notCommon1], x2[, notCommon2]);
}



#==========================================================================================================
#
# Load/build a standard set of collections for anRichment
#
#==========================================================================================================

standardCollections = function(organism,
  getHDWGCNA = TRUE,
  getAgingWGCNA = TRUE,
  getMiRNATargets = TRUE,
  getEnrichr.Aging = TRUE,
  genomicSpacings = c(5e6),
  buildExternal = TRUE,
  rootDir = NULL,
  ...
  )
{
  wd = getwd();
  if (is.null(rootDir)) 
  {
    rootDir = file.path(sub("/Work/.*", "", wd), "Work");
    printFlush(spaste("Using root path: ", rootDir));
  }
  library(anRichment)
  source(file.path(rootDir, "RLibs/inHouseGeneAnnotation/anRichmentMethods/R/enrichmentAnalysis.R"));
  source(file.path(rootDir, "RLibs/inHouseGeneAnnotation/anRichment/R/loadCollections.R"));

  collections = allCollections(organism = organism,
                     buildExternal = buildExternal,
                     genomicSpacings = genomicSpacings,
                     merge = FALSE);

  if (getHDWGCNA)
  {
    library(HuntingtonDiseaseWGCNACollection);
    collections = c(collections, list(HD.WGCNA = HuntingtonDiseaseWGCNACollection(organism)));
  }

  if (getMiRNATargets)
  {
     if (organism=="mouse")
     {
       printFlush(" ...MiRNA targets..");
       mirTargets = loadMirnaTargets(file.path(rootDir, "Data/miRNA/Mouse"));
       collections = c(collections, list(MiRNATargets = collectionOfMirnaTargets(mirTargets)));
     } else warning("MiRNA target collection is at present only available for mouse.");
  }

  if (getAgingWGCNA)
  {
    x = loadAsList(file.path(rootDir, "HuntingtonsDisease/OutputResults/200-CreateWGCNAAgingModuleCollection",
          "Collection-RObject/WGCNA.Age.PL.collection.RData"));
    x2 = loadAsList(file.path(rootDir, "HuntingtonsDisease/IndividualAnalyses/070-CHDIAllelicSeries",
                  "CommonAnalysis-040-mRNA-aging/015-StandardAnalysis/RData/AgeAssociationCollection.RData"));
    col1 = mergeCollections(x$WGCNA.Age.PL.collection, x2$assocCollection);
    collections = c(collections, 
         list(AgingWGCNA = convertCollectionToOrganism(col1, organism = organism)));
  }

  if (getEnrichr.Aging)
  {
    enrichrDir = file.path(rootDir, "Data/GeneAnnotation/Enrichr");
    enrichrFiles = file.path(enrichrDir, c("Aging_Perturbations_from_GEO_down", "Aging_Perturbations_from_GEO_up"));
    IDBases = spaste("EnrichrGEOAge", c("Down", "Up"), ".");
    nameBases = spaste(c("Down", "Up"), " with age, ");
    descriptionBases = spaste("Enrichr gene sets: ", c("Down", "Up"), " with age, ");
    enrAgeColl = do.call(mergeCollections,
           mymapply(agingCollectionFromEnrichr, enrichrFile = enrichrFiles,
                    IDBase = IDBases,
                    nameBase = nameBases,
                    descriptionBase = descriptionBases,
                 MoreArgs = list(
                    enrichrClass = "Enrichr age-associated gene sets")));
    
    collections = c(collections, 
          list(Enrichr.Aging = convertCollectionToOrganism(enrAgeColl, organism = organism)))
  }

  collections;
}


rearrangeInternalCollections = function(collections,
   separateCircadianSets = FALSE, circadianPattern = "Circadian oscillations")
{
  select = c("HDSigDB", "HDTargetDB", "internal", "MillerAIBS", "YangLiterature");
  if (any(!select %in% names(collections)))
     warning("The following standard collection names were not found in 'collections':\n",
        paste(select[!select %in% names(collections)], collapse = ", "));

  select = select[select %in% names(collections)];
  if (length(select)==0) stop("No collections to rearrage.");

  pool = do.call(mergeCollections, collections[select]);

  # ChEA/TF binding
  TFTargets = subsetCollection(pool, "Transcription factor targets", 
                    matchComponents = c("groups", "groupAlternateNames"));

  # DE genes for Htt mutation
  DEForHttMutation = subsetCollection(pool, "Effect of Htt CAG length expansion",
                                        matchComponents = c("groups", "groupAlternateNames"));

  # Cell type, tissue/region markers
  CellAndRegionMarkers = subsetCollection(pool, c("Cell type markers", "Spatial markers"),
                                        matchComponents = c("groups", "groupAlternateNames"));

  # PPIs
  PPI = subsetCollection(pool, c("Protein-protein interactions"), 
                                       matchComponents = c("groups", "groupAlternateNames"));

  # WGCNA modules
  allWGCNA = subsetCollection(
                 subsetCollection(pool, c("WGCNA"), matchComponents = "groupsAndAlternates"),
                 tags = "Cell type markers", matchComponents = "groupsAndAlternates", invertSearch = TRUE);

  proteinWGCNA = subsetCollection(allWGCNA, c("Protein WGCNA"), matchComponents = "groupsAndAlternates");
  mRNA.WGCNA = subsetCollection(allWGCNA, c("Protein WGCNA", "Methylation WGCNA"), 
                            matchComponents = "groupsAndAlternates", invertSearch = TRUE);

  # general literature (none of the above)

  rearranged = list(TFTargets = TFTargets, DEForHttMutation = DEForHttMutation, 
                   CellAndRegionMarkers = CellAndRegionMarkers, PPI = PPI,
                   mRNA.WGCNA = mRNA.WGCNA, proteinWGCNA = proteinWGCNA);

  if (separateCircadianSets)
    rearranged = c(rearranged, list(Circadian = subsetCollection(pool, tags = circadianPattern, exactMatch = FALSE)));

  pool2 = do.call(mergeCollections, rearranged);

  rearranged = lapply(c(rearranged, 
             list(generalLiterature = subsetCollection(pool, tags = as.character(dataSetIDs(pool2)),
                                  matchComponents = "ID", invertSearch = TRUE))), 
             dropUnreferencedGroups,  verbose = 0);

  c(rearranged, collections[-match(select, names(collections))]);
}

MSigDBFile = function(version = "5.0")
{
  dir = getwd();
  home = sub("Work/.*", "", dir);
  file = spaste(home, "Work/Data/GeneAnnotation/MSigDB/Version", version, "/msigdb_v", version, ".xml");
  if (!file.exists(file)) stop("Cannot find MSigDB file in its usual location:\n   ", file);
  file;
}
     

#==========================================================================================================
#
# combine groups of samples based on inter-cluster vs. intra-cluster distances
#
#==========================================================================================================

combineGroups = function(x, groups, robust = TRUE, threshold, getDetails = FALSE)
{
  x = as.matrix(x);
  nSamples = nrow(x);
  dst = as.matrix(dist(x));
  diag(dst) = NA;
  groupLevels = sort(unique(groups))
  nGroups = length(groupLevels);
  if (nGroups==1) stop("There must be at least 2 groups.");
  if (nSamples != length(groups)) stop("Length of 'groups' must equal the number of rows in 'x'.");
  if (any(is.na(groups))) stop("All entries of 'groups' must be non-missing.");

  meanDist = mean(dst, na.rm = TRUE);

  groupSampleDist = matrix(NA, nSamples, nGroups)
  for (g in 1:nGroups)
     groupSampleDist[, g] = colMeans(dst[ groups==groupLevels[g], , drop = FALSE], na.rm = TRUE);

  rownames(groupSampleDist) = make.unique(groups);
  colnames(groupSampleDist) = groupLevels;

  meanInGroupDst = mean(groupSampleDist[ cbind(1:nSamples, match(groups, groupLevels))], na.rm = TRUE)

  # This hack should force groups with just 1 element to be merged below. 

  groupSampleDist[is.na(groupSampleDist)] = max(groupSampleDist, na.rm = TRUE);

  groupGroupZ = matrix(NA, nGroups, nGroups)
  for (g1 in 1:(nGroups-1)) for (g2 in (g1+1):nGroups)
  {
    within = c(groupSampleDist[groups==groupLevels[g1], g1], groupSampleDist[groups==groupLevels[g2], g2]);
    out = c(groupSampleDist[groups==groupLevels[g1], g2]);
    mw = mean(within, na.rm = TRUE);
    mo = mean(out, na.rm = TRUE);
    groupGroupZ[g1, g2] = groupGroupZ[g2, g1] = mo/meanInGroupDst;
  }

  originalZ = groupGroupZ;

  n = nGroups;
  while (n > 1 && any(groupGroupZ < threshold, na.rm = TRUE))
  {
    which = which.min(c(groupGroupZ))-1;
    g1 = floor(which/n) + 1;
    g2 = floor(which %% n) + 1;
    
    inG1 = groups==groupLevels[g1];
    inG2 = groups==groupLevels[g2];

    groups[inG2] = groupLevels[g1];
    groupSampleDist[, g1] = colMeans(dst[ groups==groupLevels[g1], , drop = FALSE], na.rm = TRUE);
    for (g3 in c(1:n)[-c(g1, g2)])
    {
      within = c(groupSampleDist[groups==groupLevels[g1], g1], groupSampleDist[groups==groupLevels[g3], g3]);
      out = c(groupSampleDist[groups==groupLevels[g1], g3]);
      mw = mean(within, na.rm = TRUE);
      mo = mean(out, na.rm = TRUE);
      groupGroupZ[g1, g3] = groupGroupZ[g3, g1] = mo/meanInGroupDst;
    }
    groupGroupZ = groupGroupZ[-g2, -g2];
    groupLevels = groupLevels[-g2];
    n = n-1;
  }

  if (getDetails) list(groups = groups, originalZ = originalZ, finalZ = groupGroupZ) else groups;
}

combineGroupsByDistance = function(x, groups, threshold,
                                   minGroupSize, getDetails = FALSE,
                                   getIntermediateGroups = FALSE)
{
  x = as.matrix(x);
  nSamples = nrow(x);
  dst = as.matrix(dist(x));
  diag(dst) = NA;
  groupLevels = sort(unique(groups))
  nGroups = length(groupLevels);
  if (nGroups==1) stop("There must be at least 2 groups.");
  if (nSamples != length(groups)) stop("Length of 'groups' must equal the number of rows in 'x'.");
  if (any(is.na(groups))) stop("All entries of 'groups' must be non-missing.");

  groupGroupDist = matrix(NA, nGroups, nGroups)
  for (g1 in 1:(nGroups-1)) for (g2 in (g1+1):nGroups)
  {
    mo = mean(dst[groups==groupLevels[g1], groups==groupLevels[g2]], na.rm = TRUE);
    groupGroupDist[g1, g2] = groupGroupDist[g2, g1] = mo;
  }

  originalDist = groupGroupDist;
  printFlush(spaste("Using threshold ", threshold));

  if (getIntermediateGroups) imGroups = list(groups);
  n = nGroups;
  while (n > 1 && any(groupGroupDist < threshold, na.rm = TRUE))
  {
    which = which.min(c(groupGroupDist))-1;
    g1 = floor(which/n) + 1;
    g2 = floor(which %% n) + 1;
    
    inG1 = groups==groupLevels[g1];
    inG2 = groups==groupLevels[g2];

    groups[inG2] = groupLevels[g1];
    for (g3 in c(1:n)[-c(g1, g2)])
    {
      mo = mean(dst[groups==groupLevels[g1], groups==groupLevels[g3]])
      groupGroupDist[g1, g3] = groupGroupDist[g3, g1] = mo;
    }
    groupGroupDist = groupGroupDist[-g2, -g2, drop = FALSE];
    groupLevels = groupLevels[-g2];
    n = n-1;
    if (getIntermediateGroups) imGroups = c(imGroups, list(groups));
  }

  groupSizes = as.numeric(table(groups));
  smallGroups = which(groupSizes < minGroupSize);
  largeGroups = which(groupSizes >= minGroupSize);
  nSmall = length(smallGroups);
  if (nSmall > 0)
  {
    printFlush("Merging small groups with larger ones.");
    slGroupDist = groupGroupDist[smallGroups, , drop = FALSE];
    while (nSmall > 0)
    {
      which = which.min(c(slGroupDist))-1;
      lg = floor(which/nSmall) + 1;
      sgi = floor(which %% nSmall) + 1;
      sg = smallGroups[sgi];

      inLG = groups==groupLevels[lg];
      inSG = groups==groupLevels[sg];

      groups[inSG] = groupLevels[lg];
      if (getIntermediateGroups) imGroups = c(imGroups, list(groups));

      groupSizes = as.numeric(table(groups));
      groupLevels = sort(unique(groups));
      nGroups = length(groupLevels);
      smallGroups = which(groupSizes < minGroupSize);
      nSmall = length(smallGroups);
      if (nSmall > 0)
      {
        # Recalculate small to large group distances since we may have one extra large group now.
        slGroupDist = matrix(NA, nSmall, nGroups);
        for (sg in 1:nSmall) for (lg in 1:nGroups) if (smallGroups[sg]!=lg)
          slGroupDist[sg, lg] = mean(dst[groups==groupLevels[smallGroups[sg]],
                                         groups==groupLevels[lg]], na.rm = TRUE)
      }
    }
  }

  if (getDetails) list(groups = groups, originalDist = originalDist, 
        finalDistBeforeMergingSmall = groupGroupDist,
        intermediateGroups = if (getIntermediateGroups) imGroups else NULL
     ) else groups;
}




#=========================================================================================================
#
# Find peak and half-width for expression histograms
#
#=========================================================================================================

histPeakAndHalfWidth = function(h, drop = 1/2)
{
  imax = which.max(h$counts);
  max = h$counts[imax];

  i = imax;
  hv = h$counts[1:imax];
  hv1 = c(0, hv[-imax])
  
  idrop  = which(hv >= max*drop & hv1 < max*drop);
  c(peak = h$mids[imax], width = h$mids[imax] - h$mids[idrop]);
}


#=========================================================================================================
#
# loadAsList
#
#=========================================================================================================

loadAsList = function(file)
{
  env = new.env();
  load(file = file, envir = env);
  as.list(env, all.names = TRUE);
};

      
#=========================================================================================================
#
# Personal: calculator of investment growth
#
#=========================================================================================================

investmentGrowth = function(base, contribution, time, compoundRate = 0.05)
{
  C = base + contribution/compoundRate;
  C*exp(compoundRate*time) - contribution/compoundRate;
}


if (FALSE)
{
  base = 366000
  compoundRate = 1.035^(1/12)-1;
  contribution = -1643
  time = 1:(30*12);
  xx = investmentGrowth(base, contribution, time, compoundRate);
  plot(time/12, xx);
  x2 = investmentGrowth(base, contribution-1000, time, compoundRate);

  lines(time/12, x2)
}



   
#========================================================================================================
#
# histogram mode
#
#========================================================================================================

histMode = function(x, range, breakFrequency = 0.2, plot = FALSE, exclude = NULL, ...)
{
  keep = abs(x) < range;
  keep[is.na(keep)] = FALSE;
  x = x[keep];
  if (!is.null(exclude)) x = x[!x%in% exclude];
  h = hist(x, breaks = floor(length(x) * breakFrequency), plot = plot, ...);
  xx = h$mids;
  y = log(h$counts);
  xxs = xx*xx
  y[!is.finite(y)] = 0;
  fit = lm(y~xx + xxs);

  center = -fit$coefficients[2]/(2*fit$coefficients[3]);
  if (plot)
  {
     lines(xx, exp(predict(fit)));
     points(center, exp(predict(fit, newdata = data.frame(xx = center, xxs = center^2))),
              pch = 21, col = "red", bg = "red");
  }
  center
}


#========================================================================================================
#
# setColnames
#
#========================================================================================================
setColnames = function(x, newColnames)
{
  colnames(x) = newColnames;
  x;
}

setNewNames = function(x, newNames)
{
  names(x) = newNames;
  x;
}


#========================================================================================================
#
# multiHist - plot multiple histograms in one plot
#
#========================================================================================================

plotMultiHist = function(data, nBreaks = 100, col = 1:length(data), scaleBy = c("area", "max", "none"), 
                         logCounts = FALSE,
                         cumulative = FALSE, ...)
{
  if (is.atomic(data)) data = list(data);
  range = range(data, na.rm = TRUE);
  breaks = seq(from = range[1], to = range[2], length.out = nBreaks + 1);
  breaks[nBreaks + 1] = range[2] + 0.001 * (range[2] - range[1]);

  hists = lapply(data, hist, breaks = breaks, plot = FALSE);
  if (logCounts) hists = lapply(hists, function(h) 
   {
     h$counts = log10(h$counts + 1);
     h;
   });

  scaleBy = match.arg(scaleBy);
  invisible(c(plotMultiHist.fromHist(hists, col = col, scaleBy = scaleBy, cumulative = cumulative, ...),
      list(histograms = hists)));
}

plotMultiHist.fromHist = function(hists, col = 1:length(hists), scaleBy = c("area", "max", "none"),
                         cumulative = FALSE, lower.tail = TRUE, ...)
{
  scaleBy = match.arg(scaleBy);
  n = length(hists)

  if (length(cumulative)==1) cumulative = rep(cumulative, n);

  if (any(cumulative))
  {
     hists[cumulative] = lapply(hists[cumulative], function(h) 
     {
          h$counts = cumsum(h$counts)/sum(h$counts);
          if (!lower.tail) h$counts = 1-h$counts;
          h;
     })
     if (any(!cumulative))
        hists[!cumulative] = mymapply(function(h1) 
            {h1$counts = h1$counts/max(h1$counts); h1}, hists[!cumulative]);
  } else {
    if (scaleBy=="max")
    {
       scale = lapply(hists, function(h1) max(h1$counts));
       hists = mapply(function(h1, s1) {h1$counts = h1$counts/s1; h1}, hists, scale, SIMPLIFY = FALSE);
    } else if (scaleBy=="area")
    { 
       scale = lapply(hists, function(h1) sum(h1$counts));
       hists = mapply(function(h1, s1) {h1$counts = h1$counts/s1; h1}, hists, scale, SIMPLIFY = FALSE);
    } 
  }

  multiPlot(x = lapply(hists, getElement, "mids"),
            y = lapply(hists, getElement, "counts"),
            type = "l", col = col, bg = col, ...);
  invisible(list(x = lapply(hists, getElement, "mids"), y = lapply(hists, getElement, "counts")))
}


#========================================================================================================
#
# calculate Wald Z statistics from the output of LRT test. Specifically, convert the p-values to Z scores
# using the coefficient sign.
#
#========================================================================================================

replaceLRTZbyWaldZ = function(deseqResults)
{
  deseqResults[, 4] = qnorm(deseqResults[, 5]/2, lower.tail= FALSE) * sign(deseqResults[, 2])
  deseqResults;
}


#========================================================================================================
#
# Calculate Z from p-value and a sign statistic
#
#========================================================================================================

ZfromP = function(statistics, pValues)
{
 sign(statistics) * qnorm(pValues/2, lower.tail = FALSE);
}


#========================================================================================================
#
# merge technical annotation tables. Basically a form of 'merge' based on the first column, with a check that
# the first columns agree.
#
#========================================================================================================

mergeTechAnnotTables = function(t1, t2)
{
  t2.unique = match(setdiff(colnames(t2), colnames(t1)), colnames(t2));

  t1tot2 = match(t1[, 1], t2[, 1])
  if (!isTRUE(all.equal(t1[, 1], t2[, 1])))
    printFlush("Caution: table rows do not agree.");

  cbind(t1, t2[ t1tot2, t2.unique]);
}


#========================================================================================================
#
# Permutation test on DESeq.
#
#========================================================================================================

permutedDESeqResults = function(counts, colData, design, ..., permuteCol, resultName, nPermutations,
                                sampledSamples = NULL,
                                randomSeed = NULL, verbose = 0, maxContiguousErrors = 10)
{
  if (!is.null(randomSeed))
  {
    savedSeed = .Random.seed;
    on.exit({.Random.seed <<- savedSeed});
    set.seed(randomSeed);
  }
  if (is.character(permuteCol)) permuteCol = match(permuteCol, colnames(colData));
  if (is.na(permuteCol)) 
    stop("'permuteCol' must be finite; perhaps was not found among columns of 'colData'.");

  nGenes = ncol(counts);
  nSamples = nrow(counts);
  permFC = permStat = matrix(NA, nPermutations, nGenes)
  if (is.null(sampledSamples))
    sampledSamples = replicate(n = nPermutations, sample(nSamples),simplify = "array");

  if (verbose > 1) pind = initProgInd();
  p = 1;
  errCount = 0;
  while (p <= nPermutations && errCount <= maxContiguousErrors)
  {
    ss = sampledSamples[, p];
    permCD = colData;
    permCD[, permuteCol] = colData[ss, permuteCol];
    res = try({
      dss = DESeq(DESeqDataSetFromMatrix.PL(t(counts), colData = permCD, design = design), ...);
      results(dss, name = resultName, independentFiltering = FALSE, cooksCutoff = Inf);
    });
    if (inherits(res, "try-error"))
    { 
       errCount = errCount + 1;
    } else {
      permFC[p, ] = res$log2FoldChange;
      permStat[p, ] = res$stat;
      p = p+1;
      errCount = 0;
    } 
    if (verbose > 1) pind = updateProgInd(p/nPermutations, pind);
  }
  if (errCount > maxContiguousErrors)
    stop("Too many failed DESeq calls in a row.")
  if (verbose > 1) printFlush("");
  list(FC = permFC, stat = permStat, sampledSamples = sampledSamples);
}  
    

#========================================================================================================
#
# Extend table by a column with multiple entries per row.
#
#========================================================================================================

duplicateRowsAndMerge = function(tab, addList,
                          nameCol, listColName)
{
  nrow = nrow(tab);
  if (!is.numeric(nameCol)) nameCol = match(nameCol, colnames(tab))
  if (length(nameCol)==0 || is.na(nameCol))
    stop("'nameCol' must be a valid 'tab' column index or name.");
  if (is.null(names(addList)))
  {
    if (length(addList)!=nrow) 
      stop("When 'addList' does not have valid 'names', its length must equal number of rows in 'tab'.");
    names(addList) = tab[, nameCol];
  }

  out = do.call(rbind, lapply(1:nrow, function(r)
  {
    r2 = match(tab[r, nameCol], names(addList));
    n = length(addList[[r2]]);
    if (n==0)
    {
      cbind(tab[r, , drop = FALSE], xyz91827491827491="")
    } else
      cbind(do.call(rbind, listRep(tab[r, , drop = FALSE], n)),
            xyz91827491827491 = addList[[r2]])
  }));
  colnames(out)[ colnames(out)=="xyz91827491827491"] = listColName;
  out;
}

  
#========================================================================================================
#
# Drop group and set ID from an enrichment table
#
#========================================================================================================

dropSetIDandGroups = function(enrTab)
{
  enrTab[, multiGrep(c("^dataSetID$", "^inGroups$"), colnames(enrTab), invert = TRUE)];
}

#========================================================================================================
#
# Split up a raw nanostring excel file.
#
#========================================================================================================

# Nanostring files come in rather clumsy format. There are a few lines of sample information, then the
# expression data.

processNanoStringFile = function(xlsx, sheetIndex = 1,
                                 probeTypes = c("Positive", "Negative", "Endogenous", "Housekeeping"),
                                 organism = "mouse", 
                                 convertGeneNamesToEntrez = TRUE,
                                 convertProbeTypes = c("Endogenous"),
                                 dropUninformative = TRUE)
{
  library("xlsx")
  if (length(sheetIndex)>1) 
     return(lapply(sheetIndex, function(sh) processNanoStringFile(xlsx, sheetIndex = sh)));

  data0 = read.xlsx(xlsx, sheetIndex = sheetIndex, header = FALSE);

  dataRows = which(data0[, 1]%in% probeTypes)
  headerRows = 1:(min(dataRows)-1);
  header0 = data0[headerRows, -c(2:3)];
  header1 = t(header0[, -1]);
  colnames(header1) = header0[, 1];
  header2 = convertNumericColumnsToNumeric(header1);
  if (dropUninformative)
    header2 = dropConstantColumns(header2);

  expr0 = data0[dataRows, -c(1:3)];
  probeType = data0[dataRows, 1];
  geneName = data0[dataRows, 2] 
  transcriptID = data0[dataRows, 3];
  expr1 = lapply(probeTypes, function(pt) 
     setColnames(convertNumericColumnsToNumeric(t(expr0[probeType==pt,  ])), geneName[probeType==pt]) );


  geneAnnot.0 = lapply(probeTypes, function(pt) 
    setColnames(data0[dataRows[probeType==pt], c(1:3)],
                   c("ProbeType", "Symbol.NS", "Transcript")));

  
  names(expr1) = names(geneAnnot.0) = probeTypes;

  if (convertGeneNamesToEntrez) 
  {
    entrez = lapply( expr1[ convertProbeTypes ], function(x)  
                        convert2entrez(symbol = colnames(x), organism = organism));
    keep = lapply(entrez, function(x) !is.na(x));
    entrez.keep = mymapply(function(x, i) x[i], entrez, keep);
    expr1[ convertProbeTypes ] = mymapply( function(x, cn, index) setColnames(x[, index], cn[index]),
         expr1[ convertProbeTypes ], entrez, keep);
    geneAnnot1 = lapply(entrez.keep, geneAnnotationFromEntrez, organism = organism);
    geneAnnot.0[convertProbeTypes] = mymapply(function(x, i) x[i, ], 
           geneAnnot.0[convertProbeTypes], keep);
    geneAnnot = mymapply(cbind, geneAnnot1, geneAnnot.0[convertProbeTypes]);
  }

  list(expr = expr1, sampleAnnot = header2, probeInfo = geneAnnot.0,
       geneAnnot = geneAnnot, originalData = data0);
}
    
   

#========================================================================================================
#
# Create a module membership label based on kME.
#
#========================================================================================================

moduleLabelsFromKME = function(expr, MEs, namePattern = "^ME",
                              corFnc= cor, corOptions = list(use = 'p'),
                              labelsAreNumeric = TRUE,
                              minKMEtoStay = 0, ignoreLabels = c("0", "grey"),
                              unassignedLabel = 0)
{
  labelLevels = sub(namePattern, "", colnames(MEs));
  if (any(labelLevels %in% ignoreLabels))
  {
     MEs = MEs[, !labelLevels %in% ignoreLabels];
     labelLevels = sub(namePattern, "", colnames(MEs));
  }

  corFnc = match.fun(corFnc);
  corOptions = c(corOptions, list(y = expr, x = MEs));
  kme = do.call(corFnc, corOptions);
  
  indexData = minWhichMin(-kme);
  indexData[, "min"] = -indexData[, "min"];
  indexData[, "which"] [ indexData[, "min"] < minKMEtoStay] = NA;
  
  labels = labelLevels[ indexData[, "which"]];
  labels[is.na(labels)] = unassignedLabel;
  if (labelsAreNumeric) as.numeric(labels) else labels;
}

#========================================================================================================
#
# Run FastLMM-EWASher
#
#========================================================================================================
# This doesn't work (yet)

runEWASher.base = function(ewasherDir,
                          tmpDirBase = "EWASherTempDir",
                          meth, pheno, covar, keepTempDir = FALSE)
{
  keepSamples = is.finite(pheno);
  tmpDir = tempfile(tmpDirBase, tmpdir = ".");
  dir.create(tmpDir);
  #if (!keepTempDir) on.exit(unlink(tmpDir, recursive = TRUE, force = TRUE));
  ewDir = file.path("..", ewasherDir);
  setwd(tmpDir);

  nSamples = nrow(meth);
  sampleNames = rownames(meth)
  if (is.null(sampleNames)) {
    sampleNames = spaste("Sample.", 1:nSamples);
    rownames(meth) = sampleNames;
  }
  pheno.df = data.frame(pheno = pheno);
  rownames(pheno.df) = sampleNames;
  rownames(covar) = sampleNames;

  FLE(ewDir, dataForTable(meth[keepSamples, ], transpose = TRUE, IDcolName = "ID"),
                dataForTable(pheno.df[keepSamples, , drop = FALSE], transpose = FALSE, IDcolName = "Sample"),
                data.frame(aa1203984203 = rep(1, sum(keepSamples)), 
                  dataForTable(covar[keepSamples, , drop = FALSE],  transpose = FALSE, IDcolName = "Sample")));
  setwd("..");
  res.ew = read.table(spaste("results/out_ewasher.txt"), sep = "\t", header = TRUE)
  res.lm = read.table(spaste("results/out_linreg.txt"), sep = "\t", header = TRUE)
  summary = readLines(spaste("results/summary.txt"), warn = FALSE);
  setwd("..");
  if (!keepTempDir) on.exit(unlink(tmpDir, recursive = TRUE, force = TRUE));
  list(res.ew = res.ew, res.lm = res.lm, summary = summary);
}

  

runEWASher= function(ewasherDir,
                     meth, pheno, covar,
                     blockSize = 1e5,
                     randomSeed = 12345,
                     tmpDirBase = "EWASherTempDir",
                     keepTempDir = FALSE)
{

  if (!is.null(dim(pheno))) stop("'pheno' must be a single vector at this time.");

  if (!is.null(randomSeed))
  {
    seed = .Random.seed;
    on.exit({ .Random.seed<<-seed });
    set.seed(randomSeed);
  }

  source(file.path(ewasherDir, "fastlmm-ewasher.r"));
  source(file.path(ewasherDir, "utils.r"));

  nVars = ncol(meth);
   
  nBlocks = ceil(nVars/blockSize);
  blocks.0 = allocateJobs(nVars, nBlocks)

  varOrder = sample(nVars);
  blocks = lapply(blocks.0, function(x) varOrder[x]);
  blockRes = lapply(blocks, function(b) 
    runEWASher.base(ewasherDir, meth = meth[, b], pheno = pheno, covar = covar, 
                    tmpDirBase = tmpDirBase, keepTempDir = keepTempDir));

  gc();
  res.ew.all = do.call(rbind, lapply(blockRes, function(r) r$res.ew[, 1:7]));
  res.lm.all = do.call(rbind, lapply(blockRes, function(r) r$res.lm[, 1:7]));

  res.ew.all$FDR = p.adjust(res.ew.all$Pvalue, method = "fdr");
  res.lm.all$FDR = p.adjust(res.lm.all$Pvalue, method = "fdr");
  list(res.ew.all = res.ew.all, res.lm.all = res.lm.all, blockResults = blockRes);
}

#========================================================================================================
#
# Relate variables that could be binary or continuous
#
#========================================================================================================

relateXY = function(x, y = NULL,
                    cpvFnc = "corAndPvalue",
                    cpvOptions = list())
{
  if (is.null(y))
  {
    y = x;
    single = TRUE
  } else
    single = FALSE;

  cpv = do.call(match.fun(cpvFnc), c(list(x=x, y=y), cpvOptions));

  binary.x = apply(x, 2, function(x) length(unique(x[!is.na(x)]))==2);
  binary.y = apply(y, 2, function(x) length(unique(x[!is.na(x)]))==2);

  p.out = cpv$p;

  for (c1 in 1:(ncol(x)-single)) for (c2 in (if (single) c1+1 else 1):ncol(y))
  {
     zz = try({
     if (binary.x[c1] & binary.y[c2])
       p.out[c1, c2] = fisher.test(x[, c1], y[, c2])$p.value;
     if (binary.x[c1] & !binary.y[c2])
       p.out[c1, c2] = kruskal.test(tapply(y[, c2], x[, c1], identity))$p.value;
     if (!binary.x[c1] & binary.y[c2])
       p.out[c1, c2] = kruskal.test(tapply(x[, c1], y[, c2], identity))$p.value;
     })
     if (inherits(zz, "try-error")) p.out[c1, c2] = NA;
     if (single) p.out[c2, c1] = p.out[c1, c2]
  }   

  list(cor = cpv[[1]], p = p.out, p.cor = cpv[[2]]);
}


#========================================================================================================
#
# Methylation M-value
#
#========================================================================================================

methylationMValue = function(beta)
{
  A = 1/(1/beta-1);
  log2(A);
}

#========================================================================================================
#
# Turning lists into data frames by extending vector componets as necessary
#
#========================================================================================================

complete = function(x, n)
{
  c(x, rep("", n-length(x)));
}

irregularList2DataFrame = function(lst, namesFromAttr = NULL, check.names = TRUE)
{
  lens = sapply(lst, length);
  out = as.data.frame(lapply(lst, complete, n = max(lens)));
  if (!is.null(namesFromAttr))
  {
    names(out) = sapply(lst, attr, namesFromAttr);
  } else 
    names(out) = names(lst);
  out;
}

verboseBarplot.2.fromMeans = function (means, errors = NULL, cellCounts = NULL, main = "",
    xlab = "", ylab = "", cex = 1, cex.axis = 1.5, cex.lab = 1.5,
    cex.main = 1.5, color="grey", numberStandardErrors=1,
    two.sided=TRUE, addCellCounts=!is.null(cellCounts), horiz = FALSE, lim = NULL, ...) 
{
    SE = errors;
    err.bp = function(dd, error, two.sided = FALSE, numberStandardErrors, 
        horiz = FALSE) {
        if (!is.numeric(dd)) {
            stop("All arguments must be numeric")
        }
        if (is.vector(dd)) {
            xval = (cumsum(c(0.7, rep(1.2, length(dd) - 1))))
        }
        else {
            if (is.matrix(dd)) {
                xval = cumsum(array(c(1, rep(0, dim(dd)[1] - 
                  1)), dim = c(1, length(dd)))) + 0:(length(dd) - 
                  1) + 0.5
            }
            else {
                stop("First argument must either be a vector or a matrix")
            }
        }
        MW = 0.25 * (max(xval)/length(xval))
        NoStandardErrors = 1
        ERR1 = dd + numberStandardErrors * error
        ERR2 = dd - numberStandardErrors * error
        if (horiz) {
            for (i in 1:length(dd)) {
                segments(dd[i], xval[i], ERR1[i], xval[i])
                segments(ERR1[i], xval[i] - MW, ERR1[i], xval[i] + 
                  MW)
                if (two.sided) {
                  segments(dd[i], xval[i], ERR2[i], xval[i])
                  segments(ERR2[i], xval[i] - MW, ERR2[i], xval[i] + 
                    MW)
                }
            }
        }
        else {
            for (i in 1:length(dd)) {
                segments(xval[i], dd[i], xval[i], ERR1[i])
                segments(xval[i] - MW, ERR1[i], xval[i] + MW, 
                  ERR1[i])
                if (two.sided) {
                  segments(xval[i], dd[i], xval[i], ERR2[i])
                  segments(xval[i] - MW, ERR2[i], xval[i] + MW, 
                    ERR2[i])
                }
            }
        }
    }
    Means1 = means;
    maxSE = max(c(0, as.vector(SE)), na.rm = TRUE);
    if (is.null(lim)) 
    {
      lim = range(Means1,na.rm = TRUE) + c(-maxSE, maxSE) * numberStandardErrors * (numberStandardErrors>0);
      if (lim[1] > 0) lim[1] = 0;
      if (lim[2] <0) lim[2] = 0;
    }
    
    ret = barplot(Means1, main = main, col = color, xlab = xlab, 
        ylab = ylab, cex = cex, cex.axis = cex.axis, cex.lab = cex.lab, 
        cex.main = cex.main, horiz = horiz, ylim = if (horiz) NULL else lim,
        xlim = if (horiz) lim else NULL,
        ...)
    if (addCellCounts) {
       mtext(text=cellCounts,side=if(horiz) 2 else 1,outer=FALSE,at=ret, col="darkgrey",las=2,cex=.8,...)
    } # end of if (addCellCounts)
    abline(h = 0)
    if (numberStandardErrors > 0 && !is.null(SE)) {
        err.bp(as.vector(Means1), as.vector(SE), two.sided = two.sided, 
            numberStandardErrors = numberStandardErrors, horiz = horiz)
    }
    attr(ret, "height") = as.vector(Means1)
    attr(ret, "stdErr") = as.vector(SE)
    invisible(ret)
}

verboseBarplot.2 = function (x, g,  main = "",
    xlab = NA, ylab = NA, cex = 1, cex.axis = 1.5, cex.lab = 1.5,
    cex.main = 1.5, color="grey", numberStandardErrors=1,
    KruskalTest=TRUE,  AnovaTest=FALSE, two.sided=TRUE, 
    addCellCounts=FALSE, horiz = FALSE, ...) 
{
  g.factor = as.factor(g);
  stderr1 = function(x) {
       sqrt(var(x, na.rm = TRUE)/sum(!is.na(x)))
  }
  SE = tapply(x, g.factor, stderr1)
  if (is.na(ylab)) 
       ylab = as.character(match.call(expand.dots = FALSE)$x)
  if (is.na(xlab)) 
       xlab = as.character(match.call(expand.dots = FALSE)$g)
  Means1 = tapply(x, g.factor, mean, na.rm = TRUE)

  if (length(unique(x)) > 2) {
       p1 = signif(kruskal.test(x ~ g.factor)$p.value, 2)
       if (AnovaTest) 
           p1 = signif(anova(lm(x ~ g.factor))$Pr[[1]], 2)
  }
  else {
       p1 = tryCatch(signif(fisher.test(x, g, alternative = "two.sided")$p.value, 
           2), error = function(e) {
           NA
       })
   }
  if (AnovaTest | KruskalTest) 
        main = paste(main, "p =", p1)
  if (addCellCounts) {
       cellCountsF = function(x) {  sum(!is.na(x)) }
       cellCounts=tapply(x, g.factor, cellCountsF)
  } else cellCounts = NULL;

  ret = verboseBarplot.2.fromMeans(means = Means1, errors = SE, cellCounts = cellCounts, main = main, 
    xlab = xlab, ylab = ylab, cex = cex, cex.axis = cex.axis, cex.lab = cex.lab,
    cex.main = cex.main, color=color, numberStandardErrors=numberStandardErrors,
    two.sided=TRUE, 
    addCellCounts=addCellCounts, horiz = horiz, ...)  
  invisible(ret);
}


#=========================================================================================================
#
# Parallel standard deviation
#
#=========================================================================================================

psd = function(...)
{
    pars = list(...)
    nPars = length(pars)
    dn = NULL
    for (p in 1:nPars) {
        if (mode(pars[[p]]) != "numeric") 
            stop(paste("Argument number", p, " is not numeric."))
        if (p == 1) {
            dim = WGCNA:::.dimensions(pars[[p]])
        }
        else {
            if (!isTRUE(all.equal(WGCNA:::.dimensions(pars[[p]]), dim))) 
                stop("Argument dimensions are not consistent.")
        }
        if (prod(dim) == 0) 
            stop(paste("Argument has zero dimension."))
        if (is.null(dn)) 
            dn = dimnames(pars[[p]])
        pars[[p]] = as.numeric(pars[[p]])
    }
    x = as.matrix(as.data.frame(pars))
    if (any(is.na(x))) 
        warning("The input contains missing data that will be removed.")
    q = rowSds(x, na.rm = TRUE)
    if (length(dim) > 1) 
        dim(q) = dim
    if (!is.null(dn)) 
        dimnames(q) = dn
    q
}


#=========================================================================================================
#
# Pick a reduced set of representative variables from module analysis, based on (consensus) KME
#
#=========================================================================================================

moduleRepresentatives = function(multiExpr, labels,
       maxRepresentatives = 30000,
       minModuleRepresentatives = 20,
       consensusQuantile = 0, 
       signed = TRUE,
       useModules = NULL,
       metaAnalysisWeights = NULL,
       corAndPvalueFnc = corAndPvalue, corOptions = list(), corComponent = "cor",
       randomSeed = 12345,
       regularizationShift = 0)
{
  if (!is.null(randomSeed))
  {
    if (exists(".Random.seed"))
    {
      savedSeed = .Random.seed
      on.exit({.Random.seed <<- savedSeed});
    }
    set.seed(randomSeed);
  }

  printFlush("Calculating eigengenes...");
  MEs = multiSetMEs(multiExpr, universalColors = labels, excludeGrey = TRUE);

  printFlush("Calculating consensus KME...");
  conKME.lst = list();
  nProbes = checkSets(multiExpr)$nGenes;
  nKMEBlocks = ceil(nProbes/30000);
  kmeBlocks = allocateJobs(nProbes, nKMEBlocks)
  for (b in 1:nKMEBlocks)
  {
    printFlush(spaste("Working on block ", b))
    conKME.lst[[b]] = consensusKME(
                      mtd.subset(multiExpr, , kmeBlocks[[b]]),
                      labels, MEs, consensusQuantile = consensusQuantile, corAndPvalueFnc = corAndPvalueFnc,
                      signed = signed, useModules = useModules, metaAnalysisWeights = metaAnalysisWeights,
                      corOptions = corOptions, getQvalues = FALSE,
                      corComponent = corComponent ,useRankPvalue = FALSE);

    drop = multiGrep(c("weightedAverage", "^kME", "equalWeights", "\\.DoFWeights", "^p.kME",
                       "meta.p.RootDoFWeights", "^Z.kME"), colnames(conKME.lst[[b]]));
    conKME.lst[[b]] = conKME.lst[[b]][, -drop];
    gc();
  }

  consKME = do.call(rbind, conKME.lst);

  printFlush("Selecting representatives...");

  modLevels = sort(unique(labels))
  modLevels = modLevels[modLevels!=0];

  moduleSizes = table(labels[labels!=0]);
  nModules = length(modLevels);

  if (minModuleRepresentatives > min(moduleSizes)) minModuleRepresentatives = min(moduleSizes);

  modSizes.2 = moduleSizes - minModuleRepresentatives;
  sizeFactors = sqrt(modSizes.2);
  nTake = minModuleRepresentatives + 
           round((  maxRepresentatives - nModules * minModuleRepresentatives ) * sizeFactors/sum(sizeFactors));

  keepVars = numeric(0)
  for (m in modLevels)
  {
    inModule = which(labels==m);
    kme.in = consKME[[spaste("consensus.kME", m)]][inModule];
    prob = (kme.in - min(kme.in))^2 + regularizationShift;
    keep = sample(inModule, min(length(inModule), nTake[m]), prob = prob/sum(prob));
    keepVars = c(keepVars, keep);
  }

  keepVars;
}


#======================================================================================================
#
# Censored cor and p-value
#
#======================================================================================================
# Correlation from Student t
# t = sqrt(n-2) * r/sqrt(1-r^2)
# r^2 (n-2) = t^2( 1-r^2)
# r^2[ n-2 + t^2] = t^2
# r = t/sqrt(t^2 + n -2)
corFromStudentT = function(t, n)
{
  t/sqrt(t^2 + n -2);
}

# For now I assume that censoring is present only in the variable x, and that x is continuous.
# y can be continous or binary (but not ordinal with multiple classes).
# y can have missing values as well; these will be removed.

censoredCorAndPvalue = function(x, y, 
        censorDirection = c("down", "up"),
        corAndPvalueFnc = "corAndPvalue",
        corAndPvalueOptions = list(),
        corComponent = "cor", 
        binaryTest = c("Student", "Kruskal"),
        var.equal = FALSE,
        flatOutput = FALSE,
        yNames = colnames(y),
        getFDR = FALSE,
        nameSep = ".for.")
{
  x = as.matrix(x);
  y = as.matrix(y);

  nx = ncol(x);
  ny = ncol(y);

  if (is.null(yNames))
  {
    if (ny==1) {yNames = ""; nameSep = ""} else {yNames = spaste("y.", 1:ny)};
  } else {
    if (length(yNames)!=ny) stop("Length of 'yNames' must equal the number of columns in 'y'.");
    colnames(y) = yNames;
  }

  censorDirection = match.arg(censorDirection);
  censorSign = if (censorDirection=="down") -1 else 1;

  # Association with non-missing values
  cpvFnc = match.fun(corAndPvalueFnc);
  cpv = do.call(cpvFnc, c(list(x= x, y = y), corAndPvalueOptions));

  # Association with missing values
  
  missing = (is.na(x) + 0);

  nObs.missing = t(missing) %*% (!is.na(y) + 0)

  nLevels = apply(y, 2, function(x) length(unique(x[!is.na(x)])));

  yIsBinary = nLevels == 2;
  yIsContinuous = nLevels > 2;
  yIsValid = nLevels > 1;

  bt = array(NA, dim = c(ny, length(TTestStatNames()), nx));
  if (any(yIsBinary))
  {
    binInd = which(yIsBinary);
    for (ix in 1:nx) for (iy in 1:ny)
    {
      tt = table(x[, ix], y[, binInd[iy]]);
      if (prod(dim(tt))==4)
      {
        ft = fisher.test(as.matrix(tt), conf.int = FALSE);
        bt[binInd[iy], , ix] = c(ft$estimate, 1, ft$p.value, 
                         censorSign * sign(ft$estimate-1) * qnorm(ft$p.value/2, lower.tail = FALSE));
      } else
        bt[binInd[iy], , ix] = c(NA, NA, NA, NA);
    }
  } 
  #anyMissingX = colSums(missing) > 0;
  if (any(yIsContinuous))
  {
    binaryTest = match.arg(binaryTest);
    if (binaryTest=="Student")
    {
       tt = matrixTTest(y[, yIsContinuous, drop = FALSE], missing, reverseGroups = FALSE, 
                                           var.equal = var.equal);
       bt[yIsContinuous, , ] = tt;
    } else 
       bt[yIsContinuous, , ] = matrixKruskalTest(y[, yIsContinuous, drop = FALSE], 1-missing)
  }

  # meta-analysis of associations with non-missing values and with missingness

  transSpecial = function(x)
  {
    dim(x) = c(ny, nx);
    setColnames(t(x), yNames);
  }

  # Do not use the correlation Z since it does not exactly correspond to the p-value.
  # Calculate Z that corresponds exactly to the p-value.
  Z1 = sign(getElement(cpv, corComponent)) * qnorm(cpv$p/2, lower.tail = FALSE);
  Z2 = transSpecial(bt[ , 4, ]);

  Z1Present = 1-is.na(Z1);
  Z2Present = 1-is.na(Z2);

  Z1[is.na(Z1)] = 0;
  Z2[is.na(Z2)] = 0;

  nObs.present = cpv$nObs;
 
  Z = (sqrt(nObs.present) * Z1 * Z1Present + sqrt(nObs.missing) * Z2 * Z2Present) / 
          sqrt( nObs.present * Z1Present + nObs.missing * Z2Present);
  p = 2*pnorm(abs(Z), lower.tail = FALSE); 
  nObs.all = nObs.present + nObs.missing;
  t = sign(Z) * qt(p/2, df = nObs.all-2, lower.tail = FALSE);
  cor.all = corFromStudentT(t, nObs.all);
  cor.all[!is.finite(cor.all)] = NA;

  if (getFDR)
  {
    FDR = apply(p, 2, p.adjust, method = "fdr");
    FDR.present = apply(cpv$p, 2, p.adjust, method = "fdr");
    FDR.missing = apply(transSpecial(bt[,3 ,]), 2, p.adjust, method = "fdr");
  } else
    FDR = FDR.present = FDR.missing = NULL;

  if (flatOutput)
  {
    il = list(cor = cor.all,
              p = p,
              FDR = FDR,
              Z = Z,
              nObs = nObs.all,
              present.cor = getElement(cpv,  corComponent),
              present.p = cpv$p,
              present.FDR = FDR.present,
              present.Z = cpv$Z,
              present.nObs = cpv$nObs,
              missing.statistic = transSpecial(bt[,1 ,]),
              missing.p = transSpecial(bt[,3 ,]),
              missing.FDR = FDR.missing,
              missing.Z = transSpecial(bt[,4 ,]),
              missing.nObs = nObs.missing);
    il = il[ sapply(il, length) > 0];
    interleave(il, baseFirst = TRUE, sep = nameSep);
  } else if (getFDR)
    {
      list(cor = cor.all, p = p, FDR = FDR, Z = Z,
           present = c(cpv, list(FDR = FDR.present)),
           missing = list(statistic = transSpecial(bt[,1 ,]),
                          p = transSpecial(bt[,3 ,]),
                          FDR = FDR.missing,
                          Z = transSpecial(bt[,4 ,]),
                          nObs = nObs.missing));  
    } else
      list(cor = cor.all, p = p, Z = Z,
           present = cpv,
           missing = list(statistic = transSpecial(bt[,1 ,]),
                          p = transSpecial(bt[,3 ,]),
                          Z = transSpecial(bt[,4 ,]),
                          nObs = nObs.missing));  
  
}

#=================================================================================================
#
# metaAnalysis using my censored association function
#
#=================================================================================================

censoredMetaAnalysis = function(multiExpr, multiTrait, 
                        reportWeights = c("Equal", "RootDoF", "DoF"),
                        metaAnalysisWeights = NULL,
                        #alternative = c("two.sided", "less", "greater"),
                        
                        censoredCorAndPvalueFnc = censoredCorAndPvalue,

                        censoredCorAndPvalueOptions = list(
                            censorDirection = c("down", "up"),
                            corAndPvalueFnc = "corAndPvalue",
                            corAndPvalueOptions = list(),
                            corComponent = "cor",
                            binaryTest = c("Student", "Kruskal"),
                            var.equal = TRUE),

                        #getQvalues = FALSE,  ### censoredCorAndPvalue cannot calc'e q values (yet).
                        getFDR = FALSE,
                        useRankPvalue = TRUE,
                        rankPvalueOptions = list(),

                        setNames = NULL, 
                        setNameSep = ".in.",

                        varNames = NULL,
                        additionalVarInfo = NULL,
                        addIDcolumn = TRUE,

                        traitNameSep = ".for.",
                        traitNames = NULL,

                        doCollectGarbage = FALSE)
{

  weightTypes = c("Equal", "RootDoF", "DoF", "User");

  size = checkSets(multiExpr);
  nSets = size$nSets;

  nVars = size$nGenes;

  reportWeights = unique(reportWeights)
  if (!all(tolower(reportWeights) %in% tolower(weightTypes))) 
    stop("Some of the supplied 'reportWeights' are not recognized: ", 
         paste( reportWeights[!tolower(reportWeights) %in% tolower(weightTypes)], collapse = ", "));

  if (!is.null(additionalVarInfo))
  {
    d = dim(additionalVarInfo)
    if (length(d) < 2) stop("If given, 'additionalVarInfo' must be a 2-dimensional array or data frame.");
    if (d[1]!=nVars) 
      stop("If givem, the number of rows in 'additionalVarInfo' must equal the number\n",
           "  of variables (columns) in components of 'multiExpr'.");
  }

  if (is.null(varNames)) varNames = mtd.colnames(multiExpr);
  if (is.null(varNames)) varNames = spaste("Variable.", prependZeros(1:nVars))

  if (length(varNames)!=nVars)
    stop("Length of 'varNames' does not match the number of variables (columns) in 'data'.");

  for (set in 1:nSets)
    multiTrait[[set]] $ data = as.matrix(multiTrait[[set]] $ data);

  #alternative = match.arg(alternative);

  tSize = checkSets(multiTrait);
  if (tSize$nGenes > 1)
  {
    # Call self recursively for each individual trait, put results together and return
    nTraits = tSize$nGenes;
    if (is.null(traitNames)) traitNames = mtd.colnames(multiTrait);
    if (is.null(traitNames)) traitNames = spaste("Trait.", c(1:nTraits));
    traitNames= make.unique(traitNames);
    if (length(traitNames)!=nTraits)
      stop("Length of 'traitNames' must equal the number of traits in 'multiTrait'.");
    for (t in 1:nTraits)
    {
       printFlush(spaste("Working on trait ", traitNames[t], " (", t, " of ", nTraits, ")"));
       mTrait1 = mtd.subset(multiTrait, colIndex = t)
       if (doCollectGarbage) gc();
       ma1 = censoredMetaAnalysis(multiExpr, mTrait1, 
                          reportWeights = reportWeights,
                          metaAnalysisWeights = metaAnalysisWeights,
                          censoredCorAndPvalueFnc = censoredCorAndPvalueFnc,
                          censoredCorAndPvalueOptions = censoredCorAndPvalueOptions,
                          #getQvalues = getQvalues,
                          getFDR = getFDR,
                          useRankPvalue = useRankPvalue,
                          rankPvalueOptions = rankPvalueOptions,
                          setNames= setNames,
                          setNameSep = setNameSep,

                          varNames = varNames,
                          traitNames = traitNames[t],
                          traitNameSep = traitNameSep,
                          addIDcolumn = addIDcolumn && t==1,
                          additionalVarInfo = if (t==1) additionalVarInfo else NULL);
       if (t==1) {
         out = ma1;
       } else 
         out = cbind(out, ma1);
     }
     return(out);
  }
                           
  if (size$nSets!=tSize$nSets)
     stop("The number of sets in 'multiExpr' and 'multiTrait' must be the same.");

  if (!all.equal(size$nSamples, tSize$nSamples))
     stop("Numbers of samples in each set of 'multiExpr' and 'multiTrait' must be the same.");

  #if (!is.finite(consensusQuantile) || consensusQuantile < 0 || consensusQuantile > 1)
  #   stop("'consensusQuantile' must be between 0 and 1.");

  if (is.null(setNames))
     setNames = names(multiExpr);

  if (is.null(setNames))
     setNames = spaste("Set_", c(1:nSets));

  if (!is.null(metaAnalysisWeights))
  {
    if (length(metaAnalysisWeights)!=nSets)
      stop("Length of 'metaAnalysisWeights' must equal the number of sets in 'multiExpr'.")
    if (any (!is.finite(metaAnalysisWeights)) || any(metaAnalysisWeights < 0))
      stop("All weights in 'metaAnalysisWeights' must be positive.");
  }

  setResults = list();

  for (set in 1:size$nSets)
  {
      setResults[[set]] = do.call(censoredCorAndPvalueFnc,
                            c(list(x = multiExpr[[set]]$data, 
                                   y = as.vector(multiTrait[[set]]$data),
                                   getFDR = getFDR,
                                   flatOutput = TRUE,
                                   yNames = traitNames,
                                   nameSep = traitNameSep),
                              censoredCorAndPvalueOptions));
      metaStat = "Z";
  }

  #comb = NULL;
  #for (set in 1:nSets)
  #{
  #  if (set==1) 
  #  {
  #    comb = setResults[[set]];
  #    colNames= colnames(comb);
  #    nColumns = ncol(comb);
  #    colnames(comb) = spaste("X", c(1:nColumns));
  #  } else {
  #    xx = setResults[[set]];
  #    colnames(xx) = spaste("X", c(1:nColumns));
  #    comb = rbind(comb, xx);
  #  }
  #}

  # Re-arrange comb:

  #comb = matrix(as.matrix(as.data.frame(comb)), size$nGenes, nColumns * nSets);
#
#  colnames(comb) = spaste( rep( colNames, rep(nSets, nColumns)), 
 #                          if (!is.null(traitNames)) spaste(traitNameSep, traitNames[1]) else "",
 #                          setNameSep, rep(setNames, nColumns));

  comb = interleave(setResults, nameBase = setNames, sep = setNameSep, baseFirst = FALSE);

  # Find the columns from which to do meta-analysis
  statCols = grep(spaste("^", metaStat), colnames(comb));
  if (length(statCols)==0) stop("Internal error: no columns for meta-analysis found. Sorry!");
  if (length(statCols)!=nSets) 
    stop(spaste("Columns for meta-analysis are not unique. This could be due to an internal error\n",
                "  or because a trait or set name matches one of the names of the meta-analysis \n",
                "  statistics. To facilitate solving this problem, here is the relevant information:\n",
                "  meta statistic search pattern: ^", metaStat, "\n",
                "  Column names of the combined results in which meta statistic is searched for: \n",
                fixLabels(paste( colnames(comb), collapse = ", "), maxCharPerLine = 60, split = " ")));
 
  setZ = comb[, statCols, drop = FALSE];

  colnames(setZ) = spaste("Z", if (!is.null(traitNames)) spaste(traitNameSep, traitNames[1]) else "",
                          setNameSep, setNames);
  nObsCols = grep("^nObs", colnames(comb));
  nObs = comb[, nObsCols, drop = FALSE];

  useWeightIndex = match(tolower(reportWeights), tolower(weightTypes));
  powers = c(0, 0.5, 1, NA)[useWeightIndex];
  nMeta = length(useWeightIndex);

  metaNames = c("equalWeights", "RootDoFWeights", "DoFWeights", "userWeights")[useWeightIndex];
  if (nMeta==1) metaNames = "metaAnalysis";

  metaResults = NULL;
  for (m in 1:nMeta)
  {
    if (useWeightIndex[m] <=3) 
    {
      weights = nObs^powers[m]
    } else
      weights = matrix( metaAnalysisWeights, size$nGenes, nSets, byrow = TRUE);

    metaZ = rowSums( setZ * weights, na.rm = TRUE) / sqrt(rowSums(weights^2, na.rm = TRUE))
    p.meta = 2*pnorm(abs(metaZ), lower.tail = FALSE);
    #if (getQvalues)
    #{
    #  q.meta = qvalue.restricted(p.meta);
    #  meta1 = cbind(metaZ, p.meta, q.meta)
    #  colnames.1 = c("Z.", "p.", "q.");
    #} else {
      q.meta = NULL;
      meta1 = cbind(metaZ, p.meta);
      colnames.1 = c("Z.", "p.");
    #}
    if (getFDR)
    {
      fdr.meta = p.adjust(p.meta, method = "fdr");
      meta1 = cbind(meta1, fdr.meta);
      colnames.1 = c(colnames.1, "FDR.");
    }
    colnames(meta1) = spaste(colnames.1, metaNames[m], 
                             if (!is.null(traitNames)) spaste(traitNameSep, traitNames[1]) else "")
    metaResults = cbind(metaResults, meta1);
  }

  # Use rankPvalue to produce yet another meta-analysis

  rankMetaResults = NULL;
  if (useRankPvalue)
  {
    rankPvalueOptions$datS = as.data.frame(setZ);
    if (is.na(match("calculateQvalue", names(rankPvalueOptions))))
      rankPvalueOptions$calculateQvalue = getFDR# getQvalues;
    for (m in 1:nMeta)
    {
      if (useWeightIndex[m] <=3) {
        weights = nObs^powers[m]
      } else
        weights = matrix( metaAnalysisWeights, size$nGenes, nSets, byrow = TRUE);

      # rankPvalue requires a vector of weights... so compress the weights to a vector.
      # Output a warning if the compression loses information.
      nDifferent = apply(weights, 2, function(x) {length(unique(x)) });
      if (any(nDifferent)>1)
        printFlush(paste("Warning in metaAnalysis: rankPvalue requires compressed weights.\n", 
                         "Some weights may not be entirely accurate."));
      rankPvalueOptions$columnweights = colMeans(weights, na.rm = TRUE);
      rankPvalueOptions$columnweights = rankPvalueOptions$columnweights / sum(rankPvalueOptions$columnweights)
      rp = do.call(rankPvalue, rankPvalueOptions);
      colnames(rp) = spaste(colnames(rp), ".", metaNames[m],
                            if (!is.null(traitNames)) spaste(traitNameSep, traitNames[1]) else "");
      rankMetaResults = cbind(rankMetaResults, as.matrix(rp));
    }
  }

  # Put together the output

  out = list(ID = if (addIDcolumn) varNames else NULL,
             additionalVarInfo,
             metaResults,
             rankMetaResults,
             comb,
             NULL);   # The last NULL is necessary so the line below works even if nothing else is NULL

  out = do.call(data.frame, out[ -which(sapply(out,is.null),arr.ind=TRUE)])

  out;
}

#======================================================================================================
#
# format p-values to a given number of significant digits
#
#======================================================================================================

formatPvalue = function(p, nSignif)
{
  out = as.character(signif(p,nSignif));
  out[nchar(out) > nSignif + 6] = ""
  out;
}


#======================================================================================================
#
# plot points by group plus mean and error bars for each group
#
#======================================================================================================

plotPointsWithMeansAndErrorBars = function(
  x, y,
  spread,
  meanBar.width = 4*strwidth("M"),
  meanBar.lwd = 3,
  errorBar.width = 3*strwidth("M"),
  errorBar.lwd = 2,
  errorBar.height = 1,
  barCol = 1,
  verbose = FALSE,
  ...
)
{
  keep = !is.na(x) & !is.na(y);
  x = x[keep];
  y = y[keep];

  group = as.numeric(factor(x));
  group.x = sort(unique(x));
  groupMeans = tapply(y, group, mean);
  groupSds = tapply(y, group, stdErr);

  set.seed(123);
  n = length(x);
  noise = runif(n, min = -spread/2, max = spread/2);

  if (verbose) verboseScatterplot(x+noise, y, ...) else plot(x + noise, y, ...);

  ng = length(groupMeans);

  for (g in 1:ng)
  {
    lines(c(group.x[g]-meanBar.width/2, group.x[g]+meanBar.width/2),
          rep(groupMeans[g], 2), lwd = meanBar.lwd, col = barCol);
    lines(rep(group.x[g], 2),
          c(groupMeans[g] - errorBar.height * groupSds[g], groupMeans[g] + errorBar.height * groupSds[g]),
          lwd = errorBar.lwd, col = barCol);
    lines(c(group.x[g]-errorBar.width/2, group.x[g]+errorBar.width/2),
          rep(groupMeans[g] - errorBar.height * groupSds[g], 2),
          lwd = errorBar.lwd, col = barCol);
    lines(c(group.x[g]-errorBar.width/2, group.x[g]+errorBar.width/2),
          rep(groupMeans[g] + errorBar.height * groupSds[g], 2),
          lwd = errorBar.lwd, col = barCol);
  }
}

  
  

#=======================================================================================================
#
# removeListNames
#
#=======================================================================================================

removeListNames = function(lst) 
{
  names(lst) = NULL;
  lst;
}

#=======================================================================================================
#
# orthogonalize columns to a reference. Data should be a matrix, ref a single vector.
#
#=======================================================================================================

orthogonalizeToRef = function(data, ref)
{
  data = as.matrix(data);
  ref = as.matrix(ref);
  sp = as.numeric(t(data) %*% ref)/sum(ref*ref)
  data - as.numeric(ref) %o% sp;
}
  

#=======================================================================================================
#
# Get hybridization date from a set of CEL files.
#
#=======================================================================================================


getHybridizationDayFromCELFiles = function(
   files,
   maxSearchLines = 200,
   subtractMinDay = TRUE,
   clusterCut = 0,
   clusterMethod = c("average", "single", "complete"))
{
  nFiles = length(files)
  hybDate = rep("", nFiles)
  hybTime = rep("", nFiles);

  for (f in 1:nFiles)
  {
    con = gzfile(file.path(exprDir.raw, files[f]));
    open(con);
    found = FALSE; count = 0;
    while (!found && count < 200)
    {
      line = readLines(con, n=1);
      if (grepl("DatHeader=", line)) found = TRUE;
      count = count + 1;
    }
    close(con);
    if (found)
    {
      start = regexpr("[0-9][0-9]/[0-9][0-9]/[0-9][0-9] ", line);
      if (start>-1) 
      {
        hybDate[f] = substring(line, start, start+7)
        hybTime[f] = substring(line, start+9, start + 16)
      }
    }
  }
  linTime = linearizeDate(hybDate, sep = "/", subtractMin = subtractMinDay)
  if (clusterCut > 0)
  {
    tx = linTime;
    tx[!is.finite(tx)] = -1e6;
    tree = hclust(dist(tx), method = match.arg(clusterMethod));
    hybBatch = cutree(tree, h = clusterCut)
  } else 
    hybBatch = as.numeric(factor(linTime));

  out = data.frame(HybridizationDay = linTime, Batch = hybBatch);
  rownames(out)= files;
  out;
}

#============================================================================================================
#
# removeNames - remove names from a list
#
#============================================================================================================


removeNames = function(lst)
{
  names(lst) = NULL;
  lst;
}

#============================================================================================================
#
# get requested columns, or columns full of NA's if the requested columns don't exist.
#
#============================================================================================================

getMatchingColumns = function(x, colNames)
{
  out = do.call(cbind, lapply(colNames, function(n)
  {
    if (n %in% colnames(x)) {
      x[, match(n, colnames(x)), drop = FALSE]
    } else {
      out = data.frame(v = rep(NA, nrow(x)));
      setColnames(out, n);
    }
  }))
  if (is.matrix(x)) out = as.matrix(out);
  out;
}

#============================================================================================================
#
# autoLayout
#
#============================================================================================================

autoLayout = function(nPlots, wide = FALSE, maxPerPage = 10)
{
  if (nPlots > maxPerPage) 
  {
    n1 = nPlots/maxPerPage+1;
    n2 = nPlots/n1 + 1
    nPlots = n2;
  }
  if (nPlots < 6)
  {
    if (wide || nPlots < 4)
    {
       return(c(1,nPlots));
    } else
       return(c(2, ceil(nPlots/2)));
  };

  c(2, ceil(nPlots/2));
}


#============================================================================================================
#
# significantOrTop
#
#============================================================================================================

combineMainAndBackupIDs = function(mainIDs, backupIDs = NULL)
{
  if (!is.null(backupIDs))
  {
    replace = replaceMissing(mainIDs=="", TRUE);
    if (any(replace)) mainIDs[replace] = backupIDs[replace];
  }
  mainIDs;
}


significantOrTop = function(IDs, backupIDs = NULL, Z, p=NULL, pThreshold = 0.05, 
                            p.adjusted, pAdjThreshold = 0.05, zThreshold = 0, nTop = 100,
                            extraStat = NULL,
                            extraThreshold = NULL,
                            returnIDs = TRUE, warn = TRUE)
{
  IDs = combineMainAndBackupIDs(IDs, backupIDs);
  if (any(is.na(IDs)))
  {
    Z[is.na(IDs)] = p[is.na(IDs)] = p.adjusted[is.na(IDs)] = NA;
    if (warn) warning("Some elements in 'IDs' are missing. The corresponding Z and p are dropped.");
  }
  order = order(Z, na.last = TRUE);
  keep = !duplicated(IDs[order]);
  if (is.null(p)) { p = p.adjusted; pThreshold = pAdjThreshold + 1; }
  if (nTop > sum(keep)) 
  {
    nTop = sum(keep);
    if (warn) warning("'nTop' is greater than the number of present and non-duplicated 'IDs'.\n",
            " All present and non-duplicated 'IDs' will be returned.");
  }
  if (is.null(extraStat) || is.null(extraThreshold))
  { 
    extraStat = rep(1, length(Z));
    extraThreshold = 0;
  }
  if (sum(p[keep] < pThreshold & p.adjusted[keep] < pAdjThreshold & Z[keep] < zThreshold & extraStat[keep] > extraThreshold, 
          na.rm = TRUE) < nTop)
  {
    index = order[keep] [1:nTop]
  } else
    index = order[keep] [replaceMissing(p.adjusted, pAdjThreshold + 1)[order[keep]] < pAdjThreshold  & 
                         replaceMissing(p, pThreshold + 1)[order[keep]] < pThreshold  & 
                         replaceMissing(Z, zThreshold + 1)[order[keep]] < zThreshold & 
                         replaceMissing(extraStat, extraThreshold-1) [order[keep]] >= extraThreshold ];
  if (returnIDs) IDs[index] else index;
}

#=============================================================================================================
#
# addAblines
#
#=============================================================================================================

addAblines = function(at=-qnorm(0.025), at.x = at, at.y = at, col = "grey", ...)
{
  if (is.finite(at.y))
  {
    abline(h=at.y, col = col, ...);
    abline(h=-at.y, col = col, ...);
  }
  if (is.finite(at.x))
  {
    abline(v=at.x, col = col, ...);
    abline(v=-at.x, col = col, ...);
  }
}

fdrFromZ = function(z, p = NULL, alternative = c("two.sided", "greater", "less"))
{
  if (is.null(p))
  {
    alternative = match.arg(alternative);
    if (alternative == "two.sided") {
      p = 2 * pnorm(abs(z), lower.tail = FALSE);
    } else p = pnorm(z, lower.tail = alternative=="less");
  }
  p.adjust(p, method = "fdr");
}


# Return NA if none of the FDRs are below threshold; return a crudely interpolated value between the lowest
# non-significant and highest significant value.
ablinePositionFromFDR = function(Z, fdr = fdrFromZ(Z,p, alternative), threshold = 0.1, p = NULL, 
                          alternative = c("two.sided", "greater", "less"))
{
  if (is.null(fdr)) fdr = fdrFromZ(Z,p, alternative);
  if (is.null(Z) && !is.null(p)) Z = ZfromP(p, p);
  if (is.null(Z)) Z = ZfromP(fdr, fdr);
  signif = which(replaceMissing(fdr<=threshold))
  nonSignif = which(replaceMissing(fdr > threshold));
  lt = log(threshold);
  if (length(signif)>0)
  {
    i1 = signif [ which.max(fdr[signif])];
    i2 = nonSignif[ which.max(abs(Z[nonSignif]))];
    f1 = log(fdr[i1]);
    f2 = log(fdr[i2]);
    r1 = (lt-f2)/(f1-f2);
    r2 = (lt-f1)/(f2-f1);
    out = r1*abs(Z[i1]) + r2 * abs(Z[i2])
    out;
  } else NA;
}
  
#========================================================================================================
#
# find Z thresholds that correspond to FDR<cut
#
#========================================================================================================

getZThresholdsForFDR = function(Z, fdr, threshold = 0.05)
{
  pos = replaceMissing(Z > 0);
  neg = replaceMissing(Z < 0);
  c(negative = max(Z[neg] [fdr[neg] < threshold], na.rm = TRUE),
    positive = min(Z[pos] [fdr[pos] < threshold], na.rm = TRUE));
}



#=============================================================================================================
#
# grepSingleOrError
#
#=============================================================================================================

grepSingleOrError = function(pattern, x, ..., stopOnMultiple = TRUE)
{
  g = grep(pattern, x, ...)
  if (length(g)==0 || (stopOnMultiple && length(g)>1))
    stop("None or more than one match for ", pattern, " found.");
  g[1];
}


#=============================================================================================================
#
# replace d by greek delta in dN17
#
#=============================================================================================================

fix.dN17 = function(s, pattern = "dN17|deltaN17", replacement = "\U0394N17", ...) 
    gsub(pattern, replacement, s, ...)


#=============================================================================================================
#
# prependDash
#
#=============================================================================================================
prependDash = function(s)
{
  sapply(s, function(s1) if (s1=="") s1 else spaste("-", s1));
}

prependMinus = prependDash;

prependPrefix = function(p, s)
{
  sapply(s, function(s1) if (s1=="") s1 else spaste(p, s1));
}

appendSuffix = function(s, suff)
{
  sapply(s, function(s1) if (s1=="") s1 else spaste(s1, suff));
}

#=============================================================================================================
#
# prettifyNames
#
#=============================================================================================================

prettifyStrings = function(s, prettifyList, from=prettifyList[[1]], to = prettifyList[[2]], fixed = TRUE)
{
  keep = replaceMissing(from!="");
  multiGSub(from[keep], to[keep], s, fixed = fixed);
}

# Multi-step prettify
prettifyStrings.multi = function(s, prettifyLists)
{
  for (pl in 1:length(prettifyLists))
    s = prettifyStrings(s, prettifyLists[[pl]]);
  s;
}

prettifyNames = function(df, prettifyList, from=prettifyList[[1]], to = prettifyList[[2]], fixed = TRUE)
{
  names(df) = multiGSub(from, to, names(df), fixed = fixed);
  df;
}

prettifyColumns = function(df, columns, prettifyList, from=prettifyList[[1]], to = prettifyList[[2]], fixed = TRUE)
{
  if (is.character(columns)) columns = match(columns, names(df));
  columns = as.numeric(columns);
  if (any(is.na(columns))) stop("Some 'columns' are missing or could not be matched to names of 'df'.");
  df[columns] = lapply(df[columns], function(x) multiGSub(from, to, x, fixed = fixed));
  df;
}

removeDuplicatesFromPrettifyList = function(prettifyList)
{
  keep = which(!duplicated(prettifyList[[1]]));
  lapply(prettifyList, `[`, keep);
}

#=============================================================================================================
#
# Optimize existing network layout by switching nodes to minimize a penalty function
#
#=============================================================================================================

layoutPenalty = function(x, y, dstSqMat = NULL, adjacency, signed = FALSE,
     type = c("distance times adjacency", "correlation", "spearman", "squared rank differences"))
{
  if (is.null(dstSqMat))
    dstSqMat = outer(x, x, `-`)^2 + outer(y, y, `-`)^2;

  type = match.arg(type);

  if (!signed) adjacency = abs(adjacency);

  # There are multiple options here: correlation, spearman correlation, sum of rank differences, AUC, ...

  switch(type, 
      `distance times adjacency` = mean(as.numeric(as.dist(dstSqMat)) * as.numeric(as.dist(adjacency))),
      correlation = as.numeric(cor(as.numeric(as.dist(dstSqMat)), as.numeric(as.dist(adjacency)))),
      spearman = cor(as.numeric(as.dist(dstSqMat)), as.numeric(as.dist(adjacency)), method = "spearman"),
      `squared rank differences` = sum( (rank(as.dist(dstSqMat)) - rank(as.dist(-adjacency)))^2 ));
}
 
optimizeLayoutBySwitching = function(x, y, adjacencyMat, penaltyFnc = layoutPenalty,
          penaltyFncOptions = list(signed = FALSE),
          nodeClasses = NULL,
          nIterations = 1000,
          randomSeed = 817,
          temp = 1, gcInterval = 10000,
          verbose = 1, indent = 0)
{
  spaces = indentSpaces(indent);
  n = length(x);
  if (n!=length(y)) stop("Length of 'x' and 'y' must be the same.")
  checkAdjMat(adjacencyMat, min = -Inf, max = Inf);
  if (n!=ncol(adjacencyMat)) stop("Length of 'x' and number of nodes in 'adjacencyMat' must be the same.");

  if (is.null(nodeClasses)) nodeClasses = rep(1, n);
  if (length(nodeClasses)!=n) 
    stop("Length of 'nodeClasses' must equal number of nodes.");

  order = 1:n;

  if (!is.null(randomSeed))
  {
    if (exists(".Random.seed"))
    {
      saved.seed = .Random.seed;
      on.exit(.Random.seed <<- saved.seed);
    }
    set.seed(randomSeed);
  }
  
  dstSqMat = outer(x, x, `-`)^2 + outer(y, y, `-`)^2;

  bestPenalty = do.call(layoutPenalty, c(list(dstSqMat = dstSqMat, adjacency = adjacencyMat),
                                        penaltyFncOptions));
  bestLayout = list(x = x, y = y);
  bestOrder = order;

  penalties = rep(NA, nIterations);
  currentPenalty = bestPenalty;
  startPenalty = bestPenalty;

  classList = lapply(sort(unique(nodeClasses)), function(l) which(nodeClasses==l));

  if (verbose ==1) pind = initProgInd(spaces);
  for (it in 1:nIterations)
  {
    nc = 0;
    while (nc<2) {
      i = floor(runif(1, min = 1, max = n+1));
      nc = length(classList[[nodeClasses[i] ]]);
    }
    j = setdiff(classList[[nodeClasses[i] ]], i)[ floor(runif(1, min = 1, max = nc*(1-1e-13)))]
    if (i==j || nodeClasses[i]!=nodeClasses[j]) stop("Internal error in sampling: sample across classes");
    # t = x[i]; x[i] = x[j]; x[j] = t;
    # t = y[i]; y[i] = y[j]; y[j] = t;
    neworder = order;
    newDst = dstSqMat;
    neworder[i] = order[j]; neworder[j] = order[i];
    t = newDst[, i]; newDst[, i] = newDst[, j]; newDst[, j] = t;
    t = newDst[i, ]; newDst[i, ] = newDst[j, ]; newDst[j, ] = t;

    newPen = do.call(layoutPenalty, c(list(dstSqMat = newDst, adjacency = adjacencyMat),
                                        penaltyFncOptions));

    if (runif(1) < exp(-(newPen - currentPenalty)/temp))
    {
      # accept the proposed move
      order = neworder; 
      dstSqMat = newDst;

      if (newPen < bestPenalty)
      {
         bestLayout = list(x = x[order], y = y[order]);
         bestPenalty = newPen;
         bestOrder = order;
      }
    }
    if ((it %% gcInterval)==0) gc();
    penalties[it] = newPen;
    if (verbose==1) pind = updateProgInd(it/nIterations, pind);
  }

  if (verbose==1) printFlush("");

  c(bestLayout, list(penalties = penalties, startPenalty = startPenalty, bestPenalty = bestPenalty,
     order = bestOrder));
}


if (FALSE)
{
  # test
  n = 20;
  adjMat = matrix(0, n, n);
  for (i in 1:n) for (j in 1:n) adjMat[i,j] = as.numeric(abs(i-j) < 3);
  #adjMat = (max(adjMat) - adjMat);
  #adjMat = adjMat/max(adjMat);

  set.seed(123);
  #x = runif(n);
  #y = runif(n);

  set.seed(18)
  order = sample(n);
  x = cos(order/n * 2 * pi);
  y = sin(order/n * 2 * pi);
  
  nodeClasses = rep( c(1,2), each = n/2);

  layout = optimizeLayoutBySwitching(x = x, y = y, adjacencyMat = adjMat, temp = 0.2, nIterations = 40000,
                                     nodeClasses = nodeClasses, verbose = 4);


  hist(layout$penalties, breaks = 100);

  sizeGrWindow(14, 7);
  par(mfrow = c(1,2));
  networkPlot(adjMat, labels = spaste(1:n),
              pos.x = x, pos.y = y,
              alignAngles = 0, labelAngles = 0, 
              cex.labels = 1, cex.points = 1,
              main = "Random" )


  networkPlot(adjMat, labels = 1:n,
              pos.x = layout$x, pos.y = layout$y,
              alignAngles = 0, labelAngles = 0,
              cex.labels = 1, cex.points = 1,
              main = "Optimized")

}
 
#=========================================================================================================
#
# allocateObjects
#
#=========================================================================================================

# Allocate objects into groups based on pre-defined separators
# here separators are a numeric vector of indices of either first or last elements of each group

allocateObjects = function(nObjects, separators, separatorsAtStart = TRUE)
{
  separators = separators[is.finite(separators)];
  separators = separators[ separators > 0 & separators <= nObjects];

  if (nObjects==0) return(list());
  nSeparators = length(separators);
  if (nSeparators==0) return(list(1:nObjects));

  starts = separators + !separatorsAtStart;
  if (starts[1]>1) starts = c(1, starts);
  starts = starts[starts <=nObjects];

  ends = separators - separatorsAtStart;
  if (ends[nSeparators] < nObjects) ends = c(ends, nObjects);
  ends = ends[ends>0];
  
  mymapply(function(s, e) s:e, starts, ends);
} 

allocateObjects.merged = function(maxLength, nObjects, separators, separatorsAtStart = TRUE,
                                  separatorSpace = 0, tryEven = TRUE)
{
  alloc = allocateObjects(nObjects, separators, separatorsAtStart);
  out = list(list());
  len = 0; index = 1;
  if (separatorSpace >=maxLength) 
    stop("allocateObjects.merged: 'separatorSpace' must be less than 'maxLength'.")
  for (a in 1:length(alloc))
  {
    la = length(alloc[[a]]) + separatorSpace
    if (len==0 || len + la <= maxLength)
    {
      out[[index]] = c(out[[index]], list(alloc[[a]]))
      len = len + la;
    } else {
      index = index + 1;
      out[[index]] = list(alloc[[a]]);
      len = la;
    }
  }

  # Try to even out the distribution if possible.

  if (tryEven)
  {
    last = function(x) x[[length(x)]];
    n = length(out);
    if (n==1) return(lapply(out, unlist));
    trySwitch = rep(TRUE, n-1);
    done = FALSE;
    while (any(trySwitch))
    {
      lens.comp = lapply(out, function(x) sapply(x, length) + separatorSpace);
      lens.total = sapply(lens.comp, sum);
  
      differences = lens.total[-n] - lens.total[-1];
      switch = which(trySwitch) [ which.max(abs(differences[trySwitch]))];
  
      if (abs(differences[switch])<=separatorSpace) { trySwitch[switch] = FALSE; next; }
      largerFirst = differences[switch] > 0;
      larger = if (largerFirst) switch else switch + 1;
      smaller = if (largerFirst) switch+1 else switch;
      compIndex = if (largerFirst) length(out[[larger]]) else 1;
      comp = if (largerFirst) last(out[[larger]]) else out[[larger]] [[1]]
  
      if (lens.comp[[larger]][compIndex] + separatorSpace < abs(differences[switch]))
      {
        if (largerFirst) {
          out[[smaller]] = c(list(comp), out[[smaller]]) 
          out[[larger]] = out[[larger]] [-length(out[[larger]])];
        } else {
          out[[smaller]] =c(out[[smaller]], list(comp));
          out[[larger]] = out[[larger]] [-1];
        }
      } else trySwitch[switch] = FALSE;
    }
  }
  lapply(out, unlist);
}

#=========================================================================================================
#
# dropEmptyElements
#
#=========================================================================================================

dropEmptyElements = function(lst) lst[sapply(lst, length) > 0];


#=========================================================================================================
#
# standardized plot ofdata in groups 
#
#=========================================================================================================

groupBarplot = function(data, group, groupOrder, 
    errors.upper = NULL,
    errors.lower = NULL,
    col = "white", border = 1,
    lines.col = "grey", cex.text = 1,
    yOffsetForXLabels = 1.5,
    xLabel.srt = 0,
    xLabel.adj = c(0.5, 0.5),
    returnData = TRUE,
    errorBarWidth = strwidth("II"),
    errorBarLwd = 1,
    ...)
{
  numGroup = match(group, groupOrder);
  sampleOrder = order(numGroup);
  nSamples = length(data);
  mp = barplot(data[sampleOrder],
          names.arg = rep("", nSamples), #group[sampleOrder], las = 2, 
         col = col, border = border,
         ...)

  if (!is.null(errors.upper))
  {
    if (length(errors.upper)!=length(data)) 
       stop("When 'errors.upper' are given, their length must be the same as length of 'data'.")
    upper = data + errors.upper;
  } else upper = NULL;
  if (!is.null(errors.lower))
  {
    if (length(errors.lower)!=length(data)) 
       stop("When 'errors.lower' are given, their length must be the same as length of 'data'.")
    lower = data + errors.lower;
  } else lower = NULL;

  if (!is.null(errors.upper)) addErrorBars.2sided(x = mp, means = data[sampleOrder], 
      upper = if (is.null(upper)) NULL else upper[sampleOrder],
      lower= if (is.null(lower)) NULL else lower[sampleOrder], 
      width = errorBarWidth,
      lwd = errorBarLwd);

  linePositions = which(group[sampleOrder][-1]!=group[sampleOrder][-nSamples]);

  h = par("usr")[4] - par("usr")[3];
  h2 = par("plt")[3] / (par("plt")[4] - par("plt")[3]) * h;
  ybot = par("usr")[3]
  line.x = rowMeans(cbind(mp[linePositions], mp[linePositions + 1]))
  for (lp in linePositions)
    lines(x = rep(mean(mp[c(lp, lp+1)]), 2), y = c(par("usr")[3] - 2*strheight("M"), par("usr")[4]),
          col = lines.col, xpd = TRUE);

  textPos = rowMeans(cbind(c(par("usr")[1], line.x),
                           c(line.x, par("usr")[2])));

  text(x = textPos, y = rep(ybot - yOffsetForXLabels*strheight("M"), length(groupOrder)),
       labels = groupOrder, adj = xLabel.adj, xpd = TRUE, cex = cex.text, srt = xLabel.srt)

  #plotData = list(bars = data.frame(x = mp, y = data[sampleOrder]),
  #                text = data.frame(
                        

  invisible(list(barplotMP = mp, textPos.x = textPos, 
                 textPos.y = rep(ybot - 1.5*strheight("M"), length(groupOrder))));
}


#=========================================================================================================
#
# Useful calculations on plotting regions; returns various dimensions in user coordinates
#
#=========================================================================================================

plotRegionDimensions = function()
{
  box = par("usr");
  xMin = box[1];
  xMax = box[2];
  yMin = box[3];
  yMax = box[4];

  userPerInch = par("cxy")[1] / par("cin")[1];

  marginWidths = par("mai") * userPerInch;

  leftMarWidth = marginWidths[2];
  rightMarWidth = marginWidths[4];

  plotRegionWidth = par("pin")[1] * userPerInch
  plotRegionHeight = par("pin")[2] * userPerInch;
  # useful for titles...
  plotRegionCenter.x = (xMax + xMin)/2 + (rightMarWidth - leftMarWidth)/2;

  list(userPerInch = userPerInch,
       box = box,
       xMin = xMin,
       xMax = xMax,
       yMin = yMin,
       yMax = yMax,
       marginWidths = marginWidths,
       leftMarginWidth = leftMarWidth,
       rightMarginWidth = rightMarWidth,
       plotRegionWidth = plotRegionWidth,
       plotRegionHeight = plotRegionHeight);
}

  
#========================================================================================================
#
# Enrichment plot functions
#
#========================================================================================================

## This function also works for general collections, not just collections of modules.

matrixEnrichmentPlot.moduleOverlaps = function(
   enrichment,
   analysisName,
   analysisName.pretty = analysisName,

   plotDir = "Plots",
   plotDevice = "CairoPDF",
   plotFileID = "",
   plotFileBase = "enrichmentMatrix",

   width = NULL,
   height = NULL,
   baseWidth = 6.6,
   widthPerCol = 0.3,
   baseHeight = 2.2,
   heightPerRow = 0.20,
   mar = c(6, 16, 2, 0.5),

   collection, 

   maxColsPerPage = 20,
   maxRowsPerPage = 20,

   rowPThreshold = 1,
   colPThreshold = 1,

   keepColModules = NULL,
   keepRowModules = NULL,

   dropNonNumericModules = FALSE,

   rowOrder = NULL,
   colOrder = NULL,

   colPrefix = "M.", colSuffix = "",
   rowPrefix = "M.", rowSuffix = "",

   lim.logp = 20,

   ...
)
{
  counts1 = t(enrichment$countsInDataSet)
  stopifnot(isTRUE(all.equal(dataSetIDs(collection), colnames(counts1))))
  p1 = t(enrichment$pValues);

  changed = TRUE;
  keepRows = rowSums(p1<rowPThreshold) > 0;
  keepCols = colSums(p1 < colPThreshold) > 0;

  if (!any(keepRows)  || !any(keepCols))
  {
    printFlush("No rows or columns remained after restriction. Returning without generating a plot.");
    return(NULL);
  }

  p1 = p1[keepRows, keepCols, drop = FALSE];
  counts1 = counts1[keepRows, keepCols, drop = FALSE];
  rowModules = as.numeric(rownames(p1))
  if (!is.null(keepRowModules))
  {
    keepRows = rowModules %in% keepRowModules;
    if (sum(keepRows)==0)
    {
      printFlush("No rows or columns remained after restriction. Returning without generating a plot.");
      return(NULL);
    }
    p1 = p1[keepRows, , drop = FALSE];
    counts1 = counts1[keepRows, , drop = FALSE];
    rowModules = rowModules[keepRows];
    rowOrder = order(match(rowModules, keepRowModules));
  } else
    rowOrder = 1:nrow(p1);

  setNames0 = dataSetNames(collection)[keepCols];

  setNames = spaste(sub(").+", "", setNames0), ")")
  colModules = as.numeric(multiSub(c("M\\.*", " .+"), c("", ""), setNames));
  if (any(is.na(colModules)))
  {
    setNames = sub(":.+", "", setNames0);
    colModules = as.numeric(multiSub(c("M\\.*", " .+"), c("", ""), setNames));
  }
  if ( any(is.na(colModules)) && !dropNonNumericModules)
    colModules = setNames;

  if (!is.null(keepColModules)) 
  {
    keepCols = replaceMissing(colModules %in% keepColModules)
    if (sum(keepCols)==0)
    {
      printFlush("No rows or columns remained after restriction. Returning without generating a plot.");
      return(NULL);
    }
    p1 = p1[, keepCols, drop = FALSE]; 
    counts1 = counts1[, keepCols, drop = FALSE];
    colModules = colModules[keepCols];
    colOrder = order(match(colModules, keepColModules));
    setNames = setNames[keepCols];
  } else
    colOrder = 1:ncol(p1);

  if (!is.null(plotDevice))
  {
    if (is.null(height)) height = baseHeight + heightPerRow*nrow(p1);
    if (is.null(width)) width = baseWidth + widthPerCol*min(maxColsPerPage, ncol(p1));
    device = match.fun(plotDevice);
    device(file = file.path(plotDir,
                     spaste("enrichmentInOlderWGCNAModules-Matrix", prependDash(analysisName),
                            prependDash(plotFileID), ".pdf")),
           width = width, height = height);
    on.exit(dev.off());
  }

  yLabels = spaste(rowPrefix, rowModules, rowSuffix);
  xLabels = spaste(colPrefix, colModules, colSuffix);

  mat = -log10(p1);
  mat[mat > lim.logp] = lim.logp;

  par(mar = mar)

  labeledHeatmap.multiPage(mat[rowOrder, colOrder, drop = FALSE],
                 xLabels = spaste("ME", labels2colors(colModules))[colOrder],
                 xSymbols = xLabels[colOrder],
                 yLabels = spaste("ME", labels2colors(rowModules))[rowOrder],
                 ySymbols = yLabels[rowOrder],
                 textMatrix = counts1[rowOrder, colOrder, drop = FALSE],
                 colors = blueWhiteRed(100)[50:100],
                 zlim = c(0, lim.logp), setStdMargins = FALSE, 
                 maxRowsPerPage = maxRowsPerPage,
                 maxColsPerPage = maxColsPerPage, ...)
}
 
enrichmentBarplot.standardAnalysis = function(
   analysis,

   analysisName = analysis$analysisName,
   analysisName.pretty = analysis$analysisName,

   prettifyList = NULL,

   collectionName = "WGCNA",
   plotComponent = "core GNV",
   rescueOnly = TRUE,

   printGenes = TRUE,
   entrezColumn = "Entrez",
   symbolColumn = "Symbol",
   rowOrderFnc = rowOrderFnc.downUp,

   useThreshold = NULL,
   ...
)
{
  rowOrderFnc.downUp = function(enrichmentTable) 
      c(multiGrep(c("^Downregulated", "^down.for"), enrichmentTable$class), 
        multiGrep(c("^Upregulated", "^up.for"), enrichmentTable$class),
        multiGrep(c("^Differentially"), enrichmentTable$class),
        multiGrep(spaste(c("Rescue", "Exacerbation", "ntMove"), ".*down.for"), enrichmentTable$class),
        multiGrep(spaste(c("Rescue", "Exacerbation", "ntMove"), ".*up.for"), enrichmentTable$class),
        multiGrep(spaste(c("Enhancement", "Reversal", "ntModification"), ".*down.for"), enrichmentTable$class),
        multiGrep(spaste(c("Enhancement", "Reversal", "ntModification"), ".*up.for"), enrichmentTable$class));

  enrComp = make.names(spaste(analysisName, ".", plotComponent));
  if (!enrComp %in% names(analysis$enrichment)) enrComp = make.names(plotComponent);
  if (!enrComp %in% names(analysis$enrichment)) 
    stop("Invalid 'plotComponent'.");
  if (is.null(useThreshold)) useThreshold = names(analysis$enrichment[[enrComp]]);
  enrTab = do.call(rbind, lapply( analysis$enrichment[[ enrComp ]][useThreshold], 
               function(enr1) enr1[[collectionName]]$enrichmentTable));
  if (rescueOnly) enrTab = enrTab[grepl("Rescue|Exacerbation", enrTab$class), ];

  if (is.null(prettifyList)) prettifyList = analysis$prettifyList.plots else
    prettifyList = mymapply(c, analysis$prettifyList.plots, prettifyList);

  prettifyList = lapply(prettifyList, as.vector)

  enrichmentBarplot.general(
     enrichmentTable = enrTab,
     analysisName = analysisName,
     #analysisName.pretty = analysis$analysisName,

     prettifyList = prettifyList,

     classGenes = if (printGenes)  ## Need to unlist the top level corresponding to top thresholds
                      unlist(removeListNames(analysis$topEntrezForEnrichment[[plotComponent]]), recursive = FALSE) else NULL,
     geneToSymbolTranslation = analysis$combinedDESeqResults[[plotComponent]] [c(entrezColumn, symbolColumn)],

     rowOrderFnc = rowOrderFnc,

     ...);
}

orderClassesByPValue = function(enrichmentTable)
{
  bestP = tapply(enrichmentTable$pValue, factor(enrichmentTable$class, levels = unique(enrichmentTable$class)), 
                 min, na.rm = TRUE);
  classRank = rank(bestP)
  order(classRank[ match(enrichmentTable$class, names(bestP))]);
}
  
enrichmentBarplot.general = function(
   enrichmentTable,
   analysisName,
   #analysisName.pretty,  # Currently not used

   # This needs to be a named list of genes in each class. Names must match the "class"
   # component of  enrichmentTable. The classGenes could contain gene IDs but then
   # geneToSymbolTranslation must provide a translation table.
   # If NULL (or length 0), no genes will be written out.

   classGenes,
   geneToSymbolTranslation = NULL,

   rowOrderFnc = NULL,
   rowOrderArgs = list(),

   plotDir = "Plots",
   plotDevice = "CairoPDF",
   plotFileBase = "enrichmentOfMovedGenes",
   plotFileID = "",

   prettifyList = NULL,
   classPrefix = NULL,
   prettifyDataSetNames = TRUE,

   width = 12,
   height = NULL,
   maxRowsPerPage = 15,
   heightPerRow = 0.5,
   baseHeight = 1.5,

   keepClassOnOnePage = FALSE,
   maxTermsPerClass = Inf,

   mar = c(6, 25, 2, 0.5),
   mgp = c(2, 0.7, 0),

   plotColumn = "Bonferroni",  
   xlab = expression(-log[10](italic(P)[Bonferroni])),
   enrichmentThreshold = 0.01,
   thresholdColumn = "pValue",

   useModules = NULL,  ## Ignored if valid module label cannot be determined
   keepTermPattern = NULL,

   classBesideBars = TRUE,
   splitClassOnAnd = TRUE,
   separatorSpace = if (classBesideBars) 0 else 3,

   maxLines.genes = 2,
   maxLines.geneSets = 2,

   barColor = if (length(classGenes)>0) "#88FFAA" else "#33FF55",
   barBorder = "royalblue",
   barColorFnc = NULL,
   barColorArgs = list(),
   thresholdLineColor = "#FFCCCC",

   dropModuleEnrichmentLabels = FALSE,
 
   lim = NULL,  ## limit for x axis
   axisLimit = 50,
   dropClassPattern = NULL,
   keepClassPattern = NULL,

   barGap = 0.2,
   drawSeparators = TRUE,
   maxRank = NA,

   classLabelFrac = 0.45,
   termLabelFrac = 0.57,

   cex.classLabels = 0.9,
   cex.termLabels = 1,
   cex.genes = 1,

   classLabelShift = 0,

   additionalTables = NULL,
   additional.pch = 21,
   additional.col = "grey",
   additional.bg = c(1:length(additionalTables)),
   additional.pt.cex = 1,
   additional.plotLegend = FALSE,
   additional.legend = NULL,
   additional.legendArgs = list(),

   ...
)
{
  enrTab = enrichmentTable
  # match the class-data set combinations from the main table to additional tables.
  if (!is.null(additionalTables))
  {
    ids = spaste(enrTab$class, "|", if (is.null(enrTab$dataSetID)) "" else enrTab$dataSetID, "|", 
                 enrTab$dataSetName);
    if (any(duplicated(ids))) warning("Found possibly duplicated enrichment results.");
    additionalTables = lapply(additionalTables, function(tab)
    {
      ids1 = spaste(tab$class, "|", if (is.null(tab$dataSetID)) "" else tab$dataSetID, "|",
                 tab$dataSetName);
      out = tab[match(ids, ids1), ];
      #out[c("pValue", "Bonferroni", "FDR")] = as.data.frame(lapply(
      #     out[c("pValue", "Bonferroni", "FDR")], replaceMissing, 1));
      out;
    })
  }

  thresholdCol.i = match(thresholdColumn, names(enrTab));
  if (is.na(thresholdCol.i)) 
    stop("'thresholdColumn' ", thresholdColumn, " not found among columns in enrichment table.");

  keepRows = enrTab[, thresholdCol.i] <= enrichmentThreshold;
  if (length(dropClassPattern) > 0)
    keepRows = keepRows & !multiGrepl(dropClassPattern, enrTab$class);
  if (length(keepClassPattern) > 0)
    keepRows = keepRows & multiGrepl(keepClassPattern, enrTab$class);

  if (length(keepTermPattern) > 0)
    keepRows = keepRows & (multiGrepl(keepTermPattern, enrTab$dataSetName) | 
                            multiGrepl(keepTermPattern, enrTab$shortDataSetName) | 
                            multiGrepl(keepTermPattern, enrTab$inGroups))

  if (is.finite(maxRank)) keepRows = keepRows & enrTab$rank <= maxRank;

  if (!any(keepRows)) 
  {
    #printFlush("Nothing to plot (1).");
    return(NULL);
  }

  tab1 = enrTab[keepRows, , drop = FALSE];

  if (!is.null(additionalTables)) 
    additionalTables = lapply(additionalTables, function(tab) tab[keepRows, , drop = FALSE]);

  setNames = tab1$dataSetName;
  #setNames = sub("\\[WGCNA.+", "", setNames0);
  moduleLabels = suppressWarnings(as.numeric(multiSub(c("^M\\.", "[ :].+"), c("", ""), setNames)));
  setNames[is.na(moduleLabels)] = tab1$shortDataSetName[is.na(moduleLabels)];

  moduleAnalysis = multiSub(c("[^[]*\\[", "\\]"), c("", ""), setNames);
  #if (any(!is.na(moduleLabels))) browser()
  if (dropModuleEnrichmentLabels)
  {
    isModule = is.finite(moduleLabels);
    if (any(isModule)) {
      setNames[isModule] = spaste("M", moduleLabels[isModule], " [", moduleAnalysis[isModule], "]");
    } else
      warning(immediate. = TRUE,
          "'dropModuleEnrichmentLabels' is 'TRUE' but valid module labels could not be found.");
  }

  if (length(useModules)>0)
  {
    keepRows = is.na(moduleLabels) | (replaceMissing( moduleLabels) %in% useModules )
    if (!any(keepRows))
    {
      printFlush("Nothing to plot (2).");
      return(NULL);
    }
    tab1 = tab1[keepRows, , drop = FALSE];
    setNames = setNames[keepRows];
    if (!is.null(additionalTables))
      additionalTables = lapply(additionalTables, function(tab) tab[keepRows, , drop = FALSE]);
    #moduleLabels = moduleLabels[keepRows];
  }

  # Drop enrichment terms in excess of maxTermsPerClass
  rownames(tab1) = spaste("Row.", 1:nrow(tab1));
  classLevels = unique(tab1$class);
  keepRows = unlist(lapply(classLevels, function(cl)
  {
    n1 = sum(tab1$class==cl);
    rownames(tab1[which(tab1$class==cl)[1:min(n1, maxTermsPerClass)], ]);
  }));
  
  setNames = setNames[rownames(tab1)%in%keepRows];
  #moduleLabels = moduleLabels[rownames(tab1)%in%keepRows];
  keep = rownames(tab1)%in%keepRows;
  if (!is.null(additionalTables))
    additionalTables = lapply(additionalTables, function(tab) tab[keep, , drop = FALSE]);
  tab1 = tab1[keep, ];

  if (is.null(tab1$overlapGenes))
    tab1$overlapGenes = rep("", nrow(tab1));
  overlapGenes0 = strsplit(tab1$overlapGenes, split = "|", fixed = TRUE);
  # This should order the genes in the same  order as they appear in original list submitted to
  # enrichmentAnalysis.
  if (length(classGenes) > 0)
  {
    classGenes = classGenes[tab1$class];
    if (!is.null(geneToSymbolTranslation))
       classGenes = translateUsingTable(classGenes, geneToSymbolTranslation)

    overlapGenes = mymapply(intersect,  classGenes, overlapGenes0);
    overlapGenes.ps = sapply(overlapGenes, base::paste, collapse = ", ");
  } else overlapGenes.ps = NULL;

  if (!is.null(rowOrderFnc)) {
     rowOrder = do.call(rowOrderFnc, c(list(enrichmentTable = tab1), rowOrderArgs));
  } else 
     rowOrder = c(1:nrow(tab1));

  if (!is.null(barColorFnc))
    barColor = do.call(match.fun(barColorFnc), c(list(table = tab1), barColorArgs));

  yLabels = setNames[rowOrder];
  if (prettifyDataSetNames)
  {
    if (!is.null(prettifyList)) 
      yLabels = prettifyStrings(yLabels, prettifyList);
    yLabels = gsub("_", " ", yLabels);
  }

  capitalize = function(s) { substr(s,1,1) = toupper(substr(s, 1,1)); s }
  classes = unique(tab1$class[rowOrder]);
  classes.pretty = classes;
  if (!is.null(prettifyList)) 
    classes.pretty = prettifyStrings(classes.pretty, prettifyList);
  
  if (!is.null(classPrefix)) classes.pretty = spaste(classPrefix, classes.pretty);

  classes.pretty = capitalize(classes.pretty);

  if (splitClassOnAnd) classes.pretty = sub(" and ", "\nand ", classes.pretty);

  plotColi = match(plotColumn, names(tab1));
  if (is.na(plotColi)) stop("Plot column ", plotColumn, " not found among columns of 'enrichmentTable'.");
  lp = -log10(tab1[rowOrder, plotColi]);
  if (!is.null(additionalTables)) {
      additionalLP = lapply(additionalTables, function(tab) -log10(tab[rowOrder, plotColi]));
      additionalLP = lapply(additionalLP, function(x) { x[x>400] = 400; x })
  } else additionalLP = NULL;

  lp[lp>400] = 400;

  if (is.null(lim)) lim = (floor( max(c(lp, 1))) / floor(max(c(lp/4, 1))) +1) * floor(max(c(lp/4,1)));
  if (lim > axisLimit)
  {
    lim = axisLimit;
  }
  if (drawSeparators) {
    sep = match(classes, tab1$class[rowOrder])-1;
  } else 
    sep = NULL;


  if (keepClassOnOnePage) 
  {
    rowsPerPage = allocateObjects.merged(maxLength = maxRowsPerPage, nObjects = nrow(tab1),
                     separators = sep, separatorsAtStart = FALSE, 
                     separatorSpace = separatorSpace/maxLines.genes);
  } else
    rowsPerPage = allocateJobs(nrow(tab1), nWorkers = ceil(nrow(tab1)/maxRowsPerPage));

  maxRP = max(sapply(rowsPerPage, function(rows)
    length(rows) + separatorSpace/maxLines.genes * sum((sep + 1) %in% rows)));

  if (!is.null(plotDevice))
  {
    device = match.fun(plotDevice);
    if (is.null(height)) height = baseHeight + heightPerRow * maxRP;
    suppressWarnings(dir.create(plotDir, recursive = TRUE))
    device(file = file.path(plotDir,
                     spaste(plotFileBase, prependDash(analysisName), 
                            prependDash(plotFileID), ".pdf")),
           width = width, height = height);
    on.exit(dev.off());
  }
  par(mar = mar)
  par(mgp = mgp)

  plotAnnot = function(plotCoords, 
                       rows,
                       separatorIndex,
                       separatorPositions,
                       ...,

                       additionalLogPValues,
                       pt.cex, pt.pch, pt.bg, pt.col,

                       tab1,
                       overlapGenes.ps,
                       sepLabels,
                       formatSepLabels,
                       classBesideBars)
  {
    abline(v = -log10(0.05), col = thresholdLineColor, lwd =2)
    lim = plotCoords$box[2]-plotCoords$box[1];
    if (length(classGenes) > 0)
    {
      n1 = tab1$nCommonGenes[rowOrder][rows];
      genes = formatLabels(spaste(n1, ifelse (n1==1, " gene: ", " genes: "), 
                         overlapGenes.ps[rowOrder][rows]), 
                         maxWidth = 0.96 * lim, cex = cex.genes, maxLines = maxLines.genes);
      text(rep(0.02*lim, length(rows)), plotCoords$yMid,
         genes, cex = cex.genes, adj = c(0, 0.5));
    }

    #labs = formatLabels(classes.pretty, maxWidth = 0.45 * plotCoords$leftMargin);
    if (length(separatorIndex) > 0)
    {
      labs = sepLabels[separatorIndex];
      pd = plotRegionDimensions();
      if (classBesideBars)
      {
        if (formatSepLabels) labs = formatLabels(labs, maxWidth = classLabelFrac * plotCoords$leftMargin);
        labLines = sapply(gregexpr("\n", labs), length) + 1 + classLabelShift;
        labs.x = spaste( sapply(labLines, function(n1) paste(rep("\n", as.integer(n1-1)/2), collapse = "")), 
                        labs);
        text(rep(plotCoords$box[1]-0.98*plotCoords$leftMargin, length(separatorIndex)), 
           plotCoords$yMid[separatorPositions+1],
           labs.x, adj = c(0,0.5), cex = cex.classLabels, xpd = TRUE, font = 2);
      } else {
        if (formatSepLabels) labs = formatLabels(labs, maxWidth = 0.94 * pd$leftMarginWidth);
        text(rep(plotCoords$box[1]-0.98*plotCoords$leftMargin, length(separatorIndex)),
           plotCoords$yTop.separators-0.5* strheight("M"),
             labs, adj = c(0,1), cex = cex.classLabels, xpd = TRUE, font = 2);
      }
    }
    nTabs = length(additionalLogPValues);
    if (nTabs > 0)
    {
      y0 = plotCoords$ymid;
      pt.bg = .checkOrExtend(pt.bg, nTabs);
      pt.col = .checkOrExtend(pt.col, nTabs);
      pt.cex = .checkOrExtend(pt.cex, nTabs);
      pt.pch = .checkOrExtend(pt.pch, nTabs);
      #ySpace = (plotCoords$ytop-plotCoords$ybottom)[1]
      for (t in 1:nTabs)
      {
        x1 = additionalLogPValues[[t]][ rows ];
        x1[x1 > axisLimit] = axisLimit;
        points(x1, plotCoords$yMid, pch = pt.pch[t], bg = pt.bg[t], col = pt.col[t], cex = pt.cex[t]);
      }
      if (additional.plotLegend)
      {
        do.call(legend, 
            c(list(x = plotCoords$box[1]-plotCoords$leftMargin + strwidth("M"), 
                   y = plotCoords$box[3] - strheight("M"), xpd = TRUE,
                   legend = additional.legend, pch = pt.pch, pt.bg = pt.bg, pt.cex = pt.cex, col = pt.col),
              additional.legendArgs));
      }
    
    }
  }

  plotAnnotArgs = list(sepLabels = classes.pretty,
                       tab1 = tab1, 
                       overlapGenes.ps = overlapGenes.ps,
                       classBesideBars = classBesideBars,
                       formatSepLabels = TRUE,
                       additionalLogPValues = additionalLP,
                       pt.pch = additional.pch,
                       pt.col = additional.col,
                       pt.bg = additional.bg,
                       pt.cex = additional.pt.cex);
  bp = labeledBarplot3.multiPage(lp,
               border = barBorder,
               col = barColor, 
               yLabels = yLabels,
               separatorPositions = sep,
               sep.col = "grey",
               sep.ext = 1,
               barGap = barGap,
               xlab = xlab,
               lim = c(0, lim),
               formatYLabels = TRUE,
               yLabels.maxFracExtWidth = if (classBesideBars) termLabelFrac else 0.98,
               yLabels.maxLines = maxLines.geneSets,
               rowsPerPage = rowsPerPage,
               afterFnc = plotAnnot,
               afterArgs = plotAnnotArgs,
               separatorSpace = separatorSpace*strheight("M"),
               cex.lab.y = cex.termLabels,
               ...)
}
  
#=========================================================================================================
#
# Rescue/normalization plot
#
#=========================================================================================================

# Plot the DE Z for rescue vs. WT for all genes that are significantly DE between disease genotype and WT.
# Order the plot by module, then by a reference statistic.
rescuePlot.base = function(
   plotStat,
   rescueStatus,
   rescueColors,
   rescueBgColors,
   rescueLevels = sort(unique(rescueStatus)),
   class = NULL,
   referenceStat = NULL,

   plotIndex = NULL,
   excludeIndex = NULL,

   nodeNames,

   statusPlotOrder = NULL,

   classOrder = NULL,
   separator.col = "grey",
   separator.lwd = 2,

   pt.cex = 1,
   pch = 1,

   classPrefix = "M",
   classSuffix = "",
   classLabels.cex = 1.4,
   classLabels.col = "darkgrey",
   classPValue = NULL,  ## Assumed to correspond to classOrder (!!)

   plotLegend = TRUE,
   main = "",
   cex.main = 1.2,

   nConsider = 1000, nLabel = 30,
   cex.labels = 0.8,

   addTrendLines = TRUE,
   trendLine.width = 50,
   trendLine.col = "orange",
   trendLine.lwd = 3,
   trendLine.lty = 1,
   
   returnData = TRUE,
   ...
)
{
  class = replaceMissing(class);

  n0 = length(plotStat);
  if (is.null(class)) class = rep(1, n0);
  if (length(class)!=n0) stop("Length of 'plotStat' and 'class' must be the same.");
  if (is.null(referenceStat)) referenceStat = seq_len(n0);
  if (length(referenceStat)!=n0) stop("Length of 'plotStat' and 'referenceStat' must be the same.");

  if (is.null(classOrder)) classOrder = sort(unique(class));
  if (is.null(plotIndex)) plotIndex = class %in% classOrder;

  if (length(excludeIndex)>0) plotIndex[excludeIndex] = FALSE;

  plotStat = plotStat[plotIndex];
  class = class[plotIndex];
  referenceStat = referenceStat[plotIndex];
  rescueStatus = rescueStatus[plotIndex];
  nodeNames = nodeNames[plotIndex];

  class.num = match(class, classOrder);

  order = order(class.num, -abs(referenceStat))

  nRescueLevels = length(rescueLevels);
  if (is.null(rescueColors)) rescueColors = 1:nRescueLevels;
  rescueColors = checkOrExtend(rescueColors, nRescueLevels, "rescueColors");
  if (is.null(names(rescueColors))) names(rescueColors) = rescueLevels;
  n = length(plotStat);

  if (!is.null(statusPlotOrder))
  {
     plotOrder = order(match(rescueStatus[order], statusPlotOrder));
  } else plotOrder = 1:n;

  plot(1:n, plotStat, type = "n", xaxt = "none", ..., xlim = c(0, n+1), xlab = "");
  separators = which(class[order][-1] != class[order][-length(class)]) + 0.5;
  lapply(separators, function(s) abline(v = s, col = separator.col, lwd = separator.lwd));
  points(c(1:n)[plotOrder], plotStat[order][plotOrder], pch = pch, cex = pt.cex, 
         col = rescueColors[rescueStatus[order]][plotOrder],
         bg = rescueBgColors[rescueStatus[order]][plotOrder]);

  abline(h = 0, col = "grey10");
  if (addTrendLines)
  {
    cuts = c(0, separators, n+1)
    xCut = cut(1:n, breaks = cuts)
    xCut = tapply(1:n, xCut, identity);
    yCut = lapply(xCut, function(i) plotStat[order] [i]);

    nMods = length(xCut);
    trendLine.width = checkOrExtend(trendLine.width, nMods, "trendLine.width");
    ySmooth = mymapply(smoothGauss, x = xCut, y = yCut, xtest = xCut, width = trendLine.width);
    trendLine.lwd = checkOrExtend(trendLine.lwd, nMods, "trendLine.lwd");
    trendLine.lty = checkOrExtend(trendLine.lty, nMods, "trendLine.lty");
    trendLine.col = checkOrExtend(trendLine.col, nMods, "trendLine.col");
    mymapply(lines, xCut, ySmooth, lwd = trendLine.lwd, lty = trendLine.lty,
               col = trendLine.col)
  }

  title(main, adj = 0, cex.main = cex.main)


  x = c(0.5, separators, n+0.5)

  areaMids = (x[-1] + x[-length(x)])/2;

  prd = plotRegionDimensions();
  
  if (!is.null(classPValue)) stars = pValueStars2(classPValue) else stars = rep("", length(classOrder));
  yShift1 = 1.2*strheight("M", cex = classLabels.cex);
  yShift = rep( c(0:2) *yShift1, length(separators)+1)[1:(length(separators)+1)];
  text(areaMids, prd$yMax + yShift + strheight("M", cex = classLabels.cex), 
       spaste(classPrefix, classOrder, classSuffix, stars),
       cex = classLabels.cex, col = classLabels.col, xpd = TRUE);

  if (plotLegend) legendClean("auto", points.x = 1:n, points.y = plotStat[order],
                    legend = names(rescueColors), pch = pch, col = rescueColors, 
                    pt.bg = rescueBgColors, pt.cex = 1.4, nPositions = 20, tryNCol = c(1:4));

  if (!is.null(nodeNames))
    labelExtremePoints2(c(1:n), plotStat[order], nodeNames[order],
         nConsider = nConsider, nLabel = nLabel, cex = cex.labels,
         directions = c("0+", "0-"), ratio.pointToChar = 0.33 * pt.cex);

  if (returnData)
  {
    plotData = data.frame(x = c(1:n)[plotOrder], y = plotStat[order][plotOrder], 
                          Symbol = nodeNames[order][plotOrder], 
                          BorderColor = RColors2hex(rescueColors[rescueStatus[order]][plotOrder]),
                          FillColors = RColors2hex(rescueBgColors[rescueStatus[order]][plotOrder]));
    nSeparators = length(separators);
    verticalLineData = data.frame(x = separators, 
               lineColor = RColors2hex(checkOrExtend(separator.col, nSeparators, "separator.col")));
    moduleLabelData = data.frame(x = areaMids, text = spaste(classPrefix, classOrder, classSuffix),
               textColor = RColors2hex(checkOrExtend(classLabels.col, length(areaMids), "classLabels.col")));
    invisible(list(pointData = plotData, verticalLineData = verticalLineData, 
                   moduleLabelData  = moduleLabelData));
  } else NULL;
}

#========================================================================================================
#
# Filter RPKMs
#
#========================================================================================================
 
filterRPKMs = function(x, minRatio, minCount, minProportion = minCount/nrow(x),
                       minProportion.automatic = 0.8)
{
  propNonZero = colMeans(x>0, na.rm = TRUE);
  x2 = x;
  x2[x==0] = NA;
  minNonZero = colMins(x2, na.rm = TRUE);
  max = colMaxs(x2, na.rm = TRUE);
  ratio = replaceMissing(max/minNonZero);
  keep = (propNonZero >= minProportion & ratio >= minRatio) | (propNonZero >= minProportion.automatic)
  keep;
}


#========================================================================================================
#
# Legend for number plots
#
#========================================================================================================
 
legendForNumberPlots = function(
   nLeft, nRight, yMid,
   ...)
{
  n = length(nLeft);

  minLeft = min(nLeft)
  minRight = min(nRight)

  box = par("usr");
  maxRight = box[2]
  maxLeft = -box[1]

  if (maxLeft-minLeft > maxRight-minRight)
  {
    x = -maxLeft
    just = 0;
    y = yMid[which.min(nLeft)]
  } else {
    x = maxRight;
    just = 1;
    y = yMid[which.min(nRight)];
  }

  legendClean(x = x, y = y, xjust = just, yjust = 0.5, ncol = 2, ...);
}


#========================================================================================================
#
# Estimate cell type fractions
#
#========================================================================================================

estimateCellTypeFractions = function(data, markerList, dataIsLogTransformed = TRUE, logBase = 2,
                                     dataGeneID = colnames(data))
{
  markerList = lapply(markerList, intersect, dataGeneID);
  nMarkers = sapply(markerList, length);
  if (is.null(names(markerList))) stop("Marker list must have names.");
  if (any(nMarkers==0)) 
      stop("The following marker lists are empty after intersecting with 'dataGeneID':\n",
           paste(names(markerList)[nMarkers==0], collapse = ", "));
  if (dataIsLogTransformed)
  {
    linData = logBase^data;
    logData = data;
  } else {
    linData = data;
    logData = logb(data, base = logBase);
  }

  markerData.lin = lapply(markerList, function(cols) linData[, match(cols, dataGeneID), drop = FALSE]);
  markerData.log = lapply(markerList, function(cols) logData[, match(cols, dataGeneID), drop = FALSE]);

  logShift = min(logData, na.rm = TRUE);
  markerWeights = lapply(markerData.log, function(x) colMeans(x, na.rm = TRUE)-logShift);

  markerMeans = mymapply(function(data, weights)
    rowWeightedMeans(data, w = weights, na.rm = TRUE),
    markerData.log, markerWeights);

  grandMeans = sapply(markerMeans, mean, na.rm = TRUE);

  markerLogFC = do.call(cbind, mymapply(`-`, markerMeans, grandMeans));
  markerFC = logBase^markerLogFC;

  A = t(markerFC) %*% markerFC;
  B = colSums(markerFC)

  fractions0 = solve(A, B)
  nSamples = nrow(data);
  nCellTypes = length(markerList);
  fractions = matrix(fractions0, nSamples, nCellTypes, byrow = TRUE) * markerFC;

  list(cellTypeFractions = fractions, fractionBase = fractions0,
       markerFC = markerFC, markerLogFC = markerLogFC, 
       markerMeans = markerMeans, markerWeights = markerWeights);
}
   

#========================================================================================================
#
# R colors to hex
#
#========================================================================================================

color2hex.1 = function(col)
{
  rgb1 = col2rgb(col);
  rgb(rgb1[1], rgb1[2], rgb1[3], maxColorValue = 255);
}

RColors2hex = function(col)
{
  matched = col %in% colors();
  out = col;
  if (any(matched)) out[matched] = sapply(col[matched], color2hex.1);
  out;
}


#========================================================================================================
#
# group scatterplot: scatterplot in which points are grouped and represented by dots with error bars
#
#========================================================================================================

groupScatterplot = function(x, y, group, 
                     meanFnc = "mean", meanOptions = list(na.rm = TRUE),
                     errorWidth = 1, errorFnc = "stdErr", errorOptions = list(),
                     xlim = NULL, ylim = NULL,
                     barWidth = 2*strwidth("o"),
                     barHeight = 2*strheight("o"),
                     bar.col = 1, bar.lwd = 1,
                     legend = NULL,
                     legendOptions = list(),
                     ...)
{
  meanFnc = match.fun(meanFnc);
  errorFnc = match.fun(errorFnc);
  meanOptions = c(list(X = x, INDEX = group, FUN = meanFnc), meanOptions);
  errorOptions = c(list(X = x, INDEX = group, FUN = errorFnc), errorOptions);

  xMean = do.call(tapply, meanOptions);
  xSD = errorWidth * do.call(tapply, errorOptions);

  meanOptions$X = y;
  errorOptions$X = y;
  yMean = do.call(tapply, meanOptions);
  ySD = errorWidth * do.call(tapply, errorOptions);

  if (is.null(xlim)) xlim = range(c(xMean, xMean + xSD, xMean - xSD), na.rm = TRUE);
  if (is.null(ylim)) ylim = range(c(yMean, yMean + ySD, yMean - ySD), na.rm = TRUE);

  plot(xMean, yMean, xlim = xlim, ylim = ylim, ...);

  nPoints = length(xMean);
  bar.col = checkOrExtend(bar.col, nPoints, "bar.col");
  bar.lwd = checkOrExtend(bar.lwd, nPoints, "bar.col");

  for (p in 1:nPoints)
  {
    if (is.finite(xSD[p]))
    {
      segments(xMean[p] - xSD[p], yMean[p], xMean[p] + xSD[p], yMean[p],
               col = bar.col[p], lwd = bar.lwd[p]);
      segments(xMean[p] - xSD[p], yMean[p] - barHeight/2, xMean[p] - xSD[p], yMean[p] + barHeight/2,
               col = bar.col[p], lwd = bar.lwd[p]);
      segments(xMean[p] + xSD[p], yMean[p] - barHeight/2, xMean[p] + xSD[p], yMean[p] + barHeight/2,
               col = bar.col[p], lwd = bar.lwd[p]);
    }
    if (is.finite(ySD[p]))
    {
      segments(xMean[p], yMean[p]-ySD[p], xMean[p], yMean[p]+ySD[p],
               col = bar.col[p], lwd = bar.lwd[p]);
      segments(xMean[p]-barWidth/2, yMean[p]-ySD[p], xMean[p]+barWidth/2, yMean[p]-ySD[p],
               col = bar.col[p], lwd = bar.lwd[p]);
      segments(xMean[p]-barWidth/2, yMean[p]+ySD[p], xMean[p]+barWidth/2, yMean[p]+ySD[p],
               col = bar.col[p], lwd = bar.lwd[p]);
    }
  }

  if (!is.null(legend)) 
    do.call(legendClean, c(list(legend = legend, 
           points.x = c(xMean, xMean, xMean, xMean - xSD, xMean + xSD), 
           points.y = c(yMean, yMean-ySD, yMean+ySD, yMean, yMean)),
           legendOptions));
  invisible(list(xMean = xMean, yMean = yMean, xError = xSD, yError = ySD));
}



#==================================================================================================
#
# Attempt to do quantile normalization when data are probably not missing at random.
#
#==================================================================================================

# idea is to use only those variables that have no missing data to define the quantile, then interpolate
# non-missing values from variables that do have some missing data.

piecewiseLinearPrediction = function(x, y, xtest, min = -Inf, max = Inf)
{
  dup = duplicated(x);
  x = x[!dup];
  y = y[!dup];
  order = order(x);

  x.ord = x[order];
  y.ord = y[order];
  nTrain = length(x);
  intervals = cut(xtest, breaks = c(min, x.ord, max), labels = FALSE)-1;
  intervals[intervals==0] = 1;
  intervals[intervals==nTrain] = nTrain-1;

  x0 = x.ord[intervals];
  x1 = x.ord[intervals + 1];
  y0 = y.ord[intervals];
  y1 = y.ord[intervals + 1];

  prediction = (xtest - x0)/(x1-x0) * y1 + (x1-xtest)/(x1-x0) * y0;
  prediction;
}

normalizeQuantiles.interpolateMissing = function(x, splineFit = TRUE, splineDFFactor = 100)
{
  x = as.matrix(x);
  nMissing = colSums(is.na(x));

  noMissing = nMissing==0;
  nNoMissing = sum(nMissing==0);
  ncols = ncol(x);
  printFlush(spaste("Found ", nNoMissing, " (", round(100*nNoMissing/ncols), 
                    "% of all) variables without missing data."));
  x0 = x[, nMissing == 0];

  x.qn = t(normalize.quantiles(t(x0)));
  out = x;
  out[, noMissing] = x.qn;

  nSamples = nrow(x);

  for (s in 1:nSamples)
  {
    x1 = x0[s, ];
    y1 = x.qn[s, ];
    if (splineFit)
    {
      fit = lm(y1~bs(x1, df = round(nNoMissing/splineDFFactor), degree = 1, intercept = TRUE));
      predicted = predict(fit, newdata = data.frame(x1 = x[s, ]));
    } else 
      predicted = piecewiseLinearPrediction(x1, y1, x[s, ]);
    #Sanity check
    if (FALSE)
    {
      sizeGrWindow(8,8);
      plot(x[s, ], predicted)
      points(x1, y1, pch = 21, cex = 0.7, bg = "green")
    }
    # if (!isTRUE(all.equal(as.numeric(predicted[noMissing]), as.numeric(y1))))
    #   browser();
    out[s, ] = predicted;
  }
  out;
}

#=========================================================================================================
#
# retrieve Ensemble annotation data
#
#=========================================================================================================

retrieveEnsembleData = function(type = c("Gene", "Transcript", "Protein"),
                                organism,
                                version = "current",
                                rootDir = NULL,
                                pathFromRoot = "Data/GeneAnnotation/Ensemble/015-TextTables",
                                verbose = 1, indent = 0)
{
  spaces = indentSpaces(indent);
  type = match.arg(type);
  fileBases = c(Gene = "EnsembleGeneAnnotation-",
               Transcript = "EnsembleTranscriptAnnotation-",
               Protein = "EnsembleProteinAnnotation-");

  fb = fileBases[type];

  if (is.null(rootDir))
  {
    rootDir = file.path(sub("/Work/.*", "", getwd()), "Work");
    if (verbose > 0) printFlush(spaste(spaces, "Using root path: ", rootDir));
  }
  sn = tolower(organismLabels(organism)[2]);
  substring(sn, 1,1) = multiSub(letters, LETTERS, substring(sn, 1,1));
  sn = sub(" ", "_", sn);
  file = file.path(rootDir, pathFromRoot, spaste(fb, sn, "-", version, ".csv.gz"));
  read.csv(gzfile(file));
}

getRootDir = function()
{
  file.path(sub("/Work/.*", "", getwd()), "Work")
}

#=========================================================================================================
#
# retrieve STRING and HIPPIE PPI data
#
#=========================================================================================================

retrieveSTRINGPPI.raw = function(organism,
                             version = "current",
                             rootDir = NULL,
                             pathFromRoot = "Data/InteractionNetworks/STRING-PPI/015-TextTables", ...)
{
  fb = "BioGRID-physical-";

  if (is.null(rootDir))
  {
    rootDir = file.path(sub("/Work/.*", "", getwd()), "Work");
    printFlush(spaste("Using root path: ", rootDir));
  }
  sn = tolower(organismLabels(organism)[1]);
  edges = read.csv(gzfile(file.path(rootDir, pathFromRoot, spaste(fb, sn, "-", version, "-edges.csv.gz"))));
  nodes = read.csv(gzfile(file.path(rootDir, pathFromRoot, spaste(fb, sn, "-", version, "-nodes.csv.gz"))));
  names(edges) = c("index.1", "index.2", "weight");
  list(nodes = nodes, edges = edges);
}


retrieveBiogridPPI.raw= function(organism,
                             version = "current",
                             rootDir = NULL,
                             pathFromRoot = "Data/InteractionNetworks/BioGRID/015-TextTables", ...)
{
  fb = "STRING-PPI-";

  if (is.null(rootDir))
  {
    rootDir = file.path(sub("/Work/.*", "", getwd()), "Work");
    printFlush(spaste("Using root path: ", rootDir));
  }
  sn = tolower(organismLabels(organism)[1]);
  edges = read.csv(gzfile(file.path(rootDir, pathFromRoot, spaste(fb, sn, "-", version, "-edges.csv.gz"))));
  nodes = read.csv(gzfile(file.path(rootDir, pathFromRoot, spaste(fb, sn, "-", version, "-nodes.csv.gz"))));
  names(edges) = c("index.1", "index.2", "weight");
  list(nodes = nodes, edges = edges);
}


retrieveHIPPIE.raw = function(version = "current",
                             rootDir = NULL,
                             pathFromRoot = "Data/InteractionNetworks/HIPPIE/015-TextTables", ...)
{
  if (is.null(rootDir))
  {
    rootDir = file.path(sub("/Work/.*", "", getwd()), "Work");
    printFlush(spaste("Using root path: ", rootDir));
  }
  edges1 = read.csv(gzfile(file.path(rootDir, pathFromRoot, "HIPPIE-basicPPIInfo.csv.gz")));
  entrez = unique(c(edges1$Entrez.1, edges1$Entrez.2));
  n = length(entrez);
  nodes = data.frame(index = 1:n, entrez = entrez);
  edges = data.frame(index.1 = nodes$Index[match(edges1$Entrez.1, entrez)],
                     index.2 = nodes$Index[match(edges1$Entrez.2, entrez)],
                     weight = edges1[, 3]);
  list(nodes = nodes, edges = edges)
}

PPItoMatrix = function(ppiList, namesFrom = "entrez", convertToPlainMatrix = FALSE)
{
  require(Matrix);
  names = getElement(ppiList$nodes, namesFrom);
  if (is.null(names))
    warning(.immediate = TRUE, 
            "Column ", namesFrom, " was not found in 'ppiList$nodes'.\n",
            "  The resulting matrix will have empty dimnames.");
  n = nrow(ppiList$nodes);
  out = sparseMatrix(ppiList$edges$index.1, ppiList$edges$index.2, x = ppiList$edges$weight,
            dims = c(n,n), dimnames = list(names, names));
  out = (out + t(out))/2;

  if (convertToPlainMatrix)
    out = as.matrix(out);

  out;
}

#===================================================================================================
#
# subset or nothing :)
#
#===================================================================================================
subsetOrNothing = function(x, i) if (length(x)==0) x else x[[i]];


#===================================================================================================
#
# Are two files different?
#
#===================================================================================================

filesDiffer = function(f1, f2, useCheckSums = TRUE)
{
  
  require(tools);
  info = file.info(f1, f2, extra_cols = FALSE);
  if (any(is.na(info$size)) || (info$size[1]!=info$size[2])) return(TRUE);

  if (useCheckSums) 
  {
    cs = md5sum(c(f1, f2));
    if (any(is.na(cs)) || (cs[1]!=cs[2])) return(TRUE);
  }

  return(FALSE);
}
  

#===================================================================================================
#
# network neighborhod of a set of genes
#
#===================================================================================================

MTOMfor2nodes = function(adj, i1, i2)
{
  if (i1==i2) stop("The two indices cannot be the same.");

  #printFlush("  ...matrix product...");
  adjSq = adj %*% adj[, c(i1,i2)];
  
  denom = pmin(adjSq[i2, 1] - 2*adj[i1,i2] - adj[, i1] * adj[, i2],
               adjSq[, 1] - 2*adj[, i1] - adj[i1, i2] * adj[, i2],
               adjSq[, 2] - 2*adj[, i2] - adj[i1, i2] * adj[, i1]) + 3;
 
  num = t(adj[, i1] * adj[, i2]) %*% adj + adj[i1, i2]*(1-adj[, i1] - adj[, i2]) - adj[, i1]*adj[, i2] 
           + adj[, i1] + adj[, i2];

  out = as.numeric(num/denom);
  out[c(i1, i2)] = 1;
  out;
}

#===================================================================================================
#
# consensus eigengene network
#
#===================================================================================================

consensusEigengeneNetwork.simple = function(multiMEs, networkOptions, consensusTree)
{
  size = checkSets(multiMEs);
  if (isMultiData(networkOptions))
  {
    if (length(networkOptions)!=size$nSets) 
      stop("If 'networkOptions' are in multiData format, its length must equal the length of 'multiMEs'.");
  } else 
    networkOptions = list2multiData(listRep(networkOptions, size$nSets));

  networkOptions = mtd.apply(networkOptions, function(no)
  {
    no$networkType = "signed";
    no$power = 1;
    no$TOMType = "none";
    no;
  });
  adjMats = mtd.mapply(WGCNA:::.networkCalculation, multiMEs, networkOptions, MoreArgs = list(verbose =0))

  # For eigengene networks: transform back to correlations
  adjMats = mtd.apply(adjMats, function(x) 2*x-1, returnList = TRUE);

  # For now: ignore consensus tree and simply take the mean...

  consensus = pmean.fromList(adjMats);
  colnames(consensus) = rownames(consensus) = mtd.colnames(multiMEs);
  consensus;
}
  
#===================================================================================================
#
# Cross-set summary of DE/association analysis
#
#===================================================================================================
.grepSingle = function(...)
{
  out = grep(...);
  if (length(out)!=1) stop("Must have a single match.")
  out;
}

crossSetDESummary = function(DETables, setNames = names(DETables),
  traitNames,
  statPatterns = spaste("Z.for.",  traitNames),
  pPatterns = spaste("p.for.", traitNames),
  FDRPatterns = spaste("FDR.for.", traitNames),
  pThreshold = 0.05,
  FDRThreshold = 0.10,
  signif = c("less", "greater"),
  dirWords = c("Downreg", "Upreg"),
  getDetails = FALSE,
  idCols = c(1:3))
{
  nSets = length(DETables);
  directions = c(-1, 1);
  signif = match.arg(signif);
  signifDir = if (signif=="less") -1 else 1;

  table = do.call(cbind, removeListNames(mymapply(function(name, spat, ppat, fpat)
  {
     Zs.lst = lapply(DETables, function(x) x[, .grepSingle(spat, colnames(x))]);
     medianZ = pquantile.fromList(Zs.lst, 0.5);
     medianZ.lst = list(medianZ);
     names(medianZ.lst) = spaste("medianZ.for.", name);
     Zs = do.call(cbind, Zs.lst);
     ps = do.call(cbind, lapply(DETables, function(x) x[, .grepSingle(ppat, colnames(x))]));
     FDRs = do.call(cbind, lapply(DETables, function(x) x[, .grepSingle(fpat, colnames(x))]));
     commonTables = do.call(cbind, c(medianZ.lst, removeListNames(mymapply(function(dir, dirWord)
     {
       pSignif = replaceMissing(ps*signifDir > pThreshold*signifDir & Zs*dir > 0);
       FDRSignif = replaceMissing(FDRs*signifDir > FDRThreshold*signifDir & Zs * dir > 0)
       colnames(pSignif) = spaste(dirWord, ".and.pSignif.for.", name, ".in.", setNames);
       colnames(FDRSignif) = spaste(dirWord, ".and.FDRSignif.for.", name, ".in.", setNames);

       nP = rowSums(pSignif)
       nFDR = rowSums(FDRSignif);

       out.lst = list(nSetsWithFDRSignif = nFDR, nSetsWithPSignif = nP);
       
       names(out.lst) = spaste(names(out.lst), ".", dirWord, ".for.", name);
       out.lst = c(if (FDRThreshold < 1) out.lst[1] else NULL,
                  if (getDetails && FDRThreshold < 1) list(FDRSignif) else NULL,
                  if (pThreshold < 1) out.lst[2] else NULL,
                  if (getDetails && pThreshold < 1) list(pSignif) else NULL);
       as.data.frame(do.call(cbind, out.lst[sapply(out.lst, length) > 0]))
     }, directions, dirWords))));
  }, traitNames, statPatterns, pPatterns, FDRPatterns)))

  cbind(DETables[[1]][, idCols], table);
}
      

#========================================================================================================
#
# Transform rescue score for plotting
#
#========================================================================================================

rescueFractionForPlotting = function(fraction, significantRescue, significantExacerbation, scaleTo = 1)
{
  trafo = 1- abs(1-fraction)
  trafo[!is.finite(trafo)] = 0;
  #trafo[trafo<0] = 0;
  trafo[trafo > 2] = 0;
  trafo[trafo < -2] = 0;
  trafo[trafo < -1] = -2 - trafo[trafo < -1];
  trafo[significantRescue==0 & significantExacerbation==0] = 0;
  trafo * scaleTo;
}


#========================================================================================================
#
# Clip clustering tree heights
#
#========================================================================================================

clipTree = function(tree, height = NULL, p = 0.05)
{
  if (length(height)==0) height = quantile(tree$height, p = p);
  tree$height[ tree$height < height] = height;
  tree;
}

#========================================================================================================
#
# Number and a selected top genes
#
#========================================================================================================

numberAndTopGenes = function(geneLists, keep = 15, maxKeep = 17)
{
  sapply(geneLists, function(lst)
  {
    n = length(lst);
    if (n>maxKeep) s = paste(c(lst[1:keep], "..."), collapse = ", ") else s = paste(lst[1:n], collapse = ", ");
    spaste(n, ": ", s);
  });
}

#========================================================================================================
#
# Suitable margin width
#
#========================================================================================================

marginWidth = function(text, cex=1, device = NULL, vfont = NULL, font = NULL, textMargin = 1)
{
  if (is.null(device)) device = "pdf";
  device = match.fun(device);
  file = tempfile();

  device(file = file, width = 10, height = 10);
  par(mar = c(1,1,1,1));
  plot(c(0, 1), type = "n");
  marginLineWidth.in = par("mai")[1];
  maxWidth.in0 = max(strwidth(text, cex = cex, vfont = vfont, font = font, units= "inches")) + 
                strwidth("M", cex = cex, vfont = vfont, font = font, units= "inches") * textMargin;
  marginWidth.lines = ceil(maxWidth.in0/marginLineWidth.in)
  maxWidth.in = marginWidth.lines * marginLineWidth.in;
  dev.off();
  unlink(file);
  list(inch = maxWidth.in, lines = marginWidth.lines);
}

# Margin width in user coordinates
marginWidth.usr = function(side = 2)  # 1: bottom, 2: left, 3:top, 4:right
{
  marginWidth.in = par("mai")[side];
  box = par("usr");
  width = if (side %in% c(2,4)) box[2]-box[1] else box[4] - box[3];
  dir = if (side %in% c(2,4)) 1 else 2
  plotDim.in = par("pin");
  marginWidth.in / plotDim.in[dir] * width;
}
  

#=========================================================================================================
#
# Approximate variance stabilization
#
#=========================================================================================================

varianceSlope = function(lx, colWeights)
{
  vars = colVars(lx, na.rm = TRUE);
  means = colMeans(lx, na.rm = TRUE);
  out = c(WGCNA::cor(means, vars, weights.x = colWeights, weights.y = colWeights, use = 'p'));
  printFlush(spaste("Means-variance correlation: ", round(out, 2)));
  out;
}

approximateVST = function(x, min = NULL, max = 500, 
      xIsLogTransformed = FALSE, xLogBase = 2, xLogShift = 1, epsilon = 0.3,
      xWeightScale = 0)
{
  if (xIsLogTransformed) x = xLogBase^x-xLogShift;
  if (is.null(min)) min = min(x, na.rm = TRUE) + 1;

  if (xWeightScale > 0) printFlush("Using weight scale coeff ", xWeightScale);

  m = colMeans(x, na.rm = TRUE);
  colWeights = m/(m + xWeightScale);

  s1 = varianceSlope(log2(x+min), colWeights = colWeights);
  s2 = varianceSlope(log2(x + max), colWeights = colWeights);

  while(s2 < 0 && max < 1e6)
  {
    min = max;
    max = 2*max;
    s1 = s2;
    s2 = varianceSlope(log2(x + max), colWeights = colWeights);
  }

  while (max-min > epsilon)
  {
    if (s1*s2 > 0)
    {
      i = which.min(c(abs(s1), abs(s2)));
      printFlush("For your information: could not find a point of constant variance.");
      shift = c(min, max)[i];
      out = log2(x + shift);
      attr(out, "logShift") = shift;
      return(out);
    }  
    middle = (min + max)/2;
    cat(spaste("Trying ", middle, ": "));
    sm = varianceSlope(log2(x + middle), colWeights = colWeights);
    if (s1*sm > 0)
    {
      s1= sm;
      min = middle;
    } else {
      s2 = sm;
      max = middle;
    }
  }
  shift = min;
  out = log2(x + shift);
  attr(out, "logShift") = shift;
  out; 
}
    
#===========================================================================================================
#
# Experimental: enrichment of an overlap
#
#===========================================================================================================

enrichmentOfOverlaps = function(overlappingSets, 
                                finalClassNames = names(overlappingSets),
                                keepIndividualEnrichment = TRUE,
                                getBonferroniCorrection = TRUE,
                                getFDR = TRUE, ...)
{
  if (is.atomic(overlappingSets[[1]])) 
    overlappingSets = list(overlappingSets);
    
  active = lapply(overlappingSets, multiIntersect);
  names(active) = finalClassNames;
  inactive1 = lapply(overlappingSets, `[[`, 1);
  inactive2 = lapply(overlappingSets, `[[`, 2);
  inactive.lst = list(inactive1, inactive2);
  enr1 = enrichmentAnalysis(active = active, inactive = inactive1, combineEnrichmentTables = TRUE, 
                   getBonferroniCorrection = getBonferroniCorrection, getFDR = getFDR, ...);
  enr2 = enrichmentAnalysis(active = active, inactive = inactive2, combineEnrichmentTables = TRUE, 
                   getBonferroniCorrection = getBonferroniCorrection, getFDR = getFDR, ...);

  p.common = cbind(as.numeric(enr1$pValues), as.numeric(enr2$pValues));
  maxIndex = minWhichMin(-p.common, byRow = TRUE)$which;
  indexMat = cbind(1:nrow(p.common), maxIndex);
  p.max = p.common[ indexMat];
  n.max = cbind(as.numeric(enr1$countsInDataSet), as.numeric(enr2$countsInDataSet))[indexMat];

  dim(p.max) = dim(n.max) = dim(enr1$pValues);
  dimnames(p.max) = dimnames(n.max) = dimnames(enr1$pValues);

  if (getBonferroniCorrection) {
    Bonferroni = p.max * length(p.max);
    Bonferroni[Bonferroni > 1] = 1;
  }
  if (getFDR) 
  {
    fdr = p.adjust(as.numeric(p.max), method = "fdr");
    dim(fdr) = dim(p.max);
    dimnames(fdr) = dimnames(p.max);
  }

  # Intersect enrichment tables

  enr.lst = list(enr1, enr2);
  rowIDs = lapply(enr.lst, function(enr1) with(enr1$enrichmentTable, spaste(class, dataSetID)));

  keepRows = multiIntersect(rowIDs);
  nKeep = length(keepRows);
  if (nKeep > 0)
  {
    et.keep = mymapply(function(enr1, rowID1) enr1$enrichmentTable[match(keepRows, rowID1), ],
                enr.lst, rowIDs);

    IDcols = c("class", "dataSetID", "dataSetName", "inGroups", "shortDataSetName", "overlapGenes",
               "nCommonGenes");
    et.dup = mymapply(function(et1, set) 
     {
       out = et1[ , !colnames(et1)%in% IDcols];
       out = setColnames(out, spaste(colnames(out), ".in.background.", set));
       activeIndex = match(et1$class, names(active))
       out1 = data.frame(Background = names(inactive.lst[[set]])[activeIndex]);
       names(out1) = spaste("Background.", set);
       cbind(out1, et1);
     }, et.keep, 1:2);

    etPValues = do.call(cbind, lapply(et.keep, getElement, "pValue"));

    maxIndex.et = minWhichMin(-etPValues, byRow = TRUE)$which;
    indexMat = cbind(1:nKeep, maxIndex.et);
    eTab.lst = lapply(names(enr1$enrichmentTable), function(col)
      do.call(cbind, lapply(et.keep, getElement, col))[indexMat]);
    names(eTab.lst) = names(enr1$enrichmentTable);
    eTab = as.data.frame(eTab.lst);

    if (getBonferroniCorrection || getFDR)
    {
      rows = match(eTab$dataSetID, rownames(enr1$pValues));
      cols = match(eTab$class, colnames(enr1$pValues));
      for (row in 1:nKeep)
      {
        if (getBonferroniCorrection) eTab$Bonferroni[row] = Bonferroni[ rows[row], cols[row]];
        if (getFDR) eTab$FDR[row] = fdr[ rows[row], cols[row]];
      }
    }
    eTab.ext =  cbind(eTab, do.call(cbind, et.dup));
    classOrder = match(eTab.ext$class, unique(eTab.ext$class));
    order = order(classOrder, eTab.ext$pValue);
    eTab.ext = eTab.ext[order, ];
    eTab.ext$rank = unlist(lapply(unique(eTab.ext$class), function(cl) 1:sum(eTab.ext$class==cl)));
  } else {
    eTab.ext = enr1$enrichmentTable[numeric(0), ];
  }

  out = list(enrichmentTable = eTab.ext,
       pValues = p.max,
       Bonferroni = if (getBonferroniCorrection) Bonferroni else NULL,
       FDR = if (getFDR) fdr else NULL,
       countsInDataSet = n.max);

  if (keepIndividualEnrichment)
    out = c(out, list(enrichmentInSet1 = enr1, enrichmentInSet2 = enr2));

  out;
}


#=====================================================================================================
#
# convenience
#
#=====================================================================================================

write.csv.nr = function(...) write.csv(..., row.names = FALSE);
  
# Find a vector that best represents variation due to well letter
# For now this function is GNV-specific
#
#=====================================================================================================


partialScaleTargetToReference = function(x.ref, group.ref, x.target, group.target)
{
  commonGroups = intersect(group.ref, group.target);
  refx1 = list();
  targetx1 = list();
  for (g in commonGroups)
  {
    ref1 = x.ref[ g==group.ref];
    tg1 = x.target[g==group.target];
    #nr1 = length(ref1);
    nt1 = length(tg1);
    mean.ref = mean(ref1);
    #mean.tg = mean(tg1);

    refx1[[g]] = rep(mean.ref, nt1);
    targetx1[[g]] = tg1;
  }

  refx1 = unlist(refx1);
  targetx1 = unlist(targetx1);

  fit = lm(targetx1~refx1);
  fit$coefficients;
}


  
GNV.wellLetterVariable = function(data, wellLetter, genotype, nPCs=4)
{
  svd1 = svd(scale(data), nu = nPCs);
  PCs = as.data.frame(svd1$u);
  wellNumber = as.numeric(factor(wellLetter));
  fit = lm(wellNumber~., data = PCs);
  predicted = predict(fit);
  genoLevels = unique(genotype)
  n = length(wellLetter);
  predictedByGeno = tapply(predicted, genotype, identity);
  letterByGeno = tapply(wellLetter, genotype, identity);
  sampleByGeno = tapply(1:n, genotype, identity);

  nGeno = length(genoLevels);
  ref = which.max(sapply(letterByGeno, function(x) length(unique(x))));

  predictedByGeno.scaled = predictedByGeno;
  for (g in c(1:nGeno)[-ref])
  {
    coeffs = partialScaleTargetToReference(predictedByGeno[[ref]], letterByGeno[[ref]],
                                predictedByGeno[[g]], letterByGeno[[g]]);
    predictedByGeno.scaled[[g]] = predictedByGeno.scaled[[g]] * coeffs[2] + coeffs[1];
  };

  if (FALSE)
  {
    par(mfrow = c(2,2));
    mymapply(plot, predictedByGeno, predictedByGeno.scaled)
  }

  #unlist(predictedByGeno.scaled)[ match(c(1:n), unlist(sampleByGeno)) ];
  meansByGroup = tapply(unlist(predictedByGeno.scaled), unlist(letterByGeno), mean);
  meansByGroup[as.character(wellLetter)];
}


shortenEnrichmentLabels = function(eLabels, 
  patterns = list(character(0), character(0)),
  basePatterns = list(c("\\([^);]*; *", " from PMID [0-9]*"), c("(", "")))
{
  patterns = mymapply(c, patterns, basePatterns);
  if (is.data.frame(eLabels)) lab = eLabels$enrichmentLabel else lab = eLabels;
  lab = multiGSub(patterns[[1]], patterns[[2]], lab);
  if (is.data.frame(eLabels)) eLabels$enrichmentLabel = lab else eLabels = lab;
  eLabels;
}



#=====================================================================================================
#
# Mask a vector by reference direction
#
#=====================================================================================================

# more precisely, this function sets components of x where ref goes in the wrong direction to NA.

maskVectorByReferenceDirection = function(
  x,
  ref,
  name.x = "",
  name.ref = "ref",
  sep.x = if (name.x=="") "" else ".and.",
  directionWords = c("down", "up"),
  wordSep = ".for.",
  includeAll = TRUE,
  mask = NA,
  returnDataFrame = FALSE)
{
  out = list();
  if (includeAll) { out = list(x); names(out) = name.x };
  sign = c(-1, 1);
  out1 = lapply(1:2, function(dir) { x[ref*sign[dir] < 0] = mask; x; });
  names(out1) = spaste(name.x, sep.x, directionWords, wordSep, name.ref);
  if (returnDataFrame) as.data.frame(c(out, out1)) else c(out, out1);
}
  

prepChar = function (s, char) 
{
    ifelse(s == "", s, paste(char, s))
}


#=====================================================================================================
#
# Isolate candidate gene symbols from a HTML file.
#
#=====================================================================================================

getGenesFromHTML = function(file, organisms, 
   opening = "<em>",
   closing = "</em>",
   wordSep = "[^0-9A-Za-z]+", ...)
{
  getGeneSymbolsFromFile(file(file, open = "rt"), organisms = organisms, opening = opening, closing = closing,
       wordSep = wordSep, ...);
}
  

getGeneSymbolsFromFile = function(con, organisms, opening = NULL, closing = NULL, 
    wordSep = "[^0-9A-Za-z]+", 
    excludeWords = NULL, ...)
{
  lines = readLines(con);
  if (!is.null(opening))
  {
    if (is.null(closing)) stop("If 'opening' is given, 'closing' must be given as well.");
    lines2 = character(0);
    nLines = length(lines);
    state = 0; index = 1; new1 = TRUE;
    while (index <= nLines)
    {
      if (new1) { l = lines[index]; new1 = FALSE; }
      if (state==0)
      {
        i = regexpr(opening, l, ...);
        if (i==-1) { 
          index = index + 1; 
          new1 = TRUE 
        } else { 
          state = 1;
          l = substring(l, i + attr(i, "match.length"));
        }
      }
      if (state==1)
      {
        i = regexpr(closing, l, ...);
        if (i==-1) {
          index = index + 1;
          new1 = TRUE;
          lines2 = c(lines2, l);
        } else {
          state = 0;
          if (i>1) lines2 = c(lines2, substring(l, 1, i-1));
          l = substring(l, i + attr(i, "match.length"));
        }
      }
    }
  } else lines2 = lines;

  words = unlist(strsplit(lines2, split = wordSep));
  candidates = lapply(organisms, convert2entrez, symbol = words);
  tables = mymapply(function(cands, org)
  {
    keep = !is.na(cands);
    out = data.frame(Word = words[keep], Entrez = cands[keep],
               Symbol = convert2symbol(entrez = cands[keep], organism = org));
    out[!duplicated(out$Word), ];
  }, candidates, organisms);

  names(tables) = organisms;
  tables;
}

#=========================================================================================================
#
# orthogonalize a set of vectors
#
#=========================================================================================================

gramSchmidt = function(x)
{
  nc = ncol(x);
  if (nc==1) return(x);
  xx = t(x)%*%x;
  det = (det(scale(xx, center = FALSE))/max(abs(scale(xx, center = FALSE)))^nc)^(1/nc)
  if (abs(det) < 1e-10) stop("The matrix is close to sigular.");
  for (cl in 2:nc)
  {
    x1 = x[, cl];
    for (i in 1:(cl-1))
      x1 = x1 - x[, i] * sum(x[, i] * x1)/sum(x[, i]^2);
    x[, cl] = x1;
  };
  x;
}

# Scaterplot margin parameters
scpp = function(mainMar = 2)
{
  par(mar = c(3.1, 3.1, mainMar, 1));
  par(mgp = c(2, 0.7, 0));
}

#=========================================================================================================
#
# My attempt at a Venn diagram function. Keyword: venn
#
#=========================================================================================================

vennDiagram.2sets = function(
  n1, n12, n2,
  
  mar = c(1, 1, 2, 1),

  setLabels,
  cex.labels = c(1,1),
  col.labels = c(1,1),

  pos.labels = c(270, 270),

  cex.counts = 1,
  col.counts = 1,

  areaScaleFunction = identity,
  main = "",
  ...)
{
  stop("haven't even started!")
}

  

#=========================================================================================================
#
# Permutation study of inter-sample distances
#
#=========================================================================================================

inAndOutGroupDistances = function(data, sameGroupIndicator = NULL, groupLabels)
{
  if (is.null(sameGroupIndicator))
  {
    if (nrow(data)!=length(groupLabels)) stop("nrow(data) and length(labels) must be the same.");
    sameGroupIndicator = c(as.dist(outer(groupLabels, groupLabels, FUN=`==`)));
  }
  dst = dist(data);
  out = tapply(c(dst), sameGroupIndicator, mean);
  names(out) = c("out", "in");
  out;
}

intersampleDistPermutationStudy = function(data, groupLabels, design = NULL, covariates = NULL,
       nPermutations = 1000, 
       randomSeed = 132,
       verbose = 1)
{
  if (nrow(data)!=length(groupLabels)) stop("nrow(data) and length(labels) must be the same.");
  if (!is.null(design)) {
    if (is.null(covariates))
      stop("When specifying 'design', also must specify 'covariates'.");
    data.corr = residuals(lm(as.formula(spaste("data~", design)), data = covariates ));
  } else
    data.corr = data;

  sameGroupIndicator = c(as.dist(outer(groupLabels, groupLabels, FUN=`==`)))

  observed = inAndOutGroupDistances(data.corr, sameGroupIndicator = sameGroupIndicator);

  permuted = matrix(NA, nPermutations, 2);
  colnames(permuted) = names(observed);

  if (!is.null(randomSeed))
  {
    .seed = .Random.seed;
    set.seed(randomSeed);
    on.exit(.Random.seed <<- .seed);
  }
  if (verbose > 0) pind = initProgInd();
  for (p in 1:nPermutations)
  {
    data.perm = data.corr[sample(1:nrow(data)), ];
    permuted[p, ] = inAndOutGroupDistances(data.perm, sameGroupIndicator = sameGroupIndicator);
    if (verbose > 0) pind = updateProgInd(p/nPermutations, pind);
  }

  if (verbose > 0) printFlush("");

  obsDiff = observed[1] - observed[2];
  permDiff = apply(permuted, 1, function(x) x[1]-x[2]);
  list(observedDifference = obsDiff,
       p.permutation = permutationPValues(c(obsDiff, obsDiff), cbind(permDiff, permDiff))[1],
       p.parametric = parametricPermutationPValues(c(obsDiff, obsDiff), cbind(permDiff, permDiff))$p[1],
       permutationDifferences = permDiff,
       observedMeans = observed, permutedMeans = permuted);
}

#=============================================================================================================
#
# Drop levels with few (typically 1) samples
#
#=============================================================================================================

dropLevelsWithFewSamples = function(x, min = 2)
{
  t = table(x);
  drop = names(t)[t<min];
  x[x%in%drop] = NA;
  x;
}

#=============================================================================================================
#
# cbindNonNULL
#
#=============================================================================================================

cbind.nonNULL = function(...)
{
  args = list(...);
  args = args[ sapply(args, length) > 0];
  do.call(cbind, args);
}

data.frame.nonNULL = function(...)
{
  args = list(...);
  args = args[ sapply(args, length) > 0];
  do.call(data.frame, args);
}

#=============================================================================================================
#
# changePoints
#
#=============================================================================================================

changepoints = function(x, first = FALSE, last = FALSE)
{
  n = length(x);
  out = which(x[-1]!=x[-n]) + 1;
  if (first) out = c(1, out);
  if (last) out = c(out, n+1);
  out;
}

#==============================================================================================================
#
# Functions to retrieve the origin counts for the top neighbors of each gene
#
#==============================================================================================================

getConsensusAndInputData = function(consensusInfo, calibrated = FALSE)
{
  consensus = BD.getData(consensusInfo$consensusData);
  if (calibrated) {
     inputs = list2multiData(lapply(consensusInfo$calibratedIndividualData, BD.getData));
  } else
     inputs = mtd.apply(consensusInfo$inputs, function(inp)
     {
       if (is.character(inp)) BD.getData(consensusInfo$individualTOMInfo$blockwiseAdjacencies[[inp]]$data) else
         BD.getData(inp$consensusData)
     });
  list(consensus = consensus, inputs = inputs);
}

traverseConsensusInfo = function(consensusInfo)
{
  out1 = list(consensusInfo);
  indivTOMInfo = consensusInfo$individualTOMInfo;
  isInfo = !mtd.apply(consensusInfo$inputs, is.character, mdaSimplify = TRUE)
  if (any(isInfo))
     out1 = c(out1, unlist(mtd.apply(consensusInfo$inputs[isInfo], traverseConsensusInfo, returnList = TRUE), 
                           recursive = FALSE));  
  if (!is.null(indivTOMInfo))
    out1 = lapply(out1, function(.info) 
            {
               if (is.null(.info$indivTOMInfo)) .info$individualTOMInfo = indivTOMInfo;
               .info;
            });
  out1
}

getOriginCountsForTopNeighbors = function(
   consensusTOMInfoList,
   calibrated = FALSE,
   nNeighbors)
{
   consensusInfos = list();
   nAnalyses = length(consensusTOMInfoList)
   for (ana in 1:nAnalyses)
   {
     tci = traverseConsensusInfo(consensusTOMInfoList[[ana]]);
     consensusInfos = c(consensusInfos,  tci);
   }
   nConsensusInfos = length(consensusInfos);
   originCounts = list();
   for (cii in 1:nConsensusInfos)
   {
     printFlush(spaste("Working on consensus ", cii))
     ci1 = getConsensusAndInputData(consensusInfos[[cii]], calibrated = calibrated);
     cmat = as.matrix(ci1$consensus);
     nGenes = ncol(cmat);
     p = nNeighbors/nGenes;
     q = colQuantileC(cmat, p = 1-p)
     topNeighbors = mymapply(function(col, .q) which(cmat[, col] >= .q),
                             1:nGenes, q);
     topNeighborConsensus = do.call(cbind, mymapply(function(col, i) cmat[i, col], 1:nGenes, topNeighbors))
     topNeighborTOMs = mtd.apply(ci1$inputs, function(inp)
     {
       mat1 = as.matrix(inp);
       do.call(cbind, mymapply(function(col, i) mat1[i, col], 1:nGenes, topNeighbors));
     });
     gc()
     
     originCounts[[cii]] = mtd.apply(topNeighborTOMs, function(x) colSums(x<=topNeighborConsensus));
   }
   originCounts.mean = lapply(originCounts, mtd.apply, mean, mdaSimplify = TRUE);
   list(originCounts = originCounts, meanOriginCounts = originCounts.mean);
}


#===================================================================================================================
#
# partial string match: a mildly complicated wrapper for grep that is kind of an inverse grep
#
#===================================================================================================================

partialMatch = function(text, patterns, value = FALSE, ...)
{
  out = as.numeric(sapply(text, function(s) 
  {
     out = which(sapply(patterns, grepl, s, ...));
     if (length(out)==0) out = NA;
     if (length(out)>1) warning(spaste("More than one match in ", s, ": using the first one."));
     out[1];
  }));
  if (value) out = patterns[out];
  out;
}

#===================================================================================================================
#
# some more variations on grep
#
#===================================================================================================================

grepv = function(...) grep(..., value = TRUE)
multiGrepv = function(...) multiGrep(..., value = TRUE);
     
#===================================================================================================================
#
# calculate and plot module overlaps
#
#===================================================================================================================

calculateAndPlotModuleOverlaps = function(
  labels1, labels2,
  ignoreLabels = 0,
  keepLabels1 = NULL,
  keepLabels2 = NULL,

  rowLabelPrefix = "M",
  colLabelPrefix = "M",
  rowLabelTranslationForColor = NULL,
  colLabelTranslationForColor = NULL,
  includeSize = TRUE,
  ...) # Arguments to plotOverlapHeatmap
{
  if (length(ignoreLabels)==0)  ignoreLabels = -324875910
  if (!is.null(keepLabels1))
    labels1 [ ! replaceMissing(labels1 %in% keepLabels1) ] = ignoreLabels[1];
    
  if (!is.null(keepLabels2))
    labels2 [ ! replaceMissing(labels2 %in% keepLabels2) ] = ignoreLabels[1];

  present = !is.na(labels1) & !is.na(labels2);
  labels1 = labels1[present];
  labels2 = labels2[present];

  tab1 = table(labels1);
  tab2 = table(labels2);

  ot = overlapTable(labels1, labels2, ignore = ignoreLabels, levels1 = keepLabels1, levels2 = keepLabels2);

  levels1 = rownames(ot$countTable);
  levels2 = colnames(ot$countTable);
  n1 = tab1[ match(levels1, names(tab1))];
  n2 = tab2[match(levels2, names(tab2))];

  if (is.null(rowLabelTranslationForColor)) rowLabelTranslationForColor = data.frame(levels1, as.numeric(levels1));
  if (is.null(colLabelTranslationForColor)) colLabelTranslationForColor = data.frame(levels2, as.numeric(levels2));

  plotOverlapHeatmap(
    countTab = ot$countTable,
    pTab = ot$pTable,
    colLabels = spaste("ME", labels2colors(translateUsingTable(levels2, colLabelTranslationForColor))),
    colSymbols = spaste(colLabelPrefix, levels2, if (includeSize) spaste( " (", n2, ")") else ""),
    rowLabels = spaste("ME", labels2colors(translateUsingTable(levels1, rowLabelTranslationForColor))),
    rowSymbols = spaste(rowLabelPrefix, levels1, if (includeSize) spaste( " (", n1, ")") else ""),
    ...);
}

   

plotOverlapHeatmap = function(
  countTab, pTab,
  rowLabels = rowNames(countTab),
  rowSymbols = NULL,
  colLabels = colnames(countTab),
  colSymbols = NULL,

  plotDevice = "pdf",
  plotFile = NULL,
  width = NULL,
  height = NULL,
  baseWidth = 2.5,
  baseHeight = 2.1,
  widthPerCol = 0.4,
  heightPerRow = 0.2,
  minWidth = 5,

  mar = c(7,7,2,1),

  threshold = 1,
  logpLimit = 30, 
  separatorInterval = 0,
  separatorInterval.col = separatorInterval,
  separatorInterval.row = separatorInterval,
  ...)
{
  ctab = countTab;
  ptab = pTab;
  lp = -log10(ptab);
  lp[lp > logpLimit] = logpLimit;
  keepRows = rowSums(ptab<=threshold) > 0;
  keepCols = colSums(ptab<=threshold) > 0;
  if (sum(keepRows)==0 || sum(keepCols)==0) return(NULL);

  ctab = ctab[keepRows, keepCols, drop = FALSE];
  lp = lp[keepRows, keepCols, drop = FALSE];
  ptab = ptab[keepRows, keepCols, drop = FALSE];
  text = ctab;
  text[ptab > threshold] = ""
  nr = nrow(ctab);
  nc = ncol(ctab);

  colSep = rowSep = NULL;
  if (separatorInterval.col > 0 && separatorInterval.col < nc) 
     colSep = seq(from = separatorInterval.col, to = nc-1, by = separatorInterval.col);
  if (separatorInterval.row > 0 && separatorInterval.row < nr) 
     rowSep = seq(from = separatorInterval.row, to = nr-1, by = separatorInterval.row);
  if (!is.null(plotDevice) && !is.null(plotFile))
  {
    dev = match.fun(plotDevice);
    dev(file = plotFile, width = max(minWidth, baseWidth + nc * widthPerCol), height = baseHeight + nr * heightPerRow);
    on.exit(dev.off());
  }
  par(mar = mar);
  out = labeledHeatmap(lp,
    xLabels = colLabels[keepCols],
    xSymbols = colSymbols[keepCols],
    yLabels = rowLabels[keepRows],
    ySymbols = rowSymbols[keepRows],
    textMatrix = text,
    colors = blueWhiteRed(100)[50:100],
    zlim = c(0, logpLimit),
    setStdMargins = FALSE,
    horizontalSeparator.y = rowSep,
    horizontalSeparator.col = "grey",
    horizontalSeparator.ext = 0.9,
    verticalSeparator.x = colSep,
    verticalSeparator.col = "grey",
    verticalSeparator.ext = 0.9,
    ...);

  list(countTable = ctab, pTable = ptab, heatmapData = out);
}


#==========================================================================================================
#
# Graph fuctions
#
#==========================================================================================================

graphMatrix = function(graph)
{
  lst = edgeL(graph);
  n = length(lst);
  mat = matrix(0, n, n);
  colnames(mat) = rownames(mat) = names(lst);
  for (i in 1:n) if (length(lst[[i]]$edges) > 0) mat[i, lst[[i]]$edges] = 1;
  mat
}

transitiveReduction = function(graph)
{
  graphMat = graphMatrix(graph);
  n = ncol(graphMat);
  for (i in 1:n) for (j in 1:n) if (graphMat[i,j] > 0)
  {
    gm1 = graphMat;
    gm1[i,j] = 0;
    gmx = gm1
    for (p in 2:n)
    {
       gmx = gmx %*% gm1
       if (gmx[i,j] > 0) break;
    }
    if (gmx[i,j] > 0) graphMat[i,j] = 0;
  }
  newList = edgeL(graph);
  for (i in 1:n)
    newList[[i]]$edges = which(graphMat[i, ]==1);
  graphNEL(nodes = nodes(graph), edgeL = newList, edgemode = "directed");
}

# this function creates a graphNEL from a named list of downstream targets
graphFromEdgeList = function(edgeList)
{
  nodes = unique(c(names(edgeList), unlist(edgeList)));
  edgeL = lapply(nodes, function(.node)
  {
    list(edges = if (.node %in% names(edgeList)) match(edgeList[[.node]], nodes) else numeric(0));
  });
  names(edgeL) = nodes;
  graphNEL(nodes = nodes, edgeL = edgeL, edgemode = "directed");
}

graphFromEdgeData = function(edgeData, nodeColumn = "Symbol", weightColumn = NULL)
{
  edgeData = edgeData[ sapply(edgeData, length) > 0 ]
  edgeData = edgeData[ sapply(edgeData, nrow) > 0 ]
  nodes = unique(c(names(edgeData), unlist(lapply(edgeData, getElement, nodeColumn))));
  edgeL = lapply(nodes, function(.node)
  {
    list(edges = if (.node %in% names(edgeData)) match(edgeData[[.node]]$Symbol, nodes) else numeric(0));
  });
  names(edgeL) = nodes;
  graphNEL(nodes = nodes, edgeL = edgeL, edgemode = "directed");
}

igraphFromEdgeData = function(edgeData, nodeColumn = "Symbol", weightColumn = NULL)
{
  edges = mymapply(function(from, to) 
        if (length(to) > 0) as.matrix(data.frame(from, to)) else NULL,
  from = names(edgeData), to = lapply(edgeData, getElement, nodeColumn));
  edges = do.call(rbind, edges[sapply(edges, length) > 0]);
  graph_from_edgelist(el = edges, directed = TRUE);
}



#==========================================================================================================
#
# Get proper deviance from output of DESeq
#
#==========================================================================================================

devianceForDESeq = function(ds, mc = mcols(ds), counts = counts(ds), weights = assays(ds)[["weights"]])
{
  if (length(weights) > 0)
  {
    if (!isTRUE(all.equal(dim(counts), dim(weights))))
      stop("'weights' and 'counts' must have the same dimensions.");
    if (length(weights)!=length(counts)) 
      stop("'weights' and 'counts' must have the same length.");
    deviance.saturated = -2*rowSums(weights * 
                dnbinom(x = counts, mu = counts, 
                        size = matrix(1/mc$dispersion, nrow = nrow(counts), ncol = ncol(counts)),
                        log = TRUE))
  } else 
    deviance.saturated = -2*rowSums(dnbinom(x = counts, mu = counts, 
                        size = matrix(1/mc$dispersion, nrow = nrow(counts), ncol = ncol(counts)),
                        log = TRUE))
  mc$deviance-deviance.saturated;
}

#==========================================================================================================
#
# Contrast numerator and denominator for a series of levels
#
#==========================================================================================================

contrastNumeratorsAndDenominators = function(levels)
{
  levels = rev(levels)
  n = length(levels)
  if (n<2) stop("Not enough levels.");
  num = den = character(n*(n-1)/2);
  index = 1;
  for (i in 1:(n-1)) for (j in (i+1):n)
  {
    num[index] = levels[i];
    den[index] = levels[j];
    index = index + 1;
  }
  list(num = num, den = den);
}

#==========================================================================================================
#
# Variance and standard deviation of an estimate of variance
#
#==========================================================================================================

varOfVariance = function(varEstimate, n) varEstimate^2 * 2/(n-1);
sdOfVariance = function(varEstimate, n) varEstimate * sqrt(2/(n-1));

sdOfStdDeviation = function(sdEstimate, n) sdEstimate * sqrt(1/(2*(n-1)));

#==========================================================================================================
#
# Variance and standard deviation of an estimate of variance
#
#==========================================================================================================

ZWeightFunction = function(z, Z0 = 3, power = 6) 
  tanh(tanh(z/Z0) ^ power);  ## The double use of tanh flattens the tails a bit

#==========================================================================================================
#
# Clip data
#
#==========================================================================================================

clipData = function(data, max)
{
  data[replaceMissing(data) < -max] = -max
  data[replaceMissing(data) > max] = max;
  data;
}

#==========================================================================================================
#
# Save collection in plain text tables
#
#==========================================================================================================

# A convenience function to automate what I always do when a new collection is created
exportCollectionToTextTables = function(
   collection, 
   dir, 
   nameBase,
   saveComponents = c("geneSetInfo", "geneSetContent", "groupInfo"))
{
  dfs = collection2dataFrames(collection);
  n = length(dfs);
  names = names(dfs);
  dir.create(dir, recursive = TRUE);
  for (c in 1:n)
  {
    if (names[c] %in% saveComponents)
       write.csv(dfs[[c]], file = gzfile(spaste(dir, "/", nameBase, "-", names[c], ".csv.gz")),
                 quote = TRUE, row.names = FALSE);
  }
}


#==========================================================================================================
#
# String similarity
#
#==========================================================================================================

stringDissimilarity = function(strings, distMethod = "manhattan")
{
  require(gtools)
  if (length(strings)==0) return(numeric(0));
  if (length(strings)==1) return(matrix(0, 1,1));
  nStrings = length(strings);
  strings2 = gsub("[^[:alnum:]]", " ", strings);
  strings2 = gsub(" +", " ", strings2);
  ascii = asc(strings2, simplify = FALSE);
  unique = unique(unlist(ascii));
  ascii2 = lapply(ascii, match, unique); 
  nUnique = length(unique);
  coocMats = list();
  space = match(asc(" ", simplify = FALSE)[[1]][1], unique, nomatch = -1);
  for (s in 1:nStrings)
  {
    coocMat1 = array(0, dim = c(nUnique, nUnique, nUnique));
    for (i in 1:(length(ascii[[s]])-2)) if (all(ascii2[[s]][i:(i+2)]!=space))
      coocMat1[ascii2[[s]][i], ascii2[[s]][i+1], ascii2[[s]][i+2]] = 
          coocMat1[ascii2[[s]][i], ascii2[[s]][i+1], ascii2[[s]][i+2]] + 1;
    coocMats[[s]] = as.vector(coocMat1);
  }
  coocMats = as.matrix(as.data.frame(coocMats));
  colnames(coocMats) = strings;
  dist(t(coocMats), method = distMethod);
}
  

#==========================================================================================================
#
# Alternative corAndPvalue
#
#==========================================================================================================

corAndNObs = function(x, y = NULL,
                        use = "pairwise.complete.obs",
                        weights.x = NULL, weights.y = NULL,
                        ...)
{
  cor = cor(x, y, use = use, weights.x = weights.x, weights.y = weights.y, ...);
  x = as.matrix(x);
  finMat = !is.na(x)
  if (!is.null(weights.x)) 
    finMat = finMat & weights.x > 0;
  if (is.null(y))
  {
    np = t(finMat) %*% finMat;
  } else {
    y = as.matrix(y);
    finMat.y = !is.na(y);
    if (!is.null(weights.y))
      finMat.y = finMat.y & weights.y > 0;
    np = t(finMat) %*% finMat.y
  }
  list(cor = cor, nObs = np);
}


#==========================================================================================================
#
# Convenience reading and saving functions
#
#==========================================================================================================

read.tsv = function(..., sep = "\t", header = TRUE)
  read.table(..., sep = sep, header = header);

write.csvr = function(...) write.csv(..., row.names = FALSE)

sourcee = function(...) source(..., echo = TRUE)


#==========================================================================================================
#
# Column and row indices of a selected elements of a matrix
#
#==========================================================================================================

rowAndColumnIndex = function(logMat)
{
  logMat = as.matrix(logMat);
  nr = nrow(logMat);
  i = which(logMat);
  row = (i-1)%%nr + 1;
  col = floor((i-1)/nr) + 1;
  cbind(row = row, col = col)
}

spacer = function(n) paste(rep(" ", round(n)), collapse = "");

decapitalize = function(s)
{
  change = which(s!="" & grepl("^.[a-z]", s));
  substr(s[change], 1, 1) = tolower(substr(s[change], 1, 1));
  s;
}

standardizeSampleNames = function(s) multiGSub(c("^X", " ", "-", "\\."), c("", "_", "_", "_"), s);



#==========================================================================================================
#
# Christmas tree plot of numbers of significant genes
#
#==========================================================================================================

plotNSignificant = function(
  nSignif,
  rowLabels,
  rowGroups = NULL,

  plotDevice = "pdf",
  plotDir = "Plots",
  plotFile = NULL,

  width = NULL,
  height = NULL,
  baseWidth = 4,
  baseHeight = 0.5,
  heightPerRow = 0.23,

  cex.rowLabels = 1,
  cex.groupLabels = 0.9,
  font.groupLabels = 2,

  col.right = "#FF9070", col.left = "skyblue",
  border.right = "brown", border.left = "blue",

  mar = c(0.5, NA, 2, 0.5),
  barGap = 0.2,

  sep.col = "grey30",

  ...)
{
  if (is.null(height)) height = baseHeight + heightPerRow * nrow(nSignif);

  mw = marginWidth(rowLabels, cex = cex.rowLabels, device = plotDevice);
  if (!is.null(rowGroups))
  {
    mw2 = marginWidth(rowGroups, cex = cex.groupLabels, device = plotDevice, font = font.groupLabels);
    mw$inch = mw$inch + mw2$inch;
    mw$lines = mw$lines + mw2$lines;
  }

  if (is.null(width)) width = baseWidth + mw$inch;
  if (is.na(mar[2])) mar[2] = mw$lines + 0.2;

  if (!is.null(plotFile) && !is.null(plotDevice))
  {
    plotFnc = match.fun(plotDevice);
    plotFnc(file = file.path(plotDir, plotFile), width = width, height = height);
    on.exit(dev.off());
  } 

  par(mar = mar);

  if (!is.null(rowGroups))
  {
    sep.pos = changepoints(rowGroups, first = TRUE)-1;
    sepLabels = rowGroups[sep.pos + 1];
  } else 
    sep.pos = NULL;

  nSep = length(sep.pos);
  plot = twoSideBarplot(nSignif[, 1], nSignif[, 2],
       col.left = col.left, col.right = col.right,
       border.left = border.left, border.right = border.right,
       yLabels = rowLabels, barGap = barGap,
       ..., 
       separatorPositions = sep.pos, sep.ext = TRUE, sep.col = sep.col);
  
   if (!is.null(sep.pos)) 
     text(rep(plot$box[1]-plot$leftMargin + 0.3 * strwidth("M"), nSep), plot$ytop[sep.pos+1],
          sepLabels, adj = c(0,1), cex = cex.groupLabels, xpd = TRUE, font = font.groupLabels)
}


#==========================================================================================================
#
# Generate list of loaded and attached packages via sessionInfo
#
#==========================================================================================================

listAttachedPackages = function(type = c("otherPkgs", "loadedOnly"), get = c("Version", "License"))
{
  si = sessionInfo();
  
  out = lapply(type, function(.comp)
  {
    pkg1 = names(si[[.comp]]);
    data = do.call(data.frame, lapply(get, function(.info) sapply(si[[.comp]], getElement, .info)));
    names(data) = get;
    out1 = data.frame(Package = pkg1, data)
  });
  setNames(out, type);
}

#==========================================================================================================
#
# Convenience function
#
#==========================================================================================================

data.frame.ncn = function(...) data.frame(..., check.names = FALSE);
