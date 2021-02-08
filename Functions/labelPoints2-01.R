
#=========================================================================================================
#
# labelExtremePoints
#
#=========================================================================================================

labelExtremePoints = function(...)
{
  labelExtremePoints.common(fnc = "labelPoints", ...);
}

labelExtremePoints2 = function(...)
{
  labelExtremePoints.common(fnc = "labelPoints2", ...);
}

labelExtremePoints.common = function(fnc, x, y, labels, nLabel, nConsider, 
                      directions = c("0+", "++", "+0", "+-", "0-", "--", "-0", "-+"), 
                      forceLabel = NULL, scaleForSelection = FALSE, pt.cex = 1,
                      cex = 1, font = 1, ptExtension = 1.1, col = 1,
                      ...)
{
  pt.cex = .extend(pt.cex, length(x));
  cex = .extend(cex, length(x));
  font = .extend(font, length(x));
  ptExtension = .extend(ptExtension, length(x));
  col = .extend(col, length(x));
  
  if (scaleForSelection) {
    xs = scale(x);
    ys = scale(y);
  } else {
    xs = x;
    ys = y;
  }
  if (nLabel > 0)
  {
    label = unique(c(
      if ("-0" %in% directions) order(xs)[1:nLabel] else NULL, 
      if ("+0" %in% directions) order(-xs)[1:nLabel] else NULL,
      if ("0-" %in% directions) order(ys)[1:nLabel] else NULL, 
      if ("0+" %in% directions) order(-ys)[1:nLabel] else NULL,
      if ("--" %in% directions) order(xs+ys)[1:nLabel] else NULL, 
      if ("++" %in% directions) order(-xs-ys)[1:nLabel] else NULL,
      if ("-+" %in% directions) order(xs-ys)[1:nLabel] else NULL, 
      if ("+-" %in% directions) order(-xs+ys)[1:nLabel] else NULL));
  } else label = numeric(0);

  if (nConsider > 0)
  {
    consider = unique(c(
      if ("-0" %in% directions) order(xs)[1:nConsider] else NULL, 
      if ("+0" %in% directions) order(-xs)[1:nConsider] else NULL,
      if ("0-" %in% directions) order(ys)[1:nConsider] else NULL, 
      if ("0+" %in% directions) order(-ys)[1:nConsider] else NULL,
      if ("--" %in% directions) order(xs+ys)[1:nConsider] else NULL, 
      if ("++" %in% directions) order(-xs-ys)[1:nConsider] else NULL,
      if ("-+" %in% directions) order(xs-ys)[1:nConsider] else NULL, 
      if ("+-" %in% directions) order(-xs+ys)[1:nConsider] else NULL));
  } else consider = numeric(0);

  if (is.character(forceLabel)) forceLabel = labels %in% forceLabel;
  if (is.logical(forceLabel)) forceLabel = which(forceLabel);
  label = unique(c(label, forceLabel))
  consider = unique(c(consider, forceLabel));

  labels.consider = labels[consider];
  x.c = x[consider];
  y.c = y[consider];
  labels.consider[ !consider %in% label] = "";
  do.call(match.fun(fnc), list(x = x.c, y = y.c, labels = labels.consider, pt.cex = pt.cex[consider], 
               cex = cex[consider], font = font[consider], ptExtension = ptExtension[consider], col = col[consider], ...))
}

fillLabelSpaceWithPoints = function(x, y, labels, cex, ...)
{
  out.x = numeric(0);
  out.y = numeric(0);
  fin = is.finite(x) & is.finite(y);
  x = x[fin];
  y = y[fin];
  labels = labels[fin];
  h = strheight(labels, cex = cex, ...);
  w = strwidth(labels, cex = cex, ...);
  for (i in 1:length(labels)) if (nchar(labels[i]) > 0)
  {
    ny = ceil(cex) + 1;
    nx = ceil(cex * nchar(labels[i]) + 1);

    xl = x[i] - w[i]/2;
    xr = x[i] + w[i]/2;
    yb = y[i] - h[i]/2;
    yt = y[i] + h[i]/2;
    out.x = c(out.x, rep(seq(from = xl, to = xr, length.out = nx), each = ny));
    out.y = c(out.y, rep(seq(from = yb, to = yt, length.out = ny), nx));
  }
  list(x = out.x, y = out.y);
}


#=================================================================================================
#
# labelPoints2: label points in a scatterplot while trying to avoid labels overlapping with one another
# and with points.
#
#=================================================================================================

labelOverlap = function(
    x.ref, y.ref, adj.ref = c(0.5, 0.5), text.ref, 
    cex=1, font=1, 
    widths.ref = strwidth(text.ref, cex = cex, font = font),
    heights.ref = strheight(text.ref, cex = cex, font = font),
    x.test, y.test, adj.test = c(0.5, 0.5),
    text.test, cex.test = cex, font.test = font,
    widths.test = strwidth(text.test, cex = cex.test, font = font.test),
    heights.test = strheight(text.test, cex = cex.test, font = font.test)
)
{
  nRef = length(widths.ref);
  if (nRef !=length(x.ref)) stop("Lengths of 'text.ref' and 'x.ref' must be the same.");
  if (nRef !=length(y.ref)) stop("Lengths of 'text.ref' and 'y.ref' must be the same.");
  if (nRef !=length(heights.ref)) stop("Lengths of 'text.ref' and 'heights.ref' must be the same.");

  nTest = length(widths.test);
  if (nTest !=length(x.test)) stop("Lengths of 'widths.test' and 'x.test' must be the same.");
  if (nTest !=length(y.test)) stop("Lengths of 'widths.test' and 'y.test' must be the same.");
  if (nTest !=length(heights.test)) stop("Lengths of 'widths.test' and 'heights.test' must be the same.");

  xl.ref = x.ref - adj.ref[1] * widths.ref;
  xr.ref = x.ref + (1-adj.ref[1]) * widths.ref;
  yb.ref = y.ref - adj.ref[2] * heights.ref;
  yt.ref = y.ref + (1-adj.ref[2]) * heights.ref;

  xl.test = x.test - adj.test[1] * widths.test;
  xr.test = x.test  + (1 - adj.test[1]) * widths.test;
  yb.test = y.test - adj.test[2] * heights.test;
  yt.test = y.test + (1 - adj.test[2]) * heights.test;

  xl.ref.x = rep(xl.ref, each = nTest);
  xr.ref.x = rep(xr.ref, each = nTest);
  yb.ref.x = rep(yb.ref, each = nTest);
  yt.ref.x = rep(yt.ref, each = nTest);
  
  xl.test.x = rep(xl.test, nRef);
  xr.test.x = rep(xr.test, nRef);
  yb.test.x = rep(yb.test, nRef);
  yt.test.x = rep(yt.test, nRef);

  overlap.x = (xl.ref.x >= xl.test.x & xl.ref.x <= xr.test.x) | (xr.ref.x >= xl.test.x & xr.ref.x <= xr.test.x) | 
              (xl.test.x >= xl.ref.x & xl.test.x <= xr.ref.x);

  overlap.y = (yb.ref.x >= yb.test.x & yb.ref.x <= yt.test.x) | (yt.ref.x >= yb.test.x & yt.ref.x <= yt.test.x) | 
              (yb.test.x >= yb.ref.x & yb.test.x <= yt.ref.x);

  overlap = overlap.x & overlap.y;
  dim(overlap) = c(nTest, nRef)

  overlap;
}


labelDistance = function(
    x.ref, y.ref, adj.ref = c(0.5, 0.5), text.ref,
    cex=1, font=1,
    widths.ref = strwidth(text.ref, cex = cex, font = font),
    heights.ref = strheight(text.ref, cex = cex, font = font),
    x.test, y.test, adj.test = c(0.5, 0.5),
    text.test, cex.test = cex, font.test = font,
    widths.test = strwidth(text.ref, cex = cex.test, font = font.test),
    heights.test = strheight(text.ref, cex = cex.test, font = font.test)
)
{
  nRef = length(widths.ref);
  if (nRef !=length(x.ref)) stop("Lengths of 'widths.ref' and 'x.ref' must be the same.");
  if (nRef !=length(y.ref)) stop("Lengths of 'widths.ref' and 'y.ref' must be the same.");
  if (nRef !=length(heights.ref)) stop("Lengths of 'widths.ref' and 'heights.ref' must be the same.");

  nTest = length(widths.test);
  if (nTest !=length(x.test)) stop("Lengths of 'widths.test' and 'x.test' must be the same.");
  if (nTest !=length(y.test)) stop("Lengths of 'widths.test' and 'y.test' must be the same.");
  if (nTest !=length(heights.test)) stop("Lengths of 'widths.test' and 'heights.test' must be the same.");

  xl.ref = x.ref - adj.ref[1] * widths.ref;
  xr.ref = x.ref + (1-adj.ref[1]) * widths.ref;
  yb.ref = y.ref - adj.ref[2] * heights.ref;
  yt.ref = y.ref + (1-adj.ref[2]) * heights.ref;

  xl.test = x.test - adj.test[1] * widths.test;
  xr.test = x.test  + (1 - adj.test[1]) * widths.test;
  yb.test = y.test - adj.test[2] * heights.test;
  yt.test = y.test + (1 - adj.test[2]) * heights.test;

  xl.ref.x = rep(xl.ref, each = nTest);
  xr.ref.x = rep(xr.ref, each = nTest);
  yb.ref.x = rep(yb.ref, each = nTest);
  yt.ref.x = rep(yt.ref, each = nTest);

  xl.test.x = rep(xl.test, nRef);
  xr.test.x = rep(xr.test, nRef);
  yb.test.x = rep(yb.test, nRef);
  yt.test.x = rep(yt.test, nRef);

  #overlap.x = (xl.ref.x >= xl.test.x & xl.ref.x <= xr.test.x) | (xr.ref.x >= xl.test.x & xr.ref.x <= xr.test.x) | 
  #            (xl.test.x >= xl.ref.x & xl.test.x <= xr.ref.x);

  #overlap.y = (yb.ref.x >= yb.test.x & yb.ref.x <= yt.test.x) | (yt.ref.x >= yb.test.x & yt.ref.x <= yt.test.x) | 
  #            (yb.test.x >= yb.ref.x & yb.test.x <= yt.ref.x);

  #dist.x = (1-overlap.x) * pmin(abs(xl.ref.x - xr.test.x), abs(xr.ref.x - xl.test.x))
  #dist.y = (1-overlap.y) * pmin(abs(yb.ref.x - yt.test.x), abs(yt.ref.x - yb.test.x));

  # the following overlaps are positive if the intervals don't overlap, and negative if they do overlap. The negative
  # number is (up to a sign) the overlap length.
  ovrLen.x = pmax(xl.ref.x, xl.test.x) - pmin(xr.ref.x, xr.test.x);
  ovrLen.y = pmax(yb.ref.x, yb.test.x) - pmin(yt.ref.x, yt.test.x);

  # If only one interval overlaps, 
  dist = ifelse( ovrLen.x * ovrLen.y > 0, sign(ovrLen.x) * sqrt(ovrLen.x^2 + ovrLen.y^2), pmax(ovrLen.x, ovrLen.y));

  dim(dist) = c(nTest, nRef)
  dist;
}


# This function returns, for each label (position and width/height), distance to the nearest point (zero if
# there is an ovrelapping point).

labelPointDistance = function(
  # Coordinates of the label(s)
  x.ref, y.ref, 
  adj.ref = c(0.5, 0.5), cex=1, font=1, 
  # Optional label(s) 
  text.ref, 
  # width and height of labels
  widths.ref = strwidth(text.ref, cex = cex, font = font),
  heights.ref = strheight(text.ref, cex = cex, font = font),
  # point coordinates
  x.pt, y.pt, pt.cex, 
  ratio.pointToChar = 0.3,
  ratio.pointToChar.width = ratio.pointToChar,
  ratio.pointToChar.height = ratio.pointToChar)
{
  nRef = length(widths.ref);
  if (nRef !=length(x.ref)) stop("Lengths of 'text.ref' and 'x.ref' must be the same.");
  if (nRef !=length(y.ref)) stop("Lengths of 'text.ref' and 'y.ref' must be the same.");

  nTest = length(x.pt);
  if (nTest !=length(y.pt)) stop("Lengths of 'text.test' and 'x.test' must be the same.");
  pt.cex = checkOrExtend(pt.cex, nTest, "pt.cex");

  xl.ref = x.ref - adj.ref[1] * widths.ref;
  xr.ref = x.ref + (1-adj.ref[1]) * widths.ref;
  yb.ref = y.ref - adj.ref[2] * heights.ref;
  yt.ref = y.ref + (1-adj.ref[2]) * heights.ref;

  pt.width = ratio.pointToChar.width * par("cxy")[1];
  pt.height = ratio.pointToChar.height * par("cxy")[2];

  #browser()
  xl.pt = x.pt - pt.width;
  xr.pt = x.pt + pt.width;
  yb.pt = y.pt - pt.height;
  yt.pt = y.pt + pt.height;
  
  distance = mapply(function(xl, xr, yb, yt)   ## use mapply here to simplify to a vector
  {
     overlap.x = (xl < xl.pt & xr > xl.pt | xl < xr.pt & xr > xr.pt |
                   xl > xl.pt & xl < xr.pt | xr > xl.pt & xr < xr.pt);
     dist.x = (1-overlap.x) * pmin(abs(xl-xl.pt), abs(xl-xr.pt), abs(xr - xr.pt), abs(xr-xl.pt));
     overlap.y = (yb < yb.pt & yt > yb.pt | yb < yt.pt & yt > yt.pt | 
                   yb > yb.pt & yb < yt.pt | yt > yb.pt & yt < yt.pt);
     dist.y = (1-overlap.y) * pmin(abs(yb-yb.pt), abs(yb-yt.pt), abs(yt - yb.pt), abs(yt - yt.pt));
     sqrt(min(dist.x^2 + dist.y^2));
  }, xl.ref, xr.ref, yb.ref, yt.ref);

  distance
}

# from the given query points given by x, y
maxDistIndex = function( 
    x1, y1, adj1 = c(0.5, 0.5),
    width1, height1,
    x2, y2, adj2 = c(0.5, 0.5),
    widths2, heights2)
{
  n = length(x1);
  dists = labelDistance(
     x.ref = x1, y.ref = y1, 
     adj.ref = adj1,
     widths.ref = rep(width1, n), 
     heights.ref = rep(height1, n),
     x.test = x2, y.test = y2, adj.test = adj2,
     widths.test = widths2, heights.test = heights2);
  minDists = colMins(dists);
  which.max(minDists);
}

if (FALSE)
{
  labelDistance(x.ref = x1[7], y.ref = y1[7], widths.ref = width1,
             heights.ref = height1, x.test = x2[5], y.test = y2[5], widths.test = widths2[5], heights.test = heights2[5])
}


# Make sure that the plot boundaries are ordered, i.e., xmax > xmin (box[2] > box[1]) and siumilarly for y dimension.
# If an axis is reversed, so will be the corresponding elements of par("usr").

normalizePlotBoundaries = function(box)
{
  if (box[2]<box[1]) box[1:2] = box[2:1];  
  if (box[4]<box[3]) box[3:4] = box[4:3];  
  box;
}

# Within plot boundaries: which positions are entirely within plot.
withinPlotBoundaries = function(x.ref, y.ref,
          adj.ref = c(0.5, 0.5), widths.ref, heights.ref,
          box = par("usr"))
{
  box = normalizePlotBoundaries(box);
  xl.ref = x.ref - adj.ref[1] * widths.ref;
  xr.ref = x.ref + (1-adj.ref[1]) * widths.ref;
  yb.ref = y.ref - adj.ref[2] * heights.ref;
  yt.ref = y.ref + (1-adj.ref[2]) * heights.ref;

  within = xl.ref >= box[1] & xr.ref <= box[2] & yb.ref >= box[3] & yt.ref <= box[4];
  within;
}

# generate potential positions of a label around a point. Argument adj adjusts the output so that it can be
# used to place the label with the given adj argument. 

potentialPositions = function(
  x, y, # point coordinates
  ptw, pth,  # point half-width and half-height, user coordinates
  labw, labh, # label width and height, user coordinates
  step.x, step.y,
  adj = c(0.5, 0.5),
  allowDistantLabels = 0,
  penaltyFnc = NULL,
  penaltyArgs = list())
{
  if (any(is.na(c(x,y)))) return(data.frame(x = numeric(0), y = numeric(0)));
  out.x = out.y = penalty = numeric(0);
  for (ds in 0:allowDistantLabels)
  {
    xl = x - adj[1] * labw - ds * step.x - ptw;
    xr = xl + labw  + 2*ds*step.x + 2*ptw;
    yb = y - adj[2] * labh -ds * step.y - pth;
    yt = yb + labh + 2*pth + 2*ds*step.y;

    nx = max(2, round((xr - xl)/step.x));
    ny = max(2, round((xr - xl)/step.y));

    #aa = try( {
    out.x = c(out.x, seq(from = xl, to = xr, length.out = nx)[-1],
          rep(xr, ny-1),
          seq(from = xr, to = xl, length.out = nx)[-1],
          rep(xl, ny-1));
    out.y = c(out.y, rep(yb, nx-1),
          seq(from=yb, to=yt, length.out = ny)[-1],
          rep(yt, nx-1), 
          seq(from=yt, to = yb, length.out = ny)[-1]);
    penalty = c(penalty, rep(ds, 2*(nx + ny-2)));
  }
  keep = withinPlotBoundaries(out.x, out.y, adj.ref = adj, widths.ref = labw, heights.ref = labh);
  out.x = out.x[keep];
  out.y = out.y[keep];
  penalty = penalty[keep];
  #}); 
  #if (inherits(aa, "try-error")) browser()
  if (length(penaltyFnc)>0)
    penalty = do.call(match.fun(penaltyFnc), c(list(x = out.x, y = out.y, basePenalty = penalty), penaltyArgs));
  data.frame(x = out.x, y = out.y, penalty = penalty);
}



labelPoints2 = function(
   x, y, labels, 
   cex = 1, font = 1, 
   pt.cex = 1, ptExtension = 1.1, 
   col = 1, 
   xpd = TRUE,
   protectEdges = TRUE, doPlot = TRUE, 
   ratio.pointToChar = 0.3,
   ratio.pointToChar.width = ratio.pointToChar,
   ratio.pointToChar.height = ratio.pointToChar,
   testStep.x = strwidth("M", cex = cex)/2,
   testStep.y = strheight("M", cex = cex)/2,
   labelPreference = NULL,
   allowDistantLabels = 0, 
      ## Non-zero here allows the labels to be positioned away from points with maximum distance in units of testStep.x
      ## and testStep.y above
   maxTries = 10,
   penaltyFunction = NULL,  ## Not implemented yet
   penaltyArgs = list(),
   verbose = 0, indent = 0,
   ...)
{
  spaces = indentSpaces(indent);
  if (length(x)!=length(y)) stop("Lengths of x and y must be the same.");
  if (length(x)!=length(labels)) stop("Lengths of x and labels must be the same.")

  extraArgs = list(...)
  if ("offs" %in% names(extraArgs))
  {
    ratio.pointToChar = 0.3 * extraArgs$offs/0.07;
    ratio.pointToChar.width = ratio.pointToChar;
    ratio.pointToChar.height = ratio.pointToChar;
  }

  pt.cex = .extend(pt.cex, length(x));
  ptExtension = .extend(ptExtension, length(x));
  col = .extend(col, length(x));
  cex = .extend(cex, length(x));
  font = .extend(font, length(x));
  ratio.pointToChar = .extend(ratio.pointToChar, length(x));
  ratio.pointToChar.width = .extend(ratio.pointToChar.width, length(x));
  ratio.pointToChar.height = .extend(ratio.pointToChar.height, length(x));

  #cex = .extend(cex, length(x));

  fin = is.finite(x) & is.finite(y);
  x = x[fin];
  y = y[fin];
  labels = replaceMissing(labels[fin]);
  pt.cex = pt.cex[fin];
  ptExtension = ptExtension[fin];
  ratio.pointToChar = ratio.pointToChar[fin]
  ratio.pointToChar.width = ratio.pointToChar.width[fin]
  ratio.pointToChar.height = ratio.pointToChar.height[fin]
  cex = cex[fin];
  font = font[fin];

  if ((length(labels)==0) || all(labels=="", na.rm = TRUE))
    return(invisible(list(x = numeric(0), y = numeric(0), label = character(0), plotLabels = character(0),
                 initial = list(x = numeric(0), y = numeric(0),
                                labels = character(0), plotLabels = character(0)))));

  nPoints = length(labels);
  if (nPoints==0) return(0);
  if (is.null(labelPreference)) labelPreference = rep(1, nPoints);

  #box = par("usr");
  dims = par("pin");
  scaleX = scaleY = 1;

  if (par("xlog"))
  {
     xx = log10(x);
  } else
     xx = x;

  if (par("ylog"))
  {
     yy = log10(y);
  } else
     yy = y;

  pt.width = pt.cex * ratio.pointToChar.width * strwidth("o");
  pt.height = pt.cex * ratio.pointToChar.height * strheight("o");

  xx = xx * scaleX;
  yy = yy * scaleY;

  labWidth = as.numeric(mapply(strwidth, labels, cex=cex, font = font)) * scaleX; 
  labHeight = as.numeric(mapply(strheight, labels, cex=cex, font = font)) * scaleY; 
  plotLabels = labels!="";
  testPositions = list();
  testPoints = numeric(0);

  # Generate possible positions for each label: positions that do not overlap with any of the points and don't cross the
  # edges of the plot.

  if (verbose > 0) printFlush(spaste(spaces, " ..generating potential positions.."));
  for (p in which(plotLabels))
  {
    testPositions1 = potentialPositions(
        xx[p], yy[p], 
        labw = labWidth[p],
        labh = labHeight[p],
        ptw = pt.width[p] * ptExtension[p],
        pth = pt.height[p] * ptExtension[p],
        step.x = testStep.x,
        step.y = testStep.y,
        adj = c(0.5, 0.5),
        allowDistantLabels = allowDistantLabels);
    nPositions1 = nrow(testPositions1);
    if (FALSE)
    {
      plot(x, y)
      points(testPositions1$x, testPositions1$y, pch = 21, bg  = c(1,2,3,4,5));
      i = nPositions1;
      text(testPositions1$x[i], testPositions1$y[i], labels[p], cex = cex[i], font = font[i])
      rect(testPositions1$x[i] - labWidth[p]/2, testPositions1$y[i] - labHeight[p]/2, testPositions1$x[i] + labWidth[p]/2, 
           testPositions1$y[i] + labHeight[p]/2);
    }
    dist1 = labelPointDistance(
          x.ref = testPositions1$x, y.ref = testPositions1$y, 
          adj.ref = c(0.5, 0.5), 
          widths.ref = rep(labWidth[p], nPositions1), heights.ref = rep(labHeight[p], nPositions1),
          x.pt = xx[-p], y.pt = yy[-p], pt.cex = pt.cex[-p] * ptExtension[-p], 
          ratio.pointToChar.width = ratio.pointToChar.width[-p],
          ratio.pointToChar.height = ratio.pointToChar.height[-p]);
    keep = dist1 > 0;
    if (FALSE)
    {
      plot(x, y)
      points(testPositions1$x[keep], testPositions1$y[keep], pch = 21, bg  = c(1,2,3,4,5));
      for (i in which(keep))
      {
      #text(testPositions1$x[i], testPositions1$y[i], labels[p], cex = cex, font = font)
        rect(testPositions1$x[i] - labWidth[p]/2, testPositions1$y[i] - labHeight[p]/2, testPositions1$x[i] + labWidth[p]/2, 
             testPositions1$y[i] + labHeight[p]/2);
      }
    }
    #zzzzz = any(keep);
    #if (is.na(zzzzz)) browser()
    if (any(keep))
    {
      testPositions = c(testPositions, list(cbind(testPositions1[keep, ], dst = dist1[keep])));
      testPoints = c(testPoints, p);
    } else 
      plotLabels[p] = FALSE;
  }

  # Initial label positions are those that minimize penalty and then maximize distances from the nearest point of each label.
  initialPositions = do.call(rbind, lapply(testPositions, function(pos) 
  {
    mp = min(pos$penalty);
    cand = which(pos$penalty==mp);
    as.matrix(pos)[ cand[which.max(pos$dst[cand])], 1:2]
  }));
  plotLabels.initial = plotLabels;
  p.x = initialPositions[, 1];
  p.y = initialPositions[, 2];

  if (FALSE)
  {
    plot(x, y, cex = pt.cex)
    text(initialPositions[, 1], initialPositions[, 2], labels[plotLabels], cex = cex, adj = c(0.5, 0.5));
  }
  labW = labWidth[testPoints];
  labH = labHeight[testPoints];

  # Get label-label distances and see if any of the labels overlap.
  labelDists = labelDistance(
    x.ref = p.x, y.ref = p.y, 
    widths.ref = labW, heights.ref = labH,
    x.test = p.x, y.test = p.y, 
    widths.test = labW, heights.test = labH);
  #labelDists = round(labelDists, 2)
  diag(labelDists) = 1;
  overlaps = labelDists<0;
  done = !any(overlaps);

  # If there are overlaps, try moving labels around.
  tries = 1;
  while (!done)
  {
    if (verbose > 0) printFlush(spaste(spaces, " ..attempting to move labels around: step ", tries));
    run = 0;
    overlaps = labelDists<0;
    nOverlaps1 = colSums(overlaps);
    done = !any(overlaps);
    nAvailablePos1 = sapply(testPositions, nrow);
    if (verbose > 1) pind = initProgInd();
    while (run < maxTries && !done)
    {
      if (verbose > 1) pind = updateProgInd(run/maxTries, pind);
      testOrder = order(-nOverlaps1, -nAvailablePos1);
      # Try moving each point such that the distance to other labels is maximized.
      ip = 1; 
      nTestPoints = length(testOrder);
      while (ip <= nTestPoints && !done)
      {
        p = testOrder[ip];
        i = maxDistIndex(
           x1 = testPositions[[p]]$x, y1 = testPositions[[p]]$y,
           width1 = labW[p], height1 = labH[p],
           x2 = p.x[-p], y2 = p.y[-p], 
           widths2 = labW[-p], heights2 = labH[-p]);
        p.x[p] = testPositions[[p]]$x[i];
        p.y[p] = testPositions[[p]]$y[i];
        labelDists[-p, p] = labelDistance(x.ref = p.x[p], y.ref = p.y[p],
           widths.ref = labW[p], heights.ref = labH[p],
           x.test = p.x[-p], y.test = p.y[-p],
           widths.test = labW[-p], heights.test = labH[-p])
        done = all(labelDists > 0)
        ip = ip + 1;
      }
      if (FALSE)
      {
         sizeGrWindow(12, 6);
         par(mfrow = c(1,2))
         par(cex = 1)
         plot(x, y, cex = pt.cex)
         text(initialPositions[, 1], initialPositions[, 2], labels[plotLabels], cex = cex, adj = c(0.5, 0.5));
         plot(x, y, cex = pt.cex)
         text(p.x, p.y, labels[plotLabels], cex = cex, adj = c(0.5, 0.5));
       }
      run = run + 1;
    }
    if (verbose > 1) { pind = updateProgInd(1, pind); printFlush()}
    if (!done)
    {
      # need to remove one of the offending labels.
      candidates = which(colSums(labelDists<0) > 0);
      negativeSpace = colSums( (labelDists < 0) * labelDists);
      removeOrder = order(labelPreference[testPoints], negativeSpace, -labW);
      remove = removeOrder[1];
      remove.origIndex = testPoints[remove];
      testPoints = testPoints[-remove];
      plotLabels[remove.origIndex] = FALSE;
      labW = labW[-remove];
      labH = labH[-remove];
      labelDists = labelDists[-remove, -remove];
      p.x = p.x[-remove];
      p.y = p.y[-remove];
      testPositions = testPositions[-remove];
    }
    tries = tries + 1
  }

  if (par("xlog")) p.x = 10^p.x;
  if (par("ylog")) p.y = 10^p.y;

  if (doPlot && any(plotLabels))
    text(p.x, p.y, labels[plotLabels], col = col[plotLabels], cex = cex, xpd = xpd, adj = c(0.5, 0.5), font = font, ...) 

  # For compatibility with older labelPoints, the return value should contain x, y, label at top level.
  invisible(list(x = p.x, y = p.y, label = labels[plotLabels], plotLabels = plotLabels, 
                 initial = list(x = initialPositions[, 1], y = initialPositions[, 2],
                                labels = labels[plotLabels.initial], plotLabels = plotLabels.initial)));
}


#=================================================================================================================
#
# Another point label function
#
#=================================================================================================================

# this function labels a set of selected points in a larger set of points. The selected points are all labeled; the
# positions are chosen with increasing radius from the labeled points until one position fits.


### not even started on this... delete if no work done in the near future.
labelPoints.mustLabel = function(
   x, y, 
   labelIndex, labels, 
   cex = 1, font = 1, 
   pt.cex = 1, ptExtension = 1.1, 
   xpd = TRUE,
   protectEdges = TRUE, doPlot = TRUE, 
   ratio.pointToChar = 0.3,
   ratio.pointToChar.width = ratio.pointToChar,
   ratio.pointToChar.height = ratio.pointToChar,
   testStep.x = strwidth("M", cex = cex)/2,
   testStep.y = strheight("M", cex = cex)/2,
   labelPreference = NULL,
   maxTries = 10,
   ...)
{
  fin = is.finite(x) & is.finite(y);
  x = x[fin];
  y = y[fin];
  labels = replaceMissing(labels[fin]);

  if ((length(labels)==0) || all(labels=="", na.rm = TRUE))
    return(invisible(list(x = numeric(0), y = numeric(0), label = character(0), plotLabels = character(0),
                 initial = list(x = numeric(0), y = numeric(0),
                                labels = character(0), plotLabels = character(0)))));

  extraArgs = list(...)
  if ("offs" %in% names(extraArgs))
  {
    ratio.pointToChar = 0.3 * extraArgs$offs/0.07;
    ratio.pointToChar.width = ratio.pointToChar;
    ratio.pointToChar.height = ratio.pointToChar;
  }
  nPoints = length(labels);
  if (nPoints==0) return(0);
  if (is.null(labelPreference)) labelPreference = rep(1, nPoints);

  #box = par("usr");
  dims = par("pin");
  scaleX = scaleY = 1;

  if (par("xlog"))
  {
     xx = log10(x);
  } else
     xx = x;

  if (par("ylog"))
  {
     yy = log10(y);
  } else
     yy = y;

  pt.width = ratio.pointToChar.width * strwidth("o");
  pt.height = ratio.pointToChar.height * strheight("o");

  xx = xx * scaleX;
  yy = yy * scaleY;

  labWidth = strwidth(labels, cex=cex, font = font) * scaleX; 
  labHeight = strheight(labels, cex=cex, font = font) * scaleY; 
  plotLabels = labels!="";
  testPositions = list();
  testPoints = numeric(0);

  # Generate possible positions for each label: positions that do not overlap with any of the points and don't cross the
  # edges of the plot.

  for (p in which(plotLabels))
  {
    testPositions1 = potentialPositions(
        xx[p], yy[p], 
        labw = labWidth[p],
        labh = labHeight[p],
        ptw = pt.width * ptExtension,
        pth = pt.height * ptExtension,
        step.x = testStep.x,
        step.y = testStep.y,
        adj = c(0.5, 0.5),
        allowDistantLabels = allowDistantLabels);
    nPositions1 = nrow(testPositions1);
    if (FALSE)
    {
      plot(x, y)
      points(testPositions1$x, testPositions1$y, pch = 21, bg  = c(1,2,3,4,5));
      i = nPositions1;
      text(testPositions1$x[i], testPositions1$y[i], labels[p], cex = cex, font = font)
      rect(testPositions1$x[i] - labWidth[p]/2, testPositions1$y[i] - labHeight[p]/2, testPositions1$x[i] + labWidth[p]/2, 
           testPositions1$y[i] + labHeight[p]/2);
    }
    dist1 = labelPointDistance(
          x.ref = testPositions1$x, y.ref = testPositions1$y, 
          adj.ref = c(0.5, 0.5), 
          widths.ref = rep(labWidth[p], nPositions1), heights.ref = rep(labHeight[p], nPositions1),
          x.pt = xx[-p], y.pt = yy[-p], pt.cex = cex * ptExtension, 
          ratio.pointToChar.width = ratio.pointToChar.width,
          ratio.pointToChar.height = ratio.pointToChar.height);
    keep = dist1 > 0;
    if (FALSE)
    {
      plot(x, y)
      points(testPositions1$x[keep], testPositions1$y[keep], pch = 21, bg  = c(1,2,3,4,5));
      for (i in which(keep))
      {
      #text(testPositions1$x[i], testPositions1$y[i], labels[p], cex = cex, font = font)
        rect(testPositions1$x[i] - labWidth[p]/2, testPositions1$y[i] - labHeight[p]/2, testPositions1$x[i] + labWidth[p]/2, 
             testPositions1$y[i] + labHeight[p]/2);
      }
    }
    #zzzzz = any(keep);
    #if (is.na(zzzzz)) browser()
    if (any(keep))
    {
      testPositions = c(testPositions, list(cbind(testPositions1[keep, ], dst = dist1[keep])));
      testPoints = c(testPoints, p);
    } else 
      plotLabels[p] = FALSE;
  }

  # Initial label positions are those that maximize distances from the nearest point of each label.
  initialPositions = do.call(rbind, lapply(testPositions, function(pos) as.matrix(pos)[ which.max(pos$dst), 1:2]));
  plotLabels.initial = plotLabels;
  p.x = initialPositions[, 1];
  p.y = initialPositions[, 2];

  if (FALSE)
  {
    plot(x, y, cex = pt.cex)
    text(initialPositions[, 1], initialPositions[, 2], labels[plotLabels], cex = cex, adj = c(0.5, 0.5));
  }
  labW = labWidth[testPoints];
  labH = labHeight[testPoints];

  # Get label-label distances and see if any of the labels overlap.
  labelDists = labelDistance(
    x.ref = p.x, y.ref = p.y, 
    widths.ref = labW, heights.ref = labH,
    x.test = p.x, y.test = p.y, 
    widths.test = labW, heights.test = labH);
  #labelDists = round(labelDists, 2)
  diag(labelDists) = 1;
  overlaps = labelDists<0;
  done = !any(overlaps);

  # If there are overlaps, try moving labels around.
  while (!done)
  {
    run = 0;
    overlaps = labelDists<0;
    nOverlaps1 = colSums(overlaps);
    done = !any(overlaps);
    nAvailablePos1 = sapply(testPositions, nrow);
    while (run < maxTries && !done)
    {
      testOrder = order(-nOverlaps1, -nAvailablePos1);
      # Try moving each point such that the distance to other labels is maximized.
      ip = 1; 
      nTestPoints = length(testOrder);
      while (ip <= nTestPoints && !done)
      {
        p = testOrder[ip];
        i = maxDistIndex(
           x1 = testPositions[[p]]$x, y1 = testPositions[[p]]$y,
           width1 = labW[p], height1 = labH[p],
           x2 = p.x[-p], y2 = p.y[-p], 
           widths2 = labW[-p], heights2 = labH[-p]);
        p.x[p] = testPositions[[p]]$x[i];
        p.y[p] = testPositions[[p]]$y[i];
        labelDists[-p, p] = labelDistance(x.ref = p.x[p], y.ref = p.y[p],
           widths.ref = labW[p], heights.ref = labH[p],
           x.test = p.x[-p], y.test = p.y[-p],
           widths.test = labW[-p], heights.test = labH[-p])
        done = all(labelDists > 0)
        ip = ip + 1;
      }
      if (FALSE)
      {
         sizeGrWindow(12, 6);
         par(mfrow = c(1,2))
         par(cex = 1)
         plot(x, y, cex = pt.cex)
         text(initialPositions[, 1], initialPositions[, 2], labels[plotLabels], cex = cex, adj = c(0.5, 0.5));
         plot(x, y, cex = pt.cex)
         text(p.x, p.y, labels[plotLabels], cex = cex, adj = c(0.5, 0.5));
       }
      run = run + 1;
    }
    if (!done)
    {
      # need to remove one of the offending labels.
      candidates = which(colSums(labelDists<0) > 0);
      negativeSpace = colSums( (labelDists < 0) * labelDists);
      removeOrder = order(labelPreference[testPoints], negativeSpace, -labW);
      remove = removeOrder[1];
      remove.origIndex = testPoints[remove];
      testPoints = testPoints[-remove];
      plotLabels[remove.origIndex] = FALSE;
      labW = labW[-remove];
      labH = labH[-remove];
      labelDists = labelDists[-remove, -remove];
      p.x = p.x[-remove];
      p.y = p.y[-remove];
      testPositions = testPositions[-remove];
    }
  }

  if (par("xlog")) p.x = 10^p.x;
  if (par("ylog")) p.y = 10^p.y;

  if (doPlot && any(plotLabels))
    text(p.x, p.y, labels[plotLabels], cex = cex, xpd = xpd, adj = c(0.5, 0.5), font = font, ...) 

  # For compatibility with older labelPoints, the return value should contain x, y, label at top level.
  invisible(list(x = p.x, y = p.y, label = labels[plotLabels], plotLabels = plotLabels, 
                 initial = list(x = initialPositions[, 1], y = initialPositions[, 2],
                                labels = labels[plotLabels.initial], plotLabels = plotLabels.initial)));
}

