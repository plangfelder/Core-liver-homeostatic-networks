#--------------------------------------------------------------------------
#
# .reverseRows = function(Matrix)
#
#--------------------------------------------------------------------------
#


.reverseRows = function(mat)
{
  ind = seq(from=dim(mat)[1], to=1, by=-1);
  mat[ind,, drop = FALSE];
}

.extend = function(x, n)
{
  nRep = ceiling(n/length(x));
  rep(x, nRep)[1:n];
}

# Adapt a numeric index to a subset
# Aim: if 'index' is a numeric index of special entries of a vector,
#    create a new index that references 'subset' elements of the vector  
.restrictIndex = function(index, subset)
{
  out = match(index, subset);
  out[!is.na(out)];
}

  
#--------------------------------------------------------------------------
#
# heatmap.wg
#
#--------------------------------------------------------------------------

.findFnc = function(name)
{
  if (exists(name, mode = "function")) match.fun(name) else
    eval(parse(text = spaste("WGCNA:::", name)))
}

wgHeatmap = function (
  # Content of heatmap
  Matrix, 
  colorMatrix = NULL,
  colors = NULL, 
  zlim = NULL,
  invertColors = FALSE, 
  naColor = "grey",
  textMatrix = NULL, cex.text = NULL, 
  textAdj = c(0.5, 0.5),

  # Labeling of rows and columns
  colLabels.bottom = NULL, 
  colLabels.top = NULL, 
  rowLabels.left = NULL, 
  rowLabels.right = NULL, 
  colColors.bottom = NULL, 
  colColors.top = NULL, 
  rowColors.left = NULL,
  rowColors.right = NULL,
  colLabelsAngle.bottom = 45,
  colLabelsAngle.top = 45,
  colLabelsAdj.bottom = NULL,
  colLabelsAdj.top = NULL,
  rowLabelsAdj.left = 1,
  rowLabelsAdj.right = 0,

  colColorText.top = NULL,
  cex.colColorText.top = 1,
  colColorTextRows.top = NULL,
  colColorRowWidths.top = NULL,
  colColorTextAlignment.top = "center",
  
  colColorText.bottom = NULL,
  cex.colColorText.bottom = 1,
  colColorTextRows.bottom = NULL,
  colColorRowWidths.bottom = NULL,
  colColorTextAlignment.bottom = "center",

  colColorWidth1 = 0.2, ## This is in inches
  rowColorWidth1 = colColorWidth1,

  textLabelGap = 0.04,  ## In inches
  mapToColorGap = 0.04, ## in inches

  leftMarginFreeSpace = textLabelGap, ## in inches; labels wil not be positioned in this space
  rightMarginFreeSpace = textLabelGap, ## in inches

  colColorLabels.top = NULL,
  colColorLabels.bottom = NULL,
  colColorLabelsPosition.top = "left",
  colColorLabelsPosition.bottom = "right",
  rowColorLabels.left = NULL,
  rowColorLabels.right = NULL,
  colColorLabelsAngle.top = 0,
  colColorLabelsAngle.bottom = 45,
  rowColorLabelsPosition.left = "bottom",
  rowColorLabelsPosition.right = "bottom",
  rowColorLabelsAngle.left = 45,
  rowColorLabelsAngle.right = 45,
  cex.lab = NULL, 
  cex.lab.col.bottom = cex.lab,
  cex.lab.col.top = cex.lab,
  cex.lab.row.left = cex.lab,
  cex.lab.row.right = cex.lab,
  colors.lab.col.top = 1,
  colors.lab.col.bottom = 1,
  colors.lab.row.left = 1,
  colors.lab.row.right = 1,
  font.lab.col.top = 1,
  font.lab.col.bottom = 1,
  font.lab.row.left = 1,
  font.lab.row.right = 1,
  bg.lab.col = NULL,
  bg.lab.row = NULL,

  # Legend for the main heatmap
  plotLegend = TRUE, 
  keepLegendSpace = plotLegend,
  legendLabel = "",
  cex.legendLabel = 1,

  # Separator line specification                   
  verticalSeparator.x = NULL,
  verticalSeparator.col = 1,  
  verticalSeparator.lty = 1,
  verticalSeparator.lwd = 1,
  verticalSeparator.ext = 0,
  verticalSeparator.interval = 0,

  horizontalSeparator.y = NULL,
  horizontalSeparator.col = 1,  
  horizontalSeparator.lty = 1,
  horizontalSeparator.lwd = 1,
  horizontalSeparator.ext = 0,
  horizontalSeparator.interval = 0,

  # optional restrictions on which rows and columns to actually show
  showRows = NULL,
  showCols = NULL,

  # Other arguments...
  ... ) 
{
  textFnc = match.fun("text");
  
  #if (is.null(rowLabels) & (!is.null(colLabels)) & (dim(Matrix)[1]==dim(Matrix)[2])) 
  #  rowLabels = colLabels; 

  nCols = ncol(Matrix);
  nRows = nrow(Matrix);

  if (is.null(showRows)) showRows = c(1:nRows);
  if (is.null(showCols)) showCols = c(1:nCols);

  nShowCols = length(showCols);
  nShowRows = length(showRows);

  if (nShowCols==0) stop("'showCols' is empty.");
  if (nShowRows==0) stop("'showRows' is empty.");

  if (length(colColors.bottom) > 0)
  {
    colColors.bottom = as.matrix(colColors.bottom);
    if (!is.null(colColorText.bottom))
    {
      colColorText.bottom = as.matrix(colColorText.bottom);
      nColColorText.bottom = ncol(colColorText.bottom);
    } else nColColorText.bottom = 0;
    if (!is.null(colColorRowWidths.bottom))
    {
      if (length(colColorRowWidths.bottom) > nrow(colColors.bottom) + nColColorText.bottom)
      stop("When 'colColorRowWidths.bottom' are given, their length must correspond to number of columns\n",
           " of 'colColors.bottom' and 'colColorText.bottom'.");
      nColColorRows.bottom = sum(colColorRowWidths.bottom);
    } else 
      nColColorRows.bottom = ncol(colColors.bottom) + nColColorText.bottom;
  } else
    nColColorRows.bottom = 0

  if (length(colColors.top) > 0)
  {
    colColors.top = as.matrix(colColors.top);
    if (!is.null(colColorText.top))
    {
      colColorText.top = as.matrix(colColorText.top);
      nColColorText.top = ncol(colColorText.top);
    } else nColColorText.top = 0;
    if (!is.null(colColorRowWidths.top))
    {
      if (length(colColorRowWidths.top) > nrow(colColors.top) + nColColorText.top)
      stop("When 'colColorRowWidths.top' are given, their length must correspond to number of columns\n",
           " of 'colColors.top' and 'colColorText.top'.");
      nColColorRows.top = sum(colColorRowWidths.top);
    } else
      nColColorRows.top = ncol(colColors.top) + nColColorText.top;
  } else
    nColColorRows.top = 0

  if (length(rowColors.left) > 0)
  {
    rowColors.left = as.matrix(rowColors.left);
    nRowColorRows.left = ncol(rowColors.left);
  } else 
    nRowColorRows.left = 0;

  if (length(rowColors.right) > 0)
  {
    rowColors.right = as.matrix(rowColors.right);
    nRowColorRows.right = ncol(rowColors.right);
  } else 
    nRowColorRows.right = 0;

  haveColLabels.top = length(colLabels.top) > 0;
  haveColLabels.bottom = length(colLabels.bottom) > 0;
  haveRowLabels.left = length(rowLabels.left) > 0;
  haveRowLabels.right= length(rowLabels.right) > 0;

  if (length(colLabels.top)==0) colLabels.top = rep("", nCols);
  if (length(colLabels.bottom)==0) colLabels.bottom = rep("", nCols);
  if (length(rowLabels.left)==0) rowLabels.left = rep("", nRows);
  if (length(rowLabels.right)==0) rowLabels.right= rep("", nRows);

  if (length(colLabels.top)!=nCols) 
    stop("Length of 'colLabels.top' must equal the number of columns in 'Matrix.'");
  if (length(colLabels.bottom)!=nCols) 
    stop("Length of 'colLabels.bottom' must equal the number of columns in 'Matrix.'");

  if (length(rowLabels.left)!=nRows)
    stop("Length of 'rowLabels.left' must equal the number of rows in 'Matrix.'");
  if (length(rowLabels.right)!=nRows)
    stop("Length of 'rowLabels.right' must equal the number of rows in 'Matrix.'");

  colLabels.top.show = colLabels.top[showCols];
  colLabels.bottom.show = colLabels.bottom[showCols];
  rowLabels.left.show = rowLabels.left[showRows];
  rowLabels.right.show = rowLabels.right[showRows];

  
  if (nColColorRows.top > 0) {
     colColors.top.show = colColors.top[showCols, , drop = FALSE] 
  } else colColors.top.show = colColors.top

  if (nColColorRows.bottom > 0) {
     colColors.bottom.show = colColors.bottom[showCols, , drop = FALSE] 
  } else colColors.bottom.show  = colColors.bottom;

  if (nRowColorRows.left > 0) {
    rowColors.left.show = rowColors.left[showRows, , drop = FALSE]
  } else rowColors.left.show = rowColors.left;

  if (nRowColorRows.right > 0) {
    rowColors.right.show = rowColors.right[showRows, , drop = FALSE]
  } else rowColors.right.show = rowColors.right;

  if (length(rowColorLabels.left)==0) rowColorLabels.left = rep("", nRowColorRows.left);
  if (length(rowColorLabels.right)==0) rowColorLabels.right = rep("", nRowColorRows.right);

  if (length(colColorLabels.top)==0) colColorLabels.top = rep("", nColColorRows.top);
  if (length(colColorLabels.bottom)==0) colColorLabels.bottom = rep("", nColColorRows.bottom);

  if (is.null(colors)) colors = heat.colors(30);
  if (invertColors) colors = rev(colors);

  labPos = .findFnc(".heatmapWithLegend")(Matrix[showRows, showCols, drop = FALSE], 
                 signed = FALSE, colorMatrix = if (is.null(colorMatrix)) NULL else colorMatrix[showRows, showCols], 
                 colors = colors, naColor = naColor, cex.legend = cex.lab, 
                 plotLegend = plotLegend,  keepLegendSpace = keepLegendSpace, zlim = zlim, legendLabel = legendLabel,
                 cex.legendLabel = cex.legendLabel, ...)
  plotbox = labPos$box;
  xmin = plotbox[1]; xmax = plotbox[2]; ymin = plotbox[3]; yrange = plotbox[4]-ymin;
  ymax = plotbox[4]; xrange = xmax - xmin;
  # The positions below are for showCols/showRows-restricted data
  xLeft = labPos$xLeft;
  xRight = labPos$xRight;
  yTop = labPos$yTop;
  yBot = labPos$yBot;

  xspacing = labPos$xMid[2] - labPos$xMid[1];
  yspacing = abs(labPos$yMid[2] - labPos$yMid[1]);

  figureDims = par("pin");
  figureBox = par("usr");
  figXrange = figureBox[2] - figureBox[1];
  figYrange = figureBox[4] - figureBox[3];
  ratio = figureDims[1]/figureDims[2] * figYrange/figXrange;

  colColorWidth1.usr = colColorWidth1 * figYrange/figureDims[2];
  rowColorWidth1.usr = rowColorWidth1 * figXrange/figureDims[1];

  colColorWidth.top.all = colColorWidth1.usr * nColColorRows.top;
  colColorWidth.bottom.all = colColorWidth1.usr * nColColorRows.bottom;
  rowColorWidth.left.all = rowColorWidth1.usr * nRowColorRows.left;
  rowColorWidth.right.all = rowColorWidth1.usr * nRowColorRows.right;

  # Additional angle-dependent offsets for x axis labels
  textOffsetY.top = strheight("M") * cos(colLabelsAngle.top/180 * pi);
  textOffsetY.bottom = strheight("M") * cos(colLabelsAngle.bottom/180 * pi);

  textLabelGap.x = textLabelGap * figXrange/figureDims[1];
  textLabelGap.y = textLabelGap * figYrange/figureDims[2];

  leftMarginFreeSpace.usr = leftMarginFreeSpace * figXrange/figureDims[1];
  rightMarginFreeSpace.usr = rightMarginFreeSpace * figXrange/figureDims[1];

  textLabelGap.y.top = textLabelGap.y + textOffsetY.top;
  textLabelGap.y.bottom = textLabelGap.y + textOffsetY.bottom;

  mapToColorGap.x = mapToColorGap * figXrange/figureDims[1];
  mapToColorGap.y = mapToColorGap * figYrange/figureDims[2];

  mapToColorGap.x.left = mapToColorGap.x * (nRowColorRows.left > 0)
  mapToColorGap.x.right = mapToColorGap.x * (nRowColorRows.right > 0)

  mapToColorGap.y.top = mapToColorGap.y * (nColColorRows.top > 0);
  mapToColorGap.y.bottom = mapToColorGap.y * (nColColorRows.bottom > 0);

  # Create the background for column and row labels.
  extension.left = par("mai")[2] * # left margin width in inches
                      figXrange/figureDims[1];
  extension.right = par("mai")[4] * # right margin width in inches
                      figXrange/figureDims[1];
  extension.bottom = par("mai")[1] * figYrange/figureDims[2] - colColorWidth.bottom.all - textLabelGap.y.bottom - 
                        mapToColorGap.y.bottom;   
  extension.top = par("mai")[3]  * figYrange/figureDims[2] - colColorWidth.top.all - textLabelGap.y.top - mapToColorGap.y.top;

  if (!is.null(bg.lab.col))
  {
    bg.lab.col = .extend(bg.lab.col, nCols)[showCols];
    if (haveColLabels.top)
    {
      angle = colLabelsAngle.top/180*pi;
      ext.x = extension.top * 1/tan(angle)/ratio
      ext.y = extension.top * sign(sin(angle));
      sign = -1; y0 = ymax + colColorWidth.top.all + 0.5 * textLabelGap.y.top + mapToColorGap.y;
      offset = 0.5*textLabelGap.y.top;
      polyLst = list(x = list(), y = list(), col = bg.lab.col);
      for (cc in 1:nShowCols)
      {
         polyLst$x = c(polyLst$x, list(c(xLeft[cc], xLeft[cc], xLeft[cc] + ext.x, xRight[cc] + ext.x, xRight[cc], xRight[cc])));
         polyLst$y = c(polyLst$y, list(c(y0, y0-sign*offset, y0-sign*offset - ext.y, y0-sign*offset - ext.y,
                       y0-sign*offset, y0)));
         polygon(x = c(xLeft[cc], xLeft[cc], xLeft[cc] + ext.x, xRight[cc] + ext.x, xRight[cc], xRight[cc]),
                 y = c(y0, y0-sign*offset, y0-sign*offset - ext.y, y0-sign*offset - ext.y, 
                       y0-sign*offset, y0), 
                 border = bg.lab.col[cc], col = bg.lab.col[cc], xpd = TRUE);
      }
      labPos$colLabelBgPolygons.top = polyLst;
    } 
    if (haveColLabels.bottom)
    { 
      angle = colLabelsAngle.bottom/180*pi;
      ext.x = extension.bottom * 1/tan(angle)/ratio
      ext.y = extension.bottom * sign(sin(angle));
      sign = 1; y0 = ymin - colColorWidth.bottom.all - mapToColorGap.y - 0.5*textLabelGap.y.bottom;
      offset = 0.5 * textLabelGap.y.bottom 
      polyLst = list(x = list(), y = list(), col = bg.lab.col);
      for (cc in 1:nShowCols)
      {
         polygon(x = c(xLeft[cc], xLeft[cc], xLeft[cc] - ext.x, xRight[cc] - ext.x, xRight[cc], xRight[cc]),
                 y = c(y0, y0-sign*offset, y0-sign*offset - ext.y, y0-sign*offset - ext.y,
                       y0-sign*offset, y0), 
                 border = bg.lab.col[cc], col = bg.lab.col[cc], xpd = TRUE);
         polyLst$x = c(polyLst$x, list(c(xLeft[cc], xLeft[cc], xLeft[cc] - ext.x, xRight[cc] - ext.x, xRight[cc], xRight[cc])));
         polyLst$y = c(polyLst$y, list(c(y0, y0-sign*offset, y0-sign*offset - ext.y, y0-sign*offset - ext.y,
                       y0-sign*offset, y0)));
      }
      labPos$colLabelBgPolygons.bottom = polyLst;
    } 
  }

  if (!is.null(bg.lab.row))
  {
    bg.lab.row = .extend(bg.lab.row, nRows)
   # reverseRows = TRUE;
   # if (reverseRows) bg.lab.row = rev(bg.lab.row);
    bg.lab.row = rev(bg.lab.row);
    bg.lab.row = bg.lab.row[showRows];
    if (haveRowLabels.left)
    {
      xl = xmin-extension.left;
      xr = xmin - rowColorWidth.left.all - mapToColorGap.x.left # -0.5*textLabelGap.x 
      for (r in 1:nShowRows)
        rect(xl, yBot[r], xr, yTop[r], col = bg.lab.row[r], border = bg.lab.row[r], xpd = TRUE);
    } 
    if (haveRowLabels.right)
    {
      xl = xmax  + rowColorWidth.right.all + mapToColorGap.x.right # + 0.5*textLabelGap.x 
      xr = xmax + extension.right;
      for (r in 1:nShowRows)
        rect(xl, yBot[r], xr, yTop[r], col = bg.lab.row[r], border = bg.lab.row[r], xpd = TRUE);
    }
  }

  # Column labels and colors
  colors.lab.col.top = .extend(colors.lab.col.top, nCols)[showCols];
  colors.lab.col.bottom = .extend(colors.lab.col.bottom, nCols)[showCols];
  font.lab.col.top = .extend(font.lab.col.top, nCols)[showCols];
  font.lab.col.bottom = .extend(font.lab.col.bottom, nCols)[showCols];

  # Plot top colors and top labels
  if (nColColorRows.top > 0)
  {
    baseY = ymax + colColorWidth.top.all + mapToColorGap.y.top;
    labPos$topColorRectangles = .findFnc(".plotOrderedColorSubplot")(order = showCols,
       colors = colColors.top, 
       rowLabels = colColorLabels.top,
       rowText = colColorText.top,
       cex.rowText = cex.colColorText.top,
       textPositions = colColorTextRows.top,
       rowWidths = colColorRowWidths.top,
       rowTextAlignment = colColorTextAlignment.top,
       horizontal = TRUE,
       rowLabelsPosition = colColorLabelsPosition.top,
       rowLabelsAngle = colColorLabelsAngle.top,
       cex.rowLabels = cex.lab.col.top,
       plotBox = c(xmin, xmax, baseY, baseY - colColorWidth.top.all),
       align = "edge", limExpansionFactor.x = 0, limExpansionFactor.y = 0,
       checkOrderLength = FALSE);
  }

  if (nColColorRows.bottom > 0)
  {
    baseY = ymin - colColorWidth.bottom.all - mapToColorGap.y.bottom;
    labPos$bottomColorRectangles = .findFnc(".plotOrderedColorSubplot")(order = showCols,
       colors = colColors.bottom,
       rowLabels = colColorLabels.bottom,
       rowText = colColorText.bottom,
       cex.rowText = cex.colColorText.bottom,
       textPositions = colColorTextRows.bottom,
       rowWidths = colColorRowWidths.bottom,
       rowTextAlignment = colColorTextAlignment.bottom,
       horizontal = TRUE,
       rowLabelsPosition = colColorLabelsPosition.bottom,
       rowLabelsAngle = colColorLabelsAngle.bottom,
       cex.rowLabels = cex.lab.col.bottom,
       plotBox = c(xmin, xmax, baseY, baseY + colColorWidth.bottom.all),
       align = "edge", limExpansionFactor.x = 0, limExpansionFactor.y = 0,
       checkOrderLength = FALSE);
  }

  if (is.null(colLabelsAdj.top)) colLabelsAdj.top = 0;
  if (is.null(colLabelsAdj.bottom)) colLabelsAdj.bottom = 1;

  if (haveColLabels.top)
  {
    baseY = ymax + colColorWidth.top.all + mapToColorGap.y.top + textLabelGap.y.top;
    mapply(textFnc, x = labPos$xMid,
         labels = colLabels.top.show,
         col = colors.lab.col.top,
         font = font.lab.col.top,
         MoreArgs = list( adj = colLabelsAdj.top, 
           y = baseY,
           xpd = TRUE, srt = colLabelsAngle.top, cex = cex.lab.col.top));
    labPos$colLabelsPos.top = list(x = labPos$xMid, y = baseY, adj = colLabelsAdj.top, srt = colLabelsAngle.top);
  }

  if (haveColLabels.bottom)
  {
    baseY = ymin - colColorWidth.bottom.all - mapToColorGap.y.bottom - textLabelGap.y.bottom;
    mapply(textFnc, x = labPos$xMid,
         labels = colLabels.bottom.show,
         col = colors.lab.col.bottom,
         font = font.lab.col.bottom,
         MoreArgs = list( adj = colLabelsAdj.bottom,
           y = baseY, 
           xpd = TRUE, srt = colLabelsAngle.bottom, cex = cex.lab.col.bottom));
    labPos$colLabelsPos.bottom = list(x = labPos$xMid, y = baseY, adj = colLabelsAdj.bottom, srt = colLabelsAngle.bottom);
  }

  #-----------------------------------------------------------------------------------------------------
  # Row colors and labels
  rowLabelsAdj.left = .extend(rowLabelsAdj.left, nRows)[showRows]
  rowLabelsAdj.right = .extend(rowLabelsAdj.right, nRows)[showRows]

  colors.lab.row.left = .extend(colors.lab.row.left, nRows)[showRows];
  font.lab.row.left = .extend(font.lab.row.left, nRows)[showRows];
  colors.lab.row.right = .extend(colors.lab.row.right, nRows)[showRows];
  font.lab.row.right = .extend(font.lab.row.right, nRows)[showRows];

  marginWidth.left = par("mai")[2] / par("pin")[1] * xrange
  marginWidth.right = par("mai")[4] / par("pin")[1] * xrange

  # Left colors and labels
  xSpaceForRowLabels = marginWidth.left - leftMarginFreeSpace.usr - mapToColorGap.x.left - textLabelGap.x - 
                            rowColorWidth.left.all;
  xPosOfRowLabels.relative = xSpaceForRowLabels * (1-rowLabelsAdj.left)

  xr = xmin - mapToColorGap.x.left
  xl = xr-rowColorWidth.left.all; 
  xtext = xl - textLabelGap.x - xPosOfRowLabels.relative;
  if (nRowColorRows.left > 0)
  {
    .findFnc(".plotOrderedColorSubplot")(order = rev(1:nShowRows),
       colors = rowColors.left.show,
       rowLabels = rowColorLabels.left,
       horizontal = FALSE,
       rowLabelsPosition = if (rowColorLabelsPosition.left=="bottom") "left" else "right",
       cex.rowLabels = cex.lab.row.left, rowLabelsAngle = rowColorLabelsAngle.left,
       plotBox = c(xl, xr, ymin, ymax),
       align = "edge", limExpansionFactor.x = 0, limExpansionFactor.y = 0);
  }

  mapply(textFnc, y = labPos$yMid, labels = rowLabels.left.show,
             x = xtext,
             adj = lapply(rowLabelsAdj.left, c, 0.5),
             col = colors.lab.row.left,
             font = font.lab.row.left,
          MoreArgs = list(srt = 0, xpd = TRUE, cex = cex.lab.row.left));

  # Right colors and labels

  xSpaceForRowLabels = marginWidth.right- rightMarginFreeSpace.usr - mapToColorGap.x.right - textLabelGap.x - 
                            rowColorWidth.right.all;
  xPosOfRowLabels.relative = xSpaceForRowLabels * rowLabelsAdj.right
  xl = xmax + mapToColorGap.x.right;
  xr = xmax + rowColorWidth.right.all;
  xtext = xr + textLabelGap.x + xPosOfRowLabels.relative;

  if (nRowColorRows.right > 0)
  {
    .findFnc(".plotOrderedColorSubplot")(order = rev(1:nShowRows),
       colors = rowColors.right.show,
       rowLabels = rowColorLabels.right,
       horizontal = FALSE,
       rowLabelsPosition = if (rowColorLabelsPosition.right=="bottom") "right" else "right",
       cex.rowLabels = cex.lab.row.right, rowLabelsAngle = rowColorLabelsAngle.right,
       plotBox = c(xl, xr, ymin, ymax),
       align = "edge", limExpansionFactor.x = 0, limExpansionFactor.y = 0);
  }

  mapply(textFnc, y = labPos$yMid, labels = rowLabels.right.show,
             x = xtext,
             adj = lapply(rowLabelsAdj.right, c, 0.5),
             col = colors.lab.row.right,
             font = font.lab.row.right,
          MoreArgs = list(srt = 0, xpd = TRUE, cex = cex.lab.row.right));

  # Draw separator lines, if requested

  showCols.ext = c(if (1 %in% showCols) 0 else NULL, showCols);
  showCols.shift = if (0 %in% showCols.ext) 1 else 0;

  if (length(verticalSeparator.x) > 0)
  {
    if (any(verticalSeparator.x < 0 | verticalSeparator.x > nCols))
      stop("If given. 'verticalSeparator.x' must all be between 0 and the number of columns.");
    colSepShowIndex = which(verticalSeparator.x %in% showCols.ext);
    verticalSeparator.x.show = .restrictIndex(verticalSeparator.x, showCols.ext)-showCols.shift;
  } else if (verticalSeparator.interval > 0)
  {
    verticalSeparator.x.show = verticalSeparator.x = 
           seq(from = verticalSeparator.interval, by = verticalSeparator.interval,
                                    length.out = floor(length(showCols)/verticalSeparator.interval));
    colSepShowIndex = 1:length(verticalSeparator.x);
  } else 
    verticalSeparator.x.show = NULL;

  if (length(verticalSeparator.x.show) > 0)
  {
    nLines = length(verticalSeparator.x);
    vs.col = .extend(verticalSeparator.col, nLines)[colSepShowIndex];
    vs.lty = .extend(verticalSeparator.lty, nLines)[colSepShowIndex];
    vs.lwd = .extend(verticalSeparator.lwd, nLines)[colSepShowIndex];
    vs.ext = .extend(verticalSeparator.ext, nLines)[colSepShowIndex];

    x.lines = ifelse(verticalSeparator.x.show>0, labPos$xRight[verticalSeparator.x.show], labPos$xLeft[1]);
    nLines.show = length(verticalSeparator.x.show);

    for (l in 1:nLines.show)
    {
      lines.x = numeric(0);
      lines.y = numeric(0);
      offset.top = textLabelGap.y.top + (nColColorRows.top > 0) * mapToColorGap.y + colColorWidth.top.all
      offset.bottom = textLabelGap.y.bottom + (nColColorRows.bottom > 0) * mapToColorGap.y + colColorWidth.bottom.all
      if (haveColLabels.top)
      {
        angle = colLabelsAngle.top/180*pi;
        if (angle==0) angle = pi/2;
        ext.x = extension.top* 1/tan(angle)/ratio;
        ext.y = extension.top * sign(sin(angle))
        lines.x = c(lines.x, vs.ext[l] * ext.x, 0); 
        lines.y = c(lines.y, ymax + offset.top + vs.ext[l] * ext.y, ymax + offset.top);
      }
      lines.x = c(lines.x, 0, 0);
      lines.y = c(lines.y, ymax, ymin);
      if (haveColLabels.bottom)
      {
        angle = colLabelsAngle.bottom/180*pi;
        if (angle==0) angle = pi/2;
        ext.x = -extension.bottom* 1/tan(angle)/ratio;
        ext.y = extension.bottom * sign(sin(angle))
        lines.x = c(lines.x, 0, vs.ext[l] * ext.x); 
        lines.y = c(lines.y, ymin - offset.bottom, ymin - offset.bottom - vs.ext[l] * ext.y);
      }
      lines(lines.x + x.lines[l], lines.y, col = vs.col[l], lty = vs.lty[l], lwd = vs.lwd[l], xpd = TRUE);
    }
  }

  showRows.ext = c(if (1 %in% showRows) 0 else NULL, showRows);
  showRows.shift = if (0 %in% showRows.ext) 1 else 0;

  if (length(horizontalSeparator.y) >0)
  {
    if (any(horizontalSeparator.y < 0 | horizontalSeparator.y > nRows))
      stop("If given. 'horizontalSeparator.y' must all be between 0 and the number of rows.");
    rowSepShowIndex = which( horizontalSeparator.y %in% showRows.ext);
    horizontalSeparator.y.show = .restrictIndex(horizontalSeparator.y, showRows.ext)-showRows.shift;
  } else if (horizontalSeparator.interval > 0)
  {
    horizontalSeparator.y.show = horizontalSeparator.y = 
            seq(from = horizontalSeparator.interval, by = horizontalSeparator.interval,
                                    length.out = floor(length(showRows)/horizontalSeparator.interval));
    rowSepShowIndex = 1:length(horizontalSeparator.y);
  } else 
    horizontalSeparator.y.show = NULL;
  
  if (length(horizontalSeparator.y.show) > 0)
  {
    reverseRows = TRUE;
    if (reverseRows) 
    {
      horizontalSeparator.y.show = nShowRows - horizontalSeparator.y.show+1;
      y.lines = ifelse( horizontalSeparator.y.show <=nShowRows, 
                               labPos$yBot[horizontalSeparator.y.show], labPos$yTop[nShowRows]);
    } else {
      y.lines = ifelse( horizontalSeparator.y.show > 0, labPos$yBot[horizontalSeparator.y.show], labPos$yTop[1]);
    }
    nLines = length(horizontalSeparator.y);
    vs.col = .extend(horizontalSeparator.col, nLines)[rowSepShowIndex];
    vs.lty = .extend(horizontalSeparator.lty, nLines)[rowSepShowIndex];
    vs.lwd = .extend(horizontalSeparator.lwd, nLines)[rowSepShowIndex];
    vs.ext = .extend(horizontalSeparator.ext, nLines)[rowSepShowIndex];
    nLines.show = length(horizontalSeparator.y.show);
    for (l in 1:nLines.show)
    {
      xl = xmin;
      xr = xmax; 
      if (haveRowLabels.left)
        xl = xmin - vs.ext[l]*extension.left 
      if (haveRowLabels.right)
        xr = xmax + vs.ext[l]*extension.right;
      lines(c(xl, xr), rep(y.lines[l], 2), 
            col = vs.col[l], lty = vs.lty[l], lwd = vs.lwd[l], xpd = TRUE);
    }
  }

  if (!is.null(textMatrix))
  {
    if (is.null(cex.text)) cex.text = par("cex");
    if (is.null(dim(textMatrix)))
      if (length(textMatrix)==prod(dim(Matrix))) dim(textMatrix)=dim(Matrix);
    if (!isTRUE(all.equal(dim(textMatrix), dim(Matrix))))
      stop("labeledHeatmap: textMatrix was given, but has dimensions incompatible with Matrix.");
    for (rw in 1:nShowRows)
      for (cl in 1:nShowCols)
      {
        text(labPos$xMid[cl], labPos$yMid[rw],
             as.character(textMatrix[showRows[rw],showCols[cl]]), xpd = TRUE, cex = cex.text, adj = textAdj);
      }
  }
  axis(1, labels = FALSE, tick = FALSE)
  axis(2, labels = FALSE, tick = FALSE)
  axis(3, labels = FALSE, tick = FALSE)
  axis(4, labels = FALSE, tick = FALSE)
  invisible(labPos)
}

#Some basic tests.

if (FALSE)
{
   set.seed(1)
   n = 10;
   m = 13
   Matrix = matrix(rnorm(n*m), n, m);
   txt = round(Matrix, 2);

   max  = max(abs(Matrix));
   sizeGrWindow(12, 9);
   par(mar = c(7,7,8,8));
   wgHeatmap(Matrix = Matrix,
      rowLabels.left = spaste("LRow ", 1:n),
      rowLabels.right = spaste("RRow ", 1:n),
      rowColors.left = cbind(standardColors(n), numbers2colors(Matrix)[, 1]),
      rowColors.right = cbind(standardColors(n), numbers2colors(Matrix)[, 2], standardColors(n)),
      rowColorLabels.left = c("Color", "Number"),
      rowColorLabels.right = c("RColor", "RNumber", "RColor2"),
      rowColorLabelsAngle.left = 45,
      rowColorLabelsAngle.right = 60,
      colLabels.bottom = spaste("Col ", 1:m),
      colLabels.top = spaste("TCol ", 1:m),
      colColorLabelsAngle.bottom = 0,
      colColorLabelsAngle.top = 0,
      colColors.bottom = cbind(standardColors(m+n)[(n+1):(n+m)], numbers2colors(Matrix)[1, ]),
      colColors.top = cbind(standardColors(m+n)[(n+1):(n+m)], t(numbers2colors(Matrix)[c(1:2), ])),
      colColorLabels.bottom = c("Color.col", "Number.col"),
      colColorLabels.top = c("Color.col", "Number.col", "Number.col2"),
      colLabelsAngle.bottom = 45,
      colLabelsAngle.top = 90,
      
      main = "Test 1.", zlim = c(-max, max), colors = blueWhiteRed(100),
      );
}

#===================================================================================================
#
# multi-page labeled heatmap
#
#===================================================================================================


wgHeatmap.multiPage = function(
   Matrix,
   # Paging options
   rowsPerPage = NULL, maxRowsPerPage = 20,
   colsPerPage = NULL, maxColsPerPage = 10,
   addPageNumberToMain = TRUE,
   mainSep = "\n",
   mainExtStart = "(page ", mainExtEnd = ")",
   main = "",

   ...)
{

  nr = nrow(Matrix);
  nc = ncol(Matrix);

  if (is.null(rowsPerPage))
  {
    nPages.rows = ceiling(nr/maxRowsPerPage);
    rowsPerPage = allocateJobs(nr, nPages.rows);
  } else 
    nPages.rows = length(rowsPerPage);

  if (is.null(colsPerPage))
  {
    nPages.cols = ceiling(nc/maxColsPerPage);
    colsPerPage = allocateJobs(nc, nPages.cols);
  } else 
    nPages.cols = length(colsPerPage);

  page = 1;
  multiPage = (nPages.cols > 1 | nPages.rows > 1)

  for (page.col in 1:nPages.cols) for (page.row in 1:nPages.rows)
  {
    rows = rowsPerPage[[page.row]];
    cols = colsPerPage[[page.col]];
    main.1 = main;
    if (addPageNumberToMain & multiPage) main.1 = spaste(main, mainSep, mainExtStart, page, mainExtEnd);
    wgHeatmap(Matrix = Matrix,
                   main = main.1, 
                   showRows = rows, showCols = cols,
                   ...);
    page = page + 1;
  }
}
                   


