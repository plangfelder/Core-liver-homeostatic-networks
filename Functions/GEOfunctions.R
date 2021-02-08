nUnique = function(x) { length(unique(x))}

library("RCurl");

removeTrailingSlash = function(path)
{
  n = length(path);
  for (i in 1:n)
  {
    nc = nchar(path[i]);
    if (substring(path[i], nc, nc)=="/") path[i] = substring(path[i], 1, nc-1);
  }
  path;
}

createAndCheckDir = function(path)
{
  x = dir.create(path, recursive = TRUE);
  x = file.info(path);
  if (length(x)==1) stop(spaste("Directory ", path, " could not be created."));
  if (!x$isdir) stop(spaste("File ", path, " exists and is not a directory."));
  return(0);
}

saveData = function(x, fileBase, useBZip2 = TRUE, useGZip = FALSE, message = "", quote = FALSE)
{
  if (useBZip2) {
    fileBase = spaste(fileBase, ".bz2")
  } else if (useGZip)
    fileBase = spaste(fileBase, ".gz");

  if (!is.null(message))
  {
    message = gsub("%file", fileBase, message, fixed = TRUE);
    printFlush(message);
  }

  file = fileBase; 
  if (useBZip2) file = bzfile(fileBase);
  if (useGZip) file = gzfile(fileBase);

  write.table(dataForTable(x, IDcolName = "ProbeID", transpose = FALSE), 
              file = file, quote = quote, sep = "\t", row.names = FALSE);
}


# Drop columns that have median length longer than threshold characters.

dropLongColumns = function(x, threshold = 100)
{
  medLen = apply( as.matrix(apply(as.matrix(x), 2, nchar)), 2, median);
  drop = replaceMissing(medLen) > threshold;
  x[, !drop, drop = FALSE];
}

dropUninformativeColumns = function(x)
{
  x = as.data.frame(x);
  nlevels = sapply(x, nUnique)
  informative = nlevels>1
  x[informative, drop = FALSE];
}

# This function is adapted from GEOquery:::getAndParseGSEMatrices
# It attempts to download supplementary data for the given GEO identifier and put them into the destdir
# directory. It returns a vector with the names of the downloaded files or NULL if no files were
# downloaded.

downloadSupplementaryData = function (GEO, destdir) 
{
    GEO <- toupper(GEO)
    stub = gsub("\\d{1,3}$", "nnn", GEO, perl = TRUE)
    gdsurl <- "ftp://ftp.ncbi.nlm.nih.gov/geo/series/%s/%s/suppl/"
    a <- try(getURL(sprintf(gdsurl, stub, GEO)));
    if (inherits(a, "try-error"))
    {
      printFlush(spaste("GEO identifier ", GEO, " does not contain supplementary data."));
      return(NULL);
    }
    tmpcon <- textConnection(a, "r")
    b <- read.table(tmpcon)
    close(tmpcon)
    b <- as.character(b[, ncol(b)])
    message(sprintf("Found %d file(s)", length(b)))
    ret <- list()
    for (i in 1:length(b)) {
        message(b[i])
        destfile = file.path(destdir, b[i])
        if (file.exists(destfile)) {
            message(sprintf("File %s already exists", destfile))
        }
        else {
            download.file(sprintf("ftp://ftp.ncbi.nlm.nih.gov/geo/series/%s/%s/suppl/%s", 
                stub, GEO, b[i]), destfile = destfile, mode = "wb", 
                method = getOption("download.file.method.GEOquery"))
        }
    }
    return(b)
}

   

processGEOExpressionSeries = function(identifier, dirBase = "", addIDtoBase = TRUE, rawDir = "GEOdata", 
                                      sampleAnnotDir = "SampleAnnotation",
                                      expressionDir = "Expression/020-AuthorNormalized",
                                      arrayAnnotDir = "Expression/ArrayAnnotation",
                                      sampleAnnotFile = spaste(identifier, "-sampleAnnotation"),
                                      expressionFile = spaste(identifier, "-expression"),
                                      arrayAnnotFile = "arrayAnnotation-%GPLid", 
                                      dropLongColumnsInAA = TRUE,
                                      dropLongThreshold = 100,
                                      useBZip2 = TRUE, useGZip = FALSE, 
                                      quote = FALSE,
                                      getArrayAnnot = TRUE,
                                      keepExpr = TRUE, keepArrayAnnot = TRUE, 
                                      keepSampleAnnot = TRUE, 
                                      getSupplementaryData = FALSE,
                                      supplementaryDataDir = "Expression/010-RawData",
                                      verbose = 1)
{

  dirBase = removeTrailingSlash(dirBase);
  sampleAnnotDir = removeTrailingSlash(sampleAnnotDir);
  expressionDir = removeTrailingSlash(expressionDir);

  if (addIDtoBase) 
    dirBase = spaste(dirBase, "-", identifier);

  if (verbose > 0) printFlush(spaste(
       "\n====================================================================================\n\n",
       "Querying GEO for", identifier, " (", dirBase, ") and parsing data..."));

  rawDirX = file.path(dirBase, rawDir);
  createAndCheckDir(rawDirX);

  data = getGEO(identifier, destdir = rawDirX, GSEMatrix=TRUE, AnnotGPL = FALSE);

  if (length(data)==0) stop("Returned data set is not valid.");

  if (getSupplementaryData) 
  {
    out = try({
    if (verbose > 1) printFlush("Attempting to download supplementary data.");
    suppDirX = file.path(dirBase, supplementaryDataDir);
    createAndCheckDir(suppDirX);
    supplementaryFiles = downloadSupplementaryData(identifier, suppDirX);
    })
    if (inherits(out, "try-error")) supplementaryFiles = NULL;
  } else 
    supplementaryFiles = NULL

  exprDirX = spaste(dirBase, "/", expressionDir);
  createAndCheckDir(exprDirX);
  sampleDirX = spaste(dirBase, "/", sampleAnnotDir);
  createAndCheckDir(sampleDirX);

  nSets = length(data)
  expr = list();
  pheno = list();
  arrayAnnot = list();
  arrayTitle = character(0);


  if (nSets > 1) {
     setSuffix = spaste("-", 1:nSets, ".txt");
     printFlush(paste("The data consists of", nSets, "matrices.\n",  
                      "Will use suffixes -1, -2, etc when saving expression and phenotype data."));
  } else
     setSuffix = ".txt";
  message = NULL;
  for (set in 1:nSets)
  {
    expr[[set]] = exprs(data[[set]]);
    if (verbose)  message = "Saving expression data into file %file";
    saveData(expr[[set]], spaste(exprDirX, "/", expressionFile, setSuffix[set]), useBZip2 = useBZip2, 
             useGZip = useGZip, message = message);

    pheno[[set]] = dropUninformativeColumns(pData(data[[set]]));
    if (verbose)  message = "Saving sample annotation data into file %file";
    saveData(pheno[[set]], spaste(sampleDirX, "/", sampleAnnotFile, setSuffix[set]), useBZip2 = FALSE, 
             useGZip = FALSE, message = message, quote = quote);

    arrID = annotation(data[[set]]);
    if (length(arrID) > 1) 
      printFlush(paste("Something's weird: multiple array annotations:\n", paste(arrID, collapse = ", ")))

    #arrayAnnot

    if (length(arrID) == 1 && getArrayAnnot) { 
      out = try({
       if (verbose >0) printFlush(paste("Queryig GEO for array annotation", arrID, "..."));

       arrayAnnotDirX = spaste(dirBase, "/", arrayAnnotDir);
       createAndCheckDir(arrayAnnotDirX);

       # Attempt to download the annotation file as well
       arrayData = getGEO(arrID, destdir = rawDirX);
       arrayAnnot0 = arrayData@dataTable@table
       if (dropLongColumnsInAA) {
          arrayAnnot[[set]] = dropLongColumns(dropUninformativeColumns(arrayData@dataTable@table), 
                                       threshold = dropLongThreshold)
       } else 
          arrayAnnot[[set]] = dropUninformativeColumns(arrayData@dataTable@table);

       aaFile = gsub("%GPLid", arrID, arrayAnnotFile, fixed = TRUE);
       aaDFile = spaste(arrayAnnotDirX, "/", aaFile, setSuffix[set]);
       if (verbose)  message = "Saving array annotation into file %file";
       saveData(arrayAnnot[[set]], aaDFile, useBZip2 = useBZip2, 
               useGZip = useGZip, message = message);
       arrayTitle[set] = arrayData@header$title;
       });
      if (inherits(out, "try-error")) arrayAnnot[[set]] = NA;
     } else if (verbose > 0)
       printFlush("No array annotation available from GEO.");

     if (!keepExpr) expr[[set]] = NA;
     if (!keepSampleAnnot) pheno[[set]] = NA;
     if (!keepArrayAnnot) arrayAnnot[[set]] = NA;
  }

  if (nSets==1) {
    expr = expr[[1]];
    pheno = pheno[[1]];
    arrayAnnot = if (getArrayAnnot) arrayAnnot[[1]] else NULL;
  }
  
  invisible(list(expr = expr, sampleAnnot = pheno, arrayAnnot = arrayAnnot, arrayTitles = arrayTitle,
                 supplementaryFiles = supplementaryFiles));
}


#=====================================================================================================
#
# Split character string into column name and value
#
#=====================================================================================================

splitToNameAndValue = function(x, split, fixed, convertToNumeric = NULL, NAstrings = "NA",
                               usePart = c("2nd", "last", "allBut1st"), checkNames = TRUE)
{
  if (!fixed && usePart=="allBut1st") 
     warning(.immediate = TRUE,
             "Caution advised: returned values may no be exactly as in input\n",
             "  if 'split' contains special characters.");

  usePart = match.arg(usePart);
  if (inherits(x, "data.frame"))
  {
    out.list = lapply(x, splitToNameAndValue, split = split, fixed = fixed, usePart = usePart, checkNames = checkNames,
                         convertToNumeric = convertToNumeric, NAstrings = NAstrings);
    out = as.data.frame(lapply(out.list, getElement, "values"), check.names = checkNames);
    names(out) = (if (checkNames) make.names else identity) (sapply(out.list, getElement, "name"));
    out;
  } else {
    x.split = strsplit(x, split = split, fixed = fixed)
    use = !(x %in% NAstrings);
    lengths = sapply(x.split, length)
    if (usePart=="2nd") 
    {
       if (any(lengths[use]!=2)) stop("Data cannot be unambiguously split into name and value.");
       useIndex = rep(2, length(x.split));
    } else if (usePart=="allBut1st") {
       x.split = lapply(x.split, function(s) c(s[1], paste(s[-1], collapse = split)))
       useIndex = rep(2, length(x.split));
    } else
       useIndex = sapply(x.split, length);
    useIndex[!use] = 1;
    firstUse = which(use)[1];
    firstEntries = sapply(x.split[use], `[`, 1);
    if (checkNames && any(firstEntries!=firstEntries[1], na.rm = TRUE)) 
    {   
        printFlush("Error in splitToNameAndValue, dropping into browser (place 1).");
        browser();
          stop("The name is not the same for all components.");
    }
    values.char = mapply(`[`, x.split, useIndex);
    if (is.null(convertToNumeric))
    {
      values.num = as.numeric(as.character(values.char))
      bad = is.na(values.num);
      if (all(values.char[bad] %in% NAstrings))
      {
         values = values.num
      } else
         values = values.char;
    } else if (convertToNumeric) {
      values = as.numeric(as.character(values.char));
    } else 
      values = values.char;
      
    list(name = firstEntries[1], values = values);
  }
}

splitToNameAndValue.df = function(...)
{
  lst = splitToNameAndValue(...)
  out = data.frame(V = lst$values);
  names(out) = lst$name;
  out;
}

# Alternate way of fishing out sample characteristics. In some GEO soft matrices, the order of the characteristics is not
# the same for all samples. Instead of going by column, simply collect a set of attributes for all samples; then arrange
# attributes into columns based on names.

splitToNameAndValue.byName = function(x, split, fixed, convertToNumeric = NULL, NAstrings = "NA",
                               usePart = c("2nd", "last", "allBut1st"), checkNames = TRUE)
{
  if (!fixed && usePart=="allBut1st")
     warning(.immediate = TRUE,
             "Caution advised: returned values may no be exactly as in input\n",
             "  if 'split' contains special characters.");

  if (is.atomic(x)) x = data.frame(x);
  usePart = match.arg(usePart);
  attributes0 = lapply(x, function(x1)
  {
    x.split = strsplit(x1, split = split, fixed = fixed)
    use = !(x1 %in% NAstrings);
    lengths = sapply(x.split, length)
    if (usePart=="2nd")
    {
       if (any(lengths[use]!=2)) stop("Data cannot be unambiguously split into name and value.");
       useIndex = rep(2, length(x.split));
    } else if (usePart=="allBut1st") {
       x.split = lapply(x.split, function(s) c(s[1], paste(s[-1], collapse = split)))
       useIndex = rep(2, length(x.split));
    } else
       useIndex = sapply(x.split, length);
    names1 = rep(NA, length(x1));
    values1 = rep(NA, length(x1));
    names1[use] = sapply(x.split[use], `[`, 1);
    values1[use] = mapply(`[`, x.split[use], useIndex[use])
    list(characts = names1, values = values1);
  });
  columns = unique(unlist(lapply(attributes0, getElement, "characts")))
  columns = columns[!is.na(columns)]
  out = as.data.frame(lapply(columns, function(col)
  {
    values1 = do.call(cbind, lapply(attributes0, function(.attr)
      ifelse(.attr$characts==col, .attr$values, NA)));
    values.char = apply(values1, 1, function(x) 
    {
      x = x[!is.na(x)];
      if (length(x)==0) return(NA)
      if (all(x==x[1])) return(x[1]);
      stop("Have at least two different values for characteristic ", col, ": \n", paste(unique(x), collapse = " | "));
    });
    if (is.null(convertToNumeric))
    {
      values.num = as.numeric(as.character(values.char))
      bad = is.na(values.num);
      if (all(values.char[bad] %in% NAstrings))
      {
         values = values.num
      } else
         values = values.char;
    } else if (convertToNumeric) {
      values = as.numeric(as.character(values.char));
    } else
      values = values.char;
    values;
  }));
  setColnames(out, columns);
}

#=========================================================================================================
#
# Process the GEO sample annotation table
#
#=========================================================================================================

processGEOSampleAnnotation = function(sa, processColPattern = "characteristics", split = ": ",
         fixed = TRUE, convertToNumeric = NULL, NAstrings = "NA",
         usePart = c("2nd", "last", "allBut1st"), checkNames = TRUE,
         rearrangeColumns = FALSE,
         forceNumeric = TRUE)
{
  processCols = grep(processColPattern, names(sa));
  if (length(processCols)==0) return(sa);

  xx = (if (rearrangeColumns) splitToNameAndValue.byName else splitToNameAndValue) (sa[processCols], split = split, 
        fixed = fixed, convertToNumeric = convertToNumeric, NAstrings = NAstrings,
        usePart = usePart, checkNames = checkNames);

  sa[processCols] = xx;
  names(sa)[processCols] = names(xx);

  if (forceNumeric)
  {
    character = sapply(sa[processCols], mode)=="character"
    for (c in processCols[which(character)])
    {
      rest = gsub("[-0-9.eE+]", "", sa[[c]]);
      if (all(replaceMissing(rest==rest[1])))
      {
        num = suppressWarnings(as.numeric(gsub("[^-.0-9eE+]", "", sa[[c]])));
        if (!any(is.na(num))) sa[[c]] = num;
      }
    }
  }

  sa;
}
        
        
#=========================================================================================================

parseGSEMatrix = function (fname, AnnotGPL = FALSE, destdir = tempdir()) 
{
    require(Biobase)
    dat <- readLines(fname)
    nseries <- sum(grepl("^!Series_", dat))
    nsamples <- sum(grepl("^!Sample_", dat))
    con <- GEOquery:::fileOpen(fname)
    header <- read.table(con, sep = "\t", header = FALSE, 
        nrows = nseries)
    tmpdat <- read.table(con, sep = "\t", header = FALSE, 
        nrows = nsamples)
    tmptmp <- t(tmpdat)
    sampledat <- rbind(data.frame(), tmptmp[-1, ])
    colnames(sampledat) <- make.unique(sub("!Sample_", "", as.character(tmpdat[, 
        1])))
    readLines(con, 1)
    datamat <- read.delim(con, sep = "\t", header = TRUE, 
        na.strings = c("NA", "null", "NULL", "Null"), comment.char = "");
    close(con)
    tmprownames = datamat[, 1]
    datamat <- as.matrix(datamat[!is.na(tmprownames), -1])
    rownames(datamat) <- tmprownames[!is.na(tmprownames)]
    if (nrow(datamat) == 1) {
        datamat <- datamat[1:(nrow(datamat) - 1), ]
    }
    else {
        datamat <- as.matrix(datamat[1:(nrow(datamat) - 1), ])
    }
    rownames(sampledat) <- colnames(datamat)
    GPL = as.character(sampledat[1, grep("platform_id", colnames(sampledat), 
        ignore.case = TRUE)])
    printFlush(spaste("Getting GPL: ", GPL));
    x = try({
    gpl <- getGEO(GPL, AnnotGPL = AnnotGPL, destdir = destdir)
    vmd <- Columns(gpl)
    dat <- Table(gpl)
    tmpnames = as.character(dat$ID)
    tmpnames[is.na(tmpnames)] = "NA"
    rownames(dat) <- tmpnames
    dat <- dat[match(tolower(rownames(datamat)), tolower(rownames(dat))), ]
    rownames(vmd) <- make.unique(colnames(Table(gpl)))
    colnames(dat) <- rownames(vmd)
    fd <- new("AnnotatedDataFrame", data = dat, varMetadata = vmd)
    if (is.null(nrow(datamat))) {
        tmpnames <- names(datamat)
        rownames(sampledat) <- tmpnames
        datamat = matrix(nrow = 0, ncol = nrow(sampledat))
        colnames(datamat) <- tmpnames
    }
    else {
        rownames(datamat) <- rownames(dat)
    }})
    if (inherits(x, "try-error"))
    {
       rownames(datamat) <- rownames(dat);
      eset <- new("ExpressionSet", phenoData = as(sampledat, "AnnotatedDataFrame"), 
          annotation = GPL, exprs = as.matrix(datamat))
    } else 
      eset <- new("ExpressionSet", phenoData = as(sampledat, "AnnotatedDataFrame"), 
          annotation = GPL, featureData = fd, exprs = as.matrix(datamat))
    return(list(GPL = as.character(sampledat[1, grep("platform_id", 
        colnames(sampledat), ignore.case = TRUE)]), eset = eset))
}

getGEO = function (GEO = NULL, filename = NULL, destdir = tempdir(), GSElimits = NULL, 
    GSEMatrix = TRUE, AnnotGPL = FALSE) 
{
    con <- NULL
    if (!is.null(GSElimits)) {
        if (length(GSElimits) != 2) {
            stop("GSElimits should be an integer vector of length 2, like (1,10) to include GSMs 1 through 10")
        }
    }
    if (is.null(GEO) & is.null(filename)) {
        stop("You must supply either a filename of a GEO file or a GEO accession")
    }
    if (is.null(filename)) {
        GEO <- toupper(GEO)
        geotype <- toupper(substr(GEO, 1, 3))
        if (GSEMatrix & geotype == "GSE") {
            return(getAndParseGSEMatrices(GEO, destdir, AnnotGPL = AnnotGPL))
        }
        filename <- getGEOfile(GEO, destdir = destdir, AnnotGPL = AnnotGPL)
    }
    ret <- parseGEO(filename, GSElimits)
    return(ret)
}


getAndParseGSEMatrices = function (GEO, destdir, AnnotGPL) 
{
    GEO <- toupper(GEO)
    stub = gsub("\\d{1,3}$", "nnn", GEO, perl = TRUE)
    gdsurl <- "ftp://ftp.ncbi.nlm.nih.gov/geo/series/%s/%s/matrix/"
    url = sprintf(gdsurl, stub, GEO);
    printFlush(spaste("Attempting to open URL: ", url))
    a <- getURL(sprintf(gdsurl, stub, GEO))
    tmpcon <- textConnection(a, "r")
    b <- read.table(tmpcon)
    close(tmpcon)
    b <- as.character(b[, ncol(b)])
    message(sprintf("Found %d file(s)", length(b)))
    ret <- list()
    for (i in 1:length(b)) {
        message(b[i])
        destfile = file.path(destdir, b[i])
        if (file.exists(destfile)) {
            message(sprintf("Using locally cached version: %s", 
                destfile))
        }
        else {
            download.file(sprintf("ftp://ftp.ncbi.nlm.nih.gov/geo/series/%s/%s/matrix/%s", 
                stub, GEO, b[i]), destfile = destfile, mode = "wb", 
                method = getOption("download.file.method.GEOquery"))
        }
        ret[[b[i]]] <- parseGSEMatrix(destfile, destdir = destdir, 
            AnnotGPL = AnnotGPL)$eset
    }
    return(ret)
}

getGEOfile = function (GEO, destdir = tempdir(), AnnotGPL = FALSE, amount = c("full", "brief", "quick", "data")) 
{
    amount <- match.arg(amount)
    geotype <- toupper(substr(GEO, 1, 3))
    mode <- "wb"
    GEO <- toupper(GEO)
    stub = gsub("\\d{1,3}$", "nnn", GEO, perl = TRUE)
    if (geotype == "GDS") {
        gdsurl <- "ftp://ftp.ncbi.nlm.nih.gov/geo/datasets/%s/%s/soft/%s"
        myurl <- sprintf(gdsurl, stub, GEO, paste0(GEO, ".soft.gz"))
        destfile <- file.path(destdir, paste0(GEO, ".soft.gz"))
    }
    if (geotype == "GSE" & amount == "full") {
        gseurl <- "ftp://ftp.ncbi.nlm.nih.gov/geo/series/%s/%s/soft/%s"
        myurl <- sprintf(gseurl, stub, GEO, paste0(GEO, "_family.soft.gz"))
        destfile <- file.path(destdir, paste(GEO, ".soft.gz", 
            sep = ""))
    }
    if (geotype == "GSE" & amount != "full" & amount != "table") {
        gseurl <- "http://www.ncbi.nlm.nih.gov/geo/query/acc.cgi"
        myurl <- paste(gseurl, "?targ=self&acc=", GEO, "&form=text&view=", 
            amount, sep = "")
        destfile <- file.path(destdir, paste(GEO, ".soft", sep = ""))
        mode <- "w"
    }
    if (geotype == "GPL") {
        if (AnnotGPL) {
            gplurl <- "ftp://ftp.ncbi.nlm.nih.gov/geo/platforms/%s/%s/annot/%s"
            myurl <- sprintf(gplurl, stub, GEO, paste0(GEO, ".annot.gz"))
            printFlush(spaste("Attempting to download GPL URL: ", myurl))
            destfile <- file.path(destdir, paste(GEO, ".annot.gz", 
                sep = ""))
            res = try({
                if (!file.exists(destfile)) {
                  downloadRes = download.file(myurl, destfile, mode = mode, 
                    quiet = TRUE, method = getOption("download.file.method.GEOquery"))
                  if (downloadRes > 0) stop();
                  message("File stored at: ")
                  message(destfile)
                }
                else {
                  message(sprintf("Using locally cached version of %s found here:\n%s ", 
                    GEO, destfile))
                }
            }, silent = TRUE)
            if (!inherits(res, "try-error")) {
                return(invisible(destfile))
            }
            else {
                message("Annotation GPL not available, so will use submitter GPL instead")
            }
        }
        gseurl <- "http://www.ncbi.nlm.nih.gov/geo/query/acc.cgi"
        myurl <- paste(gseurl, "?targ=self&acc=", GEO, "&form=text&view=", 
            amount, sep = "")
        destfile <- file.path(destdir, paste(GEO, ".soft", sep = ""))
        mode <- "w"
        if (!file.exists(destfile)) {
            download.file(myurl, destfile, mode = mode, quiet = TRUE, 
                method = getOption("download.file.method.GEOquery"))
            message("File stored at: ")
            message(destfile)
        }
        else {
            message(sprintf("Using locally cached version of %s found here:\n%s ", 
                GEO, destfile))
        }
        return(invisible(destfile))
    }
    if (geotype == "GSM") {
        gseurl <- "http://www.ncbi.nlm.nih.gov/geo/query/acc.cgi"
        myurl <- paste(gseurl, "?targ=self&acc=", GEO, "&form=text&view=", 
            amount, sep = "")
        destfile <- file.path(destdir, paste(GEO, ".soft", sep = ""))
        mode <- "w"
    }
    if (!file.exists(destfile)) {
        download.file(myurl, destfile, mode = mode, quiet = TRUE, 
            method = getOption("download.file.method.GEOquery"))
        message("File stored at: ")
        message(destfile)
    }
    else {
        message(sprintf("Using locally cached version of %s found here:\n%s ", 
            GEO, destfile))
    }
    invisible(destfile)
}


#===============================================================================================
#
# Map SRX accessions to run numbers (SRRxxxxxxx) 
#
#===============================================================================================

getSRRaccessions = function(SRXlinks, sleep = 10)
{
  n = length(SRXlinks)
  links.http = sub("ftp", "http", SRXlinks);
  srr = rep("", n);
  for (l in 1:n)
  {
    printFlush(spaste("Working on link ", links.http[l]))
    conn = url(links.http[l]);
    lines = readLines(conn);
    close(conn);
    line = lines[grep("SRR", lines)];
    if (length(line)==0) stop("Cannot find 'SRR' in lines\n", paste(lines, collapse = "\n"));
    srr[l] = multiSub(c(".*href=\"", "/\\\".*"), c("", ""), line);
    Sys.sleep(sleep);
  }
  srr;
}
    

