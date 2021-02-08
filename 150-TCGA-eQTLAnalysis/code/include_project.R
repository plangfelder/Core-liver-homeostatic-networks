### include_project.R:  libraries, "CONSTANTS", etc. for project `eqtl-tcga`.

## PROJECT-SPECIFIC LIBRARIES

library(lcmix) # for mixture model stratification
library(survival) # for basic survical fitting
library(survminer) # for nice survival plots

## DANIEL'S STANDARD LIBRARIES

library(ARC.utils) # also imports dplyr, parallel, R.methodsS3
library(matrixStats) # we always end up needing `rowSds`, etc.
library(stringr) # `str_*` functions are often better than base R
library(readxl) # since everything is in Excel, apparently ...
library(reshape2) # `melt` and `(a)cast` should never be deprecated!
library(data.table) # `fread` is great, also own `melt` and `dcast`
library(knitr) # for `kable` to print tables usable in Markdown etc.

### "CONSTANTS"

## directories

RAWDIR <- "../raw_data"
PROCDIR <- "../processed_data"
REPTDIR <- "../reports"

## standard settings

options(stringsAsFactors       = FALSE,
        datatable.showProgress = FALSE)

### FUNCTION DEFINITIONS

## dplyr functions that can get overwritten by other methods

select <- dplyr::select
filter <- dplyr::filter
mutate <- dplyr::mutate
transmute <- dplyr::transmute
arrange <- dplyr::arrange
first <- dplyr::first
last <- dplyr::last
nth <- dplyr::nth
sym <- dplyr::sym # useful in using column names in string form for `arrange`; <https://stackoverflow.com/questions/27034655/how-to-use-dplyrarrangedesc-when-using-a-string-as-column-name?lq=1>

## project-specific settings

## project-specific functions

# Given a vector of strings `x`, convert any element of `x` but gene strings, i.e. strings containing only the letters 'A', 'G', 'C', and 'T', to blanks, and return the cleaned vector.
clean_allele_vector <- function(x)
{ 
	badidx = grep("[^ATGC]", x)
	x[badidx] = ""
	return(x)
}

# Given a vector of factors or factor-like objects `x`, calculate the entropy of their distribution.
factent <- function(x)
{
	cts = table(x, useNA="ifany")
	props = cts / sum(cts)
	-2 * sum(na.sub(props * log(props)))
}


# Get the chi-squared p-value from an `object` of class "survdiff"; see <https://stackoverflow.com/questions/36368652/how-to-get-survdiff-returned-p-value>
setMethodS3("pval", "survdiff", function(object, ...)
	pchisq(object$chisq, length(object$n)-1, lower.tail=FALSE), 
	conflict="quiet")

### CLEANUP

closeAllConnections()
invisible(gc())