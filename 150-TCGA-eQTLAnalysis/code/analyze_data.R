### analyze_data.R:  analysis file for project `eqtl-tcga`.

### INCLUDES

source("include_project.R")

### PREPROCESSING

if(FALSE) {
	
	message(notification("preprocessing eigengenes and clinical data"))
	indir = file.path(RAWDIR, "WGCNA-TCGA-liver")
	invisible(gc())
	
	message(notification("sample annotation", 1))
	raw = fread(file.path(indir, "Data", "TCGA", "SampleAnnotation", 
		"030-OutlierSamplesRemoved", "LIHC-US-sampleAnnotation.OSR.csv"))
	sampid = raw$icgc_sample_id 
	patid = str_extract(raw$submitted_sample_id, "^[^-]+-[^-]+-[^-]+")
		# extracts patient identifiers of the form "TCGA-XX-YYYY"
	spectype = tolower(str_trim(str_extract(raw$specimen_type, "[A-Za-z]+ ")))
		# "normal" = normal tissue adjacent to primary tumor, "primary" = solid tissue primary tumor, "recurrent" = solid tissue recurrent tumor
	isnorm = (spectype == "normal")
		# normal tissue adjacent to tumor if TRUE, tumor tissue if FALSE
	sampannot = data.table(sampid, patid, spectype, isnorm)
	setkey(sampannot, patid)
	invisible(gc())
	
	message(notification("module eigengenes", 1))
	maineig = fread(file.path(indir,  "130-NetworkAnalysis-TCGAPhase",
		"Results", "eigengenes-LIHC-US.main.csv"))
	names(maineig)[1] = "sampid"
	setkey(maineig, sampid)
	subeig = fread(file.path(indir,  "130-NetworkAnalysis-TCGAPhase",
		"Results", "eigengenes-LIHC-US.sub.csv"))
	names(subeig)[1] = "sampid"
	setkey(subeig, sampid)
	invisible(gc())
	
	message(notification("patient survival", 1))
	raw = fread(file.path(RAWDIR, "tcga_survival_files", "clinical_data.csv"))
	tumordict = c("TUMOR FREE"=0, "WITH TUMOR"=1)
	patsurv = raw[, .(patid=bcr_patient_barcode, type=type, time=OS.time/365.24, 
		event=OS, tumor=tumordict[raw$tumor_status])] # could add covariates
	patsurv = patsurv[!is.na(time), ]
	setkey(patsurv, patid)
	#save(patsurv, file=file.path(PROCDIR, "patient_survival.RData"))
	invisible(gc())
	
	message(notification("combining and saving", 1))
	clindat = patsurv[sampannot]
	clindat$spectype = factor(clindat$spectype,
		levels=c("primary", "recurrent", "normal"))
	setorder(clindat, patid, spectype, sampid)
		# we are currently treating multiple samples from single patients as separate samples, but that could change: for example, we eliminate duplicate patient IDs so that (according to the ordering above) primary tumor samples are retained by preference, then "recurrent" (secondary?) tumors, then normal samples (DD 2019-11-16)
	eigens = list(main = as.data.frame(as.matrix(maineig, rownames="sampid")),
	              sub  = as.data.frame(as.matrix(subeig,  rownames="sampid")))
	eigens = mlapply(eigens, function(eig) eig[clindat$sampid, ])
	save(clindat, eigens, file=file.path(PROCDIR, 
		"preprocessed_main.RData"))
	invisible(gc())
}

if(FALSE) {
	
	message(notification("preprocessing SNPs"))
	
	message(notification("includes", 1))
	library(biomaRt)
	invisible(gc())
	
	message(notification("reading significant eQTL SNPs", 1))
	snp_sig = unique(fread(
		select=c("SNP", "major", "minor", "MAF"),
		file=file.path(RAWDIR, "LIHC_eQTL.tsv")))
	names(snp_sig) = c("snp", "major", "minor", "maf")
	setkey(snp_sig, snp)
	invisible(gc())
	
	message(notification("getting new eQTL data for signifcant SNPS", 1))
	mart = useMart("ENSEMBL_MART_SNP", host="grch37.ensembl.org", 
		dataset="hsapiens_snp") # <https://support.bioconductor.org/p/65831/>
	res = getBM(filters="snp_filter", values=snp_sig$snp, mart=mart, 
		attributes=c("refsnp_id", "chr_name", "chrom_start", "chrom_end", "allele", "minor_allele", "minor_allele_freq"))
	res$major_allele = mmapply(res$minor_allele, strsplit(res$allele, "/"), 
		FUN=function(x, y) first(setdiff(y, x))) # It is surprisingly difficult to get major alleles as such from Biomart.  What this does is split all the alleles listed in `res$allele` by slash, remove `res$minor_allele`, and return the first element of what's left.  It may not be perfect, but it is the approach that lines up best with previous results.
	snp_new = with(res, data.table(snp=refsnp_id, chr=chr_name, 
		start=pmin(chrom_start, chrom_end), end=pmax(chrom_start, chrom_end), major=major_allele, minor=minor_allele, maf=minor_allele_freq))
	chroms = c(1:22, "X", "Y", "M")
	snp_new = snp_new[chr %in% chroms]
	snp_new$major = clean_allele_vector(snp_new$major)
	snp_new$minor = clean_allele_vector(snp_new$minor)
	snp_new = snp_new[major != minor]
	setkey(snp_new, snp)
	invisible(gc())
	
	message(notification("saving results", 1))
	csvwrite(snp_new, file=file.path(PROCDIR, "snp_info.csv"))
	save(snp_new, snp_sig, file=file.path(PROCDIR, "preprocessed_snp.RData"))
	invisible(gc())
}

if(0 && FALSE && !TRUE) { # this is an OLD version of the preprocessing that includes the giant all-eQTLs list rather than just the significant ones; it is preserved here _just in case_ we ever need to go back to that, but should probably never again be executed as-is, as indicated by the overkill in the `if` statement ... (DD 2019-12-09)
	
	message(notification("preprocessing SNPs"))
	
	message(notification("includes"))
	library(biomaRt)
	invisible(gc())
	
	message(notification("reading and combining current SNPs", 1))
	
	message(notification("positional data", 2))
	snp_pos = unique(fread(select=1:4,
		file=file.path(RAWDIR, "LIHC_tumor.cis_eQTL.csv")))
	names(snp_pos) = c("snp", "chr", "pos", "alleles")
	setkey(snp_pos, snp)
	majmin = do.call(rbind, strsplit(snp_pos$alleles, "/"))
	colnames(majmin) = c("major", "minor")
	snp_pos = cbind(snp_pos, majmin)
	set(snp_pos, j="alleles", value=NULL)
	invisible(gc())
	
	message(notification("cis-significance", 2))
	snp_sig = unique(fread(
		select=c("SNP", "major", "minor"),
		file=file.path(RAWDIR, "LIHC_eQTL.tsv")))
	names(snp_sig) = c("snp", "major", "minor")
	setkey(snp_sig, snp)
	invisible(gc())
	
	message(notification("combined", 2))
	snp_old = merge(snp_pos, snp_sig, all=T)
	invisible(gc())
	
	message(notification("getting new eQTL data", 1))
	
	mart = useMart("ENSEMBL_MART_SNP", host="grch37.ensembl.org", 
		dataset="hsapiens_snp") # <https://support.bioconductor.org/p/65831/>
	res = getBM(filters="snp_filter", values=snp_old$snp, mart=mart, 
		attributes=c("refsnp_id", "chr_name", "chrom_start", "chrom_end", "allele", "minor_allele", "minor_allele_freq"))
	res$major_allele = mmapply(res$minor_allele, strsplit(res$allele, "/"), 
		FUN=function(x, y) first(setdiff(y, x))) # It is surprisingly difficult to get major alleles as such from Biomart.  What this does is split all the alleles listed in `res$allele` by slash, remove `res$minor_allele`, and return the first element of what's left.  It may not be perfect, but it is the approach that lines up best with previous results.
	snp_new = with(res, data.table(snp=refsnp_id,
		chr=paste("chr", chr_name, sep=""), start=pmin(chrom_start, chrom_end), end=pmax(chrom_start, chrom_end), major=major_allele, minor=minor_allele, maf=minor_allele_freq))
	chroms = paste("chr", c(1:22, "X", "Y", "M"), sep="")
	snp_new = snp_new[chr %in% chroms]
	setkey(snp_new, snp)
	invisible(gc())
	
	message(notification("saving results", 1))
	csvwrite(snp_new, file=file.path(PROCDIR, "snp_info.csv"))
	save(snp_new, snp_old, snp_pos, snp_sig,
		file=file.path(PROCDIR, "preprocessed_snp.RData"))
	invisible(gc())
}

### FUNCTION DEFINITIONS

# Given a numerical vector `x`, return a vector of factors "high" and "low" (and "mid" if `three` is TRUE) corresponding to the elements of `x`.  The effect of `cut` depends on `method`; see below.  Argument `method` describes how the stratification is carried out:
# - "quantile":  by quantiles based on `cut`:  if `three` is TRUE, then "low" for less than or equal to the `cut`th quantile of `x`, "high" for greater than the `1-cut`th quantile, and "mid" otherwise; if `three` is FALSE, then "low" for less than or equal to the `cut`th quantile of `x` and "high" otherwise.  Note that `three=FALSE` and `cut=1/2` means "low" is less than or equal to, and "high" is greater than, the median of `x`.
# - "zscore":  like "quantile" but using the quantiles of the normal distribution based on the mean and standard deviation of `x`.  Here `three=FALSE` and `cut=1/2` means "low" is less than or equal to, and "high" is greater than, the mean of `x`.
# - "mixture":  apply a normal mixture model to `x` (3 components if `three` is TRUE, 2 if FALSE) and stratify based on component assignment.  In this case, `cut` has no effect.
stratify <- function(x, three=FALSE, cut=ifelse(three, 1/3, 1/2),
	method=c("quantile", "zscore", "mixture"))
{
	method = match.arg(method)
	retn = vector(length=length(x), mode="character")
	
	if(method == "quantile") {
		if(three)
			crit = quantile(x, c(cut, 1-cut))
		else
			crit = quantile(x, cut)
	} else if(method == "zscore") {
		mu = mean(x)
		sigma = sd(x)
		if(three)
			crit = c(qnorm(cut, mu, sigma), qnorm(1-cut, mu, sigma))
		else
			crit = qnorm(cut, mu, sigma)
	} else if(method == "mixture") {
		if(three) {
			mod = mixmod(x, 3)
			mindex = which.min(mod$params$observed$mean)
			maxdex = which.max(mod$params$observed$mean)
			crit = c(max(x[mod$assignment == mindex]),
			         min(x[mod$assignment == maxdex]) - LC_EPS)
				# `- LC_EPS` so that the least value of `x` assigned to the highest-mean component will make the `> crit[[2]]` cutoff below
		} else {
			mod = mixmod(x, 2)
			mindex = which.min(mod$params$observed$mean)
			crit = max(x[mod$assignment == mindex])
		}
	} else stop("Inconceivable!") # this should never be able to happen

	if(three) {
		retn[x <= crit[[1]]] = "low"
		retn[x > crit[[2]]] = "high"
		retn[retn == ""] = "mid"
		factor(retn, levels=c("high", "low", "mid"))
	} else {
		retn[x <= crit] = "low"
		retn[retn == ""] = "high"
		factor(retn, levels=c("high", "low"))
	}
}

### MAIN EXECUTION

if(0) {
	
	message(notification("survival analysis by eigengene"))
	
	message(notification("loading data", 1))
	u = load(file.path(PROCDIR, "preprocessed_main.RData"))
	tedat = clindat[, .(time, event)]
	invisible(gc())
	
	message(notification("unstratified tests", 1))
	raw_pvals <- function(x)
	{
		dloc = cbind(tedat, x)
		#!#pmod = survreg(Surv(time+LC_EPS, event) ~ x, data=dloc)
		#!#	# `+LC_EPS` because `survreg` doesn't handle 0 times well
		#!#para = summary(pmod)$table[2, "p"]
		cmod = coxph(Surv(time, event) ~ x, data=dloc)
		Cox = summary(cmod)$sctest[["pvalue"]]
		#!#data.frame(para, Cox)
		return(Cox)
	}
	#!#mainres = do.call(rbind, mlapply(eigens$main, raw_pvals))
	#!#subres = do.call(rbind, mlapply(eigens$sub, raw_pvals))
	#!#names(mainres) = paste("main", names(mainres))
	#!#names(subres) = paste("sub", names(subres))
	#!#mainres = name_by_rownames(round(mainres, 5), "module")
	#!#subres = name_by_rownames(round(subres, 5), "module")
	#!#rept = merge(mainres, subres, by="module", all=TRUE)
	#!#csvwrite(rept, file=file.path(REPTDIR,
	#!#	"unstratified_regression_pvals.csv"))
	mainpvals = msapply(eigens$main, raw_pvals)
	subpvals = msapply(eigens$sub, raw_pvals)
	csvwrite(name_by_rownames(data.frame(pval=mainpvals), "module"),
		file=file.path(REPTDIR, "unstratified_main_module_pvals.csv"))
	csvwrite(name_by_rownames(data.frame(pval=subpvals), "module"),
		file=file.path(REPTDIR, "unstratified_sub_module_pvals.csv"))
	invisible(gc())
	
	message(notification("stratified tests", 1))
	
	message(notification("function definitions", 2))
	.strattest_engine <- function(x, title)
	{
		# prepare data and calculate p-value for log-rank test of differences
		dloc = cbind(tedat, x)
		p = pval(survdiff(Surv(time, event) ~ x, data=dloc))
		
		# fit and plot Kaplan-Meier curves
		fit = survfit(Surv(time, event) ~ x, data=dloc)
		nms = gsub("x=", "", names(fit$strata))
		plt = ggsurvplot(fit, data=dloc, xlab="Time (years)",
			title=title, legend="top", legend.labs=nms, legend.title="",
			pval=sprintf("p = %.5f", p))
		plt$plot = plt$plot + theme( # main (title) text
			plot.title = element_text(hjust=0.5, face="bold"))
		plt$plot = plt$plot + theme( # legend text
			legend.text = element_text(size=14))
		print(plt)
		
		# return p-value for further use
		return(p)
	}
	strattest <- function(x, type=c("main", "sub"))
	{
		# set up base file information and plot
		type = match.arg(type)
		dev.new(height=4, width=6)
		fbase = paste("KM", type, "%s_%s_strata", sep="_")
			# formatting arguments in order are module, then stratification method and number of strata, e.g. "ME1" and "quantile 2"
		fbase = paste(fbase, "eps", sep=".")
		kmdir = file.path(REPTDIR, paste(type, "KM_plots", sep="_"))
		if(!dir.exists(kmdir)) dir.create(kmdir, recursive=TRUE)
		pathbase = file.path(REPTDIR, kmdir, fbase)
		mainbase = paste(type, "%s %s strata")
			# same formatting arguments as `fbase`
		
		# get data names and iterate through
		xnames = names(x)
		retn = vector(mode="list", length=length(xnames))
		names(retn) = xnames
		for(xname in names(x))
		{
			raw = x[[xname]]
			strat = data.table(
				"median 2"   = stratify(raw, FALSE, 1/2, "quantile"),
				"quantile 3" = stratify(raw, TRUE,  1/3, "quantile"),
				"qu.wide 3"  = stratify(raw,  TRUE,  1/4, "quantile"),
				"mean 2"     = stratify(raw, FALSE, 1/2, "zscore"),
				"zscore 3"   = stratify(raw, TRUE,  1/3, "zscore"),
				"zs.wide 3"  = stratify(raw, TRUE,  1/4, "zscore"),
				"mix 2"      = stratify(raw, FALSE, method="mixture"),
				"mix 3"      = stratify(raw, TRUE,  method="mixture"))
			nms = names(strat)
			res = setNames(rep(as.numeric(NA), length(nms)), nm=nms)
			for(nm in nms)
			{
				title = sprintf(mainbase, xname, nm)
				xloc = strat[[nm]]
				if(length(unique(xloc)) > 1) {
					res[[nm]] = .strattest_engine(xloc, title)
					fpath = sprintf(pathbase, xname, nm)
					fpath = gsub(" +", "_", fpath)
					dev.copy2eps(file=fpath)
				} # "else": no good strata, just skip it
			}
			retn[[xname]] = res
		}
		
		# close plot window, package and return
		dev.off()
		name_by_rownames(do.call(rbind, retn), "module")
	}
	invisible(gc()) # end local "function definitions" block
	
	message(notification("main", 2))
	x = eigens$main[names(mainpvals)[mainpvals <= 0.05]]
	mainrept = strattest(x, type="main")
	csvwrite(mainrept, file=file.path(REPTDIR, 
		"stratified_main_module_pvals.csv"))
	invisible(gc())
	
	message(notification("sub", 2))
	x = eigens$sub[names(subpvals)[subpvals <= 0.05]]
	subrept = strattest(x, type="sub")
	csvwrite(subrept, file=file.path(REPTDIR, 
		"stratified_sub_module_pvals.csv"))
	invisible(gc())
	
	message(notification("formatted tables for report", 1))
	ow = options()$width
	options(width=10000, knitr.kable.NA="")
	sink(file.path(REPTDIR, "stratified_main_module_pvals.txt"))
	print(kable(mainrept, digits=5, row.names=FALSE))
	sink()
	sink(file.path(REPTDIR, "stratified_sub_module_pvals.txt"))
	print(kable(subrept, digits=5, row.names=FALSE))
	sink()
	options(width=ow, knitr.kable.NA="NA")
	invisible(gc())
	
	message(notification("other report stuff", 1))
	catprint(colSums(numericols(mainrept) <= 0.05, na.rm=T))
	catprint(colSums(numericols(subrept) <= 0.05, na.rm=T))
	cat("\n\n")
	invisible(gc())
}

if(0) {
	
	message(notification("extracting eQTLs for modules"))
	
	message(notification("includes", 1))
	library(seqminer)
	invisible(gc())
	
	message(notification("loading data", 1))
	u1 = load(file.path(PROCDIR, "preprocessed_main.RData"))
	u2 = load(file.path(PROCDIR, "preprocessed_snp.RData"))
	invisible(gc())
	
	message(notification("genotypes", 1))
		# for the sake of speed, and keeping SNP identifiers, here we read the file in chunks and then put it together rather than relying on `tabix.read.table`, even though that would be more convenient
	
	message(notification("reading column names", 2))
	# code copied _almost_ verbatim from `tabix.read.table`
	inpath = file.path(RAWDIR, 
		"PCA.r1.TCGAbarcode.merge.tnSwapCorrected.10389.vcf.gz")
	header = .Call("readTabixHeader", inpath, PACKAGE="seqminer")
	header = header[nchar(header) > 0]
	header = header[length(header)]
	header = sub("^#", "", header)
	header = strsplit(header, "\t")[[1]]
	invisible(gc())
	
	message(notification("cleaning/subsetting column names", 2))
	tcga_idx = grep("^TCGA", header)
	header[tcga_idx] = str_extract(
		header[tcga_idx], "^TCGA-[A-Za-z0-9]+-[A-Za-z0-9]+")
	col_idx = which(header %in% clindat$patid)
	invisible(gc())
	
	message(notification("reading data body", 2))
	readdat <- function(i)
	{
		# get SNP information and format range
		snp = snp_new[i,]
		range = with(snp, sprintf("%s:%d-%d", chr, start, end))
		
		# read data and put into matrix
		res = do.call(rbind, strsplit(tabix.read(inpath, range), "\t"))
		if(is.null(res)) return(NULL) # nothing we can do here
		colnames(res) = header # assume we got through the above
		idx = which.min(abs(as.numeric(res[,"POS"]) - snp$start))
		res = res[idx, col_idx] # row closest to SNP, dropping extra columns
		
		# convert numerical allele codes ("./.", "0/1", etc.) to zygosity codes; this works because the lower number is always first in "allele1/allele2" format, i.e. "." comes before "0", "0" comes before "1", etc.
		retn = vector(mode="character", length=length(res))
		names(retn) = names(res)
		retn[grep("^\\.", res)]   = "0" # major allele homozygous
		retn[grep("^0", res)]     = "1" # minor allele heterozygous
		retn[grep("^[1-9]", res)] = "2" # minor allele homozygous
		
		# make sure "retn" has at least 10 non-major allele hits, ; if so, return with identifying information, otherwise return NULL
		if(sum(retn != "0") >= 10) c(snp=snp$snp, retn) else NULL
	}
	genotype = do.call(rbind, mlapply(seq.int(nrow(snp_new)), readdat))
	snps = genotype[, 1] # preserve SNP identifiers
	genotype = genotype[, -1] # now there is no SNP, only genotype
	genotype = apply(genotype, 2, as.integer) # "0" -> 0, etc.
	rownames(genotype) = snps # put SNP identifiers back
	invisible(gc())
	
	message(notification("saving", 1))
	save(genotype, file=file.path(PROCDIR, "module_eQTL_genotype.RData"))
	invisible(gc())
}

if(0) {
	
	message(notification("module eQTLs"))
	
	message(notification("includes", 1))
	library(MatrixEQTL)
	invisible(gc())
	
	message(notification("loading data", 1))
	u1 = load(file.path(PROCDIR, "preprocessed_main.RData"))
	u2 = load(file.path(PROCDIR, "preprocessed_snp.RData"))
	u3 = load(file.path(PROCDIR, "module_eQTL_genotype.RData"))
	cdat = copy(clindat) # local version we can manipulate
	invisible(gc())
	
	message(notification("organizing data", 1))
	eigmat = t(eigens$sub[, c("ME5.1", "ME17", "ME38")])
	common_samples = intersect(cdat$sampid, colnames(eigmat))
	common_patients = intersect(cdat$patid, colnames(genotype))
	cdat = cdat[  (sampid %in% colnames(eigmat))
	            & (patid %in% colnames(genotype))  ]
	eigmat = eigmat[, cdat$sampid]
	genmat = genotype[, cdat$patid]
	colnames(genmat) = cdat$sampid
	invisible(gc())
	
	message(notification("prepping for Matrix eQTL", 1))
	eigdat = SlicedData$new()$CreateFromMatrix(eigmat)
	gendat = SlicedData$new()$CreateFromMatrix(genmat)
	invisible(gc())
	
	message(notification("running linear model", 1))
	res = suppressMessages(Matrix_eQTL_engine(gendat, eigdat,
		pvOutputThreshold=1, verbose=FALSE, useModel=modelLINEAR,
		output_file_name=NULL)$all$eqtls)
	res = with(res, data.frame(SNP=snps, module=gene, pval=pvalue, qval=FDR,
		coef=beta, sdcoef=beta/statistic, tstat=statistic))
	csvwrite(res, file=file.path(REPTDIR, "module_eQTL_linear.csv"))
	invisible(gc())
	
	message(notification("running three-way ANOVA", 1))
	res = suppressMessages(Matrix_eQTL_engine(gendat, eigdat,
		pvOutputThreshold=1, verbose=FALSE, useModel=modelANOVA,
		output_file_name=NULL)$all$eqtls)
	res = with(res, data.frame(SNP=snps, module=gene, pval=pvalue, qval=FDR))
	csvwrite(res, file=file.path(REPTDIR, "module_eQTL_3wayANOVA.csv"))
	invisible(gc())
	
	message(notification("plotting", 1))
	
	message(notification("function definitions", 2))
	alleledict = c("0"="ref hom", "1"="alt/ref", "2"="alt hom")
	module_snp_plot <- function(module, SNP)
	{
		y = eigmat[module, ]
		x = alleledict[as.character(genmat[SNP, ])]
		x = factor(x, levels=alleledict)
		boxplot(y ~ x, xlab="SNP alleles", ylab="eigengene expression",
			main=sprintf("%s eigengene expression by SNP %s", module, SNP))
		fname = sprintf("module_%s_SNP_%s_expr.eps", module, SNP)
		dev.copy2eps(file=file.path(REPTDIR, fname))
	}
	invisible(gc())
	
	message(notification("plotting", 2))
	dev.new(width=10)
	isgood = ((res$qval <= 0.05) | (res$SNP %in% c("rs103197", "rs8076470")))
	rr = res[isgood, c("module", "SNP")]
	for(i in 1:nrow(rr)) do.call(module_snp_plot, rr[i,])
	dev.off()
	invisible(gc())
}

if(1) {
	
	message(notification("sample size report"))
	# see Grant's 2021-02-06 @ 1637 e-mail "Fwd: Editorial Decision on manuscript CELL-SYSTEMS-D-20-00235R1" et seq.
	
	message(notification("loading data", 1))
	u1 = load(file.path(PROCDIR, "preprocessed_main.RData"))
	u2 = load(file.path(PROCDIR, "preprocessed_snp.RData"))
	u3 = load(file.path(PROCDIR, "module_eQTL_genotype.RData"))
	invisible(gc())
	
	message(notification("sample and subject counts", 1))
	sample_count_report = data.table(
		"clinical samples" = nrow(clindat),
		"clinical subjects" = length(unique(clindat$patid)),
		"genotype subjects" = ncol(genotype))
	csvwrite(sample_count_report,
		file=file.path(REPTDIR, "sample_count_report.csv"))
	invisible(gc())
	
	message(notification("genotype counts", 1))
	alleledict = c("0"="ref hom", "1"="alt/ref", "2"="alt hom")
	typestab = snp_sig[, .(snp,
		"0" = paste(major, major, sep=""),
		"1" = paste(major, minor, sep=""),
		"2" = paste(minor, minor, sep=""))]
	ctfn <- function(i)
	{
		snpname = rownames(genotype)[i]
		counts = table(factor(genotype[i, ], levels=names(alleledict)))
		names(counts) = paste(alleledict, "count")
		alleles = typestab[snpname]
		names(alleles)[2:4] = paste(alleledict, "allele")
		cbind(alleles, as.data.table(t(as.matrix(counts))))
	}
	genotype_count_report = rbindlist(mlapply(seq.int(nrow(genotype)), ctfn))
	ord = order(as.numeric(gsub("rs", "", genotype_count_report$snp)))
	genotype_count_report = genotype_count_report[ord, ]
	csvwrite(genotype_count_report,
		file=file.path(REPTDIR, "genotype_count_report.csv"))
	invisible(gc())
}

### CLEANUP

message(notification("cleaning up"))
closeAllConnections()
invisible(gc())
message(notification("ALL DONE!"))
cat("\n")