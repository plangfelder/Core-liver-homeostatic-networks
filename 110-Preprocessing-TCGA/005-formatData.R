# Standard preprocessing

source("../Functions/networkFunctions-extras-19.R")
source("../Functions/labelPoints2-01.R");
source("../Functions/outlierRemovalFunctions.R")
source("../Functions/preprocessing-General-013.R")
source("../Functions/GNVFunctions-016.R")

library(anRichment)

source("../../../../RLibs/inHouseGeneAnnotation/anRichmentMethods/R/enrichmentAnalysis.R");

# Load data

baseDir = "../Data/TCGA/";
exprFiles = file.path(baseDir, "Expression/010-RawData",
    c("exp_seq.LICA-FR-V2.tsv.gz", "exp_seq.LIHC-US.tsv.gz", "exp_seq.LIRI-JP.tsv.gz"))

setNames = c("LICA-FR", "LIHC-US", "LIRI-JP");

data0 = list();
for (fi in 1:length(exprFiles))
{
  f = exprFiles[fi];
  printFlush(f);
  data0[[fi]] = read.tsv(gzfile(f));
}

# Read the old data and transfer annotation to new data.

data0.x = read.tsv(gzfile(file.path(baseDir, "Expression/010-RawData", "exp_seq.LICA-FR.tsv.gz")));

samples.new = unique(data0[[1]]$analyzed_sample_id);
samples.old = unique(data0.x$submitted_sample_id);

table(samples.new %in% samples.old)

data0[[1]] = cbind(data0[[1]], data0.x[ match(data0[[1]]$analyzed_sample_id, data0.x$submitted_sample_id),
                    c("icgc_donor_id", "icgc_specimen_id", "icgc_sample_id", "submitted_sample_id", "analysis_id")]);

data1 = lapply(data0, dropConstantColumns, checkDuplicates = FALSE);

rm(data0);

gc();


ensemblAnnot = read.csv(gzfile("../Data/Annotation/Ensembl/Homo_sapiens_GRCh37p7.txt.gz"));
NCBIAnnot = read.csv(gzfile("../Data/Annotation/NCBI/geneAnnotation-NCBIHuman37.3.csv.gz"));

if (FALSE)
{
refSeqAnnot = read.table(gzfile("../Data/Annotation/NCBI/refSeq-HomoSapiens-geneAnnotation-2019-07-22.tsv.gz"),
                 header = FALSE, sep = "\t");
names(refSeqAnnot) = c("TaxonID", "Entrez", "Symbol", "RefSeqID");
dim(refSeqAnnot)
}

NCBIAnnot2 = read.tsv(gzfile("../Data/Annotation/NCBI/Homo_sapiens.gene_info.gz"), quote = "", comment.char = "")


# Create an alias to official symbol mapping from NCBI annotation.
split = strsplit(NCBIAnnot2$Synonyms, split = "|", fixed = TRUE);
synonymIndex = anRichmentMethods:::.indexedFlattenedList(split);

synonym2symbol = data.frame(synonym = replaceMissing(synonymIndex$data), symbol = NCBIAnnot2$Symbol[synonymIndex$index]);

# Remove those synonyms that map to more than 1 symbol. This also gets rid of all missing data since those are
# duplicated...
nonuniqueSynonyms = with(synonym2symbol,
    unique(synonym[ which( duplicated(synonym) & !duplicated(spaste(synonym, symbol)))]))

length(nonuniqueSynonyms)
#[1] 3246


synonym2symbol = synonym2symbol[ !synonym2symbol$synonym %in% nonuniqueSynonyms, ]

# sanity check
with(synonym2symbol, any(duplicated(synonym) & !duplicated(spaste(synonym, symbol))))
with(synonym2symbol, any(duplicated(synonym) ))

synonym2symbol = synonym2symbol[!duplicated(synonym2symbol$synonym), ];

dim(synonym2symbol)

# We are left with 61192 synonyms mapping to unique gene symbols.

sets= c(2:3)
data2 = data1;
data2[sets] = lapply(data2[sets], function(.data)
{
  # In set 2, some gene IDs are ?, exclude those
  .data = .data[ .data$gene_id!="?", ]
  # any(duplicated(spaste(.data$gene_id, .data$icgc_sample_id)))
  # Convert synonyms to unique symbols where possible.
  .data$gene_symbol = translateUsingTable(.data$gene_id, synonym2symbol, keepUntranslated = TRUE)
  inNCBIAnnot = unique(.data$gene_symbol) %in% NCBIAnnot2$Symbol;
  inBioCAnnot = unique(.data$gene_symbol) %in% convert2symbol(entrez = allOrgGenes("human"), organism = "human");
  print(table(inNCBIAnnot, inBioCAnnot));

  table(inNCBIAnnot);
  table(inBioCAnnot);
  table(table(.data$gene_symbol))
  tt = table(.data$gene_symbol)
  if (FALSE)
  {
    print(tt[tt>length(unique(.data$icgc_sample_id))]);
    # Look at the counts of the SLC35E2
    rows = .data$gene_id=="SLC35E2";
    table(is.na(.data$raw_read_count[rows]))
    lst = with(.data[rows, ],  tapply(raw_read_count, icgc_sample_id, identity));
    lst;
    # It looks like those are different counts; will write the code to sum them.
  }
})
#### Unfortunately, there are hundreds of genes that have multiple synonyms mapping to them that were used as gene ids in
### the table, no matter which annotation I use. 
### Go with the usual way of using Bioconductor annotation, which seems to be very similar to the newest NCBI annotation
### file; for now keep the data indexed to the IDs provided with the data. Will deal with the non-uniqueness later when
### combining the results and mapping them to mouse; that way, I don't start losing genes earlier than necessary.
  

# Turn the linear format into appropriate matrices.

data1.ext = lapply(data1, function(.data)
{
  geneLevels = unique(.data$gene_id);
  geneLevels = setdiff(geneLevels, "?")  ## get rid of the missing gene IDs
  if (length(grep("^ENSG00000", geneLevels))>10000)
  {
    geneLevels = sub("\\..*", "", geneLevels);
    .data$gene_id = sub("\\..*", "", .data$gene_id);
    print(table(geneLevels%in% ensemblAnnot$Ensembl.Gene.ID));  
    ## Not too bad: 
    #FALSE  TRUE 
    # 6780 51040 
    # Still quite disappointing given the amount of time I spent chasing down a "correct" version of ensembl annotation
    geneEntrez = translateUsingTable(geneLevels, ensemblAnnot[, c("Ensembl.Gene.ID", "EntrezGene.ID")])
  } else
    geneEntrez = convert2entrez(geneLevels, organism = "human");
  .data = .data[.data$gene_id %in% geneLevels, ]
  print(table(is.na(geneEntrez)))
  sampleLevels = unique(.data$icgc_sample_id);
  rawCounts = matrix(0, length(sampleLevels), length(geneLevels));
  dimnames(rawCounts) = list(sampleLevels, geneLevels);
  sampleCol = "icgc_sample_id";
  if (!sampleCol %in% names(.data)) sampleCol = "analyzed_sample_id";
  pind = initProgInd();
  for (s in 1:length(sampleLevels))
  {
    rows = which(.data[[sampleCol]]==sampleLevels[s]);
    index.fwd = match(.data$gene_id[rows], geneLevels);
    index.fwd = index.fwd;
    index.b = tapply(1:length(index.fwd), index.fwd, identity);
    rcIndex = as.numeric(names(index.b));
    maxL = max(sapply(index.b, length));
    for (step in 1:maxL)
    {
      source = sapply(index.b, `[`, step);
      rawCounts[s, rcIndex] = rawCounts[s, rcIndex] + replaceMissing(.data$raw_read_count[rows][source]);
    }
    pind = updateProgInd(s/length(sampleLevels), pind);
  }
  printFlush("");
  sampleData = .data[match(sampleLevels, .data[[sampleCol]]), 
                       c(sampleCol, "icgc_donor_id", "icgc_specimen_id", "submitted_sample_id")];
  gc()
  list(counts = rawCounts, sampleAnnot = sampleData);
});

counts0 = lapply(data1.ext, getElement, "counts");
sampleAnnot0 = lapply(data1.ext, getElement, "sampleAnnot");

if (FALSE)
{
  # Sanity check
  trees = lapply(counts0, function(x) hclust(dist(x), method = "a"));
  gc();
  #sizeGrWindow(14, 10);
  par(mfrow = c(1,3));
  lapply(trees, function(x) {plot(x, labels = FALSE, xlab = "", sub = "", hang = 0); abline(h = 0)});
  # Looks ok except 5 samples in the third data set seem to be duplicated.
}

clinicalData1 = read.csv(file.path("../Data/TCGA/SampleAnnotation/010-AsSupplied",
                   "List_of_HCC_samples_in_TCGA_and_ICGC-ICGCLiverAndNormalTissue-Bioportal.csv"))

sampleData2 = read.csv(file.path("../Data/TCGA/SampleAnnotation/010-AsSupplied",
                   "List_of_HCC_samples_in_TCGA_and_ICGC-ICGCLiverAndNormalTissue.csv"))

sampleDataBySet.sample = lapply(setNames, function(name)
  read.tsv(gzfile(spaste("../Data/TCGA/SampleAnnotation/010-AsSupplied/sampleData.", name, ".tsv.gz"))));

sampleDataBySet.donor = lapply(setNames, function(name)
  read.tsv(gzfile(spaste("../Data/TCGA/SampleAnnotation/010-AsSupplied/donor.", name, ".tsv.gz"))));

sampleDataBySet.specimen = lapply(setNames, function(name)
  read.tsv(gzfile(spaste("../Data/TCGA/SampleAnnotation/010-AsSupplied/specimen.", name, ".tsv.gz"))));

sampleAnnot1 = mymapply(function(ref, sa.s, sa.d, sa.spec)
{
  printFlush("==============================================================");
  printFlush("to sample:")
  print(table(!is.na(match(ref$icgc_specimen_id, sa.s$icgc_specimen_id))))
  out = cbind(ref, sa.s[ match(ref$icgc_specimen_id, sa.s$icgc_specimen_id), !colnames(sa.s) %in% colnames(ref)])
  printFlush("to donor:")
  print(table(!is.na(match(ref$icgc_donor_id, sa.d$icgc_donor_id))))
  out = cbind(out, sa.d[ match(ref$icgc_donor_id, sa.d$icgc_donor_id), !colnames(sa.d) %in% colnames(out)])
  printFlush("to specimen:")
  print(table(!is.na(match(ref$icgc_specimen_id, sa.spec$icgc_specimen_id))))
  out = cbind(out, sa.spec[ match(ref$icgc_specimen_id, sa.spec$icgc_specimen_id), !colnames(sa.spec) %in% colnames(out)])
  out;
}, ref = sampleAnnot0, sa.s = sampleDataBySet.sample, sa.d = sampleDataBySet.donor, sa.spec = sampleDataBySet.specimen)


lapply(sampleAnnot1, function(sa) table(sa$icgc_specimen_id %in% sampleData2$icgc_specimen_id))

# Save the data at this point

exprDir.11 = file.path(baseDir, "Expression/011-RawData-Reformatted");
dir.create(exprDir.11, recursive = TRUE);

nSets = length(setNames)
for (set in 1:nSets)
{
  write.csv.nr(dataForTable(counts0[[set]], transpose = TRUE, IDcolName = "GeneID"),
    file = gzfile(spaste(exprDir.11, "/countMatrix-", setNames[set], ".csv.gz")),
    quote = FALSE);
}

sampleDir.11 = file.path(baseDir, "SampleAnnotation/011-Combined");
dir.create(sampleDir.11, recursive = TRUE);

for (set in 1:nSets)
{
  write.csv.nr(sampleAnnot1[[set]],
    file = gzfile(spaste(exprDir.11, "/sampleData-ICGCcombinedAnnotation-", setNames[set], ".csv.gz")),
    quote = TRUE);
}


