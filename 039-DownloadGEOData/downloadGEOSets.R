# Download a few data sets (technically Series) from GEO automatically

source("../Functions/networkFunctions-extras-18.R");
library(GEOquery);
source("../Functions/GEOfunctions.R");

identifiers = c("GSE126848")
citations = c("Suppli2019-HumanNAFLD");

nSets = length(identifiers);

geoData = list();
for (set in 1:nSets)
  geoData[[set]] = processGEOExpressionSeries(identifier = identifiers[set], 
                              spaste(dirBase = "../Data/", citations[set]), verbose = 5);

sampleAnnot0 = read.table(
  "../Data/Suppli2019-HumanNAFLD-GSE126848/SampleAnnotation/010-AsSupplied/GSE126848-sampleAnnotation.txt", sep = "\t",
  quote = "", comment.char = "", header = TRUE)

sampleAnnot = processGEOSampleAnnotation(sampleAnnot0, usePart = "2nd", checkNames = FALSE)
names(sampleAnnot) = sub("ProbeID", "SampleID", names(sampleAnnot));

write.csv(sampleAnnot, 
  file = "../Data/Suppli2019-HumanNAFLD-GSE126848/SampleAnnotation/010-AsSupplied/GSE126848-sampleAnnotation-processed.csv",
  quote = FALSE, row.names = FALSE);
