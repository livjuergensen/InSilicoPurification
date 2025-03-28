# -----------------------------------------------------------------------------
# Download of example data and preprocessing workflow as performed in the study.
#
# This script uses MNPpreprocessIllumina (Capper et al., 2018) for preprocessing.
# If this function is not available, similar methods (e.g., preprocessIllumina) may be used.
# Reference profiles should be processed the same way as tumor samples.
# -----------------------------------------------------------------------------

library(GEOquery)
library(minfi)
library(R.utils)
#library(IlluminaHumanMethylationEPICmanifest) -> not sure if necessary

dir.create("example_data")
dir.create("results")

# Download example data from GEO (Sturm et al., 2023)
geo_ids <- c( "GSM6629544", "GSM6629974", "GSM6630028")

for (gsm in geo_ids) {
  getGEOSuppFiles(GEO = gsm, makeDirectory =  FALSE, baseDir = "example_data/")
}

gz_files <- list.files("example_data", pattern = "\\.gz$", full.names = TRUE)
sapply(gz_files, gunzip, remove = TRUE)

# Preprocessing
source("/path/to/function/MNPprocessIDAT_functions.R")
rgset <- read.metharray.exp(base = "example_data", force = TRUE, extended = TRUE)
sampleNames(rgset) <- sub(".*_(\\d+_R\\d+C\\d+)$", "\\1", sampleNames(rgset))
mset  <- MNPpreprocessIllumina(rgset, bg.correct = TRUE, normalize = "controls")
beta  <- getBeta(mset, type = "Illumina")

# Save beta values
saveRDS(beta, file.path("results", "beta_values.rds"))
