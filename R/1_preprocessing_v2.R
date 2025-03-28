# -----------------------------------------------------------------------------
# Download of example data from Sturm et al. (2023)
# Preprocessing workflow as performed in the study.
#
# This script uses MNPpreprocessIllumina (Capper et al., 2018) for preprocessing.
# If this function is not available, similar methods (e.g., preprocessIllumina) may be used.
# Reference profiles should be processed the same way as tumor samples.
# -----------------------------------------------------------------------------
library(GEOquery)
library(minfi)
library(R.utils)
library(IlluminaHumanMethylationEPICmanifest)

dir.create("results")

# Download example data from GEO (Sturm et al., 2023)
# This attempts to download 6 Illumina .idat files (Green and Red channel for 3 samples) 
# If the automatic download fails you may need to download the files manually from GEO
geo_ids <- c( "GSM6629544", "GSM6629974", "GSM6630028")

for (gsm in geo_ids) {
  getGEOSuppFiles(GEO = gsm, makeDirectory =  FALSE, baseDir = "data/")
}

gz_files <- list.files("data", pattern = "\\.gz$", full.names = TRUE)
sapply(gz_files, gunzip, remove = TRUE)


# Preprocessing
rgset <- read.metharray.exp(base = "example_data", force = TRUE, extended = TRUE)
sampleNames(rgset) <- sub(".*_(\\d+_R\\d+C\\d+)$", "\\1", sampleNames(rgset))

# If MNPpreprocessIllumina is not available to you, you can use:
# mset <- preprocessIllumina(rgset, bg.correct = TRUE, normalize = "controls")
source("/path/to/function/MNPprocessIDAT_functions.R")
mset <- preprocessIllumina(rgset, bg.correct = TRUE, normalize = "controls")
beta  <- getBeta(mset, type = "Illumina")

saveRDS(beta, file.path("results", "beta_values.rds"))
