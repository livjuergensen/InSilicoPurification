# ------------------------------------------------------------------------------
# Preprocessing Workflow using example data from Sturm et al. (2023)
#
# This script uses MNPpreprocessIllumina (Capper et al., 2018) for preprocessing.
# If this function is not available, similar methods (e.g., preprocessIllumina) 
# may be used.
# Reference profiles should be processed the same way as tumor samples.
# ------------------------------------------------------------------------------
library(minfi)
library(IlluminaHumanMethylationEPICmanifest)

dir.create("results")

untar("data/idats_example.tar.gz", exdir = "data/")
rgset <- read.metharray.exp(base = "data/", force = TRUE, extended = TRUE)
sampleNames(rgset) <- sub(".*_(\\d+_R\\d+C\\d+)$", "\\1", sampleNames(rgset))

# If MNPpreprocessIllumina is not available to you, you can use:
# mset <- preprocessIllumina(rgset, bg.correct = TRUE, normalize = "controls")
source("/path/to/function/MNPprocessIDAT_functions.R")
mset <- MNPpreprocessIllumina(rgset, bg.correct = TRUE, normalize = "controls")

beta <- getBeta(mset, type = "Illumina")

saveRDS(beta, file.path("results", "beta_values.rds"))
