# -----------------------------------------------------------------------------
# This script demonstrates the framework for removing non-malignant cell-type 
# signatures from tumor methylation profiles, followed by classification of 
# the adjusted profiles.
#
# The approach is model-agnostic, but this script uses the 
# Heidelberg Classifier v11b4 (Capper et al.). You may need to adapt the 
# classification steps to fit your own classifier.
# -----------------------------------------------------------------------------

library(randomForest)
library(glmnet)

# Load classifier 
# When using the Heidelberg Classifier, you'll need: 
# rf.pred, calfit, cgroups, reflist
load("/path/to/your/classifier/rfpred.v11b4.RData")
source("/path/to/your/classifier/R/MNPpredict_betas.R")

# Load beta values (from 01_preprocessing or your own data)
beta <- readRDS(file.path("results", "beta_values.rds"))
beta_proc <- as.data.frame(beta[rownames(rf.pred$importance), ])
colnames(beta_proc) <- colnames(beta)

# Load reference signatures (use provided or replace with your own)
ref <- read.csv(file.path("data", "ref.csv"), check.names = FALSE)
rownames(ref) <- ref$CpG
ref_al <- ref[rownames(rf.pred$importance), ]
ref_list <- setNames(
  lapply(ref_al[-1], function(x) setNames(as.numeric(x), ref_al$CpG)),
  names(ref_al)[-1]
)

# Define fractions of non-malignant cell signatures to remove (adjust if needed)
fractions <- seq(0, 0.98, by = 0.02)

# Create purified profiles by removing reference signatures
beta_adj_all <- do.call(
  cbind,
  unlist(
    lapply(names(ref_list), function(ct) {
      ct_median <- ref_list[[ct]]
      lapply(fractions, function(x) {
        mat <- (beta_proc - x * ct_median) / (1 - x)
        colnames(mat) <- paste(ct, x, colnames(mat), sep = "_")
        mat
      })
    }),
    recursive = FALSE
  )
)

# Classify adjusted profiles
prob_matrix <- MNPpredict_betas(
  beta_adj_all,
  calibrate = TRUE,
  type      = "prob",
  MCF       = FALSE
)

# Group into Methylation Class Families (MCFs), if applicable to your classifier
mcf_probs <- sapply(cgroups, function(grp) rowSums(prob_matrix[, grp, drop = FALSE]))
remaining_classes <- setdiff(colnames(prob_matrix), unlist(cgroups))
prob_matrix_mcf   <- cbind(mcf_probs, prob_matrix[, remaining_classes, drop = FALSE])

# Save results
saveRDS(prob_matrix_mcf, file.path("results", "prob_matrix_mcf.RDS"))
saveRDS(prob_matrix, file.path("results", "prob_matrix.RDS"))
