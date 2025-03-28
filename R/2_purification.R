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
ref <- ref[rownames(rf.pred$importance), ]
ref_cts <- setNames(
  lapply(ref[-1], function(x) setNames(as.numeric(x), ref$CpG)),
  names(ref)[-1]
)

# Define fractions of non-malignant cell signatures to remove (adjust if needed)
frac <- seq(0, 0.98, by = 0.02)

# Create purified profiles by removing reference signatures
beta_adj <- do.call(cbind, unlist(lapply(names(ref_cts), function(ct) {
  ref_ct <- ref_cts[[ct]]
  lapply(frac, function(f) {
    mat <- (beta_proc - f * ref_ct) / (1 - f)
    colnames(mat) <- paste(ct, f, colnames(mat), sep = "_")
    mat
  })
}), recursive = FALSE))

# Classify adjusted profiles
prob <- MNPpredict_betas(
  beta_adj,
  calibrate = TRUE,
  type      = "prob",
  MCF       = FALSE
)

# Group into Methylation Class Families (MCFs), if applicable to your classifier
mcf_mat <- sapply(cgroups, function(grp) rowSums(prob[, grp, drop = FALSE]))
mcf_rest <- setdiff(colnames(prob), unlist(cgroups))
prob_mcf <- cbind(mcf_mat, prob[, mcf_rest, drop = FALSE])

# Save results
saveRDS(prob,     file.path("results", "prob.RDS"))
saveRDS(prob_mcf, file.path("results", "prob_mcf.RDS"))

# Exclude control classes
ctrl <- c("CONTR, REACT", "CONTR, CEBM", "CONTR, INFLAM",
          "CONTR, HEMI", "CONTR, HYPTHAL", "CONTR, WM",
          "CONTR, ADENOPIT", "CONTR, PINEAL", "CONTR, PONS")

prob_filt     <- prob[, !colnames(prob) %in% ctrl]
prob_mcf_filt <- prob_mcf[, !colnames(prob_mcf) %in% ctrl]

# Identify highest scoring tumor per profile, excluding controls
# Highest scoring methylation class family (MCF)
mcf_idx   <- apply(prob_mcf_filt, 1, which.max)
mcf_class <- colnames(prob_mcf_filt)[mcf_idx]
mcf_score <- prob_mcf_filt[cbind(seq_len(nrow(prob_mcf_filt)), mcf_idx)]

# Highest scoring methylation class (MC)
mc_idx   <- apply(prob_filt, 1, which.max)
mc_class <- colnames(prob_filt)[mc_idx]
mc_score <- prob_filt[cbind(seq_len(nrow(prob_filt)), mc_idx)]

# Summarize classification results per adjusted profile
results <- data.frame(
  idat   = sub("^[^_]*_[^_]*_(.*)$", "\\1", rownames(prob)),
  cell_type     = sub("^(.*?)_.*$", "\\1", rownames(prob)),
  fraction_removed   = sub("^[^_]*_(.*?)_.*$", "\\1", rownames(prob)),
  mcf    = mcf_class,
  mcf_score = mcf_score,
  mc    = mc_class,
  mc_score = mc_score,
  stringsAsFactors = FALSE
)

results$classifiable <- with(results, mcf_score >= 0.9 & mc_score >= 0.5)

# Identify the most confident classification per sample across all adjusted profiles
res_max <- do.call(rbind, lapply(unique(results$idat), function(id) {
  r <- results[results$idat == id, ]
  best <- r[which.max(r$mc_score), ]
  best$cell_type <- ifelse(best$fraction_removed == 0, "n.a.", best$cell_type)
  best
}))

# Identify the most confident classification per sample–cell type combination
res_max_ct <- do.call(rbind, lapply(unique(results$idat), function(id) {
  do.call(rbind, lapply(unique(results$cell_type), function(cell_type) {
    r <- results[results$idat == id & results$ct == cell_type, ]
    if (nrow(r) == 0) return(NULL)
    r[which.max(r$mcf_score), ]
  }))
}))

# Save results
write.csv(res_max,    file.path("results", "max_pur.csv"),    row.names = FALSE)
write.csv(res_max_ct, file.path("results", "max_ct_pur.csv"), row.names = FALSE)
