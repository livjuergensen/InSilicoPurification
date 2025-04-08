# -----------------------------------------------------------------------------
# Visualization of classification results for adjusted methylation profiles
# -----------------------------------------------------------------------------

dir.create("plots")

# Set to TRUE to plot scores for methylation class (MC), methylation class family (MCF), or both
plot_mcf <- TRUE
plot_mc  <- TRUE

colors <- readRDS(file.path("data", "colors_pur.rds"))

# General plotting function
plot_prob <- function(prob, color_vector, output_suffix) {
  prob <- as.data.frame(prob)
  prob$sample             <- sub("^[^_]*_[^_]*_(.*)$",  "\\1", rownames(prob))
  prob$subtracted_percent <- sub("^[^_]*_(.*?)_.*$",    "\\1", rownames(prob))
  prob$cell_type          <- sub("^(.*?)_.*$",           "\\1", rownames(prob))
  
  keep_cols <- c("sample", "subtracted_percent", "cell_type")
  df_long <- reshape(
    data      = prob,
    varying   = setdiff(names(prob), keep_cols),
    v.names   = "Probability",
    timevar   = "Outcome",
    times     = setdiff(names(prob), keep_cols),
    direction = "long"
  )
  
  df_long$Outcome[df_long$Probability < 0.15] <- "OTHER"
  
  for (smpl in unique(df_long$sample)) {
    df_sample   <- df_long[df_long$sample == smpl, ]
    all_sp      <- unique(df_sample$subtracted_percent)
    cell_types  <- unique(df_sample$cell_type)
    n_ct        <- length(cell_types)
    
    pdf(file = file.path("plots", paste0(smpl, "_purified_", output_suffix, ".pdf")), 
        width = 8.27, height = 11.69)
    
    par(
      mfrow    = c(n_ct, 1),
      oma      = c(5, 5, 4, 12),
      mar      = c(1, 1, 2, 1),
      xaxs     = "i",
      cex.main = 1,
      cex.lab  = 1,
      cex.axis = 0.9
    )
    
    for (i in seq_along(cell_types)) {
      ct    <- cell_types[i]
      df_ct <- df_sample[df_sample$cell_type == ct, ]
      
      sumProbs <- tapply(df_ct$Probability, df_ct$Outcome, sum)
      sumProbs[is.na(sumProbs)] <- 0
      reversed <- c(setdiff(names(sort(sumProbs, decreasing = TRUE)), "OTHER"), "OTHER")
      
      mat_ct <- tapply(
        df_ct$Probability,
        INDEX = list(df_ct$Outcome, df_ct$subtracted_percent),
        FUN   = sum
      )
      
      mat_ct <- mat_ct[, all_sp, drop = FALSE]
      mat_ct[is.na(mat_ct)] <- 0
      mat_ct <- mat_ct[reversed, , drop = FALSE]
      
      bar_colors <- rep("gray70", nrow(mat_ct))
      matched    <- intersect(reversed, names(color_vector))
      bar_colors[match(matched, reversed)] <- color_vector[matched]
      
      x_mid <- barplot(
        mat_ct,
        col       = bar_colors,
        border    = NA,
        space     = 0,
        axes      = FALSE,
        axisnames = FALSE
      )
      
      abline(h = 0.9, lty = 2)
      axis(side = 2, at = c(0, 0.5, 1), labels = c("0", "0.5", "1"), las = 1)
      
      lab_idx <- c(seq(1, length(x_mid), by = 10), length(x_mid))
      axis(side = 1, at = x_mid[lab_idx], pos = 0, labels = FALSE)
      if (i == n_ct) {
        axis(side = 1, at = x_mid[lab_idx], pos = 0, 
             labels = colnames(mat_ct)[lab_idx], las = 1)
      }
      
      mtext(ct, side = 3, line = 0.5, cex = 0.9, font = 2)
    }
    
    mtext(smpl, side = 3, outer = TRUE, line = 1, cex = 1, font = 2)
    mtext("Removed fraction of non-malignant cell type", side = 1, outer = TRUE, line = 1.7, cex = 0.9)
    mtext("Calibrated Probability Score", side = 2, outer = TRUE, line = 1.7, cex = 0.9)
    
    par(xpd = NA)
    out_smpl    <- intersect(names(color_vector), unique(df_sample$Outcome))
    legend_cols <- color_vector[out_smpl]
    y_mid       <- grconvertY(0.5, from = "ndc", to = "user")
    
    legend(
      x      = ncol(mat_ct) + 1,
      y      = y_mid,
      legend = out_smpl,
      fill   = legend_cols,
      bty    = "n",
      cex    = 1,
      yjust  = 0.5,
      xjust  = 0,
      title  = "Methylation Class Family",
      title.cex = 1
    )
    
    dev.off()
  }
}

# Generate plots
if (plot_mcf) {
  prob_mcf <- readRDS(file.path("results", "prob_mcf.RDS"))
  plot_prob(prob_mcf, colors, "mcf")
}

if (plot_mc) {
  prob <- readRDS(file.path("results", "prob.RDS"))
  plot_prob(prob, colors, "mc")
}
