library(pheatmap)
library(gridExtra)
library(grid)
library(RColorBrewer)
setwd("/Users/laura/Desktop/post/2 nmf/R")

######################################################################## 
# DATI TOMMASO # HEATMAP COMPARISON: EFA vs bayesNMF

pal <- colorRampPalette(c("navy", "white", "firebrick3"))(200)

# normalizza nmf_loadings e ordina per fattore principale
nmf_norm <- sweep(nmf_loadings, 2, apply(nmf_loadings, 2, max, na.rm = TRUE), "/")
main_factor_nmf <- apply(nmf_norm, 1, which.max)
order_items_nmf <- order(main_factor_nmf, -apply(nmf_norm, 1, max))
NMF_plot <- nmf_norm[order_items_nmf, , drop = FALSE]

# normalizza loadings EFA e ordina
efa_norm <- sweep(efa_loadings, 2, apply(efa_loadings, 2, max, na.rm = TRUE), "/")
main_factor_efa <- apply(efa_norm, 1, which.max)
order_items_efa <- order(main_factor_efa, -apply(efa_norm, 1, max))
EFA_plot <- efa_norm[order_items_efa, , drop = FALSE]


p1 <- pheatmap( NMF_plot, color = pal, cluster_rows = FALSE, cluster_cols = FALSE, scale = "none", border_color = "black", legend = TRUE, fontsize_row = 7, fontsize_col = 10, cellwidth = 20, cellheight = 6, main = "NMF loadings", angle_col = 45, show_rownames = TRUE, show_colnames = TRUE, labels_row = rownames(NMF_plot), labels_col = colnames(NMF_plot), treeheight_row = 0, treeheight_col = 0 )
p2 <- pheatmap( EFA_plot, color = pal, cluster_rows = FALSE, cluster_cols = FALSE, scale = "none", border_color = "black", legend = TRUE, fontsize_row = 7, fontsize_col = 10, cellwidth = 20, cellheight = 6, main = "EFA loadings", angle_col = 45, show_rownames = TRUE, show_colnames = TRUE, labels_row = rownames(EFA_plot), labels_col = colnames(EFA_plot), treeheight_row = 0, treeheight_col = 0 )

agg_png("images/heatmaps_EFA-NMF.png", width = 1400, height = 900, res = 150)
grid.arrange(p1[[4]], p2[[4]], ncol = 2)
dev.off()


######################################################################## 
# WORK IN PROGRESS
# SIMULAZIONE # HEATMAP COMPARISON: EFA vs bayesNMF

# Normalizza e ordina NMF
nmf_loadings_sim <- Psim  # matrice NMF normalizzata
nmf_norm <- sweep(nmf_loadings_sim, 2, apply(nmf_loadings_sim, 2, max, na.rm = TRUE), "/")
main_factor_nmf <- apply(nmf_norm, 1, which.max)
order_items_nmf <- order(main_factor_nmf, -apply(nmf_norm, 1, max))
NMF_plot <- nmf_norm[order_items_nmf, , drop = FALSE]


# Normalizza e ordina EFA
efa_loadings_sim <- as.matrix(efa$loadings)
storage.mode(efa_loadings_sim) <- "numeric"
efa_norm <- sweep(efa_loadings_sim, 2, apply(efa_loadings_sim, 2, max, na.rm = TRUE), "/")
main_factor_efa <- apply(efa_norm, 1, which.max)
order_items_efa <- order(main_factor_efa, -apply(efa_norm, 1, max))
EFA_plot <- efa_norm[order_items_efa, , drop = FALSE]

agg_png("images/heatmapsSIM_EFA-NMF.png", width = 1400, height = 900, res = 150)
par(mfrow = c(1, 2))
p1 <- pheatmap(
  NMF_plot,
  color = pal,
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  scale = "none",
  border_color = "black",
  legend = TRUE,
  fontsize_row = 7,
  fontsize_col = 10,
  cellwidth = 20,
  cellheight = 6,
  main = "NMF (simulazione)",
  angle_col = 45,
  show_rownames = TRUE,
  show_colnames = TRUE,
  labels_row = rownames(NMF_plot),
  labels_col = colnames(NMF_plot),
  treeheight_row = 0,
  treeheight_col = 0
)

p2 <- pheatmap(
  EFA_plot,
  color = pal,
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  scale = "none",
  border_color = "black",
  legend = TRUE,
  fontsize_row = 7,
  fontsize_col = 10,
  cellwidth = 20,
  cellheight = 6,
  main = "EFA (simulazione)",
  angle_col = 45,
  show_rownames = TRUE,
  show_colnames = TRUE,
  labels_row = rownames(EFA_plot),
  labels_col = colnames(EFA_plot),
  treeheight_row = 0,
  treeheight_col = 0
)

agg_png("images/heatmapSIM_EFA-NMF.png", width = 1400, height = 900, res = 150)
grid.arrange(p1[[4]], p2[[4]], ncol = 2)
dev.off()
