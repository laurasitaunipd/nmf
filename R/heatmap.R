# heatmap
setwd("/Users/laura/Desktop/post/2 nmf/R")

library(pheatmap)
library(gridExtra)
library(grid)
library(RColorBrewer)

# crea la cartella di output se non esiste
if (!dir.exists("images")) dir.create("images")

# paletta colori
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

# heatmap NMF
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
  main = "NMF loadings",
  angle_col = 45,
  show_rownames = TRUE,
  show_colnames = TRUE,
  labels_row = rownames(NMF_plot),
  labels_col = colnames(NMF_plot),
  treeheight_row = 0,
  treeheight_col = 0
)

# heatmap EFA
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
  main = "EFA loadings",
  angle_col = 45,
  show_rownames = TRUE,
  show_colnames = TRUE,
  labels_row = rownames(EFA_plot),
  labels_col = colnames(EFA_plot),
  treeheight_row = 0,
  treeheight_col = 0
)

# --- Sezione finale stabile per macOS (senza Cairo) ---
# apre una finestra grafica Quartz
quartz(width = 14, height = 9)

# disegna le due heatmap affiancate
grid.newpage()
grid.arrange(p1[[4]], p2[[4]], ncol = 2, widths = c(1, 1))

# copia la grafica nella cartella images come PNG
dev.copy(png, filename = "images/heatmaps_EFA-NMF.png", width = 1400, height = 900, res = 150)
dev.off()  # chiude il dispositivo PNG
dev.off()  # chiude la finestra Quartz

