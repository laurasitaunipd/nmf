library(bayesNMF)
library(corrplot)
library(lavaan)
library(ragg)
library(pheatmap)
library(gridExtra)
library(grid)
library(RColorBrewer)
setwd("/Users/laura/Desktop/post/2 nmf/R")

######################################################################## 
# EFA
######################################################################## 

# PUNTEGGI GREZZI (conteggi risposte corrette) INTELLIGENZA IN BAMBINI CON DSA 
df <- read.csv("df.csv", header = TRUE)

# selezione dei 10 subtest essenziali
df_w <- df[, 7:16]
colnames(df_w) <- paste0("item_", seq_len(ncol(df_w)))

df_w <- df_w[complete.cases(df_w), ] # rimuove soggetti con almeno un NA
df_w <- as.data.frame(lapply(df_w, as.numeric))
range(df_w, na.rm = TRUE)

hist(as.numeric(as.matrix(df_w)),
     main = "Distribuzione punteggi grezzi Wechsler",
     xlab = "Punteggio grezzo")

# rendiamo la matrice numerica intera
df_w <- as.matrix(df_w)
df_w <- round(df_w)
storage.mode(df_w) <- "numeric"
any(df_w %% 1 != 0)

# so che ci sono 4 fattori specifici sottostanti
efa2 <- factanal(df_w, factors = 4, rotation = "promax", scores = "regression")

print(efa2, digits = 2, cutoff = 0.3)
efa2$loadings
# prop var. (influence of each factor / total variance of data)%
# cum var. our 4 factors explain 68% of the total variance

######################################################################## NMF
# save results in the outpiut folder
if (!dir.exists("output2")) dir.create("output2")

file_result2 <- "output2/result_bayesNMF2.rds"

if (file.exists(file_result2)) {
  message("✓ Loading saved result from file...")
  result2 <- readRDS(file_result2)
} else {
  message("Running bayesNMF()")
  
  M_t2 <- t(as.matrix(df_w))
  dim(M_t2)
  
  nmf2 <- bayesNMF(
    data = M_t2,
    likelihood = "poisson", # dati = conteggi discreti 
    prior = "truncnormal", # positivi (0 – 148)
    rank = 4
  )
  
  saveRDS(nmf2, file_result2)
  message("Result saved in: ", file_result2)
}

MAP2 <- nmf2$get_MAP()

# extract the two matrices (P and E)
P2 <- MAP2$P
E2 <- MAP2$E

# avarage value of the upper and lower bound
if (is.list(P2)) {
  if (all(c("lower", "upper") %in% names(P2))) P2 <- (P2$lower + P2$upper) / 2
  else P2 <- as.matrix(P2[[1]])
}
if (is.list(E2)) {
  if (all(c("lower", "upper") %in% names(E2))) E2 <- (E2$lower + E2$upper) / 2
  else E2 <- as.matrix(E2[[1]])
}

P22 <- as.matrix(P2)
E22 <- as.matrix(E2)

######################################################################## EFA VS NMF
# loadings EFA
efa2_loadings <- as.matrix(efa2$loadings[, 1:4]) # MATRIX item x factor
rownames(efa2_loadings) <- paste0("item_", 1:nrow(efa2_loadings))
colnames(efa2_loadings) <- paste0("EFA-factor", 1:ncol(efa2_loadings))

# loadings nmf (P)
nmf2_loadings <- P22 # MATRIX item x factor
rownames(nmf2_loadings) <- paste0("item_", 1:nrow(nmf2_loadings))
colnames(nmf2_loadings) <- paste0("NMF-factor", 1:ncol(nmf2_loadings))

######################################################################## 
# correspondences factor scores individuals
tE2 = t(E22)
colnames(tE2) = paste0("nmf-FACTORscores",1:4)
nmf2_scores = tE2 # MATRIX unit x factor

efa2_scores = efa2$scores # MATRIX unit x factor
colnames(efa2_scores) = paste0("efa-FACTORscores",1:4)

# correspondences FACTOR LOADINGS !

allLoadings2 = cbind(nmf2_loadings,efa2_loadings)

agg_png("images/2nmf-efa.png", width = 1400, height = 900, res = 150)
corrplot(
  cor(allLoadings2),
  method = "color",   # fill cells with colors
  type = "full",      # show full matrix
  addCoef.col = "black",   # show correlation coefficients
  number.cex = 0.7,        # text size for numbers
  tl.cex = 0.8             # text size for labels
)
dev.off()

# CFA 
######################################################################## 

# costruisci il modello 

fit2 = cfa(model = model2, data = df_w, ordered = T)
summary(fit2, standardized = T)
fitMeasures(fit2, fit.measures=c("rmsea","srmr","cfi","nnfi"))

modificationIndices(fit2, sort.=T)[1:10,]

####

lambda_unstd2 = inspect(fit2, what = "coef")$lambda
lambda_unstd2

lambda2 = inspect(fit2, what = "std")$lambda
lambda2

corrplot(cor(efa2_loadings,lambda2))

###

cfa2_scores = predict(fit2) # MATRIX unit x factor

# correspondences factor CFA loadings
rownames(lambda2) <- paste0("item_", 1:nrow(lambda2))
lambda2 <- lambda2[order(as.numeric(gsub("item_", "", rownames(lambda2)))), ]

########################################################################

allLoadings32 = cbind(nmf2_loadings,efa2_loadings,lambda2)

agg_png("images/23nmf-efa-cfa.png", width = 1400, height = 900, res = 150)
corrplot(
  cor(allLoadings32),
  method = "color",   # fill cells with colors
  type = "full",      # show full matrix
  addCoef.col = "black",   # show correlation coefficients
  number.cex = 0.7,        # text size for numbers
  tl.cex = 0.8             # text size for labels
)
dev.off()

######################################################################## 
# HEATMAP COMPARISON: EFA vs bayesNMF

pal <- colorRampPalette(c("navy", "white", "firebrick3"))(200)

# normalize nmf_loadings and sort by leading factor
nmf2_norm <- sweep(nmf2_loadings, 2, apply(nmf2_loadings, 2, max, na.rm = TRUE), "/")
main_factor_nmf2 <- apply(nmf2_norm, 1, which.max)
order_items_nmf2 <- order(main_factor_nmf2, -apply(nmf2_norm, 1, max))
NMF_plot2 <- nmf2_norm[order_items_nmf2, , drop = FALSE]

# normalize efa_loadings and sort by leading factor
efa2_norm <- sweep(efa2_loadings, 2, apply(efa2_loadings, 2, max, na.rm = TRUE), "/")
main_factor_efa2 <- apply(efa2_norm, 1, which.max)
order_items_efa2 <- order(main_factor_efa2, -apply(efa2_norm, 1, max))
EFA_plot2 <- efa2_norm[order_items_efa2, , drop = FALSE]

# normalize efa_loadings and sort by leading factor
lambda2_norm <- sweep(lambda2, 2, apply(lambda2, 2, max, na.rm = TRUE), "/")
main_factor_cfa2 <- apply(lambda2_norm, 1, which.max)
order_items_cfa2 <- order(main_factor_cfa2, -apply(lambda2_norm, 1, max))
CFA_plot2 <- lambda2_norm[order_items_cfa2, , drop = FALSE]

p1 <- pheatmap( NMF_plot2, color = pal, cluster_rows = FALSE, cluster_cols = FALSE, scale = "none", border_color = "black", legend = TRUE, fontsize_row = 7, fontsize_col = 10, cellwidth = 20, cellheight = 6, main = "NMF loadings", angle_col = 45, show_rownames = TRUE, show_colnames = TRUE, labels_row = rownames(NMF_plot2), labels_col = colnames(NMF_plot2), treeheight_row = 0, treeheight_col = 0 )
p2 <- pheatmap( EFA_plot2, color = pal, cluster_rows = FALSE, cluster_cols = FALSE, scale = "none", border_color = "black", legend = TRUE, fontsize_row = 7, fontsize_col = 10, cellwidth = 20, cellheight = 6, main = "EFA loadings", angle_col = 45, show_rownames = TRUE, show_colnames = TRUE, labels_row = rownames(EFA_plot2), labels_col = colnames(EFA_plot2), treeheight_row = 0, treeheight_col = 0 )
p3 <- pheatmap(
  CFA_plot2,
  color = pal,
  cluster_rows = FALSE, cluster_cols = FALSE,
  scale = "none", border_color = "black",
  legend = TRUE, fontsize_row = 7, fontsize_col = 10,
  cellwidth = 20, cellheight = 6,
  main = "CFA loadings", angle_col = 45,
  show_rownames = TRUE, show_colnames = TRUE,
  labels_row = rownames(CFA_plot2), labels_col = colnames(CFA_plot2),
  treeheight_row = 0, treeheight_col = 0
)

agg_png("images/2HEATMAPnmf-efa.png", width = 1400, height = 900, res = 150)
grid.arrange(p1[[4]], p2[[4]], ncol = 2)
dev.off()

agg_png("images/2HEATMAPnmf-efa-cfa.png", width = 2000, height = 900, res = 150)
grid.arrange(p1[[4]], p2[[4]], p3[[4]], ncol = 3)
dev.off()

######################################################################## 
# FACTOR SCORES comparison

allScores = cbind(efa_scores,nmf_scores,cfa_scores)

agg_png("images/SCORESnmf-efa-cfa.png", width = 1400, height = 900, res = 150)
corrplot(
  cor(allScores),
  method = "color",   # fill cells with colors
  type = "full",      # show full matrix
  addCoef.col = "black",   # show correlation coefficients
  number.cex = 0.7,        # text size for numbers
  tl.cex = 0.8             # text size for labels
)
dev.off()