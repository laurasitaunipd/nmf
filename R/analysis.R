library(readxl)
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

# BESSI questionnaire 
data <- read_excel("data.xlsx")
item_names = paste0("bessi_", 1:45)
item_data <- data[, item_names]

item_data <- data.frame(lapply(item_data, as.numeric))
item_data_complete <- na.omit(item_data)
str(item_data_complete)

# 5 factors
risultato_fa <- factanal(item_data_complete, factors = 5, rotation = "promax", scores = "regression")

print(risultato_fa, digits = 2, cutoff = 0.3)
risultato_fa$loadings
# prop var. (influence of each factor / total variance of data)%
# cum var. our 5 factors explain 38% of the total variance

######################################################################## NMF
# save results in the outpiut folder
if (!dir.exists("output")) dir.create("output")

file_result <- "output/result_bayesNMF.rds"

if (file.exists(file_result)) {
  message("✓ Loading saved result from file...")
  result <- readRDS(file_result)
} else {
  message("Running bayesNMF()")
  
  M_t <- t(as.matrix(item_data_complete))
  dim(M_t)
  
  result <- bayesNMF(
    data = M_t,
    likelihood = "normal",
    prior = "truncnormal",
    rank = 5
  )
  
  saveRDS(result, file_result)
  message("Result saved in: ", file_result)
}

MAP <- result$get_MAP()

# extract the two matrices (P and E)
P <- MAP$P
E <- MAP$E

# avarage value of the upper and lower bound
if (is.list(P)) {
  if (all(c("lower", "upper") %in% names(P))) P <- (P$lower + P$upper) / 2
  else P <- as.matrix(P[[1]])
}
if (is.list(E)) {
  if (all(c("lower", "upper") %in% names(E))) E <- (E$lower + E$upper) / 2
  else E <- as.matrix(E[[1]])
}

P <- as.matrix(P)
E <- as.matrix(E)

######################################################################## EFA VS NMF
# loadings EFA
efa_loadings <- as.matrix(risultato_fa$loadings[, 1:5]) # MATRIX item x factor
rownames(efa_loadings) <- paste0("item_", 1:nrow(efa_loadings))
colnames(efa_loadings) <- paste0("EFA-factor", 1:ncol(efa_loadings))

# loadings nmf (P)
nmf_loadings <- P # MATRIX item x factor
rownames(nmf_loadings) <- paste0("item_", 1:nrow(nmf_loadings))
colnames(nmf_loadings) <- paste0("NMF-factor", 1:ncol(nmf_loadings))

######################################################################## 

# correspondences factor scores individuals

tE = t(E)
colnames(tE) = paste0("nmf_scores_",1:5)
efa_scores = risultato_fa$scores

allScores = cbind(tE,efa_scores)
round(cor(allScores),3)

quartz()
corrplot(
  cor(allScores),
  method = "color",   # fill cells with colors
  type = "full",      # show full matrix
  addCoef.col = "black",   # show correlation coefficients
  number.cex = 0.7,        # text size for numbers
  tl.cex = 0.8             # text size for labels
)

# correspondences FACTOR LOADINGS !

allLoadings = cbind(nmf_loadings,efa_loadings)

agg_png("images/nmf-efa.png", width = 1400, height = 900, res = 150)
corrplot(
  cor(allLoadings),
  method = "color",   # fill cells with colors
  type = "full",      # show full matrix
  addCoef.col = "black",   # show correlation coefficients
  number.cex = 0.7,        # text size for numbers
  tl.cex = 0.8             # text size for labels
)
dev.off()

########################################################################

plot(P[,3],efa_loadings[,5])

efa_loadings_diff <- efa_loadings[,4]-efa_loadings[,5]-efa_loadings[,3]
plot(P[,5],efa_loadings_diff)
cor(P[,5],efa_loadings_diff)

######################################################################## 
# CFA 
######################################################################## 

SMD = paste0("SMD =~ ",paste0("bessi_",seq(1,41,5),collapse=" + "))
IND = paste0("IND =~ ",paste0("bessi_",seq(5,45,5),collapse=" + "))
COD = paste0("COD =~ ",paste0("bessi_",seq(3,43,5),collapse=" + "))
SED = paste0("SED =~ ",paste0("bessi_",seq(2,42,5),collapse=" + "))
ESD = paste0("ESD =~ ",paste0("bessi_",seq(4,44,5),collapse=" + "))
model = paste0(c(SMD, IND, COD, SED, ESD),collapse=" \n ")

####

model = "
SMD =~ bessi_1 + bessi_6 + bessi_11 + bessi_16 + bessi_21 + bessi_26 + bessi_31 + bessi_36 + bessi_41 
IND =~ bessi_5 + bessi_10 + bessi_15 + bessi_20 + bessi_25 + bessi_30 + bessi_35 + bessi_40 + bessi_45
COD =~ bessi_3 + bessi_8 + bessi_13 + bessi_18 + bessi_23 + bessi_28 + bessi_33 + bessi_38 + bessi_43 
SED =~ bessi_2 + bessi_7 + bessi_12 + bessi_17 + bessi_22 + bessi_27 + bessi_32 + bessi_37 + bessi_42 
ESD =~ bessi_4 + bessi_9 + bessi_14 + bessi_19 + bessi_24 + bessi_29 + bessi_34 + bessi_39 + bessi_44
"
fit = cfa(model=model, data=data, ordered=T)
summary(fit, standardized=T)
fitMeasures(fit, fit.measures=c("rmsea","srmr","cfi","nnfi"))

modificationIndices(fit, sort.=T)[1:10,]

####

lambda_unstd = inspect(fit, what = "coef")$lambda
lambda_unstd

lambda = inspect(fit, what = "std")$lambda
lambda

######################################################################## 
# CFA vs NMF
# correspondences factor CFA loadings
rownames(lambda) <- paste0("item_", 1:nrow(lambda))
lambda <- lambda[order(as.numeric(gsub("item_", "", rownames(lambda)))), ]

allLoadings = cbind(P,lambda)

agg_png("images/nmf-cfa.png", width = 1400, height = 900, res = 150)
corrplot(
  cor(allLoadings),
  method = "color",   # fill cells with colors
  type = "full",      # show full matrix
  addCoef.col = "black",   # show correlation coefficients
  number.cex = 0.7,        # text size for numbers
  tl.cex = 0.8             # text size for labels
)
dev.off()

########################################################################

allLoadings3 = cbind(nmf_loadings,efa_loadings,lambda)

agg_png("images/nmf-efa-cfa.png", width = 1400, height = 900, res = 150)
corrplot(
  cor(allLoadings3),
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
nmf_norm <- sweep(nmf_loadings, 2, apply(nmf_loadings, 2, max, na.rm = TRUE), "/")
main_factor_nmf <- apply(nmf_norm, 1, which.max)
order_items_nmf <- order(main_factor_nmf, -apply(nmf_norm, 1, max))
NMF_plot <- nmf_norm[order_items_nmf, , drop = FALSE]

# normalize efa_loadings and sort by leading factor
efa_norm <- sweep(efa_loadings, 2, apply(efa_loadings, 2, max, na.rm = TRUE), "/")
main_factor_efa <- apply(efa_norm, 1, which.max)
order_items_efa <- order(main_factor_efa, -apply(efa_norm, 1, max))
EFA_plot <- efa_norm[order_items_efa, , drop = FALSE]

# normalize efa_loadings and sort by leading factor
lambda_norm <- sweep(lambda, 2, apply(lambda, 2, max, na.rm = TRUE), "/")
main_factor_cfa <- apply(lambda_norm, 1, which.max)
order_items_cfa <- order(main_factor_cfa, -apply(lambda_norm, 1, max))
CFA_plot <- lambda_norm[order_items_cfa, , drop = FALSE]

p1 <- pheatmap( NMF_plot, color = pal, cluster_rows = FALSE, cluster_cols = FALSE, scale = "none", border_color = "black", legend = TRUE, fontsize_row = 7, fontsize_col = 10, cellwidth = 20, cellheight = 6, main = "NMF loadings", angle_col = 45, show_rownames = TRUE, show_colnames = TRUE, labels_row = rownames(NMF_plot), labels_col = colnames(NMF_plot), treeheight_row = 0, treeheight_col = 0 )
p2 <- pheatmap( EFA_plot, color = pal, cluster_rows = FALSE, cluster_cols = FALSE, scale = "none", border_color = "black", legend = TRUE, fontsize_row = 7, fontsize_col = 10, cellwidth = 20, cellheight = 6, main = "EFA loadings", angle_col = 45, show_rownames = TRUE, show_colnames = TRUE, labels_row = rownames(EFA_plot), labels_col = colnames(EFA_plot), treeheight_row = 0, treeheight_col = 0 )
p3 <- pheatmap(
  CFA_plot,
  color = pal,
  cluster_rows = FALSE, cluster_cols = FALSE,
  scale = "none", border_color = "black",
  legend = TRUE, fontsize_row = 7, fontsize_col = 10,
  cellwidth = 20, cellheight = 6,
  main = "CFA loadings", angle_col = 45,
  show_rownames = TRUE, show_colnames = TRUE,
  labels_row = rownames(CFA_plot), labels_col = colnames(CFA_plot),
  treeheight_row = 0, treeheight_col = 0
)

agg_png("images/HEATMAPnmf-efa.png", width = 1400, height = 900, res = 150)
grid.arrange(p1[[4]], p2[[4]], ncol = 2)
dev.off()

agg_png("images/HEATMAPnmf-efa-cfa.png", width = 2000, height = 900, res = 150)
grid.arrange(p1[[4]], p2[[4]], p3[[4]], ncol = 3)
dev.off()

