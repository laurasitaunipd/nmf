library(lavaan)
library(MASS)
library(bayesNMF)
library(ragg)
library(corrplot)
library(psych)
setwd("/Users/laura/Desktop/post/2 nmf/R")

N = 1000
r = .20
x = mvrnorm(N, mu=c(-1,-1), Sigma=matrix(c(1,r,
                                           r,1),2,2))

Depression = rbinom(N,50,plogis(x[,1]))
Anxiety = rbinom(N,50,plogis(x[,2]))
cor(Depression,Anxiety)

Depr_1 = 1.2*Depression+rbinom(N,50,.7)
Depr_2 = 0.5*Depression+rbinom(N,80,.5)
Depr_3 = 0.8*Depression+rbinom(N,40,.5)
Depr_4 = 1.5*Depression+rbinom(N,30,.3)
Depr_5 = 1.0*Depression+rbinom(N,100,.1)

Anx_1 = 1.2*Anxiety+rbinom(N,50,.7)
Anx_2 = 0.5*Anxiety+rbinom(N,80,.5)
Anx_3 = 0.8*Anxiety+rbinom(N,40,.5)
Anx_4 = 1.5*Anxiety+rbinom(N,30,.3)
Anx_5 = 1.0*Anxiety+rbinom(N,100,.1)

df = data.frame(
  ID = 1:N,
  Depr_1, Depr_2, Depr_3, Depr_4, Depr_5,
  Anx_1, Anx_2, Anx_3, Anx_4, Anx_5
)


model = "
depr_latent =~ Depr_1+Depr_2+Depr_3+Depr_4+Depr_5
anx_latent =~ Anx_1+Anx_2+Anx_3+Anx_4+Anx_5
"
SIMcfa = cfa(model, df, std.lv=T)
summary(SIMcfa,standardized=T)
SIMcfa_scores = predict(SIMcfa)
# plot(pred_depress, Depression)
# cor(pred_depress, Depression)

######################################################################## 
# efa 
SIMefa <- fa(df[, 2:ncol(df)],
             nfactors = 2,
             fm       = "minres",
             rotate   = "promax",
             scores   = "regression")
agg_png("images/SIMbiplot.png", width = 1400, height = 900, res = 150)
biplot.psych(SIMefa, choose = c(1, 2))
dev.off()


# NMF

if (!dir.exists("output_sim")) dir.create("output_sim")

file_result_sim <- "output_sim/result_bayesNMF_sim.rds"

if (file.exists(file_result_sim)) {
  message("✓ Loading saved result from file...")
  result_sim <- readRDS(file_result_sim)   
} else {
  message("Running bayesNMF()")
  
  # The first column is the ID, so only take columns 2 to the end
  M_t_2 <- t(as.matrix(df[, -1]))  
  dim(M_t_2)
  
  result_sim <- bayesNMF(
    data = M_t_2,
    likelihood = "normal",
    prior = "truncnormal",
    rank = 2
  )
  
  saveRDS(result_sim, file_result_sim)  
  message("Result saved in: ", file_result_sim)
}

MAPsim <- result_sim$get_MAP()  

Psim <- MAPsim$P
Esim <- MAPsim$E

# Compute mean between lower and upper
if (is.list(Psim)) {
  if (all(c("lower", "upper") %in% names(Psim))) {
    Psim <- (Psim$lower + Psim$upper) / 2
  } else {
    Psim <- as.matrix(Psim[[1]])
  }
}

if (is.list(Esim)) {
  if (all(c("lower", "upper") %in% names(Esim))) {
    Esim <- (Esim$lower + Esim$upper) / 2
  } else {
    Esim <- as.matrix(Esim[[1]])
  }
}

Psim <- as.matrix(Psim)
Esim <- as.matrix(Esim)

# COMPARISON: corrplot efa cfa nmf
# factor scores individuals
SIMnmf_scores = t(Esim) # MATRIX unit x factor
colnames(SIMnmf_scores) = paste0("nmf-FACTORscores",1:2)

SIMefa_scores = SIMefa$scores # MATRIX unit x factor
colnames(SIMefa_scores) = paste0("efa-FACTORscores",1:2)


# EFA loadings
SIMefa_loadings <- as.matrix(SIMefa$loadings[, 1:2]) # MATRIX item x factor
rownames(SIMefa_loadings) <- paste0("item_", 1:nrow(SIMefa_loadings))
colnames(SIMefa_loadings) <- paste0("EFA-factor", 1:ncol(SIMefa_loadings))

# NMF loadings (Psim)
SIMnmf_loadings <- Psim # MATRIX item x factor
rownames(SIMnmf_loadings) <- paste0("item_", 1:nrow(SIMnmf_loadings))
colnames(SIMnmf_loadings) <- paste0("NMF-factor", 1:ncol(SIMnmf_loadings))

allLoadings <- cbind(SIMnmf_loadings, SIMefa_loadings)

print(allLoadings)

# COMPARISON: corrplot efa, cfa, nmf
# CFA loadings
SIMcfa_loadings <- lavaan::inspect(SIMcfa, "std")$lambda
rownames(SIMcfa_loadings) <- names(df)[-1]
#colnames(SIMcfa_loadings) <- paste0("CFA-factor", 1:ncol(SIMcfa_loadings))

allLoadings2 <- cbind(SIMnmf_loadings, SIMefa_loadings, SIMcfa_loadings)

print(allLoadings2)

# corrplot
agg_png("images/SIMnmf-efa-cfa.png", width = 1400, height = 900, res = 150)
par(mfrow = c(1, 2))
corrplot(cor(allLoadings, use = "pairwise.complete.obs"),
         method="color", type="full", addCoef.col="black",
         number.cex=0.9, tl.cex=1, mar=c(0,0,1,0), title="nmf vs efa")
corrplot(cor(allLoadings2, use = "pairwise.complete.obs"),
         method="color", type="full", addCoef.col="black",
         number.cex=0.9, tl.cex=1, mar=c(0,0,1,0), title="nmf vs efa vs cfa")
dev.off()

######################################################################## 
# FACTOR SCORES comparison
allScores3 = cbind(SIMefa_scores,SIMnmf_scores,SIMcfa_scores)

agg_png("images/SIM-SCORESnmf-efa-cfa.png", width = 1400, height = 900, res = 150)
corrplot(
  cor(allScores3),
  method = "color",   # fill cells with colors
  type = "full",      # show full matrix
  addCoef.col = "black",   # show correlation coefficients
  number.cex = 0.7,        # text size for numbers
  tl.cex = 0.8             # text size for labels
)
dev.off()
