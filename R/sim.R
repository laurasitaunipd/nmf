library(lavaan)
library(MASS)
library(bayesNMF)
library(ragg)
library(RColorBrewer)
library(corrplot)
library(psych)
setwd("/Users/laura/Desktop/post/2 nmf/R")

# solo valori positivi con approssimazione alla normale
# in modo plausibile a livello biologico/psicologico

N = 1000
r = .20
x = mvrnorm(N, mu=c(-1,-1), Sigma=matrix(c(1,r,
                                         r,1),2,2))

Depressione = rbinom(N,50,plogis(x[,1]))
Ansia = rbinom(N,50,plogis(x[,2]))
cor(Depressione,Ansia)

Depr_1 = 1.2*Depressione+rbinom(N,50,.7)
Depr_2 = 0.5*Depressione+rbinom(N,80,.5)
Depr_3 = 0.8*Depressione+rbinom(N,40,.5)
Depr_4 = 1.5*Depressione+rbinom(N,30,.3)
Depr_5 = 1.0*Depressione+rbinom(N,100,.1)

Anx_1 = 1.2*Ansia+rbinom(N,50,.7)
Anx_2 = 0.5*Ansia+rbinom(N,80,.5)
Anx_3 = 0.8*Ansia+rbinom(N,40,.5)
Anx_4 = 1.5*Ansia+rbinom(N,30,.3)
Anx_5 = 1.0*Ansia+rbinom(N,100,.1)

df = data.frame(
  ID = 1:N,
  Depr_1, Depr_2, Depr_3, Depr_4, Depr_5,
  Anx_1, Anx_2, Anx_3, Anx_4, Anx_5
)


model = "
depr_latent =~ Depr_1+Depr_2+Depr_3+Depr_4+Depr_5
anx_latent =~ Anx_1+Anx_2+Anx_3+Anx_4+Anx_5
"
fit = cfa(model, df, std.lv=T)
summary(fit,standardized=T)
# pred_depress = as.numeric( predict(fit) )
# plot(pred_depress, Depressione)
# cor(pred_depress, Depressione)

######################################################################## 
# efa 
efa <- fa(df[, 2:ncol(df)],
          nfactors = 2,
          fm       = "minres",
          rotate   = "varimax",
          scores   = "regression")
quartz(); biplot.psych(efa, choose = c(1, 2))

# ======================================
# NMF (Bayesian NMF)
# ======================================

if (!dir.exists("output_sim")) dir.create("output_sim")

file_result_sim <- "output_sim/result_bayesNMF_sim.rds"

if (file.exists(file_result_sim)) {
  message("✓ Loading saved result from file...")
  result_sim <- readRDS(file_result_sim)   
} else {
  message("Running bayesNMF()")
  
  # La prima colonna è l'ID, quindi prendo solo le colonne 2: fine
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

# Calcolo media tra lower e upper se presenti
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

# CONFRONTO: corplot efa cfa nmf
# loadings EFA
efa_loadingsSIMSIMSIMSIM <- as.matrix(efa$loadings[, 1:2]) # MATRIX item x factor
colnames(efa_loadingsSIMSIMSIMSIM) <- paste0("efa_factor", 1:ncol(efa_loadingsSIMSIMSIMSIM))

# loadings NMF (P2)
nmf_loadingsSIM <- round(Psim, 2)  # matrice item × 1
rownames(nmf_loadingsSIM) <- paste0("item_", 1:nrow(nmf_loadingsSIM))
colnames(nmf_loadingsSIM) <- paste0("nmf_factor", 1:ncol(nmf_loadingsSIM))

allLoadings <- cbind(nmf_loadingsSIM, efa_loadingsSIMSIMSIMSIM)

print(allLoadings)

# CONFRONTO: corplot efa, cfa, nmf
# loadings CFA
cfa_loadingsSIM <- lavaan::inspect(fit, "std")$lambda
round(cfa_loadingsSIM, 2)
rownames(cfa_loadingsSIM) <- names(df)[-1]
colnames(cfa_loadingsSIM) <- paste0("cfa_factor", 1:ncol(cfa_loadingsSIM))

allLoadings2 <- cbind(nmf_loadingsSIM, efa_loadingsSIMSIMSIMSIM, cfa_loadingsSIM)

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

