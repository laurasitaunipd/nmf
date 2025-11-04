library(bayesNMF)
library(corrplot)
library(lavaan)
library(psych)
library(GPArotation)
setwd("/Users/laura/Desktop/post/2 nmf/R")

# 45 item (punteggio da 1 a 5) # 100 soggetti # 5 fattori 

######################################################################## 
# SOLUZIONE DELLA SIMULAZIONE: MATRICE VERA
######################################################################## 
set.seed(123)

# ---------------------------
# 1) Matrice di loadings 45x5
# ---------------------------
n_items <- 45; n_factors <- 5; n <- 100
item_names <- paste0("item_", 1:n_items)
factor_names <- paste0("factor_", 1:n_factors)

loadings <- matrix(0, nrow = n_items, ncol = n_factors)
for (f in 1:n_factors) {
  idx <- ((f-1)*9 + 1):(f*9)
  loadings[idx, f] <- pmax(0, rnorm(9, 0.70, 0.05))
  loadings[idx, setdiff(1:n_factors, f)] <- pmax(0, rnorm(9*4, 0.20, 0.05))
}
loadings[loadings > 1] <- 1
rownames(loadings) <- item_names
colnames(loadings) <- factor_names

# ------------------------------------
# 2) Genera punteggi latenti e rumore
# ------------------------------------
# fattori ~ N(0,1) indipendenti (puoi introdurre correlazioni se vuoi)
S <- matrix(rnorm(n * n_factors), nrow = n, ncol = n_factors)
colnames(S) <- factor_names
# rumore specifico dell’item
epsilon <- matrix(rnorm(n * n_items, sd = 0.30), nrow = n, ncol = n_items)

# ----------------------------------------
# 3) Dati continui dal modello: X = S %*% L'
# ----------------------------------------
X_cont <- S %*% t(loadings) + epsilon      # n x p
colnames(X_cont) <- item_names
rownames(X_cont) <- paste0("subj_", 1:n)

# ------------------------------------
# 4) Discretizza in Likert 1–5 per item
# ------------------------------------
datasim <- sapply(1:n_items, function(j){
  br <- quantile(X_cont[, j], probs = seq(0, 1, length.out = 6), na.rm = TRUE)
  as.integer(cut(X_cont[, j], breaks = br, include.lowest = TRUE, labels = 1:5))
})
datasim <- as.matrix(datasim)
colnames(datasim) <- item_names
rownames(datasim) <- paste0("subj_", 1:n)

# ===============================
# 1) Calcola la matrice policorica
# ===============================
rho_poly <- polychoric(datasim)$rho   # matrice 45x45 di correlazioni policoriche

# Controllo base (dimensioni e struttura)
dim(rho_poly)
round(rho_poly[1:5, 1:5], 2)  # prime 5x5 per vedere l'andamento

# ===============================
# 2) EFA con 5 fattori
# ===============================
efa <- fa(
  r = rho_poly,       # matrice di correlazioni policoriche
  nfactors = 5,       # numero di fattori latenti
  fm = "minres",      # metodo di estrazione (robusto)
  rotate = "oblimin"  # rotazione obliqua (fattori correlati)
)

# ===============================
# 3) Risultati principali
# ===============================

# Matrice dei loadings ruotati
print(efa$loadings, cutoff = 0.3)

# Varianza spiegata cumulativa 
efa$Vaccounted # 76%
# Correlazioni tra fattori (perché rotazione obliqua)
efa$r.scores

# loadings EFA
efa_loadingssim <- as.matrix(efa$loadings[, 1:5]) # MATRIX item x factor
efa_loadingssim <- round(efa_loadingssim, 2)
colnames(efa_loadingssim) = paste0("factor_",1:5)

######################################################################## 
# NMF su datasim: P matrice di output su DATASIM (sub x item)
######################################################################## 
########################################################################
# BAYESIAN NMF SUI DATI SIMULATI (datasim)
########################################################################
if (!dir.exists("output")) dir.create("output")
file_result_sim <- "output/NMF_sim.rds"

if (file.exists(file_result_sim)) {
  message("✓ Carico il risultato salvato da file...")
  nmf_sim <- readRDS(file_result_sim)
} else {
  message("Eseguo bayesNMF()")
  
  # Trasponi la matrice: NMF vuole item (righe) × soggetti (colonne)
  M_t <- t(as.matrix(datasim))
  dim(M_t)  # 45 x 100
  
  # Esegui il modello
  nmf_sim <- bayesNMF(
    data = M_t,
    likelihood = "normal",
    prior = "truncnormal",
    rank = 5
  )
  
  # Salva il risultato
  saveRDS(nmf_sim, file_result_sim)
  message("✓ Risultato salvato in: ", file_result_sim)
}

MAP <- nmf_sim$get_MAP()
P <- MAP$P
E <- MAP$E

# Funzione di utilità: porta qualsiasi output di P/E a matrice numerica
to_numeric_matrix <- function(M) {
  # Se è una lista con intervalli, prendo la media (midpoint)
  if (is.list(M) && all(c("lower","upper") %in% names(M))) {
    M <- (M$lower + M$upper) / 2
  } else if (is.list(M)) {
    M <- as.matrix(M[[1]])
  }
  # Assicuro il tipo matrice
  if (!is.matrix(M)) M <- as.matrix(M)
  # Forzo modalità numerica
  storage.mode(M) <- "numeric"
  # Se per qualche motivo è ancora non-numerica, vectorizzo e ricostruisco
  if (!is.numeric(M)) {
    M <- matrix(as.numeric(M), nrow = nrow(M), dimnames = dimnames(M))
  }
  M
}

# Applico ai tuoi oggetti
Psim <- to_numeric_matrix(P)
Esim <- to_numeric_matrix(E)

# loadings nmf (P)
Psim <- round(Psim, 2) # MATRIX item x factor
rownames(Psim) <- paste0("item_", 1:nrow(Psim))
colnames(Psim) <- paste0("fattore_", 1:ncol(Psim))

# Normalizzazione colonne di P (somma=1), con denominatore sicuro
den <- colSums(Psim)
den[den == 0 | !is.finite(den)] <- 1
Psim_norm <- sweep(Psim, 2, den, "/") # paragonabili alle loadings dell’EFA (che sono standardizzate).

######################################################################## 
# heatmap # vedi file a parte heatmap.R
######################################################################## 

######################################################################## 
# matrice di correlazione tra le matrici P ed EFA
######################################################################## 

allLoadings = cbind(Psim_norm,efa_loadingssim)

# apri una finestra grafica quartz (macOS nativa)
quartz(width = 10, height = 10)

# crea il grafico
corrplot(
  cor(allLoadings),
  method = "color",      # celle colorate
  type = "full",         # matrice completa
  addCoef.col = "black", # mostra i valori di correlazione
  number.cex = 0.7,      # grandezza numeri
  tl.cex = 0.8           # grandezza etichette
)

# salva automaticamente in PNG nella cartella images
dev.copy(png, filename = "images/corrplot_sim.png", width = 1200, height = 1200, res = 150)
dev.off()  # chiude il PNG
dev.off()  # chiude la finestra quartz
