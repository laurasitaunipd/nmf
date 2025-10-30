library(readxl)
library(bayesNMF)
library(corrplot)
library(lavaan)

######################################################################## EFA
# dati tommaso
dati <- read_excel("dati.xlsx")
item_names = paste0("bessi_", 1:45)
item_data <- dati[, item_names]

item_data <- data.frame(lapply(item_data, as.numeric))
item_data_complete <- na.omit(item_data)
str(item_data_complete)

# Esegui l’analisi fattoriale esplorativa
# (puoi modificare il numero di fattori da 3 a quello che vuoi esplorare)
risultato_fa <- factanal(item_data_complete, factors = 5, rotation = "promax", scores = "regression")

# visualizza i risultati principali
print(risultato_fa, digits = 2, cutoff = 0.3)

# visualizza i carichi fattoriali
risultato_fa$loadings
# prop var. peso di ogni fattore / varianza totale dei dati IN PERCENTUALE
# cum var. i 5 fattori spieganon il 38% della varianza totale

#risultato_fa$uniquenesses

######################################################################## NMF
# mi assicuro che la cartella "output" esista
if (!dir.exists("output")) dir.create("output")

# salvo result
file_result <- "output/result_bayesNMF.rds"

# se il file esiste, lo carico
if (file.exists(file_result)) {
  message("✓ Carico il risultato salvato da file...")
  result <- readRDS(file_result)
} else {
  message("Eseguo bayesNMF()")
  
  # calcolo solo se serve
  M_t <- t(as.matrix(item_data_complete))
  dim(M_t)
  
  result <- bayesNMF(
    data = M_t,
    likelihood = "normal",
    prior = "truncnormal",
    rank = 5
  )
  
  # salvo il risultato
  saveRDS(result, file_result)
  message("Result salvato in: ", file_result)
}

MAP <- result$get_MAP()


######################################################################## 
# estrai le matrici dei fattori (P ed E)
P <- MAP$P
E <- MAP$E

# Se sono liste con intervalli, prendi la media dei limiti
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

######################################################################## CONFRONTO
# loadings EFA
efa_loadings <- as.matrix(risultato_fa$loadings[, 1:5]) # matrice item x fattore

# loadings nmf (P)
nmf_loadings <- round(P, 2) # matrice item x fattore
rownames(nmf_loadings) <- paste0("bessi_", 1:nrow(nmf_loadings))
colnames(nmf_loadings) <- paste0("Fattore_", 1:ncol(nmf_loadings))

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

# correspondences factor loadings

colnames(P) = paste0("nmf_loadings_",1:5)
efa_loadings = risultato_fa$loadings

allLoadings = cbind(P,efa_loadings)

quartz()
corrplot(
  cor(allLoadings),
  method = "color",   # fill cells with colors
  type = "full",      # show full matrix
  addCoef.col = "black",   # show correlation coefficients
  number.cex = 0.7,        # text size for numbers
  tl.cex = 0.8             # text size for labels
)

########################################################################

plot(P[,3],efa_loadings[,5])

efa_loadings_diff <- efa_loadings[,4]-efa_loadings[,5]-efa_loadings[,3]
plot(P[,5],efa_loadings_diff)
cor(P[,5],efa_loadings_diff)

######################################################################## 

# CFA

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
fit = cfa(model=model, data=dati, ordered=T)
summary(fit, standardized=T)
fitMeasures(fit, fit.measures=c("rmsea","srmr","cfi","nnfi"))

modificationIndices(fit, sort.=T)[1:10,]

####

model_1 = "
SMD =~ bessi_1 + bessi_6 + bessi_11 + bessi_16 + bessi_21 + bessi_26 + bessi_31 + bessi_36 + bessi_41 
IND =~ bessi_5 + bessi_10 + bessi_15 + bessi_20 + bessi_25 + bessi_30 + bessi_35 + bessi_40 + bessi_45
COD =~ bessi_3 + bessi_8 + bessi_13 + bessi_18 + bessi_23 + bessi_28 + bessi_33 + bessi_38 + bessi_43 
SED =~ bessi_2 + bessi_7 + bessi_12 + bessi_17 + bessi_22 + bessi_27 + bessi_32 + bessi_37 + bessi_42 
ESD =~ bessi_4 + bessi_9 + bessi_14 + bessi_19 + bessi_24 + bessi_29 + bessi_34 + bessi_39 + bessi_44

bessi_15 ~~ bessi_40 
"
fit_1 = cfa(model=model_1, data=dati, ordered=T)
summary(fit_1, standardized=T)
fitMeasures(fit_1, fit.measures=c("rmsea","srmr","cfi","nnfi"))

modificationIndices(fit_1, sort.=T)[1:10,]

####
lambda_unstd = inspect(fit, what = "coef")$lambda
lambda_unstd

lambda = inspect(fit, what = "std")$lambda
lambda

#### confronto fattori trovati tramite EFA e CFA
quartz()
allLoadings = cbind(efa_loadings[item_names,], lambda[item_names,])
corrplot(
  cor(allLoadings),
  method = "color",   # fill cells with colors
  type = "full",      # show full matrix
  addCoef.col = "black",   # show correlation coefficients
  number.cex = 0.7,        # text size for numbers
  tl.cex = 0.8             # text size for labels
)

######################################################################## 
# CONFRONTO CFA vs NMF
# correspondences factor CFA loadings

# Riordina le righe di lambda in base al numero negli item bessi_1 ... bessi_45
lambda <- lambda[order(as.numeric(gsub("bessi_", "", rownames(lambda)))), ]
#lambda_unstd <- lambda_unstd[order(as.numeric(gsub("bessi_", "", rownames(lambda_unstd)))), ]

allLoadings = cbind(P,lambda)

quartz()
corrplot(
  cor(allLoadings),
  method = "color",   # fill cells with colors
  type = "full",      # show full matrix
  addCoef.col = "black",   # show correlation coefficients
  number.cex = 0.7,        # text size for numbers
  tl.cex = 0.8             # text size for labels
)

######################################################################## 
# CONFRONTO su DUE ITEM con EFA vs NMF

# EFA
item_data_2 <- dati[, c("bessi_13", "bessi_25", "bessi_30")]
#item_data_2 <- dati[, c("bessi_15", "bessi_25", "bessi_30", "bessi_40")]

item_data_2 <- data.frame(lapply(item_data_2, as.numeric))
item_data_complete_2 <- na.omit(item_data_2)
str(item_data_complete_2)

risultato_fa <- factanal(item_data_complete_2, factors = 1, rotation = "promax", scores = "regression")
print(risultato_fa, digits = 2, cutoff = 0.3)
risultato_fa$loadings

# NMF
# mi assicuro che la cartella "output" esista
if (!dir.exists("output")) dir.create("output")

# salvo il file di output
file_result_2 <- "output/result_bayesNMF_2item.rds"

# se il file esiste, lo carico
if (file.exists(file_result_2)) {
  message("✓ Carico il risultato salvato da file...")
  result2 <- readRDS(file_result_2)
} else {
  message("Eseguo bayesNMF()")
  
  # calcolo solo se serve
  M_t_2 <- t(as.matrix(item_data_complete_2))
  dim(M_t_2)
  
  result2 <- bayesNMF(
    data = M_t_2,
    likelihood = "normal",
    prior = "truncnormal",
    rank = 1
  )
  
  # salvo il risultato
  saveRDS(result2, file_result_2)
  message("Result salvato in: ", file_result_2)
}

# estraggo la stima MAP
MAP_2 <- result2$get_MAP()

######################################################################## 
# estrai le matrici dei fattori (P ed E)
P2 <- MAP_2$P
E2 <- MAP_2$E

# Se sono liste con intervalli, prendi la media dei limiti
if (is.list(P2)) {
  if (all(c("lower", "upper") %in% names(P2))) {
    P2 <- (P2$lower + P2$upper) / 2
  } else {
    P2 <- as.matrix(P2[[1]])
  }
}

if (is.list(E2)) {
  if (all(c("lower", "upper") %in% names(E2))) {
    E2 <- (E2$lower + E2$upper) / 2
  } else {
    E2 <- as.matrix(E2[[1]])
  }
}

P2 <- as.matrix(P2)
E2 <- as.matrix(E2)

########################################################################
# CONFRONTO TRA EFA E BAYESIAN NMF (1 fattore, 3 item)
########################################################################

# loadings EFA
efa_loadings <- as.matrix(risultato_fa$loadings[, 1, drop = FALSE])  # matrice item × 1

# loadings NMF (P2)
nmf_loadings <- round(P2, 2)  # matrice item × 1
rownames(nmf_loadings) <- paste0("bessi_", 1:nrow(nmf_loadings))
colnames(nmf_loadings) <- "nmf_loading_1"

########################################################################
# Combina i loadings EFA e NMF in un’unica matrice per il confronto
########################################################################

# rinomina la colonna EFA per chiarezza
colnames(efa_loadings) <- "efa_loading_1"

# unisci le due matrici per item
allLoadings <- cbind(nmf_loadings, efa_loadings)

# stampa di controllo
print(allLoadings)

########################################################################
# MATRICE DI CORRELAZIONE E PLOT
########################################################################

quartz()

corrplot(
  cor(allLoadings, use = "pairwise.complete.obs"),
  method = "color",       # riempie le celle con colori
  type = "full",          # mostra tutta la matrice
  addCoef.col = "black",  # mostra i coefficienti
  number.cex = 0.9,       # grandezza dei numeri
  tl.cex = 1,             # grandezza delle etichette
  mar = c(0, 0, 1, 0)
)
