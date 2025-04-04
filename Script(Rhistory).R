library(SummarizedExperiment)
library(readxl)

GastricCancer_Data <- read_excel("C:/Users/Jaume/Desktop/GastricCancer_NMR.xlsx", sheet = "Data")
GastricCancer_Peak <- read_excel("C:/Users/Jaume/Desktop/GastricCancer_NMR.xlsx", sheet = "Peak")
GastricCancer_Data <- GastricCancer_Data[, -c(1:4)]
GastricCancer_Peak <- GastricCancer_Peak[, -1]
Data <- as.matrix(GastricCancer_Data)
colData <- as.data.frame(GastricCancer_Peak)
se <- SummarizedExperiment(assay=list(Data=Data),
                           colData=colData)
gorrionesNum <- scale(assays(se)$Data, center = TRUE, scale=FALSE)
apply(gorrionesNum,2, mean)
Data <- na.omit(Data)
se <- SummarizedExperiment(assay=list(Data=Data),
                           colData=colData)
gorrionesNum <- scale(assays(se)$Data, center = TRUE, scale=FALSE)
apply(gorrionesNum,2, mean)
class(assays(se)$Data)
dim(assays(se)$Data)
colnames(assays(se)$Data)
View(Data)
gorrionesNum <- scale(assays(se)$Data, use = "pairwise.complete.obs", center = TRUE, scale=FALSE)
Data <- as.matrix(GastricCancer_Data)
se <- SummarizedExperiment(assay=list(Data=Data),
                           colData=colData)
gorrionesNum <- scale(assays(se)$Data, use = "pairwise.complete.obs", center = TRUE, scale=FALSE)
gorrionesNum <- scale(assays(se)$Data, center = TRUE, scale=FALSE)
apply(gorrionesNum,2, mean)
n<- dim(gorriones)[1]
n<- dim(assays(se)$Data)[1]
S<-cov(gorrionesNum)*(n-1)/n
show(S)
R<-cor(gorrionesNum)
show(R)
#Substituit los NA por la media de cada columna
assays(se)$Data[] <- lapply(assays(se)$Data, function(x) {
  mean_val <- mean(x, na.rm = TRUE)
  x[is.na(x)] <- mean_val
  return(x)
})
#Substituit los NA por la media de cada columna
assays(se)$Data[] <- lapply(assays(se)$Data, function(x) {
  mean_val <- mean(x, na.rm = TRUE)  # Calcular la media, ignorando los NA
  x[is.na(x)] <- mean_val  # Reemplazar los NA con la media
  return(x)
})
#Substituit los NA por la media de cada columna
Data[] <- lapply(Data, function(x) {
  mean_val <- mean(x, na.rm = TRUE)
  x[is.na(x)] <- mean_val
  return(x)
})
colData <- as.data.frame(GastricCancer_Peak)
se <- SummarizedExperiment(assay=list(Data=Data),
                           colData=colData)
View(Data)
#Substituit los NA por la media de cada columna
GastricCancer_Data [] <- lapply(GastricCancer_Data , function(x) {
  mean_val <- mean(x, na.rm = TRUE)
  x[is.na(x)] <- mean_val
  return(x)
})
Data <- as.matrix(GastricCancer_Data)
colData <- as.data.frame(GastricCancer_Peak)
se <- SummarizedExperiment(assay=list(Data=Data),
                           colData=colData)
#Matriu de covariança de les dades centrades
DataNum <- scale(assays(se)$Data, center = TRUE, scale=FALSE)
apply(gorrionesNum,2, mean)
GastricCancer_Data <- read_excel("C:/Users/Jaume/Desktop/GastricCancer_NMR.xlsx", sheet = "Data")
GastricCancer_Data <- GastricCancer_Data[, -c(1:4)]
#Substituit los NA por la media de cada columna
GastricCancer_Data [] <- lapply(GastricCancer_Data , function(x) {
  mean_val <- mean(x, na.rm = TRUE)
  x[is.na(x)] <- mean_val
  return(x)
})
Data <- as.matrix(GastricCancer_Data)
colData <- as.data.frame(GastricCancer_Peak)
se <- SummarizedExperiment(assay=list(Data=Data),
                           colData=colData)
se
class(assays(se)$Data)
dim(assays(se)$Data)
#Histograma muestras
hist(assays(se)$Data[,1], main="Sample_1")
#Histograma metabòlito
hist(assays(se)$Data[1,], main="M1")
#Matriu de covariança de les dades centrades
DataNum <- scale(assays(se)$Data, center = TRUE, scale=FALSE)
apply(gorrionesNum,2, mean)
#Matriu de varances
n<- dim(assays(se)$Data)[1]
S<-cov(DataNum)*(n-1)/n
show(S)
#Matriu de correlacions
R<-cor(DataNum)
show(R)
apply(DataNum,2, mean)
#Diagonalización de la matriz de covarianzas
EIG <- eigen(S)
show(EIG)
#Multiplicam la matriu original per la matriu de vectors pròpis
eigenVecs1 <- EIG$vectors
PCAS1 <- DataNum %*% eigenVecs1
head(PCAS1)
#Representació PCA 1 i 2
plot(PCAS1[,1], PCAS1[,2], main = "Gorriones. 2 primeras PCs")
#Porcentatges de variabilitat
vars1<- EIG$values/sum(EIG$values)
round(vars1,3)
xlabel <- paste("PCA1 ", round(vars1[1]*100, 2),"%" )
ylabel <- paste("PCA2 ", round(vars1[2]*100,2),"%" )
plot(PCAS1[,1], PCAS1[,2], main = "Metabòlits 2 primeres PCs",
     xlab=xlabel, ylab=ylabel)
PCAS2 <- princomp(ays(se)$Data)
PCAS2 <- princomp(assays(se)$Data)
PCAS2 <- princomp(Data)
EIG$vectors
#Visualitzar les classes
se
class(assays(se)$Data)
class(colData(se))
dim(assays(se)$Data)
dim(colData(se))
colnames(assays(se)$Data)
colnames(colData(se))
head(colData(se))
assays(se)$Data
colData(se)
head(colData(se))
assays(se)$Data
colData(se)
head(assays(se)$Data)
head(colData(se))
summary(assays(se)$Data)
#Histograma muestres
hist(assays(se)$Data[,1], main="Sample_1")
#Histograma metabòlito
hist(assays(se)$Data[1,], main="M1")
#Matriu de covariança de les dades centrades
DataNum <- scale(assays(se)$Data, center = TRUE, scale=FALSE)
apply(DataNum,2, mean)
n<- dim(assays(se)$Data)[1]
S<-cov(DataNum)*(n-1)/n
show(S)
#Matriu de correlacions
R<-cor(DataNum)
show(R)
#Diagonalizació de la matriu de covarianças
EIG <- eigen(S)
show(EIG)
#Multiplicam la matriu original per la matriu de vectors pròpis
eigenVecs1 <- EIG$vectors
PCAS1 <- DataNum %*% eigenVecs1
head(PCAS1)
#Representació PCA 1 i 2
plot(PCAS1[,1], PCAS1[,2], main = "Metabòlits 2 primeres PCs")
#Percentatges de variabilitat
vars1<- EIG$values/sum(EIG$values)
round(vars1,3)
xlabel <- paste("PCA1 ", round(vars1[1]*100, 2),"%" )
ylabel <- paste("PCA2 ", round(vars1[2]*100,2),"%" )
plot(PCAS1[,1], PCAS1[,2], main = "Metabòlits 2 primeres PCs",
     xlab=xlabel, ylab=ylabel)
# Guardar l'objecte SummarizedExperiment a un fitxer .Rda
save(se, file = "GastricCancer_SummarizedExperiment.Rda")
