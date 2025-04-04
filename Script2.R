library(SummarizedExperiment)
library(readxl)

#Carregar els archius
GastricCancer_Data <- read_excel("C:/Users/Jaume/Desktop/GastricCancer_NMR.xlsx", sheet = "Data")
GastricCancer_Peak <- read_excel("C:/Users/Jaume/Desktop/GastricCancer_NMR.xlsx", sheet = "Peak")

#Eliminar les columnes que no interessen
GastricCancer_Data <- GastricCancer_Data[, -c(1:4)]
GastricCancer_Peak <- GastricCancer_Peak[, -1]

#Substituit els NA per la mitjana de cada metabòlit
GastricCancer_Data [] <- lapply(GastricCancer_Data , function(x) {
  mean_val <- mean(x, na.rm = TRUE) 
  x[is.na(x)] <- mean_val 
  return(x)
})

#Convertir les dades a matriu i dataframe
Data <- as.matrix(GastricCancer_Data)
colData <- as.data.frame(GastricCancer_Peak)

#Crear SummarizedExperiment
se <- SummarizedExperiment(assay=list(Data=Data),
                           colData=colData)

#Visualitzar les classes
se

class(assays(se)$Data)
class(colData(se))

dim(assays(se)$Data)
dim(colData(se))

colnames(assays(se)$Data)
colnames(colData(se))

head(assays(se)$Data)
head(colData(se))

#Resum estadísitic
summary(assays(se)$Data)

#Histograma muestres
hist(assays(se)$Data[,1], main="Sample_1")

#Histograma metabòlito
hist(assays(se)$Data[1,], main="M1")

#Estimació de la matriu de covariança
DataNum <- scale(assays(se)$Data, center = TRUE, scale=FALSE)
apply(DataNum,2, mean)
n<- dim(assays(se)$Data)[1]
S<-cov(DataNum)*(n-1)/n
show(S)

#Matriu de correlacions
R<-cor(DataNum)
show(R)

#Càlcul dels valors propis
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

# Mostrar el directori de treball actual
getwd()








