library(geoR)
library(kableExtra)
library(dplyr)
library(lattice)


geodados <- read.geodata("AMF.txt", h=T)

dados <- geodados$data
coord <- geodados$coords


# análise descritiva da variável; 
summary(dados)

plot(geodados)

hmax2 <- summary(geodados)[[3]][[2]]*0.8

vario2 <- variog(geodados, trend="1st", estimator="classical", max.dist=hmax2)

vario.env <- variog.mc.env(geodados, obj.var=vario2)

plot(vario2)
plot(vario2, envelope=vario.env)






geodados <- read.geodata("AMF.txt", h=T)

dados <- geodados$data
coord <- geodados$coords



x <- coord[,1]
y <- coord[,2]

matriz <- xtabs(dados ~ x + y)
matriz

matriz[matriz == 0] <- 0.00001

fit1 <- mxarma2d.fit(matriz, 1, 1)


imagem <- fit1$fitted

image(imagem)
image(matriz)







