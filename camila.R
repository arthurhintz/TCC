library(R.matlab) # read mat files
library(raster)
library(mmand) # morphologic operations
#library(hexbin)

# comando para enviar para o github
# git push origin main

source("fit_2d_mxarma.R")


# read in our data
#data1 <- readMat("2s1_real_A_elevDeg_015_azCenter_011_22_serial_b01.mat")

# lendo dados do github

files <- paste0("2s1_synth_A_elevDeg_015_azCenter_0",16,"_22_serial_b01.mat")
url <- "https://raw.githubusercontent.com/benjaminlewis-afrl/SAMPLE_dataset_public/1bae4156afcba7ad9c9f2eb743ed20675ed14fb1/mat_files/synth/2s1/"


data1 <- readMat(paste0(url,files))

# check out data structure
str(data1)

im_c<- data1$complex.img # complex SAR
im_a<- raster(abs(data1$complex.img)) # amplitude of SAR data
im_matrix<-ifelse(as.matrix(im_a)==0,0.0001,as.matrix(im_a)) # take off zero values

dim(im_matrix)
summary(as.vector(im_matrix))

plot(im_a, col =grey(seq(0, 1, length = 512))) # plot with raster

image(1:dim(im_matrix)[1], 1:dim(im_matrix)[2], im_matrix, xlab="", ylab="",
      main="SAR image",col =grey(seq(0, 1, length = 512)))


y <- im_matrix[20:90, 40:110]


fit <- mxarma2d.fit(y, 1, 1)

fit$fitted

image(fit$fitted, xlab="", ylab="",
      main="SAR image Fitted",col =grey(seq(0, 1, length = 512)))

im <- y[-1,-1]
eqmM <-mean((im-fit$fitted)^2, na.rm = T)
eqmM

# residuos estão como vetor -> colocar como matriz (att func)
resi <- matrix(fit$resid, nrow(im))

matbin <- ifelse(abs(resi) > 3, 1, 0)

image(resi)

sum(matbin)
image(matbin, col = c("white", "black"), axes = FALSE)

k <-matrix(1, nrow = 3, ncol = 3)

E <- erode(matbin,k)
D <-dilate(E, k)

sum(D)

image(t(D)[,nrow(D):1], col = c("white", "black"), axes = FALSE)


#==========/==========/==========/==========/==========/==========/==========/==========/
# irei treinar meu modelo sem o alvo e depois ajustar com ele

treino <- im_matrix[1:50, 1:50]
teste <- im_matrix[41:90, 41:90]

image(treino)
image(teste)

fit_1 <- mxarma2d.fit(treino, 1, 1) 


fit_1$model

alpha <- fit_1$model[1,1]
phi <- matrix(c(fit_1$model[2:4,1],0), 2)
theta <- matrix(c(fit_1$model[5:7,1],0), 2)

n = k = 50
m =1
etahat = errorhat = matrix(0, ncol = k, nrow = n)

ylog <- log(teste)

for (i in (m+1):n)
{
  for (j in (m+1):k)
  {
    ynew1 <- as.vector(t(ylog[(i-1):i,(j-1):j]))
    ynew1 <- ynew1[-(3+1)]
    
    errorhat_new <- as.vector(t(errorhat[(i-1):i,(j-1):j]))
    errorhat_new <- errorhat_new[-(3+1)]
    
    etahat[i,j]  = alpha + sum(phi * ynew1) + sum(theta * errorhat_new)
    errorhat[i,j] = ylog[i,j] - etahat[i,j]
  }
}

fit_f = exp(etahat[(m+1):n,(m+1):k])

image(fit_f)

resido <- qnorm(MxARMA::pmax(teste[-1,-1], fit_f))

matbin <- ifelse(abs(resido) > 3, 1, 0)
matbin <- matrix(matbin, nrow = n-1)

image(t(matbin)[,nrow(matbin):1])


sum(matbin)
image(matbin, col = c("black", "white"), axes = FALSE)

rmse <- sqrt(mean((teste[-1,-1]-fit_f)^2, na.rm = T))
rmse

# Don't work
# k <-matrix(1, nrow = 3, ncol = 3)
# E <-erode(matbin,k)
# D <-dilate(E, k)
# sum(D)
# image(t(D)[,nrow(D):1], col = c("black", "white"), axes = FALSE)

# Renata Cod

ker = mmand::shapeKernel(c(3,3), type="box") # kernel
img2<- mmand::closing(matbin,ker)
img3<- mmand::opening(img2,ker)
plot(raster(img2), col =grey(seq(0, 1, length = 512)))
plot(raster(img3), col =grey(seq(0, 1, length = 512)))
