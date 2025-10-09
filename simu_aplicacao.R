library(R.matlab) # read mat files
library(raster)

source("fit_2d_mxarma.R")
source("https://raw.githubusercontent.com/arthurhintz/RARMA2D/refs/heads/main/R/rarma.R")
source("https://raw.githubusercontent.com/arthurhintz/RARMA2D/refs/heads/main/R/pr.R")

# sets
rmseM <- c()
rmseR <- c()

# lendo dados do github (sintéticos)
files <- paste0("2s1_synth_A_elevDeg_015_azCenter_0",c(10:36, 38:46) ,"_22_serial_b01.mat")
base_url <- "https://raw.githubusercontent.com/benjaminlewis-afrl/SAMPLE_dataset_public/1bae4156afcba7ad9c9f2eb743ed20675ed14fb1/mat_files/synth/2s1/"

data_list <- lapply(files, function(f) {
  url <- paste0(base_url, f)
  readMat(url)
})

names(data_list) <- files


for (i in seq_along(files)) {

  print(i)
  data1 <- data_list[[i]]
  
  # check out data structure
  im_c<- data1$complex.img # complex SAR
  im_a<- raster(abs(data1$complex.img)) # amplitude of SAR data
  im_matrix<-ifelse(as.matrix(im_a)==0,0.0001,as.matrix(im_a)) # take off zero values
  
  plot(im_a, col =grey(seq(0, 1, length = 512))) # plot with raster
  
  y <- im_matrix[30:90, 40:100]
  
  fitM <- mxarma2d.fit(y, 1, 1)
  ajustM <- fitM$fitted
  
  fitR <- rarma2d(y, ar = 1, ma = 1, ne = 1)
  ajustR <- fitR$fitted
  
  im <- y[-1,-1]
  
  #RMSE
  rmseM[i] <- sqrt(mean((im-ajustM)^2, na.rm = T))
  
  rmseR[i] <- sqrt(mean((im-ajustR)^2, na.rm = T))
}

dif <- rmseR - rmseM

max(dif)
which.min(rmseM)

# perde 5, 22, 27, 33, 36
# Imagem 23 mais diferença
# Menor RMSE imagem 7