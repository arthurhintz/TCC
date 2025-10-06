library(R.matlab) # read mat files
library(raster)
#library(devtools) 
#install_github("murilosagrillo/ladistribution")
library(ladistribution) # LA distribution
library(mmand) # morphologic operations


# read in our data
data1 <- readMat("2s1_real_A_elevDeg_015_azCenter_011_22_serial_b01.mat")

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


