library(MxARMA)

# Sets
set.seed(1212)

n = 10
posicoes <- c(14, 18, 33, 42, 55, 59, 66, 72, 74, 86)
nw <- length(posicoes)

x <- c(seq(-0.3, -0.05, 0.05), seq(0.05, 0.6, 0.05)) # tirei o 0

nrep <- 30

erro1 <- matrix(NA, nrow = nrep, ncol = 4)
erro2 <- matrix(NA, nrow = nrep, ncol = 4)

rotate90 <- function(mat) {
  t(apply(mat, 2, rev))
}

# Simualação

for (i in 1:nrep){

w <- runif(nw, 20, 30)


alpha <- runif(1, max = 0.6)
phi <- matrix(c(sample(x,3), 0), ncol = 2)
theta <- matrix(c(sample(x,3), 0), ncol = 2)

y <- mxarma2d.sim(n, n, alpha, phi, theta)$y



# irei fixar os valores das posições
#n_total <- length(y)
#sample(n_total, nw)

y[posicoes] <- w

#image(y)

m_bin <- matrix(0, nrow = nrow(y), ncol = ncol(y))
m_bin[posicoes] <- 1

# Visualizar
# image(
#   t(apply(m_bin, 2, rev)),
#   col = c("gray", "black"),
#   axes = FALSE
# )

  for(j in 1:4){
    
  y <- rotate90(y)
  
  m_bin <- rotate90(m_bin)
  posicoes <- which(m_bin > 0)

  fit <- mxarma2d.fit(y, 1, 1)

  resi <- fit$resid

  resi_bin <- ifelse(abs(resi) > 2, 1, 0) # defini como 2 desvios, FAbio deixa 3
  #image(resi_bin)



# Resultados
# comparando as detecões

  resi_bin <- rbind(rep(0, n-1), resi_bin)   # vira n x (n-1)
  resi_bin <- cbind(rep(0, n), resi_bin)

  

# Acerto (detectar quando deve detectar)
#sum(resi_bin[posicoes]) / 10

# Erro tipo II(Falso negativo) (Não detectar)

  erro1[i, j] <- 1 - sum(resi_bin[posicoes]) / nw

# Erro tipo I (Falso positivo) (detectou quando não existe)

  erro2[i, j] <- sum(resi_bin[-posicoes]) / (n*n-nw)

  }
}

erro1
erro2

mean(erro1)
mean(erro2)

matplot(erro2, type = "l", lty = 1, col = 1:ncol(erro2),
        xlab = "Índice", ylab = "Erro2")
legend("topright", legend = paste("Col", 1:ncol(erro2)),
       col = 1:ncol(erro2), lty = 1, cex = 0.8)


e1 <- read.table("erro1.txt")

colSums(is.na(e1)) # 10% não convergiu

mean(rowMeans(e1, na.rm = T), na.rm = T)


e2 <- read.table("erro2.txt")

mean(rowMeans(e2, na.rm = T), na.rm = T)
