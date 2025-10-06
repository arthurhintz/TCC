library(MxARMA)
library(doParallel)

source("2D_MxARMA/fit_2d_mxarma.R")
source("2D_MxARMA/simu_2d_mxarma.R")


# Nessa simulação fixei 9 posições juntas, porém irei validar apenas quando 
# mais de 2 posições detectar a mesma posição

# Sets
set.seed(1248)

n = 10
nrep <- 100

matriz_pos <- matrix(0, 4, 9)
resu <- matrix(NA, 4, 9)
final <- numeric(nrep)
erro2 <- numeric(nrep)

matriz_pos[1,] <- c(14, 15, 16, 24, 25, 26, 34, 35, 36)

nw <- ncol(matriz_pos)

# Valores para os parametros phi e theta
x <- c(seq(-0.3, -0.05, 0.02), seq(0.05, 0.3, 0.02)) # tirei o 0

rotate90 <- function(mat) {
  t(apply(mat, 2, rev))
}

# loop para preencher a matriz das posições dos alvos
for (i in 2:4) {
  m_bin <- matrix(0, nrow = n, ncol = n)
  
  m_bin[matriz_pos[i-1,]] <- 1
  
  m_bin <- rotate90(m_bin)
  
  matriz_pos[i,] <- which(m_bin == 1)
}
matriz_pos  

# Simualação com paralelismo

n_cores <- max(1L, parallel::detectCores() - 2L)
doParallel::registerDoParallel(cores = n_cores)


for (i in 1:nrep){
  print(i)

  w <- runif(nw, 1, 40)
  
  alpha <- runif(1, max = 0.6)
  phi <- matrix(c(sample(x,3), 0), ncol = 2)
  theta <- matrix(c(sample(x,3), 0), ncol = 2)
  
  sim <- tryCatch(
    mxarma2d.sim(n, n, alpha, phi, theta),
    error = function(e) return(NULL)
  )
  
  # Verifica se houve erro na simulação
  if (is.null(sim)) {
    cat("Erro na simulação ARMA na repetição", i, "\n")
    erro2[i] <- NA
    next  # pula para a próxima repetição
  }
  
  y <- sim$y
  y[matriz_pos[1,]] <- w
  
  resu <- matrix(NA, 4, 9)
  
  for (j in 1:4) {
    if (j > 1) y <- rotate90(y)
    
    fit <- tryCatch(
      mxarma2d.fit(y, 1, 1),
      error = function(e) return(NULL)
    )
    
    # Verifica se houve erro na simulação
    if (is.null(fit) || fit$conv != 0) {
      cat("Erro no fit na repetição", i, "\n")
      erro2[i] <- NA
      next  # pula para a próxima repetição
    }
    
    resi <- fit$resid
    
    resi_bin <- ifelse(abs(resi) > 3, 1, 0) # defini como 2 desvios, FAbio deixa 3
    resi_bin <- rbind(rep(0, n-1), resi_bin)   # vira n x (n-1)
    resi_bin <- cbind(rep(0, n), resi_bin)
    
    resu[j,] <-  resi_bin[matriz_pos[j,]]
  
  }
  
  
  # Para cada posição (coluna de resu), checar quantas rotações detectaram
  detec_por_pos <- colSums(resu, na.rm = TRUE)
  
  # Posição é considerada detectada se apareceu em >= 3 rotações
  pos_detectadas <- detec_por_pos >= 3
  
  # Erro tipo II = não detectar
  erro2[i] <- 1 - sum(as.numeric(pos_detectadas)) / nw
  
  
}

write.table(erro2, "erro2_simu3.txt", row.names = FALSE, col.names = FALSE)

erro2

mean(erro2)

matplot(erro2, type = "l", lty = 1, col = 1:ncol(erro2),
        xlab = "Índice", ylab = "Erro2")
legend("topright", legend = paste("Col", 1:ncol(erro2)),
       col = 1:ncol(erro2), lty = 1, cex = 0.8)


e2 <- read.table("erro2_simu2.txt")

mean(rowMeans(e2, na.rm = T), na.rm = T)
