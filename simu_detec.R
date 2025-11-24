library(MxARMA)

source("fit_2d_mxarma.R")
source("simu_2d_mxarma.R")

set.seed(1248)

# Parâmetros verdadeiros
phi_values   <- c(0.03, 0.35, 0.25, 0)
theta_values <- c(-0.1, -0.06, -0.008, 0)
phi   <- matrix(phi_values,   ncol = 2, byrow = TRUE)
theta <- matrix(theta_values, ncol = 2, byrow = TRUE)
alpha <- -1.2

n = k = 10
nrep = 10
m = 1

#---- Posições e rotações ----#

matriz_pos <- matrix(0, 4, 9)
matriz_pos[1,] <- c(14,15,16, 24,25,26, 34,35,36)

rotate90 <- function(mat) t(apply(mat, 2, rev))

# gerar as 4 rotações
for (i in 2:4) {
  temp <- matrix(0, n, n)
  temp[matriz_pos[i-1,]] <- 1
  temp <- rotate90(temp)
  matriz_pos[i,] <- which(temp == 1)
}
nw <- ncol(matriz_pos)


coefs <- matrix(NA, nrow = 4, ncol = 7)

sim <- mxarma2d.sim(n, n, alpha, phi, theta)
y <- sim$y  
y0 <- y

for(j in 1:4){
  if (j > 1) y <- rotate90(y)
  
  fit0 <- mxarma2d.fit(y, 1, 1)
  
  alpha_hat <- fit0$alpha
  phi_hat   <- fit0$phi
  theta_hat <- fit0$theta
  
  coefs[j, ] <- c(alpha_hat, phi_hat, theta_hat)
}


coefs


erro1 <- erro2 <- numeric(nrep)

#---- Loop principal ----#

for (r in 1:nrep) {
  cat("Repetição:", r, "\n")
  
  resu <- matrix(0, 4, 9)
  
  for (rot in 1:4) {
    
    # Rotação da imagem
    y <- y0
    if (rot > 1) y <- rotate90(y)
    
    # Inserir o alvo na rotação correta
    y1 <- y
    y1[matriz_pos[rot,]] <- 2 * max(y)
    ylog <- log(y1)
    
    # Inicializar matrizes para erro e etahat
    etahat    <- matrix(0, n, n)
    errorhat  <- matrix(0, n, n)
    
    # Previsões usando parâmetros fixados (fit inicial!)
    for (i in (m+1):n) {
      for (j in (m+1):k) {
        
        y_block <- as.vector(t(ylog[(i-1):i,(j-1):j]))
        y_block <- y_block[-(3+1)]
        
        err_block <- as.vector(t(errorhat[(i-1):i,(j-1):j]))
        err_block <- err_block[-(3+1)]
        
        etahat[i,j]  <- alpha_hat + sum(phi_hat * y_block) +
          sum(theta_hat * err_block)
        
        errorhat[i,j] <- ylog[i,j] - etahat[i,j]
      }
    }
    
    fit_f <- exp(etahat[(m+1):n,(m+1):k])
    resi  <- qnorm( MxARMA::pmax(y[-1,-1], fit_f) )
    
    # binariza
    resi_bin <- ifelse(abs(resi) > 3, 1, 0)
    resi_bin <- rbind(rep(0, n-1), resi_bin)
    resi_bin <- cbind(rep(0, n), resi_bin)
    
    # guarda apenas nas posições da rotação
    resu[rot,] <- resi_bin[matriz_pos[rot,]]
  }
  
  # detecção quando pelo menos 1 rotação detectou
  detec_por_pos <- colSums(resu)
  pos_detectadas <- detec_por_pos >= 1
  
  # Erros
  posicoes <- matriz_pos[1,]   # sempre usa o alvo da rotação original
  
  erro2[r] <- 1 - sum(pos_detectadas[posicoes]) / nw
  erro1[r] <- sum(pos_detectadas[-posicoes]) / (n*n - nw)
}

erro1
erro2
