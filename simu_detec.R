library(MxARMA)

source("fit_2d_mxarma.R")
source("simu_2d_mxarma.R")

set.seed(1248)

# set pars
phi_values   <- c(0.03, 0.35, 0.25, 0)
theta_values <- c(-0.1, -0.06, -0.008, 0)
phi   <- matrix(phi_values,   ncol = 2, byrow = TRUE)
theta <- matrix(theta_values, ncol = 2, byrow = TRUE)
alpha <- -1.2

n = k = 10
nrep = 100
m = 1

# Posições e rotações
matriz_pos <- matrix(0, 4, 16)
matriz_pos[1,] <- c(14,15,16,17, 24,25,26,27, 34,35,36,37, 44,45,46,47)

rotate90 <- function(mat) t(apply(mat, 2, rev))

# gerar as 4 rotações
for (r in 2:4) {
  temp <- matrix(0, n, n)
  temp[ matriz_pos[r-1,] ] <- 1
  temp <- rotate90(temp)
  matriz_pos[r,] <- which(temp == 1)
}
nw <- ncol(matriz_pos)  # número de alvo (9)


erro2 <- numeric(nrep)  # falso negativo por réplica
erro1 <- numeric(nrep)  # falso positivo por réplica


for (i in 1:nrep) {
  
  cat("rep:", i ,"\n")
  
  sim <- mxarma2d.sim(n, n, alpha, phi, theta)
  y_orig <- sim$y  # sem alvo injetado
  
  detections_by_rotation <- matrix(0, nrow = 4, ncol = nw)
  
  fp_count_by_rotation <- numeric(4)
  
  for (rot in 1:4) {
    
    y_rot <- y_orig
    if (rot > 1) {
      for (t in 2:rot) {
        y_rot <- rotate90(y_rot)
      }
    }
    
    fit <- mxarma2d.fit(y_rot, 1, 1)
    alpha_hat <- fit$alpha
    phi_hat   <- fit$phi
    theta_hat <- fit$theta
    
    # 2) colocando os alvos
    y_inj <- y_rot
    #y_inj[ matriz_pos[rot,] ] <- 3 * max(y_rot)
    y_inj[ matriz_pos[rot,] ] <- 20 * mean(y_rot) # é maios ou menos isso a diferença da aplicacao
    
    ylog <- log(y_inj)
    
    # 3) calcula etahat e errorhat usando os coeficientes estimados
    etahat   <- matrix(0, n, n)
    errorhat <- matrix(0, n, n)
    
    for (ii in (m+1):n) {
      for (jj in (m+1):k) {
        y_block  <- as.vector(t(ylog[(ii-1):ii, (jj-1):jj]))
        y_block  <- y_block[-(3+1)]     # remover o elemento central conforme seu código
        err_block <- as.vector(t(errorhat[(ii-1):ii, (jj-1):jj]))
        err_block <- err_block[-(3+1)]
        
        etahat[ii, jj] <- alpha_hat + sum(phi_hat * y_block) + sum(theta_hat * err_block)
        errorhat[ii, jj] <- ylog[ii, jj] - etahat[ii, jj]
      }
    }
    
    fit_f <- exp(etahat[(m+1):n, (m+1):k])
    # observe que y[-1,-1] e fit_f têm dimensão (n-1)x(n-1)
    
    resi_mat <- qnorm(MxARMA::pmax(y_rot[-1, -1], fit_f))
    # garantir forma (n-1) x (n-1)
    
    resi_mat <- matrix(resi_mat, nrow = n-1)
    
    resi_bin <- ifelse(abs(resi_mat) > 3, 1, 0)
    resi_bin <- rbind(rep(0, n-1), resi_bin)   # adiciona primeira linha de zeros
    resi_bin <- cbind(rep(0, n), resi_bin)     # adiciona primeira coluna de zeros
    # agora resi_bin é n x n e pode indexar por posições lineares
    
    # 6) guarda detections nas posições verdadeiras DESSA rotação
    detections_by_rotation[rot, ] <- resi_bin[ matriz_pos[rot,]]
    
    # 7) falso positivos: contar detections fora das nw posições
    fp_count_by_rotation[rot] <- sum(resi_bin) - sum(resi_bin[matriz_pos[rot,]])
  }
  
  # --- agregação entre rotações: posição detectada se detectada em pelo menos 1 rotação ---
  detected_any_rotation <- apply(detections_by_rotation, 2, function(col) any(col == 1))
  # erro tipo II (falso negativo) = 1 - (nº detectadas / nw)
  erro2[i] <- 1 - (sum(detected_any_rotation) / nw)
  
  # erro tipo I (falso positivo)
  detected_grid_union <- matrix(0, n, n)
  for (rot in 1:4) {
    NULL
  }
  # usar a média de FP por rotação (simplificação): FP rate = mean(fp_count_by_rotation) / (n*n - nw)
  erro1[i] <- mean(fp_count_by_rotation) / (n*n - nw)
  
  #if (i %% 10 == 0) cat("Rep:", i, " FN=", round(erro2[i],3), " FP=", round(erro1[rep_idx],4), "\n")
}

# resultados
cat("Erro tipo II (FN) por réplica:\n")
print(erro2)
cat("Erro tipo I (FP) por réplica (estimativa média por rotação):\n")
print(erro1)

# médias
cat("Média Erro tipo 2 =", mean(erro2), "   Média Erro tipo 1 =", mean(erro1), "\n")

