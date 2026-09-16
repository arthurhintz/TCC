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

n = k = 30
nrep = 1000
m = 1



# Target position
#matriz_pos <- matrix(0, 4, 16)
#matriz_pos[1,] <- c(14,15,16,17, 24,25,26,27, 34,35,36,37, 44,45,46,47)

# 4 x 4 target
t_dim <- 4

# center of target
start_row <- floor((n - t_dim) / 2) + 1
start_col <- floor((k - t_dim) / 2) + 1

target_mask <- matrix(0, n, k)

target_mask[
  start_row:(start_row + t_dim - 1),
  start_col:(start_col + t_dim - 1)
] <- 1

nw <- sum(target_mask)

rotate90 <- function(mat) t(apply(mat, 2, rev))

# # gerar as 4 rotações
# for (r in 2:4) {
#   temp <- matrix(0, n, n)
#   temp[ matriz_pos[r-1,] ] <- 1
#   temp <- rotate90(temp)
#   matriz_pos[r,] <- which(temp == 1)
# }
# nw <- ncol(matriz_pos)  # number of target 16


# Error vectors
erro2 <- numeric(nrep)  # falso negativo por réplica
erro1 <- numeric(nrep)  # falso positivo por réplica

#==========/==========/==========/==========/==========/==========/==========/==========/
#-- Loop start 

for (i in 1:nrep) {
  
  cat("rep:", i ,"\n")
  
  sim <- mxarma2d.sim(n, n, alpha, phi, theta)
  y_orig <- sim$y  # sem alvo injetado
  
  #detections_by_rotation <- matrix(0, nrow = 4, ncol = nw)
  #fp_count_by_rotation <- numeric(4)
  
  detection_union <- matrix(0, n, k)
  
#==========/==========/==========/==========/==========/==========/==========/==========/  
  for (rot in 1:4) {
    
    y_rot <- y_orig
    mask_rot <- target_mask
    
    if (rot > 1) {
      for (r in 2:rot) {
        y_rot    <- rotate90(y_rot)
        mask_rot <- rotate90(mask_rot)
      }
    }
#==========/==========/==========/==========/==========/==========/==========/==========/    
    # Fit model without target
    fit <- mxarma2d.fit(y_rot, 1, 1)
    alpha_hat <- fit$alpha
    phi_hat   <- fit$phi
    theta_hat <- fit$theta
    
#==========/==========/==========/==========/==========/==========/==========/==========/
    # Inject Targets 
    y_inj <- y_rot
    
    y_inj[mask_rot == 1] <- 20 * mean(y_rot) # é maios ou menos isso a diferença da aplicacao
    
    ylog <- log(y_inj)
    
#==========/==========/==========/==========/==========/==========/==========/==========/
# Calculate fitted values and errors
    
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
    
#==========/==========/==========/==========/==========/==========/==========/==========/    
    # Quantile residuals
    
    fit_f <- exp(etahat[(m+1):n, (m+1):k])
    # observe que y[-1,-1] e fit_f têm dimensão (n-1)x(n-1)
    
    resi_mat <- qnorm(MxARMA::pmax(y_inj[-1, -1], fit_f))
    # garantir forma (n-1) x (n-1)
    
    resi_mat <- matrix(
      resi_mat,
      nrow = n - m,
      ncol = k - m
    )
    
#==========/==========/==========/==========/==========/==========/==========/==========/    
    # Detection
    
    resi_bin_small <- ifelse(abs(resi_mat) > 3, 1,0)
    
    # Put residual matrix back into n x k grid
    resi_bin <- matrix(0, n, k)
    
    resi_bin[(m + 1):n,(m + 1):k] <- resi_bin_small
    
#==========/==========/==========/==========/==========/==========/==========/==========/
    # Return detection to original orientation
        
    detection_original <- resi_bin
    
    if (rot > 1) {
      # inverse rotation:
      # rot = 2 -> rotate 3 times
      # rot = 3 -> rotate 2 times
      # rot = 4 -> rotate 1 time
      
      for (r in 1:(5 - rot)) {
        detection_original <- rotate90(detection_original)
      }
    }
    
    # Union:
    # detected if detected in at least one rotation
    detection_union <- base::pmax(
      detection_union,
      detection_original
    )
  }
    
#==========/==========/==========/==========/==========/==========/==========/==========/
  # Type I and Type II errors
  
  # True Positive
  TP <- sum(detection_union == 1 & target_mask == 1)
  
  # False negatives
  FN <- sum(detection_union == 0 & target_mask == 1)
  
  # False positives
  FP <- sum(detection_union == 1 & target_mask == 0)
  
  # True negatives
  TN <- sum(detection_union == 0 & target_mask == 0)
  
  # Type II error = false negative rate
  erro2[i] <- FN / (TP + FN)
  
  # Type I error = false positive rate
  erro1[i] <- FP / (FP + TN)
  
  cat("TP =", TP, "FN =", FN, "FP =", FP, "TN =", TN, "\n")
  
  cat(
    "Type II =", round(erro2[i], 4),
    "Type I =", round(erro1[i], 4),
    "\n"
  )
}
  
#     # 6) guarda detections nas posições verdadeiras DESSA rotação
#     detections_by_rotation[rot, ] <- resi_bin[ matriz_pos[rot,]]
#     
#     # 7) falso positivos: contar detections fora das nw posições
#     fp_count_by_rotation[rot] <- sum(resi_bin) - sum(resi_bin[matriz_pos[rot,]])
#   }
#   
#   # --- agregação entre rotações: posição detectada se detectada em pelo menos 1 rotação ---
#   detected_any_rotation <- apply(detections_by_rotation, 2, function(col) any(col == 1))
#   # erro tipo II (falso negativo) = 1 - (nº detectadas / nw)
#   erro2[i] <- 1 - (sum(detected_any_rotation) / nw)
#   
#   # erro tipo I (falso positivo)
#   detected_grid_union <- matrix(0, n, n)
#   for (rot in 1:4) {
#     NULL
#   }
#   # usar a média de FP por rotação (simplificação): FP rate = mean(fp_count_by_rotation) / (n*n - nw)
#   erro1[i] <- mean(fp_count_by_rotation) / (n*n - nw)
#   
#   #if (i %% 10 == 0) cat("Rep:", i, " FN=", round(erro2[i],3), " FP=", round(erro1[rep_idx],4), "\n")
# }

#==========/==========/==========/==========/==========/==========/==========/==========/
# Results

cat("\nType II error (FN rate):\n")
print(erro2)

cat("\nType I error (FP rate):\n")
print(erro1)

cat(
  "\nMean Type II =", mean(erro2),
  "\nMean Type I  =", mean(erro1),
  "\n"
)

