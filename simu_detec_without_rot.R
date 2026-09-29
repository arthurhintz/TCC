library(MxARMA)

source("fit_2d_mxarma.R")
source("simu_2d_mxarma.R")

set.seed(1248)

#===========================================================
# Parameters
#===========================================================

phi_values   <- c(0.03, 0.35, 0.25, 0)
theta_values <- c(-0.1, -0.06, -0.008, 0)

phi <- matrix(
  phi_values,
  ncol = 2,
  byrow = TRUE
)

theta <- matrix(
  theta_values,
  ncol = 2,
  byrow = TRUE
)

alpha <- -1.2

n = k = 50

target_dim <- 10

nrep <- 500
m <- 1

# Target position

# Target size: 4 x 4
t_dim <- 4

# Centralized target
start_row <- floor((n - t_dim) / 2) + 1
start_col <- floor((k - t_dim) / 2) + 1

target_mask <- matrix(0, n, k)

target_mask[
  start_row:(start_row + t_dim - 1),
  start_col:(start_col + t_dim - 1)
] <- 1

# Number of target pixels
nw <- sum(target_mask)

cat("Number of target pixels:", nw, "\n")
cat(
  "Rows:",
  start_row:(start_row + t_dim - 1),
  "\n"
)
cat(
  "Columns:",
  start_col:(start_col + t_dim - 1),
  "\n"
)

#===========================================================
# Storage

erro2 <- numeric(nrep)   # False negative rate
erro1 <- numeric(nrep)   # False positive rate

TP_vec <- numeric(nrep)
FN_vec <- numeric(nrep)
FP_vec <- numeric(nrep)
TN_vec <- numeric(nrep)

#===========================================================
# Simulation loop

for (i in 1:nrep) {
  
  cat("rep:", i, "\n")
  
  #---------------------------------------------------------
  # 1. Simulate image WITHOUT target

  sim <- mxarma2d.sim(n,k,alpha,phi,theta )
  
  y_orig <- sim$y
  
  #---------------------------------------------------------
  # 2. Fit model WITHOUT target

  fit <- mxarma2d.fit(y_orig,1,1)
  
  alpha_hat <- fit$alpha
  phi_hat   <- fit$phi
  theta_hat <- fit$theta
  
  #---------------------------------------------------------
  # 3. Inject target
  
  y_inj <- y_orig
  
  y_inj[target_mask == 1] <- target_dim * mean(y_orig)
  
  ylog <- log(y_inj)
  
  #---------------------------------------------------------
  # 4. Calculate fitted values and errors
  
  etahat <- matrix(0, n,k)
  
  errorhat <- matrix(0, n,k)
  
  for (ii in (m + 1):n) {
    
    for (jj in (m + 1):k) {
      
      # Previous observations
      y_block <- as.vector(
        t(
          ylog[
            (ii - 1):ii,
            (jj - 1):jj
          ]
        )
      )
      
      # Remove current observation
      y_block <- y_block[-4]
      
      # Previous errors
      err_block <- as.vector(
        t(
          errorhat[
            (ii - 1):ii,
            (jj - 1):jj
          ]
        )
      )
      
      # Remove current error
      err_block <- err_block[-4]
      
      # Linear predictor
      etahat[ii, jj] <-
        alpha_hat +
        sum(phi_hat * y_block) +
        sum(theta_hat * err_block)
      
      # Error
      errorhat[ii, jj] <-
        ylog[ii, jj] -
        etahat[ii, jj]
    }
  }
  
  #---------------------------------------------------------
  # 5. Fitted values

  fit_f <- exp(etahat[(m + 1):n,(m + 1):k])
  
  #---------------------------------------------------------
  # 6. Quantile residuals

  resi_mat <- qnorm(
    MxARMA::pmax(
      y_inj[
        (m + 1):n,
        (m + 1):k
      ],
      fit_f
    )
  )
  
  resi_mat <- matrix(
    resi_mat,
    nrow = n - m,
    ncol = k - m
  )
  
  #---------------------------------------------------------
  # 7. Detection

  # Detection if |residual| > 2
  resi_bin_small <- ifelse(abs(resi_mat) > 2, 1,0)
  
  # Return to original n x k dimension
  resi_bin <- matrix(0,n, k)
  
  resi_bin[
    (m + 1):n,
    (m + 1):k
  ] <- resi_bin_small
  
  #---------------------------------------------------------
  # 8. Compare detection with TRUE target

  # True positive:
  # target exists AND was detected
  TP <- sum(
    target_mask == 1 &
      resi_bin == 1
  )
  
  # False negative:
  # target exists BUT was not detected
  FN <- sum(
    target_mask == 1 &
      resi_bin == 0
  )
  
  # False positive:
  # no target BUT detection occurred
  FP <- sum(
    target_mask == 0 &
      resi_bin == 1
  )
  
  # True negative:
  # no target AND no detection
  TN <- sum(
    target_mask == 0 &
      resi_bin == 0
  )
  
  #---------------------------------------------------------
  # 9. Store results

  TP_vec[i] <- TP
  FN_vec[i] <- FN
  FP_vec[i] <- FP
  TN_vec[i] <- TN
  
  # Type II error = False Negative Rate
  erro2[i] <- FN / (TP + FN)
  
  # Type I error = False Positive Rate
  erro1[i] <- FP / (FP + TN)
  
  cat(
    "TP =", TP,
    "FN =", FN,
    "FP =", FP,
    "TN =", TN,
    "\n"
  )
  
  cat(
    "Type II =", round(erro2[i], 4),
    "Type I =", round(erro1[i], 4),
    "\n"
  )
}

#===========================================================
# Results
#===========================================================


file_1 <- paste0("erro_1_nk", n, "target", target_dim, ".txt")
file_2 <- paste0("erro_2_nk", n, "target", target_dim, ".txt")


write.table(erro1, file_1, row.names = FALSE, col.names = FALSE)
write.table(erro2, file_2, row.names = FALSE, col.names = FALSE)




cat("\nType II error (FN rate):\n")
print(erro2)

cat("\nType I error (FP rate):\n")
print(erro1)

cat(
  "\nMean Type II =", mean(erro2),
  "\nMean Type I  =", mean(erro1),
  "\n"
)
