# Simulação
#

source("simu_2d_mxarma.R")
source("fit_2d_mxarma.R")



# Ver os valores do artigo rarma2d
phi_values <- c(0.4, 0.3, 0.1, 0)
phi <- matrix(phi_values, ncol = 2, nrow = 2, byrow = T)
phi

theta_values <- c(0.3, 0.15, 0.05, 0)
theta <- matrix(theta_values, ncol = 2, nrow = 2, byrow = T)
theta

alpha = 0.4

p <- 1
q <- 1

p1 = (p + 1)^2 - 1
q1 = (q + 1)^2 - 1

nrep = 2

nvalores = c(30, 50, 100, 200)
kvalores = c(30, 50, 100, 200)

vp <- c(alpha, phi_values[-4], theta_values[-4])
mf <- matrix(NA, ncol = length(vp), nrow = nrep)


results <- array(NA, dim = c(length(vp), 4, length(nvalores)))

colnames(results) <- c("parâmetros", "par_estimados", "viés_relativo", "EQM")
rownames(results) <- c("alpha", paste0("phi", 1:p1), paste0("theta", 1:q1))

for (j in seq_along(nvalores)) {

  k <- kvalores[j]
  n <- nvalores[j]

  for (i in 1:nrep) {

    rasu <- mxarma2d.sim(n = n, k = k, alpha = alpha,
                         phi = phi, theta = theta)

    y <- rasu$y

    fit <- mxarma2d.fit(y = y, p = p, q = q)
    mf[i,] <- fit$coeff

    cat("n = ", j, "  ",  "nrep =", i)
  }

  par_est <- apply(mf, 2, mean)
  vviesm <- (par_est - vp)
  viest <- (par_est - vp) / vp
  vvarm <- apply(mf, 2, var)
  veqmm <- vviesm^2 + vvarm

  results[ , ,j] <- cbind(vp, round(par_est, digits = 4),
                                        round(viest, digits = 5), round(veqmm, digits = 3))

  }

results

write.table(results, "Simu_MXARMA_2D_5000nrep.txt")
