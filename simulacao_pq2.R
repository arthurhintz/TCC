# Simulação


source("simu_2d_mxarma.R")
source("fit_2d_mxarma.R")

library(doParallel)
library(doRNG)


# Sets
phi_values <- c(rep(0.1, 8), 0)
phi <- matrix(phi_values, ncol = 3, nrow = 3, byrow = T)
phi

theta = phi

alpha = 0.1
nrep = 2
p = q = 2

nvalores = c(10, 20, 40, 80, 120)
kvalores = c(10, 20, 40, 80, 120)

vp <- c(alpha, phi_values[-4], theta_values[-4])
nm_par <- c("alpha", paste0("phi", 1:8), paste0("theta", 1:8))
stopifnot(length(vp) == length(nm_par))
npar <- length(nm_par)


# paralelismo

set.seed(2025)
n_cores <- max(1L, parallel::detectCores() - 1L)
cl <- parallel::makeCluster(n_cores)
on.exit(parallel::stopCluster(cl), add = TRUE)

registerDoParallel(cl)
doRNG::registerDoRNG(123) 


parallel::clusterExport(cl, c("mxarma2d.sim", "mxarma2d.fit"), envir = environment())

message(sprintf("Paralelizando com %d núcleos.", n_cores))


for (j in seq_along(nvalores)) {
  
  k <- kvalores[j]
  n <- nvalores[j]
  
  nrep_j <- nrep
  
  mf <- foreach(i = 1:nrep_j, .combine = rbind, .inorder = TRUE,
                .export = c("alpha","phi","theta","p","q"),
                .errorhandling = "pass") %dopar% {
                  
                  out <- tryCatch({
                    rasu <- mxarma2d.sim(n = n, k = k, alpha = alpha, phi = phi, theta = theta)
                    y <- rasu$y
                    
                    fit <- mxarma2d.fit(y = y, p = p, q = q)
                    
                    # Checa se convergiu (caso sua função tenha esse campo)
                    if (!is.null(fit$conv) && fit$conv != 0) stop("Não convergiu")
                    
                    as.numeric(fit$coeff)
                    
                  }, error = function(e) {
                    cat(sprintf("[Erro] n=%d k=%d rep=%d: %s\n", n, k, i, conditionMessage(e)))
                    rep(NA_real_, npar)
                  })
                  
                  out
                }
  
  mf <- as.data.frame(mf, optional = TRUE, stringsAsFactors = FALSE)
  colnames(mf) <- nm_par
  
  
  arq_out <- paste0("simu_2DMxARMA(22)_nk", n, ".txt")
  utils::write.table(mf, file = arq_out, row.names = FALSE, col.names = TRUE, quote = FALSE)
  
  message(sprintf("Concluído: j=%d, n=%d, k=%d, réplicas=%d -> %s",
                  j, n, k, nrep_j, arq_out))
}

