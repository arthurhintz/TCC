suppressPackageStartupMessages({
  library(doParallel)
  library(foreach)
  if (!requireNamespace("doRNG", quietly = TRUE)) install.packages("doRNG")
  library(doRNG)
})


script_simu <- normalizePath("simu_2d_mxarma.R", mustWork = TRUE)
script_fit  <- normalizePath("fit_2d_mxarma.R",  mustWork = TRUE)


source(script_simu)
source(script_fit)


phi_values   <- c(0.4, 0.3, 0.1, 0)
phi          <- matrix(phi_values, ncol = 2, nrow = 2, byrow = TRUE); phi

theta_values <- c(0.3, 0.15, 0.05, 0)
theta        <- matrix(theta_values, ncol = 2, nrow = 2, byrow = TRUE); theta

alpha <- 0.4
p <- 1
q <- 1

p1 <- (p + 1)^2 - 1
q1 <- (q + 1)^2 - 1


nrep <- 2

nvalores <- c(10, 30, 50, 80, 150, 200, 400)
kvalores <- c(10, 30, 50, 80, 150, 200, 400)


vp <- c(alpha, phi_values[-4], theta_values[-4])
nm_par <- c("alpha", paste0("phi", 1:p1), paste0("theta", 1:q1))
stopifnot(length(vp) == length(nm_par))
npar <- length(nm_par)


set.seed(2025)
n_cores <- max(1L, parallel::detectCores() - 1L)
cl <- parallel::makeCluster(n_cores)
on.exit(parallel::stopCluster(cl), add = TRUE)

registerDoParallel(cl)
doRNG::registerDoRNG(123) 


parallel::clusterExport(cl, c("script_simu", "script_fit"), envir = environment())
parallel::clusterEvalQ(cl, {
  source(script_simu)
  source(script_fit)
  NULL
})

message(sprintf("Paralelizando com %d núcleos.", n_cores))


for (j in seq_along(nvalores)) {
  
  k <- kvalores[j]
  n <- nvalores[j]
  
  nrep_j <- nrep
  
  mf <- foreach(i = 1:nrep_j, .combine = rbind, .inorder = TRUE,
                .export = c("alpha","phi","theta","p","q"),
                .errorhandling = "pass") %dopar% {
                  out <- try({
                    rasu <- mxarma2d.sim(n = n, k = k, alpha = alpha, phi = phi, theta = theta)
                    y <- rasu$y
                    fit <- mxarma2d.fit(y = y, p = p, q = q)
                    as.numeric(fit$coeff)
                  }, silent = TRUE)
                  
                  if (inherits(out, "try-error")) {
                    rep(NA_real_, npar)
                  } else {
                    if (length(out) != npar) {
                      tmp <- rep(NA_real_, npar)
                      tmp[seq_len(min(npar, length(out)))] <- out[seq_len(min(npar, length(out)))]
                      tmp
                    } else {
                      out
                    }
                  }
                }
  
  mf <- as.data.frame(mf, optional = TRUE, stringsAsFactors = FALSE)
  colnames(mf) <- nm_par
  

  arq_out <- paste0("simu_2DMxARMA_nk", n, ".txt")
  utils::write.table(mf, file = arq_out, row.names = FALSE, col.names = TRUE, quote = FALSE)
  
  message(sprintf("Concluído: j=%d, n=%d, k=%d, réplicas=%d -> %s",
                  j, n, k, nrep_j, arq_out))
}

