# Simulação
#
source("simu_2d_mxarma.R")
source("fit_2d_mxarma.R")

library(doParallel)
library(doRNG)

        
        
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

nvalores = c(10, 30, 50, 80, 150, 200, 400)
kvalores = c(10, 30, 50, 80, 150, 200, 400)

vp <- c(alpha, phi_values[-4], theta_values[-4])
mf <- matrix(NA, ncol = length(vp), nrow = nrep)


#results <- array(NA, dim = c(length(vp), 4, length(nvalores)))

#colnames(results) <- c("parâmetros", "par_estimados", "viés_relativo", "EQM")
colnames(mf) <- c("alpha", paste0("phi", 1:p1), paste0("theta", 1:q1))


# paralelismo

# Infraestrutura paralela
n_cores <- parallel::detectCores() - 2L
cl <- parallel::makeCluster(n_cores)
on.exit(parallel::stopCluster(cl), add = TRUE)

doParallel::registerDoParallel(cl)
# reprodutibilidade independente do número de núcleos
doRNG::registerDoRNG(1248)   

parallel::clusterExport(cl, c("mxarma2d.fit", "mxarma2d.sim"), 
                        envir = environment())


for (j in seq_along(nvalores)) {
  
  k <- kvalores[j]
  n <- nvalores[j]
  
  start_time <- Sys.time()
  
  #if (n < 100) nrep <- 4000 else nrep <- 2000
  
foreach::foreach(i = 1:nrep, .combine = rbind) %dopar% {
    
    rasu <- mxarma2d.sim(n = n, k = k, alpha = alpha,
                         phi = phi, theta = theta)
    
    y <- rasu$y
    
    fit <- mxarma2d.fit(y = y, p = p, q = q)
    mf[i,] <- fit$coeff
    
    end_time <- Sys.time()
    elapsed <- round(difftime(end_time, start_time, units = "secs"), 2)
    
    cat(
      sprintf("[%s] n = %d (%d×%d) | rep = %d | tempo = %s s\n",
              format(Sys.time(), "%H:%M:%S"),
              j, n, k, i, elapsed)
    )
    
  }
  
  write.table(mf, paste0("simu_2DMxARMA_nk", n, ".txt"),
              row.names = FALSE, quote = FALSE)
}

cat("Simulação concluída às", format(Sys.time(), "%H:%M:%S"), "\n")



# par_est <- apply(mf, 2, mean)
# vviesm <- (par_est - vp)
# viest <- (par_est - vp) / vp
# vvarm <- apply(mf, 2, var)
# veqmm <- vviesm^2 + vvarm
# 
# results[ , ,j] <- cbind(vp, round(par_est, digits = 4),
#                         round(viest, digits = 5), round(veqmm, digits = 3))
