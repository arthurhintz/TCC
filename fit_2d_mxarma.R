#' 2D-MxARMA fitting without regressors
#'
#' Used for fitting the Maxwell auto-regressive moving averages  bidimensional model without regressors. The model is given by
#'
#' \deqn{\eta_t = \alpha + \sum_{j\in ar}\phi_j \log(y_{t-j}) +\sum_{j\in ma}\theta_jr_{t-j}}
#' Where \itemize{
#' \item{\eqn{y} are the variables}
#' \item{\eqn{\alpha} is the intercept}
#' \item{\eqn{ar} are the indices for the auto-regression}
#' \item{\eqn{ma} are the indices for the moving-averages}
#' \item{\eqn{\phi} are the auto-regression coefficients}
#' \item{\eqn{\theta} are the moving-averages coefficients}
#' \item{\eqn{r} are the errors}}
#'
#'
#' @param y The matrix of variables to fit
#' @param p are the indices for the auto-regression.
#' @param q are the indices for the moving-averages
#'
#'
#'@examples
#'
#'
#'@importFrom stats optim
#'@importFrom stats qnorm
#'@importFrom stats pnorm
#'
#'@export
mxarma2d.fit <-  function(y, p = NA, q = NA) {


  if (is.matrix(y) == F) stop("y must be a numeric matrix")

  if (nrow(y) != ncol(y)) stop('y must be a square matrix')


  maxit1 = 1000

  isar <- !is.na(p)
  isma <- !is.na(q)

  ar <- seq_along(p)
  ma <- seq_along(q) # corrigir para NA
  #p = max(ar) # AR order
  #q = max(ma) # MA order
  n = dim(y)[1] # n is the number of rows
  k = dim(y)[2] # k is the number of columns
  m = max(p,q,na.rm = T)

  ynew = log(y) # ynew é uma matriz se y é uma matriz

  # inicializacao dos parametros alpha e phi (beta)
  if (isar == T) # se tem componente ar
    {
      p1 = (p + 1)^2 - 1

      XX = c()

      for (i in (m+1):n)
      {
        for (j in (m+1):k)
        {
          #ynew <- matrix(1:100, ncol = 10)

          xx1 = as.vector(t(ynew[(i-ar):i,(j-ar):j]))
          #XX = rbind(XX,t(as.matrix(xx1[(length(xx1)):1])))
          XX = rbind(XX,xx1)

        }
      }

      # maybe to use 2:ncol(XX), to test

      #P = XX[,2:dim(XX)[2]]
      #Y = as.matrix(XX[,1])

      n_col <- ncol(XX)

      P <- as.matrix(XX[,-n_col])
      Y <- as.matrix(XX[,n_col])

      Z <- cbind(rep(1, (n - m) * (k - m)), P)

      x = as.matrix(Z)

      ajuste = lm.fit(x, Y)
      mqo = c(ajuste$coef)

    }else{

      Z <- as.matrix(rep(1, (n - m) * (k - m)))

      Y = as.vector(ynew[(m + 1):n, (m + 1):k])

      x = as.matrix(Z)
      ajuste = lm.fit(x, Y)
      mqo = c(ajuste$coef)
    }

#_____________________AR___________________________________

  if(isar == T && isma == F){

    p1 = (p+1)^2-1
    q1 = 0
    reg = c(mqo) # initializing the parameter values


    loglik = function(z)
    {
      alpha = z[1]
      phi = z[2:(p1+1)]
      eta = matrix(0, ncol = k, nrow = n)

      for (i in (m+1):n) {

        for (j in (m+1):k) {

          ynew1 <- as.vector(t(ynew[(i-ar):i,(j-ar):j]))
          ynew1 <- ynew1[-(p1+1)]

          eta[i,j]  = alpha + sum(phi * ynew1)
        }

      }

      mu = as.vector(exp((t(eta[(m+1):n,(m+1):k]))))
      y1 = as.vector(t(y[(m+1):n,(m+1):k]))

      ll <- log(32) - 2*log(pi) - 3*log(mu) + log(y1^2) + ((-4*y1^2)/(pi*mu^2))

      sum(ll)

    }

    escore <- function(z) {

      alpha = z[1]
      phi = z[2:(p1+1)]

      eta = matrix(0, ncol=k,nrow=n)

      for (i in (m+1):n)
      {
        for (j in (m+1):k)
        {
          ynew1 <- as.vector(t(ynew[(i-ar):i,(j-ar):j]))
          ynew1 <- ynew1[-(p1+1)]

          eta[i,j]  = alpha + sum(phi * ynew1)
        }

      }

      mu = exp(t(eta[(m+1):n,(m+1):k])) # tirei no formato vector, nao fazia sentido
      y1 = as.vector(t(y[(m+1):n,(m+1):k]))

      dmu <- as.vector(((8*y1^2)/(mu^3*pi))-(3/mu))

      mT = as.vector(diag(mu))

      mT = mu

      a = as.vector(rep(1,(n-m)*(k-m)))

      # fazer no formato de somatório para utilizar as matrizes dos parâmetros

      Ualpha = sum(a * mT * dmu)

      Uphi = c()
      for(i in 1:ncol(P)){

        Uphi[i] <- sum(P[,i] * mT * dmu)
      }
      #class(mT)
      #length(dmu)
      #Ualpha =  t(a) %*% mT %*% dmu
      #Uphi =    t(P) %*% mT %*% dmu

      rval = c(Ualpha,Uphi)
    }

    names_par <- c("alpha", paste0("phi", 1:p1))

    opt = optim(reg, loglik,
                escore,
                method = "BFGS",
                hessian = TRUE,
                control = list(fnscale = -1
                              #maxit = 400, abstol = 1e-6,
                              #factr = 1e20
                              ))


    if (opt$conv != 0)
      warning("FUNCTION DID NOT CONVERGE!")

    z <- c()
    z$conv <- opt$conv
    coef <- opt$par
    names(coef) <- names_par
    z$coeff <- coef

    z$loglik = opt$value  ## AQUI

    alpha = coef[1]
    phi = coef[2:(p1+1)]

    z$alpha = alpha
    z$phi = phi

    etahat = matrix(0, ncol = k, nrow = n)

    for (i in (m+1):n)
    {
      for (j in (m+1):k)
      {
        ynew1 = ynew[(i-ar):i,(j-ar):j]
        etahat[i,j]  = alpha + sum(phi*ynew1)
      }
    }


    z$fitted = exp(etahat[(m+1):n,(m+1):k])
    z$etahat = etahat[(m+1):n,(m+1):k]

    y1 = as.vector(t(y[(m+1):n,(m+1):k]))

    muhat2 = as.vector(exp(t(etahat[(m+1):n,(m+1):k])))

    W = diag(((6)/(muhat2^2))*(muhat2^2))

    a = as.vector(rep(1,(n-m)*(k-m)))

    Kaa = t(a) %*% W %*% a
    Kpa = t(P) %*% W %*% a
    Kap = t(Kpa)
    Kpp = t(P) %*% W %*% P

    K = rbind(
      cbind(Kaa,Kap),
      cbind(Kpa,Kpp)
    )

  }

#_____________________MA___________________________________

  if(isar == F && isma == T){

    p1 = 0
    q1 = (q+1)^2-1
    reg = c(mqo[1],rep(0,q1)) # initializing the parameter values

    loglik = function(z)
      {
      alpha = z[1]
      theta = z[2:(q1+1)]

      eta = error = matrix(0, ncol = k, nrow = n)

      for (i in (m+1):n)
      {
        for (j in (m+1):k)
        {
          error_new <- as.vector(t(error[(i-ma):i,(j-ma):j]))
          error_new <- error_new[-(q1+1)]

          eta[i,j]  = alpha + sum(theta*error_new)

          error[i,j] = ynew[i,j] - eta[i,j]
        }

      }

      mu = as.vector(exp((t(eta[(m+1):n,(m+1):k]))))
      y1 = as.vector(t(y[(m+1):n,(m+1):k]))

      ll <- log(32) - 2*log(pi) - 3*log(mu) + log(y1^2) + ((-4*y1^2)/(pi*mu^2))

      sum(ll)

    }

    escore = function(z)
    {
      alpha = z[1]
      theta = z[2:(q1+1)] #definir como matriz, se definir como matriz posso mudar linha 273
      #theta = matrix(c(z[2:(q1+1)]), ncol = (q + 1), byrow = T)

      eta = error = matrix(0, ncol = k, nrow = n)

      for (i in (m+1):n)
      {
        for (j in (m+1):k)
        {

          error_new <-  as.vector(t(error[(i-ma):i,(j-ma):j]))
          error_new = error_new[-(q1+1)]

          eta[i,j]  = alpha + sum(theta*error_new)
          error[i,j] = ynew[i,j] - eta[i,j]
        }

      }

      mu = exp(t(eta[(m+1):n,(m+1):k])) # tirei no formato vector, nao fazia sentido
      y1 = as.vector(t(y[(m+1):n,(m+1):k]))

      dmu <- as.vector(((8*y1^2)/(mu^3*pi))-(3/mu))

      mT = as.vector(diag(mu))

      # dalpha
      deta.dalpha <-  matrix(0, nrow = n, ncol = k)

      for(i in (m+1):n)
      {
        for(j in (m+1):k)
        {
          deta.dalpha1 = as.vector(t(deta.dalpha[(i-ma):i,(j-ma):j]))
          deta.dalpha1 = deta.dalpha1[-(q1+1)]
          deta.dalpha[i,j] = 1 - sum(theta*deta.dalpha1)
        }
      }

      a = as.vector(deta.dalpha[(m+1):n,(m+1):k])

      #dtheta
      XX = c()

      for (i in (m+1):n)
      {
        for (j in (m+1):k)
        {
          xx1 = as.vector(t(error[(i-ma):i,(j-ma):j]))
          #XX = rbind(XX,t(as.matrix(xx1[(length(xx1)):1])))
          XX = rbind(XX,xx1)
        }
      }

      #R = rbind(rep(0,q1),XX[,2:dim(XX)[2]])
      n_col <- ncol(XX)
      rownames(XX) <- NULL
      R <- rbind(rep(0,q1),XX[,-n_col])

      deta.dtheta = matrix(0, ncol = q1, nrow = dim(R)[1])
      dsum.theta = matrix(0, ncol = q1, nrow = q1)

      for(i in m:dim(R)[1])
      {
        if( i == m)
        {
          deta.dtheta[i,] <- R[i,]
        }else{
          for(j in 1:q1)
          {
            dsum.theta[j,] = theta[j] * deta.dtheta[i-ma,]
          }

          deta.dtheta[i,] <- R[i,] - apply(dsum.theta,2,mean)
        }
      }
      rR <- deta.dtheta[-m,]
#==========/==========/==========/==========/==========/==========/==========/==========/
      Ualpha = sum(a * mT * dmu)

      Utheta = c()
      for(i in 1:ncol(R)){

        Utheta[i] <- sum(rR[,i] * mT * dmu)
      }

      #Ualpha =  - a %*% mT %*% dmu
      #Utheta =  - t(rR) %*% mT %*% dmu

      rval = c(Ualpha,Utheta)
    }

    names_par <- c("alpha", paste0("theta", 1:q1))

    opt = optim(reg, loglik,
                escore,
                method = "BFGS",
                hessian = TRUE,
                control = list(fnscale = -1
                               #maxit = 400, abstol = 1e-6,
                               #factr = 1e20
                ))


    if (opt$conv != 0)
      warning("FUNCTION DID NOT CONVERGE!")

    z = c()
    z$conv = opt$conv
    coef = opt$par
    names(coef) = names_par
    z$coeff = coef
    z$loglik = opt$value  ### AQUI

    alpha = coef[1]
    theta = coef[2:(q1 + 1)]

    z$alpha = alpha
    z$theta = theta

    etahat = errorhat = matrix(0, ncol = k,nrow = n)

    for (i in (m+1):n)
    {
      for (j in (m+1):k)
      {
        errorhat_new = as.vector(t(errorhat[(i-ma):i,(j-ma):j]))
        errorhat_new = errorhat_new[-(q1+1)]

        etahat[i,j]  = alpha + sum(theta * errorhat_new)
        errorhat[i,j] = ynew[i,j] - etahat[i,j]
      }
    }


    z$fitted = exp(etahat[(m+1):n,(m+1):k])
    z$etahat = etahat[(m+1):n,(m+1):k]
    z$errorhat = errorhat[(m+1):n,(m+1):k]

    y1 = as.vector(t(y[(m+1):n,(m+1):k]))
#==========/==========/==========/==========/==========/==========/==========/==========/
    # dalpha
    deta.dalpha <-  matrix(0, nrow = n, ncol = k)

    for(i in (m+1):n)
    {
      for(j in (m+1):k)
      {
        deta.dalpha1 = as.vector(t(deta.dalpha[(i-ma):i,(j-ma):j]))
        deta.dalpha1 = deta.dalpha1[-(q1+1)]
        deta.dalpha[i,j] = 1 - sum(theta*deta.dalpha1)
      }
    }

    a = as.vector(deta.dalpha[(m+1):n,(m+1):k])

    #dtheta
    XX = c()

    for (i in (m+1):n)
    {
      for (j in (m+1):k)
      {
        xx1 = as.vector(t(errorhat[(i-ma):i,(j-ma):j]))
        #XX = rbind(XX,t(as.matrix(xx1[(length(xx1)):1])))
        XX = rbind(XX,xx1)
      }
    }


    #R = rbind(rep(0,q1),XX[,2:dim(XX)[2]])
    n_col <- ncol(XX)
    rownames(XX) <- NULL
    R <- rbind(rep(0,q1),XX[,-n_col])

    deta.dtheta = matrix(0, ncol = q1, nrow = dim(R)[1])
    dsum.theta = matrix(0, ncol = q1, nrow = q1)

    for(i in m:dim(R)[1])
    {
      if( i == m)
      {
        deta.dtheta[i,] <- R[i,]
      }else{
        for(j in 1:q1)
        {
          dsum.theta[j,] = theta[j] * deta.dtheta[i-ma,]
        }
        deta.dtheta[i,]<- R[i,] - apply(dsum.theta,2,mean)
      }
    }
    rR <- deta.dtheta[-m,]

#==========/==========/==========/==========/==========/==========/==========/==========/

    muhat2 = as.vector(exp(t(etahat[(m+1):n,(m+1):k])))

    W = diag(((6)/(muhat2^2))*(muhat2^2))

    Kaa = t(a) %*% W %*% a
    Kta = t(rR) %*% W %*% a
    Kat = t(Kta)
    Ktt = t(rR) %*% W %*% rR

    K = rbind(
      cbind(Kaa,Kat),
      cbind(Kta,Ktt)
    )

  }

#_____________________ARMA_________________________________

  if(isar == T && isma == T){

    p1 = (p+1)^2-1
    q1 = (q+1)^2-1

    reg = c(mqo, rep(0,q1)) # initializing the parameter values

    loglik = function(z){

      alpha = z[1]
      phi = z[2:(p1+1)]
      theta = z[(p1+2):(p1+q1+1)]

      eta = error = matrix(0, ncol = k, nrow = n)

      for (i in (m+1):n)
      {
        for (j in (m+1):k)
        {
          ynew1 <- as.vector(t(ynew[(i-ar):i,(j-ar):j]))
          ynew1 <- ynew1[-(p1+1)]

          error_new <- as.vector(t(error[(i-ma):i,(j-ma):j]))
          error_new <- error_new[-(q1+1)]

          eta[i,j]  = alpha + sum(phi * ynew1) + sum(theta * error_new)
          error[i,j] = ynew[i,j] - eta[i,j]
        }
      }

      mu = as.vector(exp((t(eta[(m+1):n,(m+1):k]))))
      y1 = as.vector(t(y[(m+1):n,(m+1):k]))

      ll <- log(32) - 2*log(pi) - 3*log(mu) + log(y1^2) + ((-4*y1^2)/(pi*mu^2))

      sum(ll)
    }

    escore <- function(z){

      alpha = z[1]
      phi = z[2:(p1+1)]
      theta = z[(p1+2):(p1+q1+1)]

      eta = error = matrix(0, ncol = k, nrow = n)

      for (i in (m+1):n)
      {
        for (j in (m+1):k)
        {
          ynew1 <- as.vector(t(ynew[(i-ar):i,(j-ar):j]))
          ynew1 <- ynew1[-(p1+1)]

          error_new <- as.vector(t(error[(i-ma):i,(j-ma):j]))
          error_new <- error_new[-(q1+1)]

          eta[i,j]  = alpha + sum(phi * ynew1) + sum(theta * error_new)
          error[i,j] = ynew[i,j] - eta[i,j]
        }
      }

      mu = as.vector(exp(t(eta[(m+1):n,(m+1):k])))
      y1 = as.vector(t(y[(m+1):n,(m+1):k]))

      dmu <- as.vector(((8*y1^2)/(mu^3*pi))-(3/mu))

      #mT = as.vector(diag(mu))
      mT = as.vector(mu)
#==========/==========/==========/==========/==========/==========/==========/==========/
      # dalpha
      deta.dalpha = matrix(0, nrow = n, ncol=k)

      for(i in (m+1):n)
      {
        for(j in (m+1):k)
        {
          deta.dalpha1 = as.vector(t(deta.dalpha[(i-ma):i,(j-ma):j]))
          deta.dalpha1 = deta.dalpha1[-(q1+1)]
          deta.dalpha[i,j] = 1 - sum(theta*deta.dalpha1)
        }
      }
      a = as.vector(deta.dalpha[(m + 1):n,(m + 1):k])

      # dphi
      P1 = rbind(rep(0,p1),P)
      deta.dphi = matrix(0, ncol = p1, nrow = dim(P1)[1])
      dsum.phi = matrix(0, ncol = p1, nrow = q1)

      for(i in m:dim(P1)[1]) {
        if( i == m) {

          deta.dphi[i,] <- P1[i,]
        }else{

          for(j in 1:q1) {
            dsum.phi[j,] = theta[j] * deta.dphi[i-ma,]
          }
          deta.dphi[i,] <- P1[i,] - apply(dsum.phi,2,mean)
        }
      }
      rP <- deta.dphi[-m,]

      #dtheta
      XX = c()

      for (i in (m+1):n)
      {
        for (j in (m+1):k)
        {
          xx1 = as.vector(t(error[(i-ma):i,(j-ma):j]))
          #XX = rbind(XX,t(as.matrix(xx1[(length(xx1)):1])))
          XX = rbind(XX,xx1)
        }
      }

      #R = rbind(rep(0,q1),XX[,2:dim(XX)[2]])
      n_col <- ncol(XX)
      rownames(XX) <- NULL
      R <- rbind(rep(0,q1),XX[,-n_col])

      deta.dtheta = matrix(0, ncol = q1, nrow = dim(R)[1])
      dsum.theta = matrix(0, ncol = q1, nrow = q1)

      for(i in m:dim(R)[1])
      {
        if( i == m)
        {
          deta.dtheta[i,]<- R[i,]
        }else{
          for(j in 1:q1)
          {
            dsum.theta[j,] = theta[j] * deta.dtheta[i-ma,]
          }
          deta.dtheta[i,]<- R[i,] - apply(dsum.theta,2,mean)
        }
      }
      rR <- deta.dtheta[-m,]
#==========/==========/==========/==========/==========/==========/==========/==========/
      Ualpha = sum(a * mT * dmu)

      Uphi = c()
      for(i in 1:ncol(P)) {
        Uphi[i] <- sum(P[,i] * mT * dmu)
      }

      Utheta = c()
      for(i in 1:ncol(R)) {
        Utheta[i] <- sum(rR[,i] * mT * dmu)
      }

      #Ualpha =  - a %*% mT %*% dmu
      #Uphi =    - t(rP) %*% mT %*% dmu
      #Utheta =  - t(rR) %*% mT %*% dmu

      rval = c(Ualpha,Uphi,Utheta)
    }

    names_par <- c("alpha", paste0("phi", 1:p1), paste0("theta", 1:q1))

    opt = optim(reg, loglik,
                escore,
                method = "BFGS",
                hessian = TRUE,
                control = list(fnscale = -1
                               #maxit = 400, abstol = 1e-6,
                               #factr = 1e20
                ))


    if (opt$conv != 0)
      warning("FUNCTION DID NOT CONVERGE!")

    z = c()
    z$conv = opt$conv
    coef = opt$par
    names(coef) = names_par
    z$coeff = coef
    z$loglik = opt$value

    alpha = coef[1]
    phi = coef[2:(p1+1)]
    theta = coef[(p1+2):(p1+q1+1)]


    z$alpha = alpha
    z$phi = phi
    z$theta = theta

    etahat = errorhat = matrix(0, ncol = k, nrow = n)

    for (i in (m+1):n)
    {
      for (j in (m+1):k)
      {
        ynew1 <- as.vector(t(ynew[(i-ar):i,(j-ar):j]))
        ynew1 <- ynew1[-(p1+1)]

        errorhat_new <- as.vector(t(errorhat[(i-ma):i,(j-ma):j]))
        errorhat_new <- errorhat_new[-(q1+1)]

        etahat[i,j]  = alpha + sum(phi * ynew1) + sum(theta * errorhat_new)
        errorhat[i,j] = ynew[i,j] - etahat[i,j]
      }
    }

    z$fitted = exp(etahat[(m+1):n,(m+1):k])
    z$etahat = etahat[(m+1):n,(m+1):k]
    z$errorhat = errorhat[(m+1):n,(m+1):k]

    y1 = as.vector(t(y[(m+1):n,(m+1):k]))
#==========/==========/==========/==========/==========/==========/==========/==========/
    # dalpha
    deta.dalpha = matrix(0, nrow = n, ncol=k)

    for(i in (m+1):n)
    {
      for(j in (m+1):k)
      {
        deta.dalpha1 = as.vector(t(deta.dalpha[(i-ma):i,(j-ma):j]))
        deta.dalpha1 = deta.dalpha1[-(q1+1)]
        deta.dalpha[i,j] = 1 - sum(theta*deta.dalpha1)
      }
    }
    a = as.vector(deta.dalpha[(m + 1):n,(m + 1):k])

    # dphi
    P1 = rbind(rep(0,p1),P)
    deta.dphi = matrix(0, ncol = p1, nrow = dim(P1)[1])
    dsum.phi = matrix(0, ncol = p1, nrow = q1)

    for(i in m:dim(P1)[1]) {
      if( i == m) {

        deta.dphi[i,] <- P1[i,]
      }else{

        for(j in 1:q1) {
          dsum.phi[j,] = theta[j] * deta.dphi[i-ma,]
        }
        deta.dphi[i,] <- P1[i,] - apply(dsum.phi,2,mean)
      }
    }
    rP <- deta.dphi[-m,]

    ##dtheta
    XX = c()

    for (i in (m+1):n)
    {
      for (j in (m+1):k)
      {
        xx1 = as.vector(t(errorhat[(i-ma):i,(j-ma):j]))
        #XX = rbind(XX,t(as.matrix(xx1[(length(xx1)):1])))
        XX = rbind(XX,xx1)
      }
    }

    #R = rbind(rep(0,q1),XX[,2:dim(XX)[2]])
    n_col <- ncol(XX)
    rownames(XX) <- NULL
    R <- rbind(rep(0,q1),XX[,-n_col])
    deta.dtheta = matrix(0, ncol=q1,nrow=dim(R)[1])
    dsum.theta = matrix(0, ncol=q1,nrow=q1)

    for(i in m:dim(R)[1])
    {
      if( i == m)
      {
        deta.dtheta[i,]<- R[i,]
      }else{
        for(j in 1:q1)
        {
          dsum.theta[j,] = theta[j] * deta.dtheta[i-ma,]
        }
        deta.dtheta[i,]<- R[i,] - apply(dsum.theta,2,mean)
      }
    }
    rR <- deta.dtheta[-m,]
#==========/==========/==========/==========/==========/==========/==========/==========/
    muhat2 = as.vector(exp(t(etahat[(m+1):n,(m+1):k])))
    W = diag(((6)/(muhat2^2))*(muhat2^2))

    Kaa = t(a) %*% W %*% a
    Kpa = t(rP) %*% W %*% a
    Kap = t(Kpa)
    Kta = t(rR) %*% W %*% a
    Kat = t(Kta)
    Kpp = t(rP) %*% W %*% rP
    Kpt = t(rP) %*% W %*% rR
    Ktp = t(Kpt)
    Ktt = t(rR) %*% W %*% rR

    K = rbind(
      cbind(Kaa,Kap,Kat),
      cbind(Kpa,Kpp,Kpt),
      cbind(Kta,Ktp,Ktt)
    )
    }

#==========/==========/==========/==========/==========/==========/==========/==========/

  z$image = y
  y1 = y[(m+1):n,(m+1):k]
  z$fitted = exp(z$etahat)

  # residuals
  z$resid = qnorm(MxARMA::pmax(y1,z$fitted)) # quantile residuals

  vcov = chol2inv(chol(K))
  z$vcov = vcov

  stderror = sqrt(diag(vcov))
  z$stderror = stderror

  z$zstat = abs(z$coef/stderror)
  z$pvalues = 2*(1 - pnorm(z$zstat) )

  model_presentation = cbind(round(z$coef,4),round(z$stderror,4),
                             round(z$zstat,4),round(z$pvalues,4))
  colnames(model_presentation)=c("Estimate","Std. Error","z value","Pr(>|z|)")

  z$model = model_presentation

  return(z)
}

#resu <- mxarma2d.fit(y = y, p = 1, q = 1)
#resu
#resu$coeff





