#'  Simulation from the 2D-MxARMA Model
#'
#'  Generates a bivariates time series from an ARMA-type bidimensional model whose innovations follow
#' the Maxwell distribution. The user can simulate pure autoregressive (AR),
#' pure moving average (MA), or mixed ARMA dynamics depending on which parameters
#' are provided.
#'
#' @param n Integer. Number of rows to simulate.
#' @param k Integer. Number of coluns to simulate.
#' @param alpha  Numeric. Intercept of the linear predictor.
#' @param phi Numeric matrix Autoregressive coefficients. If \code{NULL} (default),
#'   no AR structure is included.
#' @param theta Numeric matrix Moving average coefficients. If \code{NULL} (default),
#'   no MA structure is included.
#'
#' @details
#' The MxARMA(p,q) model assumes that the conditional mean of the process is linked
#' to a linear predictor via a logarithmic link function. The innovations are assumed
#' to follow a Maxwell distribution, generated internally with \code{rmax()}.
#'
#' When only \code{phi} is supplied, an MxAR(p) model is simulated. When only
#' \code{theta} is supplied, an MxMA(q) model is generated. When both are supplied,
#' an MxARMA(p,q) model is simulated.
#'
#' @return A numeric vector of length \code{n} containing the simulated time series.
#'
#' @author Arthur Hintz
#'
#' @seealso \code{\link{rmax}}
#'
#' @examples
#' # ARMA(1,1) Example
#' alpha <-  0.2
#' phi <-  matrix(c(0.4,0.3, 0.1, 0), ncol = 2, byrow = T)
#' theta <- matrix(c(0.3,0.15, 0.05, 0), ncol = 2, byrow = T)
#' n <- 10; k <- 10
#'
#' Y <- mxarma2d.sim(n, k, alpha, phi, theta)
#'
#' @export
mxarma2d.sim <- function(n, k, alpha = 0.0, phi = NULL, theta = NULL) {


  if (!is.numeric(n) || n <= 0) {
    stop("n must be a positive number.")
  }
  if (!is.numeric(phi) && !is.numeric(theta)) {
    stop("phi or theta must be numeric vectors.")
  }
  if (any(is.null(phi)) && any(is.null(theta))) {
    stop("At least phi or theta must be specified.")
  }

  if (!is.numeric(alpha) && length(alpha) != 1) {
    stop("alpha must be a single value")
  }

  isphi   <- !is.null(phi)
  istheta <- !is.null(theta)


  if (isphi == T && !is.matrix(phi)) {
    stop("phi deve ser uma matriz")
  }

  if (istheta == T && !is.matrix(theta)) {
    stop("theta deve ser uma matriz")
  }

#________________AR___________________________________
    if (isphi == T & istheta == F) {

      ar = length(phi) - 1

      p = sqrt(ar + 1) - 1

      m = p

      phi[(p + 1), (p + 1)] <- 0 #Forçar o phi do elemento max(p)=max(q) ser 0

      ynew  =  matrix(alpha, ncol = (k + m),nrow = (n + m))

      y = matrix(0, ncol=(k+m),nrow=(n+m))
      eta = mu = y

      for (i in (m+1):(n+m))
        {
        for (j in (m+1):(k+m))
          {
          ynew1 = ynew[(i-m):i,(j-m):j]  # filtrando intervalo de interesse

          eta[i,j]  = alpha + sum(phi*ynew1)
          mu[i,j] = exp(eta[i,j])
          y[i,j] = MxARMA::rmax(1, mu[i,j])
          ynew[i,j] = log(y[i,j])
        }
      }

      y = y[(m+1):(n+m),(m+1):(k+m)]

      mu = mu[(m+1):(n+m),(m+1):(k+m)]

      return(list(y = y, mu = mu))
    }

#________________MA___________________________________
  if (isphi == F & istheta == T) {

    ma = length(theta) - 1

    q = sqrt(ma + 1) - 1

    m = q

    theta[(q + 1), (q + 1)] <- 0

    ynew  = matrix(alpha, ncol = (k + m), nrow = (n + m))

    y = matrix(0, ncol = (k + m), nrow = (n + m))

    eta = error = mu = y

    for (i in (m+1):(n+m)){

      for (j in (m+1):(k+m)){

        error_new <- error[(i-q):i,(j-q):j]

        eta[i,j]  = alpha + sum(theta*error_new)
        mu[i,j] = exp(eta[i,j])
        y[i,j] = MxARMA::rmax(1, mu[i,j])
        ynew[i,j] = log(y[i,j])

        error[i,j] = ynew[i,j]-log(mu[i,j])
      }
    }

    y = y[(m+1):(n+m),(m+1):(k+m)]

    return(y)
  }

#________________ARMA___________________________________
  if(isphi == T & istheta == T){

    ar = length(phi) - 1
    p = sqrt(ar + 1) - 1

    ma = length(theta) - 1
    q = sqrt(ma + 1) - 1

    m = max(p,q)

    phi[(p + 1), (p + 1)] <- 0
    theta[(q + 1), (q + 1)] <- 0

    ynew = matrix(alpha, ncol = (k + m), nrow = (n + m))

    y = matrix(0, ncol = (k + m), nrow = (n + m))

    eta = error = mu = y

    for (i in (m+1):(n+m)) {

      for (j in (m+1):(k+m)) {

        ynew1 = ynew[(i-m):i,(j-m):j]

        error_new <- error[(i-q):i,(j-q):j]

        eta[i,j]  = alpha + sum(phi * ynew1) + sum(theta * error_new)
        mu[i,j] = exp(eta[i,j])
        y[i,j] = MxARMA::rmax(1, mu[i,j])
        ynew[i,j] = log(y[i,j])

        error[i,j] = ynew[i,j]-log(mu[i,j])
      }

    }

    y = y[(m+1):(n+m),(m+1):(k+m)]
    mu = mu[(m+1):(n+m),(m+1):(k+m)]

    return(list(y = y, mu = mu))

  }
}



#
# teste <- as.vector(t(phi))
# teste
#
# phi <- matrix(c(0.4,0.3,0.1, 0), ncol = 2, nrow = 2, byrow = T)
# phi
# theta <- matrix(c(0.3,0.15,0.1, 0), ncol = 2, nrow = 2, byrow = T)
# theta
#
# alpha = 0.4
# n = 100
# k = 100
#
#
# resu <- mxarma2d.sim(n = n, k = k, alpha = alpha, phi = phi, theta = theta)
# y = resu$y
# mu = resu$mu
# y=resu
# hist(y)
#
# plot.ts(y[1,])
# acf(y)
#
# summary(y)
#
#
#
#
# hist(y, prob = TRUE, breaks = 30, ylim = c(0,0.8), xlim = c(0,13))
# curve(dmax(x, mean(mu)), add = TRUE, col = "blue", lwd = 2)
#
#
#
#
# ynew <- matrix(1:992, ncol = (k+2*m), nrow = (n+m))
