#' Function \code{pseudo}
#' 
#' Function to calculate pseudo-observations
#' 
#' @param x An n-by-d data matrix.
#' @return An n-by-d matrix of ranks divided by n.
#' @export
#' @examples
#' (x <- matrix(runif(10), nrow = 5, ncol = 2))
#' pseudo(x)
pseudo <- function(x) {
  d <- dim(x)[2]
  if( is.null(d) ) { stop("Number of variables needs to be larger than 1.") }
  n <- dim(x)[1]
  pdata <- apply(x, 2, rank) / (n+1)
  return(pdata)
}

#' Function \code{evcPickands}
#' 
#' Function to compute the rank-based Pickands estimator of the Pickands dependence function of a bivariate extreme-value copula.
#' 
#' @param t The arguments (between 0 and 1) of the Pickands dependence function.
#' @param x The matrix of bivariate block maxima, a k-by-2 matrix.
#' @param pdata Alternatively, the data may be directly specified as the pseudo-observations, i.e. the ranks divided by k+1.
#' @param abc If TRUE, apply an additive boundary correction.
#' @return A vector of length \code{length(t)}.
#' @export
#' @references
#' C. Genest and J. Segers (2009), The Annals of Statistics.
evcPickands <- function(t, x = NULL, pdata = pseudo(x), abc = TRUE) {
  l.t <- length(t)
  A.inv <- numeric(l.t)
  lpdata <- -log(pdata)
  for(i in 1:l.t) {
    A.inv[i] <- mean(pmin(lpdata[,1]/(1-t[i]),lpdata[,2]/t[i]))
  }
  A <- 1/A.inv
  if (abc) {
    A0inv <- mean(lpdata[,1])
    A1inv <- mean(lpdata[,2])
    A <- 1 / (A.inv - (1-t) * (A0inv-1) - t * (A1inv-1))
  } else {
    A <- 1 / A.inv
  }
  return( A )
}

#' Function \code{evcCFG}
#' 
#' Function to compute the rank-based Capéràa-Fougères-Genest estimator of the Pickands dependence function of a bivariate extreme-value copula.
#' 
#' @param t The arguments (between 0 and 1) of the Pickands dependence function.
#' @param x The matrix of bivariate block maxima, a k-by-2 matrix.
#' @param pdata Alternatively, the data may be directly specified as the pseudo-observations, i.e. the ranks divided by k+1.
#' @param mbc If TRUE, apply a multiplicative boundary correction.
#' @return A vector of length \code{length(t)}.
#' @export
#' @references
#' C. Genest and J. Segers (2009), Annals of Statistics.
evcCFG <- function(t, x = NULL, pdata = pseudo(x), abc = TRUE) {
  EulerGamma <- -0.5772156
  l.t <- length(t)
  A.log <- numeric(l.t)
  lpdata <- -log(pdata)
  for(i in 1:l.t) {
    xihat <- pmin(lpdata[,1]/(1-t[i]), lpdata[,2]/t[i])
    A.log[i] <- EulerGamma - mean( log(xihat) ) 
  }
  A <- exp(A.log)
  if (mbc) {
    A0log <- EulerGamma - mean( log(lpdata[,1]) )
    A1log <- EulerGamma - mean( log(lpdata[,2]) )
    A <- A * exp(- (1-t) * Alog0 - t * Alog1)
  }
  return(A)
}

empCop <- function(u, pdata) {
  data.vec <- as.numeric(pdata)
  k <- length(pdata[,1])
  d <- length(pdata[1,])
  u.vec <- as.numeric(u)
  l.u <- length(u[,1])
  res <- .C("empcop_c",
            as.double(data.vec),
            as.double(u.vec),
            as.integer(k),
            as.integer(d),
            as.integer(l.u),
            result = double(l.u))$result
  return(res)
}

#' Function \code{evcBDV}
#' 
#' Function to compute the rank-based Capéràa-Fougères-Genest estimator of the Pickands dependence function of a bivariate extreme-value copula.
#' 
#' @param t The arguments (between 0 and 1) of the Pickands dependence function.
#' @param x The matrix of bivariate block maxima, a k-by-2 matrix.
#' @param pdata Alternatively, the data may be directly specified as the pseudo-observations, i.e. the ranks divided by k+1.
#' @param kappa The shape parameter of the weight function in the definition of the estimator.
#' @param gamma The truncation parameter in the definition of the estimator
#' @param mbc If TRUE, apply a multiplicative boundary correction.
#' @return A vector of length \code{length(t)}.
#' @export
#' @references
#' A. Buecher, H. Dette, and S. Volgushev (2010), The Annals of Statistics.
evcBDV <- function(t, x = NULL, pdata = pseudo(x), kappa = 0.5, gamma = 2/3, abc = TRUE) {
  l.t <- length(t)
  A <- numeric(l.t)
  k <- dim(pdata)[1]
  f <- function(y, t) {
    u <- cbind(y^(1-t), y^t)
    Ctilde <- pmax(k^(-gamma), empCop(u, pdata))
    return(- (kappa+1) * y^kappa * log(Ctilde))
  }
  for (i in 1:l.t) {
    A[i] <- integrate(f, lower = 0, upper = 1, t = t[i])$value
  }
  if (abc) {
    A0 <- integrate(f, lower = 0, upper = 1, t = 0)$value
    A1 <- integrate(f, lower = 0, upper = 1, t = 1)$value
    A <- A - (1-t) * (A0 - 1) - t * (A1 - 1)
  }
  return(A)
}