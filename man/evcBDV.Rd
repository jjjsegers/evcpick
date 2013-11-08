\name{evcBDV}
\alias{evcBDV}
\title{Function \code{evcBDV}}
\usage{
  evcBDV(t, x = NULL, pdata = pseudo(x), kappa = 0.5,
    gamma = 2/3, abc = TRUE)
}
\arguments{
  \item{t}{The arguments (between 0 and 1) of the Pickands
  dependence function.}

  \item{x}{The matrix of bivariate block maxima, a k-by-2
  matrix.}

  \item{pdata}{Alternatively, the data may be directly
  specified as the pseudo-observations, i.e. the ranks
  divided by k+1.}

  \item{kappa}{The shape parameter of the weight function
  in the definition of the estimator.}

  \item{gamma}{The truncation parameter in the definition
  of the estimator}

  \item{mbc}{If TRUE, apply a multiplicative boundary
  correction.}
}
\value{
  A vector of length \code{length(t)}.
}
\description{
  Function to compute the rank-based
  Capéràa-Fougères-Genest estimator of the Pickands
  dependence function of a bivariate extreme-value copula.
}
\references{
  A. Buecher, H. Dette, and S. Volgushev (2010), The Annals
  of Statistics.
}

