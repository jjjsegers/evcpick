\name{evcCFG}
\alias{evcCFG}
\title{Function \code{evcCFG}}
\usage{
  evcCFG(t, x = NULL, pdata = pseudo(x), abc = TRUE)
}
\arguments{
  \item{t}{The arguments (between 0 and 1) of the Pickands
  dependence function.}

  \item{x}{The matrix of bivariate block maxima, a k-by-2
  matrix.}

  \item{pdata}{Alternatively, the data may be directly
  specified as the pseudo-observations, i.e. the ranks
  divided by k+1.}

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
  C. Genest and J. Segers (2009), Annals of Statistics.
}

