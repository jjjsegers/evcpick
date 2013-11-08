\name{pseudo}
\alias{pseudo}
\title{Function \code{pseudo}}
\usage{
  pseudo(x)
}
\arguments{
  \item{x}{An n-by-d data matrix.}
}
\value{
  An n-by-d matrix of ranks divided by n.
}
\description{
  Function to calculate pseudo-observations
}
\examples{
(x <- matrix(runif(10), nrow = 5, ncol = 2))
pseudo(x)
}

