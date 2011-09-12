\name{ex1.data}
\alias{ex1.data}

\docType{data}
\title{
	Example data frame
}
\description{
   Simulated data from a Gaussian contamination model
}
\usage{data(ex1.data)}
\format{
  A data frame with 500 observations on the following 2 variables.
  \describe{
    \item{\code{X1}}{numeric error-free variable}
    \item{\code{Y1}}{numeric contaminated variable}
  }
}
\details{
Non contaminated data have been simulated by Gaussian contamination model with one error-free covariate (X1) and one
contaminated variable (Y1) with parameters B=(-0.26, 1.26), sigma=1.21, w=0.05, lambda=10.
}

\examples{
data(ex1.data)

}
\keyword{datasets}
