\name{ex2.data}
\alias{ex2.data}
\docType{data}
\title{
Example Data for package SeleMix
}
\description{
  Simulated data from a Gaussian contaminated model
}
\usage{data(ex2.data)}
\format{
   A data frame with 500 observations on the following 2 variables.
  \describe{
    \item{\code{Y1}}{first numeric contaminated variable}
    \item{\code{Y2}}{second numeric contaminated variable}
  }
}
\details{
Data have been simulated by Gaussian contamination model with two
contaminated variables (Y1,Y2) with parameters 
B=(1.03, 0.96), sigma=matrix(c(1.22, 1.42,1.42, 2.89),2,2), w=0.05, lambda=10.
}

\examples{
data(ex2.data)

}
\keyword{datasets}
