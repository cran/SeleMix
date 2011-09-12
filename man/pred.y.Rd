\name{pred.y}
\alias{pred.y}

\title{
    Prediction of y variables  
}
\description{
    Provides predictions of y variables according to a Gaussian contamination model
}
\usage{
    pred.y (y, x=NULL, B, sigma, lambda, w, model="LN", t.outl=0.5)
}

\arguments{
  \item{y}{
    matrix or data frame containing the response variables
  }
  
  \item{x}{
    optional matrix or data frame containing the error free covariates
  }
  
  \item{B}{
    matrix of regression coefficients
  }
  \item{sigma}{
        covariance matrix
    }
  \item{lambda}{
        variance inflation factor
    }
  \item{w}{
        proportion of erroneous data
    } 
 
 \item{model}{
        data distribution: LN = lognormal(default), N=normal
    }
 \item{t.outl}{
        threshold value for posterior probabilities of identifying outliers (default=0.5)
     }
}

\details{
This function provides expected values of a set of variables (\code{y1.p,y2.p,\dots} ) according 
to a mixture of two regression models with Gaussian residuals (see \code{\link{ml.est}}). If no covariates are available 
(\code{x} variables), a two component Gaussian mixture is used. Expected values (predictions) are computed
 on the base of a set of parameters of appropriate dimensions (\code{B, sigma, lambda,w}) and (possibly) 
 a matrix (or data frame) containing the error-free \code{x} variables. 

Missing values in the \code{x} variables are not allowed. However, robust predictions of \code{y} variables are 
also provided when these variables are not observed. A vector of missing pattern (\code{pattern}) indicates 
which item is observed and which is missing.  

For each unit in the data set the posterior probability of being erroneous (\code{tau})  is computed and a
 flag (\code{outlier}) is provided taking value 0 or 1 depending on whether \code{tau}  is greater than the user 
 specified threshold  (\code{t.outl}).
      }
   
\value{
 \code{pred.y} returns a data frame containing the following columns:\cr\cr
 \item{y1.p,y2.p,\dots}{
  predicted values for y variables
 }
 \item{tau }{ posterior probabilities of being contaminated}
 \item{outlier }{ 1 if the observation is classified as an outlier, 0 otherwise }
 \item{pattern}{ non-response patterns for y variables: 0 = missing, 1 = present value}

}
\references{
     Buglielli, M.T., Di Zio, M., Guarnera, U. (2010) \emph{"Use of Contamination Models for Selective Editing"}, 
        European Conference on Quality in Survey Statistics Q2010, Helsinki, 4-6 May 2010
}
\author{
    M. Teresa Buglielli <bugliell@istat.it>, Ugo Guarnera <guarnera@istat.it>
}
\examples{

# Parameter estimation with one contaminated variable and one covariate
  data(ex1.data)
# Parameters estimated applying ml.est to \code{ex1.data} 
  B1 <- as.matrix(c(-0.152, 1.215))
  sigma1 <- as.matrix(1.25)
  lambda1 <- 15.5
  w1 <- 0.0479

# Variable prediction
  ypred <- pred.y (y=ex1.data[,"Y1"],  x=ex1.data[,"X1"], B=B1,
          sigma=sigma1, lambda=lambda1, w=w1, model="LN", t.outl=0.5)
# Plot ypred vs Y1
  sel.pairs(cbind(ypred[,1,drop=FALSE],ex1.data[,"Y1",drop=FALSE]),
            outl=ypred[,"outlier"])
}
