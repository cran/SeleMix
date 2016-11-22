\name{ml.est}
\alias{ml.est}

\title{
    Fitting Contamination Model
}
\description{
    Provides ML estimates of a Gaussian contamination model.
}
\usage{
    ml.est (y, x=NULL, model = "LN", lambda=3,  w=0.05,
            lambda.fix=FALSE, w.fix=FALSE, eps=1e-7,  
            max.iter=500, t.outl=0.5, graph=FALSE)
}
\arguments{
  \item{y}{
    matrix or data frame containing the response variables
  }
  
  \item{x}{
    optional matrix or data frame containing the error free covariates
  }
  
  \item{model}{
    data distribution: LN = lognormal(default), N=normal
  }

  \item{lambda}{
    starting value for the variance inflation factor (default=3)
        }
    
  \item{w}{
    starting value for the proportion of contaminated data (default=0.05)
  }
  
  \item{lambda.fix}{ logical. TRUE if lambda is known
  }
  
  \item{w.fix}{ logical. TRUE if w is known
  }
  
  \item{eps}{
    epsilon : tolerance parameter for the log-likelihood convergence (default=1e-7)
  }
  
  \item{max.iter}{
    maximum number of EM iterations (default=500)
  }
  
  \item{t.outl}{
    threshold value for posterior probabilities of identifying outliers (default=0.5)
  }
  
  \item{graph}{
    logical. TRUE to display graphics (default=FALSE)
  }
}
\details{
    This function provides the parameter estimates of a contamination model where a set of \code{y}
	variables is assumed to depend on a (possibly empty) set of covariates (\code{x} variables) 
	through a mixture of two linear regressions with Gaussian residuals. The covariance matrices
	of the two mixture components are assumed to be proportional (the proportionality constant being 
	\code{lambda}). In case of no \code{x} variables a mixture of two Gaussian distribution is estimated.  
BIC and AIC scores (\code{bic.aic}) are returned corresponding to both standard Gaussian model
 and contamination model in order to help the user to avoid possible over-parametrisation.  

According to the estimated model parameters, a matrix of predictions of \sQuote{true} \code{y} values 
(\code{ypred}) is computed. To each unit in the dataset, a flag (\code{outlier}) is assigned taking 
value 0 or 1 depending on whether the posterior probability of being erroneous (\code{tau}) is 
greater than the user specified threshold  (\code{t.outl}). 

The model is estimated using complete observations. Missing values in the \code{x} variables are
 not allowed. However,  \code{y} variables can be partly observed. Robust predictions of \code{y} variables 
 are provided even when they are not observed.
 A vector of missing pattern (\code{pattern}) indicates which  item is observed and which is missing.  

In case the option \sQuote{model = LN}  is specified, each zero value is changed in 1e-7 and 
a warning is returned. 

In order to graphically monitor EM algorithm, a scatter plot is showed where outliers 
are depicted as long as they are identified. The trajectory of the \code{lambda} parameter 
is also showed until convergence. 

}

\value{
 \code{ml.est} returns a list containing the following components:\cr\cr
 \item{ypred }{matrix of predicted values for y variables}
 \item{B }{matrix of estimated regression coefficients}
 \item{sigma }{estimated covariance matrix}
 \item{lambda }{estimated variance inflation factor}
 \item{w }{estimated proportion of erroneous data}
 \item{tau }{vector of posterior probabilities of being contaminated}
 \item{outlier }{1 if the observation is classified as an outlier, 0 otherwise}
 \item{n.outlier }{total of outlier observations}
 \item{pattern }{vector of non-response patterns for y variables: 0 = missing, 1 = present value}
 \item{is.conv }{logical value: TRUE if the EM algorithm has converged}
 \item{n.iter }{number of iterations of EM algorithm}
 \item{sing }{if TRUE iteration are stopped because there is an almost perfect fit } 
 \item{bic.aic }{Bayesian Information Criterion and Akaike Information Criterion for contaminated 
 and non contaminated Gaussian models}
}
\references{
     Di Zio, M., Guarnera, U. (2013) \emph{"A Contamination Model for Selective Editing"}, 
        Journal of Official Statistics. Volume 29, Issue 4, Pages 539-555 (\url{http://www.degruyter.com/view/j/jos.2013.29.issue-4/jos-2013-0039/jos-2013-0039.xml}). \cr
        
     Buglielli, M.T., Di Zio, M., Guarnera, U. (2010) \emph{"Use of Contamination Models for Selective Editing"}, 
        European Conference on Quality in Survey Statistics Q2010, Helsinki, 4-6 May 2010
}
\author{
    M. Teresa Buglielli <bugliell@istat.it>, Ugo Guarnera <guarnera@istat.it>
}

\examples{

# Parameter estimation with one contaminated variable and one covariate
        data(ex1.data)
        ml.par <- ml.est(y=ex1.data[,"Y1"], x=ex1.data[,"X1"], graph=TRUE) 
        str(ml.par)
        sum(ml.par$outlier)  # number of outliers
# Parameter estimation with two contaminated variables and no covariates
        data(ex2.data)
        par.joint <- ml.est(y=ex2.data, x=NULL, graph=TRUE)  
        sum(par.joint$outlier)  # number of outliers           
}
