\name{sel.edit}
\alias{sel.edit}

\title{
        Influential Error Detection
}
\description{
        Computes the score function and identifies influential errors
}
\usage{
       sel.edit (y, ypred, wgt=rep(1,nrow(as.matrix(y ))), 
                 tot=colSums(ypred * wgt), t.sel=0.01) 
}

\arguments{
  \item{y}{
        matrix or data frame containing the response variables 
        }

  \item{ypred}{
        matrix of predicted values for y variables
        }

  \item{wgt}{
        optional vector of sampling weights (default=1)
        }
\item{tot}{
        optional vector containing reference estimates of totals for the y variables. 
                If omitted, it is computed as the (possibly weighted) sum of predicted values 
        }
                
  \item{t.sel}{
        threshold value for selective editing (default=0.01)   
        }

}
\details{
This function ranks observations (\code{rank}) according to the importance of their potential errors. 
The order is made with respect to the global score function values (\code{global.score}).  
The function also selects the units to be edited (\code{sel}) so that the expected residual error of
all variables is below a prefixed level of accuracy (\code{t.sel}). 
The global score (\code{global.score})  is the maximum of the local scores computed for each variable 
(\code{y1.score, y2.score,\dots}). 
The local scores are defined as a weighted (\code{weights}) absolute difference between the observed 
(\code{y1, y2,\dots}) and the predicted values (\code{y1.p, y2.p,\dots}) standardised with respect to 
the reference total estimates (\code{tot}).\cr

The selection of the units to be edited because affected by an influential error (\code{sel=1}) is 
made according to a two-step algorithm:\cr
1) order the observations with respect to the \code{global.score} (decreasing order);\cr
2) select the first k units such that, from the (k+1)th to the last observation, all the 
residual errors (\code{y1.reserr, y2.reserr,\dots}) for each variable are below \code{t.sel}.\cr
\cr
The function provides also an indicator function (\code{y1.sel, y2.sel,\dots}) reporting 
which variables contain an influential errors in a unit selected for the revision. 
 
       }
\value{
  \code{sel.edit} returns a data matrix containing the following columns:\cr

    \item{y1, y2,\dots}{observed variables}
    \item{y1.p, y2.p,\dots}{predictions of y variables}

    \item{weights}{sampling weights}
    \item{y1.score, y2.score,\dots}{local scores}
    \item{global.score}{global score}
    \item{y1.reserr, y2.reserr,\dots}{residual errors}
    \item{y1.sel, y2.sel,\dots}{influential error flags}       
    \item{rank}{rank according to global score}
    \item{sel}{1 if the observation contains an influential error, 0 otherwise} 
}
\references{
    Di Zio, M., Guarnera, U., Luzi, O. (2008) \emph{"Contamination models for the detection of outliers and influential errors in continuous multivariate data"},
	UNECE Work Session on Statistical Data Editing, Vienna, 21-23 Aprile 2008 (\url{http://www.unece.org/stats/documents/2008.04.sde.htm}).\cr
	
    Buglielli, M.T., Di Zio, M., Guarnera, U. (2010) \emph{"Use of Contamination Models for Selective Editing"}, 
    European Conference on Quality in Survey Statistics Q2010, Helsinki, 4-6 May 2010.
}
\author{
    M. Teresa Buglielli <bugliell@istat.it>, Ugo Guarnera <guarnera@istat.it>
}
\examples{
# Example 1
    data(ex1.data)
    ml.par <- ml.est(y=ex1.data[,"Y1"], x=ex1.data[,"X1"])
    sel <- sel.edit(y=ex1.data[,"Y1"], ypred=ml.par$ypred)
    head(sel)
    sum(sel[,"sel"])
# orders results for decreasing importance of score     
    sel.ord <- sel[order(sel[,"rank"]),  ] 
# adds columns to data
    ex1.data <- cbind(ex1.data, tau=ml.par$tau, outlier=ml.par$outlier,
                      sel[,c("rank", "sel")])
# plot of data with outliers and influential errors 
    sel.pairs(ex1.data[,c("X1","Y1")],outl=ml.par$outlier, sel=sel[,"sel"])
# Example 2
    data(ex2.data)
    ml.par <- ml.est(y=ex2.data)
    sel <- sel.edit(y=ex2.data, ypred=ml.par$ypred)	
    sel.pairs(ex2.data,outl=ml.par$outlier, sel=sel[,"sel"])
}
