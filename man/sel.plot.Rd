\name{sel.plot}
\alias{sel.plot}

\title{
   Scatterplot with information about outliers and influential errors
}
\description{ 
  In addition to a standard scatterplot, outliers and influential errors are highlighted.  
}
\usage{
sel.plot (data, vars=1:2, outl = rep(0, nrow(data)), sel = rep(0, nrow(data)),
          log = TRUE, n.identify=0, file=NULL, title=NULL)
}

\arguments{
  \item{data}{
    named matrix or data frame containing at least the coordinates of points
	}
	\item{vars}{
	vector with the names or column numbers of the two variables to plot
	}
  \item{outl}{
	vector identifying outliers (1 or TRUE means outlier)
	}
  \item{sel}{
	vector identifying influential errors (1 or TRUE means influential error)
	}

  \item{log}{
     if TRUE logarithm of \code{data[,vars]} are plotted
	}
   \item{n.identify}{
     number of points to be identified on the scattrerplot. Corresponding data are printed on console or file (if a file name is specified)
	} 
	\item{file}{
     name of the output file. If \code{n.identify} is equal 0 the graphic is saved in a jpeg file. If \code{n.identify} is greater than 0 data rows corresponding selected points are saved in a csv file
	} 
   \item{title}{
     an overall title for the plot
	} 	
}
\details{
  The scatterplot contains the first varaible plotted against the second. Outliers are represented as blue circles,
  influential errors as red circles and points that are both outlier and influential error as cyan circles.
  
}
\examples{

    data(ex2.data)
    par.joint <- ml.est(y=ex2.data)
    sel <- sel.edit(y=ex2.data, ypred=par.joint$ypred)	
    sel.plot(ex2.data,outl=par.joint$outlier, sel=sel[,"sel"], title="EXAMPLE 2")
  \dontrun{   
    sel.plot(ex2.data,outl=par.joint$outlier, sel=sel[,"sel"], title="EXAMPLE 2", n.identify=3)
    }
    
}    