\name{sel.pairs}
\alias{sel.pairs}

\title{
   Scatterplot Matrix
}
\description{ 
  A scatterplot matrix with outlier and influential errors is produced.
}
\usage{
sel.pairs(x, outl = rep(0, nrow(x)), sel = rep(0,nrow(x)), 
          labs = NULL, log = TRUE, legend=TRUE)
}

\arguments{
  \item{x}{
   matrix or data frame of the coordinates of points
	}
  \item{outl}{
	vector identifying outliers (1 or TRUE means outlier)
	}
  \item{sel}{
	vector identifying influential errors (1 or TRUE means influential error)
	}
  \item{labs}{
	names of the variables 
	}
  \item{log}{
     if TRUE logarithm of \code{x} are plotted
	}
   \item{legend}{
     if TRUE a legend is added to first boxplot
	} 	
}
\details{
  The ijth scatterplot contains x[,i] plotted against x[,j]. Outliers are represented as blue circles,
  influential errors as red circles and points that are both outlier and influential error as cyan circles.
}
