
check.vars <- function (y, x, model="LN", parent="pred.y") {
#print(parent.env(env) )  vedere come posso conoscere il nome della funzione chiamante (altrimenti lo passo)

   lista  <- list(ret=0, msg.err=NULL, y=NULL, x=NULL)
#   ERRORI BLOCCANTI
   if (!is.numeric(as.matrix(y))  )  {
      lista$msg.err <- "Variables must be numeric "
      lista$ret <- -9
      
   }
   else if (!inherits(y,c("data.frame", "matrix","numeric", "integer"))  )  {
      lista$msg.err <- "Variables must be supplied as a matrix or as a data frame "
      lista$ret <- -9
   }
   else if ( model == 'LN' && sum(y<0, na.rm=TRUE) > 0 )  {
        lista$msg.err <-"Negative values are not allowed if model is LN\n"
        lista$ret <- -9
   }
#------------------------------------------------------------------------------ 
#  Aggiunte il 25/02/2016  
#------------------------------------------------------------------------------    
   else if ( sum(is.nan(y), na.rm=TRUE) > 0 || sum(is.infinite(y), na.rm=TRUE) > 0)  {
        lista$msg.err <-" +/- Inf and NaN values are not allowed \n"
        lista$ret <- -9
   }


   y <- as.matrix (y)
   n <- nrow(y)
   if (!is.null(x)) {
      x<- as.matrix(x)
      if (!is.numeric(x) )  {
       lista$msg.err <- "Covariates must be numeric "
       lista$ret <- -9
     }   
     else if (!inherits(x,c("data.frame", "matrix","numeric", "integer")) )  {
      lista$msg.err <- "Covariates must be supplied as a matrix or as a data frame "
      lista$ret <- -9
      }
     else if (sum(is.na(x)) > 0)  {
        lista$msg.err <-"Covariates can not have missing values "
        lista$ret <- -9
     } else if (nrow(x) != nrow(y))  {
        lista$msg.err <-"Variables y and x must have the same number of rows"
        lista$ret <- -9
     }
     else if ( model == 'LN' &&  sum(x<0) != 0 )  {
        lista$msg.err <-"Negative values are not allowed if model is LN\n"
        lista$ret <- -9
     }
   }
    
   
   if (lista$ret == -9)
       return (lista)
#    WARNING


#   PREPARAZIONE DATI IN BASE AL MODELLO
  modificati <- 0 
  if (model == "LN") {  #lognormal model
      for (i in 1:ncol(y)) {
        #y[,i] <- sostituisci_zeri(y[,i])
        n.zeri <- length(which(y[,i] == 0))
        if (n.zeri > 0)  {
          ind0<-which(y[,i]==0)
          new.val <- min(y[,i][-ind0])/2
          y[ind0,i]<- new.val
          modificati <- modificati + n.zeri
        } 
      }
      if (parent!="ml.est" && modificati > 0)  {
          lista$ret <- lista$ret+1
          lista$msg.err <-  rbind( lista$msg.err,                 
                paste(modificati," response variable values (%", round(modificati*100/n, 2),
                       ") equal 0 are substituted by the half minimum of the corresponding variable\n",sep="") )
      }
      
      y <- as.matrix(log(y))
      
      modificati <- 0
      if (!is.null(x))  {
         x <- as.matrix (x)
         for (i in 1:ncol(x)) {
       
         n.zeri <- length(which(x[,i] == 0))
         if (n.zeri > 0)  {
          ind0<-which(x[,i]==0)
          new.val <- min(x[,i][-ind0])/2
          x[ind0,i]<- new.val
          modificati <- modificati + n.zeri
        } 
      }
  
      if (parent!="ml.est" && modificati > 0)   {                    
           lista$ret <- lista$ret+1
           lista$msg.err <- rbind( lista$msg.err,                  
                paste(modificati," covariate values (%", round(modificati*100/length(x), 2),
                     ") equal 0 are substituted by the half minimum of the corresponding variable\n", sep="") )
          }
              
        x <- log(x)
      }
      else
        x <- NULL
  } else {   # normal model
     y <- as.matrix(y)
     x <- x
  }
   lista$y <- y
   lista$x <- as.matrix(cbind(rep(1,n),x))
   return (lista)
 }
 
 
  