#' Identification of outliers based on a mixture of 2 Gaussian distributions
#'
#' @param y matrix or data.frame with the response variable(s), 
#' assumed to be affected by outliers 
#' @param x matrix or data.frame with the predictors, assumed 
#' to be free of errors (no missing values allowed)
#' @param model assumed distribution N= Gaussian while LN stands for
#' LogNormal (default)
#' @param lambda numeric indicated the variance inflation factor (default is 3)
#' @param w numeric indicating the (apriori) proportion of contaminated (outliers) data
#' (default is 0.05)
#' @param lambda.fix logical, whether lambda parameter should be kept
#' fixed or can change during the fitting procedure (default FALSE)
#' @param w.fix ogical, whether  w parameter should be kept
#' fixed or can change during the fitting procedure (default FALSE)
#' @param eps tolerance for assessing convergence of the iterative procedure
#' @param max.iter numeric indicating the maximum number of iteration of the
#' estimation procedure based on EM algorithm (default is 500)
#' @param t.outl numeric threshold, between 0 and 1, for identifying possible 
#' errors, units whose estimated probability of being in the contaminated
#' distributions are >t.outl are identified as errors
#' @param graph logical, wheter a graph summarizing the iterative fitting 
#' procedure should be plotted (default FALSE)
#'
#' @returns
#' @export
#'
#' @examples
ml.est <- function (y, x=NULL, model = "LN", lambda=3,  
                    w=0.05, lambda.fix=FALSE, w.fix=FALSE, 
                    eps=1e-7, max.iter=500, t.outl=0.5, graph=FALSE)
{
  ris<-list(
    ypred = NA,  
    B=NA, 
    sigma=NA, 
    lambda=Inf,
    w=NA, 
    tau=NA, 
    outlier = NA,
    n.outlier =0,
    pattern= NA,
    is.conv = NA,
    n.iter =NA,
    sing=NA,
    bic.aic = NA,
    msg="",
    model=model
  )
#------------------------------------------------------------------------------
# copy of input data
#------------------------------------------------------------------------------
  memo.y <- y <- as.matrix(y)         
  memo.x <- x
#------------------------------------------------------------------------------
#        Checks on input data
#------------------------------------------------------------------------------  
  ind.NA<- which(rowSums(is.na(y)) >0)

#------------------------------------------------------------------------------  
  if (length(ind.NA) > 0 )    {
      warning(paste("Input matrix y contains", length(ind.NA), " (%",length(ind.NA)*100/nrow(y) ,
      ") rows with missing values not included in parameter's estimation\n" ))
      y <- y[-ind.NA,,drop=FALSE]
      if (!is.null(x)) {
	      x <-as.matrix(x) 
        x <- x[-ind.NA,,drop=FALSE]
	  }	  
   }
#------------------------------------------------------------------------------
#        checks on the input arguments
#------------------------------------------------------------------------------
  vars <- check.vars(y,x,model,parent="ml.est")

  if (vars$ret == -9) {
     stop(vars$msg.err) 
  }
  if (vars$ret != 0) {
     warning(vars$msg.err) 
  } 
  y <- as.matrix(vars$y)
  x <- as.matrix(vars$x) 
  
  #y <- as.matrix(y, dimnames = NULL)
   
  p <- ncol(y)
  n <- nrow(y)
  omega <- rep(1,n)  
   
#  x <- as.matrix(cbind (rep(1,n),x))     
  q <- ncol(x) 
#------------------------------------------------------------------------------  
  
#---------------------       DEFine variables       ---------------------
  conv <- FALSE
  continua <- TRUE
  lik<-NA
  lik0 <- 10
  oldlik <- 0
  iter <- 0
  sing <- FALSE
  
  if (ncol(x)+ ncol(y) < 3)   
      graph=FALSE 
  if (graph )   {
   lambda_all <- lambda

   if (ncol(y) >= 2)  { 
      lab  <- names(y)[1:2]
      Var <- y[,1:2]  
   }   
   else if (ncol(y) == 1 & ncol(x) > 1) {
      lab  <- c(names(x)[2],names(y)[1])
      Var <- cbind(x[,2],y[,1])  
   }     
   
   par(mfrow=c(2,1))
  }
 
# B ha q (ncol(x)) righe e p (ncol(y)) variabili  
 B <- solve((t(x) %*% x) + (10e-8* diag(rep(1,q)))) %*% t(x) %*% y  
 # B <- try(solve(t(x) %*% x) %*% t(x) %*% y, TRUE) 
 
  B0 <- B
  sigma <- (t(y - x%*%B) %*% (y - x%*%B)) / (n-1)
  
  sigma <- sigma + (10e-8* diag(rep(1,p)))
  sigma0 <- sigma
  sigma2 <- (1 + lambda) * sigma  
  w1 <- 1-w
#------------------------------------------------------------------------------  
#   estimates BIC for comparing models
#------------------------------------------------------------------------------  
# No. parameters for Gaussina distribution
    k1 <- ncol(x) * p + (p*(p+1))/2 # p=ncol(y) 
# No. parameters for contamined model
    k2 <- k1 + 2 - w.fix - lambda.fix  
    
 if (n < k2)  {
      warning(paste("Input data are fewer than the number of model parameters\n" ))
 }
#------------------------------------------------------------------------------      
# estima lokelihood
#------------------------------------------------------------------------------  
    dati <- cbind(x, y)
    norm.mv<-function(u){dmvnorm(u[q+1:p], t(B0)%*%u[1:q], sigma0, log=TRUE)}
    lik.n <- sum(apply(dati, 1, norm.mv))  
  BIC.n <- -2*lik.n + k1*log(n)  
  
#************************  Start EM  ************************************

  while (iter < max.iter & continua == TRUE)
  {
     iter <- iter + 1
#     print(paste("E-step",iter))
#***********************    E - STEP         ************************************      
    tau1 <- post.prob(y, x, B, sigma, w1, lambda)   
    tau2 <- 1 - tau1 
#***********************    M - STEP         ************************************    
#    print(paste("M-step",iter))
#***********************    estimate weights  **********************************
     if (!w.fix)
         w1 <- sum(tau1)/n;

#***********************        omega          **********************************
  
     omega <- as.vector(tau1 + tau2 / (1+lambda))
    
#***********************        B         **********************************
     appo <- t(x) %*% (omega * x)
     appo <- solve(appo)
     B <- appo %*% t(x) %*% (omega * y)
    
#***********************        sigma         **********************************

     dif <- y - x%*%B
     
     sigma <- (t(dif) %*% (omega * dif)) / n
     if (det(sigma)  < 10e-10)  {
         ris$sing <- TRUE
         ris$is.conv <- FALSE
         ris$msg <- "Covariance matrix quasi singular: essentially perfect fit"
         warning(ris$msg)
       
      }
     s1 <- solve(sigma)
#***********************        lambda        **********************************

#     q1 <- matrix(diag(dif %*% solve(sigma1) %*% t(dif)),n,1)        ##   DIM n,1
#     q2 <- matrix(diag(dif %*% solve(sigma2) %*% t(dif)),n,1)        ##   DIM n,1
     if (!lambda.fix)  {     
        
        appo <-  t(dif) %*% (as.vector(tau2) * dif) %*% s1
        lambda <- sum(diag(as.matrix(appo))) / (p * sum(tau2)) -1
 #       if (lambda > 1e+06) {
  #         sing <- TRUE  
   #        continua <- conv <- FALSE
    #       warning (paste("lambda =" ,lambda,": iterations stopped because of essentially perfect fit", sep=""))
     #      break
      #  }
       if (lambda < 0.5) {
           continua <- conv <- FALSE
           warning (paste("lambda parameter lower than 0.5. Possible lack of model identification.", sep=""))
           
           break
        }        
     }
    
  
#***********************    convergence      **********************************

     s2 <- s1 / (1 + lambda) 
     sigma2 <- (1+lambda)* sigma 
     q1 <- matrix(tensorizza (dif, s1),n,1) 
     q2 <- matrix(tensorizza (dif, s2),n,1) 
     rm (s1,s2)

     q1 <- -0.5*q1
     q2 <- -0.5*q2
     
     ll <- w1 * exp(q1) / sqrt(2*pi*det(sigma)) + (1-w1) * (exp(q2)) / sqrt(2*pi*det(sigma2))
     lik <- sum(log(ll))
     
     if (graph)   {     
             plot(Var,  col = "lightgrey", main= "EM IN ACTION...\n Identifying outliers",  xlab=lab[1], ylab=lab[2] )
             points(Var[tau2 > t.outl, ],pch=21,col="blue",bg=paste("cyan",sample(1:4,1),sep="")) 
             
             lambda_all <- c (lambda_all, lambda)
             plot( lambda_all, xlab="n. iterations", ylab="lambda")
      }
   
     BIC.mix <- -2*lik + k2*log(n) 
         
     continua <- (abs(lik-oldlik) > eps*abs(lik-lik0) ) 
     conv <- !continua     
     if (iter > round(max.iter/5) & BIC.n < BIC.mix ) {
           continua <- conv <- FALSE
           warning (paste("EM stopped because BIC value " ,BIC.mix,"for the contamination model is greater than BIC value",
                      BIC.n, " for the Gaussian model", sep=""))
           break
     }

     #alpha <- sqrt((lambda+1) )
     oldlik <- lik
     if (iter == 1)  
        lik0 <- lik
   
  }
 #************************ end EM  ************************************

    if (iter >= max.iter)  
     conv <- FALSE
#   estimate predicted values
   
      yprev <-  pred.y(y=memo.y, x=memo.x, B, sigma, 
                     lambda, w=1-w1, model = model, t.outl=t.outl)    
 
   ###############  estimate BIC e AIC  #############
  
     BIC.n <- -2*lik.n + k1*log(n)
     BIC.mix <- -2*lik + k2*log(n)
     AIC.n <- 2*k1 - lik.n   
     AIC.mix <- 2*k2 - lik   

     ris$ypred <- as.matrix(yprev[,1:(ncol(yprev)-3)])
     ris$B <- B 
     ris$sigma <- sigma
     ris$lambda <- lambda
     ris$w <- 1-w1
     ris$tau <- yprev$tau
     ris$outlier <- yprev$outlier
     ris$n.outlier <- sum(yprev$outlier)
     ris$pattern <- yprev$pattern
     ris$is.conv <- conv
     ris$n.iter <- iter
     ris$sing <- sing
     ris$bic.aic <- c(BIC.norm=BIC.n, BIC.mix=BIC.mix, AIC.norm=AIC.n, AIC.mix=AIC.mix)

     class(ris) <- c(class(ris), "mlest" )
     if (conv == FALSE)
        warning (paste("EM algorithm failed to converge: stop after", iter, "iterations"))
    ris
 }



