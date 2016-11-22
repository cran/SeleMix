ml.est <- function (y, x=NULL, model = "LN", lambda=3,  w=0.05, lambda.fix=FALSE, w.fix=FALSE, eps=1e-7, max.iter=500, t.outl=0.5, graph=FALSE)
{
#------------------------------------------------------------------------------
#           Individuazione degli outlier basata su un modello mistura di 2 gaussiane
#------------------------------------------------------------------------------
#         PARAMETRI 
#  y  = matrice ( o data.frame) -  Variabili dipendenti (con possibili errori)
#  x  = matrice ( o data.frame) - Variabili indipendenti (dati esatti. P.e. da archivio amministrativo)
#  model = Indica se i dati osservati hanno distribuzione log-normale (LN) o normale (N).
#  w = proporzione dei dati contaminati (peso a priori) 
#  max.iter = numero massimo di iterazioni per la convergenza EM
#  eps = soglia di accettazione
#  lambda = fattore di inflazione della varianza
#  graph = visualizzazione dei grafici durante l'elaborazione
#------------------------------------------------------------------------------
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
# Copio i dati di input su aree di appoggio
#------------------------------------------------------------------------------
  memo.y <- y <- as.matrix(y)         
  memo.x <- x
#------------------------------------------------------------------------------
#        CONTROLLI SUI PARAMETRI 
#  Eliminazione dei record contenenti missing per la stima dei parametri
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
#        CONTROLLI SUI PARAMETRI 
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
  

    

 
#---------------------       DEFINIZIONE VARIABILI       ---------------------
  conv <- FALSE
  continua <- TRUE
  lik<-NA
  lik0 <- 10
  oldlik <- 0
  iter <- 0
  sing <- FALSE
  
  if (ncol(x)+ ncol(y) < 3)   # INSERIRE BOXPLOT
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
#   CALCOLO DEL BIC per il modello normale da usare per i confronti
#------------------------------------------------------------------------------  
# N. parametri per il modello normale
    k1 <- ncol(x) * p + (p*(p+1))/2 # p=ncol(y) 
# N. parametri per il modello di contaminazione
    k2 <- k1 + 2 - w.fix - lambda.fix  
    
 if (n < k2)  {
      warning(paste("Input data are fewer than the number of model parameters\n" ))
 }
#------------------------------------------------------------------------------      
# Calcolo della verisimiglianza normale
#------------------------------------------------------------------------------  
    dati<-cbind(x,y)
    norm.mv<-function(u){dmvnorm(u[q+1:p], t(B0)%*%u[1:q], sigma0, log=TRUE)}
    lik.n <- sum(apply(dati,1,norm.mv))  
  BIC.n <- -2*lik.n + k1*log(n)  
  
#************************  INIZIO CICLO EM  ************************************

  while (iter < max.iter & continua == TRUE)
  {
     iter <- iter + 1
#     print(paste("E-step",iter))
#***********************    E - STEP         ************************************      
    tau1 <- post.prob(y, x, B, sigma, w1, lambda)   
    tau2 <- 1 - tau1 
#***********************    M - STEP         ************************************    
#    print(paste("M-step",iter))
#***********************    calcolo dei pesi  **********************************
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
    
  
#***********************    CONVERGENZA      **********************************


     
     s2 <- s1 / (1 + lambda) 
     sigma2 <- (1+lambda)* sigma 
     q1 <- matrix(tensorizza (dif, s1),n,1) 
     q2 <- matrix(tensorizza (dif, s2),n,1) 
     rm (s1,s2)
    
     
     q1 <- -0.5*q1
     q2 <- -0.5*q2
     
     ll <- w1 * exp(q1)  / sqrt(2*pi*det(sigma)) + (1-w1) * (exp(q2)) / sqrt(2*pi*det(sigma2))
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
 #************************  FINE CICLO EM  ************************************

    if (iter >= max.iter)  
     conv <- FALSE
#   CALCOLO DEI VALORI PREVISTI
   
      yprev <-  pred.y(y=memo.y, x=memo.x, B, sigma, 
                     lambda, w=1-w1, model = model, t.outl=t.outl)    
 
   ###############  calcolo di BIC e AIC per i due modelli #############
  
   
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



