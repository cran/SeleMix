  
  sel.edit  <- function (y, ypred, wgt=rep(1,nrow(as.matrix(y ))), tot=colSums(ypred * wgt), t.sel=0.01)  
  {
   memo.y <- y
   y=as.matrix(y)
   id=1:nrow(y)
   ypred=as.matrix(ypred)
   n <- nrow(y)
   p <- ncol(y)
   if (is.null(colnames(y)))
       colnames(y)<-paste("y",1:p,sep="")
   if (is.null(colnames(ypred)))
       colnames(ypred)<-paste(colnames(y),".p",sep="")    
 # Assumo che nel caso di y=NA il valore sia sostituito dal valore previsto dal modello 
 #  cioe' che la previsione sia esatta
   ind <- which(is.na(y))
   y[ind] <- ypred[ind]
   
 # CALCOLO DEL PUNTEGGIO: rapporto fra l'errore atteso e la stima di riferimento (TOTALE)
   errore <- (wgt * (y - ypred) )/ matrix(rep(tot,n), n, p, byrow=TRUE)
   score <- abs(errore) 

 # definizione di un vettore per l'ordinamento: massimo fra i valori assoluti dei punteggi delle variabili 
 # score_sas <- errore^2

    global.score <- sapply (1:n, function(i) { max(score[i, ])})

#  INDICE ORDINATO
   ind_ord <- order(global.score)  #ordinamento  crescente
 
# MATRICE DEI DATI ORDINATA IN BASE AL VALORE DELLO |SCORE| 

   dati <- cbind( id[ind_ord], y[ind_ord,,drop=FALSE], ypred[ind_ord,,drop=FALSE], wgt[ind_ord], 
                  score[ind_ord,,drop=FALSE], global.score[ind_ord] )
                 
# CALCOLO DELL'ERRORE CUMULATO  
  # cumulate<-matrix ( nrow=n, ncol=p)
  #  sapply(1:p, function(i) { cumulate[ind_ord,i]<<-cumsum (errore[ind_ord,i]) } )
    cumulate <- matrix(apply (errore[ind_ord,,drop=FALSE], 2, cumsum), nrow=n, ncol=p)
#   
   sel <- rep(0,n)
   flag <- NULL
   flag <-  matrix(sapply(1:ncol(cumulate), function(i) {
                                         fl<- rep(0,n)
                                         fl[(abs(cumulate[,i]) > t.sel)] <- 1 
                                         cbind(flag,fl)
                                         }
                            ), nrow=n, ncol=p)
  
   ind <- which(rowSums (flag) > 0)
   if (length(ind) > 0) {
      k <- min(ind)                                                 
      sel[k:n] <- 1 
   }
   
   sel.y <- sapply(1:ncol(flag), function(i) {  
                                         ind <- which(flag[,i] == 1)
                                         if (length(ind) > 0) {
                                            flag[(min(ind)):n,i] <<- 1
                                            }
                                         }
                            )
                    
   
   appo <- dati[,1:ncol(dati), drop=FALSE]
   dati <- cbind(appo, cumulate,  flag, nrow(appo):1, sel )
   dati.ord<- dati[order(dati[,1]),2:ncol(dati),drop=FALSE] 
   colnames(dati.ord) <- c(colnames(y),                                
                    colnames(ypred),
                    "weights",                   
                    paste(colnames(y),".score",sep=""),
                    "global.score",
                    paste(colnames(y),".reserr",sep=""),
                    paste(colnames(y),".sel",sep=""),
                    "rank",
                    "sel" ) 
   
   dati.ord[,colnames(y)]<-as.matrix(memo.y )
   class(dati.ord) <- c(class(dati.ord), "seled" )  
   dati.ord  

 } 

  
