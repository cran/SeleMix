 sel.plot <- function(data,vars=1:2, outl = rep(0, nrow(data)), sel = rep(0, nrow(data)),
   log = TRUE,n.identify=0, file=NULL, title=NULL)
 {
 #-------------------------------------------------------------------------------------------  
   identifyPnt <- function(x,  dati, n,  filename)
   {
     xy <- xy.coords(x)
     x <- xy$x
     y <- xy$y
     vedi <- rep(FALSE, length(x))
     res <- integer(0)
     
     while(sum(vedi) < n) {
       ans <- identify(x[!vedi], y[!vedi], n=1, plot=FALSE)
       if(!length(ans)) break
       
       ans <- which(!vedi)[ans]
       res <- c(res, ans)
       text (x[ans], y[ans], labels=ans, offset=1, cex=.8)
       points(x[ans], y[ans], pch=23)
      
       vedi[ans] <- TRUE
       
     }
     if (!is.null(filename))
      write.csv(dati[res,], filename)
     else    
      print(dati[res,])                     
     res
   }
#-------------------------------------------------------------------------------------------
   if (!inherits(data,c("data.frame", "matrix"))  )  {
      warning( "Variables must be supplied as a matrix or as a data frame ")
      return
   }
   data <- as.data.frame(data)
   
   
   
    if (ncol(data) < 2 && n.identify > 0 )  {
      warning( "Identify not allowed for plot of a single variable")
      n.identify<-0
      
   }
   if (!is.null(file) && n.identify == 0)
       jpeg(file, width = 1024, height = 1024,quality = 100,  pointsize = 18)
   if (log == TRUE) {
        ldata <- data[, vars, drop = FALSE]
        data[data==0] <- 1e-07
        ldata <- log(ldata[, vars, drop = FALSE])
   }
   else ldata <- data[, vars, drop = FALSE]
   
     plot(ldata,   pch = 21, col = "lightgrey",
        sub=paste("obs",nrow(data), "outliers", sum(outl), "influential errors",sum(sel)))
     
     points(ldata[outl == 1 & sel == 0,vars , drop = FALSE], pch = 21,
                  cex = 0.9, col = "blue", bg = "blue")
     points(ldata[sel == 1 & outl == 1, vars, drop = FALSE], pch = 20,
                  cex = 0.9, col = "cyan", bg = "cyan")
     points(ldata[outl == 0 & sel == 1, vars, drop = FALSE], pch = 20,
                  cex = 0.9, col = "red",  bg = "red")
      
     if (!is.null(title))
       title (main=title)
     
     if (n.identify > 0)
        identifyPnt(ldata[,vars], n=n.identify, dati=data, filename=file)
     if (!is.null(file) && n.identify == 0) 
        dev.off()
 }
