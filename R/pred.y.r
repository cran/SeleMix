
#' Predict Y Values with Missing Data Handling and Outlier Detection
#'
#' This function performs prediction of response variables (y) with support for 
#' missing data patterns, outlier detection, and both Normal and Log-Normal models.
#' It uses a Bayesian approach with mixture modeling to handle contaminated observations.
#'
#' @param y Response variable matrix (may contain NA/Inf values)
#' @param x Predictor variable matrix (optional)
#' @param B Coefficient matrix for the linear model
#' @param sigma Covariance matrix of the response variables
#' @param lambda Precision parameter for the prior distribution
#' @param w Weight parameter for the mixture model (probability of being an outlier)
#' @param model Type of model: "N" for Normal (default) or "LN" for Log-Normal
#' @param t.outl Threshold for outlier classification (default 0.5)
#'
#' @return A data frame containing:
#'   - ypred: Predicted values for all y variables
#'   - tau: Posterior probability of being a regular observation (1 - outlier probability)
#'   - outlier: Binary indicator (1 if tau > t.outl, 0 otherwise)
#'   - pattern: Missing data pattern identifier
#'
pred.y <- function(y, x=NULL, B, sigma, lambda, w, model = "N", t.outl=0.5) {
    
    # === INPUT VALIDATION ===
    # Check and validate input variables (y, x, model type)
    vars <- check.vars(y, x, model, parent="pred.y")
    
    if (vars$ret == -9) {
        stop(vars$msg.err) 
    }
    if (vars$ret != 0) {
        warning(vars$msg.err) 
    } 
    
    # === DATA PREPARATION ===
    # Convert inputs to matrices and extract dimensions
    y <- as.matrix(vars$y)
    x <- as.matrix(vars$x) 
    nomicol.y <- colnames(y)  # Store original column names
    id <- 1:nrow(y)
    ncoly <- ncol(y)
    n <- nrow(y)
    ncolx <- ncol(x)
    
    # Create standardized variable names
    nomiy <- paste("y", 1:ncoly, sep="")
    nomix <- paste("x", 1:ncolx, sep="")
    
    # === MISSING DATA PATTERN DETECTION ===
    # Create pattern matrix: 1 for observed values, 0 for missing/infinite
    mat.pat <- (is.na(as.matrix(y))) | (is.infinite(y))   
    mat.pat <- 1 - mat.pat  # Flip: 1 = observed, 0 = missing
    
    colnames(x) <- nomix
    colnames(y) <- nomiy
    colnames(mat.pat) <- paste(nomiy, ".pat", sep="")
    
    # Convert pattern matrix to string identifier for each row
    # e.g., "111" = all observed, "101" = middle variable missing
    pattern <- apply(mat.pat, 1, function(k) { 
        appo <- ""
        for (i in 1:length(k)) { 
            appo <- paste(appo, k[i], sep="")
        }
        appo  
    })
    
    # === GROUP BY PATTERN ===
    # Combine all data and split by missing data pattern
    # This allows processing observations with the same pattern together
    appo.df <- data.frame(x, y, 1:n, mat.pat, w, pattern)
    colnames(appo.df) <- c(nomix, nomiy, "id", colnames(mat.pat), "w", "pattern")
    appo <- split(appo.df, pattern)
    
    # === PROCESS EACH PATTERN GROUP ===
    stime <- NULL  # Will accumulate results
    tau2 <- NULL
    
    for (i in 1:length(appo)) {  
        pattern.corr <- appo[[i]]    
        y.corr <- as.matrix(pattern.corr[, nomiy])
        x.corr <- as.matrix(pattern.corr[, nomix])
        
        # Identify which variables are observed vs missing in this pattern
        ind.obs <- as.logical(pattern.corr[1, colnames(mat.pat)])
        ind.mis <- !ind.obs
        
        # === COMPUTE CONDITIONAL DISTRIBUTIONS ===
        # Calculate mean and covariance for observed and missing components
        mu <- as.matrix(x.corr %*% B)
        mu.o <- mu[, ind.obs, drop=FALSE]   # Mean for observed variables
        mu.m <- mu[, ind.mis, drop=FALSE]   # Mean for missing variables
        
        # Partition covariance matrix
        sigma.oo <- sigma[ind.obs, ind.obs, drop=FALSE]  # Observed-observed block
        sigma.mm <- sigma[ind.mis, ind.mis, drop=FALSE]  # Missing-missing block
        sigma.om <- sigma[ind.obs, ind.mis, drop=FALSE]  # Observed-missing block
        
        # Posterior covariance adjustments using lambda (shrinkage parameter)
        sigma.bar.oo <- (lambda * sigma.oo) / (lambda + 1)    
        sigma.bar.om <- (lambda * sigma.om) / (lambda + 1)  # Not used directly
        
        # === HANDLE COMPLETELY MISSING CASE ===
        if (sum(ind.obs) == 0) {
            # No observed variables: use prior mean and variance
            mu.bar1 <- mu.bar2 <- mu.m
            diag.sigma.bar1 <- diag.sigma.bar2 <- diag(sigma.mm)
            tau1 <- 1 - pattern.corr[, "w"]  # tau1 is arbitrary here
        } else {
            # === CONDITIONAL PREDICTION FOR MISSING VARIABLES ===
            # Residual covariance of missing given observed
            sigma.res.mo <- sigma.mm - (t(sigma.om) %*% solve(sigma.oo) %*% sigma.om)
            
            # Adjusted covariance for Bayesian shrinkage
            sigma.bar.mm <- sigma.mm - ((t(sigma.om) %*% solve(sigma.oo) %*% sigma.om) / (lambda + 1))   
            
            # === OUTLIER DETECTION ===
            # Compute posterior probability of being a regular (non-outlier) observation
            tau1 <- post.prob(y.corr[, ind.obs, drop=FALSE], x.corr, B[, ind.obs, drop=FALSE],
                              sigma.oo, 1-w, lambda)                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                            
            
            # Shrunk observed values (toward prior mean)
            mu.tilde.obs <- (y.corr[, ind.obs, drop=FALSE] + (lambda * mu.o)) / (lambda + 1) 
            
            # === TWO PREDICTION MODES ===
            # Mode 1: Assume observation is regular (use actual observed values)
            mu.bar1 <- cbind(y.corr[, ind.obs, drop=FALSE], 
                             mu.m + (y.corr[, ind.obs, drop=FALSE] - mu.o) %*% solve(sigma.oo) %*% sigma.om)
            
            # Mode 2: Assume observation is contaminated (use shrunk values)
            mu.bar2 <- cbind(mu.tilde.obs, 
                             mu.m + (mu.tilde.obs - mu.o) %*% solve(sigma.oo) %*% sigma.om)   
            
            # Variance components for each mode
            diag.sigma.bar1 <- c(rep(0, ncol(mu.o)), as.numeric(diag(sigma.res.mo)))
            diag.sigma.bar2 <- c(diag(sigma.bar.oo), as.numeric(diag(sigma.bar.mm)))  
        }
        
        # === MIXTURE PREDICTION ===
        # Combine predictions from both modes weighted by posterior probabilities
        if (model == "LN") {
            # Log-Normal model: transform predictions back to original scale
            mu.bar <- as.vector(tau1) * (t(exp(t(mu.bar1) + (0.5 * diag.sigma.bar1)))) + 
                as.vector(1 - tau1) * (t(exp(t(mu.bar2) + (0.5 * diag.sigma.bar2))))
        } else {
            # Normal model: weighted average of predictions
            mu.bar <- as.vector(tau1) * mu.bar1 + as.vector(1 - tau1) * mu.bar2    
        }
        
        # === FORMAT RESULTS ===
        colnames(mu.bar) <- c(nomiy[ind.obs], nomiy[ind.mis])     
        stima.corr <- matrix(NA, nrow(pattern.corr), ncoly+1)
        
        # Combine predictions with outlier probabilities
        stima.corr <- cbind(mu.bar[, nomiy, drop=FALSE], 1-tau1)
        
        # Accumulate results with pattern and ID information
        stime <- rbind(stime, cbind(stima.corr, pattern.corr[, c("pattern", "id")]))
    }
    
    # === PREPARE OUTPUT ===
    # Sort results back to original order
    stime <- stime[order(stime[, "id"]), ]
    
    # Extract predictions
    ypred <- as.matrix(stime[, 1:ncoly, drop=FALSE])
    
    # Restore original column names or create new ones
    if (is.null(nomicol.y))
        colnames(ypred) <- paste("ypred", 1:ncoly, sep="") 
    else       
        colnames(ypred) <- paste(nomicol.y, "p", sep=".")
    
    # Extract outlier metrics
    tau <- stime[, ncoly+1]  # Posterior probability of being regular
    outlier <- (stime[, ncoly+1] > t.outl) + 0  # Binary outlier flag
    pattern <- stime[, "pattern"]
    
    # Return comprehensive results as data frame
    data.frame(cbind(ypred, tau, outlier), pattern)             
}
  
