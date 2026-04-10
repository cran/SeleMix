#' Selective Editing for Survey Data
#'
#' This function implements a selective editing procedure to identify observations
#' that have the largest impact on estimated totals. It prioritizes records for 
#' manual review/editing based on their contribution to cumulative estimation error.
#'
#' @param y Response variable matrix (original/reported values, may contain NA)
#' @param ypred Predicted/expected values matrix from a statistical model
#' @param wgt Vector of weights for each observation (default: equal weights of 1)
#' @param tot Vector of target totals for each variable (default: weighted sum of predictions)
#' @param t.sel Threshold(s) for selective editing (default: 0.01 = 1% error tolerance)
#'              Can be a single value (applied to all variables) or a vector (one per variable)
#'
#' @return A matrix with class "seled" containing:
#'   - Original y values
#'   - Predicted values (ypred)
#'   - Weights
#'   - Individual scores (absolute relative errors)
#'   - Global score (maximum score across variables for each observation)
#'   - Cumulative residual errors
#'   - Selection flags per variable
#'   - Overall rank
#'   - Final selection indicator (1 = needs editing, 0 = no action needed)
#'
sel.edit <- function(y, ypred, wgt=rep(1, nrow(as.matrix(y))), tot=colSums(ypred * wgt), t.sel=0.01) {
    
    # === PRESERVE ORIGINAL DATA ===
    # Store original y values (including NAs) for final output
    memo.y <- y
    
    # === DATA PREPARATION ===
    y <- as.matrix(y)
    id <- 1:nrow(y)
    ypred <- as.matrix(ypred)
    n <- nrow(y)
    p <- ncol(y)
    
    # === INPUT VALIDATION ===
    # Check that y and ypred have compatible dimensions
    if (ncol(ypred) > p) 
        stop("Response variables (y) and their predicted values (ypred) have a different number of columns")
    
    # Check threshold vector length
    if (length(t.sel) > p) 
        stop("Length of threshold vector is greater of the variables number") 
    
    # Expand threshold to vector if single value provided
    if (length(t.sel) == 1) 
        t.sel <- rep(t.sel, p)
    
    # === SET COLUMN NAMES ===
    # Create standardized names if not provided
    if (is.null(colnames(y)))
        colnames(y) <- paste("y", 1:p, sep="")
    if (is.null(colnames(ypred)))
        colnames(ypred) <- paste(colnames(y), ".p", sep="")    
    
    # === HANDLE MISSING VALUES ===
    # Replace missing values in y with model predictions
    # This allows score calculation for all observations
    ind <- which(is.na(y))
    y[ind] <- ypred[ind]
    
    # === CALCULATE EDITING SCORES ===
    # Compute weighted errors relative to estimated totals
    # This measures each observation's impact on total estimation
    errore <- (wgt * (y - ypred)) / matrix(rep(tot, n), n, p, byrow=TRUE)
    
    # Score = absolute relative error (impact on total as proportion)
    score <- abs(errore) 
    
    # === PRIORITIZATION ===
    # Global score = maximum error across all variables for each observation
    # This identifies observations with the largest impact on ANY variable
    global.score <- sapply(1:n, function(i) { max(score[i, ])})
    
    # Sort observations by increasing global score
    # (least impactful first, most impactful last)
    ind_ord <- order(global.score)  
    
    # === CREATE SORTED DATA MATRIX ===
    # Arrange all data by increasing impact score
    dati <- cbind(id[ind_ord], 
                  y[ind_ord, , drop=FALSE], 
                  ypred[ind_ord, , drop=FALSE], 
                  wgt[ind_ord], 
                  score[ind_ord, , drop=FALSE], 
                  global.score[ind_ord])
    
    # === CUMULATIVE ERROR CALCULATION ===
    # Calculate cumulative sum of errors in sorted order
    # This shows how total error accumulates as we process observations
    cumulate <- matrix(apply(errore[ind_ord, , drop=FALSE], 2, cumsum), nrow=n, ncol=p)
    
    # === SELECTIVE EDITING DECISION ===
    # Initialize selection indicator (0 = no editing needed)
    sel <- rep(0, n)
    
    # Flag observations where cumulative error exceeds threshold
    # Each variable is checked against its threshold
    flag <- abs(cumulate) > matrix(t.sel, n, p, byrow=TRUE)
    mode(flag) <- "integer"
    
    # Find first observation where ANY variable exceeds threshold
    ind <- which(rowSums(flag) > 0)
    
    if (length(ind) > 0) {
        # Mark this observation and all higher-impact ones for editing
        # Logic: if cumulative error exceeds threshold at position k,
        # then all observations from k onward need review
        k <- min(ind)                                                 
        sel[k:n] <- 1 
    }
    
    # === COMPILE RESULTS ===
    # Add cumulative errors, flags, ranking, and selection indicator
    appo <- dati[, 1:ncol(dati), drop=FALSE]
    dati <- cbind(appo, 
                  cumulate,   # Cumulative residual errors
                  flag,       # Per-variable threshold flags
                  nrow(appo):1,  # Rank (highest impact = 1)
                  sel)        # Final selection (1 = edit, 0 = keep)
    
    # === RESTORE ORIGINAL ORDER ===
    # Sort back to original observation order for output
    dati.ord <- dati[order(dati[, 1]), 2:ncol(dati), drop=FALSE] 
    
    # === ASSIGN DESCRIPTIVE COLUMN NAMES ===
    colnames(dati.ord) <- c(colnames(y),                                
                            colnames(ypred),
                            "weights",                   
                            paste(colnames(y), ".score", sep=""),      # Individual scores
                            "global.score",                             # Max score
                            paste(colnames(y), ".reserr", sep=""),     # Cumulative errors
                            paste(colnames(y), ".sel", sep=""),        # Per-variable flags
                            "rank",                                     # Priority rank
                            "sel")                                      # Overall selection
    
    # === RESTORE ORIGINAL Y VALUES ===
    # Replace imputed values with original NAs for transparency
    dati.ord[, colnames(y)] <- as.matrix(memo.y)
    
    # === SET CLASS AND RETURN ===
    # Add "seled" class for potential method dispatch
    class(dati.ord) <- c(class(dati.ord), "seled")  
    dati.ord  
}