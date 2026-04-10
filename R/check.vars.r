
#' Validate and Prepare Input Variables for Statistical Modeling
#'
#' This function performs input validation and data preparation for response variables (y)
#' and covariates (x). It checks data types, handles special cases for Log-Normal models,
#' and applies necessary transformations. Used as a preprocessing step in prediction and
#' estimation functions.
#'
#' @param y Response variable matrix (may contain NA values, numeric only)
#' @param x Covariate matrix (optional, must be complete - no NAs allowed)
#' @param model Character string specifying model type: "N" for Normal (default) or "LN" for Log-Normal
#' @param parent Character string identifying the calling function (for context-specific warnings)
#'
#' @return A list containing:
#'   \item{ret}{Return code: 0 = success, -9 = fatal error, >0 = warning count}
#'   \item{msg.err}{Character vector of error/warning messages (NULL if no issues)}
#'   \item{y}{Processed response matrix (log-transformed if model="LN")}
#'   \item{x}{Processed design matrix with intercept column prepended (log-transformed if model="LN")}
#'
#' @details
#' For Log-Normal models:
#'   - Zero values are replaced with half the minimum positive value
#'   - Log transformation is applied to both y and x
#'   - Negative values trigger a fatal error
#'
#' Fatal errors (ret=-9) occur when:
#'   - Variables are non-numeric
#'   - Data format is invalid
#'   - Negative values exist in LN model
#'   - NaN or Inf values are present
#'   - Covariates contain missing values
#'   - Dimension mismatch between y and x
#'
check.vars <- function(y, x, model="N", parent="pred.y") {
    
    # Note: Potential enhancement - automatically detect calling function name
    # using parent.env(env) instead of requiring 'parent' parameter
    
    # === INITIALIZE OUTPUT STRUCTURE ===
    # Convert y to matrix format for consistent handling
    y <- as.matrix(y)
    
    # Initialize return list with default values
    lista <- list(ret=0,        # Return code: 0=OK, -9=error, >0=warnings
                  msg.err=NULL,  # Error/warning messages
                  y=NULL,        # Processed y values
                  x=NULL)        # Processed x values (design matrix)
    
    # === VALIDATE RESPONSE VARIABLE (y) ===
    
    # Check 1: Ensure y contains only numeric values
    if (!is.numeric(y)) {
        lista$msg.err <- "Variables must be numeric "
        lista$ret <- -9
    }
    # Check 2: Ensure y is in acceptable format (matrix, data.frame, or numeric vector)
    else if (!inherits(y, c("data.frame", "matrix", "numeric", "integer"))) {
        lista$msg.err <- "Variables must be supplied as a matrix or as a dataframe "
        lista$ret <- -9
    }
    # Check 3: For Log-Normal model, negative values are mathematically invalid (can't take log)
    else if (model == 'LN' && sum(y < 0, na.rm=TRUE) > 0) {
        lista$msg.err <- "Negative values are not allowed if model is LN\n"
        lista$ret <- -9
    }
    # Check 4: NaN and Inf values cause computational issues and must be removed beforehand
    else if (sum(is.nan(y), na.rm=TRUE) > 0 || sum(is.infinite(y), na.rm=TRUE) > 0) {
        lista$msg.err <- " +/- Inf and NaN values are not allowed \n"
        lista$ret <- -9
    }
    
    # === VALIDATE COVARIATES (x) IF PROVIDED ===
    n <- nrow(y)
    
    if (!is.null(x)) {
        x <- as.matrix(x)
        
        # Check 1: Ensure covariates are numeric
        if (!is.numeric(x)) {
            lista$msg.err <- "Covariates must be numeric "
            lista$ret <- -9
        }   
        # Check 2: Ensure proper format
        else if (!inherits(x, c("data.frame", "matrix", "numeric", "integer"))) {
            lista$msg.err <- "Covariates must be supplied as a matrix or as a data frame "
            lista$ret <- -9
        }
        # Check 3: Covariates must be complete (no missing values allowed)
        # This is stricter than for y because incomplete covariates prevent model fitting
        else if (sum(is.na(x)) > 0) {
            lista$msg.err <- "Covariates can not have missing values "
            lista$ret <- -9
        } 
        # Check 4: Dimension compatibility - each y observation needs corresponding x values
        else if (nrow(x) != nrow(y)) {
            lista$msg.err <- "Variables y and x must have the same number of rows"
            lista$ret <- -9
        }
        # Check 5: For Log-Normal model, x also cannot have negative values
        else if (model == 'LN' && sum(x < 0) != 0) {
            lista$msg.err <- "Negative values are not allowed if model is LN\n"
            lista$ret <- -9
        }
    }
    
    # === EARLY RETURN ON FATAL ERROR ===
    # If any validation failed (ret=-9), return immediately with error message
    if (lista$ret == -9)
        return(lista)
    
    # === DATA PREPARATION BASED ON MODEL TYPE ===
    
    ####---------------------------------------------------------------------
    # LOG-NORMAL MODEL PREPROCESSING
    ####---------------------------------------------------------------------
    
    modificati <- 0  # Counter for modified zero values
    
    if (model == "LN") {  
        
        # === HANDLE ZERO VALUES IN RESPONSE VARIABLES ===
        # Log(0) is undefined, so zeros must be replaced before transformation
        for (i in 1:ncol(y)) {
            
            # Count zero values in this column
            n.zeri <- length(which(y[, i] == 0))
            
            if (n.zeri > 0) {
                # Identify zero positions
                ind0 <- which(y[, i] == 0)
                
                # Replace zeros with half the minimum positive value
                # This is a common approach that preserves relative scale
                new.val <- min(y[, i][-ind0]) / 2
                y[ind0, i] <- new.val
                
                modificati <- modificati + n.zeri
            } 
        }
        
        # === ISSUE WARNING ABOUT MODIFIED VALUES ===
        # Only warn if not called from ml.est (which may handle this differently)
        if (parent != "ml.est" && modificati > 0) {
            lista$ret <- lista$ret + 1  # Increment warning counter
            lista$msg.err <- rbind(lista$msg.err,                 
                                   paste(modificati, " response variable values (%", 
                                         round(modificati * 100 / n, 2),
                                         ") equal 0 are substituted by the half minimum of the corresponding variable\n", 
                                         sep=""))
        }
        
        # === APPLY LOG TRANSFORMATION TO RESPONSE ===
        y <- as.matrix(log(y))
        
        # === HANDLE ZERO VALUES IN COVARIATES ===
        modificati <- 0
        
        if (!is.null(x)) {
            x <- as.matrix(x)
            
            # Same zero-replacement procedure for covariates
            for (i in 1:ncol(x)) {
                n.zeri <- length(which(x[, i] == 0))
                
                if (n.zeri > 0) {
                    ind0 <- which(x[, i] == 0)
                    new.val <- min(x[, i][-ind0]) / 2
                    x[ind0, i] <- new.val
                    modificati <- modificati + n.zeri
                } 
            }
            
            # === WARN ABOUT MODIFIED COVARIATE VALUES ===
            if (parent != "ml.est" && modificati > 0) {                    
                lista$ret <- lista$ret + 1
                lista$msg.err <- rbind(lista$msg.err,                  
                                       paste(modificati, " covariate values (%", 
                                             round(modificati * 100 / length(x), 2),
                                             ") equal 0 are substituted by the half minimum of the corresponding variable\n", 
                                             sep=""))
            }
            
            # === APPLY LOG TRANSFORMATION TO COVARIATES ===
            x <- log(x)
        }
        else
            x <- NULL
        
    } 
    ####---------------------------------------------------------------------
    # NORMAL MODEL PREPROCESSING
    ####---------------------------------------------------------------------
    else {   
        # For Normal model, no transformation needed
        y <- as.matrix(y)
        x <- x
    }
    
    # === PREPARE FINAL OUTPUT ===
    
    # Store processed response variables
    lista$y <- y
    
    # Create design matrix by prepending intercept column
    # cbind(1, x) adds a column of ones for the intercept term in regression
    lista$x <- as.matrix(cbind(rep(1, n), x))
    
    # Return list with validated and processed data
    return(lista)
}