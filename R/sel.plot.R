#' Selective Editing Diagnostic Plot with Interactive Identification
#'
#' Creates a scatterplot or plot matrix highlighting outliers and influential observations,
#' with optional interactive identification of specific points. Useful for visual diagnosis
#' in selective editing workflows.
#'
#' @param data Data frame or matrix containing the variables to plot
#' @param vars Vector of column indices or names to plot (default: columns 1 and 2)
#' @param outl Binary vector indicating outlier status (1 = outlier, 0 = regular)
#' @param sel Binary vector indicating influential/selected status (1 = influential, 0 = not)
#' @param log Logical; if TRUE, applies log transformation to data (default: TRUE)
#' @param n.identify Number of points to interactively identify (default: 0 = no identification)
#' @param file Character string; if provided, saves plot to JPEG file (only when n.identify=0)
#' @param title Character string for plot title (optional)
#'
#' @details
#' Color coding:
#'   - Light grey: Regular observations
#'   - Blue: Only outliers (outl=1, sel=0)
#'   - Red: Only influential (outl=0, sel=1)
#'   - Cyan: Both outlier and influential (outl=1, sel=1)
#'
#' When n.identify > 0, user can click on points to identify them. Identified points
#' are labeled and their full records are printed or saved to file.
#'
#' @return Invisible vector of identified point indices (if n.identify > 0)
#'
sel.plot <- function(data, vars=1:2, outl = rep(0, nrow(data)), sel = rep(0, nrow(data)),
                     log = TRUE, n.identify=0, file=NULL, title=NULL) {
    
    #-------------------------------------------------------------------------------------------  
    #' Internal Function: Interactive Point Identification
    #'
    #' Allows user to click on plot to identify specific observations. Selected points
    #' are labeled with their row numbers and their full data records are output.
    #'
    #' @param x Coordinates object (from xy.coords) or x-y data
    #' @param dati Original data frame to extract full records from
    #' @param n Maximum number of points to identify
    #' @param filename Optional file path to save identified records as CSV
    #'
    #' @return Integer vector of identified point indices
    #'
    identifyPnt <- function(x, dati, n, filename) {
        
        # === EXTRACT COORDINATES ===
        # Convert input to standard x-y coordinate format
        xy <- xy.coords(x)
        x <- xy$x
        y <- xy$y
        
        # === INITIALIZE TRACKING VARIABLES ===
        vedi <- rep(FALSE, length(x))  # Track which points have been identified
        res <- integer(0)              # Store indices of identified points
        
        # === INTERACTIVE IDENTIFICATION LOOP ===
        # Continue until n points identified or user stops
        while(sum(vedi) < n) {
            # Let user click on one point (only show unidentified points)
            ans <- identify(x[!vedi], y[!vedi], n=1, plot=FALSE)
            
            # Break if user cancels (e.g., presses ESC or right-clicks)
            if(!length(ans)) break
            
            # Convert index back to original data position
            ans <- which(!vedi)[ans]
            
            # Add to results
            res <- c(res, ans)
            
            # === ANNOTATE IDENTIFIED POINT ===
            # Add text label with row number
            text(x[ans], y[ans], labels=ans, offset=1, cex=.8)
            # Highlight with special symbol
            points(x[ans], y[ans], pch=23)
            
            # Mark as identified
            vedi[ans] <- TRUE
        }
        
        # === OUTPUT IDENTIFIED RECORDS ===
        if (!is.null(filename))
            # Save to CSV file if filename provided
            write.csv(dati[res,], filename)
        else    
            # Otherwise print to console
            print(dati[res,])                     
        
        res  # Return vector of identified indices
    }
    #-------------------------------------------------------------------------------------------
    
    # === INPUT VALIDATION ===
    # Ensure data is in proper format
    if (!inherits(data, c("data.frame", "matrix"))) {
        warning("Variables must be supplied as a matrix or as a data frame")
        return()
    }
    data <- as.data.frame(data)
    
    # === VALIDATE IDENTIFICATION REQUEST ===
    # Interactive identification requires at least 2 variables for a scatterplot
    if (ncol(data) < 2 && n.identify > 0) {
        warning("Identify not allowed for plot of a single variable")
        n.identify <- 0
    }
    
    # === SET UP FILE OUTPUT ===
    # Open JPEG device if saving to file (only when not using interactive identification)
    if (!is.null(file) && n.identify == 0)
        jpeg(file, width = 1024, height = 1024, quality = 100, pointsize = 18)
    
    # === APPLY LOG TRANSFORMATION ===
    # Transform data if requested to handle wide-ranging values
    if (log == TRUE) {
        ldata <- data[, vars, drop = FALSE]
        # Replace zeros with small value to avoid log(0) = -Inf
        data[data==0] <- 1e-07
        ldata <- log(ldata[, vars, drop = FALSE])
    }
    else 
        ldata <- data[, vars, drop = FALSE]
    
    # === CREATE BASE PLOT ===
    # Plot all observations in light grey
    # Subtitle shows counts of total obs, outliers, and influential errors
    plot(ldata, pch = 21, col = "lightgrey",
         sub=paste("obs", nrow(data), 
                   "outliers", sum(outl), 
                   "influential errors", sum(sel)))
    
    # === OVERLAY PROBLEMATIC OBSERVATIONS ===
    # Category 1: Only outliers (blue circles)
    # Statistical outliers but low impact on estimates
    points(ldata[outl == 1 & sel == 0, vars, drop = FALSE], 
           pch = 21, cex = 0.9, 
           col = "blue", bg = "blue")
    
    # Category 2: Both outlier AND influential (cyan filled circles) - most critical
    # High priority for manual review
    points(ldata[sel == 1 & outl == 1, vars, drop = FALSE], 
           pch = 20, cex = 0.9, 
           col = "cyan", bg = "cyan")
    
    # Category 3: Only influential (red filled circles)
    # High impact on estimates but not statistical outliers
    points(ldata[outl == 0 & sel == 1, vars, drop = FALSE], 
           pch = 20, cex = 0.9, 
           col = "red", bg = "red")
    
    # === ADD TITLE IF PROVIDED ===
    if (!is.null(title))
        title(main=title)
    
    # === INTERACTIVE IDENTIFICATION ===
    # If requested, allow user to click and identify specific points
    if (n.identify > 0)
        identifyPnt(ldata[, vars], n=n.identify, dati=data, filename=file)
    
    # === CLOSE GRAPHICS DEVICE ===
    # Close JPEG file if it was opened
    if (!is.null(file) && n.identify == 0) 
        dev.off()
}