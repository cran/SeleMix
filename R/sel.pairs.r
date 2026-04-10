#' Scatterplot Matrix with Outlier and Influential Unit Highlighting
#'
#' Creates a matrix of scatterplots (pairs plot) with special highlighting for
#' outliers and influential observations. The diagonal shows boxplots instead of
#' scatterplots. Used for diagnostic visualization in selective editing workflows.
#'
#' @param x Data matrix or data frame containing the variables to plot
#' @param outl Binary vector indicating outlier status (1 = outlier, 0 = regular)
#'             Default: all zeros (no outliers)
#' @param sel Binary vector indicating influential/selected status (1 = influential, 0 = not)
#'            Default: all zeros (no influential observations)
#' @param labs Character vector of variable labels (default: column names or "x1", "x2", ...)
#' @param log Logical; if TRUE, applies log transformation to data (default: TRUE)
#'            Zero values are replaced with 1e-07 before transformation
#' @param legend Logical; if TRUE, displays legend in top-left plot (default: TRUE)
#' @param title Character string for overall plot title (default: preset message)
#'
#' @details
#' The function distinguishes three types of problematic observations:
#'   - Only Outlier (outl=1, sel=0): Blue circles - statistical outliers but low impact
#'   - Only Influential (outl=0, sel=1): Red triangles - high impact on estimates
#'   - Both (outl=1, sel=1): Cyan squares - outliers with high impact (most critical)
#'
#' Regular observations are shown in light grey.
#'
#' @return NULL (produces plot as side effect)
#'
#' @note Diagonal elements show boxplots with problematic observations overlaid.
#'       Off-diagonal elements show bivariate scatterplots.
#'
sel.pairs <- function(x, outl = rep(0, nrow(x)), sel = rep(0, nrow(x)), labs = NULL,
                      log = TRUE, legend = TRUE, title = NULL) {
    
    # === INPUT VALIDATION ===
    # Ensure data is in proper matrix or data frame format
    if (!inherits(x, c("matrix", "data.frame")))
        stop("data must be supplied as matrix or data frame ")
    
    # === DEFINE LEGEND TEXT ===
    # Three categories of problematic observations based on combinations of outl and sel flags
    testo_legenda <- c("Only Outlier           ",      # outl=1, sel=0: Blue circles
                       "Only Influential       ",       # outl=0, sel=1: Red triangles
                       "Outlier and Influential")       # outl=1, sel=1: Cyan squares (highest priority)
    
    # === SET DEFAULT TITLE ===
    # Use custom title if provided, otherwise use default
    if (is.null(title))
        title <- "Selective Editing - outliers and influential errors"
    
    # === EXTRACT DIMENSIONS ===
    nvar1 <- ncol(x)  # Number of variables determines grid size (nvar1 x nvar1)
    
    # === SET UP VARIABLE LABELS ===
    # Priority: 1) User-provided labs, 2) Column names from x, 3) Generic "x1", "x2", etc.
    if (is.null(labs)) {
        if (length(colnames(x) > 0))
            labs <- colnames(x)
        else 
            labs <- paste("x", 1:nvar1, sep = "")
    }
    
    # === APPLY LOG TRANSFORMATION ===
    # Log transformation helps visualize data with wide ranges or skewed distributions
    if (log == TRUE) {
        # Replace zeros with small positive value to avoid log(0) = -Inf
        x[x == 0] <- 1e-07
        x <- log(x)
    }
    
    # === SET UP PLOT LAYOUT ===
    # Create nvar1 x nvar1 grid of plots (scatterplot matrix)
    # mfrow: arrange plots in a matrix layout
    # mar: inner margins (bottom, left, top, right) in lines of text
    # oma: outer margins for overall title (bottom, left, top, right)
    par(mfrow = c(nvar1, nvar1), mar = c(3, 2, 2, 3), oma = c(0, 0, 3, 0))    
    
    # === CREATE SCATTERPLOT MATRIX ===
    # Double loop creates nvar1 x nvar1 grid of plots
    # Row i, Column j shows relationship between variable i (x-axis) and variable j (y-axis)
    
    for (j in 1:nvar1) {  # Loop over columns (variables on y-axis)
        for (i in 1:nvar1) {  # Loop over rows (variables on x-axis)
            
            # === DIAGONAL PLOTS: UNIVARIATE BOXPLOTS ===
            # When i==j, show distribution of single variable instead of self-scatter
            if (i == j) {
                
                # Create horizontal boxplot for variable i
                boxplot(x[, i], 
                        main = NULL,           # No individual plot title
                        col = "peachpuff",     # Fill color for box
                        horizontal = TRUE,     # Horizontal orientation
                        names = labs[i],       # Variable label
                        show.names = TRUE)
                
                # === OVERLAY PROBLEMATIC OBSERVATIONS ON BOXPLOT ===
                
                # Category 1: Only outliers (outl=1, sel=0) - Blue circles
                # Statistical outliers but low impact on estimates
                aa <- x[outl == 1 & sel == 0, i]
                points(aa, rep(1, length(aa)),    # y=1 positions points on boxplot
                       pch = 21,                   # Circle symbol
                       cex = 1.2,                  # Point size
                       col = "blue4",              # Border color
                       bg = "blue")                # Fill color
                
                # Category 2: Both outlier AND influential (outl=1, sel=1) - Cyan squares
                # HIGHEST PRIORITY: Statistical outliers with high impact
                aa <- x[outl == 1 & sel == 1, i]
                points(aa, rep(1, length(aa)), 
                       pch = 22,                   # Square symbol
                       cex = 1.2,
                       col = "cyan4", 
                       bg = "cyan")
                
                # Category 3: Only influential (outl=0, sel=1) - Red triangles
                # High impact on estimates but not statistical outliers
                aa <- x[sel == 1 & outl == 0, i]
                points(aa, rep(1, length(aa)), 
                       pch = 24,                   # Triangle symbol
                       cex = 1.2,
                       col = "red4", 
                       bg = "red") 
                
                # === ADJUST AXIS PARAMETERS ===
                par(mgp = c(2, 1, 0))  # Margin line for axis title, labels, and line
                # Commented out: Optional y-axis label customization
                # title(ylab = labs[i], col.lab = "red3", cex=1.2, font.lab=2, outer=FALSE)
                par(mgp = c(3, 1, 0))  # Reset to default
                
                # === ADD LEGEND TO TOP-LEFT PLOT ===
                # Show legend only once (first diagonal element, i.e., top-left) if requested
                if (i == 1 & j == 1 & legend == TRUE) {
                    legend("topleft", 
                           legend = testo_legenda,           # Text labels
                           pch = c(21, 24, 22),              # Symbol types (circle, triangle, square)
                           col = c("blue4", "red4", "cyan4"), # Border colors
                           pt.bg = c("blue", "red", "cyan"), # Fill colors
                           xjust = 0,                        # Horizontal justification
                           yjust = 1,                        # Vertical justification  
                           cex = 2/3)                        # Legend text size (2/3 of default)
                }
            }
            
            # === OFF-DIAGONAL PLOTS: BIVARIATE SCATTERPLOTS ===
            # When i!=j, show relationship between two different variables
            else {
                
                # Create base scatterplot with all points in light grey
                plot(x[, i], x[, j],      # x[,i] on x-axis, x[,j] on y-axis
                     xlab = labs[i],       # X-axis label
                     ylab = labs[j],       # Y-axis label
                     pch = 21,             # Circle symbol
                     col = "lightgrey")    # Color for regular observations
                
                # Extract the two variables as a matrix for subsetting
                appo <- cbind(x[, i], x[, j])
                
                # === OVERLAY PROBLEMATIC OBSERVATIONS ===
                
                # Category 1: Only outliers (outl=1, sel=0) - Blue circles
                points(appo[outl == 1 & sel == 0, , drop = FALSE], 
                       pch = 21, 
                       cex = 1.2, 
                       col = "blue4", 
                       bg = "blue")
                
                # Category 2: Both outlier AND influential (outl=1, sel=1) - Cyan squares
                # HIGHEST PRIORITY for review
                points(appo[sel == 1 & outl == 1, , drop = FALSE], 
                       pch = 22, 
                       cex = 1.2, 
                       col = "cyan4", 
                       bg = "cyan")
                
                # Category 3: Only influential (outl=0, sel=1) - Red triangles
                points(appo[outl == 0 & sel == 1, , drop = FALSE], 
                       pch = 24, 
                       cex = 1.2, 
                       col = "red4", 
                       bg = "red")
            }
        }
    }
    
    # === ADD OVERALL TITLE ===
    # Place title in outer margin at top (side=3)
    # Spans across all plots in the matrix
    mtext(title, side = 3, outer = TRUE)    
}