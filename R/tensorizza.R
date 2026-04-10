tensorizza <- function(dif, mat) {
    
    # === INPUT PREPARATION ===
    # Ensure dif is a matrix
    dif <- as.matrix(dif)
    
    # Extract dimensions
    p <- ncol(dif)  # Number of variables (dimension of each vector)
    n <- nrow(dif)  # Number of observations
    
    # === CREATE 3D ARRAYS FOR TENSOR OPERATIONS ===
    
    # Step 1: Replicate mat matrix n times (one for each observation)
    # Creates 3D array with dimensions (p, p, n)
    matrep <- array(mat, c(p, p, n))
    
    # Step 2: Replicate dif for element-wise multiplication
    # Creates 3D array with dimensions (n, p, p)
    # This will represent dif[i,j] in position [i, j, k] for all k
    difrep <- array(dif, c(n, p, p))
    
    # Step 3: Create transposed version of dif for multiplication
    # Also dimensions (n, p, p) but arranged differently
    dift <- array(dif, c(n, p, p))
    
    # === REARRANGE DIMENSIONS FOR PROPER MULTIPLICATION ===
    
    # Permute matrep from (p, p, n) to (n, p, p)
    # This aligns dimensions so mat[i] corresponds to observation i
    matrep <- aperm(matrep, c(3, 1, 2))
    
    # Permute dift to transpose the last two dimensions
    # Changes from (n, p, p) with [i, j, k] to [i, k, j]
    # This creates the transpose needed for t(dif[i,]) in the quadratic form
    dift <- aperm(dift, c(1, 3, 2))
    
    # === COMPUTE QUADRATIC FORMS VIA ELEMENT-WISE MULTIPLICATION ===
    
    # Multiply element-wise: difrep * matrep * dift
    # For each observation i, this computes all elements of:
    #   dif[i,] * mat * t(dif[i,])
    # Result has dimensions (n, p, p)
    res <- difrep * matrep * dift
    
    # === SUM ACROSS MATRIX DIMENSIONS ===
    
    # For each observation (dimension 1), sum across both matrix dimensions (2 and 3)
    # This collapses the (p x p) matrix for each observation into a single scalar
    # Equivalent to sum(dif[i,] * mat * t(dif[i,])) for each i
    resfin <- apply(res, 1, sum)
    
    # Return vector of n quadratic form values
    return(resfin)
}