#' Calculate Posterior Probability for Mixture Model Component
#'
#' Computes the posterior probability that each observation belongs to the first
#' component of a two-component Gaussian mixture model. Used in robust estimation
#' and outlier detection to identify regular observations vs. contaminated ones.
#'
#' @param y Response variable matrix (n x p)
#' @param x Design matrix including intercept (n x q)
#' @param B Coefficient matrix (q x p)
#' @param sigma Covariance matrix for component 1 (p x p)
#' @param w1 Prior probability of belonging to component 1 (regular observations)
#' @param lambda Scale parameter for component 2 variance inflation
#'              Component 2 has covariance (1+lambda)*sigma, making it more diffuse
#'
#' @return Vector of posterior probabilities (length n) that each observation
#'         belongs to component 1 (regular/non-outlier component).
#'         Values close to 1 indicate regular observations.
#'         Values close to 0 indicate potential outliers/contaminated observations.
#'
#' @details
#' The mixture model has two components:
#'   Component 1: Regular observations ~ N(x*B, sigma) with prior weight w1
#'   Component 2: Contaminated observations ~ N(x*B, (1+lambda)*sigma) with prior weight (1-w1)
#'
#' Uses numerically stable log-space calculations to avoid underflow/overflow.
#' The "max trick" (subtracting maximum log-likelihood) ensures numerical stability.
#'
post.prob <- function(y, x, B, sigma, w1, lambda) {
    
    # === EXTRACT DIMENSIONS ===
    n <- nrow(y)  # Number of observations
    
    # === INITIALIZE MAXIMUM LOG-LIKELIHOOD TRACKER ===
    # Used for numerical stability in exponential calculations
    # Starting with very large negative number ensures it will be replaced
    maxold <- matrix(-999999999, n, 1)
    
    # === CALCULATE RESIDUALS ===
    # Difference between observed and predicted values
    # dif[i,] = y[i,] - E[y[i,]] under the model
    dif <- y - x %*% B
    
    # === COMPUTE LOG PRIOR PROBABILITIES ===
    # Work in log space to avoid numerical underflow
    wp1 <- matrix(log(w1), n, 1)      # log P(component 1) - regular observations
    wp2 <- matrix(log(1 - w1), n, 1)  # log P(component 2) - contaminated observations
    
    # === CALCULATE PRECISION MATRICES ===
    # Component 1: Uses original covariance sigma
    s1 <- solve(sigma)  # Precision matrix = inverse covariance
    
    # Component 2: Uses inflated covariance (1+lambda)*sigma
    # So precision is sigma^(-1) / (1+lambda)
    s2 <- s1 / (lambda + 1)
    
    # === COMPUTE MAHALANOBIS DISTANCES ===
    # Quadratic forms measuring how far each observation is from its expected value
    # q1[i] = t(dif[i,]) %*% s1 %*% dif[i,] - distance under component 1
    # q2[i] = t(dif[i,]) %*% s2 %*% dif[i,] - distance under component 2
    q1 <- matrix(tensorizza(dif, s1), n, 1)
    q2 <- matrix(tensorizza(dif, s2), n, 1)
    
    # === MEMORY MANAGEMENT ===
    # Remove large matrices no longer needed
    rm(s1, s2)
    gc()  # Force garbage collection to free memory
    
    # === COMPUTE LOG-LIKELIHOODS ===
    # Log-likelihood for each component (unnormalized, without determinant terms yet)
    # Formula: -0.5 * Mahalanobis_distance + log(prior_weight)
    lntau1 <- -0.5 * q1 + wp1  # Log-likelihood under component 1
    lntau2 <- -0.5 * q2 + wp2  # Log-likelihood under component 2
    
    # === NUMERICAL STABILITY: MAX TRICK ===
    # Find maximum log-likelihood across both components for each observation
    # This will be subtracted to prevent overflow when exponentiating
    maxln <- pmax(maxold, lntau1, lntau2)
    maxold <- maxln
    
    # Subtract maximum for numerical stability
    # exp(lntau - max) keeps values in reasonable range
    appo1 <- lntau1 - maxold
    appo2 <- lntau2 - maxold
    
    # === CALCULATE COVARIANCE DETERMINANTS ===
    # Component 2 has inflated covariance
    sigma2 <- (1 + lambda) * sigma
    
    # Square root of determinants appear in multivariate normal density
    # |sigma|^(1/2) in the normalizing constant
    dd1 <- sqrt(det(as.matrix(sigma)))   # sqrt(|sigma|) for component 1
    dd2 <- sqrt(det(as.matrix(sigma2)))  # sqrt(|(1+lambda)*sigma|) for component 2
    
    # === COMPUTE UNNORMALIZED POSTERIOR PROBABILITIES ===
    # Numerator of Bayes' theorem for each component
    # = exp(log-likelihood - max) / sqrt(det(Sigma))
    # The division by sqrt(det) completes the multivariate normal density
    numtau1 <- exp(appo1) / dd1  # Unnormalized posterior for component 1
    numtau2 <- exp(appo2) / dd2  # Unnormalized posterior for component 2
    
    # === NORMALIZE TO GET POSTERIOR PROBABILITIES ===
    # Denominator = sum of unnormalized posteriors (total probability)
    dentau <- (numtau1 + numtau2)
    
    # Posterior probability of component 1 (regular observations)
    # P(component 1 | data) = numtau1 / (numtau1 + numtau2)
    tau1 <- numtau1 / dentau
    
    # Return posterior probabilities
    # Values close to 1: likely regular observation
    # Values close to 0: likely outlier/contaminated
    tau1
}