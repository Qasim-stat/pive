# pive
The "pive" package provides estimates for causal models using many (weak) invalid instruments and heteroscedasticity. This package uses six instrumental variables methods to robustly estimate causal effects. Additionally, it conducts bootstrap sampling to derive estimates, standard errors, and confidence intervals.

## Installation

You can install the **pive** package from GitHub
```r
# Install devtools if not already installed
install.packages("devtools")

# Install pive from GitHub
devtools::install_github("Qasim-stat/pive")
# Load the package
library(pive)

# Example
set.seed(123)
n <- 1000
L <- 30
r <- ceiling(L * 0.3)
pii <- sqrt(8 / (n * L))
pi <- rep(pii, L)
beta <- 1
delta <- c(rep(1, r), rep(0, L - r))
Z <- MASS::mvrnorm(n, rep(0, L), diag(L))
sigma_e <- (1 + 0.3 * rowSums(Z[, (r + 1):L])) * rnorm(n)
u <- 0.2 * sigma_e + sqrt(1 - 0.2^2) * rnorm(n)
X <- Z %*% pi + u
Y <- X * beta + Z %*% delta + sigma_e

# Perform causal analysis without bootstrap
results <- pive(Y, X, Z, bootstrap = FALSE)
# Perform causal analysis bootstrap
results <- pive(Y, X, Z)
# Print the results
print(results)
