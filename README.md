# pive
The "pive" package provides estimates for causal models using many (weak) invalid instruments and heteroscedasticity. This package uses six instrumental variables (IV) methods to robustly estimate causal effects. Additionally, it conducts bootstrap sampling to derive estimates, standard errors, and confidence intervals.

## Installation

You can install the **pive** package from GitHub

# Install devtools if not already installed
install.packages("devtools")

# Install pive from GitHub
devtools::install_github("Qasim-stat/pive")
# Load the package
library(pive)

## Example
```r
set.seed(123)
library(MASS)
n <- 1000
L <- 30
r <- ceiling(L * 0.3)
pii <- sqrt(8 / (n * L))
pi <- rep(pii, L)
beta <- 1
delta <- c(rep(1, r), rep(0, L - r))
Z <- mvrnorm(n, rep(0, L), diag(L))
sigma_e <- (1 + 0.3 * rowSums(Z[, (r + 1):L])) * rnorm(n)
u <- 0.2 * sigma_e + sqrt(1 - 0.2^2) * rnorm(n)
X <- Z %*% pi + u
Y <- X * beta + Z %*% delta + sigma_e

# Perform causal analysis bootstrap
# B "Number of bootstrap samples (default is 500)"
results <- pive(Y, X, Z, bootstrap = TRUE, B = 100, alpha = 0.05)
print(results)


# Lasso-type jackknife IV methods

Pz = Z %*% solve(t(Z) %*% Z) %*% t(Z); W = cbind(Y, X)
iZZ<- solve(t(Z)%*%Z)
h <- vector("numeric",length(Y))            
for(i in 1:length(Y)){
h[i] <- t(Z[i,])%*%iZZ%*%Z[i,]               
    }#n           
H = diag(h)
klimlj <- min(eigen(ginv(t(W)  %*% W) %*% (t(W) %*% (Pz - H) %*% W ))$values)
kfulj  <- (klimlj - (1-klimlj)/(n))*(1 - (1-klimlj)/(n))

cv.LJIVE(Y,X,Z,k = 1)
cv.LJIVE(Y,X,Z,k = klimlj)
cv.LJIVE(Y,X,Z,k = kfulj)
 
# Penalized K-Class IV method
cv.PKCIVE(Y,X,Z,k = 1)
