#' Lasso-type instrumental variables methods for estimating causal effects and selecting invalid instruments.
#'
#' "pive" function applies six methods to estimate coefficients and performs bootstrap sampling to provide estimates, standard errors, and confidence intervals.
#' If bootstrap is FALSE, it provides causal estimates of all methods and the number of invalid instrumental variables.
#'
# Load necessary packages
#library(boot)
#library(MASS)
#'
#' @param Y Response variable
#' @param X Exposure variable
#' @param Z Instrumental variables matrix
#' @param bootstrap Logical, whether to perform bootstrap sampling (default is TRUE)
#' @param B Number of bootstrap samples (default is 500)
#' @param alpha Significance level for confidence intervals (default is 0.05)
#' @param K K-fold cross-validation (default is 10-fold CV)
#' @importFrom stats coef predict quantile sd
#' @import boot
#'
#' @return Estimated results
#' @export
pive <- function(Y, X, Z, bootstrap = TRUE, B = 500, alpha = 0.05, K = 10) {

  # Helper function to calculate parameters
  calc_param <- function(Y, X, Z) {
    n  = nrow(Z)
    Pz = Z %*% solve(t(Z) %*% Z) %*% t(Z); W = cbind(Y, X)
    iZZ <- solve(t(Z)%*%Z);
    h   <- vector("numeric",length(Y));
    for(ii in 1:length(Y)){
      h[ii]   <- t(Z[ii,])%*%iZZ%*%Z[ii,]
    }#n
    klimln <- min(eigen(solve(t(W)  %*% W) %*% (t(W) %*% (Pz) %*% W ))$values)
    kfuln  <- (klimln - (1-klimln)/(n))*(1 - (1-klimln)/(n))

    H = diag(h)
    klimlnj <- min(eigen(solve(t(W)  %*% W) %*% (t(W) %*% (Pz - H) %*% W ))$values)
    kfulnj  <- (klimlnj - (1-klimlnj)/(n))*(1 - (1-klimlnj)/(n))

    list(klimln = klimln, kfuln = kfuln, klimlnj = klimlnj, kfulnj = kfulnj)
  }
  df <- data.frame(y = Y, X = X, Z)

  # Bootstrap function to fit models
  fit_models <- function(data, indices, alpha, params) {
    resampled_data <- data[indices, ]
    Y <- resampled_data$y
    X <- resampled_data$X
    Z <- resampled_data[, -c(1, 2)]
    Z <- cbind(apply(Z,2,as.numeric))
    params <- calc_param(Y, X, Z)

    # Perform model fits using parameters
    beta_tsls  <- cv.PKCIVE(Y, X, Z, k = 1, K = K)$beta
    beta_liml  <- cv.PKCIVE(Y, X, Z, k = params$klimln, K = K)$beta
    beta_ful   <- cv.PKCIVE(Y, X, Z, k = params$kfuln, K = K)$beta
    beta_pjtsls<- cv.LJIVE(Y, X, Z, k = 1, K = K)$beta
    beta_pjliml<- cv.LJIVE(Y, X, Z, k = params$klimlnj, K = K)$beta
    beta_pjful <- cv.LJIVE(Y, X, Z, k = params$kfulnj, K = K)$beta

    # Return coefficients of PKCIVE and LJIVE
    c(beta_tsls, beta_liml, beta_ful, beta_pjtsls, beta_pjliml, beta_pjful)
  }

  # Main function body
  params <- calc_param(Y, X, Z)

  if (bootstrap) {
    # Perform bootstrap
    set.seed(123)  # Set seed
    boot_results <- boot(data = df, statistic = fit_models, R = B, alpha = alpha, params = params)

    # Extract bootstrap results
    boot_coefs <- boot_results$t
    boot_mean  <- apply(boot_coefs, 2, mean)
    boot_se    <- apply(boot_coefs, 2, sd)
    boot_ci    <- t(apply(boot_coefs, 2, function(x) quantile(x, c(alpha/2, 1 - alpha/2))))

    # Create results data frame
    results <- data.frame(
      Method = c("PTSLS", "PLIML", "PFUL", "PJTSLS", "PJLIML", "PJFUL"),
      Estimate = c(cv.PKCIVE(Y, X, Z, k = 1, K = K)$beta,
                   cv.PKCIVE(Y, X, Z, k = params$klimln, K = K)$beta,
                   cv.PKCIVE(Y, X, Z, k = params$kfuln, K = K)$beta,
                   cv.LJIVE(Y, X, Z, k = 1, K = K)$beta,
                   cv.LJIVE(Y, X, Z, k = params$klimlnj, K = K)$beta,
                   cv.LJIVE(Y, X, Z, k = params$kfulnj, K = K)$beta),
      Invalid_Count = NA,  # To be filled later
      Invalid_Instruments = NA,  # To be filled later
      Bootstrap_Estimate = boot_mean,
      SE = boot_se,
      CI_Lower = boot_ci[, 1],
      CI_Upper = boot_ci[, 2],
      stringsAsFactors = FALSE
    )
  } else {
    # If not bootstrapping, calculate estimates without bootstrap
    beta_tsls   <- cv.PKCIVE(Y, X, Z, k = 1, K = K)$beta
    beta_liml   <- cv.PKCIVE(Y, X, Z, k = params$klimln, K = K)$beta
    beta_ful    <- cv.PKCIVE(Y, X, Z, k = params$kfuln, K = K)$beta
    beta_pjtsls <- cv.LJIVE(Y, X, Z, k = 1, K = K)$beta
    beta_pjliml <- cv.LJIVE(Y, X, Z, k = params$klimlnj, K = K)$beta
    beta_pjful  <- cv.LJIVE(Y, X, Z, k = params$kfulnj, K = K)$beta

    results <- data.frame(
      Method = c("PTSLS", "PLIML", "PFUL", "PJTSLS", "PJLIML", "PJFUL"),
      Estimate = c(beta_tsls, beta_liml, beta_ful, beta_pjtsls, beta_pjliml, beta_pjful),
      Invalid_Count = NA,        # To be filled later
      Invalid_Instruments = NA,  # To be filled later
      Bootstrap_Estimate = NA,
      SE = NA,
      CI_Lower = NA,
      CI_Upper = NA,
      stringsAsFactors = FALSE
    )
  }

  # Fill in Invalid_Count and Invalid_Instruments
  for (i in 1:nrow(results)) {
    fit <- switch(results$Method[i],
                  "PTSLS"  = cv.PKCIVE(Y, X, Z, k = 1, K = K),
                  "PLIML"  = cv.PKCIVE(Y, X, Z, k = params$klimln, K = K),
                  "PFUL"   = cv.PKCIVE(Y, X, Z, k = params$kfuln, K = K),
                  "PJTSLS" = cv.LJIVE(Y, X, Z, k = 1, K = K),
                  "PJLIML" = cv.LJIVE(Y, X, Z, k = params$klimlnj, K = K),
                  "PJFUL"  = cv.LJIVE(Y, X, Z, k = params$kfulnj, K = K)
    )

    invalid_instruments <- paste(fit$whichInvalid, collapse = ", ")
    invalid_count       <- length(unlist(strsplit(fit$whichInvalid, ",")))

    results[i, c("Invalid_Count", "Invalid_Instruments")] <- list(invalid_count, invalid_instruments)
  }

  colnames(results) <- c("Method", "Estimate", "# IVs", "Invalid instruments",
                         "Bootstrap Estimate", "Std.error", "CI Lower", "CI Upper")

  return(cbind(results))

}

#' Fit the IVs model using LJIVE
#'
#' Load necessary packages
#library(lars)
#library(MASS)
#' @param Y Response variable
#' @param X Exposure variable
#' @param Z Instrumental variables
#' @param k Parameter for K-class jackknife IVs
#' @param intercept Logical, indicating whether to include an intercept (default: TRUE).
#' @param normalize Logical, indicating whether to normalize the instruments (default: TRUE).
#' @import lars
#' @return An object of class "LJIVE" containing various outputs.
#' @export

LJIVE <- function(Y,X,Z,k,intercept=TRUE,normalize=TRUE) {
  # Error checking
  if( (!is.vector(Y) && !is.matrix(Y)) | !is.numeric(Y)) stop("Y is not in numeric vector format")
  if( (!is.vector(X) && !is.matrix(X)) | !is.numeric(X)) stop("X is not in numeric vector format")
  if( !is.matrix(Z) | !is.numeric(Z)) stop("Z is not in numeric vector format")
  if(nrow(Z) != length(Y)) stop("The dimensions of Y and the rows of Z do not match")
  if(nrow(Z) != length(X)) stop("The dimensions of X and the rows of Z do not match")
  if(length(X) != length(Y)) stop("The dimensions of X do not match those of Y")
  if(intercept) {
    if( (ncol(Z) + 1) > nrow(Z)) stop("The number of instruments is greater than the sample size.")
  } else {
    if( ncol(Z) > nrow(Z)) stop("The number of instruments is greater than the sample size.")
  }

  ## Sample size (n) and No. of instruments (L)
  n = nrow(Z); L = ncol(Z)
  if(intercept) {
    meanY = mean(Y); meanX = mean(X); meanZ = colMeans(Z);
    Y = Y - mean(Y); X = X - meanX; Z = scale(Z,center=TRUE,scale=FALSE);
  }
  if(normalize) {
    normZ = sqrt(colSums(Z^2))
    Z = scale(Z,center=FALSE,scale=TRUE) / sqrt(n-1)
  } else {
    normZ = rep(1,L)
  }
  ## JIVE procedure
  Xhatj  <- matrix(0,length(Y),1)
  iZZ    <- solve(t(Z)%*%Z)
  Ga.hat <- (iZZ)%*%(t(Z)%*%X)
  h      <- vector("numeric",length(Y))
  for(ii in 1:length(Y)){
    h[ii]      <- t(Z[ii,])%*%iZZ%*%Z[ii,];
    Xhatj[ii,] <- (t(Z[ii,])%*%Ga.hat - h[ii]*X[ii])/(1-h[ii]);
  }#n
  #PXhatj = Xhatj %*% solve(t(Xhatj) %*% Xhatj) %*% t(Xhatj); #MXhatj = diag(n) - PXhatj

  QR = qr(Z); Yhat = qr.fitted(QR,Y)
  Ztildej = Z - Xhatj %*% t(Xhatj) %*% Z / sum(Xhatj^2); #Ztildej = MXhatj%*%Z
  Ytildej = Yhat - Xhatj * (sum(Xhatj * Yhat) / sum( Xhatj^2));
  Xstarj = Xhatj + (1-k)^2*(X-Xhatj)
  klimldj = ((1-k)^2*sum(X^2) + (1-(1-k)^2)*sum(Xhatj^2))

  fit    <- lars(Ztildej,Ytildej,intercept = FALSE,normalize=TRUE)
  deltaj <- coef(fit)
  deltajSuppSize = apply(deltaj,1,function(x){sum(x != 0)}); indexProperSupp = which(deltajSuppSize < L);
  betaj  = drop(as.numeric(t(Xstarj) %*%Y) - t(Xhatj)%*% Z %*% t(deltaj) ) / klimldj
  deltaj = scale(deltaj,FALSE,normZ)
  attr(deltaj,"scaled:scale") = NULL

  # Packaging Object
  lambda = fit$lambda;
  whichInvalid = apply(deltaj,1,function(x){paste(which(abs(x) > 0),sep="",collapse=",")})
  dimnames(deltaj) = list(1:length(betaj),dimnames(Z)[[2]]);
  names(betaj) = paste(1:length(betaj));
  names(whichInvalid) = paste(1:length(betaj))

  object <- list(call = match.call(), deltaj = deltaj, betaj = betaj, whichInvalid = whichInvalid,nInvalid = deltajSuppSize,lambda = lambda,larsObject = fit, Y=Y, X=X, Z=Z, k=k, Xhatj=Xhatj, Xstarj=Xstarj, klimldj=klimldj, normZ = normZ)
  class(object) <- "LJIVE"
  return(object)
}


#' Print the summary of the model

print.LJIVE <- function(x,...) {
  object = x
  cat("\nCall:\n")
  dput(object$call)
  printOut = data.frame(object$betaj,object$nInvalid, object$whichInvalid)
  colnames(printOut) = c("Estimates of Beta",
                         "   Number of Invalid IVs","Invalid IVs")
  print(printOut,row.names=FALSE)
  invisible(object)
}

#' Predict the model results based on beta and delta parameters
#'
#' @param object A fitted model object of class "LJIVE".
#' @param lambda Numeric vector of lambda values for which predictions are desired.
#' @param type Character string specifying the type of prediction:
#'   - "coefficients" for predicting coefficients.
#'   - "instruments" for predicting instruments.
#' @param ... Additional arguments passed to methods.
#' @import lars
#' @export
predict.LJIVE <- function(object,lambda,type=c("coefficients","instruments"),...) {

  type = match.arg(type)
  if(!(type %in% c("coefficients","instruments"))) stop("For type, specify either `coefficients' or `instruments'")
  if(missing(lambda)) lambda = object$lambda
  deltaj = predict(object$larsObject,s = lambda,type="coefficients",mode="lambda",...)$coefficient
  if(is.vector(deltaj)) deltaj = matrix(deltaj,1,length(deltaj))
  if(type == "coefficients") {
    betaj =  drop(as.numeric(t(object$Xstarj) %*%object$Y) - t(object$Xhatj)%*% object$Z %*% t(deltaj) ) / sum(object$klimldj)
    deltaj = scale(deltaj,FALSE,object$normZ)
    attr(deltaj,"scaled:scale") = NULL
    return(list(lambda = lambda,deltaj = deltaj,betaj=betaj))
  }
  if(type == "instruments") {
    return(list(lambda = lambda, instruments = apply(deltaj,1,function(x){paste(which(abs(x) > 0),sep="",collapse=",")})))
  }

}

#' summary of the model
#' @export
summary.LJIVE <- function(object,...) {
  print.LJIVE(object,...)
}

#' Fit the IVs model and perform Cross-Validation for LJIVE
#'
#' @param Y Response variable
#' @param X Exposure variable
#' @param Z Instrumental variables
#' @param k Parameter for jackknife IVs
#' @param K Number of folds for cross-validation
#' @param intercept Logical, indicating whether to include an intercept (default: TRUE).
#' @param normalize Logical, indicating whether to normalize the instruments (default: TRUE).
#' @import lars
#' @return An object of class "LJIVE" containing various outputs.
#' @export


cv.LJIVE <- function(Y,X,Z,k,lambdaSeq,K = 10,intercept=TRUE,normalize=TRUE) {
  # Error checking
  if( (!is.vector(Y) && !is.matrix(Y)) | !is.numeric(Y)) stop("Y is not in numeric vector format")
  if( (!is.vector(X) && !is.matrix(X)) | !is.numeric(X)) stop("X is not in numeric vector format")
  if( !is.matrix(Z) | !is.numeric(Z)) stop("Z is not in numeric vector format")
  if(nrow(Z) != length(Y)) stop("The dimensions of Y and the rows of Z do not match")
  if(nrow(Z) != length(X)) stop("The dimensions of X and the rows of Z do not match")
  if(length(X) != length(Y)) stop("The dimensions of X do not match those of Y")

  if(intercept) {
    if( (ncol(Z) + 1) > nrow(Z)) stop("The number of instruments is greater than the sample size.")
  } else {
    if( ncol(Z) > nrow(Z)) stop("The number of instruments is greater than the sample size.")
  }
  if(!is.numeric(K) | K > length(Y)) stop("K is not a proper numeric number")
  fitall = LJIVE(Y,X,Z,k,intercept=intercept,normalize=normalize)
  if(missing(lambdaSeq) || all(is.na(lambdaSeq))) {
    #warning("Lambda sequence not provided; defaulting to using lambdas provided by LJIVE")
    lambdaSeq = c(fitall$lambda,seq(from=0,to=2*max(fitall$lambda,na.rm=TRUE),length.out = 100))
    lambdaSeq = sort(unique(lambdaSeq))
  }
  if(any(is.na(lambdaSeq))) {
    warning("Some lambda values are missing. Ignoring these lambda values for cross-validation")
    lambdaSeq = lambdaSeq[!is.na(lambdaSeq)]
  }
  if(length(lambdaSeq) < 2) stop("Only one lambda provided. Please provide multiple lambdas")
  lambdaSeq = sort(lambdaSeq,decreasing=TRUE)

  n = nrow(Z); L = ncol(Z)
  # Cross validation
  sampleIndex = sample(rep(1:K,length.out= n))
  errormat = matrix(0,K,length(lambdaSeq))
  for(i in 1:K) {
    testSet = (sampleIndex == i)
    trainfit = LJIVE(Y[!testSet], X[!testSet], Z[!testSet, , drop=FALSE],k, intercept = intercept,normalize = normalize)
    trainfitCoef = predict(trainfit, lambda = lambdaSeq, type = "coefficients")
    Y.test = Y[testSet]; X.test = X[testSet]; Z.test = Z[testSet,,drop=FALSE]
    if(intercept) {
      meanY.test = mean(Y.test); meanX.test = mean(X.test); meanZ.test = colMeans(Z.test)
      Y.test = Y.test - meanY.test; X.test = X.test - meanX.test; Z.test = scale(Z.test,center=TRUE,scale=FALSE)
    }

    QR = qr(Z.test)
    residTest = (as.numeric(Y.test) - Z.test %*% t(trainfitCoef$deltaj) - X.test %*% t(trainfitCoef$betaj))
    errormat[i,] = colSums(qr.fitted(QR, residTest)^2) #does (PZ %*% residual)^2 summed across observations
  }
  cv = colMeans(errormat)

  if(all(is.nan(cv))) {
    warning("All lambdas were invalid. Please try different values of lambda for cross-validation to work")
    return(list(lambda = rep(NA,length(lambdaSeq)),estCVError = NA, deltaj = rep(NA,ncol(Z)),betaj = NA, whichInvalid = NA))
  } else {
    stderror = apply(errormat,2,function(x){sd(x)/sqrt(K)})
    mincv.index = which.min(cv); onestderrorbound = cv[mincv.index] + stderror[mincv.index]
    onestdLambda = max(lambdaSeq[which( (cv <= onestderrorbound) & (cv >= cv[mincv.index]))]) #this is never empty vector
    returnOut = predict(fitall,onestdLambda,type="coefficients")
    return(list(lambda = onestdLambda, estCVError = cv[which(onestdLambda == lambdaSeq)], deltaj = drop(returnOut$deltaj), betaj = drop(returnOut$betaj),whichInvalid = paste(which(abs(returnOut$deltaj) > 0),sep="",collapse=",")))
  }
}

#' Fit the IVs model using PKCIVE
#'
#' @param Y Response variable
#' @param X Exposure variable
#' @param Z Instrumental variables
#' @param k Parameter for K-class jackknife IVs
#' @param intercept Logical, indicating whether to include an intercept (default: TRUE).
#' @param normalize Logical, indicating whether to normalize the instruments (default: TRUE).
#' @return An object of class "LJIVE" containing various outputs.
#' @import lars
#' @export
#'
PKCIVE <- function(Y,X,Z,k,intercept=TRUE,normalize=TRUE) {
  # Error checking
  if( (!is.vector(Y) && !is.matrix(Y)) | !is.numeric(Y)) stop("Y is not in numeric vector format")
  if( (!is.vector(X) && !is.matrix(X)) | !is.numeric(X)) stop("X is not in numeric vector format")
  if( !is.matrix(Z) | !is.numeric(Z)) stop("Z is not in numeric vector format")
  if(nrow(Z) != length(Y)) stop("The dimensions of Y and the rows of Z do not match")
  if(nrow(Z) != length(X)) stop("The dimensions of X and the rows of Z do not match")
  if(length(X) != length(Y)) stop("The dimensions of X do not match those of Y")
  if(intercept) {
    if( (ncol(Z) + 1) > nrow(Z)) stop("The number of instruments exceeds the sample size")
  } else {
    if( ncol(Z) > nrow(Z)) stop("The number of instruments exceeds the sample size")
  }

  ## Sample size (n) and No. of instruments (L)
  n = nrow(Z); L = ncol(Z)
  if(intercept) {
    meanY = mean(Y); meanX = mean(X); meanZ = colMeans(Z);
    Y = Y - mean(Y); X = X - meanX; Z = scale(Z,center=TRUE,scale=FALSE);
  }
  if(normalize) {
    normZ = sqrt(colSums(Z^2))
    Z = scale(Z,center=FALSE,scale=TRUE) / sqrt(n-1)
  } else {
    normZ = rep(1,L)
  }
  QR = qr(Z); Yhat = qr.fitted(QR,Y); Xhat = qr.fitted(QR,X);
  Ztilde = Z - Xhat %*% t(Xhat) %*% Z / sum( Xhat^2)
  Ytilde = Yhat - Xhat * (sum(Xhat * Yhat) / sum( Xhat^2))
  #Pz = Z %*% solve(t(Z) %*% Z) %*% t(Z); Mz = diag(n) - Pz
  #PXhat = Xhat %*% solve(t(Xhat) %*% Xhat) %*% t(Xhat); MXhat = diag(n) - PXhat
  #Ytilde = MXhat%*%Pz%*%Y; Ztilde = MXhat%*%Z;
  ##k is estimated as the smallest eigenvalue using the observed dataset
  Xstar  = Xhat + (1-k)^2*(X-Xhat)
  klimld = sum(Xhat^2) + (1-k)^2*(sum(X^2) - sum(Xhat^2))
  fit    = lars(Ztilde,Ytilde,intercept = FALSE,normalize=TRUE)
  delta  = coef(fit)
  deltaSuppSize = apply(delta,1,function(x){sum(x != 0)}); indexProperSupp = which(deltaSuppSize < L);
  beta = drop(as.numeric(t(Xstar) %*%Y) - t(Xhat)%*% Z %*% t(delta) ) / klimld
  delta = scale(delta,FALSE,normZ)
  attr(delta,"scaled:scale") = NULL

  # Packaging Object
  lambda = fit$lambda;
  whichInvalid = apply(delta,1,function(x){paste(which(abs(x) > 0),sep="",collapse=",")})
  dimnames(delta) = list(1:length(beta),dimnames(Z)[[2]]);
  names(beta)     = paste(1:length(beta));
  names(whichInvalid) = paste(1:length(beta))

  object <- list(call = match.call(), delta = delta, beta = beta, whichInvalid = whichInvalid,nInvalid = deltaSuppSize,lambda = lambda,larsObject = fit, Y=Y, X=X, Z=Z, k=k, Xhat=Xhat, Xstar=Xstar, klimld=klimld, normZ = normZ)
  class(object) <- "PKCIVE"
  return(object)
}

#' Print the summary of the model

print.PKCIVE  <- function(x,...) {
  object = x
  cat("\nCall:\n")
  dput(object$call)
  printOut = data.frame(object$beta,object$nInvalid, object$whichInvalid)
  colnames(printOut) = c("Estimates of Beta",
                         "   Number of Invalid IVs","Invalid IVs")
  print(printOut,row.names=FALSE)
  invisible(object)
}

#' Predict the model results based on beta and delta parameters
#'
#' @param object A fitted model object of class "PKCIVE".
#' @param lambda Numeric vector of lambda values for which predictions are desired.
#' @param type Character string specifying the type of prediction:
#'   - "coefficients" for predicting coefficients.
#'   - "instruments" for predicting instruments.
#' @param ... Additional arguments passed to methods.
#' @export

predict.PKCIVE  <- function(object,lambda,type=c("coefficients","instruments"),...) {
  type = match.arg(type)
  if(!(type %in% c("coefficients","instruments"))) stop("For type, specify either `coefficients' or `instruments'")
  if(missing(lambda)) lambda = object$lambda
  delta = predict(object$larsObject,s = lambda,type="coefficients",mode="lambda",...)$coefficient
  if(is.vector(delta)) delta = matrix(delta,1,length(delta))
  if(type == "coefficients") {
    beta  = drop(as.numeric(t(object$Xstar) %*%object$Y) - t(object$Xhat)%*% object$Z %*% t(delta) ) / sum(object$klimld)
    delta = scale(delta,FALSE,object$normZ)
    attr(delta,"scaled:scale") = NULL
    return(list(lambda = lambda,delta = delta,beta=beta))
  }
  if(type == "instruments") {
    return(list(lambda = lambda, instruments = apply(delta,1,function(x){paste(which(abs(x) > 0),sep="",collapse=",")})))
  }

}
#' summary of the model
#' @export

summary.PKCIVE  <- function(object,...) {
  print.PKCIVE (object,...)
}

#' Fit the IVs model and perform Cross-Validation for PKCIVE
#'
#' @param Y Response variable
#' @param X Exposure variable
#' @param Z Instrumental variables
#' @param k Parameter for K-class IVs
#' @param K Number of folds for cross-validation
#' @param intercept Logical, indicating whether to include an intercept (default: TRUE).
#' @param normalize Logical, indicating whether to normalize the instruments (default: TRUE).
#' @return An object of class "LJIVE" containing various outputs.
#' @export

cv.PKCIVE  <- function(Y,X,Z,k,lambdaSeq,K = 10,intercept=TRUE,normalize=TRUE) {
  # Error checking
  if( (!is.vector(Y) && !is.matrix(Y)) | !is.numeric(Y)) stop("Y is not in numeric vector format")
  if( (!is.vector(X) && !is.matrix(X)) | !is.numeric(X)) stop("X is not in numeric vector format")
  if( !is.matrix(Z) | !is.numeric(Z)) stop("Z is not in numeric vector format")
  if(nrow(Z) != length(Y)) stop("The dimensions of Y and the rows of Z do not match")
  if(nrow(Z) != length(X)) stop("The dimensions of X and the rows of Z do not match")
  if(length(X) != length(Y)) stop("The dimensions of X do not match those of Y")

  if(intercept) {
    if( (ncol(Z) + 1) > nrow(Z)) stop("The number of instruments is greater than the sample size.")
  } else {
    if( ncol(Z) > nrow(Z)) stop("The number of instruments is greater than the sample size.")
  }
  if(!is.numeric(K) | K > length(Y)) stop("K is not a proper numeric number")
  fitall = PKCIVE (Y,X,Z,k,intercept=intercept,normalize=normalize)
  if(missing(lambdaSeq) || all(is.na(lambdaSeq))) {
    #warning("Lambda sequence not provided; defaulting to using lambdas provided by PKCIVE ")
    lambdaSeq = c(fitall$lambda,seq(from=0,to=2*max(fitall$lambda,na.rm=TRUE),length.out = 100))
    lambdaSeq = sort(unique(lambdaSeq))
  }
  if(any(is.na(lambdaSeq))) {
    warning("Some lambda values are missing. Ignoring these lambda values for cross-validation")
    lambdaSeq = lambdaSeq[!is.na(lambdaSeq)]
  }
  if(length(lambdaSeq) < 2) stop("Only one lambda provided. Please provide multiple lambdas")
  lambdaSeq = sort(lambdaSeq,decreasing=TRUE)
  n = nrow(Z); L = ncol(Z);

  # Cross validation
  sampleIndex = sample(rep(1:K,length.out= n))
  errormat = matrix(0,K,length(lambdaSeq))
  for(i in 1:K) {
    testSet = (sampleIndex == i)
    trainfit = PKCIVE (Y[!testSet], X[!testSet], Z[!testSet, , drop=FALSE],k, intercept = intercept,normalize = normalize)
    trainfitCoef = predict(trainfit, lambda = lambdaSeq, type = "coefficients")
    Y.test = Y[testSet]; X.test = X[testSet]; Z.test = Z[testSet,,drop=FALSE]
    if(intercept) {
      meanY.test = mean(Y.test); meanX.test = mean(X.test); meanZ.test = colMeans(Z.test)
      Y.test = Y.test - meanY.test; X.test = X.test - meanX.test; Z.test = scale(Z.test,center=TRUE,scale=FALSE)
    }

    QR = qr(Z.test)
    residTest = (as.numeric(Y.test) - Z.test %*% t(trainfitCoef$delta) - X.test %*% t(trainfitCoef$beta))
    errormat[i,] = colSums(qr.fitted(QR, residTest)^2)
  }
  cv = colMeans(errormat)

  if(all(is.nan(cv))) {
    warning("All lambdas were invalid. Please try different values of lambda for cross-validation to work")
    return(list(lambda = rep(NA,length(lambdaSeq)),estCVError = NA, delta = rep(NA,ncol(Z)),beta = NA, whichInvalid = NA))
  } else {
    stderror = apply(errormat,2,function(x){sd(x)/sqrt(K)})
    mincv.index = which.min(cv); onestderrorbound = cv[mincv.index] + stderror[mincv.index]
    onestdLambda = max(lambdaSeq[which( (cv <= onestderrorbound) & (cv >= cv[mincv.index]))]) #this is never empty vector
    returnOut = predict(fitall,onestdLambda,type="coefficients")
    return(list(lambda = onestdLambda, estCVError = cv[which(onestdLambda == lambdaSeq)], delta = drop(returnOut$delta), beta = drop(returnOut$beta),whichInvalid = paste(which(abs(returnOut$delta) > 0),sep="",collapse=",")))
  }
}

