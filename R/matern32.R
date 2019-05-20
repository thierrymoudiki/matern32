#' Fit Matern 3/2 model
#'
#' @param x 
#' @param y 
#' @param sigma 
#' @param l 
#' @param lambda 
#' @param get_derivatives
#' @param ... 
#'
#' @return
#' @export
#'
#' @examples
#' 
#' n <- 10; p <- 4
#' 
#' set.seed(456)
#' X <- matrix(rnorm(n * p), n, p) # no intercept!
#' y <- rnorm(n)
#' 
#' lams <- 10^seq(-5, 4, length.out = 50)
#' 
#' # use matern32::find_params_matern32 to find sigma and l
#' fit_obj <- matern32::fit_matern32(x = X, y = y, lambda = lams)
#'
#' par(mfrow=c(1, 2))
#'  
#' plot(log(lams), fit_obj$GCV, type = 'l', main = "GCV", 
#' ylab = "GCV")
#' 
#' matplot(log(lams), t(fit_obj$coef), type = 'l', 
#' main = "coefficients = f(lambda)", xlab = "log(lambda)", 
#' ylab = "coefs")
#' abline(h = 0, lty = 2, lwd = 2, col = "red")
#' 
#' 
#' library(MASS)
#' X <-  longley[,-7]
#' y <- longley[, 7]
#' fit_obj <- matern32::fit_matern32(x = X, y = y, lambda = lams)
#' 
#' par(mfrow=c(1, 2))
#'  
#' plot(log(lams), fit_obj$GCV, type = 'l', main = "GCV", 
#' ylab = "GCV")
#' 
#' matplot(log(lams), t(fit_obj$coef), type = 'l', 
#' main = "coefficients = f(lambda)", xlab = "log(lambda)", 
#' ylab = "coefs")
#' abline(h = 0, lty = 2, lwd = 2, col = "red")
#' 
fit_matern32 <- function(x, y, lambda = 10^seq(-5, 4, length.out = 100),
                         l = NULL, get_derivatives = TRUE, 
                         method = c("svd", "chol", "eigen"),
                         ...)
{
  method <- match.arg(method)
  
  ## regression ----
  x <- as.matrix(x)
  y <- as.vector(y)
  nlambda <- length(lambda)
  n <- dim(x)[1]
  p <- dim(x)[2]
  stopifnot(n == length(y))
  
  # centered response
  ym <- mean(y)
  centered_y <- y - ym
  
  # construct covariance
  x_scaled <- matern32::my_scale(x)
  X <- x_scaled$res
  
  if (is.null(l))
    l <- sqrt(p)
  
  # compute kernel as a vector if necessary
  if (length(l) == 1)
    l <- rep(l, p)
  
  K <- matern32_kxx_cpp(x = X, l = l) 
  
  if (method == "svd")
  {
    if(n > 500)
      cat("Processing...", "\n")
    Xs <- La.svd(K)
    rhs <- crossprod(Xs$u, centered_y)
    d <- Xs$d
    nb_di <- length(d)
    div <- d ^ 2 + rep(lambda, rep(nb_di, nlambda))
    a <- drop(d * rhs) / div
    dim(a) <- c(nb_di, nlambda)
    coef <- crossprod(Xs$vt, a)
    colnames(coef) <- lambda
    
    centered_y_hat <- K %*% coef
    fitted_values <- drop(ym +  centered_y_hat)
    resid <- centered_y - centered_y_hat
    colnames(resid) <- lambda
    GCV <- colSums(resid^2)/(nrow(X) - colSums(matrix(d^2/div, 
                                                      nb_di)))^2
    
    if (length(lambda) > 1)
    {
      RSS <- colSums((y - fitted_values)^2)
    } else {
      RSS <- sum((y - fitted_values)^2)
    }
    
    TSS <- sum((y - ym)^2)
    R_Squared <- 1 - RSS/TSS
    names(R_Squared) <- lambda
    
    return(list(K = K, l = l, 
                coef = drop(coef), scales = x_scaled$xsd,
                ym = ym, xm = x_scaled$xm,
                fitted_values = fitted_values, resid = resid,
                GCV = GCV, R_Squared = R_Squared, 
                scaled_x = X, centered_y = centered_y))
  }
  
  if (method == "chol")
  {
    if (length(lambda) <= 1)
    {
      
      K_plus <- K + lambda*diag(n)
      invK <- chol2inv(chol(K_plus))
      coef <- invK%*%centered_y
      loocv <- drop(coef/diag(invK))
      
    } else { # length(lambda) > 1
      
      get_loocv <- function(lambda_i)
      {
        K_plus <- K + lambda_i*diag(n)
        invK <- chol2inv(chol(K_plus))
        coef <- invK%*%centered_y
        return(list(coef = coef,
                    loocv = drop(coef/diag(invK))))
      }
      
      fit_res <- lapply(lambda, function(x) get_loocv(x))
      n_fit_res <- length(fit_res)
      
      loocvs <- sqrt(colMeans(sapply(1:n_fit_res,  
                       function(i) fit_res[[i]]$loocv)^2))
      names(loocvs) <- lambda
      
      coefs <- sapply(1:n_fit_res,  function(i) fit_res[[i]]$coef)
      colnames(coefs) <- lambda
      
    }
  }
  
  if (method == "eigen")
  {
    if(length(lambda) <= 1)
    {
      eigenK <- base::eigen(K)
      eigen_values <- eigenK$values
      Q <- eigenK$vectors 
      inv_eigen <- solve_eigen(Eigenvectors = Q,
                               Eigenvalues = eigen_values,
                               y = centered_y,
                               lambda = lambda)
      
      coef <- inv_eigen$coeffs
      loocv <- inv_eigen$loocv 
      
    } else { # length(lambda) > 1
      
        get_loocv <- function(lambda_i)
        {
          eigenK <- base::eigen(K)
          eigen_values <- eigenK$values
          Q <- eigenK$vectors 
          inv_eigen <- solve_eigen(Eigenvectors = Q,
                                   Eigenvalues = eigen_values,
                                   y = centered_y,
                                   lambda = lambda_i)
          return(list(coef = inv_eigen$coef,
                      loocv = inv_eigen$loocv))
        }
        
        fit_res <- lapply(lambda, function(x) get_loocv(x))
        n_fit_res <- length(fit_res)
        
        loocvs <- sqrt(colMeans(sapply(1:n_fit_res,  
                         function(i) fit_res[[i]]$loocv)^2))
        names(loocvs) <- lambda
        
        coefs <- sapply(1:n_fit_res,  function(i) fit_res[[i]]$coef)
        colnames(coefs) <- lambda
      
      }
    
  }
  
  if (method %in% c("chol", "eigen"))
  {
    if (length(lambda) == 1)
    {
      centered_y_hat <- K %*% coef
      fitted_values <- drop(ym +  centered_y_hat)
      resid <- centered_y - centered_y_hat
      
      RSS <- sum((y - fitted_values)^2)
      TSS <- sum((y - ym)^2)
      R_Squared <- 1 - RSS/TSS
      
      return(list(K = K, l = l, 
                  coef = drop(coef), 
                  scales = x_scaled$xsd,
                  ym = ym, xm = x_scaled$xm,
                  fitted_values = fitted_values, resid = drop(resid),
                  loocv = loocvs, R_Squared = R_Squared, 
                  scaled_x = X, centered_y = centered_y))  
    } else { 
      centered_y_hat <- K %*% coefs
      fitted_values <- drop(ym +  centered_y_hat)
      resid <- centered_y - centered_y_hat
      
      RSS <- colSums((y - fitted_values)^2)
      TSS <- sum((y - ym)^2)
      R_Squared <- 1 - RSS/TSS
      names(R_Squared) <- lambda
    
      return(list(K = K, l = l, 
                  coef = drop(coefs), 
                  scales = x_scaled$xsd,
                  ym = ym, xm = x_scaled$xm,
                  fitted_values = fitted_values, resid = resid,
                  loocv = loocvs, R_Squared = R_Squared, 
                  scaled_x = X, centered_y = centered_y)) 
    }
  }
  
}


#' Predict from Matern 3/2 model
#'
#' @param fit_obj 
#' @param newx 
#' @param ci 
#'
#' @return
#' @export
#'
#' @examples
#' 
#' 
#' n <- 10 ; p <- 4
#' 
#' set.seed(456)
#' X <- matrix(rnorm(n * p), n, p) # no intercept!
#' y <- rnorm(n)
#' 
#' lams <- 10^seq(-5, 4, length.out = 50)
#' 
#' # use matern32::find_params_matern32 to find sigma and l
#' fit_obj <- fit_matern32(x = X, y = y, lambda = lams)
#' 
#' df <- data.frame(predict_matern32(fit_obj, newx = X) - y)
#' colnames(df) <- paste0(round(lams, 2))
#' summary(df)
#' boxplot(df[, c(1, 10, 25, 35, 50)], 
#' main = "distribution of bias", 
#' xlab = "lambda", ylab = "y_hat - y")
#' 
predict_matern32 <- function(fit_obj, newx, ci = NULL)
{
  if (is.vector(newx))
    newx <- t(newx)
  
  K_star <- matern32_kxstar_cpp(newx = matern32::my_scale(x = newx,
                                                          xm = as.vector(fit_obj$xm),
                                                          xsd = as.vector(fit_obj$scales)), 
                                x = fit_obj$scaled_x, 
                                l = fit_obj$l)
  
  return(drop(crossprod(K_star, fit_obj$coef)) + fit_obj$ym)
}

