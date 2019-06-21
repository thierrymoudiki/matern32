#' Fit Matern 3/2 model
#'
#' @param x 
#' @param y 
#' @param sigma 
#' @param l 
#' @param lambda 
#' @param ... 
#'
#' @return
#' @export
#'
#' @examples
#' 
#' n <- 100; p <- 4
#' 
#' set.seed(456)
#' X <- matrix(rnorm(n * p), n, p) # no intercept!
#' y <- rnorm(n)
#' 
#' lams <- 10^seq(-5, 4, length.out = 50)
#' 
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
#' fit_obj <- matern32::fit_matern32(x = X, y = y)
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
#'
#'X <- as.matrix(Boston[,-14])
#'y <- Boston[,14]
#'log_y <- log(y)
#'
#'fit_obj <- matern32::fit_matern32(x = X, y = log_y, lambda = lams, 
#'with_kmeans = TRUE, centers = 10)
#'
#'summary(fit_obj)
#'  
fit_matern32 <- function(x, y, lambda = 10^seq(-5, 4, length.out = 100),
                         l = NULL, method = c("svd", "chol", "eigen"),
                         with_kmeans = FALSE, centers = NULL, 
                         seed = 123, cl = NULL, ...)
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
  
  if (with_kmeans == FALSE) 
    K <- matern32_kxx_cpp(x = X, l = l)  
  
  if(n > 500) # can use kmeans
  {
    if (with_kmeans == TRUE)
    {
      # adjust KRR to centers = X and this new response = y
      if (is.null(centers))
      {
        # find best k in kmeans if 'cl' is NULL
        # see https://uc-r.github.io/kmeans_clustering#optimal (fastest hun?! :D) 
        stop("'centers' not provided: choose best k in kmeans")
        centers <- 2 
      } else { #is.null(centers) == FALSE
        # fitted values, residuals, etc => predict on entire dataset from reduced kernel 
        set.seed(seed)
        cclust_obj <- cclust::cclust(x = as.matrix(X), 
                                     centers = centers)
        # new training set (X_clust, y_clust)
        X_clust <- as.matrix(cclust_obj$centers)
        centered_y_clust <- sapply(1:centers, 
                                   function(i) mean(centered_y[which(cclust_obj$cluster == i)]))
        K <- matern32_kxx_cpp(x = X_clust, l = l)  
      }
      
    } else { # with_kmeans == FALSE
      cat("Processing... (try using option 'with_kmeans')", "\n") 
    }
    
  } else { # if (n <= 500)
    
    if (with_kmeans == TRUE)
    {
      warning("option 'with_kmeans' not implemented for n_obs <= 500")
    } else {
      K <- matern32_kxx_cpp(x = X, l = l)  
    }
    
  }
  
  if (with_kmeans == FALSE)
  {
    if (method == "svd")
    {
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
      Adj_R_Squared <- 1 - (1 - R_Squared)*((n - 1)/(n - p - 1))
      
      res <- list(K = K, l = l, 
                  lambda = lambda,
                  coef = drop(coef), scales = x_scaled$xsd,
                  ym = ym, xm = x_scaled$xm,
                  fitted_values = fitted_values, resid = resid,
                  GCV = GCV, R_Squared = R_Squared, 
                  Adj_R_Squared = Adj_R_Squared,
                  scaled_x = X, x = x, centered_y = centered_y, 
                  fit_method = method)
      
      class(res) <- "matern32"
      
      return(res)
    }
    
    if (method == "chol")
    {
      if (length(lambda) <= 1)
      {
        K_plus <- K + lambda*diag(n)
        invK <- chol2inv(chol(K_plus))
        coef <- invK%*%centered_y
        loocv <- sum(drop(coef/diag(invK))^2)
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
        
        loocv <- colSums(sapply(1:n_fit_res,  
                                function(i) fit_res[[i]]$loocv)^2)
        names(loocv) <- lambda
        
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
        loocv <- sum(inv_eigen$loocv^2)
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
        
        loocv <- colSums(sapply(1:n_fit_res,  
                                function(i) fit_res[[i]]$loocv)^2)
        names(loocv) <- lambda
        
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
        Adj_R_Squared <- 1 - (1 - R_Squared)*((n - 1)/(n - p - 1))
        
        res <- list(K = K, l = l, 
                    lambda = lambda,
                    coef = drop(coef), 
                    scales = x_scaled$xsd,
                    ym = ym, xm = x_scaled$xm,
                    fitted_values = fitted_values, resid = drop(resid),
                    loocv = loocv, R_Squared = R_Squared, 
                    Adj_R_Squared = Adj_R_Squared, 
                    scaled_x = X, centered_y = centered_y, 
                    fit_method = method)
        
        class(res) <- "matern32"
        
        return(res)  
      } else { 
        centered_y_hat <- K %*% coefs
        fitted_values <- drop(ym +  centered_y_hat)
        resid <- centered_y - centered_y_hat
        
        RSS <- colSums((y - fitted_values)^2)
        TSS <- sum((y - ym)^2)
        R_Squared <- 1 - RSS/TSS
        names(R_Squared) <- lambda
        Adj_R_Squared <- 1 - (1 - R_Squared)*((n - 1)/(n - p - 1))
        
        res <- list(K = K, l = l, 
                    lambda = lambda,
                    coef = drop(coefs), 
                    scales = x_scaled$xsd,
                    ym = ym, xm = x_scaled$xm,
                    fitted_values = fitted_values, resid = resid,
                    loocv = loocv, R_Squared = R_Squared, 
                    Adj_R_Squared = Adj_R_Squared, 
                    scaled_x = X, centered_y = centered_y, 
                    fit_method = method)
        
        class(res) <- "matern32"
        
        return(res) 
      }
    } 
    
  } else { # with_kmeans == TRUE
    
    if (n > 500)
    {
      
      if (method %in% c("chol", "eigen"))
        stop("'method' not implemented")
      
      if (method == "svd")
      {
        Xs <- La.svd(K) # K is based on clustered X 
        rhs <- crossprod(Xs$u, centered_y_clust)
        d <- Xs$d
        nb_di <- length(d)
        div <- d ^ 2 + rep(lambda, rep(nb_di, nlambda))
        a <- drop(d * rhs) / div
        dim(a) <- c(nb_di, nlambda)
        coef <- crossprod(Xs$vt, a)
        colnames(coef) <- lambda
        scales <- x_scaled$xsd
        xm <- x_scaled$xm
        
        K_star <- matern32_kxstar_cpp(newx = X, # X is already scaled
                                      x = X_clust, 
                                      l = l)
        
        centered_y_hat <- K_star%*%coef 
        fitted_values <- drop(ym + centered_y_hat)
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
        Adj_R_Squared <- 1 - (1 - R_Squared)*((n - 1)/(n - p - 1))
        names(R_Squared) <- lambda
        
        res <- list(K = K, l = l, 
                    lambda = lambda,
                    coef = drop(coef), scales = scales,
                    ym = ym, xm = xm,
                    fitted_values = fitted_values, resid = resid,
                    GCV = GCV, R_Squared = R_Squared, 
                    Adj_R_Squared = Adj_R_Squared,
                    scaled_x = X, x = x, 
                    with_kmeans = TRUE,
                    cclust_obj = cclust_obj, 
                    scaled_x_clust = X_clust, 
                    centered_y = centered_y, 
                    fit_method = method)
        
        class(res) <- "matern32"
        
        return(res)
      }
      
    } else {
      warning("option 'with_kmeans' not implemented for n_obs <= 500")
    }
    
  }
}
fit_matern32 <- memoise::memoize(fit_matern32)