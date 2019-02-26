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
#' 
#' n <- 10 ; p <- 4
#' 
#' 
#' X <- matrix(rnorm(n * p), n, p) # no intercept!
#' y <- rnorm(n)
#' lams <- 10^seq(-5, 4, length.out = 50)
#' 
#' fit_obj <- fit_matern32(x = X, y = y, lambda = lams)
#' 
#' plot(log(lams), fit_obj$GCV, type = 'l', main = "GCV", 
#' ylab = "GCV")
#' 
#' matplot(log(lams), t(fit_obj$coef), type = 'l', 
#' main = "coefficients = f(lambda)", xlab = "log(lambda)", ylab = "coefs")
#' abline(h = 0, lty = 2, lwd = 2, col = "red")
#' 
#' 
fit_matern32 <- function(x, y, sigma = 2, l = 0.1, ...)
{
  ## regression ----
  x <- as.matrix(x)
  y <- as.vector(y)
  lambda = 10^seq(-5, 4, length.out = 100)
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
  
  # compute kernel
  if (length(l) == 1)
    l <- rep(l, p)
  
  K <- matern32_kxx_cpp(x = X, sigma = sigma, 
                        l = l) 
  
  Xs <- La.svd(K)
  rhs <- crossprod(Xs$u, centered_y)
  d <- Xs$d
  nb_di <- length(d)
  div <- d ^ 2 + rep(lambda, rep(nb_di, nlambda))
  a <- drop(d * rhs) / div
  dim(a) <- c(nb_di, nlambda)
  coef <- crossprod(Xs$vt, a)
  
  centered_y_hat <- K %*% coef
  fitted_values <- drop(ym +  centered_y_hat)
  resid <- centered_y - centered_y_hat
  GCV <- colSums(resid^2)/(nrow(X) - colSums(matrix(d^2/div, 
                                                    nb_di)))^2
  
  return(list(K = K, sigma = sigma, l = l, 
              coef = drop(coef), scales = x_scaled$xsd,
              ym = ym, xm = x_scaled$xm,
              fitted_values = fitted_values, resid = resid,
              GCV = GCV, scaled_x = X, centered_y = centered_y))
    
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
#' X <- matrix(rnorm(n * p), n, p) # no intercept!
#' y <- rnorm(n)
#' 
#' lams <- 10^seq(-5, 4, length.out = 50)
#' fit_obj <- fit_matern32(x = X, y = y, lambda = lams)
#' 
#' df <- data.frame(predict_matern32(fit_obj, newx = X) - y)
#' summary(df)
#' boxplot(df[, c(1, 10, 25, 35, 50)], 
#' main = "distribution of bias")
#' 
#' 
predict_matern32 <- function(fit_obj, newx, ci = NULL)
{
  if (is.vector(newx))
    newx <- t(newx)
  
  scaled_newx <- matern32::my_scale(x = newx,
                                    xm = as.vector(fit_obj$xm),
                                    xsd = as.vector(fit_obj$scales))
  
  K_star <- matern32_kxstar_cpp(newx = scaled_newx, 
                                x = fit_obj$scaled_x, 
                                sigma = fit_obj$sigma, 
                                l = fit_obj$l)
  
  return(drop(crossprod(K_star, fit_obj$coef)) + fit_obj$ym)
}


#' Find params by minimizing GCV
#'
#' @param x 
#' @param y 
#' @param ... 
#'
#' @return
#' @export
#'
#' @examples
find_params_matern32 <- function(x, y, cl = NULL, ...)
{
  
  OF <- function(xx)
  {
    lams <- 10^seq(-5, 4, length.out = 100)
    GCVs <- fit_matern32(x = x, y = y, 
                         lambda = lams, 
                         sigma = xx[1], 
                         l = xx[2])$GCV
    min(GCVs)
  }
  
  lower_bound <- 10^c(-5, -5)
  upper_bound <- 10^c(4, 4)
  nb_iter <- 100
  
  if (is.null(cl))
  {
    `%op%` <- foreach::`%do%`
    
    pb <- txtProgressBar(min = 0, max = nb_iter, style = 3)
    
    out <- foreach::foreach(i = 1:nb_iter, .errorhandling = "remove")%op%{
      
      set.seed(i*100)
      
      res <- stats::nlminb(start = runif(n = 2, min = lower_bound, 
                                         max = upper_bound), objective = OF, 
                           lower = lower_bound, upper = upper_bound, ...)
      
      setTxtProgressBar(pb, i)
      
      res
    }
    close(pb)
    
  } else {
    
    `%op%` <- foreach::`%dopar%`
    
    cl_SOCK <- parallel::makeCluster(cl, type = "SOCK")
    doSNOW::registerDoSNOW(cl_SOCK)
    
    pb <- txtProgressBar(min = 0, max = nb_iter,
                         style = 3)
    
    progress <- function(n) utils::setTxtProgressBar(pb, n)
    opts <- list(progress = progress)
    
    i <- j <- NULL
    out <- foreach::foreach(i = 1:nb_iter,
                            .packages = "doSNOW",
                            .options.snow = opts,
                            .errorhandling = "remove", ...)%op%{ 
                              
                              set.seed(i*100)
                              
                              stats::nlminb(start = runif(n = 2, min = lower_bound, 
                                                          max = upper_bound), objective = OF, 
                                            lower = lower_bound, upper = upper_bound, ...)
                              
                            }
    
    parallel::stopCluster(cl_SOCK)
  }
  
  index_opt <- which.min(sapply(1:length(out),
                                function (i)
                                  out[[i]]$objective))
  
  return(out[[index_opt]])
  
}