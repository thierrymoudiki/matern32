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
#' matplot(log(lams), t(fit_obj$coef), type = 'l', 
#' main = "coefficients = f(lambda)", xlab = "log(lambda)", ylab = "coefs")
#' abline(h = 0, lty = 2, lwd = 2, col = "red")
#' 
#' 
fit_matern32 <- function(x, y, sigma = 2, l = 0.1, 
                         lambda = 10^seq(-5, 4, length.out = 50), ...)
{
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
  
  return(list(K = K, lambda = lambda, sigma = sigma, l = l, 
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