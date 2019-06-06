

#' Title
#'
#' @param fit_obj 
#'
#' @return
#' @export
#'
#' @examples
plot_coeffs1 <- function(fit_obj)
{
  stopifnot(length(fit_obj$lambda) > 1)
  matplot(log(lams), t(fit_obj$coef), type = 'l', 
  main = "coefficients = f(lambda)", xlab = "log(lambda)", 
  ylab = "coefs")
  abline(h = 0, lty = 2, lwd = 2, col = "red")
}

#' Title
#'
#' @param fit_obj 
#'
#' @return
#' @export
#'
#' @examples
plot_GCV <- function(fit_obj){
  plot(fit_obj$GCV, type = 'l')
}


#' Title
#'
#' @param fit_obj 
#' @param var 
#' @param ... 
#'
#' @return
#' @export
#'
#' @examples
plot_heterogen1 <- function(fit_obj, var = 1, ...)
{
  d1 <- derivatives(fit_obj)$`1D`
  col_names <- colnames(fit_obj$scaled_x)
  
  if (!is.null(col_names))
  {
    h <- hist(d1[, var], 
              main = paste("Heterogeneity of \n 1st order effects for ", col_names[var]), 
              xlab = "1st order effect", ...) 
  } else {
    h <- hist(d1[, var], 
              main = paste("Heterogeneity of \n 1st order effects for covariate", var), 
              xlab = "1st order effect", ...)
  }
    
  return(list(breaks = h$breaks, 
              counts = h$count))
}


plot_interactions1 <- function(fit_obj, var1 = 1, var2 = 2, ...)
{
  stopifnot(var1 != var2)
  col_names <- colnames(fit_obj$scaled_x)
  n <- nrow(fit_obj$scaled_x)
  
  if (!is.null(dim(fit_obj$coef)))
  {
    i_best <- switch(fit_obj$fit_method, 
                     "svd" = which.min(fit_obj$GCV),
                     "eigen" = which.min(fit_obj$loocv),
                     "chol" = which.min(fit_obj$loocv))
    
    res <- inters2(x = fit_obj$scaled_x, j1 = var1, j2 = var2, 
                  c = fit_obj$coef[, i_best], l = fit_obj$l[1])
    
  } else {
    
    res <- inters2(x = fit_obj$scaled_x, j1 = var1, j2 = var2, 
                   c = fit_obj$coef, l = fit_obj$l[1])
  }
  
  grid <- cbind.data.frame(x = rep(fit_obj$x[, var1], n), 
                           y = rep(fit_obj$x[, var2], n))
  grid$z <- res 
  
  grid <- grid[!is.na(grid$z), ]
  
  if (!is.null(col_names))
  {
    levelplot(z~x*y, data = grid, col.regions = terrain.colors,
              main = paste("Interactions of \n", col_names[var1], 
                                     "and", col_names[var2]), 
              xlab = col_names[var1], ylab = col_names[var2], ...)
  } else {
    levelplot(z~x*y, data = grid, col.regions = terrain.colors,
              main = paste("Interactions of covariate \n", var1, 
                           "and covariate", var2), 
              xlab = paste("covariate", var1), 
              ylab =  paste("covariate", var2), ...)
  }

  #return(grid)
}



#' Title
#'
#' @param fit_obj 
#' @param var 
#'
#' @return
#' @export
#'
#' @examples
plot_heterogen2 <- function(fit_obj, var = 1)
{
  n <- nrow(fit_obj$scaled_x)
  p <- ncol(fit_obj$scaled_x)
  i <- 2
  index <- seq(from = i, to = n^2, by = n)
  
  # result for each observation
  #res <- matrix(0, ncol = 6, nrow = n)
  res <- matrix(0, ncol = n, nrow = n)
  
  if (!is.null(dim(fit_obj$coef)))
  {
    i_best <- switch(fit_obj$fit_method, 
                     "svd" = which.min(fit_obj$GCV),
                     "eigen" = which.min(fit_obj$loocv),
                     "chol" = which.min(fit_obj$loocv))
    
    fitted_values <- fit_obj$fitted_values[,i_best]
    residuals <- fit_obj$resid[,i_best]
  } else {
    fitted_values <- fit_obj$fitted_values
    residuals <- fit_obj$resid
  }
  
  d1 <- derivatives(fit_obj)$`1D`
  
  # cpp loop
  while (i <= n)
  {
    res[i, ] <- as.numeric(d1[index, var])
    i <- i + 1
    index <- seq(from = i, to = n^2, by = n)
  }
  
  x_var <- colnames(fit_obj$scaled_x)[var]
  
  matplot(fit_obj$x[ , var], res, type = 'p', pch = 20,
          main = paste("1st order effects of",  x_var, "\n on response"), 
          xlab = x_var,
          ylab = "change for an increase of 1")
  lines(lowess(x = fit_obj$x[ , var], 
               y = apply(res, 1, median)))
}




#' Plot residuals
#'
#' @param fit_obj 
#'
#' @return
#' @export
#'
#' @examples
plot_residuals <- function(fit_obj)
{
  
  if (!is.null(dim(fit_obj$coef)))
  {
    i_best <- switch(fit_obj$fit_method, 
                     "svd" = which.min(fit_obj$GCV),
                     "eigen" = which.min(fit_obj$loocv),
                     "chol" = which.min(fit_obj$loocv))
    
    fitted_values <- fit_obj$fitted_values[,i_best]
    residuals <- fit_obj$resid[,i_best]
    
  } else {
    
    fitted_values <- fit_obj$fitted_values
    residuals <- fit_obj$resid
    
  }
  
  plot(x = fitted_values, 
       y = residuals, type = 'p', 
       col = "gray60", main = "residuals vs fitted values",
       xlab = "fitted values", ylab = "residuals")
  abline(h = 0, col = "red")
  
}