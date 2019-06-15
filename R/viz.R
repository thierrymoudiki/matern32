

#' Title
#'
#' @param fit_obj 
#'
#' @return
#' @export
#'
#' @examples
plot_coeffs <- function(fit_obj)
{
  stopifnot(length(fit_obj$lambda) > 1)
  matplot(log(fit_obj$lambda), t(fit_obj$coef), type = 'l', 
  main = "coefficients = f(lambda)", xlab = "log(lambda)", 
  ylab = "coefs")
  abline(h = 0, lty = 2, lwd = 2, col = "red")
  
  invisible(list(lambda = fit_obj$lambda, 
                 coef = fit_obj$coef))
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
  plot(log(fit_obj$lambda),  fit_obj$GCV, type = 'l', 
       main = "GCV error", 
       xlab = "log(lambda)",
       ylab = "GCV")
  
  invisible(list(lambda=fit_obj$lambda, 
                 GCV=fit_obj$GCV))
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
plot_heterogeneity <- function(fit_obj, var = 1, ...)
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


#' Title
#'
#' @param fit_obj 
#' @param var1 
#' @param var2 
#' @param ... 
#'
#' @return
#' @export
#'
#' @examples
plot_interactions <- function(fit_obj, var1 = 1, var2 = 2, 
                              degree = 1, ...)
{
  stopifnot(var1 != var2)
  col_names <- colnames(fit_obj$scaled_x)
  n <- nrow(fit_obj$scaled_x)
  
  if (!is.null(dim(fit_obj$coef))) # matrix of coeffs 
  {
    i_best <- switch(fit_obj$fit_method, 
                     "svd" = which.min(fit_obj$GCV),
                     "eigen" = which.min(fit_obj$loocv),
                     "chol" = which.min(fit_obj$loocv))
    
    res <- inters2(x = fit_obj$scaled_x, 
                   j1 = var1, j2 = var2, 
                   c = fit_obj$coef[, i_best], 
                   l = fit_obj$l[1])
    
  } else { # 1 vector
    
    res <- inters2(x = fit_obj$scaled_x, 
                   j1 = var1, j2 = var2, 
                   c = fit_obj$coef, 
                   l = fit_obj$l[1])
  }
  
  # smooth the interactions
  x_var1 <- fit_obj$x[, var1]
  x_var2 <- fit_obj$x[, var2]
  x_var1_rep <- rep(x_var1, n)
  x_var2_rep <- rep(x_var2, n)
  grid <- cbind.data.frame(x = x_var1_rep, 
                           y = x_var2_rep)
  grid$z <- res 
  index_na <- !is.na(res)
  
  # loess_obj <- loess(z ~ ., data = grid, degree = 2)
  if (degree == 3)
  {
    smooth_obj <- matern32::fit_poly(X = cbind(x_var1_rep[index_na], 
                                                x_var2_rep[index_na]), 
                                      y = grid$z[index_na], 
                                     degree = degree) 
  } 
  
  if (degree == 2)
  {
    smooth_obj <- matern32::fit_poly(X = cbind(x_var1_rep[index_na], 
                                                x_var2_rep[index_na]), 
                                      y = grid$z[index_na], 
                                     degree = degree) 
  }
  
  if (degree == 1)
  {
    smooth_obj <- matern32::fit_poly(X = cbind(x_var1_rep[index_na], 
                                                x_var2_rep[index_na]), 
                                      y = grid$z[index_na], 
                                     degree = degree) 
  }
  
  # predict smoothed interactions on a grid
  min_x <- min(x_var1); max_x <- max(x_var1)
  min_y <- min(x_var2); max_y <- max(x_var2)
  n_out <- 10
  x_new <- seq(min_x, max_x, length.out = n_out)
  y_new <- seq(min_y, max_y, length.out = n_out)
  
  new_grid <- data.frame(x = rep(x_new, each = n_out), 
                         y = rep(y_new, n_out))
  
  preds <- matern32::predict_poly(object = smooth_obj, 
                                     newx = new_grid)
  
  # output matrix of smoothed interactions
  z <- matrix(preds, 
              nrow = n_out, ncol = n_out, 
              byrow = TRUE)
  colnames(z) <- round(y_new, 2)
  rownames(z) <- round(x_new, 2)
  
  # plot stuff
  if (!is.null(col_names))
  {
    filled.contour(x = x_new, y = y_new, z = z, 
                   color = terrain.colors, 
                   main = paste("smoothed interactions effects \n", 
                                col_names[var1], "x", col_names[var2]),
                   xlab = col_names[var1], ylab = col_names[var2]
    )
  } else {
    filled.contour(x = x_new, y = y_new, z = z, 
                   color = terrain.colors, 
                   main = paste("smoothed interactions effects \n", 
                                "covariate", var1, "x covariate", var2),
                   xlab = paste("covariate", var1), 
                   ylab = paste("covariate", var2))
  }
  
  # return smoothed grid
  invisible(z)
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

#' Title
#'
#' @param fit_obj 
#' @param var 
#'
#' @return
#' @export
#'
#' @examples
# plot_heterogen2 <- function(fit_obj, var = 1)
# {
#   n <- nrow(fit_obj$scaled_x)
#   p <- ncol(fit_obj$scaled_x)
#   i <- 2
#   index <- seq(from = i, to = n^2, by = n)
#   
#   # result for each observation
#   #res <- matrix(0, ncol = 6, nrow = n)
#   res <- matrix(0, ncol = n, nrow = n)
#   
#   if (!is.null(dim(fit_obj$coef)))
#   {
#     i_best <- switch(fit_obj$fit_method, 
#                      "svd" = which.min(fit_obj$GCV),
#                      "eigen" = which.min(fit_obj$loocv),
#                      "chol" = which.min(fit_obj$loocv))
#     
#     fitted_values <- fit_obj$fitted_values[,i_best]
#     residuals <- fit_obj$resid[,i_best]
#   } else {
#     fitted_values <- fit_obj$fitted_values
#     residuals <- fit_obj$resid
#   }
#   
#   d1 <- derivatives(fit_obj)$`1D`
#   
#   while (i <= n)
#   {
#     res[i, ] <- as.numeric(d1[index, var])
#     i <- i + 1
#     index <- seq(from = i, to = n^2, by = n)
#   }
#   
#   x_var <- colnames(fit_obj$scaled_x)[var]
#   
#   matplot(fit_obj$x[ , var], res, type = 'p', pch = 20,
#           main = paste("1st order effects of",  x_var, "\n on response"), 
#           xlab = x_var,
#           ylab = "change for an increase of 1")
#   lines(lowess(x = fit_obj$x[ , var], 
#                y = apply(res, 1, median)))
# }