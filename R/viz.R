



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
#'
#' @return
#' @export
#'
#' @examples
plot_coeffs <- function(fit_obj, var = 1)
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
      # cat("i = ", i, "\n")
      # print(d1[index, var])
      # cat("\n")
      
      #res[i, ] <- as.numeric(summary(d1[index, var]))
      res[i, ] <- as.numeric(d1[index, var])
      i <- i + 1
      index <- seq(from = i, to = n^2, by = n)
    }
  
  rownames(res) <- paste0("obs", 1:n)
  #colnames(res) <- c("Min.", "1st Qu.", "Median", "Mean", "3rd Qu.", "Max.")
  
  x_var <- colnames(fit_obj$scaled_x)[var]
  
  matplot(fit_obj$x[ , var], res, type = 'p', pch = 20,
          main = paste("effects of",  x_var, "on response"), 
          xlab = x_var,
          ylab = "delta_1")
  lines(lowess(x = fit_obj$x[ , var], 
               y = apply(res, 1, median)))
}