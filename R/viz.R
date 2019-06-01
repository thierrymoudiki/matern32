



#' Title
#'
#' @param fit_obj 
#' @param verbose 
#'
#' @return
#' @export
#'
#' @examples
plot_residuals <- function(fit_obj, verbose = TRUE)
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
  
  test1 <- Box.test(residuals)
  test2 <- shapiro.test(residuals)
  
    if (verbose)
    {
      print(test1)
      cat("\n")
      print(test2)
    }
  
  return(invisible(list(Box_test = test1, shapiro_test = test2)))
}