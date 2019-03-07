
# colMeans(der$`2.31012970008316`, na.rm = TRUE)

derivs_matern32 <- function(fit_obj)
{
  n_lambda <- ncol(fit_obj$coef)
  n <- fit_obj$scaled_x
  
  `%op%` <- foreach::`%do%`
  
    res <- foreach::foreach(i = 1:n_lambda)%op%{
      derivs(fit_obj$scaled_x, fit_obj$coef[, i], 
             fit_obj$sigma, fit_obj$l[1])$deriv1
    }
  
  names(res) <- colnames(fit_obj$coef)
  
  res
}