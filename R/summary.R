#' Title
#'
#' @param fit_obj 
#' @param obs 
#'
#' @return
#' @export
#'
#' @examples
summary.matern32 <- function(fit_obj, obs=NULL, verbose=TRUE)
{
  if(class(fit_obj) != "matern32")
  {
    warning("`object` must be of class 'matern32'")
    UseMethod("summary")
    return(invisible(NULL))
  }
  
  n <- nrow(fit_obj$scaled_x)
  p <- ncol(fit_obj$scaled_x)
  
  index_lam <- switch(fit_obj$fit_method, 
                      "svd" = which.min(fit_obj$GCV),
                      "eigen" = which.min(fit_obj$loocv),
                      "chol" = which.min(fit_obj$loocv))
  
  if(!is.null(fit_obj$GCV))
  {
    GCV <- fit_obj$GCV[index_lam]
  } else {
    loocv <- fit_obj$loocv[index_lam]
  }
  
  
  if (!is.null(ncol(fit_obj$coef))) # multiple lambdas 
  {
    R_Squared <- fit_obj$R_Squared[index_lam]
    Adj_R_Squared <- fit_obj$Adj_R_Squared[index_lam]
    residuals <- fit_obj$resid[ , index_lam]
    best_lam <- fit_obj$lambda[index_lam]
    
  } else { # one lambda
    
    R_Squared <- fit_obj$R_Squared
    Adj_R_Squared <- fit_obj$Adj_R_Squared
    residuals <- fit_obj$resid
    best_lam <- fit_obj$lambda
  }
  
  summary_resid <- quantile(residuals, 
                            probs = c(0, 25, 50, 75, 100)/100)
  names(summary_resid) <- c("Min.", "1st Qu.",  "Median", 
                            "3rd Qu.", "Max.")
  
  if (is.null(obs))
  {
    derivatives <- matern32::derivatives(fit_obj)[[1]]
  } else {
    stopifnot(obs >= 1 && obs <= n && floor(obs) == obs)
    upper_bound <- n*obs
    lower_bound <- n*obs - n + 1
    derivatives <- matern32::derivatives(fit_obj)[[1]][lower_bound:upper_bound, ]
  }
  
  n_derivatives <- nrow(derivatives)
  estimate <- colMeans(derivatives)
  distro_effects <- summary(derivatives)
  empirical_var <- ((n_derivatives-1)/n_derivatives)*apply(derivatives, 2, var)
  std_error <- sqrt(empirical_var/n_derivatives)
  t_value <- estimate/std_error
  
  get_code_pval <- function(pval)
  {
    
    stopifnot(pval >= 0 && pval <= 1)
    
    if(pval >= 0 && pval < 0.001)
    {
      return ("***")
    }
    if(pval >= 0.001 && pval < 0.01)
    {
      return ("**")
    }
    if(pval >= 0.01 && pval < 0.05)
    {
      return ("*")
    }
    if(pval >= 0.05 && pval < 0.1)
    {
      return (".")
    }
    if(pval >= 0.1)
    {
      return (" ")
    }
  }
  
  p_value <- 2*pt(abs(t_value), 
                  n - p, lower.tail = FALSE)
  
  signif_codes <- sapply(1:p, function (j) get_code_pval(p_value[j]))
  
  if (verbose == TRUE)
  {
    cat("Residuals:", "\n")
    print(drop(summary_resid))
    cat("\n")
    cat("Coefficients:", "\n") 
  }
  coefficients <- cbind.data.frame(estimate, std_error, 
                                   t_value, p_value, 
                                   signif_codes)
  
  colnames(coefficients) <- c("Est", "Std. Error", "t value", 
                              "Pr(>|t|)", "") 
  
  if (verbose == TRUE)
  {
    print(coefficients)
    cat("---", "")
    cat("\n")
    cat("Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1", "\n")
    cat("\n")
    cat("Distribution of marginal effects: ", "\n")
    print(distro_effects)
    cat("\n")
    cat("Multiple R-squared:  ", R_Squared, "	Adjusted R-squared:", Adj_R_Squared,"\n")
  }
  
  if(!is.null(fit_obj$GCV))
  {
    if (verbose == TRUE)
      cat("\n")
      cat("GCV error:  ", GCV, "lambda: ", best_lam, "\n")
    return(invisible(list(coefficients = coefficients, 
                          distro_effects = distro_effects, 
                          R_Squared = as.numeric(R_Squared),
                          Adj_R_Squared = as.numeric(fit_obj$Adj_R_Squared),
                          GCV = as.numeric(GCV), 
                          best_lam = as.numeric(best_lam))))
  } else {
    if (verbose == TRUE)
      cat("\n")
      cat("LOOCV error:  ", loocv, "lambda: ", best_lam, "\n")
    return(invisible(list(coefficients = coefficients, 
                          distro_effects = distro_effects, 
                          R_Squared = as.numeric(R_Squared),
                          Adj_R_Squared = as.numeric(fit_obj$Adj_R_Squared),
                          loocv = as.numeric(loocv), 
                          best_lam = as.numeric(best_lam))))
  }
}