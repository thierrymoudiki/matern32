

#' Title
#'
#' @param fit_obj 
#' @param index_col 
#' @param h 
#'
#' @return
#' @export
#'
#' @examples
sensi1D <- function(fit_obj, index_col, h)
{
  n <- nrow(fit_obj$scaled_x)
  
  if(n > 500)
    cat("Processing...", "\n")
  
  # d1 & d2 
  if (!is.null(dim(fit_obj$coef)))
  {
    i_best <- switch(fit_obj$fit_method, 
                     "svd" = which.min(fit_obj$GCV),
                     "eigen" = which.min(fit_obj$loocv),
                     "chol" = which.min(fit_obj$loocv))
    
    ans_derivs <- derivs(x = fit_obj$scaled_x, 
                  c = fit_obj$coef[, i_best], 
                  l = fit_obj$l[1])
  } else {
    ans_derivs <- derivs(x = fit_obj$scaled_x, 
                  c = fit_obj$coef, 
                  l = fit_obj$l[1])
  }
  
  res <- ans_derivs$deriv1[, index_col]*h + 0.5*ans_derivs$deriv2[, index_col]*(h^2)
  names(res) <- paste(rep(1:n, each = n),
                      rep(1:n), sep = ".")
  
  print(summary(res[!is.na(res)]))
  
  return(invisible(res))
}


#' Title
#'
#' @param fit_obj 
#' @param index_col1 
#' @param index_col2 
#' @param h 
#'
#' @return
#' @export
#'
#' @examples
sensi2D <- function(fit_obj, index_col1, index_col2, h1, h2)
{
  n <- nrow(fit_obj$scaled_x)
  
  if(n > 500)
    cat("Processing...", "\n")
  
  # d1 & d2 
  if (!is.null(dim(fit_obj$coef)))
  {
    i_best <- switch(fit_obj$fit_method, 
                     "svd" = which.min(fit_obj$GCV),
                     "eigen" = which.min(fit_obj$loocv),
                     "chol" = which.min(fit_obj$loocv))
    
    # d1 & d2
    ans_derivs <- derivs(x = fit_obj$scaled_x, 
                  c = fit_obj$coef[, i_best], 
                  l = fit_obj$l[1])
    
    # inters
    interactions <- inters2(fit_obj$scaled_x, j1 = index_col1, j2 = index_col2, 
                            c = fit_obj$coef[, i_best], l = fit_obj$l[1])
    
  } else {
    
    # d1 & d2
    ans_derivs <- derivs(x = fit_obj$scaled_x, 
                  c = fit_obj$coef, 
                  l = fit_obj$l[1])
    
    # inters
    interactions <- inters2(fit_obj$scaled_x, j1 = index_col1, j2 = index_col2, 
                            c = fit_obj$coef, l = fit_obj$l[1])
  }
  
  series1 <- ans_derivs$deriv1[, index_col1]*h1 + 0.5*ans_derivs$deriv2[, index_col1]*(h1^2)
  series2 <- ans_derivs$deriv1[, index_col2]*h2 + 0.5*ans_derivs$deriv2[, index_col2]*(h2^2)
  res <-  series1 + series2 + interactions
  names(res) <- paste(rep(1:n, each = n),
                      rep(1:n), sep = ".")
  
  print(summary(res[!is.na(res)]))
  
  return(invisible(res))
}