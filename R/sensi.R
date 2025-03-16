
#' Title
#'
#' @param fit_obj object fitted by using \code{fit_matern32}
#' @param index_col an integer, column index 
#' @param h (must be small)
#'
#' @return
#' @export
#'
#' @examples
sensi1D <- function(fit_obj, index_col, h, verbose=TRUE)
{
  
  if (fit_obj$with_kmeans)
  {
    scaled_x <- fit_obj$scaled_x_clust
  } else {
    scaled_x <- fit_obj$scaled_x
  }
  
  n <- nrow(scaled_x)
  
  if(n > 500)
    cat("Processing...", "\n")
  
  # d1 & d2 
  if (!is.null(dim(fit_obj$coef)))
  {
    i_best <- switch(fit_obj$fit_method, 
                     "svd" = which.min(fit_obj$GCV),
                     "eigen" = which.min(fit_obj$loocv),
                     "chol" = which.min(fit_obj$loocv), 
                     "solve" = which.min(fit_obj$loocv))
    
    ans_derivs <- derivs(x = scaled_x, 
                         c = fit_obj$coef[, i_best], 
                         l = fit_obj$l[1])
  } else {
    ans_derivs <- derivs(x = scaled_x, 
                         c = fit_obj$coef, 
                         l = fit_obj$l[1])
  }
  
  res <- ans_derivs$deriv1[, index_col]*h + 0.5*ans_derivs$deriv2[, index_col]*(h^2)
  names(res) <- paste(rep(1:n, each = n),
                      rep(1:n), sep = ".")
  
  
  if (verbose)
  {
    col_names <- colnames(scaled_x)
    if (!is.null(col_names))
      cat("Heterogeneity of changes in the response for", 
          col_names[index_col], 
          ifelse(h > 0, "+", "-"), abs(h),": \n")
    print(summary(res[!is.na(res)]))
  }
  
  return(invisible(res))
}
#sensi1D <- memoise::memoize(sensi1D)

#' Title
#'
#' @param fit_obj object fitted by using \code{fit_matern32}
#' @param index_col1 an integer, column index 
#' @param index_col2 an integer, second column index 
#' @param h (must be small)
#'
#' @return
#' @export
#'
#' @examples
sensi2D <- function(fit_obj, index_col1, index_col2, h1, h2, verbose=TRUE)
{
  
  if (fit_obj$with_kmeans)
  {
    scaled_x <- fit_obj$scaled_x_clust
  } else {
    scaled_x <- fit_obj$scaled_x
  }
  
  n <- nrow(scaled_x)
  
  if(n > 500)
    cat("Processing...", "\n")
  
  # d1 & d2 
  if (!is.null(dim(fit_obj$coef)))
  {
    i_best <- switch(fit_obj$fit_method, 
                     "svd" = which.min(fit_obj$GCV),
                     "eigen" = which.min(fit_obj$loocv),
                     "chol" = which.min(fit_obj$loocv), 
                     "solve" = which.min(fit_obj$loocv))
    
    # d1 & d2
    ans_derivs <- derivs(x = scaled_x, 
                         c = fit_obj$coef[, i_best], 
                         l = fit_obj$l[1])
    
    # inters
    term_interactions <- inters2(scaled_x, j1 = index_col1, j2 = index_col2, 
                            c = fit_obj$coef[, i_best], l = fit_obj$l[1])
    
  } else {
    
    # d1 & d2
    ans_derivs <- derivs(x = scaled_x, 
                  c = fit_obj$coef, 
                  l = fit_obj$l[1])
    
    # inters
    term_interactions <- inters2(scaled_x, j1 = index_col1, j2 = index_col2, 
                            c = fit_obj$coef, l = fit_obj$l[1])
  }
  
  term1 <- ans_derivs$deriv1[, index_col1]*h1 + 0.5*ans_derivs$deriv2[, index_col1]*(h1^2)
  term2 <- ans_derivs$deriv1[, index_col2]*h2 + 0.5*ans_derivs$deriv2[, index_col2]*(h2^2)
  res <-  term1 + term2 + term_interactions
  names(res) <- paste(rep(1:n, each = n),
                      rep(1:n), sep = ".")
  
  if (verbose)
  {
    col_names <- colnames(scaled_x)
    if (!is.null(col_names))
      cat("Heterogeneity of changes in the response for \n", col_names[index_col1], 
          ifelse(h1 > 0, "+", "-"), abs(h1), " and ", col_names[index_col2], 
          ifelse(h2 > 0, "+", "-"), abs(h2), ": \n")
    print(summary(res[!is.na(res)])) 
  }
  
  return(invisible(res))
}
#sensi2D <- memoise::memoize(sensi2D)