

#' Title
#'
#' @param fit_obj 
#'
#' @return
#' @export
#'
#' @examples
derivs_matern32 <- function(fit_obj)
{
  n <- nrow(fit_obj$scaled_x)
  
  if(n > 500)
    cat("Processing...", "\n")
  
    if (!is.null(dim(fit_obj$coef)))
    {
      i_best <- switch(fit_obj$fit_method, 
                       "svd" = which.min(fit_obj$GCV),
                       "eigen" = which.min(fit_obj$loocv),
                       "chol" = which.min(fit_obj$loocv))
      
      res <- derivs(x = fit_obj$scaled_x, 
                    c = fit_obj$coef[, i_best], 
                    l = fit_obj$l[1])
    } else {
      res <- derivs(x = fit_obj$scaled_x, 
                    c = fit_obj$coef, 
                    l = fit_obj$l[1])
    }
  
  col_names_x <- colnames(fit_obj$scaled_x)
  
  if (!is.null(col_names_x))  
    colnames(res[[1]]) <- colnames(res[[2]]) <- col_names_x

  rownames(res[[1]]) <- rownames(res[[2]]) <- paste(rep(1:n, each = n),
                                                    rep(1:n), sep = ".")
  
  return(res)
}


#' Title
#'
#' @param fit_obj 
#' @param index_col1 
#' @param index_col2 
#'
#' @return
#' @export
#'
#' @examples
inters_matern32 <- function(fit_obj, index_col1, index_col2)
{
  n <- nrow(fit_obj$scaled_x)
  p <- ncol(fit_obj$scaled_x)
  l <- sqrt(p)
  
    if (!is.null(dim(fit_obj$coef)))
    {
      i_best <- switch(fit_obj$fit_method, 
                       "svd" = which.min(fit_obj$GCV),
                       "eigen" = which.min(fit_obj$loocv),
                       "chol" = which.min(fit_obj$loocv))
      
      res_inters <- inters(x = fit_obj$scaled_x, j1 = index_col1,
                           j2 = index_col2, c = fit_obj$coef[ , i_best], 
                           l = l)
    } else {
      res_inters <- inters(x = fit_obj$scaled_x, j1 = index_col1,
                           j2 = index_col2, c = fit_obj$coef, 
                           l = l)
    }
  
  res <- apply(res_inters, 1, summary)[-7,]
  colnames(res) <- paste0("obs", 1:n) 
  
  if (!is.null(colnames(fit_obj$scaled_x)))
  {
    col_names <- colnames(fit_obj$scaled_x)
    cat("Interaction effects between", col_names[index_col1]," and ", col_names[index_col2], ":\n") 
  }
  
  return(res)
}