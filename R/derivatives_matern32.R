

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
    i_best <- which.min(fit_obj$GCV)
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