

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

# derivs_matern32R <- function(fit_obj)
# {
#   X <- fit_obj$scaled_x
#   i_best <- which.min(fit_obj$GCV)
#   c <- fit_obj$coef[ , i_best]
#   temp <- sqrt(3)/fit_obj$l[1]
#   const_mult <- temp^2
#   
#   n <- nrow(X)
#   p <- ncol(X)
#   deriv1_mat <- deriv2_mat <- matrix(0, nrow = n^2, ncol = p)
#   
#   i <- 1
#   
#     for (i0 in 1:n)
#     {
#       for (k in 1:n)
#       {
#         vec <- X[i0,] - X[k,]
#         r <- sqrt(sum(vec^2))
#         temp2 <- c[k]*exp(-temp*r)
#         deriv1_mat[i, ] <- temp2*as.numeric(vec)
#         #cat("The value of deriv1: ", deriv1_mat[i0*k, 1], "\n")
#         deriv2_mat[i, ] <- temp2*((temp/r)*(vec^2) - 1)
#         i <- i + 1
#       }
#     }
#   
#   return(list(deriv1 = const_mult*deriv1_mat, 
#               deriv2 = const_mult*deriv2_mat))
# }