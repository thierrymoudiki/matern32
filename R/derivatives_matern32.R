

#' Title
#'
#' @param fit_obj 
#'
#' @return
#' @export
#'
#' @examples
derivatives <- function(fit_obj, obs = NULL)
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
  
  if (is.null(obs))
  {
    col_names_x <- colnames(fit_obj$scaled_x)
    
    if (!is.null(col_names_x))  
      colnames(res[[1]]) <- colnames(res[[2]]) <- col_names_x
  
    rownames(res[[1]]) <- rownames(res[[2]]) <- paste(rep(1:n, each = n),
                                                      rep(1:n), sep = ".")
    
    names(res) <- c("1D", "2D")
    
    return(res)
    
  } else {
    
    stopifnot(obs >= 1 && obs <= n && floor(obs) == obs)
    upper_bound <- n*obs
    lower_bound <- n*obs - n + 1
    res <- list(res[[1]][lower_bound:upper_bound, ], 
                res[[2]][lower_bound:upper_bound, ])
    
    col_names_x <- colnames(fit_obj$scaled_x)
    
    if (!is.null(col_names_x))  
      colnames(res[[1]]) <- colnames(res[[2]]) <- col_names_x
    
    rownames(res[[1]]) <- rownames(res[[2]]) <- seq(1, n, by = 1)
    
    names(res) <- c("1D", "2D")
    
    return(res)
    
  }
  
}


#' Title
#'
#' @param fit_obj 
#' @param index_col1 
#' @param index_col2 
#' @param obs 
#'
#' @return
#' @export
#'
#' @examples
interactions <- function(fit_obj, index_col1, index_col2, obs = NULL)
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
      
      res_inters <- inters(x = as.matrix(fit_obj$scaled_x), j1 = index_col1,
                           j2 = index_col2, c = fit_obj$coef[ , i_best], 
                           l = l)
    } else {
      
      res_inters <- inters(x = as.matrix(fit_obj$scaled_x), j1 = index_col1,
                           j2 = index_col2, c = fit_obj$coef, 
                           l = l)
    }
    colnames(res_inters) <- paste(rep(1, n),
                                  seq(1, n), sep = ".")
    rownames(res_inters) <- paste(seq(1, n),
                                  rep(1, n), sep = ".")
  
    if (is.null(obs)) # for all the observations 
    {
      if (!is.null(colnames(fit_obj$scaled_x)))
      {
        col_names <- colnames(fit_obj$scaled_x)
        cat("Interaction effects between", col_names[index_col1],
            " and ", col_names[index_col2], ":\n") 
      }
      
      print(summary(as.vector(res_inters)))
      return(invisible(res_inters)) 
      
    } else { # for one observation
      
      stopifnot(obs >= 1 && obs <= n && floor(obs) == obs)
      if (!is.null(colnames(fit_obj$scaled_x)))
      {
        col_names <- colnames(fit_obj$scaled_x)
        cat("Interaction effects between", col_names[index_col1],
            " and ", col_names[index_col2], "for observation #", obs, ":\n") 
      }
      
      print(summary(as.vector(res_inters[obs, ])))
      return(invisible(res_inters[obs, ])) 
    }
  
}