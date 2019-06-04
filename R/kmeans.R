
# # 0 - datasets ------------------------------------------------------------
# 
# # for a  large n only!
# n <- 1000 ; p <- 5
# 
# # dataset 1
#  set.seed(7116342)
#  X <- matrix(rlnorm(n * p), n, p) # no intercept!
# 
# # dataset 2
# #X <- as.matrix(Boston)
# 
# # 1 - BIC ------------------------------------------------------------
# 
# kmeansBIC <- function(X, k){
#   set.seed(123)
#   fit_kmeans <- stats::kmeans(X, centers = k)
#   fit_kmeans$tot.withinss + log(length(fit_kmeans$cluster))*ncol(fit_kmeans$centers)*k
# }
# 
# OF <- function(k)
# {
#   kmeansBIC(X, as.integer(k))
# }
# 
# nb_k <- min(50, nrow(X)-1)
# res <- rep(0, nb_k)
# pb <- txtProgressBar(min = 1, max = nb_k, style = 3)
# for (i in 1:nb_k){
#   res[i] <- OF(i)
#   setTxtProgressBar(pb, i)
#   close(pb)
# }
# 
# plot(res, type = 'l', main = "BIC = \n f(number of clusters)",
#      ylab = "BIC", xlab = "number of clusters")
# 
# # 2 - find k ------------------------------------------------------------
# 
# # do this in cpp
# n <- length(res)
# n_obs <- 10
# alpha <- 0.05
# rel_tol <- diff(res)/res[-n]
# check_val <- FALSE
# i <- 1
# ub <- (i + n_obs- 1)
#   while(check_val == FALSE && (ub <= length(rel_tol)))
#   {
#     rel_tols <- rel_tol[i:ub]
#     check_val <- (sum(abs(rel_tols) <= alpha) == n_obs)
#     i <- i + 1
#     ub <- (i + n_obs- 1)
#   }
# print(i)
# 
# abline(v = i, lty = 2, col = "red")
