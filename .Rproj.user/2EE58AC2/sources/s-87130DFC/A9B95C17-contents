#include <Rcpp.h>
#include <Math.h>
using namespace Rcpp;


/* 0 - utils */

// [[Rcpp::export]]
double weighted_l2_norm(NumericVector x, NumericVector l)
{
  unsigned long int n = x.size();
  if (l.size() != n) {
    ::Rf_error("you must have x.size() == l.size()");
  }
  double res = 0;
  for(int i = 0; i < n; i++) {
    res += pow(x(i), 2)/pow(l(i), 2);
  }
  return(sqrt(res));
}


/* 1 - Matérn 3/2 kernel */

// [[Rcpp::export]]
NumericMatrix matern32_kxx_cpp(NumericMatrix x, double sigma, 
                               NumericVector l)
{
  unsigned long int n = x.nrow();
  NumericMatrix res(n, n);
  double sqrt3 = sqrt(3);
  double temp = 0;
  
  for(int i = 0; i < n; i++) {
    for(int j = i; j < n; j++) {
      temp = sqrt3*weighted_l2_norm(x(i, _) - x(j, _), l);
      res(i , j) = (1 + temp)*exp(-temp);
      res(j , i) = res(i , j);
    }
  }
  
  return(pow(sigma, 2)*res);
}


// [[Rcpp::export]]
NumericMatrix matern32_kxstar_cpp(NumericMatrix newx, NumericMatrix x,
                                  double sigma, NumericVector l)
{
  unsigned long int m = newx.nrow();
  unsigned long int n = x.nrow();
  NumericMatrix res(m, n);
  double temp = 0;
  double sqrt3 = sqrt(3);
  
  for(int i = 0; i < m; i++) {
    for(int j = 0; j < n; j++) {
      temp = sqrt3*sqrt(weighted_l2_norm(newx(i, _) - x(j, _), l));
      res(i , j) = (1 + temp)*exp(-temp);
    }
  }
  
  return(pow(sigma, 2)*res);
}

// You can include R code blocks in C++ files processed with sourceCpp
// (useful for testing and development). The R code will be automatically 
// run after the compilation.
//

/*** R
n <- 7 ; p <- 2
X <- matrix(rnorm(n * p), n, p) # no intercept!
y <- rnorm(n)

matern32_kxx_cpp(X, sigma = 0.1, l = rep(0.1, p))
*/
