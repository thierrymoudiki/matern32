#include <Rcpp.h>
#include <Math.h>
using namespace Rcpp;


/* 0 - utils */

// [[Rcpp::export]]
double l2_norm(NumericVector x)
{
  unsigned long int n = x.size();
  
  double res = 0;
  for(int i = 0; i < n; i++) {
    res += pow(x(i), 2);
  }
  return(sqrt(res));
}

// [[Rcpp::export]]
NumericMatrix na_matrix(unsigned int n, unsigned int p){
  NumericMatrix m(n,p) ;
  std::fill( m.begin(), m.end(), NumericVector::get_na() ) ;
  return m ;
}

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

// [[Rcpp::export]]
List derivs(NumericMatrix x, 
            NumericVector c, 
            double sigma, double l)
{
  
  unsigned long int n = x.nrow();
  unsigned long int p = x.ncol();
  
  if (c.size() != n) {
    ::Rf_error("you must have c.size() == x.nrow()");
  }
  
  double r;
  double n2 = pow(n, 2);
  NumericVector vec(n);
  NumericMatrix res = na_matrix(n2, p);
  NumericMatrix res2 = na_matrix(n2, p);
  double temp = sqrt(3)/l;
  double temp2 = 0;
  double const_mult = pow(sigma, 2)*pow(temp, 2);
  
  for(unsigned long int i0 = 0; i0 < n; i0++){
    for(unsigned long int k = 0; k < n; k++){ // There is something to optimize here
      vec = x(i0, _) - x(k, _);
      r = l2_norm(vec);
      temp2 = c(k)*exp(-temp*r); // temp = sqrt(3)/l;
      res(i0*k, _) = temp2*vec; // first derivative
      res2(i0*k, _) = temp2*((temp/r)*pow(vec, 2) - 1); // temp = sqrt(3)/l; // second derivative
    }
  }
  
  return List::create(Rcpp::Named("deriv1") = const_mult*res,
                      Rcpp::Named("deriv2") = const_mult*res2);
}


// [[Rcpp::export]]
NumericVector inters(NumericMatrix x, 
                     NumericVector c, 
                     unsigned long int i0, 
                     double sigma, double l)
{
  
  unsigned long int n = x.nrow();
  unsigned long int p = x.ncol();
  
  if (c.size() != n) {
    ::Rf_error("you must have c.size() == x.nrow()");
  }
  
  if (i0 >= n) {
    ::Rf_error("you must have i0 < n");
  }
  
  double r;
  NumericVector vec(n);
  NumericMatrix res(n, p), res2(n, p);
  double temp = sqrt(3)/l;
  double temp2 = 0;
  double const_mult = pow(sigma, 2)*pow(temp, 2);
  
  for(unsigned long int k = 0; k < n; k++){
    vec = x(i0, _) - x(k, _);
    r = l2_norm(vec);
    temp2 = c(k)*exp(-temp*r);
    res2(k, _) = temp2*((temp/r)*pow(vec, 2) - 1); // second derivative
  }
  
  return (const_mult*res2);
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
