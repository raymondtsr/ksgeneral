#include <Rcpp.h>
#include "crossprob_new.h"
using namespace Rcpp;

// This is a simple example of exporting a C++ function to R. You can
// source this function into an R session using the Rcpp::sourceCpp
// function (or via the Source button on the editor toolbar). Learn
// more about Rcpp at:
//
//   http://www.rcpp.org/
//   http://adv-r.had.co.nz/Rcpp.html
//   http://gallery.rcpp.org/
//

//// [[Rcpp::export]]
//int up_low_bounds(double n, double x){
//
//  int result = 0;
//  upperlowerbounds(n, x);
//
//  return result;
//}

// [[Rcpp::export]]
double cont_ks_distribution_Rcpp(double n){
    
//    upperlowerbounds(n, x);
    
    return cont_ks_distribution(n);
    
//    return 0;
    
}

// [[Rcpp::export]]
double cont_ks_distribution_Rcpp_alternative(double n, double x){
    
    upperlowerbounds_alternative(n, x);
    
    return cont_ks_distribution(n);
    
    //    return 0;
    
}


//// [[Rcpp::export]]
//double disc_ks_distribution_Rcpp_alternative(double n, NumericVector lbs, NumericVector ubs){
//    return disc_ks_distribution(n, lbs, ubs);
//}

// You can include R code blocks in C++ files processed with sourceCpp
// (useful for testing and development). The R code will be automatically
// run after the compilation.
//

/*** R
up_low_bounds(10, 0.1)
*/
