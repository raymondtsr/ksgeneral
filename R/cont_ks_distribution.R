cont_ks_cdf <- function(q, n)
{

  PVAL <- NULL

  if ((n <= 140) && (n*(q)^2 > 12)){
    PVAL <- 1
  }
  else if ((n > 140) && (n <= 100000) && (n*(q)^(3/2) >= 1.4) && (n*(q)^2 >= 11)){
    PVAL <- 1
  }
  else if (n*(q)^2 >= 18){
    PVAL <- 1
  }
  else{
    PVAL <- 1 - ksgeneral::cont_ks_distribution_Rcpp_alternative(n, q)
  }

  PVAL <- min(1.0, max(0.0, PVAL))

  return(PVAL)

}
####################################################
cont_ks_c_cdf <- function(q, n)
{

  PVAL <- 1 - cont_ks_cdf(q, n)


  return(PVAL)

}
