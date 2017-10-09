# Mixed_cdf <- function(x)
# {
#   result <- 0
#   if (x < 0){
#     result <- 0
#   }
#   else if (x == 0){
#     result <- 0.5
#   }
#   else if (x < log(2.5)){
#     result <- 1 - 0.5 * exp(-x)
#   }
#   else{
#     result <- 1
#   }
#
#   return (result)
# }

#####################################################

# fun1 <- function(x)
# {
#   result <- 0
#   if (x < 0){
#     result <- 0
#   }
#   else if (x == 0){
#     result <- 0.2
#   }
#   else if (x < 0.2){
#     result <- x + 0.2
#   }
#   else if (x == 0.2){
#     result <- 0.6
#   }
#   else if (x < 0.4){
#     result <- x + 0.4
#   }
#   else{
#     result <- 1
#   }
#
#   return (result)
# }


##############################################
# Vec_Mixed_cdf_1 <- function(x)
# {
#   return (sapply(x,Mixed_cdf))
# }
##############################################
mixed_ks_test <- function(x, jump_points, Mixed_dist, ..., tol = 1e-10){

  Vec_Mixed_cdf <- function(x)
  {
    return (sapply(x,Mixed_dist))
  }

  DNAME <- deparse(substitute(x))
  x <- x[!is.na(x)]
  n <- length(x)
  if (n < 1L){
    stop("not enough 'x' data")
  }
  PVAL <- NULL

  METHOD <- "One-sample Kolmogorov-Smirnov test"

  dev_up <- max(abs(c(0, ecdf(x)(x) - Vec_Mixed_cdf(x))))
  dev_low <- max(abs(c(0, ecdf(x)(x-tol) - Vec_Mixed_cdf(x-tol))))
  STATISTIC <- max(dev_up, dev_low)
  ############################################
  a <- rep(NA,n)
  f_a <- a
  b <- a
  f_b <- a

  temp_a <- a
  temp_b <- a

  temp_length <- length(jump_points) - 1

  temp <- rep(0, (2*temp_length))

  for (i in 1:temp_length){
    temp[2*i-1] <- Vec_Mixed_cdf(jump_points[i])
    temp[2*i] <- Vec_Mixed_cdf(jump_points[i+1] - tol)

    i <- i + 1
  }

  for (i in 1:n){

    a[i] <- min(c(jump_points[which(Vec_Mixed_cdf(jump_points) + STATISTIC >= i/n + tol)[1]], Inf), na.rm = TRUE)
    b[i] <- min(c(jump_points[which(Vec_Mixed_cdf(jump_points) - STATISTIC > (i-1)/n - tol)[1]], Inf), na.rm = TRUE)


    f_a[i] <- ifelse(!(a[i] %in% jump_points), Vec_Mixed_cdf(a[i]), Vec_Mixed_cdf(a[i] - tol))

    f_b[i] <- Vec_Mixed_cdf(b[i])


    temp_a[i] <- i/n - STATISTIC + tol
    temp_b[i] <- (i-1)/n + STATISTIC - tol

    for (j in 1:temp_length){


      if ((temp_a[i] > temp[2*j-1]) && (temp_a[i] <= temp[2*j])){

        f_a[i] <- i/n - STATISTIC
      }

      if ((temp_b[i] >= temp[2*j-1]) && (temp_b[i] < temp[2*j])){

        f_b[i] <- (i-1)/n + STATISTIC
      }



      j <- j + 1
    }

    if (f_a[i] < 0){
      f_a[i] <- 0
    }
    if (f_b[i] > 1){
      f_b[i] <- 1
    }

    i <- i + 1

  }

  df <- data.frame(rbind(f_b, f_a))
  write.table(df,"Boundary_Crossing_Time.txt", sep = ", ", row.names = FALSE, col.names = FALSE)

#  if(exact){
    PVAL <- ksgeneral::cont_ks_distribution_Rcpp(n)
#  }

  nm_alternative <- "two-sided"

  names(STATISTIC) <- "D"

#  if (is.null(PVAL)){
#    PVAL <- 0
#  }

  PVAL <- min(1.0, max(0.0, PVAL))
  RVAL <- list(statistic = STATISTIC, p.value = PVAL, alternative = nm_alternative, method = METHOD, data.name = DNAME)

  class(RVAL) <- "htest"
  return(RVAL)

}


