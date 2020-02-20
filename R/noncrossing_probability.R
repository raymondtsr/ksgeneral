#' Compute the joint probability of ordered random vaiables
#' 
#' Compute the probability P(l_i<X_i<h_i for all i=1,...,n) 
#' where X_1,..., X_n are ordered continuous random vaiables
#' 
#' @param h numeric(n), the higher bound of the ordered random vaiables
#' @param l numeric(n), the lower bound of ordered random vaiables
#' @param CDF function, a CDF function that specify the distribution of
#' the random variable. The density of the random variable at x must be
#' able to be computed via CDF(x).
compute_noncrossing_prob<-function(h,l,CDF = punif,FFT= TRUE){
    if(length(l)!=length(h))
        stop("The length of upper and lower bounds must be the same")
    if(any(l>h)) 
        return(0)
    n=length(l)
    ## make sure l and h are increasing sequences
    for(i in seq_len(n-1)){
        if(l[i]>l[i+1]) l[i+1] <- l[i]
        j <- n-i
        if(h[j] > h[j+1]) h[j] <- h[j+1]
    }
    h<- CDF(h)
    l<- CDF(l)
    ecdf_noncrossing_probability_wrapper(n,h,l,FFT)
}