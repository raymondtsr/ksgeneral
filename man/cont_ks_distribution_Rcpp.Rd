\name{cont_ks_distribution_Rcpp}
\alias{cont_ks_distribution_Rcpp}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Function calling the C++ routines that compute the complementary cumulative distribution function for the one-sample two-sided Kolmogorov-Smirnov statistic
%%  ~~function to do ... ~~
}
\description{
Function calling the C++ routines that compute the complementary cumulative distribution function for the one-sample two-sided Kolmogorov-Smirnov statistic, given the sample size $n$ and the file "Boundary_Crossing_Time.txt".
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
cont_ks_distribution_Rcpp(n)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{n}{
  the sample size
%%     ~~Describe \code{n} here~~
}
}
\details{
The FFT-based algorithm requires calculation in loops, which is more efficiently carried out in C++.
Hence, the FFT-based algorithm is conducted in C++ and the result is then called by R.
%%  ~~ If necessary, more details than the description above ~~
}
\value{
the complementary cdf for the KS statistic, given particular $n$ and file "Boundary_Crossing_Time.txt"
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
Moscovich A, Nadler B (2016). "Fast Calculation of Boundary Crossing Probabilities for Poisson Processes" arXiv preprint at arXiv:1503.04363v3.

Dimitrova D, Kaishev V, Tan S (2017) "Computing the Two-Sided Kolmogorov-Smirnov Distribution When the Underlying cdf is (Dis)continuous".
%% ~put references to the literature/web site here ~
}
\author{
Senren Tan
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{


## The function is currently defined as
function (n)
{
    .Call(`_ksgeneral_cont_ks_distribution_Rcpp`, n)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
% \keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
% \keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
