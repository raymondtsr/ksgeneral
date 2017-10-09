\name{mixed_ks_test}
\alias{mixed_ks_test}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
KS Tests when the underlying cdf is mixed
%%  ~~function to do ... ~~
}
\description{
Perform a one-sample two-sided exact Kolmogorov-Smirnov test when the underlying cdf is mixed, using the FFT-based algorithm proposed in Dimitrova, Kaishev and Tan (2017).
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
mixed_ks_test(x, jump_points, Mixed_dist, ..., tol = 1e-10)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
    a numeric vector of data values.
%%     ~~Describe \code{x} here~~
}
  \item{jump_points}{
    a vector containing the points where the underlying cdf has jump(s)
%%     ~~Describe \code{jump_points} here~~
}
  \item{Mixed_dist}{
      \code{Mixed_dist} is the pre-specified mixed distribution under the null hypothesis, defined by the users.
%%     ~~Describe \code{Mixed_dist} here~~
}
  \item{\dots}{
        parameters of the distribution specified (as a character string) by \code{Mixed_dist}.
%%     ~~Describe \code{\dots} here~~
}
  \item{tol}{
   epsilon in equation (3) of Dimitrova, Kaishev and Tan (2017) when calculating the exact complementary cdf of $D_\{n\}$, with default value \code{tol = 1e-08}.
      Value of NA or 0 can risk making errors due to numerical imprecisions.
%%     ~~Describe \code{tol} here~~
}
}
\details{
To the best of our knowledge, no methods have been proposed in the literature to compute the p-value of one-sample two-sided KS test statistic, when F(x) is mixed.
The approach in this function, described in Section 2.1 of Dimitrova, Kaishev and Tan (2017), is to express $P(D_\{n\} >= q)$ as an appropriate rectangle probability for uniform order statistics, as noted by Gleser (1985), and to compute the latter probability using the Fast Fourier Transform (FFT) method.
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
A list with class "htest" containing the following components:
%%  If it is a LIST, use
  \item{statistic }{the value of the statistic.}
  \item{p.value }{the p-value of the test.}
  \item{alternative }{"two-sided".}
  \item{method }{a character string indicating what type of test was performed.}
  \item{data.name }{a character string giving the name of the data.}
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
Dimitrova D, Kaishev V, Tan S (2017) "Computing the Two-Sided Kolmogorov-Smirnov Distribution When the Underlying cdf is (Dis)continuous".

Gleser LJ (1985). "Exact Power of Goodness-of-Fit Tests of Kolmogorov Type for Discontinuous Distributions". Journal of the American Statistical Association, 80(392), 954-958.
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
# Example: Computing the p-value of the one-sample two-sided KS test,
# when the underlying distribution is a mixed distribution
# with two jumps at 0 and log(2.5):
#

Mixed_cdf_example <- function(x)
{
     result <- 0
     if (x < 0){
         result <- 0
     }
     else if (x == 0){
         result <- 0.5
     }
     else if (x < log(2.5)){
         result <- 1 - 0.5 * exp(-x)
     }
     else{
         result <- 1
     }

     return (result)
}
test_data <- c(0,0,0,0,0,0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,log(2.5),log(2.5),log(2.5),log(2.5),log(2.5),log(2.5))
ksgeneral::mixed_ks_test(test_data, c(0, log(2.5)), Mixed_cdf_example)


}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%%\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
%%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
