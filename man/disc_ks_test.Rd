\name{disc_ks_test}
\alias{disc_ks_test}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
KS Tests when the underlying cdf is purely discrete
%%  ~~function to do ... ~~
}
\description{
Perform a one-sample two-sided exact Kolmogorov-Smirnov test when the underlying cdf is purely discrete, similarly to \code{ks.test} from package \pkg{dgof}, but using the FFT-based algorithm to improve the speed and accuracy.
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
disc_ks_test(x, y, ..., exact = NULL, tol = 1e-08, sim.size = 1e+06, num.sim = 10)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
  a numeric vector of data values.
%%     ~~Describe \code{x} here~~
}
  \item{y}{
  \code{y} is the pre-specified discrete distribution under the null hypothesis.
  Note that \code{y} should be a step function within the class: \code{stepfun}, of which \code{ecdf} is a subclass.
%%     ~~Describe \code{y} here~~
}
  \item{\dots}{
      parameters of the distribution specified (as a character string) by \code{y}.
%%     ~~Describe \code{\dots} here~~
}
  \item{exact}{
  Whether the exact p-value of the KS test will be calculated.
  The default value is \code{NULL}.

  With \code{exact = TRUE}, the exact p-value is calculated.

  With \code{exact = FALSE}, the asymptotic p-value is applied to approximate the exact complementary, based on Monte Carlo simulation approach proposed by Wood and Altavela (1978).

  With \code{exact = NULL}, when the sample size \code{n} <= 100000, the exact p-value is obtained. Otherwise, the asymptotic p-value is obtained.
%%     ~~Describe \code{exact} here~~
}
  \item{tol}{
    epsilon in equation (3) of Dimitrova, Kaishev and Tan (2017) when calculating the exact complementary cdf of $D_\{n\}$, with default value \code{tol = 1e-08}.
    Value of NA or 0 can risk making errors due to numerical imprecisions.
%%     ~~Describe \code{tol} here~~
}
  \item{sim.size}{
   the number of simulated trajectories in the Monte Carlo simulation to approximate the asymptotic complementary cdf using Wood and Altavela (1978)'s approach.
  The default value is \code{sim.size = 1e+06}.
%%     ~~Describe \code{sim.size} here~~
}
  \item{num.sim}{
  the ample size in the Monte Carlo simulation to approximate the asymptotic complementary cdf using Wood and Altavela (1978)'s approach.
  The default value is \code{num.sim = 10}.
%%     ~~Describe \code{num.sim} here~~
}
}
\details{
Suppose we have an ordered set of $n$ random variates following a pre-specified discrete distribution with cdf $F(x)$, $x_{1} <= ... <= x_{n}$, Kolmogorov suggested $D_{n} = sup_\{x\}|F(x) - F_\{n\}(x)|$ as a goodness-of-fit measure, where: $F_\{n\}(x)$ denotes the empirical cdf of the random sample.
It is well known that when $F(x)$ is not continuous, the p-value of the one sample two-sided KS test depends on $F(x)$ and it is not easy to be calculated.
Arnold and Emerson (2011) provide the function \code{ks.test} as part of the package \pkg{dgof} that calculates exact p-values of the K-S test assuming $F(x)$ is purely discrete.
In \code{ks.test}, a one-sided K-S p-value is calculated by combining the approaches of Conover (1972) and Niederhausen (1981), while two-sided K-S p-values are calculated by combining the approaches of Gleser (1985) and Niederhausen (1981). However, the \code{ks.test} function due to Arnold and Emerson (2011) only provides exact p-values for \code{n} <= 30, since as noted by the authors, when \code{n} is large, numerical instabilities may occur. In the latter case, Arnold and Emerson (2011) suggest using simulation to approximate p-values, which may be rather slow and inaccurate.

Making use of FFT-based algorithm, the function \code{disc_ks_test} is more efficient in terms of both speed and accuracy.
When the sample size \code{n} is large i.e., n >= 100000, the asymptotic p-value is estimated through Monte Carlo simulation proposed by Wood and Altavela (1978).
The latter has been incorporated into \code{disc_ks_test}.
However, if the sample size \code{n} is small, the asymptotic p-value may not be accurate enough. So, the input \code{exact = FALSE} is not suggested when the sample size is small.
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
%% ...
}
\references{
Arnold TA, Emerson JW (2011). Nonparametric Goodness-of-Fit Tests for Discrete Null Distributions. The R Journal, 3(2), 34-39.

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
\code{\link{dgof::ks.test}}
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
# Comparison of results obtained from dgof::ks.test and ksgeneral::disc_ks_test:

x3 <- sample(1:10, 25, replace = TRUE)
ksgeneral::disc_ks_test(x3, ecdf(1:10), exact = TRUE)
dgof::ks.test(x3, ecdf(1:10), exact = TRUE)
ksgeneral::disc_ks_test(x3, ecdf(1:10), exact = TRUE)$p - dgof::ks.test(x3, ecdf(1:10), exact = TRUE)$p

x4 <- sample(1:10, 500, replace = TRUE)
ksgeneral::disc_ks_test(x4, ecdf(1:10), exact = TRUE)
dgof::ks.test(x4, ecdf(1:10), exact = TRUE)
ksgeneral::disc_ks_test(x4, ecdf(1:10), exact = TRUE)$p - dgof::ks.test(x4, ecdf(1:10), exact = TRUE)$p

# Using stepfun() to specify the same discrete distribution as defined by ecdf():

steps <- stepfun(1:10, cumsum(c(0, rep(0.1, 10))))
ksgeneral::disc_ks_test(x3, steps, exact = TRUE)

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%%\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
%%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
