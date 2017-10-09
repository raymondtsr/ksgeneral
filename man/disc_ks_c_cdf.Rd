\name{disc_ks_c_cdf}
\alias{disc_ks_c_cdf}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Complementary cdf of KS test statistic $D_\{n\}$ when the underlying cdf is purely discrete
%%  ~~function to do ... ~~
}
\description{
Complementary cdf for one-sample two-sided Kolmogorov-Smirnov statistic $D_\{n\}$, $P(D_\{n\} >= q)$, when the underlying cdf is purely discrete.
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
disc_ks_c_cdf(q, n, y, ..., exact = NULL, tol = 1e-08, sim.size = 1e+06, num.sim = 10)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{q}{
    the argument for the complementary cdf of $D_\{n\}$
%%     ~~Describe \code{STATISTIC} here~~
}
  \item{n}{
    the sample size
%%     ~~Describe \code{n} here~~
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
  Whether the exact result of the complementary cdf will be calculated.
  The default value is \code{NULL}.

  With \code{exact = TRUE}, the exact complementary cdf is calculated.

  With \code{exact = FALSE}, the asymptotic complementary cdf is applied to approximate the exact complementary, based on Monte Carlo simulation approach proposed by Wood and Altavela (1978).

  With \code{exact = NULL}, when the sample size \code{n} <= 100000, the exact complementary cdf is obtained. Otherwise, the asymptotic complementary cdf is obtained.
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
It is well known that when $F(x)$ is not continuous, the (complementary) cdf of $D_\{n\}$ depends on $F(x)$ and it is not easy to be calculated.
Arnold and Emerson (2011) provide the function \code{ks.test} as part of the package \pkg{dgof} that calculates exact p-values of the K-S test assuming $F(x)$ is purely discrete.
In \code{ks.test}, a one-sided K-S p-value is calculated by combining the approaches of Conover (1972) and Niederhausen (1981), while two-sided K-S p-values are calculated by combining the approaches of Gleser (1985) and Niederhausen (1981). However, the \code{ks.test} function due to Arnold and Emerson (2011) only provides exact p-values for \code{n} <= 30, since as noted by the authors, when \code{n} is large, numerical instabilities may occur. In the latter case, Arnold and Emerson (2011) suggest using simulation to approximate p-values, which may be rather slow and inaccurate.

Making use of FFT-based algorithm, our function \code{disc_ks_c_cdf} is more efficient in terms of both speed and accuracy.
When the sample size \code{n} is large i.e., n >= 100000, the asymptotic complementary cdf of $D_\{n\}$ is estimated through Monte Carlo simulation proposed by Wood and Altavela (1978).
The latter has been incorporated into \code{disc_ks_c_cdf}.
%%  ~~ If necessary, more details than the description above ~~
}
\value{
Returns $P(D_\{n\} >= q)$.
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
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
## Example to compute the complementary cdf for $D_\{n\}$ when the underlying cdf $F(x)$ is a binomial distribution

binom_3 <- stepfun(c(0:3), c(0,pbinom(0:3,3,0.5)))
ksgeneral::disc_ks_c_cdf(0.05, 400, binom_3)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
% \keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
% \keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
