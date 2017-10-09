\name{cont_ks_c_cdf}
\alias{cont_ks_c_cdf}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Kolmogorov $D_\{n\}$ complementary distribution when the underlying cdf is continuous
%%  ~~function to do ... ~~
}
\description{
Complementary cumulate distribution function (cdf) for one-sample two-sided Kolmogorov-Smirnov statistic, $P(D_\{n\} >= q), when the underlying cdf is continuous.
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
cont_ks_c_cdf(q, n)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{q}{
    the argument for the cdf of $D_\{n\}$
%%     ~~Describe \code{q} here~~
}
  \item{n}{
    the sample size
%%     ~~Describe \code{n} here~~
}
}
\details{
Given an ordered set of $n$ standard uniform variates, $u_{1} < ... < u_{n}$, Kolmogorov suggested $D_{n} = max\{D_\{n\}^-, D_\{n\}^+\} as a goodness-of-fit measure, where: $D_\{n\}^- = max_\{i = 1, ..., n\}[u_\{i\} - (i-1)/n]$ and $D_\{n\}^+ = max_\{i = 1, ..., n\}[i/n - u_\{i\}]$.

cont_ks_cdf provides the FFT-based algorithm proposed by Moscovich and Nadler (2015) to compute the cdf of $D_\{n\}$, $P(D_\{n\} < q)$, when the underlying cdf in the null hypothesis is continuous.
This algorithm is more efficient and numerically stable than the algorithm proposed by Wang, Tsang, and Marsaglia (2003), which was used by many existing packages computing the cdf of $D_\{n\}$, e.g., function \code{ks.test} in package \pkg{stat}, function \code{ks.test} in package \pkg{dgof}.
More precisely, in these packages, exact p-value is computed if the sample size $n$ is less than 100, and the computation time is $O(n^\{3\})$.
The FFT-based algorithm, on the other hand, can provide results with at least 10 digits of accuracy after the decimal points, for sample sizes $n$ up to 100000.
For $n > 100000$, results can still be computed, but with a cost of time.
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
\source{
The FFT-based algorithm to compute the distribution of $D_{n}$ when the underlying cdf is continuous, comes via Moscovich and Nadler (2016).
Also see Dimitrova, Kaishev, Tan (2017) for more details.
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

ksgeneral::cont_ks_c_cdf(0.05, 100)

n<-100
x<-1:100/500
plot(x, sapply(x, function(x) ksgeneral::cont_ks_c_cdf(x, n)), type='l')

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
% \keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
% \keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
