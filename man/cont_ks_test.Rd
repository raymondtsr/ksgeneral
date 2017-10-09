\name{cont_ks_test}
\alias{cont_ks_test}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Kolmogorov-Smirnov Tests when the underlying cdf is continuous
%%  ~~function to do ... ~~
}
\description{
Perform a one-sample two-sided exact Kolmogorov-Smirnov test, similarly to \code{ks.test} from package \pkg{stats}, but using the FFT-based algorithm to improve the speed and accuracy.
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
cont_ks_test(x, y, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
  a numeric vector of data values.
%%     ~~Describe \code{x} here~~
}
  \item{y}{
  either a numeric vector of data values, or a character string naming a continuous cumulative distribution function such a \code{pexp}, \code{pnorm}.
  Only continuous cdfs are valid!
%%     ~~Describe \code{y} here~~
}
  \item{\dots}{
  parameters of the distribution specified (as a character string) by \code{y}.
%%     ~~Describe \code{\dots} here~~
}
}
\details{
This routine is equivalent to \code{stats::ks.test(x, y, ..., exact = TRUE)}, but uses the FFT-based method based on \code{ksgeneral::cont_ks_distribution}.
For more details about the arguments, please refer to the documentation for \code{stats::ks.test}.
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
\source{
The FFT-based algorithm to compute the distribution of $D_{n}$ when the underlying cdf is continuous, comes via Moscovich and Nadler (2016).
Also see Dimitrova, Kaishev, Tan (2017) for more details.
}

\references{
Dimitrova D, Kaishev V, Tan S (2017) "Computing the Two-Sided Kolmogorov-Smirnov Distribution When the Underlying cdf is (Dis)continuous".

Moscovich A, Nadler B (2016). "Fast Calculation of Boundary Crossing Probabilities for Poisson Processes" arXiv preprint at arXiv:1503.04363v3.
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
x<-abs(rnorm(100))
p.kt <- ks.test(x, "pexp", exact = TRUE)$p
p.kt_fft <- ksgeneral::cont_ks_test(x, "pexp")$p
abs(p.kt-p.kt_fft)


}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
% \keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
% \keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
