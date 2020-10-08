![downloads](https://cranlogs.r-pkg.org/badges/grand-total/KSgeneral)
![downloads](https://cranlogs.r-pkg.org/badges/KSgeneral)
![downloads](https://cranlogs.r-pkg.org/badges/last-week/KSgeneral)
[![Rdoc](http://www.rdocumentation.org/badges/version/KSgeneral)](http://www.rdocumentation.org/packages/KSgeneral)
[![CRAN\_Status\_Badge](https://www.r-pkg.org/badges/version/KSgeneral)](https://cran.r-project.org/package=KSgeneral)

# KSgeneral

Computes a p-value of the one-sample two-sided (or one-sided, as a special case) Kolmogorov-Smirnov (KS) statistic, for any fixed critical level, and an arbitrary, possibly large sample size for a pre-specified purely discrete, mixed or continuous cumulative distribution function (cdf) under the null hypothesis.

If a data sample is supplied, 'KSgeneral' (also available from  https://CRAN.R-project.org/package=KSgeneral) computes the p-value corresponding to the value of the KS test statistic computed based on the user provided data sample.

The package 'KSgeneral' implements a novel, accurate and efficient method named Exact-KS-FFT, developed by Dimitrova, Kaishev, Tan (2017), available together with the underlying C++ code from http://openaccess.city.ac.uk/18541.

**To cite this package in publication use: Dimitrina S. Dimitrova, Vladimir K. Kaishev, and Senren Tan. Computing the Kolmogorov-Smirnov Distribution When the Underlying CDF is Purely Discrete, Mixed, or Continuous. *Journal of Statistical Software*. 2020, 95(10): 1–42. <doi:10.18637/jss.v095.i10>.**

The p-value is expressed as a double-boundary non-crossing probability for a homogeneous Poisson process, which is then efficiently computed using Fast Fourier Transform (FFT). 

The package can also be used to compute and plot the complementary cdf of the KS statistic which is known to depend on the hypothesized distribution when the latter is discontinuous (i.e. purely discrete or mixed).


# Installation
In order to build the KSgeneral package from source, a C++ compiler is required. 

The latter is contained in the Windows Rtools, available from https://cran.r-project.org/bin/windows/Rtools/, or under MacOS in Xcode, downloadable from the App Store.

The package KSgeneral uses Rcpp in R, and utilizes the C++ code that efficiently computes the complementary cdf using the Exact-KS-FFT method developed by Dimitrova, Kaishev, Tan (2017), available together with the underlying C++ code from http://openaccess.city.ac.uk/18541.

Since the latter requires computation of Fast Fourier Transform (FFT), the FFTW3 library developed by Matteo Frigo and Steven G.Johnson needs to be installed from http://www.fftw.org/index.html. 

It should be noted that the Rtools and FFTW3 should be installed in the system PATH.

For Windows users, The FFTW3 library (static library, with a ".a" extension) for Windows (32-bit or 64-bit) can be found in the local323.zip file, available from http://www.stats.ox.ac.uk/pub/Rtools/libs.html.

For Mac or Unix users, it is straightforward to install the FFTW3 library from the command line, following the instructions from http://www.fftw.org/index.html.
