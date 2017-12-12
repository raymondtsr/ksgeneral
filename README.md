![downloads](https://cranlogs.r-pkg.org/badges/grand-total/KSgeneral)
![downloads](https://cranlogs.r-pkg.org/badges/KSgeneral)
![downloads](https://cranlogs.r-pkg.org/badges/last-week/KSgeneral)

# KSgeneral
Computing (complementary)cdf for one-sample Kolmogorov-Smirnov statistic when the underlying distribution is purely discrete, mixed, or continuous.

# Installation
In order to build the KSgeneral package from source, a C++ compiler is required. The latter is contained in the Windows Rtools, available from https://cran.r-project.org/bin/windows/Rtools/, or under MacOS in Xcode, downloadable from the App Store.
The package KSgeneral uses Rcpp in R, and utilizes the C++ code that efficiently computes the complementary cdf using the Exact-KS-FFT method developed in Dimitrova, Kaishev, Tan (2017), available from http://openaccess.city.ac.uk/18541.
Since the latter requires computation of Fast Fourier Transform (FFT), the FFTW3 library developed by Matteo Frigo and Steven G.Johnson needs to be installed from http://www.fftw.org/index.html. 
It should be noted that the Rtools and FFTW3 should be installed in the system PATH.
The recommended way to building the FFTW3 library on Windows is to install the free MinGW Unix environment so that one can use the GNU C compiler. (Rtools use MinGW compiler as well). If you are using R on MacOS, the FFTW3 library is typically installed in /usr/local.

After installng the FFTW3 library and C++ compilers, you can download the source file from Github and build the packages KSgeneral from source, by running install.packages("KSgeneral_0.1.0.tar.gz", repos = NULL, type = "source", INSTALL_opts = "--no-multiarch") in RStudio, when the source package is located in the current working directory; or running R CMD INSTALL --no-multiarch KSgeneral_0.1.0.tar.gz from the command line.      
