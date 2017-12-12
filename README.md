![downloads](https://cranlogs.r-pkg.org/badges/grand-total/KSgeneral)
![downloads](https://cranlogs.r-pkg.org/badges/KSgeneral)
![downloads](https://cranlogs.r-pkg.org/badges/last-week/KSgeneral)

# ksgeneral
Computing (complementary)cdf for one-sample Kolmogorov-Smirnov statistic when the underlying distribution is continuous, mixed, or purely discrete.

# Installation
To install the R package, one has to first install the library FFTW3 available from http://www.fftw.org/index.html. And make sure that the library is installed in the systems path. If you are building the 32-bit package on Windows, you can install the 32-bit FFTW3 library, and if you are building the 64-bit package on Windows, you can install the 64-bit FFTW3 library. The recommended way to building the FFTW3 library on Windows is to install the free MinGW Unix environment so that one can use the GNU C compiler. (Rtools use MinGW compiler as well). If you are using R on MacOS, the FFTW3 library is typically installed in /usr/local.

Also, since it is required to build the R package from source, if you are working on Windows, you can download Rtools available from https://cran.r-project.org/bin/windows/Rtools/. And make sure that the Rtools is installed in the systems path. If you are working on MacOS, you can download Xcode from the App Store.

After installng the FFTW3 library and C++ compilers, you can download the source file from Github and build the packages ksgeneral from source, by running install.packages("ksgeneral_0.1.0.tar.gz", repos = NULL, type = "source", INSTALL_opts = "--no-multiarch") in RStudio, when the source file is located in the current working directory; or running R CMD INSTALL --no-multiarch ksgeneral_0.1.0.tar.gz from the command line.      
