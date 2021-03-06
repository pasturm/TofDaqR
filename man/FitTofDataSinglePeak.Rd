% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/SupplementaryFunctions.R
\name{FitTofDataSinglePeak}
\alias{FitTofDataSinglePeak}
\title{Wrapper function for FitSinglePeak().}
\usage{
FitTofDataSinglePeak(PeakTofData, peakType, SpecNo = NULL)
}
\arguments{
\item{PeakTofData}{List from the output of \code{\link{GetTofDataSinglePeak}}.}

\item{peakType}{Peak model to use.}

\item{SpecNo}{Spectrum to use. If \code{SpecNo=NULL} (default) then the averaged
spectrum is used, with \code{SpecNo=i}, the i-th spectrum from TofData is
used.}
}
\value{
A list with fitted parameters.
}
\description{
\code{FitTofDataSinglePeak} is a wrapper for \code{\link{FitSinglePeak}}
and also reports fitted resolving power, area and count rates. It takes the
output of \code{\link{GetTofDataSinglePeak}} as input parameter and calculates
a peak fit of the averaged spectrum or of a single spectrum.
}
\details{
Peak types:
\tabular{cccc}{
peak type index \tab peak function \tab symmetric \tab baseline \cr
-1 \tab Gaussian through highest point and its neighbours \tab yes \tab no \cr
0 \tab Gaussian \tab yes \tab no  \cr
1 \tab Lorentzian \tab yes \tab no \cr
2 \tab Pseudo-Voigt \tab yes \tab no \cr
3 \tab Gaussian \tab yes \tab yes \cr
4 \tab Lorentzian \tab yes \tab yes \cr
5 \tab Pseudo-Voigt \tab yes \tab yes \cr
6 \tab Gaussian \tab no \tab no \cr
7 \tab Lorentzian \tab no \tab no \cr
8 \tab Pseudo-Voigt \tab no \tab no \cr
9 \tab Gaussian \tab no \tab yes \cr
10 \tab Lorentzian \tab no \tab yes \cr
11 \tab Pseudo-Voigt \tab no \tab yes
}

Note: this function is not part of the TofDaq API, but is included in the
package for convenience.
}
\seealso{
Other wrapper functions: 
\code{\link{GetTofDataSinglePeak}()}
}
\concept{wrapper functions}
