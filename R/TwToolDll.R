# FitSinglePeak --------------------------------------------------------------
#' Performs a peak fit.
#'
#' \code{FitSinglePeak} performs a peak fit. Gaussian, Lorentzian and
#' Pseudo-Voigt peak shapes are supported.
#'
#' Peak parameters are optimized using a Levenberg-Marquardt algorithm as
#' implemented by the library \code{lmfit}. The following peak types are available:
#' \tabular{cccc}{
#' peak type index \tab peak function \tab symmetric \tab baseline \cr
#' -1 \tab Gaussian through highest point and its neighbours \tab yes \tab no \cr
#' 0 \tab Gaussian \tab yes \tab no  \cr
#' 1 \tab Lorentzian \tab yes \tab no \cr
#' 2 \tab Pseudo-Voigt \tab yes \tab no \cr
#' 3 \tab Gaussian \tab yes \tab yes \cr
#' 4 \tab Lorentzian \tab yes \tab yes \cr
#' 5 \tab Pseudo-Voigt \tab yes \tab yes \cr
#' 6 \tab Gaussian \tab no \tab no \cr
#' 7 \tab Lorentzian \tab no \tab no \cr
#' 8 \tab Pseudo-Voigt \tab no \tab no \cr
#' 9 \tab Gaussian \tab no \tab yes \cr
#' 10 \tab Lorentzian \tab no \tab yes \cr
#' 11 \tab Pseudo-Voigt \tab no \tab yes
#' } The function for peakType = -1 is not an actual peak fit but rather the
#' Gaussian function (symmetric, no baseline) defined by the maximum intensity
#' point and its neighbors. It is used to automatically generate guess values
#' for position, width and amplitude if all peak parameters are set to 0.
#'
#' @param yVals y axis data
#' @param xVals x axis data
#' @param peakType peak model to use.
#' @param blOffset Initial value of baseline offset at first data point.
#' @param blSlope Initial value of slope of baseline.
#' @param amplitude Initial value of peak amplitude.
#' @param fwhmLo Initial value of left peak width (full width at half maximum).
#' @param fwhmHi Initial value of right peak width (full width at half maximum).
#' @param peakPos Initial value of peak position.
#' @param mu Initial value of Gaussian/Lorentzian contribution.
#' @return List with peak parameters.
#'
#' @family Peak fitting functions
FitSinglePeak <- function(yVals, xVals, peakType = 0, blOffset = 0, blSlope = 0,
                          amplitude = 0, fwhmLo = 0, fwhmHi = 0, peakPos = 0, mu = 0) {
  .Call("FitSinglePeak", yVals, xVals, peakType, blOffset, blSlope,
             amplitude, fwhmLo, fwhmHi, peakPos, mu)
}

# FitSinglePeak2 --------------------------------------------------------------
#' Performs a peak fit (Initial values as vector).
#'
#' \code{FitSinglePeak2} performs a peak fit. Gaussian, Lorentzian and
#' Pseudo-Voigt peak shapes are supported.
#'
#' Same as \code{\link{FitSinglePeak}}, but takes the initial values of the
#' fit parameter as a vector argument instead of explicit parameters.
#' Peak parameters are optimized using a Levenberg-Marquardt algorithm as
#' implemented by the library \code{lmfit}. The following peak types are available:
#' \tabular{cccc}{
#' peak type index \tab peak function \tab symmetric \tab baseline \cr
#' -1 \tab Gaussian through highest point and its neighbours \tab yes \tab no \cr
#' 0 \tab Gaussian \tab yes \tab no  \cr
#' 1 \tab Lorentzian \tab yes \tab no \cr
#' 2 \tab Pseudo-Voigt \tab yes \tab no \cr
#' 3 \tab Gaussian \tab yes \tab yes \cr
#' 4 \tab Lorentzian \tab yes \tab yes \cr
#' 5 \tab Pseudo-Voigt \tab yes \tab yes \cr
#' 6 \tab Gaussian \tab no \tab no \cr
#' 7 \tab Lorentzian \tab no \tab no \cr
#' 8 \tab Pseudo-Voigt \tab no \tab no \cr
#' 9 \tab Gaussian \tab no \tab yes \cr
#' 10 \tab Lorentzian \tab no \tab yes \cr
#' 11 \tab Pseudo-Voigt \tab no \tab yes
#' } The function for peakType = -1 is not an actual peak fit but rather the
#' Gaussian function (symmetric, no baseline) defined by the maximum intensity
#' point and its neighbors. It is used to automatically generate guess values
#' for position, width and amplitude if all peak parameters are set to 0.
#'
#' @param yVals y axis data
#' @param xVals x axis data
#' @param peakType peak model to use.
#' @param param Vector of initial values (blOffset, blSlope, amplitude, fwhmLo,
#' fwhmHi, peakPos, mu).
#' @return List with peak parameters.
#'
#' @family Peak fitting functions
FitSinglePeak2 <- function(yVals, xVals, peakType = 0, param = rep(0, 7)) {
  .Call("FitSinglePeak2", yVals, xVals, peakType, param)
}

# EvalSinglePeak -------------------------------------------------------------
#' Evaluates a peak fit.
#'
#' \code{EvalSinglePeak} calculates the y-axis values for a given set of peak
#' parameters.
#'
#' @param xVals x axis data
#' @param blOffset Baseline offset at first data point.
#' @param blSlope Slope of baseline.
#' @param amplitude Peak amplitude.
#' @param fwhmLo Left peak width (full width at half maximum).
#' @param fwhmHi Right peak width (full width at half maximum).
#' @param peakPos Peak position.
#' @param mu Gaussian/Lorentzian contribution.
#' @return Vector with y axis data.
#'
#' @family Peak fitting functions
EvalSinglePeak <- function(xVals, blOffset = 0, blSlope = 0, amplitude, fwhmLo, fwhmHi = fwhmLo, peakPos, mu = 0) {
  .Call("EvalSinglePeak", xVals, blOffset, blSlope, amplitude, fwhmLo, fwhmHi, peakPos, mu)
}

# GetMoleculeMass ------------------------------------------------------------
#' Calculates the mass/charge ratio of a molecule.
#'
#' \code{GetMoleculeMass} parses the molecular formula and returns the
#' mass/charge ratio of the molecule/atom/ion.
#'
#' @param molecule Molecule string.
#' @return Mass/charge ratio.
#'
#' @family Chemistry functions
#'
#' @examples
#' GetMoleculeMass("CO2")
GetMoleculeMass <- function(molecule) {
  .Call("GetMoleculeMass", molecule)
}

# GetIsotopePattern ----------------------------------------------------------
#' Calculates the isotope pattern of a molecule.
#'
#' \code{GetIsotopePattern} parses the molecular formula and returns the
#' isotope pattern (mass and abundance).
#'
#' @param molecule Molecule string.
#' @param abundanceLimit Absolute abundance limit for the generated pattern.
#' @return List with mass and abundace vectors.
#'
#' @family Chemistry functions
#'
#' @examples
#' GetIsotopePattern("CO2", 1e-5)
GetIsotopePattern <- function(molecule, abundanceLimit) {
  rv <- .Call("GetIsotopePattern", molecule, abundanceLimit)
  if (rv$TwRetVal!=4) {warning(TwRetVal[rv$TwRetVal+1])} else {rv$TwRetVal <- NULL}
  return(rv)
}

# GetIsotopePattern2 ---------------------------------------------------------
#' Calculates the isotope pattern of a molecule.
#'
#' \code{GetIsotopePattern2} parses the molecular formula and returns the
#' isotope pattern (mass and abundance).
#'
#' Same as \code{GetIsotopePattern} but using an exact algorithm and is
#' therefore suitable only for rather small molecules.
#'
#' @param molecule Molecule string.
#' @param abundanceLimit Absolute abundance limit for the generated pattern.
#' @return List with mass and abundace vectors.
#'
#' @family Chemistry functions
#'
#' @examples
#' GetIsotopePattern2("CO2", 1e-5)
GetIsotopePattern2 <- function(molecule, abundanceLimit) {
  rv <- .Call("GetIsotopePattern2", molecule, abundanceLimit)
  if (rv$TwRetVal!=4) {warning(TwRetVal[rv$TwRetVal+1])} else {rv$TwRetVal <- NULL}
  return(rv)
}

# Tof2Mass -------------------------------------------------------------------
#' Converts from sample index to mass/charge.
#'
#' \code{Tof2Mass} converts from sample index to mass/charge.
#'
#' \tabular{cl}{
#' massCalibMode \tab Mass calibration function \cr
#' 0 \tab \eqn{i = p_1 \sqrt(m) + p_2} \cr
#' 1 \tab \eqn{i = p_1/\sqrt(m) + p_2} \cr
#' 2 \tab \eqn{i = p_1 m^{p_3} + p_2} \cr
#' 3 \tab \eqn{i = p_1 \sqrt(m) + p_2 + p_3 (m - p_4)^2} \cr
#' 4 \tab \eqn{i = p_1 \sqrt(m) + p_2 + p_3 m^2 + p_4 m + p_5} \cr
#' 5 \tab \eqn{m = p_1 i^2 + p_2 i + p_3}
#' }
#' Note: Modes 3 and 4 are flawed. Don't use them. In mode 3 the fit does not
#' converge well, because of a bug (parameters not correctly initialized).
#' Mode 4 is two sequential fits, first mode 0, then a quadratic fit to the
#' residuals, which is an inferior implementation of mode 3. Mode 1 is for FTMS
#' data.
#'
#' @param tofSample Vector of sample indices to convert.
#' @param massCalibMode Mass calibration function to use. See below.
#' @param p Vector containing the calibration parameters (number depends on
#' \code{MassCalibMode}, see below).
#' @return Mass/charge values.
#'
#' @seealso \code{\link{Mass2Tof}}
#'
#' @examples
#' Tof2Mass(100000, massCalibMode = 0, p = c(3,5))
Tof2Mass <- function(tofSample, massCalibMode = 0, p) {
  .Call("Tof2Mass", as.numeric(tofSample), massCalibMode, p)
}

# Mass2Tof -------------------------------------------------------------------
#' Converts from mass/charge to sample index.
#'
#' \code{Mass2Tof} converts from mass/charge to sample index.
#'
#' \tabular{cl}{
#' massCalibMode \tab Mass calibration function \cr
#' 0 \tab \eqn{i = p_1 \sqrt(m) + p_2} \cr
#' 1 \tab \eqn{i = p_1/\sqrt(m) + p_2} \cr
#' 2 \tab \eqn{i = p_1 m^{p_3} + p_2} \cr
#' 3 \tab \eqn{i = p_1 \sqrt(m) + p_2 + p_3 (m - p_4)^2} \cr
#' 4 \tab \eqn{i = p_1 \sqrt(m) + p_2 + p_3 m^2 + p_4 m + p_5} \cr
#' 5 \tab \eqn{m = p_1 i^2 + p_2 i + p_3}
#' }
#' Note: Modes 3 and 4 are flawed. Don't use them. In mode 3 the fit does not
#' converge well, because of a bug (parameters not correctly initialized).
#' Mode 4 is two sequential fits, first mode 0, then a quadratic fit to the
#' residuals, which is an inferior implementation of mode 3. Mode 1 is for FTMS
#' data.
#'
#' @param mass Vector of mass/charge values to convert.
#' @param massCalibMode Mass calibration function to use. See below.
#' @param p Vector containing the calibration parameters (number depends on
#' \code{MassCalibMode}, see below).
#' @return Sample indices.
#'
#' @seealso \code{\link{Tof2Mass}}
#'
#' @examples
#' Mass2Tof(100, massCalibMode = 0, p = c(3,5))
Mass2Tof <- function(mass, massCalibMode = 0, p) {
  .Call("Mass2Tof", as.numeric(mass), massCalibMode, p)
}

# MassCalibrate -------------------------------------------------------------------
#'  Performs a mass calibration.
#'
#' \code{MassCalibrate} performs a mass calibration for a list of
#' mass/sample index/weight values and for a given calibration mode.
#'
#' \tabular{cl}{
#' massCalibMode \tab Mass calibration function \cr
#' 0 \tab \eqn{i = p_1 \sqrt(m) + p_2} \cr
#' 1 \tab \eqn{i = p_1/\sqrt(m) + p_2} \cr
#' 2 \tab \eqn{i = p_1 m^{p_3} + p_2} \cr
#' 3 \tab \eqn{i = p_1 \sqrt(m) + p_2 + p_3 (m - p_4)^2} \cr
#' 4 \tab \eqn{i = p_1 \sqrt(m) + p_2 + p_3 m^2 + p_4 m + p_5} \cr
#' 5 \tab \eqn{m = p_1 i^2 + p_2 i + p_3}
#' }
#' Note: Modes 3 and 4 are flawed. Don't use them. In mode 3 the fit does not
#' converge well, because of a bug (parameters not correctly initialized).
#' Mode 4 is two sequential fits, first mode 0, then a quadratic fit to the
#' residuals, which is an inferior implementation of mode 3. Mode 1 is for FTMS
#' data.
#'
#' @param massCalibMode Mass calibration function to use. See below.
#' @param mass Vector of mass/charge values.
#' @param idx Vector of sample indices (of same length as \code{mass}).
#' @param weight Vector of weights (if \code{NULL} (default) all weights are set equal).
#' @return Vector of calibration parameters.
MassCalibrate <- function(massCalibMode = 0, mass, idx, weight = NULL) {
  if (is.null(weight)) {
    weight <- rep(1, length(mass))
  }
  .Call("MassCalibrate", massCalibMode, mass, idx, weight)
}

# SiInitializeHistograms -------------------------------------------------------
#' Initializes the single ion histogramming.
#'
#' \code{SiInitializeHistograms} initializes the single ion histogramming.
#'
#' This function must be called before calling the functions
#' \code{\link{SiSetProcessingOptions}}, \code{\link{SiProcessSpectrum}}, \code{\link{SiGetHistogram}},
#' \code{\link{SiResetHistograms}} or \code{\link{SiCleanup}}. Each mass range specified by
#' \code{loMass} and \code{hiMass} elements is associated with a spectrum type that allows to
#' get separate statistics for multi-spectrum acquisitions (bipolar or pulsed
#' experiments).
#'
#' @param loMass Vector of the lower borders of the mass ranges.
#' @param hiMass Vector of the upper borders of the mass ranges.
#' @param specType Vector of spectrum type indices (non-negative integers). If
#' specType is \code{NULL}, all mass ranges get a default spectrum type of 0.
#'
#' @family Single ion histogramming functions.
SiInitializeHistograms <- function(loMass, hiMass, specType = NULL) {
  if (is.null(specType)) {specType <- rep(0, length(loMass))}
  # Note: it is easier to do this here than to pass NULL to the API.
  rv <- .Call("SiInitializeHistograms", loMass, hiMass, as.integer(specType))
  if (rv!=4) {return(warning(TwRetVal[rv+1]))} else {return(TwRetVal[5])}
}

# SiSetProcessingOptions -------------------------------------------------------
#' Sets processing options for each spectrum type.
#'
#' \code{SiSetProcessingOptions} sets processing options for each spectrum type.
#'
#' Options:
#' \tabular{llc}{
#' Name \tab Description \tab Default value \cr
#' \code{MassCalibMode} \tab Mass calibration mode in use (see \code{\link{MassCalibrate}}). \tab 0 \cr
#' \code{MassCalibParam\emph{n}} \tab Mass calibration parameters \emph{n} =
#' 1..number of calibration parameters for the given mass calibration mode. \tab \code{c(1000, 0)} \cr
#' \code{FullScale} \tab Full digitizing range of the ADC in the same units as
#' the spectra to be analyzed (typically mV). \tab 500 \cr
#' \code{NbrBits} \tab ADC resolution (8 for AP240, 14 for ADQ114 etc.). \tab 8 \cr
#' \code{SampleInterval} \tab Sample interval in ns. \tab 1 \cr
#' \code{PreampGain} \tab Gain of external preamp. \tab 1 \cr
#' \code{PreSamples} \tab Number of samples before a threshold crosser taken into account. \tab 0 \cr
#' \code{PostSamples} \tab Number of samples after a negative threshold crosser taken into account. \tab 0 \cr
#' \code{BaselineAndThresholdFromData} \tab If >0 the baseline and threshold
#' values will be determined based on \code{NbrStdDevs} for every processed
#' spectrum. If >1.5 baseline noise is determined from a fit to a histogram of
#' all data (instead of from the standard deviation of all data). This makes the
#' noise determination more robust when real peaks are present in the spectrum. \tab 0 \cr
#' \code{NbrStdDevs} \tab Number of standard deviations of baseline noise that
#' defines the threshold. Only relevant if \code{BaselineAndThresholdFromData>0}. \tab 6 \cr
#' \code{Baseline} \tab Baseline value used for calculation of intensities. Has
#' no meaning if \code{BaselineAndThresholdFromData>0}. \tab 5 \cr
#' \code{Threshold} \tab Threshold value used for calculation of intensities.
#' Has no meaning if \code{BaselineAndThresholdFromData!=0}. \tab 8 \cr
#' \code{NegativeSignal} \tab Indicates peak polarity with respect to the baseline. \tab 0 (=FALSE) \cr
#' \code{BaselineAndThresholdInCodes} \tab Indicates whether the values in \code{Baseline}
#' and \code{Threshold} are interpreted as ADC codes or mV. \tab 1 (=TRUE)
#' }
#'
#' @param option Option to set (see below for a list of valid option strings).
#' @param value Value to set the for the given option.
#' @param specType Spectrum type index. -1 is a wildcard for all spectrum types.
#'
#' @family Single ion histogramming functions.
SiSetProcessingOptions <- function(option, value, specType = -1) {
  rv <- .Call("SiSetProcessingOptions", option, value, specType)
  if (rv!=4) {return(warning(TwRetVal[rv+1]))} else {return(TwRetVal[5])}
}

# SiProcessSpectrum ------------------------------------------------------------
#' Processes a spectrum.
#'
#' \code{SiProcessSpectrum} processes a spectrum according to the options set for
#' it's spectrum type.
#'
#' See \code{\link{SiProcessSpectrumFromShMem}} for another variant of this
#' function, where the spectrum is directly read from shared memory.
#'
#' @param spectrum Vector holding the spectrum to process.
#' @param specType Spectrum type index (non-negative integer).
#' @return A list with the baseline and threshold value.
#'
#' @family Single ion histogramming functions.
SiProcessSpectrum <- function(spectrum, specType = 0) {
  .Call("SiProcessSpectrum", spectrum, specType)
}

# SiGetHistogram ---------------------------------------------------------------
#' Gets a histogram of the single ion intensities.
#'
#' \code{SiGetHistogram} gets a histogram of the single ion intensities for a
#' mass range defined with \code{\link{SiInitializeHistograms}}.
#'
#' Note: R crashes if \code{histogramIndex} is set to max(histogramIndex)+1
#' (API bug).
#'
#' @param histogramIndex Index (0-based) of the histogram. It corresponds to the
#' mass range defined with \code{\link{SiInitializeHistograms}}.
#' @return A list with the intensities (histogram x-values), counts (histogram
#' y-values), the number of spectra that were processed for this histogram and
#' the mean histogram value i.e. sum(intensity[i]*counts[i])/sum(counts[i]).
#'
#' @family Single ion histogramming functions.
SiGetHistogram <- function(histogramIndex) {
  .Call("SiGetHistogram", histogramIndex)
}

# SiGetSumHistogram ------------------------------------------------------------
#' Gets a sum histogram of the single ion intensities.
#'
#' \code{SiGetSumHistogram} gets a histogram of the single ion intensities, which
#' is a sum over all histograms of a given \code{specType} within the rate and
#' mass range as specified by \code{minMass}, \code{maxMass}, \code{minRate} and \code{maxRate}.
#'
#' @param specType Spectrum type index (non-negative integer).
#' @param minMass Minimum mass for histogram filtering.
#' @param maxMass Maximum mass for histogram filtering.
#' @param minRate Minimum event count rate for histogram filtering.
#' @param maxRate Maximum event count rate for histogram filtering.
#' @return A list with the intensities (histogram x-values), counts (histogram
#' y-values), the number of spectra that were processed for this histogram and
#' the mean histogram value i.e. sum(intensity[i]*counts[i])/sum(counts[i]).
#'
#' @family Single ion histogramming functions.
SiGetSumHistogram <- function(specType, minMass, maxMass, minRate, maxRate) {
  .Call("SiGetSumHistogram", specType, minMass, maxMass, minRate, maxRate)
}

# SiResetHistograms ------------------------------------------------------------
#' Resets all histograms and spectrum counters to zero.
#'
#' \code{SiResetHistograms} resets all histograms and spectrum counters to zero.
#'
#' @family Single ion histogramming functions.
SiResetHistograms <- function() {
  rv <- .Call("SiResetHistograms")
  if (rv!=4) {return(warning(TwRetVal[rv+1]))} else {return(TwRetVal[5])}
}

# SiCleanup --------------------------------------------------------------------
#' Cleans up the state in the DLL.
#'
#' \code{SiCleanup} cleans up the state in the DLL. After a call to this
#' function \code{\link{SiInitializeHistograms}} can be called again.
#'
#' @family Single ion histogramming functions.
SiCleanup <- function() {
  rv <- .Call("SiCleanup")
  if (rv!=4) {return(warning(TwRetVal[rv+1]))} else {return(TwRetVal[5])}
}

# SiFitPhd ---------------------------------------------------------------------
#' Fits a (slightly modified) log-normal distribution to the histogram.
#'
#' \code{SiFitPhd} fits a (slightly modified) log-normal distribution to the
#' histogram data.
#'
#' The following equation with the parameters A, w, x0 and xc is fitted: \cr
#' \eqn{A/(\sqrt(2*pi)*w*(x-x0))*exp(-(log((x-x0)/xc)^2)/(2*w^2))}
#'
#' @param intensity Vector of intensities (histogram x-axis data).
#' @param counts Vector of counts (histogram y-axis data).
#' @return A list with the FWHM of the distribution, the position of the fitted
#' maximum and a vector with the values of the four fitting parameters.
#'
#' @family Single ion histogramming functions.
SiFitPhd <- function(intensity, counts) {
  .Call("SiFitPhd", intensity, counts)
}

# SiEvalPhd --------------------------------------------------------------------
#' Evaluates the fitted single ion distribution.
#'
#' \code{SiEvalPhd} evaluates the fitted single ion distribution.
#'
#' @param par Vector of fitted parameter values.
#' @param intensity Vector of intensities (histogram x-axis data).
#' @return Vector with y-axis data.
#'
#' @family Single ion histogramming functions.
SiEvalPhd <- function(par, intensity) {
  .Call("SiEvalPhd", par, intensity)
}

# SiFitRateFromPhd -------------------------------------------------------------
#' Fits an event rate to a multi-ion histogram.
#'
#' \code{SiFitRateFromPhd} takes a fitted single ion distribution as input
#' and fits an event rate to a multi-ion histogram (intensity and counts)
#' assuming poisson distribution of n-ion events.
#'
#' @param intensity Vector of intensities (histogram x-axis data).
#' @param counts Vector of counts (histogram y-axis data).
#' @param siPar Vector with fitted parameter values.
#' @return A list with the fitted rate (ions/extraction) and vector with the
#' fitted multi-ion distribution.
#'
#' @family Single ion histogramming functions.
SiFitRateFromPhd <- function(intensity, counts, siPar) {
  .Call("SiFitRateFromPhd", intensity, counts, siPar)
}

# FindTpsIp ------------------------------------------------------------------
#' Gets IP address of TPS2.
#'
#' \code{FindTpsIp} listens for UDP packets that TPS2 broadcast and returns
#' the IP of the TPS2.
#'
#' Note that executing this function makes your program to a UDP server and
#' Windows firewall (or other personal firewall software) will query for
#' permission.
#'
#' @param TpsSerial Serial number of TPS2.
#' @param timeout Timeout in ms to wait for the UDP packet.
#' @return String of IP address.
#'
#' @examples
#' FindTpsIp("910.33.0316", 500)
FindTpsIp <- function(TpsSerial, timeout) {
  rv <- .Call("FindTpsIp", TpsSerial, timeout)
  if (rv$TwRetVal!=4) {warning(TwRetVal[rv$TwRetVal+1])} else {rv$TwRetVal <- NULL}
  return(rv)
}

# Not implemented --------------------------------------------------------------
# TwTranslateReturnValue
# TwMultiPeakFit
# TwEvalMultiPeak
# TwFitResolution
# TwEvalResolution
# TwDecomposeMass
# TwGetComposition
# TwNistLibrarySearch
# TwNistLibraryQueryResult
# TwBruteForceCalibrate
# TwGetMassCalibInfo
# TwEncImsCorrelateProfile
# TwEncImsSharpenProfile (removed in TofDaq API 1.98)
# TwEncImsDenoiseProfile (removed in TofDaq API 1.98)
# TwEncImsCorrelateMultiProfiles
# TwEncImsCleanup
