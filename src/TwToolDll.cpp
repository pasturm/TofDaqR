#include <Rcpp.h>
using namespace Rcpp;
#include <TwToolDll.h>
#include "TofDaqR.h"

// convert R string to C string
char* RtoCstring(SEXP rstring) {
  std::string str = Rcpp::as<std::string>(rstring);
  char *cstring;
  cstring = R_alloc(str.length() + 1, sizeof(char));  // stringlength + 1 to account for null termination
  strcpy(cstring, str.c_str());  // copy rstring to cstring
  cstring[str.length()] = '\0';  // null terminate for safety
  return cstring;
}

// convert TwRetVal to String
StringVector TwRetValString(TwRetVal rv) {
  StringVector string(13);
  string[0] = "TwDaqRecNotRunning";
  string[1] = "TwAcquisitionActive";
  string[2] = "TwNoActiveAcquisition";
  string[3] = "TwFileNotFound";
  string[4] = "TwSuccess";
  string[5] = "TwError";
  string[6] = "TwOutOfBounds";
  string[7] = "TwNoData";
  string[8] = "TwTimeout";
  string[9] = "TwValueAdjusted";
  string[10] = "TwInvalidParameter";
  string[11] = "TwInvalidValue";
  string[12] = "TwAborted";
  StringVector result(1, string[(int)rv]);
  return result;
}

// FitSinglePeak ---------------------------------------------------------------
//' Performs a peak fit.
//'
//' \code{FitSinglePeak} performs a peak fit. Gaussian, Lorentzian and
//' Pseudo-Voigt peak shapes are supported.
//'
//' Peak parameters are optimized using a Levenberg-Marquardt algorithm as
//' implemented by the library \code{lmfit}. The following peak types are available:
//' \tabular{cccc}{
//' peak type index \tab peak function \tab symmetric \tab baseline \cr
//' -1 \tab Gaussian through highest point and its neighbours \tab yes \tab no \cr
//' 0 \tab Gaussian \tab yes \tab no  \cr
//' 1 \tab Lorentzian \tab yes \tab no \cr
//' 2 \tab Pseudo-Voigt \tab yes \tab no \cr
//' 3 \tab Gaussian \tab yes \tab yes \cr
//' 4 \tab Lorentzian \tab yes \tab yes \cr
//' 5 \tab Pseudo-Voigt \tab yes \tab yes \cr
//' 6 \tab Gaussian \tab no \tab no \cr
//' 7 \tab Lorentzian \tab no \tab no \cr
//' 8 \tab Pseudo-Voigt \tab no \tab no \cr
//' 9 \tab Gaussian \tab no \tab yes \cr
//' 10 \tab Lorentzian \tab no \tab yes \cr
//' 11 \tab Pseudo-Voigt \tab no \tab yes
//' } The function for peakType = -1 is not an actual peak fit but rather the
//' Gaussian function (symmetric, no baseline) defined by the maximum intensity
//' point and its neighbors. It is used to automatically generate guess values
//' for position, width and amplitude if all peak parameters are set to 0.
//'
//' @param yVals y axis data
//' @param xVals x axis data
//' @param peakType peak model to use.
//' @param blOffset Initial value of baseline offset at first data point.
//' @param blSlope Initial value of slope of baseline.
//' @param amplitude Initial value of peak amplitude.
//' @param fwhmLo Initial value of left peak width (full width at half maximum).
//' @param fwhmHi Initial value of right peak width (full width at half maximum).
//' @param peakPos Initial value of peak position.
//' @param mu Initial value of Gaussian/Lorentzian contribution.
//' @return List with peak parameters.
//'
//' @family Peak fitting functions
//' @export
// [[Rcpp::export]]
SEXP FitSinglePeak(NumericVector yVals, NumericVector xVals, int peakType = 0,
                   double blOffset = 0, double blSlope = 0,
                   double amplitude = 0, double fwhmLo = 0, double fwhmHi = 0,
                   double peakPos = 0, double mu = 0) {

  int nbrDataPoints = xVals.size();

  TwRetVal rv = TwFitSinglePeak(nbrDataPoints, &yVals[0], &xVals[0], peakType,
                                &blOffset, &blSlope, &amplitude, &fwhmLo,
                                &fwhmHi, &peakPos, &mu);

  if (rv != TwSuccess) {
    return TwRetValString(rv);
  }

  List result;
  result["blOffset"] = blOffset;
  result["blSlope"] = blSlope;
  result["amplitude"] = amplitude;
  result["fwhmLo"] = fwhmLo;
  result["fwhmHi"] = fwhmHi;
  result["peakPos"] = peakPos;
  result["mu"] = mu;
  return result;
}

// FitSinglePeak2 --------------------------------------------------------------
//' Performs a peak fit (Initial values as vector).
//'
//' \code{FitSinglePeak2} performs a peak fit. Gaussian, Lorentzian and
//' Pseudo-Voigt peak shapes are supported.
//'
//' Same as \code{\link{FitSinglePeak}}, but takes the initial values of the
//' fit parameter as a vector argument instead of explicit parameters.
//' Peak parameters are optimized using a Levenberg-Marquardt algorithm as
//' implemented by the library \code{lmfit}. The following peak types are available:
//' \tabular{cccc}{
//' peak type index \tab peak function \tab symmetric \tab baseline \cr
//' -1 \tab Gaussian through highest point and its neighbours \tab yes \tab no \cr
//' 0 \tab Gaussian \tab yes \tab no  \cr
//' 1 \tab Lorentzian \tab yes \tab no \cr
//' 2 \tab Pseudo-Voigt \tab yes \tab no \cr
//' 3 \tab Gaussian \tab yes \tab yes \cr
//' 4 \tab Lorentzian \tab yes \tab yes \cr
//' 5 \tab Pseudo-Voigt \tab yes \tab yes \cr
//' 6 \tab Gaussian \tab no \tab no \cr
//' 7 \tab Lorentzian \tab no \tab no \cr
//' 8 \tab Pseudo-Voigt \tab no \tab no \cr
//' 9 \tab Gaussian \tab no \tab yes \cr
//' 10 \tab Lorentzian \tab no \tab yes \cr
//' 11 \tab Pseudo-Voigt \tab no \tab yes
//' } The function for peakType = -1 is not an actual peak fit but rather the
//' Gaussian function (symmetric, no baseline) defined by the maximum intensity
//' point and its neighbors. It is used to automatically generate guess values
//' for position, width and amplitude if all peak parameters are set to 0.
//'
//' @param yVals y axis data
//' @param xVals x axis data
//' @param peakType peak model to use.
//' @param param Vector of initial values (blOffset, blSlope, amplitude, fwhmLo,
//' fwhmHi, peakPos, mu).
//' @return List with peak parameters.
//'
//' @family Peak fitting functions
//' @export
// [[Rcpp::export]]
SEXP FitSinglePeak2(NumericVector yVals, NumericVector xVals, int peakType = 0,
                    NumericVector param = NumericVector::create(0,0,0,0,0,0,0)) {

  int nbrDataPoints = xVals.size();

  TwRetVal rv = TwFitSinglePeak2(nbrDataPoints, &yVals[0], &xVals[0], peakType,
                                 &param[0]);
  if (rv != TwSuccess) {
    return TwRetValString(rv);
  }

  List result;
  result["blOffset"] = param[0];
  result["blSlope"] = param[1];
  result["amplitude"] = param[2];
  result["fwhmLo"] = param[3];
  result["fwhmHi"] = param[4];
  result["peakPos"] = param[5];
  result["mu"] = param[6];
  return result;
}

// EvalSinglePeak --------------------------------------------------------------
//' Evaluates a peak fit.
//'
//' \code{EvalSinglePeak} calculates the y-axis values for a given set of peak
//' parameters.
//'
//' @param xVals x axis data
//' @param blOffset Baseline offset at first data point.
//' @param blSlope Slope of baseline.
//' @param amplitude Peak amplitude.
//' @param fwhmLo Left peak width (full width at half maximum).
//' @param fwhmHi Right peak width (full width at half maximum).
//' @param peakPos Peak position.
//' @param mu Gaussian/Lorentzian contribution.
//' @return Vector with y axis data.
//'
//' @family Peak fitting functions
//' @export
// [[Rcpp::export]]
NumericVector EvalSinglePeak(NumericVector xVals, double blOffset = 0,
                             double blSlope = 0, double amplitude = 0,
                             double fwhmLo = 0, double fwhmHi = 0,
                             double peakPos = 0, double mu = 0) {

  int nbrDataPoints = xVals.size();

  double *param = new double[7];
  param[0] = blOffset;
  param[1] = blSlope;
  param[2] = amplitude;
  param[3] = fwhmLo;
  param[4] = fwhmHi;
  param[5] = peakPos;
  param[6] = mu;

  NumericVector yValsFit(nbrDataPoints);

  for (int j = 0; j<nbrDataPoints; ++j) {
    yValsFit[j] = TwEvalSinglePeak(xVals[j], param);
  }

  return yValsFit;
}

// GetMoleculeMass -------------------------------------------------------------
//' Calculates the mass/charge ratio of a molecule.
//'
//' \code{GetMoleculeMass} parses the molecular formula and returns the
//' mass/charge ratio of the molecule/atom/ion.
//'
//' The formula parser for the molecule string understands the (case sensitive)
//' element labels and round and curly brackets (square brackets are reserved for
//' isotope specification). Numbers (multipliers) have to follow either directly
//' the element symbol or a closing bracket. Specific isoptopes are specified in
//' square brackets with the mass number before the element label (e.g. [2H],
//' [14C] or [235U]). Charge indicators (+, -) have to be the last characters in
//' the formula string, multiple charges require multiple charge symbols (e.g.
//' doubly charged calcium is Ca++ and not Ca2+ which is the correct syntax for a
//' singly charged calcium dimer).
//'
//' @param molecule Molecule string.
//' @return Mass/charge ratio.
//'
//' @family Chemistry functions
//'
//' @examples
//' GetMoleculeMass("CO2")
//' GetMoleculeMass("CO2++")
//' GetMoleculeMass("[13C]O2+")
//' GetMoleculeMass("[13C][18O]2")
//' @export
// [[Rcpp::export]]
SEXP GetMoleculeMass(SEXP molecule) {

  char * cMolecule = RtoCstring(molecule);
  NumericVector mass(1);

  TwRetVal rv = TwGetMoleculeMass(cMolecule, &mass[0]);

  if (rv != TwSuccess) {
    return TwRetValString(rv);
  }
  return mass;
}

// GetIsotopePattern -----------------------------------------------------------
//' Calculates the isotope pattern of a molecule.
//'
//' \code{GetIsotopePattern} parses the molecular formula and returns the
//' isotope pattern (mass and abundance).
//'
//' The formula parser for the molecule string understands the (case sensitive)
//' element labels and round and curly brackets (square brackets are reserved for
//' isotope specification). Numbers (multipliers) have to follow either directly
//' the element symbol or a closing bracket. Specific isoptopes are specified in
//' square brackets with the mass number before the element label (e.g. [2H],
//' [14C] or [235U]). Charge indicators (+, -) have to be the last characters in
//' the formula string, multiple charges require multiple charge symbols (e.g.
//' doubly charged calcium is Ca++ and not Ca2+ which is the correct syntax for a
//' singly charged calcium dimer).
//'
//' @param molecule Molecule string.
//' @param abundanceLimit Absolute abundance limit for the generated pattern.
//' @return List with mass and abundace vectors.
//'
//' @family Chemistry functions
//'
//' @examples
//' GetIsotopePattern("CO2", 1e-5)
//' @export
// [[Rcpp::export]]
SEXP GetIsotopePattern(SEXP molecule, double abundanceLimit) {

  char * cMolecule = RtoCstring(molecule);

  int nbrIsotopes = 0;

  // get number of isotopes
  TwRetVal rv = TwGetIsotopePattern(cMolecule, abundanceLimit, &nbrIsotopes,
                                    NULL, NULL);

  if (rv != TwValueAdjusted) {
    return TwRetValString(rv);
  }

  NumericVector isoMass(nbrIsotopes);
  NumericVector isoAbundance(nbrIsotopes);

  rv = TwGetIsotopePattern(cMolecule, abundanceLimit, &nbrIsotopes, &isoMass[0],
                           &isoAbundance[0]);
  if (rv != TwSuccess) {
    return TwRetValString(rv);
  }

  List result;
  result["isoMass"] = isoMass;
  result["isoAbundance"] = isoAbundance;
  return result;
}

// GetIsotopePattern2 ----------------------------------------------------------
//' Calculates the isotope pattern of a molecule.
//'
//' \code{GetIsotopePattern2} parses the molecular formula and returns the
//' isotope pattern (mass and abundance).
//'
//' Same as \code{GetIsotopePattern} but using an exact algorithm and is
//' therefore suitable only for rather small molecules.
//'
//' The formula parser for the molecule string understands the (case sensitive)
//' element labels and round and curly brackets (square brackets are reserved for
//' isotope specification). Numbers (multipliers) have to follow either directly
//' the element symbol or a closing bracket. Specific isoptopes are specified in
//' square brackets with the mass number before the element label (e.g. [2H],
//' [14C] or [235U]). Charge indicators (+, -) have to be the last characters in
//' the formula string, multiple charges require multiple charge symbols (e.g.
//' doubly charged calcium is Ca++ and not Ca2+ which is the correct syntax for a
//' singly charged calcium dimer).
//'
//' @param molecule Molecule string.
//' @param abundanceLimit Absolute abundance limit for the generated pattern.
//' @return List with mass and abundace vectors.
//'
//' @family Chemistry functions
//'
//' @examples
//' GetIsotopePattern2("CO2", 1e-5)
//' @export
// [[Rcpp::export]]
SEXP GetIsotopePattern2(SEXP molecule, double abundanceLimit) {

  char * cMolecule = RtoCstring(molecule);

  int nbrIsotopes = 0;

  // get number of isotopes
  TwRetVal rv = TwGetIsotopePattern2(cMolecule, abundanceLimit, &nbrIsotopes,
                                    NULL, NULL);

  if (rv != TwValueAdjusted) {
    return TwRetValString(rv);
  }

  NumericVector isoMass(nbrIsotopes);
  NumericVector isoAbundance(nbrIsotopes);

  rv = TwGetIsotopePattern2(cMolecule, abundanceLimit, &nbrIsotopes, &isoMass[0],
                            &isoAbundance[0]);
  if (rv != TwSuccess) {
    return TwRetValString(rv);
  }

  List result;
  result["isoMass"] = isoMass;
  result["isoAbundance"] = isoAbundance;
  return result;
}

// Tof2Mass --------------------------------------------------------------------
//' Converts from sample index to mass/charge.
//'
//' \code{Tof2Mass} converts from sample index to mass/charge.
//'
//' \tabular{cl}{
//' massCalibMode \tab Mass calibration function \cr
//' 0 \tab \eqn{i = p_1 \sqrt(m) + p_2} \cr
//' 1 \tab \eqn{i = p_1/\sqrt(m) + p_2} \cr
//' 2 \tab \eqn{i = p_1 m^{p_3} + p_2} \cr
//' 3 \tab \eqn{i = p_1 \sqrt(m) + p_2 + p_3 (m - p_4)^2} \cr
//' 4 \tab \eqn{i = p_1 \sqrt(m) + p_2 + p_3 m^2 + p_4 m + p_5} \cr
//' 5 \tab \eqn{m = p_1 i^2 + p_2 i + p_3}
//' }
//' Note: Modes 3 and 4 are flawed. Don't use them. In mode 3 the fit does not
//' converge well, because of a bug (parameters not correctly initialized).
//' Mode 4 is two sequential fits, first mode 0, then a quadratic fit to the
//' residuals, which is an inferior implementation of mode 3. Mode 1 is for FTMS
//' data.
//'
//' @param tofSample Vector of sample indices to convert.
//' @param massCalibMode Mass calibration function to use. See below.
//' @param p Vector containing the calibration parameters (number depends on
//' \code{MassCalibMode}, see below).
//' @return Mass/charge values.
//'
//' @seealso \code{\link{Mass2Tof}}
//'
//' @examples
//' Tof2Mass(100000, massCalibMode = 0, p = c(3,5))
//' @export
// [[Rcpp::export]]
NumericVector Tof2Mass(NumericVector tofSample, int massCalibMode, NumericVector p) {

  int nbrSamples = tofSample.size();
  NumericVector mass(nbrSamples);

  for (int i = 0; i < nbrSamples; ++i) {
    mass[i] = TwTof2Mass(tofSample[i], massCalibMode, &p[0]);
  }

  return mass;
}

// Mass2Tof --------------------------------------------------------------------
//' Converts from mass/charge to sample index.
//'
//' \code{Mass2Tof} converts from mass/charge to sample index.
//'
//' \tabular{cl}{
//' massCalibMode \tab Mass calibration function \cr
//' 0 \tab \eqn{i = p_1 \sqrt(m) + p_2} \cr
//' 1 \tab \eqn{i = p_1/\sqrt(m) + p_2} \cr
//' 2 \tab \eqn{i = p_1 m^{p_3} + p_2} \cr
//' 3 \tab \eqn{i = p_1 \sqrt(m) + p_2 + p_3 (m - p_4)^2} \cr
//' 4 \tab \eqn{i = p_1 \sqrt(m) + p_2 + p_3 m^2 + p_4 m + p_5} \cr
//' 5 \tab \eqn{m = p_1 i^2 + p_2 i + p_3}
//' }
//' Note: Modes 3 and 4 are flawed. Don't use them. In mode 3 the fit does not
//' converge well, because of a bug (parameters not correctly initialized).
//' Mode 4 is two sequential fits, first mode 0, then a quadratic fit to the
//' residuals, which is an inferior implementation of mode 3. Mode 1 is for FTMS
//' data.
//'
//' @param mass Vector of mass/charge values to convert.
//' @param massCalibMode Mass calibration function to use. See below.
//' @param p Vector containing the calibration parameters (number depends on
//' \code{MassCalibMode}, see below).
//' @return Sample indices.
//'
//' @seealso \code{\link{Tof2Mass}}
//'
//' @examples
//' Mass2Tof(100, massCalibMode = 0, p = c(3,5))
//' @export
// [[Rcpp::export]]
NumericVector Mass2Tof(NumericVector mass, int massCalibMode, NumericVector p) {

  int nbrSamples = mass.size();
  NumericVector tofSample(nbrSamples);

  for (int i = 0; i < nbrSamples; ++i) {
    tofSample[i] = TwMass2Tof(mass[i], massCalibMode,  &p[0]);
  }

  return tofSample;
}

// MassCalibrate ---------------------------------------------------------------
//'  Performs a mass calibration.
//'
//' \code{MassCalibrate} performs a mass calibration for a list of
//' mass/sample index/weight values and for a given calibration mode.
//'
//' \tabular{cl}{
//' massCalibMode \tab Mass calibration function \cr
//' 0 \tab \eqn{i = p_1 \sqrt(m) + p_2} \cr
//' 1 \tab \eqn{i = p_1/\sqrt(m) + p_2} \cr
//' 2 \tab \eqn{i = p_1 m^{p_3} + p_2} \cr
//' 3 \tab \eqn{i = p_1 \sqrt(m) + p_2 + p_3 (m - p_4)^2} \cr
//' 4 \tab \eqn{i = p_1 \sqrt(m) + p_2 + p_3 m^2 + p_4 m + p_5} \cr
//' 5 \tab \eqn{m = p_1 i^2 + p_2 i + p_3}
//' }
//' Note: Modes 3 and 4 are flawed. Don't use them. In mode 3 the fit does not
//' converge well, because of a bug (parameters not correctly initialized).
//' Mode 4 is two sequential fits, first mode 0, then a quadratic fit to the
//' residuals, which is an inferior implementation of mode 3. Mode 1 is for FTMS
//' data.
//'
//' @param massCalibMode Mass calibration function to use. See below.
//' @param mass Vector of mass/charge values.
//' @param idx Vector of sample indices (of same length as \code{mass}).
//' @param weight Vector of weights (if \code{NULL} (default) all weights are set equal).
//' @return Vector of calibration parameters.
//' @export
// [[Rcpp::export]]
SEXP MassCalibrate(int massCalibMode, NumericVector mass, NumericVector tof,
                   Nullable<NumericVector> weight = R_NilValue) {

  int nbrPoints = mass.size();

  double * p_weight;
  if (weight.isNotNull()) {
    NumericVector x(weight);
    p_weight = &x[0];
  } else {
    p_weight = nullptr;
  }

  int nbrParams;

  // get the number of parameters
  char * description = new char[64];
  TwRetVal rv = TwGetMassCalibInfo(massCalibMode, description, &nbrParams);
  if (rv != TwSuccess) {
    return TwRetValString(rv);
  }

  NumericVector param(nbrParams);

  TwMassCalibrate(massCalibMode, nbrPoints, &mass[0], &tof[0], p_weight, &nbrParams, &param[0], NULL, NULL);

  return param;
}

// SiInitializeHistograms ------------------------------------------------------
//' Initializes the single ion histogramming.
//'
//' \code{SiInitializeHistograms} initializes the single ion histogramming.
//'
//' This function must be called before calling the functions
//' \code{\link{SiSetProcessingOptions}}, \code{\link{SiProcessSpectrum}}, \code{\link{SiGetHistogram}},
//' \code{\link{SiResetHistograms}} or \code{\link{SiCleanup}}. Each mass range specified by
//' \code{loMass} and \code{hiMass} elements is associated with a spectrum type that allows to
//' get separate statistics for multi-spectrum acquisitions (bipolar or pulsed
//' experiments).
//'
//' @param loMass Vector of the lower borders of the mass ranges.
//' @param hiMass Vector of the upper borders of the mass ranges.
//' @param specType Vector of spectrum type indices (non-negative integers). If
//' specType is \code{NULL}, all mass ranges get a default spectrum type of 0.
//'
//' @family Single ion histogramming functions.
//' @export
// [[Rcpp::export]]
StringVector SiInitializeHistograms(NumericVector loMass, NumericVector hiMass,
                            Nullable<IntegerVector> specType = R_NilValue) {

  int nbrHist = loMass.size();
  int * p_specType;
  if (specType.isNotNull()) {
    IntegerVector x(specType);
    p_specType = &x[0];
  } else {
    p_specType = nullptr;
  }

  TwRetVal rv = TwSiInitializeHistograms(nbrHist, &loMass[0], &hiMass[0], p_specType);

  return TwRetValString(rv);
}

// SiSetProcessingOptions ------------------------------------------------------
//' Sets processing options for each spectrum type.
//'
//' \code{SiSetProcessingOptions} sets processing options for each spectrum type.
//'
//' Options:
//' \tabular{llc}{
//' Name \tab Description \tab Default value \cr
//' \code{MassCalibMode} \tab Mass calibration mode in use (see \code{\link{MassCalibrate}}). \tab 0 \cr
//' \code{MassCalibParamn} \tab Mass calibration parameters n =
//' 1..number of calibration parameters for the given mass calibration mode. \tab \code{c(1000, 0)} \cr
//' \code{FullScale} \tab Full digitizing range of the ADC in the same units as
//' the spectra to be analyzed (typically mV). \tab 500 \cr
//' \code{NbrBits} \tab ADC resolution (8 for AP240, 14 for ADQ114 etc.). \tab 8 \cr
//' \code{SampleInterval} \tab Sample interval in ns. \tab 1 \cr
//' \code{PreampGain} \tab Gain of external preamp. \tab 1 \cr
//' \code{PreSamples} \tab Number of samples before a threshold crosser taken into account. \tab 0 \cr
//' \code{PostSamples} \tab Number of samples after a negative threshold crosser taken into account. \tab 0 \cr
//' \code{BaselineAndThresholdFromData} \tab If >0 the baseline and threshold
//' values will be determined based on \code{NbrStdDevs} for every processed
//' spectrum. If >1.5 baseline noise is determined from a fit to a histogram of
//' all data (instead of from the standard deviation of all data). This makes the
//' noise determination more robust when real peaks are present in the spectrum. \tab 0 \cr
//' \code{NbrStdDevs} \tab Number of standard deviations of baseline noise that
//' defines the threshold. Only relevant if \code{BaselineAndThresholdFromData>0}. \tab 6 \cr
//' \code{Baseline} \tab Baseline value used for calculation of intensities. Has
//' no meaning if \code{BaselineAndThresholdFromData>0}. \tab 5 \cr
//' \code{Threshold} \tab Threshold value used for calculation of intensities.
//' Has no meaning if \code{BaselineAndThresholdFromData!=0}. \tab 8 \cr
//' \code{NegativeSignal} \tab Indicates peak polarity with respect to the baseline. \tab 0 (=FALSE) \cr
//' \code{BaselineAndThresholdInCodes} \tab Indicates whether the values in \code{Baseline}
//' and \code{Threshold} are interpreted as ADC codes or mV. \tab 1 (=TRUE)
//' }
//'
//' @param option Option to set (see below for a list of valid option strings).
//' @param value Value to set the for the given option.
//' @param specType Spectrum type index. -1 is a wildcard for all spectrum types.
//'
//' @family Single ion histogramming functions.
//' @export
// [[Rcpp::export]]
StringVector SiSetProcessingOptions(SEXP option, double value, int specType) {

  char * coption = RtoCstring(option);

  TwRetVal rv = TwSiSetProcessingOptions(coption, value, specType);

  return TwRetValString(rv);
}

// SiProcessSpectrum -----------------------------------------------------------
//' Processes a spectrum.
//'
//' \code{SiProcessSpectrum} processes a spectrum according to the options set for
//' it's spectrum type.
//'
//' See \code{\link{SiProcessSpectrumFromShMem}} for another variant of this
//' function, where the spectrum is directly read from shared memory.
//'
//' @param spectrum Vector holding the spectrum to process.
//' @param specType Spectrum type index (non-negative integer).
//' @return A list with the baseline and threshold value.
//'
//' @family Single ion histogramming functions.
//' @export
// [[Rcpp::export]]
SEXP SiProcessSpectrum(NumericVector spectrum, int specType) {

  int nbrSamples = spectrum.size();

  std::vector<float> fspectrum = Rcpp::as<std::vector<float> >(spectrum);

  float blFromData;
  float thrFromData;

  TwRetVal rv = TwSiProcessSpectrum(&fspectrum[0], nbrSamples, specType,
                                    &blFromData, &thrFromData);

  if (rv != TwSuccess) {
    return TwRetValString(rv);
  }

  List result;
  result["baseline"] = wrap(blFromData);
  result["threshold"] = wrap(thrFromData);

  return result;
}

// SiGetHistogram --------------------------------------------------------------
//' Gets a histogram of the single ion intensities.
//'
//' \code{SiGetHistogram} gets a histogram of the single ion intensities for a
//' mass range defined with \code{\link{SiInitializeHistograms}}.
//'
//' Note: R crashes if \code{histogramIndex} is set to max(histogramIndex)+1
//' (API bug).
//'
//' @param histogramIndex Index (0-based) of the histogram. It corresponds to the
//' mass range defined with \code{\link{SiInitializeHistograms}}.
//' @return A list with the intensities (histogram x-values), counts (histogram
//' y-values), the number of spectra that were processed for this histogram and
//' the mean histogram value i.e. sum(intensity[i]*counts[i])/sum(counts[i]).
//'
//' @family Single ion histogramming functions.
//' @export
// [[Rcpp::export]]
SEXP SiGetHistogram(int histogramIndex) {

  unsigned int arrayLength;
  unsigned int spectrumCount;
  double meanValue;

  //get length of array
  TwRetVal rv = TwSiGetHistogram(histogramIndex, NULL, NULL, &arrayLength,
                                 &spectrumCount, &meanValue);
  if (rv != TwValueAdjusted) {
    return TwRetValString(rv);
  }

  NumericVector intensity(arrayLength);
  NumericVector counts(arrayLength);

  std::vector<float> fintensity = Rcpp::as<std::vector<float> >(intensity);
  std::vector<unsigned int> fcounts = Rcpp::as<std::vector<unsigned int> >(counts);

  rv = TwSiGetHistogram(histogramIndex, &fintensity[0], &fcounts[0],
                        &arrayLength, &spectrumCount, &meanValue);
  if (rv != TwSuccess) {
    return TwRetValString(rv);
  }

  List result;
  result["intensity"] = wrap(fintensity);
  result["counts"] = wrap(fcounts);
  result["spectrumCount"] = spectrumCount;
  result["meanValue"] = meanValue;

  return result;
}

// SiGetSumHistogram -----------------------------------------------------------
//' Gets a sum histogram of the single ion intensities.
//'
//' \code{SiGetSumHistogram} gets a histogram of the single ion intensities, which
//' is a sum over all histograms of a given \code{specType} within the rate and
//' mass range as specified by \code{minMass}, \code{maxMass}, \code{minRate} and \code{maxRate}.
//'
//' @param specType Spectrum type index (non-negative integer).
//' @param minMass Minimum mass for histogram filtering.
//' @param maxMass Maximum mass for histogram filtering.
//' @param minRate Minimum event count rate for histogram filtering.
//' @param maxRate Maximum event count rate for histogram filtering.
//' @return A list with the intensities (histogram x-values), counts (histogram
//' y-values), the number of spectra that were processed for this histogram and
//' the mean histogram value i.e. sum(intensity[i]*counts[i])/sum(counts[i]).
//'
//' @family Single ion histogramming functions.
//' @export
// [[Rcpp::export]]
SEXP SiGetSumHistogram(int specType, double minMass, double maxMass,
                       double minRate, double maxRate) {

  unsigned int arrayLength;
  unsigned int spectrumCount;
  double meanValue;

  //get length of array
  TwRetVal rv = TwSiGetSumHistogram(specType, NULL, NULL, &arrayLength,
                                    &spectrumCount, &meanValue, minMass, maxMass,
                                    minRate, maxRate);
  if (rv != TwValueAdjusted) {
    return TwRetValString(rv);
  }

  NumericVector intensity(arrayLength);
  NumericVector counts(arrayLength);

  std::vector<float> fintensity = Rcpp::as<std::vector<float> >(intensity);
  std::vector<unsigned int> fcounts = Rcpp::as<std::vector<unsigned int> >(counts);

  rv = TwSiGetSumHistogram(specType, &fintensity[0], &fcounts[0], &arrayLength,
                           &spectrumCount, &meanValue, minMass, maxMass,
                           minRate, maxRate);
  if (rv != TwSuccess) {
    return TwRetValString(rv);
  }

  List result;
  result["intensity"] = wrap(fintensity);
  result["counts"] = wrap(fcounts);
  result["spectrumCount"] = spectrumCount;
  result["meanValue"] = meanValue;

  return result;
}

// SiResetHistograms -----------------------------------------------------------
//' Resets all histograms and spectrum counters to zero.
//'
//' \code{SiResetHistograms} resets all histograms and spectrum counters to zero.
//'
//' @family Single ion histogramming functions.
//' @export
// [[Rcpp::export]]
StringVector SiResetHistograms() {

  TwRetVal rv = TwSiResetHistograms();

  return TwRetValString(rv);
}

// SiCleanup -------------------------------------------------------------------
//' Cleans up the state in the DLL.
//'
//' \code{SiCleanup} cleans up the state in the DLL. After a call to this
//' function \code{\link{SiInitializeHistograms}} can be called again.
//'
//' @family Single ion histogramming functions.
//' @export
// [[Rcpp::export]]
StringVector SiCleanup() {

  TwRetVal rv = TwSiCleanup();

  return TwRetValString(rv);
}

// SiFitPhd --------------------------------------------------------------------
//' Fits a (slightly modified) log-normal distribution to the histogram.
//'
//' \code{SiFitPhd} fits a (slightly modified) log-normal distribution to the
//' histogram data.
//'
//' The following equation with the parameters A, w, x0 and xc is fitted: \cr
//' \eqn{A/(\sqrt(2*pi)*w*(x-x0))*exp(-(log((x-x0)/xc)^2)/(2*w^2))}
//'
//' @param intensity Vector of intensities (histogram x-axis data).
//' @param counts Vector of counts (histogram y-axis data).
//' @return A list with the FWHM of the distribution, the position of the fitted
//' maximum and a vector with the values of the four fitting parameters.
//'
//' @family Single ion histogramming functions.
//' @export
// [[Rcpp::export]]
SEXP SiFitPhd(NumericVector intensity, NumericVector counts) {

  int nbrPoints = intensity.size();

  NumericVector fwhm(nbrPoints);
  NumericVector a(nbrPoints);
  NumericVector par(nbrPoints);

  TwRetVal rv = TwSiFitPhd(nbrPoints, &intensity[0], &counts[0], &fwhm[0],
                           &a[0], &par[0]);
  if (rv != TwSuccess) {
    return TwRetValString(rv);
  }

  List result;
  result["fwhm"] = fwhm;
  result["a"] = a;
  result["par"] = par;

  return result;
}

// SiEvalPhd -------------------------------------------------------------------
//' Evaluates the fitted single ion distribution.
//'
//' \code{SiEvalPhd} evaluates the fitted single ion distribution.
//'
//' @param par Vector of fitted parameter values.
//' @param intensity Vector of intensities (histogram x-axis data).
//' @return Vector with y-axis data.
//'
//' @family Single ion histogramming functions.
//' @export
// [[Rcpp::export]]
NumericVector SiEvalPhd(NumericVector par, NumericVector intensity) {

  int nbrPoints = intensity.size();

  NumericVector yValsFit(nbrPoints);

  for (int j = 0; j<nbrPoints; ++j) {
    yValsFit[j] = TwSiEvalPhd(&par[0], intensity[j]);
  }

  return yValsFit;
}

// SiFitRateFromPhd ------------------------------------------------------------
//' Fits an event rate to a multi-ion histogram.
//'
//' \code{SiFitRateFromPhd} takes a fitted single ion distribution as input
//' and fits an event rate to a multi-ion histogram (intensity and counts)
//' assuming poisson distribution of n-ion events.
//'
//' @param intensity Vector of intensities (histogram x-axis data).
//' @param counts Vector of counts (histogram y-axis data).
//' @param siPar Vector with fitted parameter values.
//' @return A list with the fitted rate (ions/extraction) and vector with the
//' fitted multi-ion distribution.
//'
//' @family Single ion histogramming functions.
//' @export
// [[Rcpp::export]]
SEXP SiFitRateFromPhd(NumericVector intensity, NumericVector counts,
                      NumericVector siPar) {

  int nbrPoints = intensity.size();

  double rate;
  NumericVector fitCounts(nbrPoints);

  TwRetVal rv = TwSiFitRateFromPhd(nbrPoints, &intensity[0], &counts[0],
                                   &siPar[0], &rate, &fitCounts[0], 0, NULL);
  if (rv != TwSuccess) {
    return TwRetValString(rv);
  }

  List result;
  result["rate"] = rate;
  result["fitCounts"] = fitCounts;

  return result;
}

#ifdef _WIN32
// FindTpsIp -------------------------------------------------------------------
//' Gets IP address of TPS2.
//'
//' \code{FindTpsIp} listens for UDP packets that TPS2 broadcast and returns
//' the IP of the TPS2.
//'
//' Note that executing this function makes your program to a UDP server and
//' Windows firewall (or other personal firewall software) will query for
//' permission.
//'
//' @param TpsSerial Serial number of TPS2.
//' @param timeout Timeout in ms to wait for the UDP packet.
//' @return String of IP address.
//'
//' @examples
//' FindTpsIp("910.33.0316", 500)
//' @export
// [[Rcpp::export]]
SEXP FindTpsIp(SEXP TpsSerial, int timeout) {

  char * cTpsSerial = RtoCstring(TpsSerial);

  int hostStrLen = 15;

  char *buffer = new char[15];

  TwRetVal rv = TwFindTpsIp(cTpsSerial, timeout, &hostStrLen, buffer);

  if (rv != TwSuccess) {
    return TwRetValString(rv);
  }

  std::string str(buffer);

  return wrap(str);
}
#endif

// Not implemented: TwTranslateReturnValue -------------------------------------
// Not implemented: TwMultiPeakFit ---------------------------------------------
// Not implemented: TwEvalMultiPeak --------------------------------------------
// Not implemented: TwFitResolution --------------------------------------------
// Not implemented: TwEvalResolution -------------------------------------------
// Not implemented: TwDecomposeMass --------------------------------------------
// Not implemented: TwGetComposition -------------------------------------------
// Not implemented: TwNistLibrarySearch ----------------------------------------
// Not implemented: TwNistLibraryQueryResult -----------------------------------
// Not implemented: TwBruteForceCalibrate --------------------------------------
// Not implemented: TwGetMassCalibInfo -----------------------------------------
// Not implemented: TwEncImsCorrelateProfile -----------------------------------
// Not implemented: TwEncImsCorrelateMultiProfiles -----------------------------
// Not implemented: TwEncImsCleanup --------------------------------------------
// Not implemented: TwGenerateImageData ----------------------------------------
// Not implemented: TwImagePaletteOffset2Value ---------------------------------
// Not implemented: TwImageValue2PaletteOffset ---------------------------------
// Not implemented: TwImageGetPaletteRGB ---------------------------------------
