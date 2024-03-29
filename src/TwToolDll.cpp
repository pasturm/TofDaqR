#include <Rcpp.h>
using namespace Rcpp;
#include <TwToolDll.h>
#include "TofDaqR.h"

// convert std::string to char*
char* StringToChar(std::string str) {
  char *cstring;
  cstring = R_alloc(str.length() + 1, sizeof(char));  // + 1 to account for null termination
  strcpy(cstring, str.c_str());  // copy string to cstring
  return cstring;
}

// convert TwRetVal to std::string
std::string TranslateReturnValue(TwRetVal rv) {
  char *retval = TwTranslateReturnValue(rv);
  std::string str(retval);
  return str;
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
//' @param yVals y axis data.
//' @param xVals x axis data.
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
List FitSinglePeak(NumericVector yVals, NumericVector xVals, int peakType = 0,
                   double blOffset = 0, double blSlope = 0,
                   double amplitude = 0, double fwhmLo = 0, double fwhmHi = 0,
                   double peakPos = 0, double mu = 0) {

  int nbrDataPoints = xVals.size();

  TwRetVal rv = TwFitSinglePeak(nbrDataPoints, &yVals[0], &xVals[0], peakType,
                                &blOffset, &blSlope, &amplitude, &fwhmLo,
                                &fwhmHi, &peakPos, &mu);

  if (rv != TwSuccess) {
    stop(TranslateReturnValue(rv));
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
//' @param yVals y axis data.
//' @param xVals x axis data.
//' @param peakType peak model to use.
//' @param param Vector of initial values (blOffset, blSlope, amplitude, fwhmLo,
//' fwhmHi, peakPos, mu).
//' @return List with peak parameters.
//'
//' @family Peak fitting functions
//' @export
// [[Rcpp::export]]
List FitSinglePeak2(NumericVector yVals, NumericVector xVals, int peakType = 0,
                    NumericVector param = NumericVector::create(0,0,0,0,0,0,0)) {

  int nbrDataPoints = xVals.size();

  TwRetVal rv = TwFitSinglePeak2(nbrDataPoints, &yVals[0], &xVals[0], peakType,
                                 &param[0]);
  if (rv != TwSuccess) {
    stop(TranslateReturnValue(rv));
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
//' Calculates the y-axis values for a given set of peak parameters.
//'
//' \code{EvalSinglePeak} calculates the y-axis values for a given set of peak
//' parameters.
//'
//' @param xVals x axis data.
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
  // auto *param{ new double[7]{ blOffset, blSlope, amplitude, fwhmLo, fwhmHi, peakPos, mu } };  // the C++11 way

  NumericVector yValsFit(nbrDataPoints);

  for (int j = 0; j<nbrDataPoints; ++j) {
    yValsFit[j] = TwEvalSinglePeak(xVals[j], param);
  }

  delete[] param;

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
double GetMoleculeMass(std::string molecule) {

  char *cMolecule = StringToChar(molecule);

  double mass;

  TwRetVal rv = TwGetMoleculeMass(cMolecule, &mass);

  if (rv != TwSuccess) {
    stop(TranslateReturnValue(rv));
  }
  return mass;
}

// MultiPeakFit ----------------------------------------------------------------
//' Performs a multi-peak fit.
//'
//' \code{MultiPeakFit} performs a multi-peak fit for (partially) overlapping
//' peaks. All peaks in a multi-peak cluster share common peak parameters except
//' amplitude and position. Various options allow for a more or less constrained
//' fit (see description of options below).
//'
//' The peak positions are not optimized, but the common mass shift parameter in
//' \code{commonPar} allows for minimal refinement of the peak positions due to
//' imperfect mass calibration.
//'
//' \code{options} is a list containing:
//' \tabular{rcl}{
//' peakModel \tab = \tab 0 (Gauss) or 1 (Lorentz) or 2 (Pseudo-Voigt) \cr
//' asymmetric \tab = \tab 0 (symmetric peak model) or 1 (asymmetric peak model) \cr
//' baseline \tab = \tab 0 (no optimization of baseline parameters) or 1 (optimize baseline parameters) \cr
//' width \tab = \tab 0 (no optimization of width parameters) or 1 (optimize width parameters) \cr
//' peakShape \tab = \tab 0 (no optimization of peak shape (mu) parameter) or 1 (optimize peak shape (mu) parameter (applies only to pseudo-Voigt)) \cr
//' massShift \tab = \tab 0 (no optimization of common mass shift parameters) or 1 (optimize common mass shift parameters) \cr
//' amplitude \tab = \tab 0 (no constraint on amplitudes) or 1 (constrain sum of baseline and all peaks to total intensity in spectrum)
//' } Note that if a given parameter is not activated for optimization, the
//' supplied (guess) values are still used (e.g. specify a known baseline
//' without optimizing the parameters).
//'
//' @param dataX x axis data.
//' @param dataY y axis data.
//' @param mass  Vector of the known positions of the peaks.
//' @param intensity  Vector of initial guesses for intensities. If all values
//' are 0 the guess values for peak intensities are generated automatically.
//' @param commonPar Vector of guess values of common peak parameters. \code{commonPar}
//' has 6 elements: [1] offset of common baseline, [2] slope of common baseline,
//' [3] left FWHM, [4] right FWHM (same as left FWHM for symmetric peaks),
//' [5] shape parameter mu (applies only for peak model Pseudo-Voigt),
//' [6] common mass shift of peaks. \code{options} determines which of the
//' common parameters are optimized.
//' @param options List of peak model and optimization options (see Details).
//'
//' @return List with the optimized intensities and common peak parameters.
//'
//' @family Peak fitting functions
//' @export
// [[Rcpp::export]]
List MultiPeakFit(NumericVector dataX, NumericVector dataY, NumericVector mass,
                  NumericVector intensity, NumericVector commonPar, List options) {

  int nbrDataPoints = dataX.size();
  int nbrPeaks = mass.size();
  if (nbrPeaks != intensity.size()) {
    stop("mass and intensity must be the same length.");
  }
  int opt = 0;
  opt = (int)options["peakModel"] | ((int)options["asymmetric"] << 2) |
    ((int)options["baseline"] << 3) | ((int)options["width"] << 4) |
    ((int)options["peakShape"] << 5) | ((int)options["massShift"] << 6) |
    ((int)options["amplitude"] << 7);

  TwRetVal rv = TwMultiPeakFit(nbrDataPoints, &dataX[0], &dataY[0], nbrPeaks,
                               &mass[0], &intensity[0], &commonPar[0], opt);

  if (rv != TwSuccess) {
    stop(TranslateReturnValue(rv));
  }

  List result;
  result["intensity"] = intensity;
  result["commonPar"] = commonPar;
  return result;
}

// EvalMultiPeak ---------------------------------------------------------------
//' Calculates the y-axis values for a given set of multi-peak parameters.
//'
//' \code{EvalMultiPeak} calculates the y-axis values for a given set of
//' multi-peak parameters.
//'
//' @param dataX x axis data.
//' @param mass  Vector of the peak positions.
//' @param intensity  Vector of (fitted) intensities.
//' @param commonPar Vector of (fitted) values of common peak parameters. \code{commonPar}
//' has 6 elements: [1] offset of common baseline, [2] slope of common baseline,
//' [3] left FWHM, [4] right FWHM, [5] shape parameter mu,
//' [6] common mass shift of peaks.
//'
//' @return Vector with y axis data.
//'
//' @family Peak fitting functions
//' @export
// [[Rcpp::export]]
NumericVector EvalMultiPeak(NumericVector dataX, NumericVector mass,
                            NumericVector intensity, NumericVector commonPar) {

  int nbrDataPoints = dataX.size();
  int nbrPeaks = mass.size();
  if (nbrPeaks != intensity.size()) {
    stop("mass and intensity must be the same length.");
  }

  NumericVector yValsFit(nbrDataPoints);

  for (int j = 0; j<nbrDataPoints; ++j) {
    yValsFit[j] = TwEvalMultiPeak(dataX[j], nbrPeaks, &mass[0], &intensity[0],
                                  &commonPar[0]);
  }

  return yValsFit;
}

// FitResolution ---------------------------------------------------------------
//' Fits mass versus resolution values to an empirical function.
//'
//' \code{FitResolution} fits mass versus resolution values to an empirical
//' function (see Details).
//'
//' An empirical function that describes the mass resolution as a function of
//' mass is: R(m) = R0 - R0/(1 + exp((m -m0)/dm)), where R0 is the nominal mass
//' resolution, m0 is the mass at which the resolution is R0/2 and dm is a slope
//' parameter.
//'
//' Reasonable initial guess values for R0, m0 and dm need to be provided.
//'
//' @param mass Vector of mass values.
//' @param resolution Vector of resolution values.
//' @param R0 Initial guess value of nominal mass resolution R0.
//' @param m0 Initial guess value of m0 (mass where R(m) = 0.5*R0).
//' @param dm Initial guess value of slope parameter dm.
//' @return List with the fitted parameters R0, m0 and dm.
//'
//' @seealso \code{\link{EvalResolution}}
//'
//' @export
// [[Rcpp::export]]
List FitResolution(NumericVector mass, NumericVector resolution, double R0,
                   double m0, double dm) {

  int nbrPoints = mass.size();

  TwRetVal rv = TwFitResolution(nbrPoints, &mass[0], &resolution[0], &R0, &m0, &dm);

  if (rv != TwSuccess) {
    stop(TranslateReturnValue(rv));
  }

  List result;
  result["R0"] = R0;
  result["m0"] = m0;
  result["dm"] = dm;
  return result;
}

// EvalResolution --------------------------------------------------------------
//' Evaluates the fitted resolution function for given mass values.
//'
//' \code{EvalResolution} evaluates the fitted resolution function for given
//' mass values.
//'
//' An empirical function that describes the mass resolution of TOFs as a
//' function of mass is: R(m) = R0 - R0/(1 + exp((m -m0)/dm)), where
//' R0 is the nominal mass resolution, m0 is the mass at which the resolution
//' is R0/2 and dm is a slope parameter.
//'
//' @param R0 Nominal mass resolution.
//' @param m0 Mass where R(m) = 0.5*R0.
//' @param dm Slope parameter.
//' @param mass Vector of mass values.
//' @return Vector with resolution values.
//'
//' @seealso \code{\link{FitResolution}}
//'
//' @export
// [[Rcpp::export]]
NumericVector EvalResolution(double R0, double m0, double dm, NumericVector mass) {

  int nbrPoints = mass.size();

  NumericVector yValsFit(nbrPoints);

  for (int j = 0; j<nbrPoints; ++j) {
    yValsFit[j] = TwEvalResolution(R0, m0, dm, mass[j]);
  }

  return yValsFit;
}

// GetIsotopePattern -----------------------------------------------------------
//' Calculates the isotope pattern of a molecule.
//'
//' \code{GetIsotopePattern} parses the molecular formula and returns the
//' isotope pattern (mass and abundance).
//'
//' The isotope pattern returned by this function may be quite inaccurate. A more
//' accurate pattern is obtained with \code{GetIsotopePattern2}, which however
//' only works for small molecules. See also
//' \code{\link[enviPat:isopattern]{enviPat::isopattern()}}
//' for an alternative function to calculate accurate isotope patterns.
//'
//' Note that abundances are returned relative to the most abundant peak of the isotope
//' pattern. The \code{abundanceLimit}, however, is an absolute abundance.
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
//' @return List with mass and abundace vectors. Abundances are normalized
//' relative to the most abundant peak of the isotope pattern.
//'
//' @family Chemistry functions
//'
//' @examples
//' GetIsotopePattern("CO2", 1e-5)
//' @export
// [[Rcpp::export]]
List GetIsotopePattern(std::string molecule, double abundanceLimit) {

  char *cMolecule = StringToChar(molecule);

  int nbrIsotopes = 0;

  // get number of isotopes
  TwRetVal rv = TwGetIsotopePattern(cMolecule, abundanceLimit, &nbrIsotopes,
                                    NULL, NULL);

  if (rv != TwValueAdjusted) {
    stop(TranslateReturnValue(rv));
  }

  NumericVector isoMass(nbrIsotopes);
  NumericVector isoAbundance(nbrIsotopes);

  rv = TwGetIsotopePattern(cMolecule, abundanceLimit, &nbrIsotopes, &isoMass[0],
                           &isoAbundance[0]);
  if (rv != TwSuccess) {
    stop(TranslateReturnValue(rv));
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
//' \code{GetIsotopePattern2} is the same as \code{GetIsotopePattern} but uses an
//' exact algorithm and is suitable only for small molecules.
//' See also \code{\link[enviPat:isopattern]{enviPat::isopattern()}}
//' for an alternative function to calculate accurate isotope patterns.
//'
//' Note that abundances are returned relative to the most abundant peak of the isotope
//' pattern. The \code{abundanceLimit}, however, is an absolute abundance.
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
//' @return List with mass and abundace vectors. Abundances are normalized
//' relative to the most abundant peak of the isotope pattern.
//'
//' @family Chemistry functions
//'
//' @examples
//' GetIsotopePattern2("CO2", 1e-5)
//' @export
// [[Rcpp::export]]
List GetIsotopePattern2(std::string molecule, double abundanceLimit) {

  char *cMolecule = StringToChar(molecule);

  int nbrIsotopes = 0;

  // get number of isotopes
  TwRetVal rv = TwGetIsotopePattern2(cMolecule, abundanceLimit, &nbrIsotopes,
                                    NULL, NULL);

  if (rv != TwValueAdjusted) {
    stop(TranslateReturnValue(rv));
  }

  NumericVector isoMass(nbrIsotopes);
  NumericVector isoAbundance(nbrIsotopes);

  rv = TwGetIsotopePattern2(cMolecule, abundanceLimit, &nbrIsotopes, &isoMass[0],
                            &isoAbundance[0]);
  if (rv != TwSuccess) {
    stop(TranslateReturnValue(rv));
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
NumericVector Tof2Mass(NumericVector tofSample, int massCalibMode,
                       NumericVector p) {

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
//' @param tof Vector of sample indices (of same length as \code{mass}).
//' @param weight Vector of weights (if \code{NULL} (default) all weights are set equal).
//' @return Vector of calibration parameters.
//' @export
// [[Rcpp::export]]
NumericVector MassCalibrate(int massCalibMode, NumericVector mass,
                            NumericVector tof,
                            Nullable<NumericVector> weight = R_NilValue) {

  int nbrPoints = mass.size();

  double *p_weight;
  if (weight.isNotNull()) {
    NumericVector x(weight);
    p_weight = &x[0];
  } else {
    p_weight = nullptr;
  }

  int nbrParams;

  // get the number of parameters
  char *description = new char[64];
  TwRetVal rv = TwGetMassCalibInfo(massCalibMode, description, &nbrParams);
  delete[] description;
  if (rv != TwSuccess) {
    stop(TranslateReturnValue(rv));
  }

  NumericVector param(nbrParams);

  rv = TwMassCalibrate(massCalibMode, nbrPoints, &mass[0], &tof[0], p_weight,
                       &nbrParams, &param[0], NULL, NULL);

  if (rv != TwSuccess) {
    stop(TranslateReturnValue(rv));
  }

  return param;
}

// GetMassCalibInfo ------------------------------------------------------------
//' Gets the description and number of parameters of the available mass
//' calibration functions.
//'
//' \code{GetMassCalibInfo} gets the description and number of parameters of the
//' available mass calibration functions.
//'
//' Note: Modes 3 and 4 are flawed. Don't use them. In mode 3 the fit does not
//' converge well, because of a bug (parameters not correctly initialized).
//' Mode 4 is two sequential fits, first mode 0, then a quadratic fit to the
//' residuals, which is an inferior implementation of mode 3. Mode 1 is for FTMS
//' data.
//'
//' @param massCalibMode Mass calibration mode (0 to 5).
//' @return List with the description and number of calibration parameters for
//' the given \code{massCalibMode}.
//' @export
// [[Rcpp::export]]
List GetMassCalibInfo(int massCalibMode) {

  int nbrParams;
  char *description = new char[64];

  TwRetVal rv = TwGetMassCalibInfo(massCalibMode, description, &nbrParams);
  if (rv != TwSuccess) {
    delete[] description;
    stop(TranslateReturnValue(rv));
  }

  std::string str(description);
  delete[] description;

  List result;
  result["description"] = wrap(str);
  result["nbrParams"] = nbrParams;

  return result;
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
//' @export
// [[Rcpp::export]]
void SiInitializeHistograms(NumericVector loMass, NumericVector hiMass,
                            Nullable<IntegerVector> specType = R_NilValue) {

  int nbrHist = loMass.size();
  int *p_specType;
  if (specType.isNotNull()) {
    IntegerVector x(specType);
    p_specType = &x[0];
  } else {
    p_specType = nullptr;
  }

  TwRetVal rv = TwSiInitializeHistograms(nbrHist, &loMass[0], &hiMass[0], p_specType);

  if (rv != TwSuccess) {
    stop(TranslateReturnValue(rv));
  }
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
//' @export
// [[Rcpp::export]]
void SiSetProcessingOptions(std::string option, double value, int specType) {

  char *coption = StringToChar(option);

  TwRetVal rv = TwSiSetProcessingOptions(coption, value, specType);

  if (rv != TwSuccess) {
    stop(TranslateReturnValue(rv));
  }
}

// SiProcessSpectrum -----------------------------------------------------------
//' Processes a spectrum.
//'
//' \code{SiProcessSpectrum} processes a spectrum according to the options set for
//' it's spectrum type.
//'
//' @param spectrum Vector holding the spectrum to process.
//' @param specType Spectrum type index (non-negative integer).
//' @return A list with the baseline and threshold value.
//'
//' @export
// [[Rcpp::export]]
List SiProcessSpectrum(NumericVector spectrum, int specType) {

  int nbrSamples = spectrum.size();

  std::vector<float> fspectrum = as<std::vector<float> >(spectrum);

  float blFromData;
  float thrFromData;

  TwRetVal rv = TwSiProcessSpectrum(&fspectrum[0], nbrSamples, specType,
                                    &blFromData, &thrFromData);

  if (rv != TwSuccess) {
    stop(TranslateReturnValue(rv));
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
//' @param histogramIndex Index (zero-based numbering) of the histogram. It
//' corresponds to the mass range defined with \code{\link{SiInitializeHistograms}}.
//' @return A list with the intensities (histogram x-values), counts (histogram
//' y-values), the number of spectra that were processed for this histogram and
//' the mean histogram value i.e. sum(intensity[i]*counts[i])/sum(counts[i]).
//'
//' @export
// [[Rcpp::export]]
List SiGetHistogram(int histogramIndex) {

  unsigned int arrayLength;
  unsigned int spectrumCount;
  double meanValue;

  //get length of array
  TwRetVal rv = TwSiGetHistogram(histogramIndex, NULL, NULL, &arrayLength,
                                 &spectrumCount, &meanValue);
  if (rv != TwValueAdjusted) {
    stop(TranslateReturnValue(rv));
  }

  NumericVector intensity(arrayLength);
  NumericVector counts(arrayLength);

  std::vector<float> fintensity = as<std::vector<float> >(intensity);
  std::vector<unsigned int> fcounts = as<std::vector<unsigned int> >(counts);

  rv = TwSiGetHistogram(histogramIndex, &fintensity[0], &fcounts[0],
                        &arrayLength, &spectrumCount, &meanValue);
  if (rv != TwSuccess) {
    stop(TranslateReturnValue(rv));
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
//' @export
// [[Rcpp::export]]
List SiGetSumHistogram(int specType, double minMass, double maxMass,
                       double minRate, double maxRate) {

  unsigned int arrayLength;
  unsigned int spectrumCount;
  double meanValue;

  //get length of array
  TwRetVal rv = TwSiGetSumHistogram(specType, NULL, NULL, &arrayLength,
                                    &spectrumCount, &meanValue, minMass, maxMass,
                                    minRate, maxRate);
  if (rv != TwValueAdjusted) {
    stop(TranslateReturnValue(rv));
  }

  NumericVector intensity(arrayLength);
  NumericVector counts(arrayLength);

  std::vector<float> fintensity = as<std::vector<float> >(intensity);
  std::vector<unsigned int> fcounts = as<std::vector<unsigned int> >(counts);

  rv = TwSiGetSumHistogram(specType, &fintensity[0], &fcounts[0], &arrayLength,
                           &spectrumCount, &meanValue, minMass, maxMass,
                           minRate, maxRate);
  if (rv != TwSuccess) {
    stop(TranslateReturnValue(rv));
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
//' @export
// [[Rcpp::export]]
void SiResetHistograms() {

  TwRetVal rv = TwSiResetHistograms();

  if (rv != TwSuccess) {
    stop(TranslateReturnValue(rv));
  }
}

// SiCleanup -------------------------------------------------------------------
//' Cleans up the state in the DLL.
//'
//' \code{SiCleanup} cleans up the state in the DLL. After a call to this
//' function \code{\link{SiInitializeHistograms}} can be called again.
//'
//' @export
// [[Rcpp::export]]
void SiCleanup() {

  TwRetVal rv = TwSiCleanup();

  if (rv != TwSuccess) {
    stop(TranslateReturnValue(rv));
  }
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
//' @export
// [[Rcpp::export]]
List SiFitPhd(NumericVector intensity, NumericVector counts) {

  int nbrPoints = intensity.size();

  NumericVector fwhm(nbrPoints);
  NumericVector a(nbrPoints);
  NumericVector par(nbrPoints);

  TwRetVal rv = TwSiFitPhd(nbrPoints, &intensity[0], &counts[0], &fwhm[0],
                           &a[0], &par[0]);
  if (rv != TwSuccess) {
    stop(TranslateReturnValue(rv));
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
//' @export
// [[Rcpp::export]]
List SiFitRateFromPhd(NumericVector intensity, NumericVector counts,
                      NumericVector siPar) {

  int nbrPoints = intensity.size();

  double rate;
  NumericVector fitCounts(nbrPoints);

  TwRetVal rv = TwSiFitRateFromPhd(nbrPoints, &intensity[0], &counts[0],
                                   &siPar[0], &rate, &fitCounts[0], 0, NULL);
  if (rv != TwSuccess) {
    stop(TranslateReturnValue(rv));
  }

  List result;
  result["rate"] = rate;
  result["fitCounts"] = fitCounts;

  return result;
}

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
//' \dontrun{
//' FindTpsIp("910.33.0316", 500)
//' }
//' @export
// [[Rcpp::export]]
String FindTpsIp(std::string TpsSerial, int timeout) {
#ifdef _WIN32
  char *cTpsSerial = StringToChar(TpsSerial);

  int hostStrLen = 15;

  char *buffer = new char[15];

  TwRetVal rv = TwFindTpsIp(cTpsSerial, timeout, &hostStrLen, buffer);

  if (rv != TwSuccess) {
    delete[] buffer;
    stop(TranslateReturnValue(rv));
  }

  std::string str(buffer);
  delete[] buffer;

  return wrap(str);
#else
  stop("This function is only implemented on Windows.");
#endif
}

// DecomposeMass ---------------------------------------------------------------
//' Calculates possible sum formulas from a list of fragment masses.
//'
//' \code{DecomposeMass} calculates possible sum formulas from a list of fragment
//' masses.
//'
//' Calculates possible sum formulas that amount to target mass +/- tolerance
//' given a target mass, mass tolerance and a list of fragment masses. Filters
//' specifying min/max values for absolute counts of a fragment or for ratios
//' between fragments can be specified in order to reduce the amount of results
//' and restrict hits to chemically reasonable sum formulae. Typically atomMass
//' and atomLabels are the masses/labels of elements but you are free to use
//' whatever you like (isotopes, amino acids, common fragments etc.).
//'
//' @param targetMass Target mass.
//' @param tolerance Tolerance of target mass.
//' @param atomMass Numeric vector with fragment masses.
//' @param atomLabel String Vector of fragment labels.
//' @param elementIndex1 Element index for count filters, first element index for ratio filters.
//' @param elementIndex2 -1 for count filters, second element for ratio filters.
//' @param filterMinVal Counts or ratios that are smaller than this value are filtered out.
//' @param filterMaxVal Counts or ratios that are larger than this value are filtered out.
//'
//' @return List with sum formula strings and the mass and mass errors of the sum
//' formulas.
//'
//' @family Chemistry functions
//'
//' @examples
//' targetMass = 314
//' tolerance = 0.5
//' atomLabel = c("C", "H", "O")
//' n = length(atomLabel)
//' atomMass = rep(0, n)
//'for (i in 1:n) {
//'  atomMass[i] = GetMoleculeMass(atomLabel[i])
//' }
//' elementIndex1 = seq(along.with = atomLabel)-1
//' elementIndex2 = rep(-1, n)
//' filterMinVal = c(20, 20, 0)
//' filterMaxVal = c(22, 40, 5)
//' DecomposeMass(targetMass, tolerance, atomMass, atomLabel, elementIndex1,
//'               elementIndex2, filterMinVal, filterMaxVal)
//' @export
// [[Rcpp::export]]
List DecomposeMass(double targetMass, double tolerance, NumericVector atomMass,
                   StringVector atomLabel, IntegerVector elementIndex1,
                   IntegerVector elementIndex2, NumericVector filterMinVal,
                   NumericVector filterMaxVal) {

  int nbrAtoms = atomMass.size();

  // get size of atomLabel
  int sizeLabel = 0;
  for( int i=0; i < nbrAtoms; i++ ) {
    for(int j=0; j < atomLabel[i].size(); j++){
      sizeLabel += 1;
    }
    sizeLabel += 1;  // null termination
  }

  char *catomLabel = new char[sizeLabel];

  int pos = 0;
  for( int i=0; i < nbrAtoms; i++ ) {
    std::string str(atomLabel[i]);
    strncpy(&catomLabel[pos], str.c_str(), atomLabel[i].size());
    pos += atomLabel[i].size() + 1;
    catomLabel[pos-1] = '\0';
  }

  int nbrCompomers;
  int nbrFilters = elementIndex1.size();

  TwRetVal rv = TwDecomposeMass(targetMass, tolerance, nbrAtoms, &atomMass[0],
                       catomLabel, nbrFilters, &elementIndex1[0], &elementIndex2[0],
                       &filterMinVal[0], &filterMaxVal[0], &nbrCompomers);
  if (rv != TwSuccess) {
    delete[] catomLabel;
    stop(TranslateReturnValue(rv));
  }
  delete[] catomLabel;

  // get composition
  int sumFormulaLength = 256; // assuming all formulas are <256 characters long
  double mass;
  double massError;
  NumericVector massVector(nbrCompomers);
  NumericVector massErrorVector(nbrCompomers);
  StringVector sumFormulaVector(nbrCompomers);

  for( int i=0; i < nbrCompomers; i++ ) {
    char *sumFormula = new char[sumFormulaLength];
    rv = TwGetComposition(i, sumFormula, &sumFormulaLength, &mass, &massError);
    if (rv != TwSuccess) {
      delete[] sumFormula;
      stop(TranslateReturnValue(rv));
    }
    std::string str(sumFormula);
    delete[] sumFormula;
    massVector[i] = mass;
    massErrorVector[i] = massError;
    sumFormulaVector[i] = str;
  }

  List result;
  result["sumFormula"] = sumFormulaVector;
  result["mass"] = massVector;
  result["massError"] = massErrorVector;
  return result;
}

// DecomposeMass2 --------------------------------------------------------------
//' Calculates possible sum formulas from a list of fragment masses.
//'
//' \code{DecomposeMass2} calculates possible sum formulas from a list of fragment
//' masses.
//'
//' Same as \code{\link{DecomposeMass}} but with additional parameters \code{maxHits}
//' and \code{maxSearch} to better fine-tune when to abort a compomer search.
//' Calculates possible sum formulas that amount to target mass +/- tolerance
//' given a target mass, mass tolerance and a list of fragment masses. Filters
//' specifying min/max values for absolute counts of a fragment or for ratios
//' between fragments can be specified in order to reduce the amount of results
//' and restrict hits to chemically reasonable sum formulae. Typically atomMass
//' and atomLabels are the masses/labels of elements but you are free to use
//' whatever you like (isotopes, amino acids, common fragments etc.).
//'
//' @param targetMass Target mass.
//' @param tolerance Tolerance of target mass.
//' @param atomMass Numeric vector with fragment masses.
//' @param atomLabel String Vector of fragment labels.
//' @param elementIndex1 Element index for count filters, first element index for ratio filters.
//' @param elementIndex2 -1 for count filters, second element for ratio filters.
//' @param filterMinVal Counts or ratios that are smaller than this value are filtered out.
//' @param filterMaxVal Counts or ratios that are larger than this value are filtered out.
//' @param maxHits Maximum number of candidate hits to report.
//' @param maxSearch Maximum number of candidate formula to search.
//'
//' @return List with sum formula strings and the mass and mass errors of the sum
//' formulas.
//'
//' @family Chemistry functions
//'
//' @examples
//' targetMass = 314
//' tolerance = 0.5
//' atomLabel = c("C", "H", "O")
//' n = length(atomLabel)
//' atomMass = rep(0, n)
//'for (i in 1:n) {
//'  atomMass[i] = GetMoleculeMass(atomLabel[i])
//' }
//' elementIndex1 = seq(along.with = atomLabel)-1
//' elementIndex2 = rep(-1, n)
//' filterMinVal = c(20, 20, 0)
//' filterMaxVal = c(22, 40, 5)
//' maxHits = 10
//' maxSearch = 1000
//' DecomposeMass2(targetMass, tolerance, atomMass, atomLabel, elementIndex1,
//'               elementIndex2, filterMinVal, filterMaxVal, maxHits, maxSearch)
//' @export
// [[Rcpp::export]]
List DecomposeMass2(double targetMass, double tolerance, NumericVector atomMass,
                  StringVector atomLabel, IntegerVector elementIndex1,
                  IntegerVector elementIndex2, NumericVector filterMinVal,
                  NumericVector filterMaxVal, int maxHits, int maxSearch) {

  int nbrAtoms = atomMass.size();

  // get size of atomLabel
  int sizeLabel = 0;
  for( int i=0; i < nbrAtoms; i++ ) {
   for(int j=0; j < atomLabel[i].size(); j++){
     sizeLabel += 1;
   }
   sizeLabel += 1;  // null termination
  }

  char *catomLabel = new char[sizeLabel];

  int pos = 0;
  for( int i=0; i < nbrAtoms; i++ ) {
   std::string str(atomLabel[i]);
   strncpy(&catomLabel[pos], str.c_str(), atomLabel[i].size());
   pos += atomLabel[i].size() + 1;
   catomLabel[pos-1] = '\0';
  }

  int nbrCompomers;
  int nbrFilters = elementIndex1.size();

  TwRetVal rv = TwDecomposeMass2(targetMass, tolerance, nbrAtoms, &atomMass[0],
                                catomLabel, nbrFilters, &elementIndex1[0],
                                &elementIndex2[0], &filterMinVal[0],
                                &filterMaxVal[0], &nbrCompomers, maxHits, maxSearch);
  if (rv != TwSuccess) {
   delete[] catomLabel;
   stop(TranslateReturnValue(rv));
  }
  delete[] catomLabel;

  // get composition
  int sumFormulaLength = 256; // assuming all formulas are <256 characters long
  double mass;
  double massError;
  NumericVector massVector(nbrCompomers);
  NumericVector massErrorVector(nbrCompomers);
  StringVector sumFormulaVector(nbrCompomers);

  for( int i=0; i < nbrCompomers; i++ ) {
   char *sumFormula = new char[sumFormulaLength];
   rv = TwGetComposition(i, sumFormula, &sumFormulaLength, &mass, &massError);
   if (rv != TwSuccess) {
     delete[] sumFormula;
     stop(TranslateReturnValue(rv));
   }
   std::string str(sumFormula);
   delete[] sumFormula;
   massVector[i] = mass;
   massErrorVector[i] = massError;
   sumFormulaVector[i] = str;
  }

  List result;
  result["sumFormula"] = sumFormulaVector;
  result["mass"] = massVector;
  result["massError"] = massErrorVector;
  return result;
}

// MatchSpectra ----------------------------------------------------------------
//' Checks two spectra for similarity.
//'
//' \code{MatchSpectra} matches two spectra and returns a similarity score
//' (0 - 100 %).
//'
//' Values < 0.0 in spec1 or spec2 will be set to 0.0. Currently, only
//' matchMethod 0 is implemented and relies on the euclidean distance between
//' the spectra.
//'
//' @param spec1 First spectrum to match.
//' @param spec2 Second spectrum to match.
//' @param matchMethod Method to use. Currently only method 0 is implemented.
//' @return Match score in units of percent.
//'
//' @export
// [[Rcpp::export]]
double MatchSpectra(NumericVector spec1, NumericVector spec2,
                    int matchMethod = 0) {

  int nbrPoints = spec1.size();
  if (nbrPoints != spec2.size()) {
    stop("spec1 and spec2 must be the same length.");
  }
  double matchScore;

  TwRetVal rv = TwMatchSpectra(&spec1[0], &spec2[0], nbrPoints, matchMethod, &matchScore);
  if (rv != TwSuccess) {
    stop(TranslateReturnValue(rv));
  }
  return matchScore;
}

// Not implemented -------------------------------------------------------------
// TwNistLibrarySearch
// TwNistLibraryQueryResult
// TwBruteForceCalibrate
// TwEncImsCorrelateProfile
// TwEncImsCorrelateMultiProfiles
// TwEncImsCleanup
// TwSiGetHistogramAmp
// TwSiGetSumHistogramAmp
// TwGenerateImageData
// TwImagePaletteOffset2Value
// TwImageValue2PaletteOffset
// TwImageGetPaletteRGB
// TwIntegrateTofSpectra
// TwIntegrateTofSpectrum
// TwMakeMqAxis
