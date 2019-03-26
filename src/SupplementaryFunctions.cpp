#include <Rcpp.h>
using namespace Rcpp;
#include <TwToolDll.h>
#ifdef _WIN32
#include <TofDaqDll.h>
#endif
#include <queue>  // std::queue (FIFO)
#include "TofDaqR.h"

// tof -------------------------------------------------------------------------
//' Time-of-flight calculation.
//'
//' \code{tof} calculates the time-of-flight of ions in an ideal two-stage
//' extraction/two-stage reflection TOFMS (with second-order space focusing).
//'
//' \code{"HTOF-W"} is the W-mode and assumes hard mirror voltage = 1.5*pulse voltage.
//' \code{x} and \code{v} can be vectors, but must be the same length.
//'
//' Note: this function is not part of the TofDaq API, but is included in the
//' package for convenience.
//'
//' @param toftype TOF type (\code{"CTOF"}, \code{"HTOF"}, \code{"HTOF-W"} or \code{"LTOF"})
//' @param drift Drift tube voltage (V).
//' @param pulse Extraction pulse voltage (V).
//' @param massToCharge Mass-to-charge ratio (Th).
//' @param x Initial position deviation of the ion(s) from the extraction plane (m).
//' @param v Initial velocity of the ion(s) in extraction direction (m/s).
//'
//' @references Laiko, V.V. and A.F. Dodonov (1994), Resolution and Spectral-line
//' Shapes in the Reflecting Time-of-flight Mass Spectrometer with Orthogonally
//' Injected Ions, \emph{Rapid Communications in Mass Spectrometry}, \strong{8},
//' 720-726.
//'
//' @return Time-of-flight (s)
//'
//' @keywords internal
//' @export
// [[Rcpp::export]]
NumericVector tof(std::string toftype = "LTOF", double drift = 6000,
                  double pulse = 1000, double massToCharge = 100,
                  NumericVector x = NumericVector::create(0),
                  NumericVector v = NumericVector::create(0)) {

  const double amu = 1.660538921e-27;  // atomic mass unit (kg)
  const double e = 1.60217657e-19;  // elementary charge (C)

  char *cToftype = StringToChar(toftype);

  int nbr = x.size();
  if (nbr != v.size()) {
    stop("length(x) != length(v)");
  }

  double Vdrift = -std::abs(drift);
  double Vpush = std::abs(pulse);
  double Vpull = -Vpush;

  double d1 = 0.004 + 0.0035;  // pull to push distance

  double u1 = Vpush - Vpull;
  double u3 = Vpull - Vdrift;
  double d3, d4, d5, d6, d7;
  double u5, u6;
  double u7 = 0;  // initialized to prevent -Wmaybe-uninitialized warning.

  if (strcmp(cToftype, "HTOF") == 0) {
    d3 = 0.0065;
    d4 = 0.507 + 0.511;
    d5 = 0.0165;
    d6 = 0.0713 - d5;
  } else if (strcmp(cToftype, "HTOF-W") == 0) {
    d3 = 0.0065;
    d4 = 0.507 + 0.511 + 2 * 0.5125;
    d5 = 0.0165;
    d6 = 0.0713 - d5;
    d7 = 0.0085;
    u7 = 1.5 * Vpush - Vdrift;
  } else if (strcmp(cToftype, "LTOF") == 0) {
    d3 = 0.014;
    d4 = 1.0465 + 1.052;
    d5 = 0.0385;
    d6 = 0.1593 - d5;
  } else if (strcmp(cToftype, "CTOF") == 0) {
    d3 = 0.006;
    d4 = 2 * 0.1435;
    d5 = 0.017;
    d6 = 0.0165;
  } else if (strcmp(cToftype, "NTOF") == 0) {
    d3 = 0.014;
    d4 = 3*0.8;
    d5 = 0.0385;
    d6 = 0.1593 - d5;
    d7 = d5;
    u7 = 1.5*Vpush - Vdrift;
  } else {
    stop("invalid toftype");
  }

  // calculate u5 and u6
  double x0 = d1-0.001;
  double k0 = x0/d1;
  double p0 = k0+u3/u1;
  double a, b;
  if (strcmp(cToftype, "HTOF-W") == 0) {
    a = d1/u1*pow(k0, -0.5) + d3/u3*(pow(p0, -0.5) - pow(k0, -0.5)) - d4/u1/2*pow(p0, -1.5) + 2*d7/u7*pow(p0, -0.5);
    b = d1/u1/2*pow(k0, -1.5) + d3/u3/2*(pow(p0, -1.5) - pow(k0, -1.5)) - d4/u1*3/4*pow(p0, -2.5) + d7/u7*pow(p0, -1.5);
    u5 = (a-2*p0*b + 4*d5/u1*pow(p0, -1.5))/(-2*b/u1);
    u6 = (-4*d6*pow(p0-u5/u1, -0.5))/(a+4*d5/u5*(pow(p0, -0.5) - pow(p0-u5/u1, -0.5)));
  } else if (strcmp(cToftype, "NTOF") == 0) {
    a = d1/u1*pow(k0, -0.5) + d3/u3*(pow(p0, -0.5) - pow(k0, -0.5)) - d4/u1/2*pow(p0, -1.5) + 2*d7/u7*pow(p0, -0.5);
    b = d1/u1/2*pow(k0, -1.5) + d3/u3/2*(pow(p0, -1.5) - pow(k0, -1.5)) - d4/u1*3/4*pow(p0, -2.5) + d7/u7*pow(p0, -1.5);
    u5 = (a-2*p0*b + 2*d5/u1*pow(p0, -1.5))/(-2*b/u1);
    u6 = (-2*d6*pow(p0-u5/u1, -0.5))/(a+2*d5/u5*(pow(p0, -0.5) - pow(p0-u5/u1, -0.5)));
  } else {
    a = d1/u1*pow(k0, -0.5) + d3/u3*(pow(p0, -0.5) - pow(k0, -0.5)) - d4/u1/2*pow(p0, -1.5);
    b = d1/u1/2*pow(k0, -1.5) + d3/u3/2*(pow(p0, -1.5) - pow(k0, -1.5)) - d4/u1*3/4*pow(p0, -2.5);
    u5 = (a-2*p0*b + 2*d5/u1*pow(p0, -1.5))/(-2*b/u1);
    u6 = (-2*d6*pow(p0-u5/u1, -0.5))/(a+2*d5/u5*(pow(p0, -0.5) - pow(p0-u5/u1, -0.5)));
  }

  NumericVector xi(nbr);
  NumericVector k(nbr);

  xi = v * sqrt(massToCharge*amu / (2 * e*u1));
  k = (x0 - x) / d1;

  NumericVector timeOfFlight(nbr);

  // calculate time-of-flight
  if (strcmp(cToftype, "HTOF-W") == 0) {

    timeOfFlight = sqrt(massToCharge*amu/(2*e)) *
      (2*d1/sqrt(u1)*(sqrt(xi*xi + k) - xi) +
      2*d3/u3*(sqrt(u1*xi*xi + k*u1 + u3) - sqrt(u1*xi*xi + k*u1)) +
      d4/sqrt(u1*xi*xi + k*u1 + u3) +
      2*4*d5/u5*(sqrt(u1*xi*xi + k*u1 + u3) - sqrt(u1*xi*xi + k*u1 + u3 - u5)) +
      2*4*d6/u6*sqrt(u1*xi*xi + k*u1 + u3 - u5) +
      4*d7/u7*sqrt(u1*xi*xi + k*u1 + u3));
  } else if (strcmp(cToftype, "NTOF") == 0) {
    timeOfFlight = sqrt(massToCharge*amu/(2*e)) *
      (2*d1/sqrt(u1)*(sqrt(xi*xi + k) - xi) +
      2*d3/u3*(sqrt(u1*xi*xi + k*u1 + u3) - sqrt(u1*xi*xi + k*u1)) +
      d4/sqrt(u1*xi*xi + k*u1 + u3) +
      4*d5/u5*(sqrt(u1*xi*xi + k*u1 + u3) - sqrt(u1*xi*xi + k*u1 + u3 - u5)) +
      4*d6/u6*sqrt(u1*xi*xi + k*u1 + u3 - u5) +
      4*d7/u7*sqrt(u1*xi*xi + k*u1 + u3));
  } else {
    timeOfFlight = sqrt(massToCharge*amu/(2*e)) *
      (2*d1/sqrt(u1)*(sqrt(xi*xi + k) - xi) +
      2*d3/u3*(sqrt(u1*xi*xi + k*u1 + u3) - sqrt(u1*xi*xi + k*u1)) +
      d4/sqrt(u1*xi*xi + k*u1 + u3) +
      4*d5/u5*(sqrt(u1*xi*xi + k*u1 + u3) - sqrt(u1*xi*xi + k*u1 + u3 - u5)) +
      4*d6/u6*sqrt(u1*xi*xi + k*u1 + u3 - u5));
  }

  return timeOfFlight;
}

// EventList2TofSpec -----------------------------------------------------------
//' Converts events of an event list into a spectrum.
//'
//' \code{EventList2TofSpec} converts events of an event list (read with
//' \code{GetEventList...FromH5} functions) into a spectrum.
//'
//' Note: this function is not part of the TofDaq API, but is included in the
//' package for convenience.
//'
//' @param events Event data, e.g. from \code{GetEventList...FromH5}.
//' @param clockPeriod Clock period (in s or ns), e.g. from
//' \code{GetFloatAttributeFromH5(filename, "FullSpectra", "ClockPeriod")}.
//' @param sampleInterval Sampling interval (in same units as \code{clockPeriod}),
//'  e.g. from \code{\link{GetH5Descriptor}}.
//' @param nbrSamples Number of samples, e.g. from \code{\link{GetH5Descriptor}}.
//' @return A vector containing the spectrum.
//' @export
// [[Rcpp::export]]
NumericVector EventList2TofSpec(NumericVector events, double clockPeriod,
                       double sampleInterval, int nbrSamples) {

  unsigned int n = events.size();

  NumericVector spectrum(nbrSamples);

  for (unsigned int i = 0; i < n; ++i) {
    const unsigned int timestamp = (unsigned int)events[i] & 0xFFFFFF;
    //convert timestamp to a sample index
    const unsigned int sampleIndex = (unsigned int)(timestamp / (sampleInterval / clockPeriod) + 0.5);
    const unsigned int dataLength = (unsigned int)events[i] >> 24;
    if (dataLength == 0) { //TDC data -> each event is 1 count
      spectrum[sampleIndex] += 1.0;
    }
    else { //ADC data (timestamp is time of first sample in packet)
      for (unsigned int j = 0; j < dataLength; ++j) {
        ++i;
        unsigned int mask = (unsigned int)events[i];
        float* adcData = (float*)& mask;
        spectrum[sampleIndex + j] += *adcData;
      }
    }
  }

  return spectrum;
}

// DecodeEventList -------------------------------------------------------------
//' Decodes an event list.
//'
//' \code{DecodeEventList} decodes an event list read with \code{GetEventList...FromH5}
//' functions into a time stamp vector and a data vector.
//'
//' Note: this function is not part of the TofDaq API, but is included in the
//' package for convenience.
//'
//' @param events Event data, e.g. from \code{GetEventList...FromH5}.
//' @param clockPeriod Clock period (in s or ns), e.g. from
//' \code{GetFloatAttributeFromH5(filename, "FullSpectra", "ClockPeriod")}.
//' @param sampleInterval Sampling interval (in same units as \code{clockPeriod}),
//'  e.g. from \code{\link{GetH5Descriptor}}.
//' @return A list with sample indices and data values (in mV).
//' @export
// [[Rcpp::export]]
List DecodeEventList(NumericVector events, double clockPeriod, double sampleInterval) {

  int n = events.size();

  IntegerVector sampleindex(n);
  NumericVector value(n);

  int k = 0;
  int i = 0;
  while (i < n) {
    unsigned int timestamp = (unsigned int)events[i] & 0xFFFFFF;
    unsigned int dataLength = (unsigned int)events[i] >> 24;
    if (dataLength == 0) {  // TDC data -> each event is 1 count
      // convert timestamp to a sample index
      sampleindex[k] = (int)(timestamp * clockPeriod / sampleInterval + 0.5);
      value[k] = 1;
      k++;
    } else { // ADC data or "Ndigo TDC" data (timestamp is time of first sample in packet)
      for (unsigned int j = 0; j < dataLength; ++j) {
        i++;
        if (i < n) {  // required if sample size of event > 256
          unsigned int mask = (unsigned int)events[i];
          float* adcData = (float*)& mask;
          sampleindex[k] = (int)(timestamp * clockPeriod / sampleInterval + 0.5) + j;
          value[k] = *adcData;
          k++;
        }
      }
    }
    i++;
  }

  IntegerVector idx(k);
  for (int i = 0; i < k; i++) {
    idx[i] = i;
  }

  List result;
  result["sampleIndex"] = sampleindex[idx];
  result["value"] = value[idx];

  return result;
}

// DecodeEventListThreshold ----------------------------------------------------
//' Decodes an event list using thresholding.
//'
//' \code{DecodeEventListThreshold} decodes an event list read with \code{GetEventList...FromH5}
//' functions into a time stamp vector and a data vector. Only event data which is
//' above the threshold (plus pre-trigger and post-trigger samples) are returned.
//'
//' Note: this function is not part of the TofDaq API, but is included in the
//' package for convenience.
//'
//' @param events Event data, e.g. from \code{GetEventList...FromH5}.
//' @param clockPeriod Clock period (in s or ns), e.g. from
//' \code{GetFloatAttributeFromH5(filename, "FullSpectra", "ClockPeriod")}.
//' @param sampleInterval Sampling interval (in same units as \code{clockPeriod}),
//'  e.g. from \code{\link{GetH5Descriptor}}.
//' @param threshold Threshold value (mV).
//' @param presamples Number of pre-trigger samples.
//' @param postsamples Number of post-trigger samples.
//' @return A list with sample indices and data values (in mV).
//'
//' @keywords internal
//' @export
// [[Rcpp::export]]
List DecodeEventListThreshold(NumericVector events, double clockPeriod,
                              double sampleInterval, double threshold,
                              int presamples, int postsamples) {

  const int n = events.size();
  IntegerVector sampleindex(n);
  NumericVector value(n);
  int k = 0;
  bool peakdetected = false;
  int m = 0;  // post samples iterator
  std::queue<int> pre_index;  // pre samples queue
  std::queue<double> pre_value;  // pre samples queue

  for (int i = 0; i < n; ++i) {
    const unsigned int timestamp = (unsigned int)events[i] & 0xFFFFFF;
    const unsigned int dataLength = (unsigned int)events[i] >> 24;
    if (dataLength == 0) {  // TDC data -> each event is 1 count
      // convert timestamp to a sample index
      sampleindex[k] = (int)(timestamp * clockPeriod / sampleInterval + 0.5);
      value[k] = 1;
      k++;
    }
    else {  // ADC data or "Ndigo TDC" data (timestamp is time of first sample in packet)
      pre_index = std::queue<int>();  // empty pre samples queue
      pre_value = std::queue<double>();  // empty pre samples queue
      m = 0;
      peakdetected = false;

      for (int j = 0; j < (int)dataLength; ++j) {
        ++i;
        const unsigned int mask = events[i];
        float* adcData = (float*) & mask;

        // buffer presamples data
        if (presamples > 0 and *adcData < threshold and !peakdetected) {
          pre_index.push( (int)(timestamp * clockPeriod  / sampleInterval + 0.5) + j);  // insert element at the end
          pre_value.push( *adcData );  // insert element at the end
          if (pre_index.size() > (unsigned int)presamples) {
            pre_index.pop();  // remove first element;
            pre_value.pop();  // remove first element;
          }
        }

        // copy presamples at threshold crosser
        if (*adcData >= threshold and !peakdetected) {
          peakdetected = true;
          if (presamples > 0) {
            // copy buffered data
            unsigned int presize = pre_index.size();
            for (unsigned int l = 0; l < presize; ++l) {
              sampleindex[k] = pre_index.front();  // first element
              pre_index.pop();  // remove first element
              value[k] = pre_value.front();  // first element
              pre_value.pop();  // remove first element
              k++;
            }
            pre_index = std::queue<int>();  // empty queue
            pre_value = std::queue<double>();  // empty queue
          }
        }

        // above threshold
        if (*adcData >= threshold and peakdetected) {
          sampleindex[k] = (int)(timestamp * clockPeriod  / sampleInterval + 0.5) + j;
          value[k] = *adcData;
          k++;
        }

        // post samples
        if (*adcData < threshold and peakdetected) {
          if (m < postsamples) {
            sampleindex[k] = (int)(timestamp * clockPeriod  / sampleInterval + 0.5) + j;
            value[k] = *adcData;
            k++;
            m++;
          } else {
            peakdetected = false;
            m = 0;
          }
        }
      }
    }
  }

  IntegerVector idx(k);
  for (int i = 0; i < k; i++) {
    idx[i] = i;
  }

  List result;
  result["sampleIndex"] = sampleindex[idx];
  result["value"] = value[idx];
  return result;
}

// SiProcessSpectrumFromShMem --------------------------------------------------
//' Processes a spectrum taken from shared memory.
//'
//' \code{SiProcessSpectrumFromShMem} processes a spectrum taken from shared
//' memory according to the options set for it's spectrum type.
//'
//' This function is a variant of the original TwToolDll function \code{\link{SiProcessSpectrum}}.
//'
//' @param specType Spectrum type index (non-negative integer).
//' @param BufIndex Buf index of data to fetch.
//' @return A list with the baseline and threshold value.
//'
//' @export
// [[Rcpp::export]]
List SiProcessSpectrumFromShMem(int specType, int BufIndex) {
#ifdef _WIN32
  //get descriptor of file
  TSharedMemoryDesc desc;
  TwRetVal rv = TwGetDescriptor(&desc);
  if (rv != TwSuccess) {
    stop(TranslateReturnValue(rv));
  }

  int specLen = desc.NbrSamples;

  std::vector<float> spectrum(specLen);

  float blFromData;
  float thrFromData;

  rv = TwSiProcessSpectrum(&spectrum[0], specLen, specType, &blFromData,
                           &thrFromData);
  if (rv != TwSuccess) {
    stop(TranslateReturnValue(rv));
  }

  List result;
  result["baseline"] = wrap(blFromData);
  result["threshold"] = wrap(thrFromData);

  return result;
#else
  stop("This function is only implemented on Windows.");
#endif
}

// KeepSharedMemMapped ---------------------------------------------------------
//' Keeps the shared memory acquisition buffers mapped.
//'
//' \code{KeepSharedMemMapped} Keeps the shared memory acquisition buffers mapped.
//'
//' The DLL periodically unmaps the shared memory to give the recorder
//' application the possibility to (re)allocate the shared buffers. Call this
//' function if you want to make sure that the shared memory pointers stay
//' valid while you work with them. In this case you must call
//' \code{\link{ReleaseSharedMemory}} explicitly when finished with your
//' processing operation.
//'
//' @export
// [[Rcpp::export]]
void KeepSharedMemMapped() {
#ifdef _WIN32
  // Note: This is part of TwGetSharedMemory, but here extracted as an
  // independent function.

  TSharedMemoryPointer pShMem;
  bool keepMapped = true;

  TwRetVal rv = TwGetSharedMemory(&pShMem, keepMapped);

  if (rv != TwSuccess) {
    stop(TranslateReturnValue(rv));
  }
#else
  stop("This function is only implemented on Windows.");
#endif
}
