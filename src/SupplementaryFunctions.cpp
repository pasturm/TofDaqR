#include <Rcpp.h>
using namespace Rcpp;
#include <TwToolDll.h>
#ifdef _WIN32
#include <TofDaqDll.h>
#endif
#include <queue>  // std::queue (FIFO)
#include "TofDaqR.h"

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
