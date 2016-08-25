#include <Rcpp.h>
using namespace Rcpp;

//' Decodes an event list.
//'
//' \code{DecodeEventList} decodes an event list read with \code{GetEventList...FromH5}
//' functions into a time stamp vector and a data vector.
//'
//' Note: this function is not part of the TofDaq API, but is included in the
//' package for convenience.
//'
//' @param events Event data, e.g. from \code{GetEventList...FromH5}.
//' @param clockPeriod Clock period in ns, e.g. from
//' \code{GetFloatAttributeFromH5(filename, "FullSpectra", "ClockPeriod")}.
//' @param sampleInterval Sampling interval in ns, e.g. from \code{\link{GetH5Descriptor}}.
//' @return A list with sample indices and data values (in mV).
//' @export
// [[Rcpp::export]]
List DecodeEventList(NumericVector events, double clockPeriod, double sampleInterval) {

  const unsigned int n = events.size();
  List out;
  NumericVector sampleindex(n);
  NumericVector value(n);
  unsigned int k = 0;

  for (unsigned int i = 0; i < n; ++i) {
    const unsigned int timestamp = (unsigned int)events[i] & 0xFFFFFF;
    const unsigned int dataLength = (unsigned int)events[i] >> 24;
    if (dataLength == 0) {  // TDC data -> each event is 1 count
      // convert timestamp to a sample index
      sampleindex[k] = (int)(timestamp * clockPeriod / sampleInterval + 0.5);
      value[k] = 1;
      k++;
    }
    else { // ADC data or "Ndigo TDC" data (timestamp is time of first sample in packet)
      for (unsigned int j = 0; j < dataLength; ++j) {
        ++i;
        const unsigned int mask = events[i];
        float* adcData = (float*) & mask;
        sampleindex[k] = (int)(timestamp * clockPeriod  / sampleInterval + 0.5) + j;
        value[k] = *adcData;
        k++;
      }
    }
  }

  IntegerVector idx(k);
  for (unsigned int i = 0; i < k; i++) {
    idx[i] = i;
  }

  out["sampleindex"] = sampleindex[idx];
  out["value"] = value[idx];
  return out;
}
