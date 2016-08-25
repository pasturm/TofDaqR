#include <Rcpp.h>
using namespace Rcpp;

//' Converts events of an event list into a spectrum.
//'
//' \code{EventList2TofSpec} converts events of an event list (read with
//' \code{GetEventList...FromH5} functions) into a spectrum.
//'
//' Note: this function is not part of the TofDaq API, but is included in the
//' package for convenience.
//'
//' @param events Event data, e.g. from \code{GetEventList...FromH5}.
//' @param clockPeriod Clock period in ns, e.g. from
//' \code{GetFloatAttributeFromH5(filename, "FullSpectra", "ClockPeriod")}.
//' @param sampleInterval Sampling interval in ns, e.g. from \code{\link{GetH5Descriptor}}.
//' @param nbrSamples Number of samples, e.g. from \code{\link{GetH5Descriptor}}.
//' @return A vector containing the spectrum.
//' @export
// [[Rcpp::export]]
NumericVector EventList2TofSpec(NumericVector events, double clockPeriod, double sampleInterval, double nbrSamples) {

  unsigned int n = events.size();
  NumericVector spectrum(nbrSamples);

  for (unsigned int i = 0; i < n; ++i) {
    const unsigned int timestamp = (unsigned int)events[i] & 0xFFFFFF;
    //convert timestamp to a sample index
  	const unsigned int sampleIndex = (unsigned int)(timestamp/(sampleInterval/clockPeriod) + 0.5);
		const unsigned int dataLength = (unsigned int)events[i] >> 24;
		if (dataLength == 0) { //TDC data -> each event is 1 count
			spectrum[sampleIndex] += 1.0f;
		}
		else { //ADC data (timestamp is time of first sample in packet)
			for (unsigned int j = 0; j < dataLength; ++j) {
        ++i;
        unsigned int mask = events[i];
  		  float* adcData = (float*) & mask;
				spectrum[sampleIndex + j] += *adcData;
			}
		}
  }

  return spectrum;
}
