#include <Rcpp.h>
using namespace Rcpp;
#include <TwH5Dll.h>
#include "TofDaqR.h"

// GetH5Descriptor -------------------------------------------------------------
//' Descriptor structure of Tofwerk HDF5 data file.
//'
//' \code{GetH5Descriptor} returns a descriptor structure for the Tofwerk HDF5 file.
//'
//' The \emph{TwH5Desc} structure contains information about data dimensions,
//' available datasets and mass calibration. Additional attributes, which are not
//' available in the structure can be read using \code{Get...AttributeFromH5}
//' functions.
//' See
//' \href{http://htmlpreview.github.io/?https://github.com/pasturm/TofDaqR/blob/master/tools/doc/TwH5Dll.htm}{TofDaq API documentation}
//' for more details.
//'
//' @param Filename Path/filename of the HDF5 file.
//' @return A list containing the \emph{TwH5Desc} structure.
//'
//' @examples
//' \dontrun{
//' GetH5Descriptor("path/to/file.h5")
//' }
//' @export
// [[Rcpp::export]]
List GetH5Descriptor(std::string Filename) {

  char *cFilename = StringToChar(Filename);

  TwH5Desc descriptor;
  TwRetVal rv = TwGetH5Descriptor(cFilename, &descriptor);

  TwCloseH5(cFilename);

  if (rv != TwSuccess) {
    stop(TranslateReturnValue(rv));
  }

  List result;
  result["nbrSamples"] = descriptor.nbrSamples;
  result["nbrPeaks"] = descriptor.nbrPeaks;
  result["nbrWaveforms"] = descriptor.nbrWaveforms;
  result["nbrSegments"] = descriptor.nbrSegments;
  result["nbrBlocks"] = descriptor.nbrBlocks;
  result["nbrMemories"] = descriptor.nbrMemories;
  result["nbrBufs"] = descriptor.nbrBufs;
  result["nbrWrites"] = descriptor.nbrWrites;
  result["nbrLogEntries"] = descriptor.nbrLogEntries;
  result["secondTof"] = descriptor.secondTof;
  result["hasSumSpectrum"] = descriptor.hasSumSpectrum;
  result["hasSumSpectrum2"] = descriptor.hasSumSpectrum2;
  result["hasBufTimes"] = descriptor.hasBufTimes;
  result["hasTofData"] = descriptor.hasTofData;
  result["hasTofData2"] = descriptor.hasTofData2;
  result["hasPeakData"] = descriptor.hasPeakData;
  result["hasPeakData2"] = descriptor.hasPeakData2;
  result["hasTpsData"] = descriptor.hasTpsData;
  result["hasNbrMemories"] = descriptor.hasNbrMemories;
  result["hasPressureData"] = descriptor.hasPressureData;
  result["hasLogData"] = descriptor.hasLogData;
  result["hasMassCalibData"] = descriptor.hasMassCalibData;
  result["hasMassCalib2Data"] = descriptor.hasMassCalib2Data;
  result["hasCh1RawData"] = descriptor.hasCh1RawData;
  result["hasCh2RawData"] = descriptor.hasCh2RawData;
  result["hasRawDataDesc"] = descriptor.hasRawDataDesc;
  result["hasEventList"] = descriptor.hasEventList;
  result["segIlf"] = descriptor.segIlf;
  result["eventListMaxElementLength"] = descriptor.eventListMaxElementLength;
  result["daqMode"] = descriptor.daqMode;
  result["acquisitionMode"] = descriptor.acquisitionMode;
  result["massCalibMode"] = descriptor.massCalibMode;
  result["massCalibMode2"] = descriptor.massCalibMode2;
  result["nbrCalibParams"] = descriptor.nbrCalibParams;
  result["nbrCalibParams2"] = descriptor.nbrCalibParams2;
  result["nbrCubes"] = descriptor.nbrCubes;
  NumericVector p(descriptor.p, descriptor.p + 16);
  result["p"] = p;
  NumericVector p2(descriptor.p2, descriptor.p2 + 16);
  result["p2"] = p2;
  result["tofPeriod"] = descriptor.tofPeriod;
  result["blockPeriod"] = descriptor.blockPeriod;
  result["sampleInterval"] = descriptor.sampleInterval;
  result["singleIonSignal"] = descriptor.singleIonSignal;
  result["singleIonSignal2"] = descriptor.singleIonSignal2;

  return result;
}

// CloseH5 ---------------------------------------------------------------------
//' Closes an open HDF5 file.
//'
//' \code{CloseH5} closes an open HDF5 file.
//'
//' This function is called internally by all \code{Get..FromH5} functions, so
//' it is usually not necessary to call \code{CloseH5} explicitely.
//'
//' @param Filename Path/filename of the HDF5 file.
//'
//' @examples
//' \dontrun{
//' CloseH5("path/to/file.h5")
//' }
//' @export
// [[Rcpp::export]]
void CloseH5(std::string Filename) {

  char *cFilename = StringToChar(Filename);

  TwRetVal rv = TwCloseH5(cFilename);

  if (rv != TwSuccess) {
    stop(TranslateReturnValue(rv));
  }
}

// CloseAll --------------------------------------------------------------------
//' Closes all open HDF5 files.
//'
//' \code{CloseAll} closes all open HDF5 files. It is a good idea to call this
//' function once before your program exits.
//'
//' @examples
//' \dontrun{
//' CloseAll()
//' }
//' @export
// [[Rcpp::export]]
void CloseAll() {

  TwRetVal rv = TwCloseAll();

  if (rv != TwSuccess) {
    stop(TranslateReturnValue(rv));
  }
}

// GetSumSpectrumFromH5 --------------------------------------------------------
//' Sum spectrum from HDF5 data file.
//'
//' \code{GetSumSpectrumFromH5} reads the sum spectrum (or average spectrum
//' depending on \code{Normalize} flag) from the HDF5 file.
//'
//' @param Filename Path/filename of the HDF5 file.
//' @param Normalize If \code{FALSE} (default) the spectrum is reported as sum,
//' if \code{TRUE} the spectrum is normalized to counts per extraction.
//' @return A vector containing the sum spectrum.
//'
//' @examples
//' \dontrun{
//' GetSumSpectrumFromH5("path/to/file.h5")
//' }
//' @export
// [[Rcpp::export]]
NumericVector GetSumSpectrumFromH5(std::string Filename, bool Normalize = false) {

  char *cFilename = StringToChar(Filename);

  TwH5Desc descriptor;
  TwRetVal rv = TwGetH5Descriptor(cFilename, &descriptor);
  if (rv != TwSuccess) {
    TwCloseH5(cFilename);
    stop(TranslateReturnValue(rv));
  }

  NumericVector Spectrum(descriptor.nbrSamples);

  rv = TwGetSumSpectrumFromH5(cFilename, &Spectrum[0], Normalize);

  TwCloseH5(cFilename);

  if (rv != TwSuccess) {
    stop(TranslateReturnValue(rv));
  }

  return Spectrum;
}

// GetTofSpectrumFromH5 --------------------------------------------------------
//' Single (averaged) TOF spectrum from HDF5 data file.
//'
//' \code{GetTofSpectrumFromH5} reads a single mass spectrum (or an
//' averaged/summed hyperslab) from the HDF5 file.
//'
//' If \code{SegmentIndex == SegmentEndIndex} and \code{BufIndex == BufEndIndex} and
//' \code{WriteIndex == WriteEndIndex} and \code{Normalize == FALSE} no
//' averaging/summing of spectra is done and the spectrum is reported as stored
//' in the dataset.
//'
//' @param Filename Path/filename of the HDF5 file.
//' @param SegmentIndex Segment start index of data to fetch.
//' @param SegmentEndIndex Segment end index of data to fetch.
//' @param BufIndex Buf start index of data to fetch.
//' @param BufEndIndex Buf end index of data to fetch.
//' @param WriteIndex Write start index of data to fetch.
//' @param WriteEndIndex Write end index of data to fetch.
//' @param BufWriteLinked Indicating whether the buf and write dimension should
//' be considered linked or treated as independent dimensions (relevant only if
//' \code{WriteIndex != WriteEndIndex}). Default is \code{FALSE}.
//' @param Normalize If \code{FALSE} the spectrum is reported as sum,
//' if \code{TRUE} (default) the spectrum is normalized to counts per extraction.
//' @return A vector containing the mass spectrum.
//'
//' @examples
//' \dontrun{
//' GetTofSpectrumFromH5("path/to/file.h5")
//' }
//' @export
// [[Rcpp::export]]
SEXP GetTofSpectrumFromH5(std::string Filename, int SegmentIndex,
                          int SegmentEndIndex, int BufIndex, int BufEndIndex,
                          int WriteIndex, int WriteEndIndex,
                          bool BufWriteLinked = false, bool Normalize = true) {

  char *cFilename = StringToChar(Filename);

  TwH5Desc descriptor;
  TwRetVal rv = TwGetH5Descriptor(cFilename, &descriptor);
  if (rv != TwSuccess) {
    TwCloseH5(cFilename);
    stop(TranslateReturnValue(rv));
  }

  std::vector<float> Spectrum(descriptor.nbrSamples);

  rv = TwGetTofSpectrumFromH5(cFilename, &Spectrum[0], SegmentIndex,
                              SegmentEndIndex, BufIndex, BufEndIndex,
                              WriteIndex, WriteEndIndex, BufWriteLinked,
                              Normalize);
  TwCloseH5(cFilename);

  if (rv != TwSuccess) {
    stop(TranslateReturnValue(rv));
  }

  return wrap(Spectrum);
}

// GetTofSpectrum2FromH5 -------------------------------------------------------
//' Single (averaged) TOF spectrum from HDF5 data file.
//'
//' \code{GetTofSpectrum2FromH5} reads a single mass spectrum (or an
//' averaged/summed hyperslab) from the HDF5 file.
//'
//' If \code{SegmentIndex == SegmentEndIndex} and \code{BufIndex == BufEndIndex} and
//' \code{WriteIndex == WriteEndIndex} and \code{Normalize == FALSE} no
//' averaging/summing of spectra is done and the spectrum is reported as stored
//' in the dataset.
//'
//' @param Filename Path/filename of the HDF5 file.
//' @param SegmentIndex Segment start index of data to fetch.
//' @param SegmentEndIndex Segment end index of data to fetch.
//' @param BufIndex Buf start index of data to fetch.
//' @param BufEndIndex Buf end index of data to fetch.
//' @param WriteIndex Write start index of data to fetch.
//' @param WriteEndIndex Write end index of data to fetch.
//' @param BufWriteLinked Indicating whether the buf and write dimension should
//' be considered linked or treated as independent dimensions (relevant only if
//' \code{WriteIndex != WriteEndIndex}). Default is \code{FALSE}.
//' @param Normalize If \code{FALSE} the spectrum is reported as sum,
//' if \code{TRUE} (default) the spectrum is normalized to counts per extraction.
//' @return A vector containing the mass spectrum.
//'
//' @examples
//' \dontrun{
//' GetTofSpectrum2FromH5("path/to/file.h5")
//' }
//' @export
// [[Rcpp::export]]
SEXP GetTofSpectrum2FromH5(std::string Filename, int SegmentIndex,
                           int SegmentEndIndex, int BufIndex, int BufEndIndex,
                           int WriteIndex, int WriteEndIndex,
                           bool BufWriteLinked = false, bool Normalize = true) {

  char *cFilename = StringToChar(Filename);

  TwH5Desc descriptor;
  TwRetVal rv = TwGetH5Descriptor(cFilename, &descriptor);
  if (rv != TwSuccess) {
    TwCloseH5(cFilename);
    stop(TranslateReturnValue(rv));
  }

  std::vector<float> Spectrum(descriptor.nbrSamples);

  rv = TwGetTofSpectrum2FromH5(cFilename, &Spectrum[0], SegmentIndex,
                               SegmentEndIndex, BufIndex, BufEndIndex,
                               WriteIndex, WriteEndIndex, BufWriteLinked,
                               Normalize);
  TwCloseH5(cFilename);

  if (rv != TwSuccess) {
    stop(TranslateReturnValue(rv));
  }

  return wrap(Spectrum);
}

// GetStickSpectrumFromH5 ------------------------------------------------------
//' Single (averaged) stick spectrum from HDF5 data file.
//'
//' \code{GetStickSpectrumFromH5} reads a single stick spectrum (or an
//' averaged/summed hyperslab) from the HDF5 file.
//'
//' If \code{SegmentIndex == SegmentEndIndex} and \code{BufIndex == BufEndIndex} and
//' \code{WriteIndex == WriteEndIndex} and \code{Normalize == FALSE} no
//' averaging/summing of spectra is done and the spectrum is reported as stored
//' in the dataset.
//'
//' @param Filename Path/filename of the HDF5 file.
//' @param SegmentIndex Segment start index of data to fetch.
//' @param SegmentEndIndex Segment end index of data to fetch.
//' @param BufIndex Buf start index of data to fetch.
//' @param BufEndIndex Buf end index of data to fetch.
//' @param WriteIndex Write start index of data to fetch.
//' @param WriteEndIndex Write end index of data to fetch.
//' @param BufWriteLinked Indicating whether the buf and write dimension should
//' be considered linked or treated as independent dimensions (relevant only if
//' \code{WriteIndex != WriteEndIndex}). Default is \code{FALSE}.
//' @param Normalize If \code{FALSE} the spectrum is reported as sum,
//' if \code{TRUE} (default) the spectrum is normalized to counts per extraction.
//' @return A vector containing the stick spectrum.
//'
//' @examples
//' \dontrun{
//' GetStickSpectrumFromH5("path/to/file.h5")
//' }
//' @export
// [[Rcpp::export]]
SEXP GetStickSpectrumFromH5(std::string Filename, int SegmentIndex,
                            int SegmentEndIndex, int BufIndex, int BufEndIndex,
                            int WriteIndex, int WriteEndIndex,
                            bool BufWriteLinked = false, bool Normalize = true) {

  char *cFilename = StringToChar(Filename);

  TwH5Desc descriptor;
  TwRetVal rv = TwGetH5Descriptor(cFilename, &descriptor);
  if (rv != TwSuccess) {
    TwCloseH5(cFilename);
    stop(TranslateReturnValue(rv));
  }

  std::vector<float> Spectrum(descriptor.nbrPeaks);

  rv = TwGetStickSpectrumFromH5(cFilename, &Spectrum[0], SegmentIndex,
                                SegmentEndIndex, BufIndex, BufEndIndex,
                                WriteIndex, WriteEndIndex, BufWriteLinked,
                                Normalize);
  TwCloseH5(cFilename);

  if (rv != TwSuccess) {
    stop(TranslateReturnValue(rv));
  }

  return wrap(Spectrum);
}

// GetStickSpectrum2FromH5 -----------------------------------------------------
//' Single (averaged) stick spectrum from HDF5 data file.
//'
//' \code{GetStickSpectrum2FromH5} reads a single stick spectrum (or an
//' averaged/summed hyperslab) from the HDF5 file.
//'
//' If \code{SegmentIndex == SegmentEndIndex} and \code{BufIndex == BufEndIndex} and
//' \code{WriteIndex == WriteEndIndex} and \code{Normalize == FALSE} no
//' averaging/summing of spectra is done and the spectrum is reported as stored
//' in the dataset.
//'
//' @param Filename Path/filename of the HDF5 file.
//' @param SegmentIndex Segment start index of data to fetch.
//' @param SegmentEndIndex Segment end index of data to fetch.
//' @param BufIndex Buf start index of data to fetch.
//' @param BufEndIndex Buf end index of data to fetch.
//' @param WriteIndex Write start index of data to fetch.
//' @param WriteEndIndex Write end index of data to fetch.
//' @param BufWriteLinked Indicating whether the buf and write dimension should
//' be considered linked or treated as independent dimensions (relevant only if
//' \code{WriteIndex != WriteEndIndex}). Default is \code{FALSE}.
//' @param Normalize If \code{FALSE} the spectrum is reported as sum,
//' if \code{TRUE} (default) the spectrum is normalized to counts per extraction.
//' @return A vector containing the stick spectrum.
//'
//' @examples
//' \dontrun{
//' GetStickSpectrum2FromH5("path/to/file.h5")
//' }
//' @export
// [[Rcpp::export]]
SEXP GetStickSpectrum2FromH5(std::string Filename, int SegmentIndex,
                             int SegmentEndIndex, int BufIndex, int BufEndIndex,
                             int WriteIndex, int WriteEndIndex,
                             bool BufWriteLinked = false, bool Normalize = true) {

  char *cFilename = StringToChar(Filename);

  TwH5Desc descriptor;
  TwRetVal rv = TwGetH5Descriptor(cFilename, &descriptor);
  if (rv != TwSuccess) {
    TwCloseH5(cFilename);
    stop(TranslateReturnValue(rv));
  }

  std::vector<float> Spectrum(descriptor.nbrPeaks);

  rv = TwGetStickSpectrum2FromH5(cFilename, &Spectrum[0], SegmentIndex,
                                SegmentEndIndex, BufIndex, BufEndIndex,
                                WriteIndex, WriteEndIndex, BufWriteLinked,
                                Normalize);
  TwCloseH5(cFilename);

  if (rv != TwSuccess) {
    stop(TranslateReturnValue(rv));
  }

  return wrap(Spectrum);
}

// GetPeakParametersFromH5 -----------------------------------------------------
//' Peak parameters from HDF5 data file.
//'
//' \code{GetPeakParametersFromH5} reads peak parameters from the data file.
//'
//' @param Filename Path/filename of the HDF5 file.
//' @param PeakIndex Index of peak (zero-based numbering). If index is -1
//' (default), peak parameters of all peaks are read.
//' @return A list with the peak paramters \emph{label}, \emph{mass}, \emph{loMass} and \emph{hiMass}.
//'
//' @examples
//' \dontrun{
//' GetPeakParametersFromH5("path/to/file.h5")
//' }
//' @export
// [[Rcpp::export]]
List GetPeakParametersFromH5(std::string Filename, int PeakIndex = -1) {

  char *cFilename = StringToChar(Filename);

  if (PeakIndex == -1) {

    TwH5Desc descriptor;
    TwRetVal rv = TwGetH5Descriptor(cFilename, &descriptor);
    if (rv != TwSuccess) {
      TwCloseH5(cFilename);
      stop(TranslateReturnValue(rv));
    }
    TPeakPar *PeakPar = new TPeakPar[descriptor.nbrPeaks];

    rv = TwGetPeakParametersFromH5(cFilename, PeakPar, PeakIndex);
    TwCloseH5(cFilename);
    if (rv != TwSuccess) {
      delete[] PeakPar;
      stop(TranslateReturnValue(rv));
    }

    List result;
    CharacterVector label(descriptor.nbrPeaks);
    NumericVector mass(descriptor.nbrPeaks);
    NumericVector loMass(descriptor.nbrPeaks);
    NumericVector hiMass(descriptor.nbrPeaks);

    for (int i=0; i<descriptor.nbrPeaks; ++i) {
      label[i] = PeakPar[i].label;
      mass[i] = PeakPar[i].mass;
      loMass[i] = PeakPar[i].loMass;
      hiMass[i] = PeakPar[i].hiMass;
    }

    delete[] PeakPar;

    result["label"] = label;
    result["mass"] = mass;
    result["loMass"] = loMass;
    result["hiMass"] = hiMass;

    return result;

  } else {

    TPeakPar PeakPar;
    TwRetVal rv = TwGetPeakParametersFromH5(cFilename, &PeakPar, PeakIndex);

    TwCloseH5(cFilename);

    if (rv != TwSuccess) {
      stop(TranslateReturnValue(rv));
    }

    List result;
    result["label"] = PeakPar.label;
    result["mass"] = PeakPar.mass;
    result["loMass"] = PeakPar.loMass;
    result["hiMass"] = PeakPar.hiMass;

    return result;
  }
}

// GetBufTimeFromH5 ------------------------------------------------------------
//' Single buf timestamp from the data file.
//'
//' \code{GetBufTimeFromH5} reads a single buf time stamp from HDF5 file.
//'
//' @param Filename Path/filename of the HDF5 file.
//' @param BufIndex Buf index.
//' @param WriteIndex Write index.
//' @return A time stamp (in seconds relative to acquisition start).
//'
//' @examples
//' \dontrun{
//' GetBufTimeFromH5("path/to/file.h5", BufIndex = 0, WriteIndex = 0)
//' }
//' @export
// [[Rcpp::export]]
double GetBufTimeFromH5(std::string Filename, int BufIndex, int WriteIndex) {

  char *cFilename = StringToChar(Filename);

  double BufTime;

  TwRetVal rv = TwGetBufTimeFromH5(cFilename, &BufTime, BufIndex, WriteIndex);

  TwCloseH5(cFilename);

  if (rv != TwSuccess) {
    stop(TranslateReturnValue(rv));
  }

  return BufTime;
}

// GetSpecXaxisFromH5 ----------------------------------------------------------
//' x-axis values of mass spectrum.
//'
//' \code{GetSpecXaxisFromH5} returns an array of x-axis values of the mass
//' spectrum.
//'
//' @param Filename Path/filename of the HDF5 file.
//' @param Type x-axis type (0: sample index, 1: mass/charge [Th] (default),
//' -1: mass/charge [Th] (2nd TOF), 2: time of flight [microsec],
//' -2: time of flight [microsec] (2nd TOF), 3: frequency [kHz]).
//' @param writeIndex Write index to use for mass calibration (relevant only for
//' \code{abs(Type)== 1 or 2}). If the data file has no \emph{/TofData/MassCalibration}
//' dataset the standard mass calibration parameters are used (same for all
//' values of writeIndex). Default is 0.
//' @return A vector containing the x-axis values.
//'
//' @examples
//' \dontrun{
//' GetSpecXaxisFromH5("path/to/file.h5")
//' }
//' @export
// [[Rcpp::export]]
NumericVector GetSpecXaxisFromH5(std::string Filename, int Type = 1,
                                 int writeIndex = 0) {

  char *cFilename = StringToChar(Filename);

  TwH5Desc descriptor;
  TwRetVal rv = TwGetH5Descriptor(cFilename, &descriptor);
  if (rv != TwSuccess) {
    TwCloseH5(cFilename);
    stop(TranslateReturnValue(rv));
  }

  NumericVector SpecAxis(descriptor.nbrSamples);

  rv = TwGetSpecXaxisFromH5(cFilename, &SpecAxis[0], Type, NULL, 0.0, writeIndex);

  TwCloseH5(cFilename);

  if (rv != TwSuccess) {
    stop(TranslateReturnValue(rv));
  }

  return SpecAxis;
}

// GetSegmentProfileFromH5 -----------------------------------------------------
//' Segment profile from HDF5 data file.
//'
//' \code{GetSegmentProfileFromH5} reads a segment profile for a given peak (or
//' all peaks) and averaged over a given buf and write range.
//'
//' @param Filename Path/filename of the HDF5 file.
//' @param PeakIndex Index of peak to fetch segment profile from. All peaks are
//' read if \code{PeakIndex = -1}.
//' @param BufStartIndex Buf start index of data to fetch.
//' @param BufEndIndex Buf end index of data to fetch.
//' @param WriteStartIndex Write start index of data to fetch.
//' @param WriteEndIndex Write end index of data to fetch.
//' @param BufWriteLinked Indicating whether the buf and write dimension should
//' be considered linked or treated as independent dimensions (relevant only if
//' \code{WriteStartIndex != WriteEndIndex}). Default is \code{FALSE}.
//' @return A vector containing the segment profile(s).
//'
//' @examples
//' \dontrun{
//' GetSegmentProfileFromH5("path/to/file.h5", PeakIndex = -1, BufStartIndex = 0,
//' BufEndIndex = 0, WriteStartIndex = 0, WriteEndIndex = 0)
//' }
//' @export
// [[Rcpp::export]]
SEXP GetSegmentProfileFromH5(std::string Filename, int PeakIndex,
                             int BufStartIndex, int BufEndIndex,
                             int WriteStartIndex, int WriteEndIndex,
                             bool BufWriteLinked = false) {

  char *cFilename = StringToChar(Filename);

  TwH5Desc descriptor;
  TwRetVal rv = TwGetH5Descriptor(cFilename, &descriptor);
  if (rv != TwSuccess) {
    TwCloseH5(cFilename);
    stop(TranslateReturnValue(rv));
  }

  if (PeakIndex != -1) {

    std::vector<float> SegmentProfile(descriptor.nbrSegments);

    rv = TwGetSegmentProfileFromH5(cFilename, &SegmentProfile[0], PeakIndex,
                                   BufStartIndex, BufEndIndex, WriteStartIndex,
                                   WriteEndIndex, BufWriteLinked);
    TwCloseH5(cFilename);
    if (rv != TwSuccess) {
      stop(TranslateReturnValue(rv));
    }

    return wrap(SegmentProfile);

  } else {

    std::vector<float> SegmentProfile(descriptor.nbrSegments*descriptor.nbrPeaks);

    rv = TwGetSegmentProfileFromH5(cFilename, &SegmentProfile[0], PeakIndex,
                                   BufStartIndex, BufEndIndex, WriteStartIndex,
                                   WriteEndIndex, BufWriteLinked);
    TwCloseH5(cFilename);
    if (rv != TwSuccess) {
      stop(TranslateReturnValue(rv));
    }

    return wrap(SegmentProfile);
  }
}

// GetSegmentProfile2FromH5 ----------------------------------------------------
//' Segment profile from HDF5 data file.
//'
//' \code{GetSegmentProfile2FromH5} reads a segment profile for a given peak (or
//' all peaks) and averaged over a given buf and write range.
//'
//' @param Filename Path/filename of the HDF5 file.
//' @param PeakIndex Index of peak to fetch segment profile from. All peaks are
//' read if \code{PeakIndex = -1}.
//' @param BufStartIndex Buf start index of data to fetch.
//' @param BufEndIndex Buf end index of data to fetch.
//' @param WriteStartIndex Write start index of data to fetch.
//' @param WriteEndIndex Write end index of data to fetch.
//' @param BufWriteLinked Indicating whether the buf and write dimension should
//' be considered linked or treated as independent dimensions (relevant only if
//' \code{WriteStartIndex != WriteEndIndex}). Default is \code{FALSE}.
//' @return A vector containing the segment profile(s).
//'
//' @examples
//' \dontrun{
//' GetSegmentProfile2FromH5("path/to/file.h5", PeakIndex = -1, BufStartIndex = 0,
//' BufEndIndex = 0, WriteStartIndex = 0, WriteEndIndex = 0)
//' }
//' @export
// [[Rcpp::export]]
SEXP GetSegmentProfile2FromH5(std::string Filename, int PeakIndex,
                              int BufStartIndex, int BufEndIndex,
                              int WriteStartIndex, int WriteEndIndex,
                              bool BufWriteLinked = false) {

  char *cFilename = StringToChar(Filename);

  TwH5Desc descriptor;
  TwRetVal rv = TwGetH5Descriptor(cFilename, &descriptor);
  if (rv != TwSuccess) {
    TwCloseH5(cFilename);
    stop(TranslateReturnValue(rv));
  }

  if (PeakIndex != -1) {

    std::vector<float> SegmentProfile(descriptor.nbrSegments);

    rv = TwGetSegmentProfile2FromH5(cFilename, &SegmentProfile[0], PeakIndex,
                                    BufStartIndex, BufEndIndex, WriteStartIndex,
                                    WriteEndIndex, BufWriteLinked);
    TwCloseH5(cFilename);
    if (rv != TwSuccess) {
      stop(TranslateReturnValue(rv));
    }

    return wrap(SegmentProfile);

  } else {

    std::vector<float> SegmentProfile(descriptor.nbrSegments*descriptor.nbrPeaks);

    rv = TwGetSegmentProfile2FromH5(cFilename, &SegmentProfile[0], PeakIndex,
                                    BufStartIndex, BufEndIndex, WriteStartIndex,
                                    WriteEndIndex, BufWriteLinked);
    TwCloseH5(cFilename);
    if (rv != TwSuccess) {
      stop(TranslateReturnValue(rv));
    }

    return wrap(SegmentProfile);
  }
}

// GetBufWriteProfileFromH5 ----------------------------------------------------
//' Gets a linked buf/write profile.
//'
//' \code{GetBufWriteProfileFromH5} gets a linked buf/write profile for a given
//' peak (or all peaks) and segment slice. If your data is not linked, use
//' \code{\link{GetPeakData}} to get the buf and/or write profiles.
//'
//' @param Filename Path/filename of the HDF5 file.
//' @param PeakIndex Index of peak to fetch buf/write profile from. All peaks are
//' read if \code{PeakIndex = -1}.
//' @param SegmentStartIndex Segment start index of data to fetch.
//' @param SegmentEndIndex Segment end index of data to fetch.
//' @return A vector containing the buf/write profile(s).
//'
//' @examples
//' \dontrun{
//' GetBufWriteProfileFromH5("path/to/file.h5", PeakIndex = -1,
//' SegmentStartIndex = 0, SegmentEndIndex = 0)
//' }
//' @export
// [[Rcpp::export]]
SEXP GetBufWriteProfileFromH5(std::string Filename, int PeakIndex,
                              int SegmentStartIndex, int SegmentEndIndex) {

  char *cFilename = StringToChar(Filename);

  TwH5Desc descriptor;
  TwRetVal rv = TwGetH5Descriptor(cFilename, &descriptor);
  if (rv != TwSuccess) {
    TwCloseH5(cFilename);
    stop(TranslateReturnValue(rv));
  }

  if (PeakIndex != -1) {

    std::vector<float> Profile(descriptor.nbrBufs*descriptor.nbrWrites);

    rv = TwGetBufWriteProfileFromH5(cFilename, &Profile[0], PeakIndex,
                                    SegmentStartIndex, SegmentEndIndex);
    TwCloseH5(cFilename);
    if (rv != TwSuccess) {
      stop(TranslateReturnValue(rv));
    }

    return wrap(Profile);

  } else {

    std::vector<float> Profile(descriptor.nbrBufs*descriptor.nbrWrites*
      descriptor.nbrPeaks);

    rv = TwGetBufWriteProfileFromH5(cFilename, &Profile[0], PeakIndex,
                                    SegmentStartIndex, SegmentEndIndex);
    TwCloseH5(cFilename);
    if (rv != TwSuccess) {
      stop(TranslateReturnValue(rv));
    }

    return wrap(Profile);
  }
}

// GetBufWriteProfile2FromH5 ----------------------------------------------------
//' Gets a linked buf/write profile.
//'
//' \code{GetBufWriteProfile2FromH5} gets a linked buf/write profile for a given
//' peak (or all peaks) and segment slice. If your data is not linked, use
//' \code{\link{GetPeakData2}} to get the buf and/or write profiles.
//'
//' @param Filename Path/filename of the HDF5 file.
//' @param PeakIndex Index of peak to fetch buf/write profile from. All peaks are
//' read if \code{PeakIndex = -1}.
//' @param SegmentStartIndex Segment start index of data to fetch.
//' @param SegmentEndIndex Segment end index of data to fetch.
//' @return A vector containing the buf/write profile(s).
//'
//' @examples
//' \dontrun{
//' GetBufWriteProfile2FromH5("path/to/file.h5", PeakIndex = -1,
//' SegmentStartIndex = 0, SegmentEndIndex = 0)
//' }
//' @export
// [[Rcpp::export]]
SEXP GetBufWriteProfile2FromH5(std::string Filename, int PeakIndex,
                               int SegmentStartIndex, int SegmentEndIndex) {

  char *cFilename = StringToChar(Filename);

  TwH5Desc descriptor;
  TwRetVal rv = TwGetH5Descriptor(cFilename, &descriptor);
  if (rv != TwSuccess) {
    TwCloseH5(cFilename);
    stop(TranslateReturnValue(rv));
  }

  if (PeakIndex != -1) {

    std::vector<float> Profile(descriptor.nbrBufs*descriptor.nbrWrites);

    rv = TwGetBufWriteProfile2FromH5(cFilename, &Profile[0], PeakIndex,
                                    SegmentStartIndex, SegmentEndIndex);
    TwCloseH5(cFilename);
    if (rv != TwSuccess) {
      stop(TranslateReturnValue(rv));
    }

    return wrap(Profile);

  } else {

    std::vector<float> Profile(descriptor.nbrBufs*descriptor.nbrWrites*
      descriptor.nbrPeaks);

    rv = TwGetBufWriteProfile2FromH5(cFilename, &Profile[0], PeakIndex,
                                    SegmentStartIndex, SegmentEndIndex);
    TwCloseH5(cFilename);
    if (rv != TwSuccess) {
      stop(TranslateReturnValue(rv));
    }

    return wrap(Profile);
  }
}

// GetRegUserDataSourcesFromH5 -------------------------------------------------
//' Lists all registered user datasets available in the data file.
//'
//' \code{GetRegUserDataSourcesFromH5} lists all registered user data sets
//' available in the data file. Registered data sources can originate from data
//' source plugins, TofDaq recorder (e.g. DAQ temperatures) or data registered
//' through \code{RegisterUserData...} functions.
//'
//' @param Filename Path/filename of the HDF5 file.
//' @return A list containing the location, the length, whether is has a description
//' and the type of the data source dataset. type 1: data source values are
//' written to disk for every write, type 2: data source values are written to
//' disk for every buf.
//'
//' @examples
//' \dontrun{
//' GetRegUserDataSourcesFromH5("path/to/file.h5")
//' }
//' @export
// [[Rcpp::export]]
List GetRegUserDataSourcesFromH5(std::string Filename) {

  char *cFilename = StringToChar(Filename);

  // get nbrSources
  int nbrSources = 0;
  TwRetVal rv = TwGetRegUserDataSourcesFromH5(cFilename, &nbrSources, NULL,
                                              NULL, NULL, NULL);
  if (rv != TwValueAdjusted) {
    TwCloseH5(cFilename);
    stop(TranslateReturnValue(rv));
  }

  char *sourceLocation = new char[256 * nbrSources];

  IntegerVector sourceLength(nbrSources);
  // LogicalVector hasDesc(nbrSources);  // does not work
  // std::vector<bool> hasDesc(nbrSources);  // does not work either
  bool *hasDesc = new bool[nbrSources];
  IntegerVector type(nbrSources);

  rv = TwGetRegUserDataSourcesFromH5(cFilename, &nbrSources, sourceLocation,
                                     &sourceLength[0], hasDesc, &type[0]);

  TwCloseH5(cFilename);
  if (rv != TwSuccess) {
    delete[] sourceLocation;
    delete[] hasDesc;
    stop(TranslateReturnValue(rv));
  }

  CharacterVector locationArray(nbrSources);
  std::string str(sourceLocation, 256 * nbrSources);
  LogicalVector hasDescArray(nbrSources);

  for (int i = 0; i < nbrSources; ++i) {
    locationArray[i] = str.substr(i*256, 256);
    hasDescArray[i] = hasDesc[i];
  }
  delete[] sourceLocation;
  delete[] hasDesc;

  List result;
  result["sourceLocation"] = locationArray;
  result["sourceLength"] = sourceLength;
  result["hasDesc"] = hasDescArray;
  result["type"] = type;

  return result;
}

// GetRegUserDataFromH5 --------------------------------------------------------
//' Reads entries from a registered data source dataset.
//'
//' \code{GetRegUserDataFromH5} reads an entry (or all entries) from a
//' registered data source dataset (created by \code{RegisterUserData...} functions).
//'
//' If \code{bufIndex = writeIndex = -1}, all data are read.
//'
//' @param Filename Path/filename of the HDF5 file.
//' @param location Location of the group or dataset where the attribute is attached to.
//' @param bufIndex Buf index.
//' @param writeIndex Write index.
//' @param readDescription If \code{TRUE} (default) the data descripton is read,
//' if \code{false} the data description is not read.
//' @return A list containing the registered user data and description, or a
//' numeric vector containing the data only if \code{readDescription = FALSE}.
//'
//' @examples
//' \dontrun{
//' GetRegUserDataFromH5("path/to/file.h5", location = )
//' }
//' @export
// [[Rcpp::export]]
SEXP GetRegUserDataFromH5(std::string Filename, std::string location,
                          int bufIndex, int writeIndex,
                          bool readDescription = true) {

  char *cFilename = StringToChar(Filename);
  char *cLocation = StringToChar(location);

  // get bufLength
  int bufLength = 0;
  TwRetVal rv = TwGetRegUserDataFromH5(cFilename, cLocation, bufIndex,
                                       writeIndex, &bufLength, NULL, NULL);
  if (rv != TwValueAdjusted) {
    TwCloseH5(cFilename);
    stop(TranslateReturnValue(rv));
  }

  NumericVector buffer(bufLength);

  if (readDescription) {

    int descLength = 0;
    if (bufIndex == -1 && writeIndex == -1) {
      // get descLength
      TwRetVal rv = TwGetRegUserDataFromH5(cFilename, cLocation, 0,
                                           0, &descLength, NULL, NULL);
      if (rv != TwValueAdjusted) {
        TwCloseH5(cFilename);
        stop(TranslateReturnValue(rv));
      }
    } else {
      descLength = bufLength;
    }


    char *description = new char[256 * descLength];

    rv = TwGetRegUserDataFromH5(cFilename, cLocation, bufIndex, writeIndex,
                                &bufLength, &buffer[0], description);
    TwCloseH5(cFilename);
    if (rv != TwSuccess) {
      delete[] description;
      stop(TranslateReturnValue(rv));
    }

    CharacterVector descriptionArray(descLength);
    std::string str(description, 256 * descLength);
    delete[] description;

    for (int i = 0; i < descLength; ++i) {
      descriptionArray[i] = str.substr(i*256, 256);
    }

    List result;
    result["data"] = buffer;
    result["description"] = descriptionArray;

    return result;
  } else {

    rv = TwGetRegUserDataFromH5(cFilename, cLocation, bufIndex, writeIndex,
                                &bufLength, &buffer[0], NULL);
    TwCloseH5(cFilename);
    if (rv != TwSuccess) {
      stop(TranslateReturnValue(rv));
    }

    return buffer;
  }
}

// GetTofData ------------------------------------------------------------------
//' Gets data stored in /FullSpectra/TofData from HDF5 data file.
//'
//' \code{GetTofData} gets data stored in \code{/FullSpectra/TofData} from HDF5 data file.
//'
//' @param Filename Path/filename of the HDF5 file.
//' @param sampleOffset Sample offset.
//' @param sampleCount Sample count.
//' @param segOffset Segment offset.
//' @param segCount Segment count.
//' @param bufOffset Buf offset.
//' @param bufCount Buf count.
//' @param writeOffset Write offset.
//' @param writeCount Write count.
//' @return A vector of length sampleCount*segCount*bufCount*writeCount.
//' @export
// [[Rcpp::export]]
SEXP GetTofData(std::string Filename, int sampleOffset, int sampleCount,
                int segOffset, int segCount, int bufOffset, int bufCount,
                int writeOffset, int writeCount) {

  char *cFilename = StringToChar(Filename);

  int n = sampleCount*segCount*bufCount*writeCount;

  std::vector<float> dataBuffer(n);

  TwRetVal rv = TwGetTofData(cFilename, sampleOffset, sampleCount, segOffset,
                             segCount, bufOffset, bufCount, writeOffset,
                             writeCount, &dataBuffer[0]);
  TwCloseH5(cFilename);

  if (rv != TwSuccess) {
    stop(TranslateReturnValue(rv));
  }

  return wrap(dataBuffer);
}


// GetTofData2 -----------------------------------------------------------------
//' Gets data stored in /FullSpectra2/TofData from HDF5 data file.
//'
//' \code{GetTofData2} gets data stored in \code{/FullSpectra2/TofData} from HDF5 data file.
//'
//' @param Filename Path/filename of the HDF5 file.
//' @param sampleOffset Sample offset.
//' @param sampleCount Sample count.
//' @param segOffset Segment offset.
//' @param segCount Segment count.
//' @param bufOffset Buf offset.
//' @param bufCount Buf count.
//' @param writeOffset Write offset.
//' @param writeCount Write count.
//' @return A vector of length sampleCount*segCount*bufCount*writeCount.
//' @export
// [[Rcpp::export]]
SEXP GetTofData2(std::string Filename, int sampleOffset, int sampleCount,
                 int segOffset, int segCount, int bufOffset, int bufCount,
                 int writeOffset, int writeCount) {

  char *cFilename = StringToChar(Filename);

  int n = sampleCount*segCount*bufCount*writeCount;

  std::vector<float> dataBuffer(n);

  TwRetVal rv = TwGetTofData2(cFilename, sampleOffset, sampleCount, segOffset,
                             segCount, bufOffset, bufCount, writeOffset,
                             writeCount, &dataBuffer[0]);
  TwCloseH5(cFilename);

  if (rv != TwSuccess) {
    stop(TranslateReturnValue(rv));
  }

  return wrap(dataBuffer);
}

// GetPeakData -----------------------------------------------------------------
//' Gets data stored in /PeakData/PeakData from HDF5 data file.
//'
//' \code{GetPeakData} gets data stored in \code{/PeakData/PeakData} from HDF5 data file.
//'
//' @param Filename Path/filename of the HDF5 file.
//' @param peakOffset Peak offset.
//' @param peakCount Peak count.
//' @param segOffset Segment offset.
//' @param segCount Segment count.
//' @param bufOffset Buf offset.
//' @param bufCount Buf count.
//' @param writeOffset Write offset.
//' @param writeCount Write count.
//' @return A vector of length peakCount*segCount*bufCount*writeCount.
//' @export
// [[Rcpp::export]]
SEXP GetPeakData(std::string Filename, int peakOffset, int peakCount,
                 int segOffset, int segCount, int bufOffset, int bufCount,
                 int writeOffset, int writeCount) {

  char *cFilename = StringToChar(Filename);

  int n = peakCount*segCount*bufCount*writeCount;

  std::vector<float> dataBuffer(n);

  TwRetVal rv = TwGetPeakData(cFilename, peakOffset, peakCount, segOffset,
                              segCount, bufOffset, bufCount, writeOffset,
                              writeCount, &dataBuffer[0]);
  TwCloseH5(cFilename);

  if (rv != TwSuccess) {
    stop(TranslateReturnValue(rv));
  }

  return wrap(dataBuffer);
}

// GetPeakData2 ----------------------------------------------------------------
//' Gets data stored in /PeakData2/PeakData from HDF5 data file.
//'
//' \code{GetPeakData2} gets data stored in \code{/PeakData2/PeakData} from HDF5 data file.
//'
//' @param Filename Path/filename of the HDF5 file.
//' @param peakOffset Peak offset.
//' @param peakCount Peak count.
//' @param segOffset Segment offset.
//' @param segCount Segment count.
//' @param bufOffset Buf offset.
//' @param bufCount Buf count.
//' @param writeOffset Write offset.
//' @param writeCount Write count.
//' @return A vector of length peakCount*segCount*bufCount*writeCount.
//' @export
// [[Rcpp::export]]
SEXP GetPeakData2(std::string Filename, int peakOffset, int peakCount,
                  int segOffset, int segCount, int bufOffset, int bufCount,
                  int writeOffset, int writeCount) {

  char *cFilename = StringToChar(Filename);

  int n = peakCount*segCount*bufCount*writeCount;

  std::vector<float> dataBuffer(n);

  TwRetVal rv = TwGetPeakData2(cFilename, peakOffset, peakCount, segOffset,
                              segCount, bufOffset, bufCount, writeOffset,
                              writeCount, &dataBuffer[0]);
  TwCloseH5(cFilename);

  if (rv != TwSuccess) {
    stop(TranslateReturnValue(rv));
  }

  return wrap(dataBuffer);
}

// GetTimingData ---------------------------------------------------------------
//' Gets data stored in /Timing/BufTimes from HDF5 data file.
//'
//' \code{GetTimingData} gets data stored in \code{/Timing/BufTimes} from HDF5 data file.
//'
//' @param Filename Path/filename of the HDF5 file.
//' @param bufOffset Buf offset.
//' @param bufCount Buf count.
//' @param writeOffset Write offset.
//' @param writeCount Write count.
//' @return A vector of length bufCount*writeCount.
//' @export
// [[Rcpp::export]]
NumericVector GetTimingData(std::string Filename, int bufOffset, int bufCount,
                            int writeOffset, int writeCount) {

  char *cFilename = StringToChar(Filename);

  int n = bufCount*writeCount;

  NumericVector dataBuffer(n);

  TwRetVal rv = TwGetTimingData(cFilename, bufOffset, bufCount, writeOffset,
                              writeCount, &dataBuffer[0]);
  TwCloseH5(cFilename);

  if (rv != TwSuccess) {
    stop(TranslateReturnValue(rv));
  }

  return dataBuffer;
}

// GetIntAttributeFromH5 -------------------------------------------------------
//' Reads an integer attribute from the HDF5 file.
//'
//' \code{GetIntAttributeFromH5} reads an integer attribute from the HDF5 file.
//'
//' Used to read attributes not available from \code{GetH5Descriptor}.
//'
//' @param Filename Path/filename of the HDF5 file.
//' @param location Location of the group or dataset where the attribute is attached to.
//' @param name Attribute name.
//' @return An integer attribute.
//'
//' @examples
//' \dontrun{
//' GetIntAttributeFromH5("path/to/file.h5", location = , name = )
//' }
//' @export
// [[Rcpp::export]]
int GetIntAttributeFromH5(std::string Filename, std::string location,
                          std::string name) {

  char *cFilename = StringToChar(Filename);
  char *cLocation = StringToChar(location);
  char *cName = StringToChar(name);

  int value;
  TwRetVal rv = TwGetIntAttributeFromH5(cFilename, cLocation, cName, &value);

  TwCloseH5(cFilename);
  if (rv != TwSuccess) {
    stop(TranslateReturnValue(rv));
  }

  return value;
}

// GetUintAttributeFromH5 ------------------------------------------------------
//' Reads an unsigned integer attribute from the HDF5 file.
//'
//' \code{GetUintAttributeFromH5} reads an unsigned integer attribute from the
//' HDF5 file.
//'
//' Used to read attributes not available from \code{GetH5Descriptor}. Unsigned
//' integers are returned as numeric values in R.
//'
//' @param Filename Path/filename of the HDF5 file.
//' @param location Location of the group or dataset where the attribute is attached to.
//' @param name Attribute name.
//' @return An numeric attribute.
//'
//' @examples
//' \dontrun{
//' GetUintAttributeFromH5("path/to/file.h5", location = , name = )
//' }
//' @export
// [[Rcpp::export]]
unsigned int GetUintAttributeFromH5(std::string Filename, std::string location,
                                    std::string name) {

  char *cFilename = StringToChar(Filename);
  char *cLocation = StringToChar(location);
  char *cName = StringToChar(name);

  unsigned int value;
  TwRetVal rv = TwGetUintAttributeFromH5(cFilename, cLocation, cName, &value);

  TwCloseH5(cFilename);
  if (rv != TwSuccess) {
    stop(TranslateReturnValue(rv));
  }

  return value;
}

// GetInt64AttributeFromH5 -----------------------------------------------------
//' Reads a 64-bit integer attribute from the HDF5 file.
//'
//' \code{GetInt64AttributeFromH5} reads a 64-bit integer attribute from the
//' HDF5 file.
//'
//' Used to read attributes not available from \code{GetH5Descriptor}. int64
//' parameters are returned as string. They can be converted to integer64
//' using \code{\link[bit64:as.integer64.character]{bit64::as.integer64()}}.
//'
//' @param Filename Path/filename of the HDF5 file.
//' @param location Location of the group or dataset where the attribute is attached to.
//' @param name Attribute name.
//' @return An string attribute.
//'
//' @examples
//' \dontrun{
//' GetInt64AttributeFromH5("path/to/file.h5", location = , name = )
//' }
//' @export
// [[Rcpp::export]]
CharacterVector GetInt64AttributeFromH5(std::string Filename,
                                        std::string location, std::string name) {

  char *cFilename = StringToChar(Filename);
  char *cLocation = StringToChar(location);
  char *cName = StringToChar(name);

  int64_t value;
  TwRetVal rv = TwGetInt64AttributeFromH5(cFilename, cLocation, cName, &value);

  TwCloseH5(cFilename);
  if (rv != TwSuccess) {
    stop(TranslateReturnValue(rv));
  }

  CharacterVector out(1);
  out[0] = std::to_string(value);
  return out;
}

// GetUint64AttributeFromH5 ----------------------------------------------------
//' Reads an unsigned 64-bit integer attribute from the HDF5 file.
//'
//' \code{GetUint64AttributeFromH5} reads an unsigned 64-bit integer attribute
//' from the HDF5 file.
//'
//' Used to read attributes not available from \code{GetH5Descriptor}. Unsigned
//' int64 parameters are returned as string.
//'
//' @param Filename Path/filename of the HDF5 file.
//' @param location Location of the group or dataset where the attribute is attached to.
//' @param name Attribute name.
//' @return An int64 attribute.
//'
//' @examples
//' \dontrun{
//' GetUint64AttributeFromH5("path/to/file.h5", location = , name = )
//' }
//' @export
// [[Rcpp::export]]
CharacterVector GetUint64AttributeFromH5(std::string Filename,
                                         std::string location, std::string name) {

  char *cFilename = StringToChar(Filename);
  char *cLocation = StringToChar(location);
  char *cName = StringToChar(name);

  uint64_t value;
  TwRetVal rv = TwGetUint64AttributeFromH5(cFilename, cLocation, cName, &value);

  TwCloseH5(cFilename);
  if (rv != TwSuccess) {
    stop(TranslateReturnValue(rv));
  }

  CharacterVector out(1);
  out[0] = std::to_string(value);
  return out;
}

// GetFloatAttributeFromH5 -----------------------------------------------------
//' Reads a float attribute from the HDF5 file.
//'
//' \code{GetFloatAttributeFromH5} reads a float attribute from the HDF5 file.
//'
//' Used to read attributes not available from \code{GetH5Descriptor}.
//'
//' @param Filename Path/filename of the HDF5 file.
//' @param location Location of the group or dataset where the attribute is attached to.
//' @param name Attribute name.
//' @return A numeric attribute.
//'
//' @examples
//' \dontrun{
//' GetFloatAttributeFromH5("path/to/file.h5", location = , name = )
//' }
//' @export
// [[Rcpp::export]]
float GetFloatAttributeFromH5(std::string Filename, std::string location,
                              std::string name) {

  char *cFilename = StringToChar(Filename);
  char *cLocation = StringToChar(location);
  char *cName = StringToChar(name);

  float value;
  TwRetVal rv = TwGetFloatAttributeFromH5(cFilename, cLocation, cName, &value);

  TwCloseH5(cFilename);
  if (rv != TwSuccess) {
    stop(TranslateReturnValue(rv));
  }

  return value;
}

// GetDoubleAttributeFromH5 ----------------------------------------------------
//' Reads a double attribute from the HDF5 file.
//'
//' \code{GetDoubleAttributeFromH5} reads a double attribute from the HDF5 file.
//'
//' Used to read attributes not available from \code{GetH5Descriptor}.
//'
//' @param Filename Path/filename of the HDF5 file.
//' @param location Location of the group or dataset where the attribute is attached to.
//' @param name Attribute name.
//' @return A numeric attribute.
//'
//' @examples
//' \dontrun{
//' GetDoubleAttributeFromH5("path/to/file.h5", location = , name = )
//' }
//' @export
// [[Rcpp::export]]
double GetDoubleAttributeFromH5(std::string Filename, std::string location,
                                std::string name) {

  char *cFilename = StringToChar(Filename);
  char *cLocation = StringToChar(location);
  char *cName = StringToChar(name);

  double value;
  TwRetVal rv = TwGetDoubleAttributeFromH5(cFilename, cLocation, cName, &value);

  TwCloseH5(cFilename);
  if (rv != TwSuccess) {
    stop(TranslateReturnValue(rv));
  }

  return value;
}

// GetStringAttributeFromH5 ----------------------------------------------------
//' Reads a string attribute from the HDF5 file.
//'
//' \code{GetStringAttributeFromH5} reads a string attribute from the HDF5 file.
//'
//' Used to read attributes not available from \code{GetH5Descriptor}.
//'
//' @param Filename Path/filename of the HDF5 file.
//' @param location Location of the group or dataset where the attribute is attached to.
//' @param name Attribute name.
//' @return A string attribute.
//'
//' @examples
//' \dontrun{
//' GetStringAttributeFromH5("path/to/file.h5", location = , name = )
//' }
//' @export
// [[Rcpp::export]]
String GetStringAttributeFromH5(std::string Filename, std::string location,
                                std::string name) {

  char *cFilename = StringToChar(Filename);
  char *cLocation = StringToChar(location);
  char *cName = StringToChar(name);
  std::vector<char> value(8192);

  TwRetVal rv = TwGetStringAttributeFromH5(cFilename, cLocation, cName, &value[0]);

  TwCloseH5(cFilename);
  if (rv != TwSuccess) {
    stop(TranslateReturnValue(rv));
  }

  std::string str(value.begin(), value.end());

  return wrap(str);
}

// SetIntAttributeInH5 ---------------------------------------------------------
//' Writes an integer attribute to the HDF5 file.
//'
//' \code{SetIntAttributeInH5} writes an integer attribute to the HDF5 file.
//'
//' @param Filename Path/filename of the HDF5 file.
//' @param location Location of the group or dataset where the attribute is attached to.
//' @param name Attribute name.
//' @param attribute Integer attribute.
//'
//' @examples
//' \dontrun{
//' SetIntAttributeInH5("path/to/file.h5", location = , name = , attribute = )
//' }
//' @export
// [[Rcpp::export]]
void SetIntAttributeInH5(std::string Filename, std::string location,
                         std::string name, int attribute) {

  char *cFilename = StringToChar(Filename);
  char *cLocation = StringToChar(location);
  char *cName = StringToChar(name);

  TwRetVal rv = TwSetIntAttributeInH5(cFilename, cLocation, cName, attribute);
  TwCloseH5(cFilename);

  if (rv != TwSuccess) {
    stop(TranslateReturnValue(rv));
  }
}

// SetUintAttributeInH5 --------------------------------------------------------
//' Writes an unsigned integer attribute to the HDF5 file.
//'
//' \code{SetUintAttributeInH5} writes an unsigned integer attribute to the HDF5 file.
//'
//' @param Filename Path/filename of the HDF5 file.
//' @param location Location of the group or dataset where the attribute is attached to.
//' @param name Attribute name.
//' @param attribute Unsigned integer attribute (passed as numeric value).
//'
//' @examples
//' \dontrun{
//' SetUintAttributeInH5("path/to/file.h5", location = , name = , attribute = )
//' }
//' @export
// [[Rcpp::export]]
void SetUintAttributeInH5(std::string Filename, std::string location,
                          std::string name, double attribute) {

  char *cFilename = StringToChar(Filename);
  char *cLocation = StringToChar(location);
  char *cName = StringToChar(name);

  TwRetVal rv = TwSetUintAttributeInH5(cFilename, cLocation, cName,
                                       (unsigned int)attribute);
  TwCloseH5(cFilename);

  if (rv != TwSuccess) {
    stop(TranslateReturnValue(rv));
  }
}

// SetInt64AttributeInH5 -------------------------------------------------------
//' Writes an int64 attribute to the HDF5 file.
//'
//' \code{SetInt64AttributeInH5} writes an int64 attribute to the HDF5 file.
//'
//' @param Filename Path/filename of the HDF5 file.
//' @param location Location of the group or dataset where the attribute is attached to.
//' @param name Attribute name.
//' @param attribute Int64 attribute passed as a string.
//'
//' @examples
//' \dontrun{
//' SetInt64AttributeInH5("path/to/file.h5", location = , name = , attribute = )
//' }
//' @export
// [[Rcpp::export]]
void SetInt64AttributeInH5(std::string Filename, std::string location,
                           std::string name, std::string attribute) {

  char *cFilename = StringToChar(Filename);
  char *cLocation = StringToChar(location);
  char *cName = StringToChar(name);

  std::stringstream ss(attribute);
  int64_t int64value;
  ss >> int64value;

  TwRetVal rv = TwSetInt64AttributeInH5(cFilename, cLocation, cName, int64value);
  TwCloseH5(cFilename);

  if (rv != TwSuccess) {
    stop(TranslateReturnValue(rv));
  }
}

// SetUint64AttributeInH5 ------------------------------------------------------
//' Writes an unsigned int64 attribute to the HDF5 file.
//'
//' \code{SetUint64AttributeInH5} writes an unsigned int64 attribute to the HDF5 file.
//'
//' @param Filename Path/filename of the HDF5 file.
//' @param location Location of the group or dataset where the attribute is attached to.
//' @param name Attribute name.
//' @param attribute Unsigned int64 attribute passed as a string.
//'
//' @examples
//' \dontrun{
//' SetUint64AttributeInH5("path/to/file.h5", location = , name = , attribute = )
//' }
//' @export
// [[Rcpp::export]]
void SetUint64AttributeInH5(std::string Filename, std::string location,
                            std::string name, std::string attribute) {

  char *cFilename = StringToChar(Filename);
  char *cLocation = StringToChar(location);
  char *cName = StringToChar(name);

  std::stringstream ss(attribute);
  uint64_t uint64value;
  ss >> uint64value;

  TwRetVal rv = TwSetUint64AttributeInH5(cFilename, cLocation, cName,
                                         uint64value);
  TwCloseH5(cFilename);

  if (rv != TwSuccess) {
    stop(TranslateReturnValue(rv));
  }
}

// SetFloatAttributeInH5 -------------------------------------------------------
//' Writes a float attribute to the HDF5 file.
//'
//' \code{SetFloatAttributeInH5} writes a float attribute to the HDF5 file.
//'
//' @param Filename Path/filename of the HDF5 file.
//' @param location Location of the group or dataset where the attribute is attached to.
//' @param name Attribute name.
//' @param attribute Float attribute.
//'
//' @examples
//' \dontrun{
//' SetFloatAttributeInH5("path/to/file.h5", location = , name = , attribute = )
//' }
//' @export
// [[Rcpp::export]]
void SetFloatAttributeInH5(std::string Filename, std::string location,
                           std::string name, double attribute) {

  char *cFilename = StringToChar(Filename);
  char *cLocation = StringToChar(location);
  char *cName = StringToChar(name);

  TwRetVal rv = TwSetFloatAttributeInH5(cFilename, cLocation, cName,
                                        (float)attribute);
  TwCloseH5(cFilename);

  if (rv != TwSuccess) {
    stop(TranslateReturnValue(rv));
  }
}

// SetDoubleAttributeInH5 ------------------------------------------------------
//' Writes a double attribute to the HDF5 file.
//'
//' \code{SetDoubleAttributeInH5} writes a double attribute to the HDF5 file.
//'
//' @param Filename Path/filename of the HDF5 file.
//' @param location Location of the group or dataset where the attribute is attached to.
//' @param name Attribute name.
//' @param attribute Double attribute.
//'
//' @examples
//' \dontrun{
//' SetDoubleAttributeInH5("path/to/file.h5", location = , name = , attribute = )
//' }
//' @export
// [[Rcpp::export]]
void SetDoubleAttributeInH5(std::string Filename, std::string location,
                            std::string name, double attribute) {

  char *cFilename = StringToChar(Filename);
  char *cLocation = StringToChar(location);
  char *cName = StringToChar(name);

  TwRetVal rv = TwSetDoubleAttributeInH5(cFilename, cLocation, cName, attribute);
  TwCloseH5(cFilename);

  if (rv != TwSuccess) {
    stop(TranslateReturnValue(rv));
  }
}

// SetStringAttributeInH5 ------------------------------------------------------
//' Writes a string attribute to the HDF5 file.
//'
//' \code{SetStringAttributeInH5} writes a string attribute to the HDF5 file.
//'
//' @param Filename Path/filename of the HDF5 file.
//' @param location Location of the group or dataset where the attribute is attached to.
//' @param name Attribute name.
//' @param attribute String attribute (max. 256 characters).
//'
//' @examples
//' \dontrun{
//' SetIntAttributeInH5("path/to/file.h5", location = , name = , attribute = )
//' }
//' @export
// [[Rcpp::export]]
void SetStringAttributeInH5(std::string Filename, std::string location,
                            std::string name, std::string attribute) {

  char *cFilename = StringToChar(Filename);
  char *cLocation = StringToChar(location);
  char *cName = StringToChar(name);
  char *cAttribute = StringToChar(attribute);

  TwRetVal rv = TwSetStringAttributeInH5(cFilename, cLocation, cName, cAttribute);
  TwCloseH5(cFilename);

  if (rv != TwSuccess) {
    stop(TranslateReturnValue(rv));
  }
}

// GetUserDataFromH5 -----------------------------------------------------------
//' Reads user data from the HDF5 file.
//'
//' \code{GetUserDataFromH5} reads a row of user data added to the file using
//' the \code{AddUserData} function.
//'
//'  If you want to access data saved by a registered data source (or a data
//'  source plugin, which uses the same mechanism) use the \code{\link{GetRegUserDataFromH5}}
//'  function.
//'
//' @param Filename Path/filename of the HDF5 file.
//' @param location Location of the group or dataset where the attribute is attached to.
//' @param rowIndex Index of row to read. If index is -1 (default), all rows are read.
//' @return A list containing the user data and data description.
//'
//' @examples
//' \dontrun{
//' GetUserDataFromH5("path/to/file.h5", "/ImageData/ScanData")
//' }
//' @export
// [[Rcpp::export]]
List GetUserDataFromH5(std::string Filename, std::string location, int rowIndex = -1) {

  char *cFilename = StringToChar(Filename);
  char *cLocation = StringToChar(location);

  // get nbrElements
  int nbrElements = 0;
  int rows = 0;
  TwRetVal rv = TwGetUserDataFromH5(cFilename, cLocation, &rows, &nbrElements,
                                    NULL, NULL);
  if (rv != TwValueAdjusted) {
    TwCloseH5(cFilename);
    stop(TranslateReturnValue(rv));
  }

  NumericVector buffer(nbrElements);
  // CharacterVector elementDescription(nbrElements);
  char *elementDescription = new char[256 * nbrElements];

  List result;

  if (rowIndex == -1) {

    NumericMatrix array(rows, nbrElements);

    for (int j = 0; j < rows; j++) {
      rv = TwGetUserDataFromH5(cFilename, cLocation, &j, &nbrElements,
                               &buffer[0], elementDescription);

      if (rv != TwSuccess) {
        delete[] elementDescription;
        TwCloseH5(cFilename);
        stop(TranslateReturnValue(rv));
      }

      for (int i = 0; i < nbrElements; i++) {
        array(j,i) = buffer[i];
      }
    }

    result["data"] = array;

  } else {

    rv = TwGetUserDataFromH5(cFilename, cLocation, &rowIndex, &nbrElements,
                             &buffer[0], elementDescription);

    if (rv != TwSuccess) {
      delete[] elementDescription;
      TwCloseH5(cFilename);
      stop(TranslateReturnValue(rv));
    }

    result["data"] = buffer;

  }

  TwCloseH5(cFilename);

  CharacterVector descriptionArray(nbrElements);
  std::string str(elementDescription, 256 * nbrElements);
  delete[] elementDescription;

  for (int i = 0; i < nbrElements; i++) {
    descriptionArray[i] = str.substr(i*256, 256);
  }

  result["elementDescription"] = descriptionArray;

  return result;
}

// GetAcquisitionLogFromH5 -----------------------------------------------------
//' Reads a single acquisition log entry.
//'
//' \code{GetAcquisitionLogFromH5} reads a single acquisition log entry and
//' returns the timestamp and the log text.
//'
//' The timestamp is the number of 100-nanosecond intervals since January 1,
//' 1601 (UTC). Use \code{bit64::as.integer64(timestamp)} to convert the
//' timestamp string into a 64-bit integer value.
//'
//' @param Filename Path/filename of the HDF5 file.
//' @param index Index of log entry.
//' @return A list containing the timestamp (as a string) and log text.
//'
//' @examples
//' \dontrun{
//' GetAcquisitionLogFromH5("path/to/file.h5", index = )
//' }
//' @export
// [[Rcpp::export]]
List GetAcquisitionLogFromH5(std::string Filename, int index) {

  char *cFilename = StringToChar(Filename);

  int64_t timestamp;
  char *logText = new char[256];

  TwRetVal rv = TwGetAcquisitionLogFromH5(cFilename, index, &timestamp, logText);
  TwCloseH5(cFilename);

  if (rv != TwSuccess) {
    delete[] logText;
    stop(TranslateReturnValue(rv));
  }

  CharacterVector time(1);
  time[0] = std::to_string(timestamp);

  std::string str(logText);

  delete[] logText;

  List result;
  result["timestamp"] = time;
  result["logText"] = str;

  return result;
}

// GetEventListSpectrumFromH5 --------------------------------------------------
//' Reads the events of a spectrum from HDF5 data file.
//'
//' \code{GetEventListSpectrumFromH5} reads the events of a single spectrum
//' given by segment, buf and write indices from HDF5 data file.
//'
//' Note: only uncompressed event lists are supported.
//'
//' @param Filename Path/filename of the HDF5 file.
//' @param segmentIndex Segment index.
//' @param bufIndex Buf index.
//' @param writeIndex Write index.
//' @return A vector containing the event data.
//'
//' @examples
//' \dontrun{
//' GetEventListSpectrumFromH5("path/to/file.h5", segmentIndex = 0,
//' bufIndex = 0, writeIndex = 0)
//' }
//' @export
// [[Rcpp::export]]
SEXP GetEventListSpectrumFromH5(std::string Filename, int segmentIndex,
                                int bufIndex, int writeIndex) {

  char *cFilename = StringToChar(Filename);

  // get bufferSize
  int bufferSize = 0;
  TwRetVal rv = TwGetEventListSpectrumFromH5(cFilename, segmentIndex, bufIndex,
                                             writeIndex, &bufferSize, NULL);
  if ((rv != TwValueAdjusted) & (rv != TwSuccess)) {
    TwCloseH5(cFilename);
    stop(TranslateReturnValue(rv));
  }

  std::vector<unsigned int> buffer(bufferSize);

  rv = TwGetEventListSpectrumFromH5(cFilename, segmentIndex, bufIndex,
                                    writeIndex, &bufferSize, &buffer[0]);
  TwCloseH5(cFilename);
  if (rv != TwSuccess) {
    stop(TranslateReturnValue(rv));
  }

  return wrap(buffer);
}

// H5GetMassCalibPar -----------------------------------------------------------
//' Gets mass calibration parameters from the data file.
//'
//' \code{H5GetMassCalibPar} gets mass calibration parameters from the data file.
//'
//' @param Filename Path/filename of the HDF5 file.
//' @param writeIndex Write index.
//' @return A list containing the calibraion mode and calibration parameters.
//'
//' @examples
//' \dontrun{
//' H5GetMassCalibPar("path/to/file.h5", writeIndex = 0)
//' }
//' @export
// [[Rcpp::export]]
List H5GetMassCalibPar(std::string Filename, int writeIndex) {

  char *cFilename = StringToChar(Filename);

  int segmentIndex = 0;
  int bufIndex = 0;

  // get nbrParams
  int nbrParams = 0;
  TwRetVal rv = TwH5GetMassCalibPar(cFilename, segmentIndex, bufIndex,
                                    writeIndex, NULL, &nbrParams, NULL);
  if (rv != TwValueAdjusted) {
    TwCloseH5(cFilename);
    stop(TranslateReturnValue(rv));
  }

  int mode;
  NumericVector p(nbrParams);

  rv = TwH5GetMassCalibPar(cFilename, segmentIndex, bufIndex,
                           writeIndex, &mode, &nbrParams, &p[0]);
  TwCloseH5(cFilename);
  if (rv != TwSuccess) {
    stop(TranslateReturnValue(rv));
  }

  List result;
  result["mode"] = mode;
  result["p"] = p;

  return result;
}

// H5AddLogEntry ---------------------------------------------------------------
//' Adds an entry to an existing data file.
//'
//' \code{H5AddLogEntry} adds an entry to an existing data file. To add
//' acquisition log entries during a running acquisition use \code{AddLogEntry}.
//'
//' @param Filename Path/filename of the HDF5 file.
//' @param LogEntryText Log text (max. 255 characters).
//' @param LogEntryTime Log entry time (number of 100-nanosecond intervals since
//' January 1, 1601 UTC, Windows FILETIME) passed as a string. Set it to "0" for
//' "now".
//'
//' @export
// [[Rcpp::export]]
void H5AddLogEntry(std::string Filename, std::string LogEntryText,
                   std::string LogEntryTime) {

  char *cFilename = StringToChar(Filename);
  char *cLogEntryText = StringToChar(LogEntryText);

  std::stringstream ss(LogEntryTime);
  uint64_t cTime;
  ss >> cTime;

  TwRetVal rv = TwH5AddLogEntry(cFilename, cLogEntryText, cTime);

  TwCloseH5(cFilename);
  if (rv != TwSuccess) {
    stop(TranslateReturnValue(rv));
  }
}

// H5AddUserDataMultiRow -------------------------------------------------------
//' Adds user data to a data file.
//'
//' \code{H5AddUserDataMultiRow} adds user data to a data file. Creates datasets
//' "Data" and "Info" at \code{Location}.
//'
//' @param filename Path/filename of the HDF5 file.
//' @param location Location of group in HDF5 file where the datasets are created.
//' @param nbrElements Number of elements to store per row (if the dataset
//' already exists this value must be the same as in the file).
//' @param nbrRows Number of rows to store per call to this function (each row
//' contains \code{NbrElements} entries).
//' @param data Vector of length \code{nbrElements*nbrRows} containing the data to be
//' stored in dataset "Data".
//' @param elementDescription Vector of length \code{nbrElements} containing the
//' text description of elements. If \code{ElementDescription} is \code{NULL}
//' the dataset "Info" is not created.
//' @param compressionLevel ZLIB compression level (0-9) for dataset creation.
//' If the dataset at Location already exists this parameter has no effect.
//'
//' @export
// [[Rcpp::export]]
void H5AddUserDataMultiRow(std::string filename, std::string location,
                           int nbrElements, int nbrRows, NumericVector data,
                           Nullable<Rcpp::StringVector> elementDescription = R_NilValue,
                           int compressionLevel = 0) {

  char *cFilename = StringToChar(filename);
  char *cLocation = StringToChar(location);
  char *cElementDescription = new char[256 * nbrElements];

  if (elementDescription.isNotNull()) {
    StringVector strvec(elementDescription); // https://stackoverflow.com/questions/43388698/rcpp-how-can-i-get-the-size-of-a-rcppnullable-numericvector
    for( int i=0; i < nbrElements; i++ ) {
      std::string str(strvec[i]);
      strncpy(&cElementDescription[i*256], str.c_str(), 256);
    }
  } else {
    cElementDescription = NULL;
  }

  TwRetVal rv = TwH5AddUserDataMultiRow(cFilename, cLocation, nbrElements, nbrRows,
                                        cElementDescription, &data[0],
                                        compressionLevel);
  delete[] cElementDescription;

  TwCloseH5(cFilename);
  if (rv != TwSuccess) {
    stop(TranslateReturnValue(rv));
  }
}

// DeleteAttributeInH5 ---------------------------------------------------------
//' Deletes an attribute.
//'
//' \code{DeleteAttributeInH5} deletes an attribute.
//'
//' WARNING: no sanity checking is performed! You can delete attributes that
//' render the data file unusable.
//'
//' @param filename Path/filename of the HDF5 file.
//' @param location Location of the group or dataset the attribute is deleted from.
//' @param name Attribute name.
//'
//' @export
// [[Rcpp::export]]
void DeleteAttributeInH5(std::string filename, std::string location,
                         std::string name) {
#ifdef _WIN32
  char *cFilename = StringToChar(filename);
  char *cLocation = StringToChar(location);
  char *cName = StringToChar(name);

  TwRetVal rv = TwDeleteAttributeInH5(cFilename, cLocation, cName);

  TwCloseH5(cFilename);
  if (rv != TwSuccess) {
    stop(TranslateReturnValue(rv));
  }
#else
  stop("This function is only implemented on Windows.");
#endif
}

// WaitForExclusiveFileAccess --------------------------------------------------
//' Checks whether a file can be opened with exclusive access rights.
//'
//' \code{WaitForExclusiveFileAccess} checks whether a file can be opened with
//' exclusive access rights. This function can be used when opening a data file
//' that just finished recording in order to make sure that the recording
//' application as well as the OS have finished writing to the file. Available
//' only under Windows, returns TwError on other platforms.
//'
//' @param filename Path/filename of the HDF5 file.
//' @param timeoutMs Timeout (in ms) after which the function returns.
//'
//' @export
// [[Rcpp::export]]
void WaitForExclusiveFileAccess(std::string filename, int timeoutMs) {
#ifdef _WIN32
  char *cFilename = StringToChar(filename);

  TwRetVal rv = TwWaitForExclusiveFileAccess(cFilename, timeoutMs);

  TwCloseH5(cFilename);
  if (rv != TwSuccess) {
    stop(TranslateReturnValue(rv));
  }
#else
  stop("This function is only implemented on Windows.");
#endif
}

// WriteNetCdfTimeSeriesFile ---------------------------------------------------
//' Writes a ANDI chromatography file.
//'
//' \code{WriteNetCdfTimeSeriesFile} writes a ANDI chromatography file.
//'
//' @param filename Path/filename of the HDF5 file.
//' @param inject_ts Injection timestamp (number of 100-nanosecond intervals since
//' January 1, 1601 UTC, Windows FILETIME) passed as a string.
//' @param expTitle Experiment title.
//' @param operator_name Operator name.
//' @param company_method_name Company method name.
//' @param source_file_reference Source file reference.
//' @param retention_unit Retention unit, e.g. "[s]".
//' @param detector_unit Detector unit, e.g. "[cps]".
//' @param sample_name Sample name.
//' @param raw_data_table_name Dataset name.
//' @param retention Time axis data.
//' @param ordinate Intensity axis data.
//'
//' @export
// [[Rcpp::export]]
void WriteNetCdfTimeSeriesFile(std::string filename,
                               std::string inject_ts,
                               std::string expTitle,
                               std::string operator_name,
                               std::string company_method_name,
                               std::string source_file_reference,
                               std::string retention_unit,
                               std::string detector_unit,
                               std::string sample_name,
                               std::string raw_data_table_name,
                               NumericVector retention,
                               NumericVector ordinate) {
#ifdef _WIN32
  char *cFilename = StringToChar(filename);
  std::stringstream ss(inject_ts);
  uint64_t cTime;
  ss >> cTime;
  char *cExpTitle = StringToChar(expTitle);
  char *cOperator_name = StringToChar(operator_name);
  char *cCompany_method_name = StringToChar(company_method_name);
  char *cSource_file_reference = StringToChar(source_file_reference);
  char *cRetention_unit = StringToChar(retention_unit);
  char *cDetector_unit = StringToChar(detector_unit);
  char *cSample_name = StringToChar(sample_name);
  char *cRaw_data_table_name = StringToChar(raw_data_table_name);
  int nbrPoints = retention.size();
  if (nbrPoints != ordinate.size()) {
    stop("retention and ordinate must be the same length.");
  }
  std::vector<float> fretention = as<std::vector<float> >(retention);
  std::vector<float> fordinate = as<std::vector<float> >(ordinate);

  TwRetVal rv = TwWriteNetCdfTimeSeriesFile(cFilename, cTime, cExpTitle,
                                            cOperator_name, cCompany_method_name,
                                            cSource_file_reference,
                                            cRetention_unit, cDetector_unit,
                                            cSample_name, cRaw_data_table_name,
                                            nbrPoints, &fretention[0], &fordinate[0]);

  TwCloseH5(cFilename);
  if (rv != TwSuccess) {
    stop(TranslateReturnValue(rv));
  }
#else
  stop("This function is only implemented on Windows.");
#endif
}

// H5SetMassCalib ----------------------------------------------------------------
//' Changes the global mass calibration in the data file.
//'
//' \code{H5SetMassCalib} changes the global mass calibration in the data file.
//' If \code{nbrParams} is 0, the calibration parameters will be determined from
//' \code{mass}, \code{tof} and \code{weight}. If calibration parameters and
//' calibration points are provided, the calibration is given by the parameters
//' (no sanity check is performed whether the points yield the same parameters).
//'
//' \tabular{cl}{
//' mode \tab Mass calibration function \cr
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
//' @param Filename Path/filename of the HDF5 file.
//' @param mode Mass calibration function to use.
//' @param nbrParams Number of mass calibration parameters.
//' @param p Vector with mass calibration parameters.
//' @param mass Vector with mass of the calibration points.
//' @param tof Vector with TOF sample index of the calibration points.
//' @param weight Vector with weight of the calibration points.
//'
//' @export
// [[Rcpp::export]]
void H5SetMassCalib(std::string Filename, int mode, int nbrParams,
                    NumericVector p, NumericVector mass, NumericVector tof,
                    NumericVector weight) {

  char *cFilename = StringToChar(Filename);
  int nbrPoints = mass.size();

  TwRetVal rv = TwH5SetMassCalib(cFilename, mode, nbrParams, &p[0], nbrPoints,
                                 &mass[0], &tof[0], &weight[0]);

  if (rv != TwSuccess) {
    stop(TranslateReturnValue(rv));
  }
}

// H5SetMassCalib2 ----------------------------------------------------------------
//' Changes the global mass calibration in the data file.
//'
//' \code{H5SetMassCalib2} changes the global mass calibration in the data file.
//' If \code{nbrParams} is 0, the calibration parameters will be determined from
//' \code{mass}, \code{tof} and \code{weight}. If calibration parameters and
//' calibration points are provided, the calibration is given by the parameters
//' (no sanity check is performed whether the points yield the same parameters).
//'
//' \tabular{cl}{
//' mode \tab Mass calibration function \cr
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
//' @param Filename Path/filename of the HDF5 file.
//' @param mode Mass calibration function to use.
//' @param nbrParams Number of mass calibration parameters.
//' @param p Vector with mass calibration parameters.
//' @param mass Vector with mass of the calibration points.
//' @param tof Vector with TOF sample index of the calibration points.
//' @param weight Vector with weight of the calibration points.
//'
//' @export
// [[Rcpp::export]]
void H5SetMassCalib2(std::string Filename, int mode, int nbrParams,
                    NumericVector p, NumericVector mass, NumericVector tof,
                    NumericVector weight) {

  char *cFilename = StringToChar(Filename);
  int nbrPoints = mass.size();

  TwRetVal rv = TwH5SetMassCalib2(cFilename, mode, nbrParams, &p[0], nbrPoints,
                                 &mass[0], &tof[0], &weight[0]);

  if (rv != TwSuccess) {
    stop(TranslateReturnValue(rv));
  }
}

// H5SetMassCalibEx ----------------------------------------------------------------
//' Changes the global mass calibration in the data file.
//'
//' \code{H5SetMassCalibEx} changes the global mass calibration in the data file.
//' If \code{nbrParams} is 0, the calibration parameters will be determined from
//' \code{mass}, \code{tof} and \code{weight}. If calibration parameters and
//' calibration points are provided, the calibration is given by the parameters
//' (no sanity check is performed whether the points yield the same parameters).
//' Labels to identify compound names/formulas used for calibration have a
//' maximum length of 255 characters.
//'
//' \tabular{cl}{
//' mode \tab Mass calibration function \cr
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
//' @param Filename Path/filename of the HDF5 file.
//' @param mode Mass calibration function to use.
//' @param nbrParams Number of mass calibration parameters.
//' @param p Vector with mass calibration parameters.
//' @param mass Vector with mass of the calibration points.
//' @param tof Vector with TOF sample index of the calibration points.
//' @param weight Vector with weight of the calibration points.
//' @param label Vector with labels/names/sum formula of the calibration points.
//'
//' @export
// [[Rcpp::export]]
void H5SetMassCalibEx(std::string Filename, int mode, int nbrParams,
                    NumericVector p, NumericVector mass, NumericVector tof,
                    NumericVector weight, StringVector label) {

  char *cFilename = StringToChar(Filename);
  int nbrPoints = mass.size();

  if (nbrPoints != label.size()) {
    stop("mass, tof, weight and label must be the same length.");
  }
  char *cLabel = new char[256 * nbrPoints];
  for( int i=0; i < nbrPoints; i++ ) {
    std::string str(label[i]);
    strncpy(&cLabel[i*256], str.c_str(), 256);
  }

  TwRetVal rv = TwH5SetMassCalibEx(cFilename, mode, nbrParams, &p[0], nbrPoints,
                                 &mass[0], &tof[0], &weight[0], cLabel);
  delete[] cLabel;

  if (rv != TwSuccess) {
    stop(TranslateReturnValue(rv));
  }
}

// H5SetMassCalib2Ex ----------------------------------------------------------------
//' Changes the global mass calibration in the data file.
//'
//' \code{H5SetMassCalib2Ex} changes the global mass calibration in the data file.
//' If \code{nbrParams} is 0, the calibration parameters will be determined from
//' \code{mass}, \code{tof} and \code{weight}. If calibration parameters and
//' calibration points are provided, the calibration is given by the parameters
//' (no sanity check is performed whether the points yield the same parameters).
//' Labels to identify compound names/formulas used for calibration have a
//' maximum length of 255 characters.
//'
//' \tabular{cl}{
//' mode \tab Mass calibration function \cr
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
//' @param Filename Path/filename of the HDF5 file.
//' @param mode Mass calibration function to use.
//' @param nbrParams Number of mass calibration parameters.
//' @param p Vector with mass calibration parameters.
//' @param mass Vector with mass of the calibration points.
//' @param tof Vector with TOF sample index of the calibration points.
//' @param weight Vector with weight of the calibration points.
//' @param label Vector with labels/names/sum formula of the calibration points.
//'
//' @export
// [[Rcpp::export]]
void H5SetMassCalib2Ex(std::string Filename, int mode, int nbrParams,
                      NumericVector p, NumericVector mass, NumericVector tof,
                      NumericVector weight, StringVector label) {

  char *cFilename = StringToChar(Filename);
  int nbrPoints = mass.size();

  if (nbrPoints != label.size()) {
    stop("mass, tof, weight and label must be the same length.");
  }
  char *cLabel = new char[256 * nbrPoints];
  for( int i=0; i < nbrPoints; i++ ) {
    std::string str(label[i]);
    strncpy(&cLabel[i*256], str.c_str(), 256);
  }

  TwRetVal rv = TwH5SetMassCalib2Ex(cFilename, mode, nbrParams, &p[0], nbrPoints,
                                   &mass[0], &tof[0], &weight[0], cLabel);
  delete[] cLabel;

  if (rv != TwSuccess) {
    stop(TranslateReturnValue(rv));
  }
}

// H5SetMassCalibDynamic -------------------------------------------------------
//' Stores dynamic mass calibration for a given spectrum in the data file.
//'
//' \code{H5SetMassCalibDynamic} stores dynamic mass calibration for a given
//' spectrum in the data file.
//'
//' This function can be used to delete the dynamic calibration information by
//' passing the special parameter set: \code{writeIndex} = -1 and \code{par} =
//' \code{NULL}.
//'
//' @param filename Path/filename of the HDF5 file.
//' @param writeIndex Write index of the calibration to store.
//' @param par Numeric vector holding the actual calibration values. It
//' must be the same number of parameters as the global mass calibration.
//'
//' @export
// [[Rcpp::export]]
void H5SetMassCalibDynamic(std::string filename, int writeIndex,
                           Nullable<Rcpp::NumericVector> par) {
  char *cFilename = StringToChar(filename);
  int segmentIndex = -1;
  int bufIndex = -1;
  int nbrParams;
  double *p_par;
  if (par.isNotNull()) {
    NumericVector par_(par);
    nbrParams = par_.size();
    p_par = &par_[0];
  } else {
    nbrParams = 0;
    p_par = nullptr;
  }

  TwRetVal rv = TwH5SetMassCalibDynamic(cFilename, segmentIndex, bufIndex,
                                        writeIndex, nbrParams, p_par, nbrParams,
                                        p_par);
  if (rv != TwSuccess) {
    stop(TranslateReturnValue(rv));
  }
}

// ChangePeakTable -------------------------------------------------------------
//' Changes the PeakTable and recomputes the PeakData.
//'
//' \code{ChangePeakTable} changes the entries in the dataset \code{/PeakData/PeakTable}
//' and recomputes \code{/PeakData/PeakData} accordingly.
//'
//' Depending on data file size recomputing the peak data can take a long time.
//'
//' @param Filename Path/filename of the HDF5 file.
//' @param PeakPar List with the new peak paramters \emph{label}, \emph{mass},
//' \emph{loMass} and \emph{hiMass}.
//' @param compressionLevel Compression Level (1..9).
//'
//' @export
// [[Rcpp::export]]
void ChangePeakTable(std::string Filename, List PeakPar, int compressionLevel) {

  char *cFilename = StringToChar(Filename);

  CharacterVector label = PeakPar["label"];
  NumericVector mass = PeakPar["mass"];
  NumericVector loMass = PeakPar["loMass"];
  NumericVector hiMass = PeakPar["hiMass"];

  int nbrNewPeakPar = mass.size();
  TPeakPar *newPeakPar = new TPeakPar[nbrNewPeakPar];

  for (int i=0; i<nbrNewPeakPar; ++i) {
    std::string str(label[i]);
    strncpy(newPeakPar[i].label, str.c_str(), 64);
    newPeakPar[i].mass = (float)mass[i];
    newPeakPar[i].loMass = (float)loMass[i];
    newPeakPar[i].hiMass = (float)hiMass[i];
  }

  TwRetVal rv = TwChangePeakTable(cFilename, newPeakPar, nbrNewPeakPar,
                                  compressionLevel, NULL);
  delete[] newPeakPar;

  if (rv != TwSuccess) {
    stop(TranslateReturnValue(rv));
  }
}

// ChangePeakTableFromFile -----------------------------------------------------
//' Changes the PeakTable and recomputes the PeakData.
//'
//' \code{ChangePeakTableFromFile} changes the entries in the dataset \code{/PeakData/PeakTable}
//' and recomputes \code{/PeakData/PeakData} accordingly. This is a convenience
//' function for \code{\link{ChangePeakTable}}. Instead of specifiying directly
//' peak parameters the information is taken from a mass table (text file) or a
//' Tofwerk HDF5 file.
//'
//' Depending on data file size  recomputing the peak data can take a long time.
//'
//' @param Filename Path/filename of the HDF5 file.
//' @param massTable Path/filename to mass table or HDF5 file.
//' @param compressionLevel Compression Level (1..9).
//'
//' @export
// [[Rcpp::export]]
void ChangePeakTableFromFile(std::string Filename, std::string massTable,
                             int compressionLevel) {

  char *cFilename = StringToChar(Filename);
  char *cMassTable = StringToChar(massTable);

  TwRetVal rv = TwChangePeakTableFromFile(cFilename, cMassTable,
                                          compressionLevel, NULL);

  if (rv != TwSuccess) {
    stop(TranslateReturnValue(rv));
  }
}

// Not implemented -------------------------------------------------------------
// TwGetBufWriteProfileFromH5_2
// TwChangePeakTable2
// TwChangePeakFromFile2
// TwProgressCallback
// TwProgressCallback2
// TwChangePeakDataInit
// TwChangePeakDataWrite
// TwChangePeakDataFinalize
// TwReadRawData
// TwGetEventListDataFromH5
// TwGetEventListBlobFromH5
// TwFreeEventListData
// TwFreeEventListData2
// TwMultiPeakFitIntegration
// TwGenerateSegmentProfilesFromEventList
// TwH5MakePaletteImage
// TwH5MakeTrueColorImage
