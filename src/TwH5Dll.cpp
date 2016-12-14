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
//' See \emph{/doc/TwH5Dll.htm} for more details.
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
SEXP GetH5Descriptor(SEXP Filename) {

  char *cFilename = RtoCstring(Filename);

  TwH5Desc descriptor;
  TwRetVal rv = TwGetH5Descriptor(cFilename, &descriptor);

  TwCloseH5(cFilename);

  if (rv != TwSuccess) {
    return TwRetValString(rv);
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
//' @param Filename Path/filename of the HDF5 file.
//'
//' This function is called internally by all \code{Get..FromH5} functions, so
//' it is usually not necessary to call \code{CloseH5} explicitely.
//'
//' @examples
//' \dontrun{
//' CloseH5("path/to/file.h5")
//' }
//' @export
// [[Rcpp::export]]
SEXP CloseH5(SEXP Filename) {

  char *cFilename = RtoCstring(Filename);

  TwRetVal rv = TwCloseH5(cFilename);

  return TwRetValString(rv);
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
SEXP CloseAll() {

  TwRetVal rv = TwCloseAll();

  return TwRetValString(rv);
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
SEXP GetSumSpectrumFromH5(SEXP Filename, bool Normalize = false) {

  char *cFilename = RtoCstring(Filename);

  TwH5Desc descriptor;
  TwRetVal rv = TwGetH5Descriptor(cFilename, &descriptor);
  if (rv != TwSuccess) {
    TwCloseH5(cFilename);
    return TwRetValString(rv);
  }

  NumericVector Spectrum(descriptor.nbrSamples);

  rv = TwGetSumSpectrumFromH5(cFilename, &Spectrum[0], Normalize);

  TwCloseH5(cFilename);

  if (rv != TwSuccess) {
    return TwRetValString(rv);
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
SEXP GetTofSpectrumFromH5(SEXP Filename, int SegmentIndex, int SegmentEndIndex,
                          int BufIndex, int BufEndIndex, int WriteIndex,
                          int WriteEndIndex, bool BufWriteLinked = false,
                          bool Normalize = true) {

  char *cFilename = RtoCstring(Filename);

  TwH5Desc descriptor;
  TwRetVal rv = TwGetH5Descriptor(cFilename, &descriptor);
  if (rv != TwSuccess) {
    TwCloseH5(cFilename);
    return TwRetValString(rv);
  }

  std::vector<float> Spectrum(descriptor.nbrSamples);

  rv = TwGetTofSpectrumFromH5(cFilename, &Spectrum[0], SegmentIndex,
                              SegmentEndIndex, BufIndex, BufEndIndex,
                              WriteIndex, WriteEndIndex, BufWriteLinked,
                              Normalize);
  TwCloseH5(cFilename);

  if (rv != TwSuccess) {
    return TwRetValString(rv);
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
SEXP GetTofSpectrum2FromH5(SEXP Filename, int SegmentIndex, int SegmentEndIndex,
                          int BufIndex, int BufEndIndex, int WriteIndex,
                          int WriteEndIndex, bool BufWriteLinked = false,
                          bool Normalize = true) {

  char *cFilename = RtoCstring(Filename);

  TwH5Desc descriptor;
  TwRetVal rv = TwGetH5Descriptor(cFilename, &descriptor);
  if (rv != TwSuccess) {
    TwCloseH5(cFilename);
    return TwRetValString(rv);
  }

  std::vector<float> Spectrum(descriptor.nbrSamples);

  rv = TwGetTofSpectrum2FromH5(cFilename, &Spectrum[0], SegmentIndex,
                               SegmentEndIndex, BufIndex, BufEndIndex,
                               WriteIndex, WriteEndIndex, BufWriteLinked,
                               Normalize);
  TwCloseH5(cFilename);

  if (rv != TwSuccess) {
    return TwRetValString(rv);
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
SEXP GetStickSpectrumFromH5(SEXP Filename, int SegmentIndex, int SegmentEndIndex,
                           int BufIndex, int BufEndIndex, int WriteIndex,
                           int WriteEndIndex, bool BufWriteLinked = false,
                           bool Normalize = true) {

  char *cFilename = RtoCstring(Filename);

  TwH5Desc descriptor;
  TwRetVal rv = TwGetH5Descriptor(cFilename, &descriptor);
  if (rv != TwSuccess) {
    TwCloseH5(cFilename);
    return TwRetValString(rv);
  }

  std::vector<float> Spectrum(descriptor.nbrPeaks);

  rv = TwGetStickSpectrumFromH5(cFilename, &Spectrum[0], SegmentIndex,
                                SegmentEndIndex, BufIndex, BufEndIndex,
                                WriteIndex, WriteEndIndex, BufWriteLinked,
                                Normalize);
  TwCloseH5(cFilename);

  if (rv != TwSuccess) {
    return TwRetValString(rv);
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
SEXP GetStickSpectrum2FromH5(SEXP Filename, int SegmentIndex, int SegmentEndIndex,
                            int BufIndex, int BufEndIndex, int WriteIndex,
                            int WriteEndIndex, bool BufWriteLinked = false,
                            bool Normalize = true) {

  char *cFilename = RtoCstring(Filename);

  TwH5Desc descriptor;
  TwRetVal rv = TwGetH5Descriptor(cFilename, &descriptor);
  if (rv != TwSuccess) {
    TwCloseH5(cFilename);
    return TwRetValString(rv);
  }

  std::vector<float> Spectrum(descriptor.nbrPeaks);

  rv = TwGetStickSpectrum2FromH5(cFilename, &Spectrum[0], SegmentIndex,
                                SegmentEndIndex, BufIndex, BufEndIndex,
                                WriteIndex, WriteEndIndex, BufWriteLinked,
                                Normalize);
  TwCloseH5(cFilename);

  if (rv != TwSuccess) {
    return TwRetValString(rv);
  }

  return wrap(Spectrum);
}

// GetPeakParametersFromH5 -----------------------------------------------------
//' Peak parameters from HDF5 data file.
//'
//' \code{GetPeakParametersFromH5} reads peak parameters from the data file.
//'
//' @param Filename Path/filename of the HDF5 file.
//' @param PeakIndex Index of peak. If index is -1 (default), peak parameters of all
//' peaks are read.
//' @return A list with the peak paramters \emph{label}, \emph{mass}, \emph{loMass} and \emph{hiMass}.
//'
//' @examples
//' \dontrun{
//' GetPeakParametersFromH5("path/to/file.h5")
//' }
//' @export
// [[Rcpp::export]]
SEXP GetPeakParametersFromH5(SEXP Filename, int PeakIndex = -1) {

  char *cFilename = RtoCstring(Filename);

  if (PeakIndex == -1) {

    TwH5Desc descriptor;
    TwRetVal rv = TwGetH5Descriptor(cFilename, &descriptor);
    if (rv != TwSuccess) {
      TwCloseH5(cFilename);
      return TwRetValString(rv);
    }
    TPeakPar *PeakPar = new TPeakPar[descriptor.nbrPeaks];

    rv = TwGetPeakParametersFromH5(cFilename, PeakPar, PeakIndex);
    TwCloseH5(cFilename);
    if (rv != TwSuccess) {
      delete[] PeakPar;
      return TwRetValString(rv);
    }

    List result;
    CharacterVector label(descriptor.nbrPeaks-1);
    NumericVector mass(descriptor.nbrPeaks-1);
    NumericVector loMass(descriptor.nbrPeaks-1);
    NumericVector hiMass(descriptor.nbrPeaks-1);

    for (int i=0; i<descriptor.nbrPeaks-1; ++i) {
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
      return TwRetValString(rv);
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
SEXP GetBufTimeFromH5(SEXP Filename, int BufIndex, int WriteIndex) {

  char *cFilename = RtoCstring(Filename);

  NumericVector BufTime(1);

  TwRetVal rv = TwGetBufTimeFromH5(cFilename, &BufTime[0], BufIndex, WriteIndex);

  TwCloseH5(cFilename);

  if (rv != TwSuccess) {
    return TwRetValString(rv);
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
//' values of writeIndex).
//' @return A vector containing the x-axis values.
//'
//' @examples
//' \dontrun{
//' GetSpecXaxisFromH5("path/to/file.h5")
//' }
//' @export
// [[Rcpp::export]]
SEXP GetSpecXaxisFromH5(SEXP Filename, int Type, int writeIndex) {

  char *cFilename = RtoCstring(Filename);

  TwH5Desc descriptor;
  TwRetVal rv = TwGetH5Descriptor(cFilename, &descriptor);
  if (rv != TwSuccess) {
    TwCloseH5(cFilename);
    return TwRetValString(rv);
  }

  NumericVector SpecAxis(descriptor.nbrSamples);

  rv = TwGetSpecXaxisFromH5(cFilename, &SpecAxis[0], Type, NULL, 0.0, writeIndex);

  TwCloseH5(cFilename);

  if (rv != TwSuccess) {
    return TwRetValString(rv);
  }

  return SpecAxis;
}

// GetSegmentProfileFromH5 -----------------------------------------------------
//' Segment profile from HDF5 data file.
//'
//' \code{GetSegmentProfileFromH5} reads a segment profile for a given peak (or
//' all peaks) and buf and write slice.
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
SEXP GetSegmentProfileFromH5(SEXP Filename, int PeakIndex, int BufStartIndex,
                             int BufEndIndex, int WriteStartIndex,
                             int WriteEndIndex, bool BufWriteLinked = false) {

  char *cFilename = RtoCstring(Filename);

  TwH5Desc descriptor;
  TwRetVal rv = TwGetH5Descriptor(cFilename, &descriptor);
  if (rv != TwSuccess) {
    TwCloseH5(cFilename);
    return TwRetValString(rv);
  }

  if (PeakIndex != -1) {

    std::vector<float> SegmentProfile(descriptor.nbrSegments);

    rv = TwGetSegmentProfileFromH5(cFilename, &SegmentProfile[0], PeakIndex,
                                   BufStartIndex, BufEndIndex, WriteStartIndex,
                                   WriteEndIndex, BufWriteLinked);
    TwCloseH5(cFilename);
    if (rv != TwSuccess) {
      return TwRetValString(rv);
    }

    return wrap(SegmentProfile);

  } else {

    std::vector<float> SegmentProfile(descriptor.nbrSegments*descriptor.nbrPeaks);

    rv = TwGetSegmentProfileFromH5(cFilename, &SegmentProfile[0], PeakIndex,
                                   BufStartIndex, BufEndIndex, WriteStartIndex,
                                   WriteEndIndex, BufWriteLinked);
    TwCloseH5(cFilename);
    if (rv != TwSuccess) {
      return TwRetValString(rv);
    }

    return wrap(SegmentProfile);
  }
}

// GetSegmentProfile2FromH5 ----------------------------------------------------
//' Segment profile from HDF5 data file.
//'
//' \code{GetSegmentProfile2FromH5} reads a segment profile for a given peak (or
//' all peaks) and buf and write slice.
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
SEXP GetSegmentProfile2FromH5(SEXP Filename, int PeakIndex, int BufStartIndex,
                             int BufEndIndex, int WriteStartIndex,
                             int WriteEndIndex, bool BufWriteLinked = false) {

  char *cFilename = RtoCstring(Filename);

  TwH5Desc descriptor;
  TwRetVal rv = TwGetH5Descriptor(cFilename, &descriptor);
  if (rv != TwSuccess) {
    TwCloseH5(cFilename);
    return TwRetValString(rv);
  }

  if (PeakIndex != -1) {

    std::vector<float> SegmentProfile(descriptor.nbrSegments);

    rv = TwGetSegmentProfile2FromH5(cFilename, &SegmentProfile[0], PeakIndex,
                                    BufStartIndex, BufEndIndex, WriteStartIndex,
                                    WriteEndIndex, BufWriteLinked);
    TwCloseH5(cFilename);
    if (rv != TwSuccess) {
      return TwRetValString(rv);
    }

    return wrap(SegmentProfile);

  } else {

    std::vector<float> SegmentProfile(descriptor.nbrSegments*descriptor.nbrPeaks);

    rv = TwGetSegmentProfile2FromH5(cFilename, &SegmentProfile[0], PeakIndex,
                                    BufStartIndex, BufEndIndex, WriteStartIndex,
                                    WriteEndIndex, BufWriteLinked);
    TwCloseH5(cFilename);
    if (rv != TwSuccess) {
      return TwRetValString(rv);
    }

    return wrap(SegmentProfile);
  }
}

// GetBufWriteProfileFromH5 ----------------------------------------------------
//' Gets a linked buf/write profile.
//'
//' \code{GetBufWriteProfileFromH5} gets a linked buf/write profile for a given
//' peak (or all peaks) and segment slice.
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
SEXP GetBufWriteProfileFromH5(SEXP Filename, int PeakIndex, int SegmentStartIndex,
                              int SegmentEndIndex) {

  char *cFilename = RtoCstring(Filename);

  TwH5Desc descriptor;
  TwRetVal rv = TwGetH5Descriptor(cFilename, &descriptor);
  if (rv != TwSuccess) {
    TwCloseH5(cFilename);
    return TwRetValString(rv);
  }

  if (PeakIndex != -1) {

    std::vector<float> Profile(descriptor.nbrBufs*descriptor.nbrWrites);

    rv = TwGetBufWriteProfileFromH5(cFilename, &Profile[0], PeakIndex,
                                    SegmentStartIndex, SegmentEndIndex);
    TwCloseH5(cFilename);
    if (rv != TwSuccess) {
      return TwRetValString(rv);
    }

    return wrap(Profile);

  } else {

    std::vector<float> Profile(descriptor.nbrBufs*descriptor.nbrWrites*
      descriptor.nbrPeaks);

    rv = TwGetBufWriteProfileFromH5(cFilename, &Profile[0], PeakIndex,
                                    SegmentStartIndex, SegmentEndIndex);
    TwCloseH5(cFilename);
    if (rv != TwSuccess) {
      return TwRetValString(rv);
    }

    return wrap(Profile);
  }
}

// GetBufWriteProfile2FromH5 ----------------------------------------------------
//' Gets a linked buf/write profile.
//'
//' \code{GetBufWriteProfile2FromH5} gets a linked buf/write profile for a given
//' peak (or all peaks) and segment slice.
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
SEXP GetBufWriteProfile2FromH5(SEXP Filename, int PeakIndex, int SegmentStartIndex,
                              int SegmentEndIndex) {

  char *cFilename = RtoCstring(Filename);

  TwH5Desc descriptor;
  TwRetVal rv = TwGetH5Descriptor(cFilename, &descriptor);
  if (rv != TwSuccess) {
    TwCloseH5(cFilename);
    return TwRetValString(rv);
  }

  if (PeakIndex != -1) {

    std::vector<float> Profile(descriptor.nbrBufs*descriptor.nbrWrites);

    rv = TwGetBufWriteProfile2FromH5(cFilename, &Profile[0], PeakIndex,
                                    SegmentStartIndex, SegmentEndIndex);
    TwCloseH5(cFilename);
    if (rv != TwSuccess) {
      return TwRetValString(rv);
    }

    return wrap(Profile);

  } else {

    std::vector<float> Profile(descriptor.nbrBufs*descriptor.nbrWrites*
      descriptor.nbrPeaks);

    rv = TwGetBufWriteProfile2FromH5(cFilename, &Profile[0], PeakIndex,
                                    SegmentStartIndex, SegmentEndIndex);
    TwCloseH5(cFilename);
    if (rv != TwSuccess) {
      return TwRetValString(rv);
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
SEXP GetRegUserDataSourcesFromH5(SEXP Filename) {

  char *cFilename = RtoCstring(Filename);

  // get nbrSources
  int nbrSources = 0;
  TwRetVal rv = TwGetRegUserDataSourcesFromH5(cFilename, &nbrSources, NULL,
                                              NULL, NULL, NULL);
  if (rv != TwValueAdjusted) {
    TwCloseH5(cFilename);
    return TwRetValString(rv);
  }

  char *sourceLocation = new char[256 * nbrSources];
  memset(sourceLocation, 0, 256 * nbrSources);

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
    return TwRetValString(rv);
  }

  CharacterVector locationArray(nbrSources);
  std::string str(sourceLocation, sourceLocation + 256 * nbrSources);

  for (int i = 0; i < nbrSources; ++i) {
    locationArray[i] = str.substr(i*256, 256);
  }
  delete[] sourceLocation;

  LogicalVector hasDescArray(hasDesc, hasDesc + nbrSources);
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
SEXP GetRegUserDataFromH5(SEXP Filename, SEXP location, int bufIndex,
                          int writeIndex, bool readDescription = true) {

  char *cFilename = RtoCstring(Filename);
  char *cLocation = RtoCstring(location);

  // get bufLength
  int bufLength = 0;
  TwRetVal rv = TwGetRegUserDataFromH5(cFilename, cLocation, bufIndex,
                                       writeIndex, &bufLength, NULL, NULL);
  if (rv != TwValueAdjusted) {
    TwCloseH5(cFilename);
    return TwRetValString(rv);
  }

  NumericVector buffer(bufLength);

  if (readDescription) {

    char *description = new char[256 * bufLength];
    memset(description, 0, 256 * bufLength);

    rv = TwGetRegUserDataFromH5(cFilename, cLocation, bufIndex, writeIndex,
                                &bufLength, &buffer[0], description);
    TwCloseH5(cFilename);
    if (rv != TwSuccess) {
      delete[] description;
      return TwRetValString(rv);
    }

    CharacterVector descriptionArray(bufLength);
    std::string str(description, description + 256 * bufLength);
    delete[] description;

    for (int i = 0; i < bufLength; ++i) {
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
      return TwRetValString(rv);
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
//' @export
// [[Rcpp::export]]
SEXP GetTofData(SEXP Filename, int sampleOffset, int sampleCount, int segOffset,
                int segCount, int bufOffset, int bufCount, int writeOffset,
                int writeCount) {

  char *cFilename = RtoCstring(Filename);

  int n = sampleCount*segCount*bufCount*writeCount;

  std::vector<float> dataBuffer(n);

  TwRetVal rv = TwGetTofData(cFilename, sampleOffset, sampleCount, segOffset,
                             segCount, bufOffset, bufCount, writeOffset,
                             writeCount, &dataBuffer[0]);
  TwCloseH5(cFilename);

  if (rv != TwSuccess) {
    return TwRetValString(rv);
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
//' @export
// [[Rcpp::export]]
SEXP GetTofData2(SEXP Filename, int sampleOffset, int sampleCount, int segOffset,
                int segCount, int bufOffset, int bufCount, int writeOffset,
                int writeCount) {

  char *cFilename = RtoCstring(Filename);

  int n = sampleCount*segCount*bufCount*writeCount;

  std::vector<float> dataBuffer(n);

  TwRetVal rv = TwGetTofData2(cFilename, sampleOffset, sampleCount, segOffset,
                             segCount, bufOffset, bufCount, writeOffset,
                             writeCount, &dataBuffer[0]);
  TwCloseH5(cFilename);

  if (rv != TwSuccess) {
    return TwRetValString(rv);
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
//' @export
// [[Rcpp::export]]
SEXP GetPeakData(SEXP Filename, int peakOffset, int peakCount, int segOffset,
                 int segCount, int bufOffset, int bufCount, int writeOffset,
                 int writeCount) {

  char *cFilename = RtoCstring(Filename);

  int n = peakCount*segCount*bufCount*writeCount;

  std::vector<float> dataBuffer(n);

  TwRetVal rv = TwGetPeakData(cFilename, peakOffset, peakCount, segOffset,
                              segCount, bufOffset, bufCount, writeOffset,
                              writeCount, &dataBuffer[0]);
  TwCloseH5(cFilename);

  if (rv != TwSuccess) {
    return TwRetValString(rv);
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
//' @export
// [[Rcpp::export]]
SEXP GetPeakData2(SEXP Filename, int peakOffset, int peakCount, int segOffset,
                 int segCount, int bufOffset, int bufCount, int writeOffset,
                 int writeCount) {

  char *cFilename = RtoCstring(Filename);

  int n = peakCount*segCount*bufCount*writeCount;

  std::vector<float> dataBuffer(n);

  TwRetVal rv = TwGetPeakData2(cFilename, peakOffset, peakCount, segOffset,
                              segCount, bufOffset, bufCount, writeOffset,
                              writeCount, &dataBuffer[0]);
  TwCloseH5(cFilename);

  if (rv != TwSuccess) {
    return TwRetValString(rv);
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
//' @export
// [[Rcpp::export]]
SEXP GetTimingData(SEXP Filename, int bufOffset, int bufCount, int writeOffset,
                 int writeCount) {

  char *cFilename = RtoCstring(Filename);

  int n = bufCount*writeCount;

  NumericVector dataBuffer(n);

  TwRetVal rv = TwGetTimingData(cFilename, bufOffset, bufCount, writeOffset,
                              writeCount, &dataBuffer[0]);
  TwCloseH5(cFilename);

  if (rv != TwSuccess) {
    return TwRetValString(rv);
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
SEXP GetIntAttributeFromH5(SEXP Filename, SEXP location, SEXP name) {

  char *cFilename = RtoCstring(Filename);
  char *cLocation = RtoCstring(location);
  char *cName = RtoCstring(name);

 int value;
  TwRetVal rv = TwGetIntAttributeFromH5(cFilename, cLocation, cName, &value);

  TwCloseH5(cFilename);
  if (rv != TwSuccess) {
    return TwRetValString(rv);
  }

  return wrap(value);
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
SEXP GetUintAttributeFromH5(SEXP Filename, SEXP location, SEXP name) {

  char *cFilename = RtoCstring(Filename);
  char *cLocation = RtoCstring(location);
  char *cName = RtoCstring(name);

  unsigned int value;
  TwRetVal rv = TwGetUintAttributeFromH5(cFilename, cLocation, cName, &value);

  TwCloseH5(cFilename);
  if (rv != TwSuccess) {
    return TwRetValString(rv);
  }

  return wrap(value);
}

// GetInt64AttributeFromH5 -----------------------------------------------------
//' Reads a 64-bit integer attribute from the HDF5 file.
//'
//' \code{GetInt64AttributeFromH5} reads a 64-bit integer attribute from the
//' HDF5 file.
//'
//' Used to read attributes not available from \code{GetH5Descriptor}. int64
//' parameters are returned as string. They can be converted to integer64
//' using \code{\link[bit64]{as.integer64}}.
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
SEXP GetInt64AttributeFromH5(SEXP Filename, SEXP location, SEXP name) {

  char *cFilename = RtoCstring(Filename);
  char *cLocation = RtoCstring(location);
  char *cName = RtoCstring(name);

  int64_t value;
  TwRetVal rv = TwGetInt64AttributeFromH5(cFilename, cLocation, cName, &value);

  TwCloseH5(cFilename);
  if (rv != TwSuccess) {
    return TwRetValString(rv);
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
SEXP GetUint64AttributeFromH5(SEXP Filename, SEXP location, SEXP name) {

  char *cFilename = RtoCstring(Filename);
  char *cLocation = RtoCstring(location);
  char *cName = RtoCstring(name);

  uint64_t value;
  TwRetVal rv = TwGetUint64AttributeFromH5(cFilename, cLocation, cName, &value);

  TwCloseH5(cFilename);
  if (rv != TwSuccess) {
    return TwRetValString(rv);
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
SEXP GetFloatAttributeFromH5(SEXP Filename, SEXP location, SEXP name) {

  char *cFilename = RtoCstring(Filename);
  char *cLocation = RtoCstring(location);
  char *cName = RtoCstring(name);

  float value;
  TwRetVal rv = TwGetFloatAttributeFromH5(cFilename, cLocation, cName, &value);

  TwCloseH5(cFilename);
  if (rv != TwSuccess) {
    return TwRetValString(rv);
  }

  return wrap(value);
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
SEXP GetDoubleAttributeFromH5(SEXP Filename, SEXP location, SEXP name) {

  char *cFilename = RtoCstring(Filename);
  char *cLocation = RtoCstring(location);
  char *cName = RtoCstring(name);

  double value;
  TwRetVal rv = TwGetDoubleAttributeFromH5(cFilename, cLocation, cName, &value);

  TwCloseH5(cFilename);
  if (rv != TwSuccess) {
    return TwRetValString(rv);
  }

  return wrap(value);
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
SEXP GetStringAttributeFromH5(SEXP Filename, SEXP location, SEXP name) {

  char *cFilename = RtoCstring(Filename);
  char *cLocation = RtoCstring(location);
  char *cName = RtoCstring(name);
  char *value = new char[256];

  TwRetVal rv = TwGetStringAttributeFromH5(cFilename, cLocation, cName, value);

  TwCloseH5(cFilename);
  if (rv != TwSuccess) {
    return TwRetValString(rv);
  }

  std::string str(value);

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
SEXP SetIntAttributeInH5(SEXP Filename, SEXP location, SEXP name, int attribute) {

  char *cFilename = RtoCstring(Filename);
  char *cLocation = RtoCstring(location);
  char *cName = RtoCstring(name);

  TwRetVal rv = TwSetIntAttributeInH5(cFilename, cLocation, cName, attribute);
  TwCloseH5(cFilename);

  return TwRetValString(rv);
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
SEXP SetUintAttributeInH5(SEXP Filename, SEXP location, SEXP name, double attribute) {

  char *cFilename = RtoCstring(Filename);
  char *cLocation = RtoCstring(location);
  char *cName = RtoCstring(name);

  TwRetVal rv = TwSetUintAttributeInH5(cFilename, cLocation, cName, (unsigned int)attribute);
  TwCloseH5(cFilename);

  return TwRetValString(rv);
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
SEXP SetInt64AttributeInH5(SEXP Filename, SEXP location, SEXP name, SEXP attribute) {

  char *cFilename = RtoCstring(Filename);
  char *cLocation = RtoCstring(location);
  char *cName = RtoCstring(name);

  std::string str = Rcpp::as<std::string>(attribute);
  std::stringstream ss(str);
  int64_t int64value;
  ss >> int64value;

  TwRetVal rv = TwSetInt64AttributeInH5(cFilename, cLocation, cName, int64value);
  TwCloseH5(cFilename);

  return TwRetValString(rv);
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
SEXP SetUint64AttributeInH5(SEXP Filename, SEXP location, SEXP name, SEXP attribute) {

  char *cFilename = RtoCstring(Filename);
  char *cLocation = RtoCstring(location);
  char *cName = RtoCstring(name);

  std::string str = Rcpp::as<std::string>(attribute);
  std::stringstream ss(str);
  uint64_t uint64value;
  ss >> uint64value;

  TwRetVal rv = TwSetUint64AttributeInH5(cFilename, cLocation, cName, uint64value);
  TwCloseH5(cFilename);

  return TwRetValString(rv);
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
SEXP SetFloatAttributeInH5(SEXP Filename, SEXP location, SEXP name, double attribute) {

  char *cFilename = RtoCstring(Filename);
  char *cLocation = RtoCstring(location);
  char *cName = RtoCstring(name);

  TwRetVal rv = TwSetFloatAttributeInH5(cFilename, cLocation, cName, (float)attribute);
  TwCloseH5(cFilename);

  return TwRetValString(rv);
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
SEXP SetDoubleAttributeInH5(SEXP Filename, SEXP location, SEXP name, double attribute) {

  char *cFilename = RtoCstring(Filename);
  char *cLocation = RtoCstring(location);
  char *cName = RtoCstring(name);

  TwRetVal rv = TwSetDoubleAttributeInH5(cFilename, cLocation, cName, attribute);
  TwCloseH5(cFilename);

  return TwRetValString(rv);
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
SEXP SetStringAttributeInH5(SEXP Filename, SEXP location, SEXP name, SEXP attribute) {

  char *cFilename = RtoCstring(Filename);
  char *cLocation = RtoCstring(location);
  char *cName = RtoCstring(name);
  char *cAttribute = RtoCstring(attribute);

  TwRetVal rv = TwSetStringAttributeInH5(cFilename, cLocation, cName, cAttribute);
  TwCloseH5(cFilename);

  return TwRetValString(rv);
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
//' @param rowIndex Row index.
//' @return A list containing the row of user data and data description.
//'
//' @examples
//' \dontrun{
//' GetUserDataFromH5("path/to/file.h5", location = , rowIndex = )
//' }
//' @export
// [[Rcpp::export]]
SEXP GetUserDataFromH5(SEXP Filename, SEXP location, int rowIndex) {

  char *cFilename = RtoCstring(Filename);
  char *cLocation = RtoCstring(location);

  // get nbrElements
  int nbrElements = 0;
  int tmp = 0;
  TwRetVal rv = TwGetUserDataFromH5(cFilename, cLocation, &tmp, &nbrElements,
                                    NULL, NULL);
  if (rv != TwValueAdjusted) {
    TwCloseH5(cFilename);
    return TwRetValString(rv);
  }

  NumericVector buffer(nbrElements);
  // CharacterVector elementDescription(nbrElements);
  char *elementDescription = new char[256 * nbrElements];
  memset(elementDescription, 0, 256 * nbrElements);

  rv = TwGetUserDataFromH5(cFilename, cLocation, &rowIndex, &nbrElements,
                           &buffer[0], elementDescription);
  TwCloseH5(cFilename);

  if (rv != TwSuccess) {
    delete[] elementDescription;
    return TwRetValString(rv);
  }

  CharacterVector descriptionArray(nbrElements);
  std::string str(elementDescription);
  delete[] elementDescription;

  for (int i = 0; i < nbrElements; ++i) {
    descriptionArray[i] = str.substr(i*256, 256);
  }

  List result;
  result["data"] = buffer;
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
SEXP GetAcquisitionLogFromH5(SEXP Filename, int index) {

  char *cFilename = RtoCstring(Filename);

  int64_t timestamp;
  char * logText = new char[256];

  TwRetVal rv = TwGetAcquisitionLogFromH5(cFilename, index, &timestamp, logText);
  TwCloseH5(cFilename);

  if (rv != TwSuccess) {
    return TwRetValString(rv);
  }

  CharacterVector time(1);
  time[0] = std::to_string(timestamp);

  std::string str(logText);

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
SEXP GetEventListSpectrumFromH5(SEXP Filename, int segmentIndex, int bufIndex,
                                int writeIndex) {

  char *cFilename = RtoCstring(Filename);

  // get bufferSize
  int bufferSize = 0;
  TwRetVal rv = TwGetEventListSpectrumFromH5(cFilename, segmentIndex, bufIndex,
                                             writeIndex, &bufferSize, NULL);
  if (rv != TwValueAdjusted) {
    TwCloseH5(cFilename);
    return TwRetValString(rv);
  }

  std::vector<unsigned int> buffer(bufferSize);

  rv = TwGetEventListSpectrumFromH5(cFilename, segmentIndex, bufIndex,
                                    writeIndex, &bufferSize, &buffer[0]);
  TwCloseH5(cFilename);
  if (rv != TwSuccess) {
    return TwRetValString(rv);
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
SEXP H5GetMassCalibPar(SEXP Filename, int writeIndex) {

  char *cFilename = RtoCstring(Filename);

  int segmentIndex = 0;
  int bufIndex = 0;

  // get nbrParams
  int nbrParams = 0;
  TwRetVal rv = TwH5GetMassCalibPar(cFilename, segmentIndex, bufIndex,
                                    writeIndex, NULL, &nbrParams, NULL);
  if (rv != TwValueAdjusted) {
    TwCloseH5(cFilename);
    return TwRetValString(rv);
  }

  int mode;
  NumericVector p(nbrParams);

  rv = TwH5GetMassCalibPar(cFilename, segmentIndex, bufIndex,
                           writeIndex, &mode, &nbrParams, &p[0]);
  TwCloseH5(cFilename);
  if (rv != TwSuccess) {
    return TwRetValString(rv);
  }

  List result;
  result["mode"] = mode;
  result["p"] = p;

  return result;
}

// H5AddLogEntry -----------------------------------------------------------------
//' Adds an entry to an existing data file.
//'
//' \code{H5AddLogEntry} adds an entry to an existing data file. To add
//' acquisition log entries during a running acquisition use \code{AddLogEntry}.
//'
//' @param Filename Path/filename of the HDF5 file.
//' @param LogEntryText Log text (max. 255 characters).
//' @param LogEntryTime Log entry time (number of 100-nanosecond intervals since
//' January 1, 1601 UTC) passed as a string. Set it to "0" for "now".
//'
//' @export
// [[Rcpp::export]]
SEXP H5AddLogEntry(SEXP Filename, SEXP LogEntryText, SEXP LogEntryTime) {

  char *cFilename = RtoCstring(Filename);
  char *cLogEntryText = RtoCstring(LogEntryText);

  std::string str = Rcpp::as<std::string>(LogEntryTime);
  std::stringstream ss(str);
  uint64_t cTime;
  ss >> cTime;

  TwRetVal rv = TwH5AddLogEntry(cFilename, cLogEntryText, cTime);

  return TwRetValString(rv);
}

// Not implemented -------------------------------------------------------------
// TwGetBufWriteProfileFromH5_2
// TwChangePeakTable
// TwChangePeakTable2
// TwChangePeakFromFile
// TwChangePeakFromFile2
// TwProgressCallback
// TwProgressCallback2
// TwChangePeakDataInit
// TwChangePeakDataWrite
// TwChangePeakDataFinalize
// TwReadRawData
// TwH5SetMassCalib
// TwH5SetMassCalibEx
// TwGetEventListDataFromH5
// TwGetEventListBlobFromH5
// TwFreeEventListData
// TwFreeEventListData2
// TwMultiPeakFitIntegration
// TwH5AddUserDataMultiRow
// TwH5SetMassCalibDynamic
// TwGenerateSegmentProfilesFromEventList
// TwH5MakePaletteImage
// TwH5MakeTrueColorImage
