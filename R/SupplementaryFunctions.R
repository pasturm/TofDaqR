# Tps1rc ----------------------------------------------------------------
#' TPS1 RC codes and names.
#'
#' Data frame of TPS1 codes (as used in tps1rc.cfg to communicate with the TPS
#' through the API) and corresponding names.
#' @export
Tps1rc = c(
  1, "RBP",
  2, "RG",
  3, "LENS",
  4, "HM",
  5, "PA",
  6, "MCP",
  9, "DRIFT",
  10, "U-low",
  11, "U-high",
  12, "U+low",
  13, "U+high",
  14, "L2",
  15, "DEFL",
  16, "DEFLFL",
  17, "IONEX",
  18, "L1",
  19, "EP",
  20, "IONCH",
  21, "FIL",
  22, "FILEM",
  23, "IFIL",
  24, "ILIM",
  27, "FILNUM",
  28, "SKIMMER",
  49, "Q1F",
  50, "Q1B",
  51, "Q2F",
  52, "Q2B",
  53, "LSK",
  55, "IMSLENS",
  56, "NOZZLE",
  57, "Q1EP",
  58, "IMSEND",
  62, "RF1",
  63, "RFAMP1",
  64, "RF2",
  65, "RFAMP2",
  66, "LFFREQ1",
  67, "LFAMP1",
  68, "LFFREQ2",
  69, "LFAMP2",
  70, "LFFREQ3",
  71, "LFAMP3",
  72, "LFFREQ4",
  73, "LFAMP4",
  83, "Grid+",
  84, "Grid-",
  85, "NEEDLE",
  86, "IMSHV",
  87, "IMSGRIDSTATE",
  88, "IMS_TR_FREQ",
  90, "IMSGateHV",
  116, "SKIM2",
  117, "REFERENCE",
  121, "CORONA",
  124, "IMR",
  162, "QUADA_COIL_SET",
  163, "QUADA_RFDAC_SET",
  164, "QUADB_COIL_SET",
  165, "QUADB_RFDAC_SET",
  200, "TOF_EXTR2",
  201, "TOF_EXTR1",
  202, "TOF_REF",
  203, "TOF_PULSE",
  600, "IONMODE",
  601, "INTERLOCK",
  602, "HVSUPPLY",
  603, "HVPOS",
  604, "HVNEG",
  605, "SPARE1",
  606, "SPARE2",
  607, "SPARE3",
  608, "MODE",
  609, "EXT24V",
  1016,	"TPSBODYTEMP",
  1017,	"TPSAIRTEMP",
  1018,	"HVATEMP",
  1019, "HVATEMP"
)
Tps1rc = as.data.frame(matrix(Tps1rc, ncol = 2, byrow = TRUE), stringsAsFactors = FALSE)
names(Tps1rc) = c("codes", "names")
Tps1rc$codes = as.integer(Tps1rc$codes)

# SaveMassTableToFile ----------------------------------------------------------
#' Saves the current peak parameters to a file.
#'
#' \code{SaveMassTableToFile} saves the current peak parameters to a file.
#'
#' Note: this function is not part of the TofDaq API, but is included in the
#' package for convenience.
#'
#' @param filename Path/filename. If not specified, "TmpMassTable.txt" in the
#' current working directory will be used.
#' @export
SaveMassTableToFile = function(filename = "TmpMassTable.txt") {
  if (file.exists(filename)) {
    file.remove(filename)
    file.create(filename)
  }
  desc = GetDescriptor()
  utils::write.table(desc$NbrPeaks, file = filename, quote = FALSE, row.names = FALSE,
              col.names = FALSE)
  for (i in 1:desc$NbrPeaks) {
    rv = GetPeakParameters(i-1)
    utils::write.table(t(unlist(rv)), file = filename, append = TRUE,
                quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)
  }
}

# GetTofDataSinglePeak ---------------------------------------------------------
#' Gets TofData of a single peak from a HDF5 data file.
#'
#' \code{GetTofDataSinglePeak} gets TofData of a single peak or a specified mass
#' range from a HDF5 data file.
#'
#' Note: this function is not part of the TofDaq API, but is included in the
#' package for convenience.
#'
#' @param filename Path/filename.
#' @param peakIndex Peak index (one-based).
#' @param massRange Mass range in Th. Specified as a two-element vector
#' containing the low and high limit of the mass range. If \code{NULL} the mass
#' range is taken from the peak parameter list.
#' @param secondTOF If \code{TRUE} the data are read from /FullSpectra2/TofData
#' @return A list containing the mass axis, time axis, TofData and averaged
#' spectrum of the peak, and the descriptor of the data file.
#'
#' @family wrapper functions
#'
#' @export
GetTofDataSinglePeak = function(filename, peakIndex, massRange = NULL,
                                secondTOF = FALSE) {

  # read mass axis data
  if (secondTOF) {
    MassAxis = GetSpecXaxisFromH5(filename, Type = -1, writeIndex = 0)
    TimeAxis = GetSpecXaxisFromH5(filename, Type = -2, writeIndex = 0)
  } else {
    MassAxis = GetSpecXaxisFromH5(filename, Type = 1, writeIndex = 0)
    TimeAxis = GetSpecXaxisFromH5(filename, Type = 2, writeIndex = 0)
  }

  if (is.null(massRange)) {
    peakpar = GetPeakParametersFromH5(filename, PeakIndex = peakIndex - 1)
  } else {
    peakpar = data.frame(loMass = massRange[1], hiMass = massRange[2])
  }

  idx = which(MassAxis>(peakpar$loMass)&MassAxis<(peakpar$hiMass))  # sample indices

  # Read H5descriptor
  desc = GetH5Descriptor(filename)

  # read FullSpectra data of one mass
  if (secondTOF) {
    TofData = GetTofData2(filename,  idx[1]-1, length(idx), 0, desc$nbrSegments, 0, desc$nbrBufs, 0, desc$nbrWrites)
  } else {
    TofData = GetTofData(filename,  idx[1]-1, length(idx), 0, desc$nbrSegments, 0, desc$nbrBufs, 0, desc$nbrWrites)
  }
  # convert units to mV
  TofData = TofData/desc$nbrWaveforms

  # average over buf and write dimension
  TofData = matrix(TofData, ncol = desc$nbrBufs*desc$nbrWrites)
  AverageSpectrum = rowSums(TofData)/(desc$nbrBufs*desc$nbrWrites)

  return(list(MassAxis = MassAxis[idx], TimeAxis = TimeAxis[idx],
              TofData = TofData, AverageSpectrum = AverageSpectrum,
              desc = desc))
}

# FitTofDataSinglePeak ---------------------------------------------------------
#' Wrapper function for FitSinglePeak().
#'
#' \code{FitTofDataSinglePeak} is a wrapper for \code{\link{FitSinglePeak}}
#' and also reports fitted resolving power, area and count rates. It takes the
#' output of \code{\link{GetTofDataSinglePeak}} as input parameter.
#'
#' Peak types:
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
#' }
#'
#' Note: this function is not part of the TofDaq API, but is included in the
#' package for convenience.
#'
#' @param PeakTofData List from the output of \code{\link{GetTofDataSinglePeak}}.
#' @param peakType Peak model to use.
#' @return A list with fitted parameters.
#'
#' @family wrapper functions
#'
#' @export
FitTofDataSinglePeak = function(PeakTofData, peakType) {

  # initial values
  mu = 0.5 # mixing parameter
  amplitude = max(PeakTofData$AverageSpectrum) # amplitude

  # TimeAxis
  peakPos = PeakTofData$TimeAxis[PeakTofData$AverageSpectrum==amplitude] # position in tof (ns)
  fwhm = diff(range(PeakTofData$TimeAxis[PeakTofData$AverageSpectrum>amplitude/2], na.rm=TRUE)) # FWHM in tof (ns)
  fwhm = max(fwhm, min(diff(PeakTofData$TimeAxis))) # make it > 0

  fit1 = FitSinglePeak(yVals = PeakTofData$AverageSpectrum, xVals = PeakTofData$TimeAxis, peakType = peakType,
                       amplitude = amplitude, fwhmLo = fwhm, fwhmHi = fwhm,
                       peakPos = peakPos, mu = mu)

  # resolution
  fit1$resolution = fit1$peakPos/(fit1$fwhmLo+fit1$fwhmHi)

  # integrate mV over time axis (microsec) to get mV*ns
  if (any(peakType == c(0, 3, 6, 9))) {
    fit1$area = fit1$amplitude*(fit1$fwhmLo+fit1$fwhmHi)/2*sqrt(pi/log(2))/2*1000  # Gaussian
  } else if (any(peakType == c(1, 4, 7, 10))) {
    fit1$area = fit1$amplitude*(fit1$fwhmLo+fit1$fwhmHi)/2*pi/2*1000  # Lorentzian
  } else if (any(peakType == c(2, 5, 8, 11))) {
    fit1$area = fit1$amplitude*(fit1$fwhmLo+fit1$fwhmHi)/2/2*(fit1$mu*pi+(1-fit1$mu)*sqrt(pi/log(2)))*1000  # Pseudo-Voigt
  }

  # normalize with single ion signal to get counts per extraction
  fit1$cpe = fit1$area/PeakTofData$desc$singleIonSignal
  # counts per second
  fit1$cps = fit1$cpe/PeakTofData$desc$tofPeriod

  # MassAxis
  peakPos = PeakTofData$MassAxis[PeakTofData$AverageSpectrum==amplitude] # position in m/Q (Th)
  fwhm = diff(range(PeakTofData$MassAxis[PeakTofData$AverageSpectrum>amplitude/2], na.rm=TRUE)) # FWHM in m/Q (Th)
  fwhm = max(fwhm, min(diff(PeakTofData$MassAxis))) # make it > 0

  fit2 = FitSinglePeak(yVals = PeakTofData$AverageSpectrum, xVals = PeakTofData$MassAxis, peakType = peakType,
                       amplitude = amplitude, fwhmLo = fwhm, fwhmHi = fwhm,
                       peakPos = peakPos, mu = mu)

  # resolution
  fit2$resolution = fit2$peakPos/(fit2$fwhmLo+fit2$fwhmHi)*2

  return(list(fit_time = fit1, fit_mass = fit2))
}

