#' TofDaqR: R Interface to TOFWERK TofDaq API.
#'
#' The TofDaqR package is a wrapper for the TofDaq API functions provided by the
#' shared libraries TofDaqDll, TwH5Dll and TwToolDll (Windows) and libtwh5 and
#' libtwtool (macOS and Linux).
#'
#' @section TofDaqR functions: TofDaqR functions can be used to configure the TofDaq recorder
#' application, control running acquisitions, access the data in realtime and
#' store arbitrary data synchronously or asynchronously to the recorded mass
#' spectra. Further, there are functions for controling the TOF Power
#' Supply and to access and modify recorded data stored in HDF5 files. The
#' package also contains additional "tool" functions, e.g. for peak fitting.
#'
#' @docType package
#' @name TofDaqR
#' @useDynLib TofDaqR
NULL
