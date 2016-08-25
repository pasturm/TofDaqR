#' TofDaqR: R wrapper of TOFWERK's TofDaq API functions.
#'
#' The TofDaqR package is a wrapper for the TofDaq API functions provided by
#' TofDaqDll.dll, TwH5Dll.dll and TwToolDll.dll.
#'
#' @section TofDaqR functions: TofDaqR functions can be used to configure the
#'   TofDaq recorder application, control running acquisitions, access the data
#'   being recorded and store arbitrary data synchronously or asynchronously to
#'   the recorded mass spectra. \code{TPS...} functions allow to control the TOF Power
#'   Supply. Further, there are functions to access and modify recorded data
#'   stored in Tofwerk HDF5 files and additional "tool" functions, e.g. for peak
#'   fitting.
#'
#' @docType package
#' @name TofDaqR
#'
#' @useDynLib TofDaqR
#' @importFrom Rcpp sourceCpp
NULL
