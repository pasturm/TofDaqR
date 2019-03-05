#' TofDaqR: R Interface to the TOFWERK TofDaq API.
#'
#' The TofDaqR package provides a R interface to the TOFWERK TofDaq API, which
#' consists of libraries for communication with the TofDAQ recorder application,
#' data file access and general (time-of-flight) mass spectrometry related
#' utility functions.
#'
#' Features:
#' \itemize{
#'   \item Acquisition setup and control
#'   \item Real-time data access
#'   \item Adding (structured) custom data to the data file alongside the TOF data
#'   \item Control of the TOF Power Supply
#'   \item Quick access to TOFWERK HDF5 data files without the need to study the file format details
#'   \item Add additional data to existing data files (e.g. post-processing results)
#'   \item Peak fitting functions
#'   \item Mass calibration
#'   \item Chemistry functions (molecular mass and isotope pattern calculation)
#'   \item Single ion analysis functions
#' }
#'
#' @seealso \link[=API-documentation]{TofDaq API documentation}.
#'
#' @docType package
#' @author Patrick Sturm <sturm@tofwerk.com>
#' @importFrom Rcpp evalCpp
#' @name TofDaqR
#' @useDynLib TofDaqR, .registration = TRUE
NULL
