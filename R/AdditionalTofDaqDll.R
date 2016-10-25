# SiProcessSpectrumFromShMem ---------------------------------------------------
#' Processes a spectrum taken from shared memory.
#'
#' \code{SiProcessSpectrumFromShMem} processes a spectrum taken from shared
#' memory according to the options set for it's spectrum type.
#'
#' This function is a variant of the original TwToolDll function \code{\link{SiProcessSpectrum}}.
#'
#' @param specType Spectrum type index (non-negative integer).
#' @param BufIndex Buf index of data to fetch.
#' @return A list with the baseline and threshold value.
#'
#' @family Single ion histogramming functions.
#' @export
SiProcessSpectrumFromShMem <- function(specType = 0, BufIndex) {
  .Call("SiProcessSpectrumFromShMem", specType, BufIndex)
}

# KeepSharedMemMapped ----------------------------------------------------------
#' Keeps the shared memory acquisition buffers mapped.
#'
#' \code{KeepSharedMemMapped} Keeps the shared memory acquisition buffers mapped.
#'
#' The DLL periodically unmaps the shared memory to give the recorder
#' application the possibility to (re)allocate the shared buffers. Set this
#' option to true if you want to make sure that the shared memory pointers stay
#' valid while you work with them. In this case you must call
#' \code{\link{ReleaseSharedMemory}} explicitly when finished with your
#' processing operation.
#'
#' @param keepMapped \code{TRUE} or \code{FALSE}.
#' @export
KeepSharedMemMapped <- function(keepMapped) {
  rv <- .Call("KeepSharedMemMapped", keepMapped)
  if (rv!=4) {return(warning(TwRetVal[rv+1]))} else {return(TwRetVal[5])}
}
