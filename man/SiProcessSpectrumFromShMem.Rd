% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/RcppExports.R
\name{SiProcessSpectrumFromShMem}
\alias{SiProcessSpectrumFromShMem}
\title{Processes a spectrum taken from shared memory.}
\usage{
SiProcessSpectrumFromShMem(specType, BufIndex)
}
\arguments{
\item{specType}{Spectrum type index (non-negative integer).}

\item{BufIndex}{Buf index of data to fetch.}
}
\value{
A list with the baseline and threshold value.
}
\description{
\code{SiProcessSpectrumFromShMem} processes a spectrum taken from shared
memory according to the options set for it's spectrum type.
}
\details{
This function is a variant of the original TwToolDll function \code{\link{SiProcessSpectrum}}.
}
