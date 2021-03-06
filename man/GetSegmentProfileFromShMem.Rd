% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/RcppExports.R
\name{GetSegmentProfileFromShMem}
\alias{GetSegmentProfileFromShMem}
\title{Segment profile for a given peak and buf index from shared memory.}
\usage{
GetSegmentProfileFromShMem(PeakIndex, BufIndex)
}
\arguments{
\item{PeakIndex}{Index of peak to fetch segment profile from. All peaks are
read if \code{PeakIndex = -1}.}

\item{BufIndex}{Buf index of data to fetch.}
}
\value{
A vector containing the segment profile(s).
}
\description{
\code{GetSegmentProfileFromShMem} reads the segment profile for a given
peak and buf index from shared memory. Use -1 for \code{PeakIndex} to get
segment profiles of all peaks.
}
