% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/RcppExports.R
\name{GetStickSpectrumFromShMem2}
\alias{GetStickSpectrumFromShMem2}
\title{Single stick spectrum from shared memory.}
\usage{
GetStickSpectrumFromShMem2(SegmentIndex, SegmentEndIndex, BufIndex)
}
\arguments{
\item{SegmentIndex}{Segment start index of data to fetch.}

\item{SegmentEndIndex}{Segment end index of data to fetch.}

\item{BufIndex}{Buf index of data to fetch.}
}
\value{
A list containing the stick spectrum and corresponding masses.
}
\description{
\code{GetStickSpectrumFromShMem2} reads a single stick spectrum from shared
memory.
}
