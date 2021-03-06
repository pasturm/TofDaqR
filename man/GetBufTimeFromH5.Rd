% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/RcppExports.R
\name{GetBufTimeFromH5}
\alias{GetBufTimeFromH5}
\title{Single buf timestamp from the data file.}
\usage{
GetBufTimeFromH5(Filename, BufIndex, WriteIndex)
}
\arguments{
\item{Filename}{Path/filename of the HDF5 file.}

\item{BufIndex}{Buf index.}

\item{WriteIndex}{Write index.}
}
\value{
A time stamp (in seconds relative to acquisition start).
}
\description{
\code{GetBufTimeFromH5} reads a single buf time stamp from HDF5 file.
}
\examples{
\dontrun{
GetBufTimeFromH5("path/to/file.h5", BufIndex = 0, WriteIndex = 0)
}
}
