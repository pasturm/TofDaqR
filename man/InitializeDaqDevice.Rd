% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/RcppExports.R
\name{InitializeDaqDevice}
\alias{InitializeDaqDevice}
\title{Initializes the DAQ board.}
\usage{
InitializeDaqDevice()
}
\description{
\code{InitializeDaqDevice} initializes the DAQ board (this is also done at
startup of TofDaqRec.exe). This can take up to 8 seconds depending on the
actual DAQ hardware.
}
