% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/RcppExports.R
\name{GetMassCalib2Ex}
\alias{GetMassCalib2Ex}
\title{Returns information about the current mass calibration.}
\usage{
GetMassCalib2Ex()
}
\value{
List with calibration parameters, calibration points and labels.
}
\description{
\code{GetMassCalib2Ex} returns information about the mass calibration currently
used in TofDaq recorder.
}
\details{
This is the same as \code{\link{GetMassCalib2}}, but additionally also returns
the labels of the calibration points.
}
