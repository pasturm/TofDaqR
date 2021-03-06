% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/RcppExports.R
\name{GetMassCalibEx}
\alias{GetMassCalibEx}
\title{Returns information about the current mass calibration.}
\usage{
GetMassCalibEx()
}
\value{
List with calibration parameters, calibration points and labels.
}
\description{
\code{GetMassCalibEx} returns information about the mass calibration currently
used in TofDaq recorder.
}
\details{
This is the same as \code{\link{GetMassCalib}}, but additionally also returns
the labels of the calibration points.
}
