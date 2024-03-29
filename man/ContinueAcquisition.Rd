% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/RcppExports.R
\name{ContinueAcquisition}
\alias{ContinueAcquisition}
\title{Signals to the TofDaq recorder to continue an acquisition.}
\usage{
ContinueAcquisition()
}
\description{
\code{ContinueAcquisition} signals to the TofDaq recorder to continue an
acquisition.
}
\details{
This is a legacy function that was used with some Acqiris DAQ cards,
where every block was armed by software. The feature is enabled by setting
the parameter ManualContinueEveryNMemories to a value > 0. All latest DAQ
devices operate in a streaming mode in order to achieve 100 \% duty cycle and
TofDaq recorder no longer has per block control of the DAQ progress.
}
