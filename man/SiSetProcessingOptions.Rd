% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/RcppExports.R
\name{SiSetProcessingOptions}
\alias{SiSetProcessingOptions}
\title{Sets processing options for each spectrum type.}
\usage{
SiSetProcessingOptions(option, value, specType)
}
\arguments{
\item{option}{Option to set (see below for a list of valid option strings).}

\item{value}{Value to set the for the given option.}

\item{specType}{Spectrum type index. -1 is a wildcard for all spectrum types.}
}
\description{
\code{SiSetProcessingOptions} sets processing options for each spectrum type.
}
\details{
Options:
\tabular{llc}{
Name \tab Description \tab Default value \cr
\code{MassCalibMode} \tab Mass calibration mode in use (see \code{\link{MassCalibrate}}). \tab 0 \cr
\code{MassCalibParamn} \tab Mass calibration parameters n =
1..number of calibration parameters for the given mass calibration mode. \tab \code{c(1000, 0)} \cr
\code{FullScale} \tab Full digitizing range of the ADC in the same units as
the spectra to be analyzed (typically mV). \tab 500 \cr
\code{NbrBits} \tab ADC resolution (8 for AP240, 14 for ADQ114 etc.). \tab 8 \cr
\code{SampleInterval} \tab Sample interval in ns. \tab 1 \cr
\code{PreampGain} \tab Gain of external preamp. \tab 1 \cr
\code{PreSamples} \tab Number of samples before a threshold crosser taken into account. \tab 0 \cr
\code{PostSamples} \tab Number of samples after a negative threshold crosser taken into account. \tab 0 \cr
\code{BaselineAndThresholdFromData} \tab If >0 the baseline and threshold
values will be determined based on \code{NbrStdDevs} for every processed
spectrum. If >1.5 baseline noise is determined from a fit to a histogram of
all data (instead of from the standard deviation of all data). This makes the
noise determination more robust when real peaks are present in the spectrum. \tab 0 \cr
\code{NbrStdDevs} \tab Number of standard deviations of baseline noise that
defines the threshold. Only relevant if \code{BaselineAndThresholdFromData>0}. \tab 6 \cr
\code{Baseline} \tab Baseline value used for calculation of intensities. Has
no meaning if \code{BaselineAndThresholdFromData>0}. \tab 5 \cr
\code{Threshold} \tab Threshold value used for calculation of intensities.
Has no meaning if \code{BaselineAndThresholdFromData!=0}. \tab 8 \cr
\code{NegativeSignal} \tab Indicates peak polarity with respect to the baseline. \tab 0 (=FALSE) \cr
\code{BaselineAndThresholdInCodes} \tab Indicates whether the values in \code{Baseline}
and \code{Threshold} are interpreted as ADC codes or mV. \tab 1 (=TRUE)
}
}
