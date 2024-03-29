% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/RcppExports.R
\name{TpsSaveSetFileRc}
\alias{TpsSaveSetFileRc}
\title{Saves TPS set values with a RC code to a file.}
\usage{
TpsSaveSetFileRc(setFile)
}
\arguments{
\item{setFile}{Path/filename of the set file to save.}
}
\description{
\code{TpsSaveSetFileRc} saves the current TPS settings to a file. Only
values with assigned RC codes will be saved.
}
\details{
Note: set files saved with this function can not be loaded through the TPS
web GUI, only \code{\link{TpsLoadSetFile}} and \code{\link{TpsLoadSetFile2}}
understand this format.
}
\seealso{
Other TPS functions: 
\code{\link{TpsChangeIonMode}()},
\code{\link{TpsConnect2}()},
\code{\link{TpsConnect}()},
\code{\link{TpsDisconnect}()},
\code{\link{TpsGetActiveFilament}()},
\code{\link{TpsGetLastSetValue}()},
\code{\link{TpsGetModuleCodes}()},
\code{\link{TpsGetModuleLimits}()},
\code{\link{TpsGetModuleProperties}()},
\code{\link{TpsGetMonitorValue}()},
\code{\link{TpsGetNbrModules}()},
\code{\link{TpsGetNmtState}()},
\code{\link{TpsGetStatus}()},
\code{\link{TpsGetTargetValue}()},
\code{\link{TpsInitialize}()},
\code{\link{TpsLoadSetFile2}()},
\code{\link{TpsLoadSetFile}()},
\code{\link{TpsSaveSetFile}()},
\code{\link{TpsSetActiveFilament}()},
\code{\link{TpsSetAllVoltages}()},
\code{\link{TpsSetNmtCmd}()},
\code{\link{TpsSetTargetValue}()},
\code{\link{TpsShutdown}()}
}
\concept{TPS functions}
