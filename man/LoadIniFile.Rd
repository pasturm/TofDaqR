% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/RcppExports.R
\name{LoadIniFile}
\alias{LoadIniFile}
\title{Loads a configuration file.}
\usage{
LoadIniFile(IniFile = NULL)
}
\arguments{
\item{IniFile}{Path/filename of the configuration file. If no path is
specified, the TofDaq recorder directory will be used. If \code{IniFile} is
an empty string or \code{NULL}, "TwApiTmpIni.ini" will be loaded.}
}
\description{
\code{LoadIniFile} loads a configuration file (*.ini) from disk.
}
