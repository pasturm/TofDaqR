% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/RcppExports.R
\name{GetUint64AttributeFromH5}
\alias{GetUint64AttributeFromH5}
\title{Reads an unsigned 64-bit integer attribute from the HDF5 file.}
\usage{
GetUint64AttributeFromH5(Filename, location, name)
}
\arguments{
\item{Filename}{Path/filename of the HDF5 file.}

\item{location}{Location of the group or dataset where the attribute is attached to.}

\item{name}{Attribute name.}
}
\value{
An int64 attribute.
}
\description{
\code{GetUint64AttributeFromH5} reads an unsigned 64-bit integer attribute
from the HDF5 file.
}
\details{
Used to read attributes not available from \code{GetH5Descriptor}. Unsigned
int64 parameters are returned as string.
}
\examples{
\dontrun{
GetUint64AttributeFromH5("path/to/file.h5", location = , name = )
}
}
