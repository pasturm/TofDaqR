# TwRetVal ---------------------------------------------------------------------
#' Return value.
#'
#' Translation of the integer \code{TwRetVal} return value into a string.
#'
#' @keywords internal
TwRetVal <- c("TwDaqRecNotRunning",
              "TwAcquisitionActive",
              "TwNoActiveAcquisition",
              "TwFileNotFound",
              "TwSuccess",
              "TwError",
              "TwOutOfBounds",
              "TwNoData",
              "TwTimeout",
              "TwValueAdjusted",
              "TwInvalidParameter",
              "TwInvalidValue",
              "TwAborted")

# GetH5Descriptor ------------------------------------------------------------
#' Descriptor structure of Tofwerk HDF5 data file.
#'
#' \code{GetH5Descriptor} returns a descriptor structure for the Tofwerk HDF5 file.
#'
#' The \emph{TwH5Desc} structure contains information about data dimensions,
#' available datasets and mass calibration. Additional attributes, which are not
#' available in the structure can be read using \code{Get...AttributeFromH5}
#' functions.
#' See \emph{/doc/TwH5Dll.htm} for more details.
#'
#' @param filename Path/filename of the HDF5 file.
#' @return A list containing the \emph{TwH5Desc} structure.
#'
#' @examples
#' \dontrun{
#' GetH5Descriptor("path/to/file.h5")
#' }
GetH5Descriptor <- function(filename) {
  .Call("GetH5Descriptor", filename)
}

# GetSumSpectrumFromH5 -------------------------------------------------------
#' Sum spectrum from HDF5 data file.
#'
#' \code{GetSumSpectrumFromH5} reads the sum spectrum (or average spectrum
#' depending on \code{Normalize} flag) from the HDF5 file.
#'
#' @param filename Path/filename of the HDF5 file.
#' @param normalize If \code{FALSE} (default) the spectrum is reported as sum,
#' if \code{TRUE} the spectrum is normalized to counts per extraction.
#' @return A vector containing the sum spectrum.
#'
#' @examples
#' \dontrun{
#' GetSumSpectrumFromH5("path/to/file.h5")
#' }
GetSumSpectrumFromH5 <- function(filename, normalize = FALSE) {
  .Call("GetSumSpectrumFromH5", filename, normalize)
}

# GetTofSpectrumFromH5 -------------------------------------------------------
#' Single (averaged) TOF spectrum from HDF5 data file.
#'
#' \code{GetTofSpectrumFromH5} reads a single mass spectrum (or an
#' averaged/summed hyperslab) from the HDF5 file.
#'
#' If \code{segment.start == segment.end} and \code{buf.start == buf.end} and
#' \code{write.start == write.end} and \code{normalize == FALSE} no
#' averaging/summing of spectra is done and the spectrum is reported as stored
#' in the dataset.
#'
#' @param filename Path/filename of the HDF5 file.
#' @param segment.start Segment start index of data to fetch.
#' @param segment.end Segment end index of data to fetch.
#' @param buf.start Buf start index of data to fetch.
#' @param buf.end Buf end index of data to fetch.
#' @param write.start Write start index of data to fetch.
#' @param write.end Write end index of data to fetch.
#' @param buf.write.linked Indicating whether the buf and write dimension should
#' be considered linked or treated as independent dimensions (relevant only if
#' \code{write.start != write.end}). Default is \code{FALSE}.
#' @param normalize If \code{FALSE} the spectrum is reported as sum,
#' if \code{TRUE} (default) the spectrum is normalized to counts per extraction.
#' @return A vector containing the mass spectrum.
#'
#' @examples
#' \dontrun{
#' GetTofSpectrumFromH5("path/to/file.h5")
#' }
GetTofSpectrumFromH5 <- function(filename, segment.start = 0, segment.end = 0,
                                   buf.start = 0, buf.end = 0, write.start = 0,
                                   write.end = 0, buf.write.linked = FALSE,
                                   normalize = TRUE) {
  .Call("GetTofSpectrumFromH5", filename, segment.start, segment.end, buf.start,
        buf.end, write.start, write.end, buf.write.linked, normalize)
}

# GetStickSpectrumFromH5 -----------------------------------------------------
#' Single (averaged) stick spectrum from HDF5 data file.
#'
#' \code{GetStickSpectrumFromH5} reads a single stick spectrum (or an
#' averaged/summed hyperslab) from the HDF5 file.
#'
#' If \code{segment.start == segment.end} and \code{buf.start == buf.end} and
#' \code{write.start == write.end} and \code{normalize == FALSE} no
#' averaging/summing of spectra is done and the spectrum is reported as stored
#' in the dataset.
#'
#' @param filename Path/filename of the HDF5 file.
#' @param segment.start Segment start index of data to fetch.
#' @param segment.end Segment end index of data to fetch.
#' @param buf.start Buf start index of data to fetch.
#' @param buf.end Buf end index of data to fetch.
#' @param write.start Write start index of data to fetch.
#' @param write.end Write end index of data to fetch.
#' @param buf.write.linked Indicating whether the buf and write dimension should
#' be considered linked or treated as independent dimensions (relevant only if
#' \code{write.start != write.end}). Default is \code{FALSE}.
#' @param normalize If \code{FALSE} the spectrum is reported as sum,
#' if \code{TRUE} (default) the spectrum is normalized to counts per extraction.
#' @return A vector containing the stick spectrum.
#'
#' @examples
#' \dontrun{
#' GetStickSpectrumFromH5("path/to/file.h5")
#' }
GetStickSpectrumFromH5 <- function(filename, segment.start = 0,
                                     segment.end = 0,  buf.start = 0,
                                     buf.end =0, write.start = 0, write.end = 0,
                                     buf.write.linked = FALSE, normalize = TRUE) {
  .Call("GetStickSpectrumFromH5", filename, segment.start, segment.end,
        buf.start, buf.end, write.start, write.end, buf.write.linked, normalize)
}

# GetPeakParametersFromH5 ----------------------------------------------------
#' Peak parameters from HDF5 data file.
#'
#' \code{GetPeakParametersFromH5} reads peak parameters from the data file.
#'
#' @param filename Path/filename of the HDF5 file.
#' @param index Index of peak. If index is -1 (default), peak parameters of all
#' peaks are read.
#' @return A list with the peak paramters \emph{label}, \emph{mass}, \emph{loMass} and \emph{hiMass}.
#'
#' @examples
#' \dontrun{
#' GetPeakParametersFromH5("path/to/file.h5")
#' }
GetPeakParametersFromH5 <- function(filename, index = -1) {
  .Call("GetPeakParametersFromH5", filename, index)
}

# GetBufTimeFromH5 -----------------------------------------------------------
#' Single buf timestamp from the data file.
#'
#' \code{GetBufTimeFromH5} reads a single buf time stamp from HDF5 file.
#'
#' @param filename Path/filename of the HDF5 file.
#' @param bufIndex Buf index.
#' @param writeIndex Write index.
#' @return A time stamp (in seconds relative to acquisition start).
#'
#' @examples
#' \dontrun{
#' GetBufTimeFromH5("path/to/file.h5", bufIndex = 0, writeIndex = 0)
#' }
GetBufTimeFromH5 <- function(filename, bufIndex, writeIndex) {
  .Call("GetBufTimeFromH5", filename, bufIndex, writeIndex)
}

# GetSpecXaxisFromH5 ---------------------------------------------------------
#' x-axis values of mass spectrum.
#'
#' \code{GetSpecXaxisFromH5} returns an array of x-axis values of the mass
#' spectrum.
#'
#' @param filename Path/filename of the HDF5 file.
#' @param type x-axis type (0: sample index, 1: mass/charge [Th] (default),
#' -1: mass/charge [Th] (2nd TOF), 2: time of flight [microsec],
#' -2: time of flight [microsec] (2nd TOF), 3: frequency [kHz]).
#' @param writeIndex Write index to use for mass calibration (relevant only for
#' \code{abs(type)== 1 or 2}). If the data file has no \emph{/TofData/MassCalibration}
#' dataset the standard mass calibration parameters are used (same for all
#' values of writeIndex).
#' @return A vector containing the x-axis values.
#'
#' @examples
#' \dontrun{
#' GetSpecXaxisFromH5("path/to/file.h5")
#' }
# x-axis values, type 0: sample index, 1: mass/charge [Th], 2: time of flight [microsec], 3: frequency [kHz]
GetSpecXaxisFromH5 <- function(filename, type = 1, writeIndex = 0) {
  .Call("GetSpecXaxisFromH5", filename, type, writeIndex)
}

# GetSegmentProfileFromH5 ----------------------------------------------------
#' Segment profile from HDF5 data file.
#'
#' \code{GetSegmentProfileFromH5} reads a segment profile for a given peak (or
#' all peaks) and buf and write slice.
#'
#' @param filename Path/filename of the HDF5 file.
#' @param PeakIndex Index of peak to fetch segment profile from. All peaks are
#' read if \code{PeakIndex = -1}.
#' @param buf.start Buf start index of data to fetch.
#' @param buf.end Buf end index of data to fetch.
#' @param write.start Write start index of data to fetch.
#' @param write.end Write end index of data to fetch.
#' @param buf.write.linked Indicating whether the buf and write dimension should
#' be considered linked or treated as independent dimensions (relevant only if
#' \code{write.start != write.end}). Default is \code{FALSE}.
#' @return A vector containing the segment profile(s).
#'
#' @examples
#' \dontrun{
#' GetSegmentProfileFromH5("path/to/file.h5", PeakIndex = -1, buf.start = 0,
#' buf.end = 0, write.start = 0, write.end = 0)
#' }
GetSegmentProfileFromH5 <- function(filename, PeakIndex, buf.start, buf.end,
                                      write.start, write.end, buf.write.linked = FALSE) {
  .Call("GetSegmentProfileFromH5", filename, PeakIndex, buf.start, buf.end,
        write.start, write.end, buf.write.linked)
}

# GetRegUserDataSourcesFromH5 --------------------------------------------------
#' Lists all registered user datasets available in the data file.
#'
#' \code{GetRegUserDataSourcesFromH5} lists all registered user data sets
#' available in the data file. Registered data sources can originate from data
#' source plugins, TofDaq recorder (e.g. DAQ temperatures) or data registered
#' through \code{RegisterUserData...} functions.
#'
#' @param filename Path/filename of the HDF5 file.
#' @return A list containing the location, the length, whether is has a description
#' and the type of the data source dataset. type 1: data source values are
#' written to disk for every write, type 2: data source values are written to
#' disk for every buf.
#'
#' @examples
#' \dontrun{
#' GetRegUserDataSourcesFromH5("path/to/file.h5")
#' }
GetRegUserDataSourcesFromH5 <- function(filename) {
  .Call("GetRegUserDataSourcesFromH5", filename)
}

# GetRegUserDataFromH5 ---------------------------------------------------------
#' Reads entries from a registered data source dataset.
#'
#' \code{GetRegUserDataFromH5} reads an entry (or all entries) from a
#' registered data source dataset (created by \code{RegisterUserData...} functions).
#'
#' If \code{bufIndex = writeIndex = -1}, all data are read.
#'
#' @param filename Path/filename of the HDF5 file.
#' @param location Location of the group or dataset where the attribute is attached to.
#' @param bufIndex Buf index.
#' @param writeIndex Write index.
#' @param readDescription If \code{TRUE} (default) the data descripton is read,
#' if \code{false} the data description is not read.
#' @return A list containing the registered user data and description.
#'
#' @examples
#' \dontrun{
#' GetRegUserDataFromH5("path/to/file.h5", location = )
#' }
GetRegUserDataFromH5 <- function(filename, location, bufIndex, writeIndex, readDescription = TRUE) {
  .Call("GetRegUserDataFromH5", filename, location, bufIndex, writeIndex, readDescription)
}

# GetTofData -----------------------------------------------------------------
#' Gets data stored in /FullSpectra/TofData from HDF5 data file.
#'
#' \code{GetTofData} gets data stored in \code{/FullSpectra/TofData} from HDF5 data file.
#'
#' @param filename Path/filename of the HDF5 file.
#' @param sampleOffset Sample offset.
#' @param sampleCount Sample count.
#' @param segOffset Segment offset.
#' @param segCount Segment count.
#' @param bufOffset Buf offset.
#' @param bufCount Buf count.
#' @param writeOffset Write offset.
#' @param writeCount Write count.
GetTofData <- function(filename, sampleOffset, sampleCount, segOffset, segCount,
                         bufOffset, bufCount, writeOffset, writeCount) {
  .Call("GetTofData", filename, sampleOffset, sampleCount, segOffset, segCount,
        bufOffset, bufCount, writeOffset, writeCount)
}

# GetTofData2 ----------------------------------------------------------------
#' Gets data stored in /FullSpectra2/TofData from HDF5 data file.
#'
#' \code{GetTofData2} gets data stored in \code{/FullSpectra2/TofData} from HDF5 data file.
#'
#' @param filename Path/filename of the HDF5 file.
#' @param sampleOffset Sample offset.
#' @param sampleCount Sample count.
#' @param segOffset Segment offset.
#' @param segCount Segment count.
#' @param bufOffset Buf offset.
#' @param bufCount Buf count.
#' @param writeOffset Write offset.
#' @param writeCount Write count.
GetTofData2 <- function(filename, sampleOffset, sampleCount, segOffset, segCount,
                         bufOffset, bufCount, writeOffset, writeCount) {
  .Call("GetTofData2", filename, sampleOffset, sampleCount, segOffset, segCount,
        bufOffset, bufCount, writeOffset, writeCount)
}

# GetPeakData ----------------------------------------------------------------
#' Gets data stored in /PeakData/PeakData from HDF5 data file.
#'
#' \code{GetPeakData} gets data stored in \code{/PeakData/PeakData} from HDF5 data file.
#'
#' @param filename Path/filename of the HDF5 file.
#' @param peakOffset Peak offset.
#' @param peakCount Peak count.
#' @param segOffset Segment offset.
#' @param segCount Segment count.
#' @param bufOffset Buf offset.
#' @param bufCount Buf count.
#' @param writeOffset Write offset.
#' @param writeCount Write count.
GetPeakData <- function(filename, peakOffset, peakCount, segOffset, segCount,
                           bufOffset, bufCount, writeOffset, writeCount) {
  .Call("GetPeakData", filename, peakOffset, peakCount, segOffset, segCount,
        bufOffset, bufCount, writeOffset, writeCount)
}

# GetPeakData2 ---------------------------------------------------------------
#' Gets data stored in /PeakData2/PeakData from HDF5 data file.
#'
#' \code{GetPeakData2} gets data stored in \code{/PeakData2/PeakData} from HDF5 data file.
#'
#' @param filename Path/filename of the HDF5 file.
#' @param peakOffset Peak offset.
#' @param peakCount Peak count.
#' @param segOffset Segment offset.
#' @param segCount Segment count.
#' @param bufOffset Buf offset.
#' @param bufCount Buf count.
#' @param writeOffset Write offset.
#' @param writeCount Write count.
GetPeakData2 <- function(filename, peakOffset, peakCount, segOffset, segCount,
                          bufOffset, bufCount, writeOffset, writeCount) {
  .Call("GetPeakData2", filename, peakOffset, peakCount, segOffset, segCount,
        bufOffset, bufCount, writeOffset, writeCount)
}

# GetTimingData --------------------------------------------------------------
#' Gets data stored in /Timing/BufTimes from HDF5 data file.
#'
#' \code{GetTimingData} gets data stored in \code{/Timing/BufTimes} from HDF5 data file.
#'
#' @param filename Path/filename of the HDF5 file.
#' @param bufOffset Buf offset.
#' @param bufCount Buf count.
#' @param writeOffset Write offset.
#' @param writeCount Write count.
GetTimingData <- function(filename, bufOffset, bufCount, writeOffset, writeCount) {
  .Call("GetTimingData", filename, bufOffset, bufCount, writeOffset, writeCount)
}

# GetIntAttributeFromH5 ------------------------------------------------------
#' Reads an integer attribute from the HDF5 file.
#'
#' \code{GetIntAttributeFromH5} reads an integer attribute from the HDF5 file.
#'
#' Used to read attributes not available from \code{GetH5Descriptor}.
#'
#' @param filename Path/filename of the HDF5 file.
#' @param location Location of the group or dataset where the attribute is attached to.
#' @param name Attribute name.
#' @return An integer attribute.
#'
#' @examples
#' \dontrun{
#' GetIntAttributeFromH5("path/to/file.h5", location = , name = )
#' }
GetIntAttributeFromH5 <- function(filename, location, name) {
  .Call("GetIntAttributeFromH5", filename, location, name)
}

# GetUintAttributeFromH5 -----------------------------------------------------
#' Reads an unsigned integer attribute from the HDF5 file.
#'
#' \code{GetUintAttributeFromH5} reads an unsigned integer attribute from the
#' HDF5 file.
#'
#' Used to read attributes not available from \code{GetH5Descriptor}. Unsigned
#' integers are returned as numeric values in R.
#'
#' @param filename Path/filename of the HDF5 file.
#' @param location Location of the group or dataset where the attribute is attached to.
#' @param name Attribute name.
#' @return An numeric attribute.
#'
#' @examples
#' \dontrun{
#' GetUintAttributeFromH5("path/to/file.h5", location = , name = )
#' }
GetUintAttributeFromH5 <- function(filename, location, name) {
  .Call("GetUintAttributeFromH5", filename, location, name)
}

# GetInt64AttributeFromH5 ----------------------------------------------------
#' Reads a 64-bit integer attribute from the HDF5 file.
#'
#' \code{GetInt64AttributeFromH5} reads a 64-bit integer attribute from the
#' HDF5 file.
#'
#' Used to read attributes not available from \code{GetH5Descriptor}. int64
#' parameters are returned as string by the dll and then converted to integer64
#' using \code{\link[bit64]{as.integer64}}.
#'
#' @param filename Path/filename of the HDF5 file.
#' @param location Location of the group or dataset where the attribute is attached to.
#' @param name Attribute name.
#' @return An int64 attribute.
#'
#' @examples
#' \dontrun{
#' GetInt64AttributeFromH5("path/to/file.h5", location = , name = )
#' }
GetInt64AttributeFromH5 <- function(filename, location, name) {
  bit64::as.integer64(.Call("GetInt64AttributeFromH5", filename, location, name))
}

# GetUint64AttributeFromH5 ---------------------------------------------------
#' Reads an unsigned 64-bit integer attribute from the HDF5 file.
#'
#' \code{GetUint64AttributeFromH5} reads an unsigned 64-bit integer attribute
#' from the HDF5 file.
#'
#' Used to read attributes not available from \code{GetH5Descriptor}. Unsigned
#' int64 parameters are returned as string by the dll and then converted to integer64
#' using \code{\link[bit64]{as.integer64}} (should work most of the time).
#'
#' @param filename Path/filename of the HDF5 file.
#' @param location Location of the group or dataset where the attribute is attached to.
#' @param name Attribute name.
#' @return An int64 attribute.
#'
#' @examples
#' \dontrun{
#' GetUint64AttributeFromH5("path/to/file.h5", location = , name = )
#' }
GetUint64AttributeFromH5 <- function(filename, location, name) {
  bit64::as.integer64(.Call("GetUint64AttributeFromH5", filename, location, name))
}

# GetFloatAttributeFromH5 ----------------------------------------------------
#' Reads a float attribute from the HDF5 file.
#'
#' \code{GetFloatAttributeFromH5} reads a float attribute from the HDF5 file.
#'
#' Used to read attributes not available from \code{GetH5Descriptor}.
#'
#' @param filename Path/filename of the HDF5 file.
#' @param location Location of the group or dataset where the attribute is attached to.
#' @param name Attribute name.
#' @return A numeric attribute.
#'
#' @examples
#' \dontrun{
#' GetFloatAttributeFromH5("path/to/file.h5", location = , name = )
#' }
GetFloatAttributeFromH5 <- function(filename, location, name) {
  .Call("GetFloatAttributeFromH5", filename, location, name)
}

# GetDoubleAttributeFromH5 ---------------------------------------------------
#' Reads a double attribute from the HDF5 file.
#'
#' \code{GetDoubleAttributeFromH5} reads a double attribute from the HDF5 file.
#'
#' Used to read attributes not available from \code{GetH5Descriptor}.
#'
#' @param filename Path/filename of the HDF5 file.
#' @param location Location of the group or dataset where the attribute is attached to.
#' @param name Attribute name.
#' @return A numeric attribute.
#'
#' @examples
#' \dontrun{
#' GetDoubleAttributeFromH5("path/to/file.h5", location = , name = )
#' }
GetDoubleAttributeFromH5 <- function(filename, location, name) {
  .Call("GetDoubleAttributeFromH5", filename, location, name)
}

# GetStringAttributeFromH5 ---------------------------------------------------
#' Reads a string attribute from the HDF5 file.
#'
#' \code{GetStringAttributeFromH5} reads a string attribute from the HDF5 file.
#'
#' Used to read attributes not available from \code{GetH5Descriptor}.
#'
#' @param filename Path/filename of the HDF5 file.
#' @param location Location of the group or dataset where the attribute is attached to.
#' @param name Attribute name.
#' @return A string attribute.
#'
#' @examples
#' \dontrun{
#' GetStringAttributeFromH5("path/to/file.h5", location = , name = )
#' }
GetStringAttributeFromH5 <- function(filename, location, name) {
  .Call("GetStringAttributeFromH5", filename, location, name)
}

# SetIntAttributeInH5 --------------------------------------------------------
#' Writes an integer attribute to the HDF5 file.
#'
#' \code{SetIntAttributeInH5} writes an integer attribute to the HDF5 file.
#'
#' @param filename Path/filename of the HDF5 file.
#' @param location Location of the group or dataset where the attribute is attached to.
#' @param name Attribute name.
#' @param attribute Integer attribute.
#'
#' @examples
#' \dontrun{
#' SetIntAttributeInH5("path/to/file.h5", location = , name = , attribute = )
#' }
SetIntAttributeInH5 <- function(filename, location, name, attribute) {
  .Call("SetIntAttributeInH5", filename, location, name, attribute)
}

# SetDoubleAttributeInH5 -----------------------------------------------------
#' Writes a numeric attribute to the HDF5 file.
#'
#' \code{SetDoubleAttributeInH5} writes a numeric attribute to the HDF5 file.
#'
#' @param filename Path/filename of the HDF5 file.
#' @param location Location of the group or dataset where the attribute is attached to.
#' @param name Attribute name.
#' @param attribute Numeric attribute.
#'
#' @examples
#' \dontrun{
#' SetDoubleAttributeInH5("path/to/file.h5", location = , name = , attribute = )
#' }
SetDoubleAttributeInH5 <- function(filename, location, name, attribute) {
  .Call("SetDoubleAttributeInH5", filename, location, name, attribute)
}

# SetStringAttributeInH5 -----------------------------------------------------
#' Writes a string attribute to the HDF5 file.
#'
#' \code{SetStringAttributeInH5} writes a string attribute to the HDF5 file.
#'
#' @param filename Path/filename of the HDF5 file.
#' @param location Location of the group or dataset where the attribute is attached to.
#' @param name Attribute name.
#' @param attribute String attribute (max. 256 characters).
#'
#' @examples
#' \dontrun{
#' SetIntAttributeInH5("path/to/file.h5", location = , name = , attribute = )
#' }
SetStringAttributeInH5 <- function(filename, location, name, attribute) {
  .Call("SetStringAttributeInH5", filename, location, name, attribute)
}

# GetUserDataFromH5 ----------------------------------------------------------
#' Reads user data from the HDF5 file.
#'
#' \code{GetUserDataFromH5} reads a row of user data added to the file using
#' the \code{AddUserData} function.
#'
#'  If you want to access data saved by a registered data source (or a data
#'  source plugin, which uses the same mechanism) use the \code{\link{GetRegUserDataFromH5}}
#'  function.
#'
#' @param filename Path/filename of the HDF5 file.
#' @param location Location of the group or dataset where the attribute is attached to.
#' @param rowIndex Row index.
#' @return A list containing the row of user data and data description.
#'
#' @examples
#' \dontrun{
#' GetUserDataFromH5("path/to/file.h5", location = , rowIndex = )
#' }
GetUserDataFromH5 <- function(filename, location, rowIndex) {
  .Call("GetUserDataFromH5", filename, location, rowIndex)
}

# GetAcquisitionLogFromH5 ----------------------------------------------------
#' Reads a single acquisition log entry.
#'
#' \code{GetAcquisitionLogFromH5} reads a single acquisition log entry and
#' returns the timestamp and the log text.
#'
#' @param filename Path/filename of the HDF5 file.
#' @param index Index of log entry.
#' @return A list containing the timestamp (as a string) and log text. Use
#' \code{bit64::as.integer64(timestamp)} to convert the timestamp into a 64-bit
#' integer value.
#'
#' @examples
#' \dontrun{
#' GetAcquisitionLogFromH5("path/to/file.h5", index = )
#' }
GetAcquisitionLogFromH5 <- function(filename, index) {
  .Call("GetAcquisitionLogFromH5", filename, index)
}

# GetEventListSpectrumFromH5 -------------------------------------------------
#' Reads the events of a spectrum from HDF5 data file.
#'
#' \code{GetEventListSpectrumFromH5} reads the events of a single spectrum
#' given by segment, buf and write indices from HDF5 data file.
#'
#' Note: only uncompressed event lists are supported.
#'
#' @param filename Path/filename of the HDF5 file.
#' @param segmentIndex Segment index.
#' @param bufIndex Buf index.
#' @param writeIndex Write index.
#' @return A vector containing the event data.
#'
#' @examples
#' \dontrun{
#' GetEventListSpectrumFromH5("path/to/file.h5", segmentIndex = 0,
#' bufIndex = 0, writeIndex = 0)
#' }
GetEventListSpectrumFromH5 <- function(filename, segmentIndex, bufIndex,
                                         writeIndex) {
  .Call("GetEventListSpectrumFromH5", filename, segmentIndex, bufIndex, writeIndex)
}

# Not implemented --------------------------------------------------------------
# TwGetBufWriteProfileFromH5
# TwGetBufWriteProfileFromH5_2
# TwSetUintAttributeInH5
# TwSetInt64AttributeInH5
# TwSetUInt64AttributeInH5
# TwSetFloatAttributeInH5
# TwChangePeakTable
# TwChangePeakTable2
# TwChangePeakFromFile
# TwChangePeakFromFile2
# TwProgressCallback
# TwProgressCallback2
# TwReadRawData
# TwH5SetMassCalib
# TwH5SetMassCalibEx
# TwGetEventListDataFromH5
# TwGetEventListBlobFromH5
# TwFreeEventListData
# TwFreeEventListData2
# TwH5GetMassCalibPar
# TwMultiPeakFitIntegration
# TwH5AddLogEntry
# TwH5AddUserDataMultiRow
# TwH5SetMassCalibDynamic
# TwGenerateSegmentProfilesFromEventList
