#### Control functions ---------------------------------------------------------

# InitializeDll -----------------------------------------------------------------
#' Initializes the TofDaqDll.dll.
#'
#' \code{InitializeDll} initializes the TofDaqDll.dll. It is usually not necessary
#' to call \code{InitializeDll} explicitely, as it is called automatically by
#' functions that need the DLL to be in an initialized state.
InitializeDll <- function() {
  rv <- .Call("InitializeDll")
  if (rv!=4) {return(warning(TwRetVal[rv+1]))} else {return(TwRetVal[5])}
}

# CleanupDll -----------------------------------------------------------------
#' Deinitializes the TofDaqDll.dll.
#'
#' \code{CleanupDll} deinitializes the TofDaqDll.dll (frees allocated memory,
#' releases mapped shared memory and closes open files). This function is
#' automatically called when the TofDaqR packages is unloaded.
CleanupDll <- function() {
  invisible(.Call("CleanupDll"))
}

# GetDllVersion --------------------------------------------------------------
#' Gets the version number of the TofDaq API.
#'
#' \code{GetDllVersion} gets the version number of the TofDaq API.
GetDllVersion <- function() {
  .Call("GetDllVersion")
}

# TofDaqRunning --------------------------------------------------------------
#' Checks if TofDaq recorder application is running.
#'
#' \code{TofDaqRunning} checks if TofDaq recorder application is running.
#'
#' @return \code{TRUE} or \code{FALSE}.
TofDaqRunning <- function() {
  .Call("TofDaqRunning")
}

# DaqActive ------------------------------------------------------------------
#' Checks if TofDaq recorder is currently acquiring data.
#'
#' \code{DaqActive} checks if TofDaq recorder is currently acquiring data.
#'
#' @return \code{TRUE} or \code{FALSE}.
DaqActive <- function() {
  .Call("DaqActive")
}

# StartAcquisition -----------------------------------------------------------
#' Starts an acquisition.
#'
#' \code{StartAcquisition} starts an acquisition.
StartAcquisition <- function() {
  rv <- .Call("StartAcquisition")
  if (rv!=4) {return(warning(TwRetVal[rv+1]))} else {return(TwRetVal[5])}
}

# StopAcquisition ------------------------------------------------------------
#' Stops the current acquisition.
#'
#' \code{StopAcquisition} stops the current acquisition.
StopAcquisition <- function() {
  rv <- .Call("StopAcquisition")
  if (rv!=4) {return(warning(TwRetVal[rv+1]))} else {return(TwRetVal[5])}
}

# CloseTofDaqRec -------------------------------------------------------------
#' Closes the TofDaq recorder application.
#'
#' \code{CloseTofDaqRec} closes the TofDaq recorder application.
CloseTofDaqRec <- function() {
  rv <- .Call("CloseTofDaqRec")
  if (rv!=4) {return(warning(TwRetVal[rv+1]))} else {return(TwRetVal[5])}
}

# InitializeDaqDevice --------------------------------------------------------
#' Initializes the DAQ board.
#'
#' \code{InitializeDaqDevice} initializes the DAQ board (this is also done at
#' startup of TofDaqRec.exe). This can take up to 8 seconds depending on the
#' actual DAQ hardware.
InitializeDaqDevice <- function() {
  rv <- .Call("InitializeDaqDevice")
  if (rv!=4) {return(warning(TwRetVal[rv+1]))} else {return(TwRetVal[5])}
}

# SetTimeout -----------------------------------------------------------------
#' Sets the timeout.
#'
#' \code{SetTimeout} sets the global timeout for all functions that can time
#' out. Default is 500 ms.
#'
#' @param timeout Timeout in ms. Default is 500 ms.
SetTimeout <- function(timeout = 500) {
  invisible(.Call("SetTimeout", timeout))
}

# GetTimeout -----------------------------------------------------------------
#' Gets the timeout.
#'
#' \code{GetTimeout} gets the current timeout value (in ms).
GetTimeout <- function() {
  .Call("GetTimeout")
}

# AutoSetupDaqDevice ---------------------------------------------------------
#' Auto setup routine for the DAQ device.
#'
#' Auto setup routine for the DAQ device. Currently implemented only for
#' AP240 averager and ndigo5G.
AutoSetupDaqDevice <- function() {
  rv <- .Call("AutoSetupDaqDevice")
  if (rv!=4) {return(warning(TwRetVal[rv+1]))} else {return(TwRetVal[5])}
}

# OnDemandMassCalibration ----------------------------------------------------
#' Arms/executes on demand mass calibration.
#'
#' \code{OnDemandMassCalibration} arms/executes on demand mass calibration.
#' Requires that parameter \code{ReCalibFreq} is set to 3 (on demand) in order to work.
#'
#' @param action 0: arms mass calibration, 1: updates mass calibration using
#' data acquired since previous arm.
OnDemandMassCalibration <- function(action) {
  .Call("OnDemandMassCalibration", action)
}


#### Configuration functions ---------------------------------------------------

# ShowConfigWindow -----------------------------------------------------------
#' Shows the TofDaq recorder configuration windows.
#'
#' \code{ShowConfigWindow} shows the different tabs of the TofDaq recorder
#' configuration window.
#'
#' @param ConfigWindowIndex Index of configuration tab to show (valid range: 0-6)
ShowConfigWindow <- function(ConfigWindowIndex = 0) {
  rv <- .Call("ShowConfigWindow", ConfigWindowIndex)
  if (rv!=4) {return(warning(TwRetVal[rv+1]))} else {return(TwRetVal[5])}
}

# LoadIniFile ----------------------------------------------------------------
#' Loads a configuration file.
#'
#' \code{LoadIniFile} loads a configuration file (*.ini) from disk.
#'
#' @param IniFile Path/filename of the configuration file. If \code{IniFile} is
#' an empty string, "TwApiTmpIni.ini" in the TofDaq recorder directory will be
#' used.
LoadIniFile <- function(IniFile) {
  rv <- .Call("LoadIniFile", IniFile)
  if (rv!=4) {return(warning(TwRetVal[rv+1]))} else {return(TwRetVal[5])}
}

# SaveIniFile ----------------------------------------------------------------
#' Saves the current configuration (*.ini) to disk.
#'
#' \code{SaveIniFile} saves the current configuration (*.ini) to disk.
#'
#' @param IniFile Path/filename of the configuration file. If \code{IniFile} is
#' an empty string, "TwApiTmpIni.ini" in the TofDaq recorder directory will be
#' used.
SaveIniFile <- function(IniFile) {
  rv <- .Call("SaveIniFile", IniFile)
  if (rv!=4) {return(warning(TwRetVal[rv+1]))} else {return(TwRetVal[5])}
}

# GetDaqParameter ------------------------------------------------------------
#' Gets a single parameter as a string.
#'
#' \code{GetDaqParameter} gets a single parameter as a string.
#'
#' @param Parameter Parameter name as a string. See
#' \emph{/doc/TofDaqDll.htm#parameter_list} for a list of all available parameters.
GetDaqParameter <- function(Parameter) {
  .Call("GetDaqParameter", Parameter)
}

# SetDaqParameter ------------------------------------------------------------
#' Sets a single parameter.
#'
#' \code{SetDaqParameter} sets a single parameter.
#'
#' @param Parameter Parameter name as a string. See
#' \emph{/doc/TofDaqDll.htm#parameter_list} for a list of all available parameters.
#' @param ValueString Value as a string.
SetDaqParameter <- function(Parameter, ValueString) {
  rv <- .Call("SetDaqParameter", Parameter, ValueString)
  if (rv!=4) {return(warning(TwRetVal[rv+1]))} else {return(TwRetVal[5])}
}


#### Data access functions -----------------------------------------------------

# GetDescriptor --------------------------------------------------------------
#' Gets various information about the active acquisition.
#'
#' \code{GetDescriptor} retrieves the current TSharedMemoryDesc structure.
#' TSharedMemoryDesc contains various static information about the active
#' acquisition that can be retrieved by \code{GetDaqParameter} functions but
#' also information of DAQ progress.
#' See \emph{/doc/TofDaqDll.htm} for more details.
#'
#' Unsigned int parameters are returned as numeric. int64 and unsigned int64
#' parameters are returned as string by the dll and then converted to integer64
#' using \code{\link[bit64]{as.integer64}}.
#'
#' @return A list containing the TSharedMemoryDesc structure
GetDescriptor <- function() {
  tmp <- .Call("GetDescriptor")
  if (is.null(tmp)) {
    return()
  } else {
    tmp$BlockPeriod <- bit64::as.integer64(tmp$BlockPeriod)
    tmp$BlockPulseDelay <- bit64::as.integer64(tmp$BlockPulseDelay)
    tmp$BlockDelay <- bit64::as.integer64(tmp$BlockDelay)
    tmp$AcquisitionLogTime <- bit64::as.integer64(tmp$AcquisitionLogTime)
    tmp$TimeZero <- bit64::as.integer64(tmp$TimeZero)
    return(tmp)
  }
}

# GetPeakParameters ----------------------------------------------------------
#' Gets parameters for a given peak.
#'
#' \code{GetPeakParameters} gets parameters for a given peak.
#'
#' @param PeakIndex Index of peak.
#' @return A list with the peak paramters \emph{label}, \emph{mass}, \emph{loMass} and \emph{hiMass}.
GetPeakParameters <- function(PeakIndex) {
  rv <- .Call("GetPeakParameters", PeakIndex)
  if (rv$TwRetVal!=4) {
    warning(TwRetVal[rv$TwRetVal+1])
    return()
  } else {
    rv$TwRetVal <- NULL
    return(rv)
  }
}

# ReleaseSharedMemory ----------------------------------------------------------
#' Manually releases the shared memory acquisition buffers.
#'
#' \code{ReleaseSharedMemory} manually releases the shared memory acquisition
#' buffers. This is needed if \code{\link{KeepSharedMemMapped}} has been set to
#' \code{TRUE}.
ReleaseSharedMemory <- function() {
  rv <- .Call("ReleaseSharedMemory")
  if (rv!=4) {return(warning(TwRetVal[rv+1]))} else {return(TwRetVal[5])}
}

# WaitForNewData -------------------------------------------------------------
#' Waits for new data.
#'
#' \code{WaitForNewData} waits for new data. Returns when new data is
#' available or when timed out.
#'
#' @param timeout Timeout in ms.
#' @param WaitForEventReset If \code{TRUE} (default) waits for application to
#' reset data available event before returning.
WaitForNewData <- function(timeout, WaitForEventReset = TRUE) {
  rv <- .Call("WaitForNewData", timeout, WaitForEventReset)
  if (rv!=4) {return(warning(TwRetVal[rv+1]))} else {return(TwRetVal[5])}
}

# WaitForEndOfAcquisition ----------------------------------------------------
#' Waits for the end of the current acquisition.
#'
#' \code{WaitForEndOfAcquisition} waits for the end of the current acquisition.
#' If \code{NbrRuns > 1} this function waits for the end of the last acquisition.
#'
#' @param timeout Timeout in ms.
WaitForEndOfAcquisition <- function(timeout) {
  rv <- .Call("WaitForEndOfAcquisition", timeout)
  if (rv!=4) {return(warning(TwRetVal[rv+1]))} else {return(TwRetVal[5])}
}

# GetSumSpectrumFromShMem ----------------------------------------------------
#' Sum spectrum from shared memory.
#'
#' \code{GetSumSpectrumFromShMem} gets the sum spectrum from shared memory.
#'
#' @param Normalize If \code{FALSE} the spectrum is reported as sum,
#' if \code{TRUE} (default) the spectrum is normalized to counts per extraction.
GetSumSpectrumFromShMem <- function(Normalize = TRUE) {
  .Call("GetSumSpectrumFromShMem", Normalize)
}

# GetTofSpectrumFromShMem ----------------------------------------------------
#' Single TOF spectrum from shared memory.
#'
#' \code{GetTofSpectrumFromShMem} reads a single TOF spectrum (possibly
#' averaged/summed over segment dimension) from shared memory. If
#' \code{SegmentIndex = SegmentEndIndex = -1} the complete block of data is
#' copied and the \code{Normalize} flag is ignored.
#'
#' @param SegmentIndex Segment start index of data to fetch (or -1 for complete
#' block copy).
#' @param SegmentEndIndex Segment end index of data to fetch (or -1 for complete
#' block copy).
#' @param BufIndex Buf index of data to fetch.
#' @param Normalize If \code{FALSE} the spectrum is reported as sum,
#' if \code{TRUE} (default) the spectrum is normalized to counts per extraction.
#' @return A vector containing the mass spectrum or an array containing the
#' block of mass spectra if \code{SegmentIndex = SegmentEndIndex = -1}.
GetTofSpectrumFromShMem <- function(SegmentIndex = 0, SegmentEndIndex = 0,
                                      BufIndex, Normalize = TRUE) {
  .Call("GetTofSpectrumFromShMem", SegmentIndex, SegmentEndIndex, BufIndex, Normalize)
}

# GetSpecXaxisFromShMem ------------------------------------------------------
#' X-axis values of mass spectrum.
#'
#' \code{GetSpecXaxisFromShMem} returns an array of x-axis values of the mass
#' spectrum.
#'
#' @param Type x-axis type (0: sample index, 1: mass/charge [Th],
#' -1: mass/charge [Th] (2nd TOF), 2: time of flight [microsec],
#' -2: time of flight [microsec] (2nd TOF), 3: frequency [kHz]).
#' @return A vector containing the x-axis values.
GetSpecXaxisFromShMem <- function(Type) {
  .Call("GetSpecXaxisFromShMem", Type)
}

# GetStickSpectrumFromShMem --------------------------------------------------
#' Single stick spectrum from shared memory.
#'
#' \code{GetStickSpectrumFromShMem} reads a single stick spectrum from shared
#' memory.
#'
#' @param SegmentIndex Segment start index of data to fetch.
#' @param SegmentEndIndex Segment end index of data to fetch.
#' @param BufIndex Buf index of data to fetch.
#' @return A list containing the stick spectrum and corresponding masses.
GetStickSpectrumFromShMem <- function(SegmentIndex = 0, SegmentEndIndex = 0, BufIndex) {
  .Call("GetStickSpectrumFromShMem", SegmentIndex, SegmentEndIndex, BufIndex)
}

# GetSegmentProfileFromShMem -------------------------------------------------
#' Segment profile for a given peak and buf index from shared memory.
#'
#' \code{GetSegmentProfileFromShMem} reads the segment profile for a given
#' peak and buf index from shared memory. Use -1 for \code{PeakIndex} to get
#' segment profiles of all peaks.
#'
#' @param PeakIndex Index of peak to fetch segment profile from. All peaks are
#' read if \code{PeakIndex = -1}.
#' @param BufIndex Buf index of data to fetch.
#' @return A vector containing the segment profile(s).
#'
GetSegmentProfileFromShMem <- function(PeakIndex, BufIndex) {
  .Call("GetSegmentProfileFromShMem", PeakIndex, BufIndex)
}

# GetBufTimeFromShMem --------------------------------------------------------
#' Time stamp for a given buf and write.
#'
#' \code{GetBufTimeFromShMem} reads the time stamp for a given buf and write
#' from shared memory.
#'
#' @param BufIndex Buf index.
#' @param WriteIndex Write index.
#' @return A time stamp (in seconds relative to acquisition start).
GetBufTimeFromShMem <- function(BufIndex, WriteIndex) {
  .Call("GetBufTimeFromShMem", BufIndex, WriteIndex)
}


#### Data storage functions ----------------------------------------------------

# AddLogEntry ----------------------------------------------------------------
#' Adds an entry to the acquisition log.
#'
#' \code{AddLogEntry} adds an entry to the acquisition log.
#'
#' Note: specifying the log entry time (for adding entries in the past) is not
#' supported.
#'
#'  @param LogEntryText Log text (max. 255 characters).
#'
#'  @family Data storage functions
AddLogEntry <- function(LogEntryText) {
  rv <- .Call("AddLogEntry", LogEntryText)
  if (rv!=4) {return(warning(TwRetVal[rv+1]))} else {return(TwRetVal[5])}
}

# AddAttributeInt ------------------------------------------------------------
#' Attaches an integer attribute to the current HDF5 file.
#'
#' \code{AddAttributeInt} attaches an integer attribute to the current HDF5 file.
#'
#' @param Object HDF5 object (group or dataset) to attach attribute (max. 255
#' characters).
#' @param AttributeName Attribute name (max. 127 characters).
#' @param Value Attribute value (integer type).
#'
#' @family Data storage functions
AddAttributeInt <- function(Object, AttributeName, Value) {
  rv <- .Call("AddAttributeInt", Object, AttributeName, Value)
  if (rv!=4) {return(warning(TwRetVal[rv+1]))} else {return(TwRetVal[5])}
}

# AddAttributeDouble ---------------------------------------------------------
#' Attaches a numeric attribute to the current HDF5 file.
#'
#' \code{AddAttributeDouble} attaches a numeric attribute to the current HDF5 file.
#'
#' @param Object HDF5 object (group or dataset) to attach attribute (max. 255
#' characters).
#' @param AttributeName Attribute name (max. 127 characters).
#' @param Value Attribute value (numeric type).
#'
#' @family Data storage functions
AddAttributeDouble <- function(Object, AttributeName, Value) {
  rv <- .Call("AddAttributeDouble", Object, AttributeName, Value)
  if (rv!=4) {return(warning(TwRetVal[rv+1]))} else {return(TwRetVal[5])}
}

# AddAttributeString ---------------------------------------------------------
#' Attaches a string attribute to the current HDF5 file.
#'
#' \code{AddAttributeString} attaches a string attribute to the current HDF5 file.
#'
#' @param Object HDF5 object (group or dataset) to attach attribute (max. 255
#' characters).
#' @param AttributeName Attribute name (max. 127 characters).
#' @param Value Attribute string value (max. 255 characters).
#'
#' @family Data storage functions
AddAttributeString <- function(Object, AttributeName, Value) {
  rv <- .Call("AddAttributeString", Object, AttributeName, Value)
  if (rv!=4) {return(warning(TwRetVal[rv+1]))} else {return(TwRetVal[5])}
}

# AddUserData ----------------------------------------------------------------
#' Stores (asynchronous) user supplied data.
#'
#' \code{AddUserData} stores user supplied data asynchronously to the TOF data
#' acquistion into the current data file. Creates datasets "Data" and "Info" at
#' \code{Location}.
#'
#' @param Location Location of group in HDF5 file where the datasets are created.
#' @param NbrElements Number of elements to store (per call to this function),
#' maximum is 1048575.
#' @param ElementDescription Vector of length \code{NbrElements} containing the
#' text description of elements. If \code{ElementDescription} is an empty string
#' the dataset "Info" is not created.
#' @param Data Vector of length \code{NbrElements} containing the data to be
#' stored in dataset "Data".
#' @param CompressionLevel ZLIB compression level (0-9) for dataset creation.
#' If the dataset at Location already exists this parameter has no effect.
#'
#' @family Data storage functions
AddUserData <- function(Location, NbrElements, ElementDescription, Data, CompressionLevel) {
  label <- ""
  for (i in 1:NbrElements) {
    tmp <- paste(rep(" ", 256), collapse = "")
    substr(tmp, 1, nchar(ElementDescription[i])) <- ElementDescription[i]
    label <- paste0(label, tmp)
  }
  rv <- .Call("AddUserData", Location, NbrElements, ElementDescription = label, Data, CompressionLevel)
  if (rv!=4) {return(warning(TwRetVal[rv+1]))} else {return(TwRetVal[5])}
}

# AddUserDataMultiRow --------------------------------------------------------
#' Stores (asynchronous) user supplied data.
#'
#' \code{AddUserDataMultiRow} stores user supplied data asynchronously to the TOF data
#' acquistion into the current data file. Creates datasets "Data" and "Info" at
#' \code{Location}.
#'
#' Same as \code{AddUserData}, but adds argument \code{NbrRows} to add several
#' lines of user data at once.
#'
#' @param Location Location of group in HDF5 file where the datasets are created.
#' @param NbrElements Number of elements to store (per call to this function),
#' maximum is 1048575.
#' @param NbrRows Number of rows to store per call to this function (each row
#' contains \code{NbrElements} entries), maximum is 2047.
#' @param ElementDescription Vector of length \code{NbrElements} containing the
#' text description of elements. If \code{ElementDescription} is an empty string
#' the dataset "Info" is not created.
#' @param Data Vector of length \code{NbrElements} containing the data to be
#' stored in dataset "Data".
#' @param CompressionLevel ZLIB compression level (0-9) for dataset creation.
#' If the dataset at Location already exists this parameter has no effect.
#'
#' @family Data storage functions
AddUserDataMultiRow <- function(Location, NbrElements, NbrRows, ElementDescription, Data, CompressionLevel) {
  label <- ""
  for (i in 1:NbrElements) {
    tmp <- paste(rep(" ", 256), collapse = "")
    substr(tmp, 1, nchar(ElementDescription[i])) <- ElementDescription[i]
    label <- paste0(label, tmp)
  }
  rv <- .Call("AddUserDataMultiRow", Location, NbrElements, NbrRows, ElementDescription = label, Data, CompressionLevel)
  if (rv!=4) {return(warning(TwRetVal[rv+1]))} else {return(TwRetVal[5])}
}

# RegisterUserDataBuf --------------------------------------------------------
#' Registers a data source to store (synchronous) user supplied data.
#'
#' \code{RegisterUserDataBuf} registers a data source to store user supplied
#' data synchronously to the TOF data acquistion (every buf) into the data file
#' being currently recorded. Creates datasets "TwData" and "TwInfo" at \code{Location}.
#'
#' Needs to be executed before starting the acquisition.
#' Use \code{\link{UpdateUserData}} to actually store the data.
#'
#' @param Location Location of group in HDF5 file where the datasets are created.
#' @param NbrElements Number of elements to store per buf.
#' @param ElementDescription Vector of length \code{NbrElements} containing the
#' text description of elements. If \code{ElementDescription} is an empty string
#' the dataset "TwInfo" is not created.
#' @param CompressionLevel Compression level used for data storage (0: no
#' compression, 1-9: increasing levels of compression (and CPU load)).
#'
#' @family Data storage functions
RegisterUserDataBuf <- function(Location, NbrElements, ElementDescription, CompressionLevel) {
  label <- ""
  for (i in 1:NbrElements) {
    tmp <- paste(rep(" ", 256), collapse = "")
    substr(tmp, 1, nchar(ElementDescription[i])) <- ElementDescription[i]
    label <- paste0(label, tmp)
  }
  rv <- .Call("RegisterUserDataBuf", Location, NbrElements, ElementDescription = label, CompressionLevel)
  if (rv!=4) {return(warning(TwRetVal[rv+1]))} else {return(TwRetVal[5])}
}

# RegisterUserDataWrite ------------------------------------------------------
#' Registers a data source to store (synchronous) user supplied data.
#'
#' \code{RegisterUserDataWrite} registers a data source to store user supplied
#' data synchronously to the TOF data acquistion (every write) into the data
#' file being currently recorded. Creates datasets "TwData" and "TwInfo" at
#' \code{Location}.
#'
#' Needs to be executed before starting the acquisition.
#' Use \code{\link{UpdateUserData}} to actually store the data.
#'
#' @param Location Location of group in HDF5 file where the datasets are created.
#' @param NbrElements Number of elements to store per write.
#' @param ElementDescription Vector of length \code{NbrElements} containing the
#' text description of elements. If \code{ElementDescription} is an empty string
#' the dataset "TwInfo" is not created.
#' @param CompressionLevel Compression level used for data storage (0: no
#' compression, 1-9: increasing levels of compression (and CPU load)).
#'
#' @family Data storage functions
RegisterUserDataWrite <- function(Location, NbrElements, ElementDescription, CompressionLevel) {
  label <- ""
  for (i in 1:NbrElements) {
    tmp <- paste(rep(" ", 256), collapse = "")
    substr(tmp, 1, nchar(ElementDescription[i])) <- ElementDescription[i]
    label <- paste0(label, tmp)
  }
  rv <- .Call("RegisterUserDataWrite", Location, NbrElements, ElementDescription = label, CompressionLevel)
  if (rv!=4) {return(warning(TwRetVal[rv+1]))} else {return(TwRetVal[5])}
}

# RegisterUserDataNoStore ------------------------------------------------------
#' Registers a data source for (synchronous) user supplied data.
#'
#' \code{RegisterUserDataNoStore} registers a data for user supplied
#' data (synchronous to the TOF data acquistion) but the data is not stored
#' in the data file.
#'
#' Needs to be executed before starting the acquisition.
#' Use \code{\link{UpdateUserData}} to update the data.
#'
#' @param Location Location of group in HDF5 file where the datasets are created.
#' @param NbrElements Number of elements to store per write.
#' @param ElementDescription Vector of length \code{NbrElements} containing the
#' text description of elements. If \code{ElementDescription} is an empty string
#' the dataset "TwInfo" is not created.
#'
#' @family Data storage functions
RegisterUserDataNoStore <- function(Location, NbrElements, ElementDescription) {
  label <- ""
  for (i in 1:NbrElements) {
    tmp <- paste(rep(" ", 256), collapse = "")
    substr(tmp, 1, nchar(ElementDescription[i])) <- ElementDescription[i]
    label <- paste0(label, tmp)
  }
  rv <- .Call("RegisterUserDataNoStore", Location, NbrElements, ElementDescription = label)
  if (rv!=4) {return(warning(TwRetVal[rv+1]))} else {return(TwRetVal[5])}
}

# UnregisterUserData ---------------------------------------------------------
#' Unregisters a data source.
#'
#' \code{UnregisterUserData} unregisters a data source previously registered
#' with \code{\link{RegisterUserDataBuf}} or \code{\link{RegisterUserDataWrite}}.
#'
#' @param Location Location of group in HDF5 file identifying the user data.
#'
#' @family Data storage functions
UnregisterUserData <- function(Location) {
  rv <- .Call("UnregisterUserData", Location)
  if (rv!=4) {return(warning(TwRetVal[rv+1]))} else {return(TwRetVal[5])}
}

# UpdateUserData -------------------------------------------------------------
#' Updates the values for a registered data source.
#'
#' \code{UpdateUserData} updates the values for a registered data source.
#'
#' @param Location Location of group in HDF5 file where the datasets are created.
#' @param NbrElements Number of elements to update.
#' @param Data Vector of length \code{NbrElements} containing the new data.
#'
#' @family Data storage functions
UpdateUserData <- function(Location, NbrElements, Data) {
  rv <- .Call("UpdateUserData", Location, NbrElements, Data)
  if (rv!=4) {return(warning(TwRetVal[rv+1]))} else {return(TwRetVal[5])}
}

# ReadRegUserData ------------------------------------------------------------
#' Reads the current values of a registered data source.
#'
#' \code{ReadRegUserData} reads the current values of a registered data source.
#'
#' @param Location Location of group in HDF5 file where the datasets are created.
#' @param NbrElements Number of elements to read.
#' @return Vector containing the registered user data.
#'
#' @family Data storage functions
ReadRegUserData <- function(Location, NbrElements) {
  .Call("ReadRegUserData", Location, NbrElements)
}

# QueryRegUserDataSize -------------------------------------------------------
#' Queries the size of a registered data source.
#'
#' \code{QueryRegUserDataSize} queries the size (nbrElements) of a registered data source.
#'
#' @param Location Location of group in HDF5 file identifying the registered data.
#'
#' @family Data storage functions
QueryRegUserDataSize <- function(Location) {
  .Call("QueryRegUserDataSize", Location)
}

# GetRegUserDataSources -------------------------------------------------------
#' Queries names, dimensions and types of all registered data sources.
#'
#' \code{GetRegUserDataSources} queries names, dimensions and types of all
#' data sources currently registered in TofDaq recorder.
#'
#' @return List with the location, nbrElements and type of the data sources.
#' type 1: data source values are written to disk for every write,
#' type 2: data source values are written to disk for every buf.
#'
#' @family Data storage functions
GetRegUserDataSources <- function() {
  .Call("GetRegUserDataSources")
}

# GetRegUserDataDesc -------------------------------------------------------
#' Reads the element descriptions of a registered data source.
#'
#' \code{GetRegUserDataDesc} reads the element descriptions of a registered data source.
#'
#' @param Location Location of group in HDF5 file identifying the registered data.
#'
#' @family Data storage functions
GetRegUserDataDesc <- function(Location) {
  .Call("GetRegUserDataDesc", Location)
}

# KeepFileOpen ---------------------------------------------------------------
#' Allows to keep the data file open at the end of an acquisition.
#'
#' \code{KeepFileOpen} allows to keep the data file open at the end of an acquisition.
#'
#' @param keepOpen Issue \code{TRUE} (after DAQ start) to signal to recorder to
#' keep the data file open. When done adding data, issue \code{FALSE} to allow
#' the recorder to close the file.
#'
#' @family Data storage functions
KeepFileOpen <- function(keepOpen) {
  .Call("KeepFileOpen", keepOpen)
}


#### TPS control functions -----------------------------------------------------

# TpsConnect -----------------------------------------------------------------
#' Connects to a remote control enabled TPSController software.
#'
#' \code{TpsConnect} connects to a remote control enabled TPSController
#' software running on the same PC.
#'
#' @family TPS functions
TpsConnect <- function() {
  rv <- .Call("TpsConnect")
  if (rv!=4) {return(warning(TwRetVal[rv+1]))} else {return(TwRetVal[5])}
}

# TpsConnect2 ----------------------------------------------------------------
#' Connects to a local or remote TPS.
#'
#' \code{TpsConnect2} connects to a local or remote TPS (TPS1: type = 0,
#' TPS2: type = 1)
#'
#' @param ip TPS2 host name or IP.
#' @param type TPS type (0: 1st generation TPS, 1: 2nd generation TPS).
#'
#' @family TPS functions
TpsConnect2 <- function(ip = "192.168.168.2", type = 1) {
  rv <- .Call("TpsConnect2", ip, type)
  if (rv!=4) {return(warning(TwRetVal[rv+1]))} else {return(TwRetVal[5])}
}

# TpsDisconnect --------------------------------------------------------------
#' Disconnects from a remote control enabled TPSController software.
#'
#' \code{TpsDisconnect} disconnects from a remote control enabled TPSController
#' software.
#'
#' @family TPS functions
TpsDisconnect <- function() {
  rv <- .Call("TpsDisconnect")
  if (rv!=4) {return(warning(TwRetVal[rv+1]))} else {return(TwRetVal[5])}
}

# TpsGetMonitorVal -----------------------------------------------------------
#' Gets the last reported monitor value for a given module.
#'
#' \code{TpsGetMonitorValue} gets the last reported monitor value for a given
#' module.
#'
#' @param moduleCode Module code.
#'
#' @family TPS functions
TpsGetMonitorValue <- function(moduleCode) {
  .Call("TpsGetMonitorValue", moduleCode)
}

# TpsGetTargetValue ----------------------------------------------------------
#' Gets the last reported target value for a given module.
#'
#' \code{TpsGetTargetValue} gets the last reported target value for a given
#' module.
#'
#' @param moduleCode Module code.
#'
#' @family TPS functions
TpsGetTargetValue <- function(moduleCode) {
  .Call("TpsGetTargetValue", moduleCode)
}

# TpsGetLastSetValue ---------------------------------------------------------
#' Gets the last reported "last set" value for a given module.
#'
#' \code{TpsGetLastSetValue} gets the last reported "last set" value for a given
#' module.
#'
#' @param moduleCode Module code.
#'
#' @family TPS functions
TpsGetLastSetValue <- function(moduleCode) {
  .Call("TpsGetLastSetValue", moduleCode)
}

# TpsSetTargetValue ----------------------------------------------------------
#' Sets the target value for a given module.
#'
#' \code{TpsSetTargetValue} sets the target value for a given module.
#'
#' @param moduleCode Module code.
#' @param value Value to set.
#'
#' @family TPS functions
TpsSetTargetValue <- function(moduleCode, value = 0) {
  rv <- .Call("TpsSetTargetValue", moduleCode, value)
  if (rv!=4) {return(warning(TwRetVal[rv+1]))} else {return(TwRetVal[5])}
}

# TpsGetNbrModules -----------------------------------------------------------
#' Gets the number of controllable modules.
#'
#' \code{TpsGetNbrModules} gets the number of controllable modules.
#'
#' @family TPS functions
TpsGetNbrModules <- function() {
  .Call("TpsGetNbrModules")
}

# TpsGetModuleCodes ----------------------------------------------------------
#' Gets the module codes of all controllable TPS modules.
#'
#' \code{TpsGetModuleCodes} gets the module codes of all controllable TPS modules.
#'
#' @family TPS functions
TpsGetModuleCodes <- function() {
  .Call("TpsGetModuleCodes")
}

# TpsInitialize --------------------------------------------------------------
#' Initializes TPS.
#'
#' \code{TpsInitialize} initializes the TPS.
#'
#' @family TPS functions
TpsInitialize <- function() {
  rv <- .Call("TpsInitialize")
  if (rv!=4) {return(warning(TwRetVal[rv+1]))} else {return(TwRetVal[5])}
}

# TpsSetAllVoltages ----------------------------------------------------------
#' Sets all voltages.
#'
#' \code{TpsSetAllVoltages} sets all voltages.
#'
#' @family TPS functions
TpsSetAllVoltages <- function() {
  rv <- .Call("TpsSetAllVoltages")
  if (rv!=4) {return(warning(TwRetVal[rv+1]))} else {return(TwRetVal[5])}
}

# TpsShutdown ----------------------------------------------------------------
#' Shuts down TPS.
#'
#' \code{TpsShutdown} shuts down the TPS.
#'
#' @family TPS functions
TpsShutdown <- function() {
  rv <- .Call("TpsShutdown")
  if (rv!=4) {return(warning(TwRetVal[rv+1]))} else {return(TwRetVal[5])}
}

# TpsGetStatus ---------------------------------------------------------------
#' Gets the status of the TPS.
#'
#' \code{TpsGetStatus} gets the status of the TPS.
#'
#' @return List with the status information: connected?, initialized?, shutdown?,
#' ion mode changable?, ion mode supported?, current ion mode?.
#'
#' @family TPS functions
TpsGetStatus <- function() {
  if (is.null(.Call("TpsGetStatus"))) {rv <-  NULL} else {
    tmp <- intToBits(.Call("TpsGetStatus"))
    connected <- as.logical(tmp)[1]
    initialized <- as.logical(tmp)[2]
    shutdown <- as.logical(tmp)[3]
    ionmode_changable <- as.logical(tmp)[4]
    ionmode_supported <- as.logical(tmp)[5]
    current_ionmode <- strtoi(paste(rev(as.integer(tmp)[6:7]), collapse=""), base = 2)
    current_ionmode_string <- c("Inconsistent", "Undetermined", "Positive", "Negative")
    rv <- list(connected = connected, initialized = initialized, shutdown = shutdown,
               ionmode_changable = ionmode_changable, ionmode_supported = ionmode_supported,
               current_ionmode = current_ionmode_string[current_ionmode+1])
  }
  return(rv)
}

# TpsLoadSetFile -------------------------------------------------------------
#' Loads a TPS set file from disk and sets all values.
#'
#' \code{TpsLoadSetFile} loads a TPS set file from disk and sets all values.
#'
#' @param setfile Path/filename of the set file to load.
#' @section Warning:
#' This does not just load the file (as the function name might suggest), but
#' also immediately sets all values.
#'
#' @family TPS functions
TpsLoadSetFile <- function(setfile) {
  rv <- .Call("TpsLoadSetFile", setfile)
  if (rv!=4) {return(warning(TwRetVal[rv+1]))} else {return(TwRetVal[5])}
}

# TpsSaveSetFile -------------------------------------------------------------
#' Saves the current TPS settings to a file.
#'
#' \code{TpsSaveSetFile} saves the current TPS settings to a file.
#'
#' @param setfile Path/filename of the set file to save.
#'
#' @family TPS functions
TpsSaveSetFile <- function(setfile) {
  rv <- .Call("TpsSaveSetFile", setfile)
  if (rv!=4) {return(warning(TwRetVal[rv+1]))} else {return(TwRetVal[5])}
}

# TpsGetActiveFilament -------------------------------------------------------
#' Gets the currently active filament.
#'
#' \code{TpsGetActiveFilament} gets the currently active filament.
#'
#' Note that \code{TpsGetMonitorValue} does not work to query the filament
#' number. Use this function instead.
#'
#' @return Returns 0 for Filament 1, 1 for Filament 2.
#'
#' @family TPS functions
TpsGetActiveFilament <- function() {
  .Call("TpsGetActiveFilament")
}

# TpsSetActiveFilament -------------------------------------------------------
#' Sets the active filament.
#'
#' \code{TpsSetActiveFilament} sets the active filament.
#'
#' Note that \code{TpsSetTargetValue} does not work to set the filament
#' number. Use this function instead.
#'
#' @param activeFilament 0 for Filament 1, 1 for Filament 2.
#'
#' @family TPS functions
TpsSetActiveFilament <- function(activeFilament = 0) {
  rv <- .Call("TpsSetActiveFilament", activeFilament)
  if (rv!=4) {return(warning(TwRetVal[rv+1]))} else {return(TwRetVal[5])}
}

# TpsGetModuleLimits ---------------------------------------------------------
#' Gets the limits for a given TPS module.
#'
#' \code{TpsGetModuleLimits} gets the (ion mode dependent) limits for a given TPS module. Only
#' works for TPS2.
#'
#' @param moduleCode Module code.
#'
#' @family TPS functions
TpsGetModuleLimits <- function(moduleCode) {
  .Call("TpsGetModuleLimits", moduleCode)
}

# TpsChangeIonMode -----------------------------------------------------------
#' Changes ion mode and sets target values to 0.
#'
#' \code{TpsChangeIonMode} changes ion mode (and sets target values to 0).
#'
#' Note: this is an undocumented function of TofDaqDll.dll.
#'
#' @param ionMode 0: positive ion mode, 1: negative ion mode
#'
#' @family TPS functions
TpsChangeIonMode <- function(ionMode) {
  rv <- .Call("TpsChangeIonMode", ionMode)
  if (rv!=4) {return(warning(TwRetVal[rv+1]))} else {return(TwRetVal[5])}
}


# Not implemented --------------------------------------------------------------
# TwContinueAcquisition
# TwManualContinueNeeded
# TwLockBuf
# TwUnLockBuf
# TwIssueDio4Pulse
# TwSetDio4State
# TwConfigVarNbrMemories
# TwGetDaqParameterInt
# TwGetDaqParameterBool
# TwGetDaqParameterFloat
# TwGetDaqParameterInt64
# TwGetDaqParameterDouble
# TwGetDaqParameterIntRef
# TwGetDaqParameterBoolRef
# TwGetDaqParameterFloatRef
# TwGetDaqParameterInt64Ref
# TwGetDaqParameterDoubleRef
# TwGetDaqParameterStringRef
# TwSetDaqParameterInt
# TwSetDaqParameterBool
# TwSetDaqParameterFloat
# TwSetDaqParameterInt64
# TwSetDaqParameterDouble
# TwSetMassCalib
# TwSetMassCalibEx
# TwGetSharedMemory
# TwGetMassCalib
# TwGetMassCalibEx
