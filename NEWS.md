# Changes in version 0.3.11.9002

* Removed `SiProcessSpectrumFromShMem()`.

* Fixed a bug in `GetDaqParameter()`.

* Fixed a bug when unloading the package.


# Changes in version 0.3.11

* Fixed a bug in `GetTofDataSinglePeak()`.

* `FitTofDataSinglePeak()` can also fit single spectra (and not only the 
  averaged spectrum).
  
* Added TwToolDll functions `DecomposeMass()`, `DecomposeMass2()`, 
  `MatchSpectra()`.

* Fixed a bug in `H5SetMassCalibDynamic()`.

* Added sample index axis to the output of `GetTofDataSinglePeak()` and 
  `FitTofDataSinglePeak()`.
  
* Fixed reading the description in `GetRegUserDataFromH5()`.

* Updated help pages of `GetSegmentProfileFromH5()`, 
  `GetBufWriteProfileFromH5()`, `ContinueAcquisition()` and 
  `ManualContinueNeeded()`.
  
* Replaced HTTP links with HTTPS.

* New TofDaq API versions are no longer hosted on https://soft.tofwerk.com/.
  Therefore, the API libraries are downloaded from 
  https://github.com/pasturm/TofDaqR/blob/master/TofDaqAPI/.
  
* Updated to TofDaq API 1.99r1369.

* Removed GitHub Actions workflows.

* Removed the binary version for macOS. The package can still be installed from 
  source on Intel-based Macs. Unfortunately, TofDaq API libraries are not 
  available for ARM-based processors.
  
* Added TofDaqDll functions `DioStartDelayActive()`, `SendDioStartSignal()`, 
  `WaitingForDioStartSignal()`, `SaturationWarning()`, 
  `ConfigureForSingleIonMeasurement()`, `TpsLoadSetFile2()`, 
  `TpsSaveSetFileRc()`, `TpsGetNmtState()`, `TpsSetNmtState()`, 
  `TpsGetModuleProperties()`.


# Changes in version 0.3.10

* Fixed a bug, where the initialization of char arrays could cause R to crash in
  some versions of R (regression from TofDaqR version 0.3.9).
  
* Using GitHub Actions instead of Travis CI for continuous integration.


# Changes in version 0.3.9

* `DeleteAttributeInH5()`, `WaitForExclusiveFileAccess()` and 
  `WriteNetCdfTimeSeriesFile()` are missing in libtwh5 of TofDaq API 1.99r759
  (Mac and Linux), causing package installation on non-Windows systems to fail.
  Removed these functions for non-Windows installations.

* The TofDaq API libraries are not included any more in the source package, but 
  downloaded from https://soft.tofwerk.com/ during installation.
  
* Added a binary version for Mac.

* Updated the source installation instructions to account for the changes in 
  devtools >= 2.0.0.
  
* Internal changes for simpler source installations on Mac and Linux.

* Using Travis CI.

* Moved the links to the API documentation from a vignette to a help file.

* Simplified the cross-platform maintenance: the TofDaqDll functions are always
  exported, but they just give an error message on non-Windows platforms.

* Turned off staged installation (a feature of R >= 3.6.0) to be able to
  dynamically link to the TofDaq libraries.
  
* Updated help pages of `GetIsotopePattern()` and `GetIsotopePattern2()`.

* `GetStringAttributeFromH5()` can also read strings which are longer than 256
  characters.


# Changes in version 0.3.8

* `GetEventListSpectrumFromH5()` returns a zero length vector instead of an
  error if there is no event in the spectrum.

* `DecodeEventList()` also works for events with more than 256 data samples.

* Updated to TofDaq API 1.99r759.

* Added TofDaqDll functions `ConfigVarNbrMemories()`, `ContinueAcquisition()`,
  `ManualContinueNeeded()`, `IssueDio4Pulse()`, `SetDio4State()`,
  `GetSumSpectrumFromShMem2()`, `GetTofSpectrumFromShMem2()`,
  `GetStickSpectrumFromShMem2()`, `GetSegmentProfileFromShMem2()`,
  `GetMassCalib()`, `GetMassCalib2()`, `GetMassCalibEx()`, `GetMassCalib2Ex()`,
  `SetMassCalib()`, `SetMassCalib2()`, `SetMassCalibEx()`, `SetMassCalib2Ex()`.

* Added TwToolDll functions `MultiPeakFit()`, `EvalMultiPeak()`,
  `GetMassCalibInfo()`, `FitResolution()`, `EvalResolution()`.

* Added TwH5Dll functions `H5AddUserDataMultiRow()`, `DeleteAttributeInH5()`,
  `WaitForExclusiveFileAccess()`, `WriteNetCdfTimeSeriesFile()`,
  `H5SetMassCalib()`, `H5SetMassCalib2()`, `H5SetMassCalibEx()`,
  `H5SetMassCalib2Ex()`, `H5SetMassCalibDynamic()`, `ChangePeakTable()`,
  `ChangePeakTableFromFile()`.

* Corrected bug in functions which have an array of strings as an input
  parameter.
  
* Using NEWS.md format for this NEWS file.


# Changes in version 0.3.7

* Corrected bug in `DecodeEventList()`.

* Added Authors@R, URL and BugReports fields in DESCRIPTION, removed Date field.

* Corrected and improved documentation.

* Added vignette that links to the API documentation files.

* Removed the macos_i386 TofDaq libraries.


# Changes in version 0.3.6

* Corrected cstring to std::string conversion bug in `GetRegUserData...()`
  functions.

* Corrected off-by-one bug in `GetPeakParametersFromH5()`.

* Corrected bug in `configure` script (Mac and Linux version).


# Changes in version 0.3.5

* Updated `configure` script to work with Rcpp >=0.12.12.

* Updated README.md with easier installation instructions.

* Updated to TofDaq API 1.99r601.

* Corrected bug in `GetUserDataFromH5()` and added the option to read all rows.


# Changes in version 0.3.4

* Corrected bug in `configure.win` script which made source package installation
  on 32-bit Windows to fail.


# Changes in version 0.3.3

* Corrected bugs in `LoadIniFile()` and `SaveIniFile()`.

* Corrected bug in `configure` script (Mac and Linux version).

* Added wrapper functions `GetTofDataSinglePeak()` and `FitTofDataSinglePeak()`.


# Changes in version 0.3.2

* Corrected bug in internal std::string to char* conversion (introduced in
  TofDaqR 0.3.0).

* Using `Rcpp::stop()` to signal errors if TwRetVal != TwSuccess.

* Corrected bug in `CleanupDll()`.

* Updated some help pages.


# Changes in version 0.3.1

* `GetDescriptor()` returns int64 and uint64 values as strings.

* `AddLogEntry()` also accepts the `LogEntryTime` parameter.

* Added additional API functions (all `GetDaqParameter...()`,
  `SetDaqParameter...()` and `Set...AttributeInH5()` functions,
  `H5AddLogEntry()`, `H5GetMassCalibPar()`, `GetBufWriteProfileFromH5()`)


# Changes in version 0.3.0

* Another major code refactoring: Rcpp is used for all C/C++ code.

* Added all `Get...2FromH5()` functions.


# Changes in version 0.2.0

* Updated to TofDaq API 1.99r450.

* NAMESPACE generated by roxygen2::roxygenise().

* Major refactoring of the package: TofDaqWrapR.dll is not any more compiled
  with Visual Studio. Instead, the C++ source files of the C++/R wrapper are
  included in the package source and compiled with the GCC compiler. This makes
  package maintenance easier and allows the same workflow for all platforms.

* TofDaqR can also be installed on Linux and Mac (for using the twtool and
  twh5 functions)!


# Changes in version 0.1.7

* Removed dependency on Rcpp package by moving `DecodeEventList()` and
  `EventList2TofSpec()` to TofDaqWrapR.dll.


# Changes in version 0.1.6

* Corrected bug when package was unloaded.

* Some general code improvements in TofDaqWrapR.dll.

* `DecodeEventList()` now returns sample indices instead of time stamps.

* Updated some help pages.

* Published TofDaqR to GitHub: https://github.com/pasturm/TofDaqR


# Changes in version 0.1.5

* Updated to TofDaq API 1.98.

* Added `GetRegUserDataSourcesFromH5()`.


# Changes in version 0.1.4

* Extended `Tof2Mass()` and `Mass2Tof()` to work with vectors as first argument.

* Added `InitializeDll()`.

* Updated some help pages.

* Corrected minor bug in `CleanupDll()` and `SetTimeout()` (no return value
  instead of NULL).

* Added undocumented option of histogram-based noise determination in
  `SiSetProcessingOptions()`.

* Added `ReleaseSharedMemory()` and `KeepSharedMemMapped()`.

* Added additional TPS1 codes in `Tps1rc`.


# Changes in version 0.1.3

* Added all single ion histogramming functions (`Si...`) of TwToolDll.
  Additionally, `SiProcessSpectrumFromShMem()` takes the spectrum directly from
  shared memory.

* Updated some help pages.

* Corrected bug in `GetPeakParameters()` (label was always empty).

* Added `SaveMassTableToFile()` function, which saves all peak parameters to a
  mass table file.


# Changes in version 0.1.2

* Renamed all functions, i.e. the `Tw` (and `Tx`) has been removed. (This is a
  breaking change, but because this package has currently no users (except the
  author), it doesn't matter ;-).

* Using NEWS (so that `news(package = "TofDaqR")` works) instead of markdown
  format NEWS/NEWS.md.


# Changes in version 0.1.1

* Made DESCRIPTION to work in RStudio (description field is longer).

* Corrected some bugs in documentation.

* Does not unload TofDaqWrapR.dll when package is unloaded. This prevents some R
  crashes when package is unloaded.

* Added `TwMassCalibrate()`.

* Added this NEWS file.

* Updated `TwLoadIniFile()` and `TwSaveIniFile()` to accept an empty string as
  an argument (-> "TwApiTmpIni.ini").

* Corrected bug in `TwGetDescriptor()` (character parameters were not read
  correctly).

* Added `TxTps1rc`, a data frame giving the codes and names of some TPS modules.

* Added `TwGetRegUserDataSources()` and `TwGetRegUserDataDesc()`.

* Corrected the reading of the data description in `TwGetRegUserDataFromH5()`
  and `TwGetUserDataFromH5()`.


# Changes in version 0.1.0

* First release of TofDaqR as a package.

* TofDaq API 1.97.
