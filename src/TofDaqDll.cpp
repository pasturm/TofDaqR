#ifdef _WIN32

#include <Rcpp.h>
using namespace Rcpp;
#include <TofDaqDll.h>
#include "TofDaqR.h"

// InitializeDll ---------------------------------------------------------------
//' Initializes the TofDaqDll.dll.
//'
//' \code{InitializeDll} initializes the TofDaqDll.dll. It is usually not necessary
//' to call \code{InitializeDll} explicitely, as it is called automatically by
//' functions that need the DLL to be in an initialized state.
//' @export
// [[Rcpp::export]]
SEXP InitializeDll() {

  TwRetVal rv = TwInitializeDll();

  return TwRetValString(rv);
}

// CleanupDll ------------------------------------------------------------------
//' Deinitializes the TofDaqDll.dll.
//'
//' \code{CleanupDll} deinitializes the TofDaqDll.dll (frees allocated memory,
//' releases mapped shared memory and closes open files). This function is
//' automatically called when the TofDaqR package is unloaded.
//' @export
// [[Rcpp::export]]
void CleanupDll() {

  void TwCleanupDll(void);
}

// GetDllVersion ---------------------------------------------------------------
//' Gets the version number of the TofDaq API.
//'
//' \code{GetDllVersion} gets the version number of the TofDaq API.
//' @export
// [[Rcpp::export]]
double GetDllVersion() {

  double rv = TwGetDllVersion();
  return rv;
}

// TofDaqRunning ---------------------------------------------------------------
//' Checks if TofDaq recorder application is running.
//'
//' \code{TofDaqRunning} checks if TofDaq recorder application is running.
//'
//' @return \code{TRUE} or \code{FALSE}.
//' @export
// [[Rcpp::export]]
bool TofDaqRunning() {

  bool rv = TwTofDaqRunning();

  return rv;
}

// DaqActive -------------------------------------------------------------------
//' Checks if TofDaq recorder is currently acquiring data.
//'
//' \code{DaqActive} checks if TofDaq recorder is currently acquiring data.
//'
//' @return \code{TRUE} or \code{FALSE}.
//' @export
// [[Rcpp::export]]
bool DaqActive() {

  bool rv = TwDaqActive();

  return rv;
}

// StartAcquisition ------------------------------------------------------------
//' Starts an acquisition.
//'
//' \code{StartAcquisition} starts an acquisition.
//' @export
// [[Rcpp::export]]
SEXP StartAcquisition() {

  TwRetVal rv = TwStartAcquisition();

  return TwRetValString(rv);
}

// StopAcquisition -------------------------------------------------------------
//' Stops the current acquisition.
//'
//' \code{StopAcquisition} stops the current acquisition.
//' @export
// [[Rcpp::export]]
SEXP StopAcquisition() {

  TwRetVal rv = TwStopAcquisition();

  return TwRetValString(rv);
}

// CloseTofDaqRec --------------------------------------------------------------
//' Closes the TofDaq recorder application.
//'
//' \code{CloseTofDaqRec} closes the TofDaq recorder application.
//' @export
// [[Rcpp::export]]
SEXP CloseTofDaqRec() {

  TwRetVal rv = TwCloseTofDaqRec();

  return TwRetValString(rv);
}

// InitializeDaqDevice ---------------------------------------------------------
//' Initializes the DAQ board.
//'
//' \code{InitializeDaqDevice} initializes the DAQ board (this is also done at
//' startup of TofDaqRec.exe). This can take up to 8 seconds depending on the
//' actual DAQ hardware.
//' @export
// [[Rcpp::export]]
SEXP InitializeDaqDevice() {

  TwRetVal rv = TwInitializeDaqDevice();

  return TwRetValString(rv);
}

// SetTimeout ------------------------------------------------------------------
//' Sets the timeout.
//'
//' \code{SetTimeout} sets the global timeout for all functions that can time
//' out. Default is 500 ms.
//'
//' @param timeout Timeout in ms. Default is 500 ms.
//' @export
// [[Rcpp::export]]
void SetTimeout(int timeout) {

  TwSetTimeout(timeout);
}

// GetTimeout ------------------------------------------------------------------
//' Gets the timeout.
//'
//' \code{GetTimeout} gets the current timeout value (in ms).
//' @export
// [[Rcpp::export]]
int GetTimeout() {

  int rv = TwGetTimeout();

  return rv;
}

// AutoSetupDaqDevice ----------------------------------------------------------
//' Auto setup routine for the DAQ device.
//'
//' Auto setup routine for the DAQ device. Currently implemented only for
//' AP240 averager and ndigo5G.
//' @export
// [[Rcpp::export]]
SEXP AutoSetupDaqDevice() {

  TwRetVal rv = TwAutoSetupDaqDevice();

  return TwRetValString(rv);
}

// OnDemandMassCalibration -----------------------------------------------------
//' Arms/executes on demand mass calibration.
//'
//' \code{OnDemandMassCalibration} arms/executes on demand mass calibration.
//' Requires that parameter \code{ReCalibFreq} is set to 3 (on demand) in order to work.
//'
//' @param action 0: arms mass calibration, 1: updates mass calibration using
//' data acquired since previous arm.
//' @export
// [[Rcpp::export]]
SEXP OnDemandMassCalibration(int action) {

  TwRetVal rv = TwOnDemandMassCalibration(action);

  return TwRetValString(rv);
}

// ShowConfigWindow ------------------------------------------------------------
//' Shows the TofDaq recorder configuration windows.
//'
//' \code{ShowConfigWindow} shows the different tabs of the TofDaq recorder
//' configuration window.
//'
//' @param ConfigWindowIndex Index of configuration tab to show (valid range: 0-6)
//' @export
// [[Rcpp::export]]
SEXP ShowConfigWindow(int ConfigWindowIndex) {

  TwRetVal rv = TwShowConfigWindow(ConfigWindowIndex);

  return TwRetValString(rv);
}

// LoadIniFile -----------------------------------------------------------------
//' Loads a configuration file.
//'
//' \code{LoadIniFile} loads a configuration file (*.ini) from disk.
//'
//' @param IniFile Path/filename of the configuration file. If \code{IniFile} is
//' an empty string or NULL, "TwApiTmpIni.ini" in the TofDaq recorder directory
//' will be used.
//' @export
// [[Rcpp::export]]
SEXP LoadIniFile(Nullable<Rcpp::String> IniFile = R_NilValue) {

  std::string tmp = Rcpp::as<std::string>(IniFile);
  char *cFilename;
  if (IniFile.isNotNull() && tmp.empty()) {
    cFilename = RtoCstring(IniFile);
  } else {
    cFilename = NULL;
  }
  TwRetVal rv = TwLoadIniFile(cFilename);

  return TwRetValString(rv);
}

// SaveIniFile -----------------------------------------------------------------
//' Saves the current configuration (*.ini) to disk.
//'
//' \code{SaveIniFile} saves the current configuration (*.ini) to disk.
//'
//' @param IniFile Path/filename of the configuration file. If \code{IniFile} is
//' an empty string or NULL, "TwApiTmpIni.ini" in the TofDaq recorder directory
//' will be used.
//' @export
// [[Rcpp::export]]
SEXP SaveIniFile(Nullable<Rcpp::String> IniFile = R_NilValue) {

  std::string tmp = Rcpp::as<std::string>(IniFile);
  char *cFilename;
  if (IniFile.isNotNull() && tmp.empty()) {
    cFilename = RtoCstring(IniFile);
  } else {
    cFilename = NULL;
  }
  TwRetVal rv = TwLoadIniFile(cFilename);

  return TwRetValString(rv);
}

// GetDaqParameter -------------------------------------------------------------
//' Gets a single parameter as a string.
//'
//' \code{GetDaqParameter} gets a single parameter as a string.
//'
//' @param Parameter Parameter name as a string. See
//' \emph{/doc/TofDaqDll.htm#parameter_list} for a list of all available parameters.
//' @export
// [[Rcpp::export]]
String GetDaqParameter(SEXP Parameter) {

  char *cParameter = RtoCstring(Parameter);

  char *buffer = new char[256];

  buffer = TwGetDaqParameter(cParameter);

  std::string str(buffer);

  return str;
}

// GetDaqParameterInt ----------------------------------------------------------
//' Gets a single integer parameter.
//'
//' \code{GetDaqParameterInt} gets a single integer parameter.
//'
//' @param Parameter Parameter name as a string. See
//' \emph{/doc/TofDaqDll.htm#parameter_list} for a list of all available parameters.
//' @export
// [[Rcpp::export]]
int GetDaqParameterInt(SEXP Parameter) {

  char *cParameter = RtoCstring(Parameter);

  int result;

  result = TwGetDaqParameterInt(cParameter);

  return result;
}

// GetDaqParameterBool ---------------------------------------------------------
//' Gets a single boolean parameter.
//'
//' \code{GetDaqParameterBool} gets a single boolean parameter.
//'
//' @param Parameter Parameter name as a string. See
//' \emph{/doc/TofDaqDll.htm#parameter_list} for a list of all available parameters.
//' @export
// [[Rcpp::export]]
bool GetDaqParameterBool(SEXP Parameter) {

  char *cParameter = RtoCstring(Parameter);

  bool result;

  result = TwGetDaqParameterBool(cParameter);

  return result;
}

// GetDaqParameterFloat --------------------------------------------------------
//' Gets a single float parameter.
//'
//' \code{GetDaqParameterFloat} gets a single float parameter.
//'
//' @param Parameter Parameter name as a string. See
//' \emph{/doc/TofDaqDll.htm#parameter_list} for a list of all available parameters.
//' @export
// [[Rcpp::export]]
double GetDaqParameterFloat(SEXP Parameter) {

  char *cParameter = RtoCstring(Parameter);

  float result;

  result = TwGetDaqParameterFloat(cParameter);

  return (double)result;
}

// GetDaqParameterInt64 --------------------------------------------------------
//' Gets a single int64 parameter as a string.
//'
//' \code{GetDaqParameterInt64} gets a single int64 parameter as a string.
//'
//' The return string can be converted to integer64 using \code{\link[bit64]{bit64::as.integer64}}.
//'
//' @param Parameter Parameter name as a string. See
//' \emph{/doc/TofDaqDll.htm#parameter_list} for a list of all available parameters.
//' @export
// [[Rcpp::export]]
String GetDaqParameterInt64(SEXP Parameter) {

  char *cParameter = RtoCstring(Parameter);

  __int64 result;

  result = TwGetDaqParameterInt64(cParameter);

  std::stringstream ss;
  ss << result;

  return ss.str();
}

// GetDaqParameterDouble -------------------------------------------------------
//' Gets a single double parameter.
//'
//' \code{GetDaqParameterDouble} gets a single double parameter.
//'
//' @param Parameter Parameter name as a string. See
//' \emph{/doc/TofDaqDll.htm#parameter_list} for a list of all available parameters.
//' @export
// [[Rcpp::export]]
double GetDaqParameterDouble(SEXP Parameter) {

  char *cParameter = RtoCstring(Parameter);

  double result;

  result = TwGetDaqParameterDouble(cParameter);

  return result;
}

// GetDaqParameterIntRef -------------------------------------------------------
//' Gets a single integer parameter.
//'
//' \code{GetDaqParameterIntRef} gets a single integer parameter.
//'
//' This is the same as \code{GetDaqParameterInt}, but additionally it checks
//' for success and a TwRetVal string is returned if it is not sucessful.
//'
//' @param Parameter Parameter name as a string. See
//' \emph{/doc/TofDaqDll.htm#parameter_list} for a list of all available parameters.
//' @export
// [[Rcpp::export]]
SEXP GetDaqParameterIntRef(SEXP Parameter) {

  char *cParameter = RtoCstring(Parameter);

  int Value;

  TwRetVal rv = TwGetDaqParameterIntRef(cParameter, &Value);
  if (rv != TwSuccess) {
    return TwRetValString(rv);
  }

  return wrap(Value);
}

// GetDaqParameterBoolRef ------------------------------------------------------
//' Gets a single boolean parameter.
//'
//' \code{GetDaqParameterBoolRef} gets a single boolean parameter.
//'
//' This is the same as \code{GetDaqParameterBool}, but additionally it checks
//' for success and a TwRetVal string is returned if it is not sucessful.
//'
//' @param Parameter Parameter name as a string. See
//' \emph{/doc/TofDaqDll.htm#parameter_list} for a list of all available parameters.
//' @export
// [[Rcpp::export]]
SEXP GetDaqParameterBoolRef(SEXP Parameter) {

  char *cParameter = RtoCstring(Parameter);

  bool Value;

  TwRetVal rv = TwGetDaqParameterBoolRef(cParameter, &Value);
  if (rv != TwSuccess) {
    return TwRetValString(rv);
  }

  return wrap(Value);
}

// GetDaqParameterFloatRef -----------------------------------------------------
//' Gets a single float parameter.
//'
//' \code{GetDaqParameterFloatRef} gets a single float parameter.
//'
//' This is the same as \code{GetDaqParameterFloat}, but additionally it checks
//' for success and a TwRetVal string is returned if it is not sucessful.
//'
//' @param Parameter Parameter name as a string. See
//' \emph{/doc/TofDaqDll.htm#parameter_list} for a list of all available parameters.
//' @export
// [[Rcpp::export]]
SEXP GetDaqParameterFloatRef(SEXP Parameter) {

  char *cParameter = RtoCstring(Parameter);

  float Value;

  TwRetVal rv = TwGetDaqParameterFloatRef(cParameter, &Value);
  if (rv != TwSuccess) {
    return TwRetValString(rv);
  }

  return wrap(Value);
}

// GetDaqParameterInt64Ref -----------------------------------------------------
//' Gets a single int64 parameter as a string.
//'
//' \code{GetDaqParameterInt64Ref} gets a single int64 parameter as a string.
//'
//' This is the same as \code{GetDaqParameterInt64}, but additionally it checks
//' for success and a TwRetVal string is returned if it is not sucessful.
//'
//' @param Parameter Parameter name as a string. See
//' \emph{/doc/TofDaqDll.htm#parameter_list} for a list of all available parameters.
//' @export
// [[Rcpp::export]]
SEXP GetDaqParameterInt64Ref(SEXP Parameter) {

  char *cParameter = RtoCstring(Parameter);

  __int64 Value;

  TwRetVal rv = TwGetDaqParameterInt64Ref(cParameter, &Value);
  if (rv != TwSuccess) {
    return TwRetValString(rv);
  }

  std::stringstream ss;
  ss << Value;

  return wrap(ss.str());
}

// GetDaqParameterDoubleRef ----------------------------------------------------
//' Gets a single double parameter.
//'
//' \code{GetDaqParameterDoubleRef} gets a single double parameter.
//'
//' This is the same as \code{GetDaqParameterDouble}, but additionally it checks
//' for success and a TwRetVal string is returned if it is not sucessful.
//'
//' @param Parameter Parameter name as a string. See
//' \emph{/doc/TofDaqDll.htm#parameter_list} for a list of all available parameters.
//' @export
// [[Rcpp::export]]
SEXP GetDaqParameterDoubleRef(SEXP Parameter) {

  char *cParameter = RtoCstring(Parameter);

  double Value;

  TwRetVal rv = TwGetDaqParameterDoubleRef(cParameter, &Value);
  if (rv != TwSuccess) {
    return TwRetValString(rv);
  }

  return wrap(Value);
}

// GetDaqParameterStringRef -----------------------------------------------------
//' Gets a single string parameter.
//'
//' \code{GetDaqParameterStringRef} gets a single string parameter.
//'
//' This is the same as \code{GetDaqParameter}, but returns \code{"TwInvalidValue"}
//' if the type of the parameter is not a string.
//'
//' @param Parameter Parameter name as a string. See
//' \emph{/doc/TofDaqDll.htm#parameter_list} for a list of all available parameters.
//' @export
// [[Rcpp::export]]
SEXP GetDaqParameterStringRef(SEXP Parameter) {

  char *cParameter = RtoCstring(Parameter);

  char *Value = new char[256];

  TwRetVal rv = TwGetDaqParameterStringRef(cParameter, Value);
  if (rv != TwSuccess) {
    return TwRetValString(rv);
  }

  std::string str(Value);

  return wrap(str);
}

// SetDaqParameter -------------------------------------------------------------
//' Sets a single parameter.
//'
//' \code{SetDaqParameter} sets a single parameter.
//'
//' @param Parameter Parameter name as a string. See
//' \emph{/doc/TofDaqDll.htm#parameter_list} for a list of all available parameters.
//' @param ValueString Value as a string.
//' @export
// [[Rcpp::export]]
SEXP SetDaqParameter(SEXP Parameter, SEXP ValueString) {

  char *cParameter = RtoCstring(Parameter);
  char *cValueString = RtoCstring(ValueString);

  TwRetVal rv = TwSetDaqParameter(cParameter, cValueString);

  return TwRetValString(rv);
}

// SetDaqParameterInt ----------------------------------------------------------
//' Sets a single parameter with an integer value.
//'
//' \code{SetDaqParameterInt} sets a single parameter with an integer value.
//'
//' @param Parameter Parameter name as a string. See
//' \emph{/doc/TofDaqDll.htm#parameter_list} for a list of all available parameters.
//' @param Value Integer value.
//' @export
// [[Rcpp::export]]
SEXP SetDaqParameterInt(SEXP Parameter, int Value) {

  char *cParameter = RtoCstring(Parameter);

  TwRetVal rv = TwSetDaqParameterInt(cParameter, Value);

  return TwRetValString(rv);
}

// SetDaqParameterBool ---------------------------------------------------------
//' Sets a single parameter with a boolean value.
//'
//' \code{SetDaqParameterBool} sets a single parameter with a boolean value.
//'
//' @param Parameter Parameter name as a string. See
//' \emph{/doc/TofDaqDll.htm#parameter_list} for a list of all available parameters.
//' @param Value \code{TRUE} or \code{FALSE}.
//' @export
// [[Rcpp::export]]
SEXP SetDaqParameterBool(SEXP Parameter, bool Value) {

  char *cParameter = RtoCstring(Parameter);

  TwRetVal rv = TwSetDaqParameterBool(cParameter, Value);

  return TwRetValString(rv);
}

// SetDaqParameterFloat --------------------------------------------------------
//' Sets a single parameter with a float value.
//'
//' \code{SetDaqParameterFloat} sets a single parameter with a float value.
//'
//' @param Parameter Parameter name as a string. See
//' \emph{/doc/TofDaqDll.htm#parameter_list} for a list of all available parameters.
//' @param Value Numeric value.
//' @export
// [[Rcpp::export]]
SEXP SetDaqParameterFloat(SEXP Parameter, double Value) {

  char *cParameter = RtoCstring(Parameter);

  TwRetVal rv = TwSetDaqParameterFloat(cParameter, (float)Value);

  return TwRetValString(rv);
}

// SetDaqParameterInt64 ----------------------------------------------------------
//' Sets a single parameter with an int64 value.
//'
//' \code{SetDaqParameterInt64} sets a single parameter with an int64 value
//' (passed as a string).
//'
//' @param Parameter Parameter name as a string. See
//' \emph{/doc/TofDaqDll.htm#parameter_list} for a list of all available parameters.
//' @param Value int64 value passed as a string.
//' @export
// [[Rcpp::export]]
SEXP SetDaqParameterInt64(SEXP Parameter, SEXP Value) {

  char *cParameter = RtoCstring(Parameter);

  std::string str = Rcpp::as<std::string>(Value);
  std::stringstream ss(str);
  __int64 int64value;
  ss >> int64value;

  TwRetVal rv = TwSetDaqParameterInt64(cParameter, int64value);

  return TwRetValString(rv);
}

// SetDaqParameterDouble -------------------------------------------------------
//' Sets a single parameter with a double value.
//'
//' \code{SetDaqParameterDouble} sets a single parameter with a double value.
//'
//' @param Parameter Parameter name as a string. See
//' \emph{/doc/TofDaqDll.htm#parameter_list} for a list of all available parameters.
//' @param Value Numeric value.
//' @export
// [[Rcpp::export]]
SEXP SetDaqParameterDouble(SEXP Parameter, double Value) {

  char *cParameter = RtoCstring(Parameter);

  TwRetVal rv = TwSetDaqParameterDouble(cParameter, Value);

  return TwRetValString(rv);
}

// GetDescriptor --------------------------------------------------------------
//' Gets various information about the active acquisition.
//'
//' \code{GetDescriptor} retrieves the current TSharedMemoryDesc structure.
//' TSharedMemoryDesc contains various static information about the active
//' acquisition that can be retrieved by \code{GetDaqParameter} functions but
//' also information of DAQ progress.
//' See \emph{/doc/TofDaqDll.htm} for more details.
//'
//' int64 and unsigned int64 parameters are returned as string. They can be
//' converted to integer64 using \code{\link[bit64]{bit64::as.integer64}}.
//'
//' @return A list containing the TSharedMemoryDesc structure
//' @export
// [[Rcpp::export]]
SEXP GetDescriptor() {

  //get the current TSharedMemoryDesc structure.
  TSharedMemoryDesc pBufDesc;
  TwRetVal rv = TwGetDescriptor(&pBufDesc);
  if (rv != TwSuccess) {
    return TwRetValString(rv);
  }

  std::stringstream ss;  // for converting int64 to string

  List result;
  result["NbrSamples"] = pBufDesc.NbrSamples;
  result["NbrRawSamples"] = pBufDesc.NbrRawSamples;
  result["NbrPeaks"] = pBufDesc.NbrPeaks;
  result["NbrWaveforms"] = pBufDesc.NbrWaveforms;
  result["NbrSegments"] = pBufDesc.NbrSegments;
  result["NbrBlocks"] = pBufDesc.NbrBlocks;
  result["NbrMemories"] = pBufDesc.NbrMemories;
  result["NbrBufs"] = pBufDesc.NbrBufs;
  result["NbrWrites"] = pBufDesc.NbrWrites;
  result["NbrRuns"] = pBufDesc.NbrRuns;
  result["iWaveform"] = pBufDesc.iWaveform;
  result["iSegment"] = pBufDesc.iSegment;
  result["iBlock"] = pBufDesc.iBlock;
  result["iMemory"] = pBufDesc.iMemory;
  result["iBuf"] = pBufDesc.iBuf;
  result["iWrite"] = pBufDesc.iWrite;
  result["iRun"] = pBufDesc.iRun;
  result["TotalBufsRecorded"] = pBufDesc.TotalBufsRecorded;
  result["TotalBufsProcessed"] = pBufDesc.TotalBufsProcessed;
  result["TotalBufsWritten"] = pBufDesc.TotalBufsWritten;
  result["OverallBufsProcessed"] = pBufDesc.OverallBufsProcessed;
  result["TotalNbrMemories"] = pBufDesc.TotalNbrMemories;
  result["TotalMemoriesProcessed"] = pBufDesc.TotalMemoriesProcessed;
  result["RawDataRecordedBuf1"] = pBufDesc.RawDataRecordedBuf1;
  result["RawDataRecordedBuf2"] = pBufDesc.RawDataRecordedBuf2;
  result["RawDataLastElementInBuffer1"] = pBufDesc.RawDataLastElementInBuffer1;
  result["RawDataLastElementInBuffer2"] = pBufDesc.RawDataLastElementInBuffer2;
  result["RawDataProcessedBuf1"] = pBufDesc.RawDataProcessedBuf1;
  result["RawDataProcessedBuf2"] = pBufDesc.RawDataProcessedBuf2;
  result["RawDataWrittenBuf1"] = pBufDesc.RawDataWrittenBuf1;
  result["RawDataWrittenBuf2"] = pBufDesc.RawDataWrittenBuf2;
  result["SampleInterval"] = pBufDesc.SampleInterval;
  result["TofPeriod"] = pBufDesc.TofPeriod;
  result["NbrCubes"] = pBufDesc.NbrCubes;
  ss << pBufDesc.BlockPeriod;
  result["BlockPeriod"] = ss.str();
  ss.str(std::string());
  ss << pBufDesc.BlockPulseDelay;
  result["BlockPulseDelay"] = ss.str();
  ss.str(std::string());
  ss << pBufDesc.BlockDelay;
  result["BlockDelay"] = ss.str();
  result["SingleIonSignal"] = pBufDesc.SingleIonSignal;
  result["SingleIonSignal2"] = pBufDesc.SingleIonSignal2;
  result["MassCalibMode"] = pBufDesc.MassCalibMode;
  result["MassCalibMode2"] = pBufDesc.MassCalibMode2;
  result["NbrMassCalibParams"] = pBufDesc.NbrMassCalibParams;
  result["NbrMassCalibParams2"] = pBufDesc.NbrMassCalibParams2;
  NumericVector p(pBufDesc.p, pBufDesc.p + 16);
  result["p"] = p;
  NumericVector p2(pBufDesc.p2, pBufDesc.p2 + 16);
  result["p2"] = p2;
  result["R0"] = pBufDesc.R0;
  result["dm"] = pBufDesc.dm;
  result["m0"] = pBufDesc.m0;
  result["SecondTof"] = pBufDesc.SecondTof;
  result["chIniFileName"] = pBufDesc.chIniFileName;
  result["CurrentDataFileName"] = pBufDesc.CurrentDataFileName;
  result["segIlf"] = pBufDesc.segIlf;
  result["iCube"] = pBufDesc.iCube;
  result["DaqMode"] = pBufDesc.DaqMode;
  result["AcquisitionMode"] = pBufDesc.AcquisitionMode;
  result["CombineMode"] = pBufDesc.CombineMode;
  result["RecalibFreq"] = pBufDesc.RecalibFreq;
  result["AcquisitionLogText"] = pBufDesc.AcquisitionLogText;
  ss.str(std::string());
  ss << pBufDesc.AcquisitionLogTime;
  result["AcquisitionLogTime"] = ss.str();
  ss.str(std::string());
  ss << pBufDesc.TimeZero;
  result["TimeZero"] = ss.str();
  result["ExternalLock"] = pBufDesc.ExternalLock;
  result["ProcessingLevel"] = pBufDesc.ProcessingLevel;
  result["AttributeType"] = pBufDesc.AttributeType;
  result["AttributeObject"] = pBufDesc.AttributeObject;
  result["AttributeName"] = pBufDesc.AttributeName;
  result["AttributeInt"] = pBufDesc.AttributeInt;
  result["AttributeDouble"] = pBufDesc.AttributeDouble;
  result["EnableVarNbrMemories"] = pBufDesc.EnableVarNbrMemories;
  result["NbrSteps"] = pBufDesc.NbrSteps;
  result["CurrentStepAtBuf"] = pBufDesc.CurrentStepAtBuf;
  result["NbrMemoriesForCurrentStep"] = pBufDesc.NbrMemoriesForCurrentStep;

  return result;
}

// GetPeakParameters -----------------------------------------------------------
//' Gets parameters for a given peak.
//'
//' \code{GetPeakParameters} gets parameters for a given peak.
//'
//' @param PeakIndex Index of peak.
//' @return A list with the peak paramters \emph{label}, \emph{mass}, \emph{loMass} and \emph{hiMass}.
//' @export
// [[Rcpp::export]]
SEXP GetPeakParameters(int PeakIndex) {

  TPeakPar PeakPar;
  TwRetVal rv = TwGetPeakParameters(&PeakPar, PeakIndex);

  if (rv != TwSuccess) {
    return TwRetValString(rv);
  }

  List result;
  result["label"] = PeakPar.label;
  result["mass"] = PeakPar.mass;
  result["loMass"] = PeakPar.loMass;
  result["hiMass"] = PeakPar.hiMass;

  return result;
}

// ReleaseSharedMemory ---------------------------------------------------------
//' Manually releases the shared memory acquisition buffers.
//'
//' \code{ReleaseSharedMemory} manually releases the shared memory acquisition
//' buffers. This is needed if \code{\link{KeepSharedMemMapped}} has been set to
//' \code{TRUE}.
//' @export
// [[Rcpp::export]]
SEXP ReleaseSharedMemory() {

  TwRetVal rv = TwReleaseSharedMemory();

  return TwRetValString(rv);
}

// WaitForNewData --------------------------------------------------------------
//' Waits for new data.
//'
//' \code{WaitForNewData} waits for new data. Returns when new data is
//' available or when timed out.
//'
//' @param timeout Timeout in ms.
//' @param WaitForEventReset If \code{TRUE} (default) waits for application to
//' reset data available event before returning.
//' @export
// [[Rcpp::export]]
SEXP WaitForNewData(int timeout, bool WaitForEventReset) {

  TSharedMemoryDesc pBufDesc;
  TSharedMemoryPointer pShMem;

  TwRetVal rv = TwWaitForNewData(timeout, &pBufDesc, &pShMem, WaitForEventReset);

  return TwRetValString(rv);
}

// WaitForEndOfAcquisition -----------------------------------------------------
//' Waits for the end of the current acquisition.
//'
//' \code{WaitForEndOfAcquisition} waits for the end of the current acquisition.
//' If \code{NbrRuns > 1} this function waits for the end of the last acquisition.
//'
//' @param timeout Timeout in ms.
//' @export
// [[Rcpp::export]]
SEXP WaitForEndOfAcquisition(int timeout) {

  TwRetVal rv = TwWaitForEndOfAcquisition(timeout);

  return TwRetValString(rv);
}

// GetSumSpectrumFromShMem -----------------------------------------------------
//' Sum spectrum from shared memory.
//'
//' \code{GetSumSpectrumFromShMem} gets the sum spectrum from shared memory.
//'
//' @param Normalize If \code{FALSE} the spectrum is reported as sum,
//' if \code{TRUE} (default) the spectrum is normalized to counts per extraction.
//' @export
// [[Rcpp::export]]
SEXP GetSumSpectrumFromShMem(bool Normalize) {

  TSharedMemoryDesc desc;
  TwRetVal rv = TwGetDescriptor(&desc);
  if (rv != TwSuccess) {
    return TwRetValString(rv);
  }

  NumericVector Spectrum(desc.NbrSamples);

  rv = TwGetSumSpectrumFromShMem(&Spectrum[0], Normalize);
  if (rv != TwSuccess) {
    return TwRetValString(rv);
  }

  return Spectrum;
}

// GetTofSpectrumFromShMem -----------------------------------------------------
//' Single TOF spectrum from shared memory.
//'
//' \code{GetTofSpectrumFromShMem} reads a single TOF spectrum (possibly
//' averaged/summed over segment dimension) from shared memory. If
//' \code{SegmentIndex = SegmentEndIndex = -1} the complete block of data is
//' copied and the \code{Normalize} flag is ignored.
//'
//' @param SegmentIndex Segment start index of data to fetch (or -1 for complete
//' block copy).
//' @param SegmentEndIndex Segment end index of data to fetch (or -1 for complete
//' block copy).
//' @param BufIndex Buf index of data to fetch.
//' @param Normalize If \code{FALSE} the spectrum is reported as sum,
//' if \code{TRUE} (default) the spectrum is normalized to counts per extraction.
//' @return A vector containing the mass spectrum or an array containing the
//' block of mass spectra if \code{SegmentIndex = SegmentEndIndex = -1}.
//' @export
// [[Rcpp::export]]
SEXP GetTofSpectrumFromShMem(int SegmentIndex, int SegmentEndIndex,
                             int BufIndex, bool Normalize) {

  TSharedMemoryDesc desc;
  TwRetVal rv = TwGetDescriptor(&desc);
  if (rv != TwSuccess) {
    return TwRetValString(rv);
  }

  int specLen;
  if (SegmentIndex == -1 && SegmentEndIndex == -1) {
    specLen = desc.NbrSamples*desc.NbrSegments;
  } else {
    specLen = desc.NbrSamples;
  }

  std::vector<float> Spectrum(specLen);
  rv = TwGetTofSpectrumFromShMem(&Spectrum[0], SegmentIndex, SegmentEndIndex,
                                 BufIndex, Normalize);

  if (rv != TwSuccess) {
    return TwRetValString(rv);
  }

  return wrap(Spectrum);
}

// GetSpecXaxisFromShMem -------------------------------------------------------
//' X-axis values of mass spectrum.
//'
//' \code{GetSpecXaxisFromShMem} returns an array of x-axis values of the mass
//' spectrum.
//'
//' @param Type x-axis type (0: sample index, 1: mass/charge [Th],
//' -1: mass/charge [Th] (2nd TOF), 2: time of flight [microsec],
//' -2: time of flight [microsec] (2nd TOF), 3: frequency [kHz]).
//' @return A vector containing the x-axis values.
//' @export
// [[Rcpp::export]]
SEXP GetSpecXaxisFromShMem(int Type) {

  //get descriptor of file
  TSharedMemoryDesc desc;
  TwRetVal rv = TwGetDescriptor(&desc);
  if (rv != TwSuccess) {
    return TwRetValString(rv);
  }

  NumericVector SpecAxis(desc.NbrSamples);
  double maxMass = 0.0;
  rv = TwGetSpecXaxisFromShMem(&SpecAxis[0], Type, NULL, maxMass);
  if (rv != TwSuccess) {
    return TwRetValString(rv);
  }

  return SpecAxis;
}

// GetStickSpectrumFromShMem ---------------------------------------------------
//' Single stick spectrum from shared memory.
//'
//' \code{GetStickSpectrumFromShMem} reads a single stick spectrum from shared
//' memory.
//'
//' @param SegmentIndex Segment start index of data to fetch.
//' @param SegmentEndIndex Segment end index of data to fetch.
//' @param BufIndex Buf index of data to fetch.
//' @return A list containing the stick spectrum and corresponding masses.
//' @export
// [[Rcpp::export]]
SEXP GetStickSpectrumFromShMem(int SegmentIndex, int SegmentEndIndex,
                               int BufIndex) {

  TSharedMemoryDesc desc;
  TwRetVal rv = TwGetDescriptor(&desc);
  if (rv != TwSuccess) {
    return TwRetValString(rv);
  }

  std::vector<float> Spectrum(desc.NbrPeaks);
  std::vector<float> Masses(desc.NbrPeaks);

  rv = TwGetStickSpectrumFromShMem(&Spectrum[0], &Masses[0], SegmentIndex,
                                   SegmentEndIndex, BufIndex);
  if (rv != TwSuccess) {
    return TwRetValString(rv);
  }

  List result;
  result["Spectrum"] = Spectrum;
  result["Masses"] = Masses;

  return result;
}

// GetSegmentProfileFromShMem --------------------------------------------------
//' Segment profile for a given peak and buf index from shared memory.
//'
//' \code{GetSegmentProfileFromShMem} reads the segment profile for a given
//' peak and buf index from shared memory. Use -1 for \code{PeakIndex} to get
//' segment profiles of all peaks.
//'
//' @param PeakIndex Index of peak to fetch segment profile from. All peaks are
//' read if \code{PeakIndex = -1}.
//' @param BufIndex Buf index of data to fetch.
//' @return A vector containing the segment profile(s).
//' @export
// [[Rcpp::export]]
SEXP GetSegmentProfileFromShMem(int PeakIndex, int BufIndex) {

  TSharedMemoryDesc desc;
  TwRetVal rv = TwGetDescriptor(&desc);
  if (rv != TwSuccess) {
    return TwRetValString(rv);
  }

  int profileLen;
  if (PeakIndex == -1) {
    profileLen = desc.NbrSegments*desc.NbrPeaks;
  } else {
    profileLen = desc.NbrSegments;
  }

  std::vector<float> SegmentProfile(profileLen);
  rv = TwGetSegmentProfileFromShMem(&SegmentProfile[0], PeakIndex, BufIndex);
  if (rv != TwSuccess) {
    return TwRetValString(rv);
  }

  return wrap(SegmentProfile);
}

// GetBufTimeFromShMem ---------------------------------------------------------
//' Time stamp for a given buf and write.
//'
//' \code{GetBufTimeFromShMem} reads the time stamp for a given buf and write
//' from shared memory.
//'
//' @param BufIndex Buf index.
//' @param WriteIndex Write index.
//' @return A time stamp (in seconds relative to acquisition start).
//' @export
// [[Rcpp::export]]
SEXP GetBufTimeFromShMem(int BufIndex, int WriteIndex) {

  NumericVector BufTime(1);

  TwRetVal rv = TwGetBufTimeFromShMem(&BufTime[0], BufIndex, WriteIndex);
  if (rv != TwSuccess) {
    return TwRetValString(rv);
  }

  return BufTime;
}

// Data storage functions ------------------------------------------------------

// AddLogEntry -----------------------------------------------------------------
//' Adds an entry to the acquisition log.
//'
//' \code{AddLogEntry} adds an entry to the acquisition log.
//'
//' @param LogEntryText Log text (max. 255 characters).
//' @param LogEntryTime Log entry time (number of 100-nanosecond intervals since
//' January 1, 1601 UTC) passed as a string. Set it to "0" for "now".
//'
//' @family Data storage functions
//' @export
// [[Rcpp::export]]
SEXP AddLogEntry(SEXP LogEntryText, SEXP LogEntryTime) {

  char *cLogEntryText = RtoCstring(LogEntryText);

  std::string str = Rcpp::as<std::string>(LogEntryTime);
  std::stringstream ss(str);
  unsigned __int64 cTime;
  ss >> cTime;

  TwRetVal rv = TwAddLogEntry(cLogEntryText, cTime);

  return TwRetValString(rv);
}

// AddAttributeInt -------------------------------------------------------------
//' Attaches an integer attribute to the current HDF5 file.
//'
//' \code{AddAttributeInt} attaches an integer attribute to the current HDF5 file.
//'
//' @param Object HDF5 object (group or dataset) to attach attribute (max. 255
//' characters).
//' @param AttributeName Attribute name (max. 127 characters).
//' @param Value Attribute value (integer type).
//'
//' @family Data storage functions
//' @export
// [[Rcpp::export]]
SEXP AddAttributeInt(SEXP Object, SEXP AttributeName, int Value) {

  char *cObject = RtoCstring(Object);
  char *cAttributeName = RtoCstring(AttributeName);

  TwRetVal rv = TwAddAttributeInt(cObject, cAttributeName, Value);

  return TwRetValString(rv);
}

// AddAttributeDouble ----------------------------------------------------------
//' Attaches a numeric attribute to the current HDF5 file.
//'
//' \code{AddAttributeDouble} attaches a numeric attribute to the current HDF5 file.
//'
//' @param Object HDF5 object (group or dataset) to attach attribute (max. 255
//' characters).
//' @param AttributeName Attribute name (max. 127 characters).
//' @param Value Attribute value (numeric type).
//'
//' @family Data storage functions
//' @export
// [[Rcpp::export]]
SEXP AddAttributeDouble(SEXP Object, SEXP AttributeName, double Value) {

  char *cObject = RtoCstring(Object);
  char *cAttributeName = RtoCstring(AttributeName);

  TwRetVal rv = TwAddAttributeDouble(cObject, cAttributeName, Value);

  return TwRetValString(rv);
}

// AddAttributeString ----------------------------------------------------------
//' Attaches a string attribute to the current HDF5 file.
//'
//' \code{AddAttributeString} attaches a string attribute to the current HDF5 file.
//'
//' @param Object HDF5 object (group or dataset) to attach attribute (max. 255
//' characters).
//' @param AttributeName Attribute name (max. 127 characters).
//' @param Value Attribute string value (max. 255 characters).
//'
//' @family Data storage functions
//' @export
// [[Rcpp::export]]
SEXP AddAttributeString(SEXP Object, SEXP AttributeName, SEXP Value) {

  char *cObject = RtoCstring(Object);
  char *cAttributeName = RtoCstring(AttributeName);
  char *cValue = RtoCstring(Value);

  TwRetVal rv = TwAddAttributeString(cObject, cAttributeName, cValue);

  return TwRetValString(rv);
}

// AddUserData -----------------------------------------------------------------
//' Stores (asynchronous) user supplied data.
//'
//' \code{AddUserData} stores user supplied data asynchronously to the TOF data
//' acquistion into the current data file. Creates datasets "Data" and "Info" at
//' \code{Location}.
//'
//' @param Location Location of group in HDF5 file where the datasets are created.
//' @param NbrElements Number of elements to store (per call to this function),
//' maximum is 1048575.
//' @param Data Vector of length \code{NbrElements} containing the data to be
//' stored in dataset "Data".
//' @param ElementDescription Vector of length \code{NbrElements} containing the
//' text description of elements. If \code{ElementDescription} is an empty string
//' the dataset "Info" is not created.
//' @param CompressionLevel ZLIB compression level (0-9) for dataset creation.
//' If the dataset at Location already exists this parameter has no effect.
//'
//' @family Data storage functions
//' @export
// [[Rcpp::export]]
SEXP AddUserData(SEXP Location, int NbrElements, NumericVector Data,
                 Nullable<Rcpp::String> ElementDescription = R_NilValue,
                 int CompressionLevel = 0) {

  char *cLocation = RtoCstring(Location);
  if (ElementDescription.isNotNull()) {
    char * cElementDescription = RtoCstring(ElementDescription);
    TwRetVal rv = TwAddUserData(cLocation, NbrElements, cElementDescription,
                                &Data[0], CompressionLevel);
    return TwRetValString(rv);
  } else {
    TwRetVal rv = TwAddUserData(cLocation, NbrElements, NULL,
                                &Data[0], CompressionLevel);
    return TwRetValString(rv);
  }
}

// AddUserDataMultiRow ---------------------------------------------------------
//' Stores (asynchronous) user supplied data.
//'
//' \code{AddUserDataMultiRow} stores user supplied data asynchronously to the TOF data
//' acquistion into the current data file. Creates datasets "Data" and "Info" at
//' \code{Location}.
//'
//' Same as \code{AddUserData}, but adds argument \code{NbrRows} to add several
//' lines of user data at once.
//'
//' @param Location Location of group in HDF5 file where the datasets are created.
//' @param NbrElements Number of elements to store (per call to this function),
//' maximum is 1048575.
//' @param NbrRows Number of rows to store per call to this function (each row
//' contains \code{NbrElements} entries), maximum is 2047.
//' @param Data Vector of length \code{NbrElements} containing the data to be
//' stored in dataset "Data".
//' @param ElementDescription Vector of length \code{NbrElements} containing the
//' text description of elements. If \code{ElementDescription} is an empty string
//' the dataset "Info" is not created.
//' @param CompressionLevel ZLIB compression level (0-9) for dataset creation.
//' If the dataset at Location already exists this parameter has no effect.
//'
//' @family Data storage functions
//' @export
// [[Rcpp::export]]
SEXP AddUserDataMultiRow(SEXP Location, int NbrElements, int NbrRows, NumericVector Data,
                         Nullable<Rcpp::String> ElementDescription = R_NilValue,
                         int CompressionLevel = 0) {

  char *cLocation = RtoCstring(Location);
  char *cElementDescription;
  if (ElementDescription.isNotNull()) {
    cElementDescription = RtoCstring(ElementDescription);
  } else {
    cElementDescription = NULL;
  }

  TwRetVal rv = TwAddUserDataMultiRow(cLocation, NbrElements, NbrRows,
                                      cElementDescription, &Data[0],
                                      CompressionLevel);

 return TwRetValString(rv);
}

// RegisterUserDataBuf ---------------------------------------------------------
//' Registers a data source to store (synchronous) user supplied data.
//'
//' \code{RegisterUserDataBuf} registers a data source to store user supplied
//' data synchronously to the TOF data acquistion (every buf) into the data file
//' being currently recorded. Creates datasets "TwData" and "TwInfo" at \code{Location}.
//'
//' Needs to be executed before starting the acquisition.
//' Use \code{\link{UpdateUserData}} to actually store the data.
//'
//' @param Location Location of group in HDF5 file where the datasets are created.
//' @param NbrElements Number of elements to store per buf.
//' @param ElementDescription Vector of length \code{NbrElements} containing the
//' text description of elements. If \code{ElementDescription} is an empty string
//' the dataset "TwInfo" is not created.
//' @param CompressionLevel Compression level used for data storage (0: no
//' compression, 1-9: increasing levels of compression (and CPU load)).
//'
//' @family Data storage functions
//' @export
// [[Rcpp::export]]
SEXP RegisterUserDataBuf(SEXP Location, int NbrElements,
                         Nullable<Rcpp::String> ElementDescription = R_NilValue,
                         int CompressionLevel = 0) {

  char *cLocation = RtoCstring(Location);

  char *cElementDescription;
  if (ElementDescription.isNotNull()) {
    cElementDescription = RtoCstring(ElementDescription);
  } else {
    cElementDescription = NULL;
  }

  TwRetVal rv = TwRegisterUserDataBuf(cLocation, NbrElements,
                                      cElementDescription, CompressionLevel);

  return TwRetValString(rv);
}

// RegisterUserDataWrite -------------------------------------------------------
//' Registers a data source to store (synchronous) user supplied data.
//'
//' \code{RegisterUserDataWrite} registers a data source to store user supplied
//' data synchronously to the TOF data acquistion (every write) into the data
//' file being currently recorded. Creates datasets "TwData" and "TwInfo" at
//' \code{Location}.
//'
//' Needs to be executed before starting the acquisition.
//' Use \code{\link{UpdateUserData}} to actually store the data.
//'
//' @param Location Location of group in HDF5 file where the datasets are created.
//' @param NbrElements Number of elements to store per write.
//' @param ElementDescription Vector of length \code{NbrElements} containing the
//' text description of elements. If \code{ElementDescription} is an empty string
//' the dataset "TwInfo" is not created.
//' @param CompressionLevel Compression level used for data storage (0: no
//' compression, 1-9: increasing levels of compression (and CPU load)).
//'
//' @family Data storage functions
//' @export
// [[Rcpp::export]]
SEXP RegisterUserDataWrite(SEXP Location, int NbrElements,
                           Nullable<Rcpp::String> ElementDescription = R_NilValue,
                           int CompressionLevel = 0) {

  char *cLocation = RtoCstring(Location);

  char *cElementDescription;
  if (ElementDescription.isNotNull()) {
    cElementDescription = RtoCstring(ElementDescription);
  } else {
    cElementDescription = NULL;
  }

  TwRetVal rv = TwRegisterUserDataWrite(cLocation, NbrElements,
                                        cElementDescription, CompressionLevel);

  return TwRetValString(rv);
}

// RegisterUserDataNoStore -----------------------------------------------------
//' Registers a data source for (synchronous) user supplied data.
//'
//' \code{RegisterUserDataNoStore} registers a data for user supplied
//' data (synchronous to the TOF data acquistion) but the data is not stored
//' in the data file.
//'
//' Needs to be executed before starting the acquisition.
//' Use \code{\link{UpdateUserData}} to update the data.
//'
//' @param Location Location of group in HDF5 file where the datasets are created.
//' @param NbrElements Number of elements to store per write.
//' @param ElementDescription Vector of length \code{NbrElements} containing the
//' text description of elements. If \code{ElementDescription} is an empty string
//' the dataset "TwInfo" is not created.
//'
//' @family Data storage functions
//' @export
// [[Rcpp::export]]
SEXP RegisterUserDataNoStore(SEXP Location, int NbrElements,
                             Nullable<Rcpp::String> ElementDescription = R_NilValue) {

  char *cLocation = RtoCstring(Location);

  char *cElementDescription;
  if (ElementDescription.isNotNull()) {
    cElementDescription = RtoCstring(ElementDescription);
  } else {
    cElementDescription = NULL;
  }

  TwRetVal rv = TwRegisterUserDataNoStore(cLocation, NbrElements,
                                        cElementDescription);

  return TwRetValString(rv);
}

// UnregisterUserData ----------------------------------------------------------
//' Unregisters a data source.
//'
//' \code{UnregisterUserData} unregisters a data source previously registered
//' with \code{\link{RegisterUserDataBuf}} or \code{\link{RegisterUserDataWrite}}.
//'
//' @param Location Location of group in HDF5 file identifying the user data.
//'
//' @family Data storage functions
//' @export
// [[Rcpp::export]]
SEXP UnregisterUserData(SEXP Location) {

  char *cLocation = RtoCstring(Location);

  TwRetVal rv = TwUnregisterUserData(cLocation);

  return TwRetValString(rv);
}

// UpdateUserData --------------------------------------------------------------
//' Updates the values for a registered data source.
//'
//' \code{UpdateUserData} updates the values for a registered data source.
//'
//' @param Location Location of group in HDF5 file where the datasets are created.
//' @param NbrElements Number of elements to update.
//' @param Data Vector of length \code{NbrElements} containing the new data.
//'
//' @family Data storage functions
//' @export
// [[Rcpp::export]]
SEXP UpdateUserData(SEXP Location, int NbrElements, NumericVector Data) {

  char *cLocation = RtoCstring(Location);

  TwRetVal rv = TwUpdateUserData(cLocation, NbrElements, &Data[0]);

  return TwRetValString(rv);
}

// ReadRegUserData -------------------------------------------------------------
//' Reads the current values of a registered data source.
//'
//' \code{ReadRegUserData} reads the current values of a registered data source.
//'
//' @param Location Location of group in HDF5 file where the datasets are created.
//' @param NbrElements Number of elements to read.
//' @return Vector containing the registered user data.
//'
//' @family Data storage functions
//' @export
// [[Rcpp::export]]
SEXP ReadRegUserData(SEXP Location, int NbrElements) {

  char *cLocation = RtoCstring(Location);

  NumericVector Data(NbrElements);

  TwRetVal rv = TwReadRegUserData(cLocation, NbrElements, &Data[0]);

  if (rv != TwSuccess) {
    return TwRetValString(rv);
  }

  return Data;

}

// QueryRegUserDataSize --------------------------------------------------------
//' Queries the size of a registered data source.
//'
//' \code{QueryRegUserDataSize} queries the size (nbrElements) of a registered data source.
//'
//' @param Location Location of group in HDF5 file identifying the registered data.
//'
//' @family Data storage functions
//' @export
// [[Rcpp::export]]
SEXP QueryRegUserDataSize(SEXP Location) {

  char *cLocation = RtoCstring(Location);

  IntegerVector NbrElements(1);

  TwRetVal rv = TwQueryRegUserDataSize(cLocation, &NbrElements[0]);

  if (rv != TwSuccess) {
    return TwRetValString(rv);
  }

  return NbrElements;
}

// GetRegUserDataSources -------------------------------------------------------
//' Queries names, dimensions and types of all registered data sources.
//'
//' \code{GetRegUserDataSources} queries names, dimensions and types of all
//' data sources currently registered in TofDaq recorder.
//'
//' @return List with the location, nbrElements and type of the data sources.
//' type 1: data source values are written to disk for every write,
//' type 2: data source values are written to disk for every buf.
//'
//' @family Data storage functions
//' @export
// [[Rcpp::export]]
SEXP GetRegUserDataSources() {

  int arrayLength = 0;

  TwRetVal rv = TwGetRegUserDataSources(&arrayLength, NULL, NULL, NULL);
  if (rv != TwValueAdjusted) {
    return TwRetValString(rv);
  }

  char *location = new char[256 * arrayLength];
  memset(location, 0, 256 * arrayLength);

  IntegerVector nbrElements(arrayLength);
  IntegerVector type(arrayLength);

  rv = TwGetRegUserDataSources(&arrayLength, location, &nbrElements[0], &type[0]);

  if (rv != TwSuccess) {
    delete[] location;
    return TwRetValString(rv);
  }

  CharacterVector locationArray(arrayLength);
  std::string str(location);
  delete[] location;

  for (int i = 0; i < arrayLength; ++i) {
    locationArray[i] = str.substr(i*256, 256);
  }

  List result;
  result["location"] = locationArray;
  result["nbrElements"] = nbrElements;
  result["type"] = type;

  return result;
}

// GetRegUserDataDesc ----------------------------------------------------------
//' Reads the element descriptions of a registered data source.
//'
//' \code{GetRegUserDataDesc} reads the element descriptions of a registered data source.
//'
//' @param Location Location of group in HDF5 file identifying the registered data.
//'
//' @family Data storage functions
//' @export
// [[Rcpp::export]]
SEXP GetRegUserDataDesc(SEXP Location) {

  char *cLocation = RtoCstring(Location);

  int nbrElements = 0;

  TwRetVal rv = TwGetRegUserDataDesc(cLocation, &nbrElements, NULL);
  if (rv != TwValueAdjusted) {
    return TwRetValString(rv);
  }

  char *elementDescription = new char[256 * nbrElements];
  memset(elementDescription, 0, 256 * nbrElements);

  rv = TwGetRegUserDataDesc(cLocation, &nbrElements, elementDescription);
  if (rv != TwSuccess) {
    delete[] elementDescription;
    return TwRetValString(rv);
  }

  CharacterVector descriptionArray(nbrElements);
  std::string str(elementDescription);
  delete[] elementDescription;

  for (int i = 0; i < nbrElements; ++i) {
    descriptionArray[i] = str.substr(i*256, 256);
  }

  return descriptionArray;
}

//' Allows to keep the data file open at the end of an acquisition.
//'
//' \code{KeepFileOpen} allows to keep the data file open at the end of an acquisition.
//'
//' @param keepOpen Issue \code{TRUE} (after DAQ start) to signal to recorder to
//' keep the data file open. When done adding data, issue \code{FALSE} to allow
//' the recorder to close the file.
//'
//' @family Data storage functions
//' @export
// [[Rcpp::export]]
SEXP KeepFileOpen(bool keepOpen) {

  TwRetVal rv = TwKeepFileOpen(keepOpen);

  return TwRetValString(rv);
}

// TpsConnect ------------------------------------------------------------------
//' Connects to a remote control enabled TPSController software.
//'
//' \code{TpsConnect} connects to a remote control enabled TPSController
//' software running on the same PC.
//'
//' @family TPS functions
//' @export
// [[Rcpp::export]]
SEXP TpsConnect() {

  TwRetVal rv = TwTpsConnect();

  return TwRetValString(rv);
}

// TpsConnect2 -----------------------------------------------------------------
//' Connects to a local or remote TPS.
//'
//' \code{TpsConnect2} connects to a local or remote TPS (TPS1: type = 0,
//' TPS2: type = 1)
//'
//' @param ip TPS2 host name or IP.
//' @param type TPS type (0: 1st generation TPS, 1: 2nd generation TPS).
//'
//' @family TPS functions
//' @export
// [[Rcpp::export]]
SEXP TpsConnect2(SEXP ip, int type) {

  if (type == 1) {
    char *cFilename = RtoCstring(ip);

    TwRetVal rv = TwTpsConnect2(cFilename, type);

    return TwRetValString(rv);

  } else if (type == 0) {

    TwRetVal rv = TwTpsConnect();

    return TwRetValString(rv);
  } else {
    return R_NilValue;
  }
}

// TpsDisconnect ---------------------------------------------------------------
//' Disconnects from a remote control enabled TPSController software.
//'
//' \code{TpsDisconnect} disconnects from a remote control enabled TPSController
//' software.
//'
//' @family TPS functions
//' @export
// [[Rcpp::export]]
SEXP TpsDisconnect() {

  TwRetVal rv = TwTpsDisconnect();

  return TwRetValString(rv);
}

// TpsGetMonitorVal ------------------------------------------------------------
//' Gets the last reported monitor value for a given module.
//'
//' \code{TpsGetMonitorValue} gets the last reported monitor value for a given
//' module.
//'
//' @param moduleCode Module code.
//'
//' @family TPS functions
//' @export
// [[Rcpp::export]]
SEXP TpsGetMonitorValue(int moduleCode) {

  NumericVector value(1);

  TwRetVal rv = TwTpsGetMonitorValue(moduleCode, &value[0]);
  if (rv != TwSuccess) {
    return TwRetValString(rv);
  }

  return value;
}

// TpsGetTargetValue -----------------------------------------------------------
//' Gets the last reported target value for a given module.
//'
//' \code{TpsGetTargetValue} gets the last reported target value for a given
//' module.
//'
//' @param moduleCode Module code.
//'
//' @family TPS functions
//' @export
// [[Rcpp::export]]
SEXP TpsGetTargetValue(int moduleCode) {

  NumericVector value(1);

  TwRetVal rv = TwTpsGetTargetValue(moduleCode, &value[0]);

  if (rv != TwSuccess) {
    return TwRetValString(rv);
  }

  return value;
}

// TpsGetLastSetValue ----------------------------------------------------------
//' Gets the last reported "last set" value for a given module.
//'
//' \code{TpsGetLastSetValue} gets the last reported "last set" value for a given
//' module.
//'
//' @param moduleCode Module code.
//'
//' @family TPS functions
//' @export
// [[Rcpp::export]]
SEXP TpsGetLastSetValue(int moduleCode) {

  NumericVector value(1);

  TwRetVal rv = TwTpsGetLastSetValue(moduleCode, &value[0]);

  if (rv != TwSuccess) {
    return TwRetValString(rv);
  }

  return value;
}

// TpsSetTargetValue -----------------------------------------------------------
//' Sets the target value for a given module.
//'
//' \code{TpsSetTargetValue} sets the target value for a given module.
//'
//' @param moduleCode Module code.
//' @param value Value to set.
//'
//' @family TPS functions
//' @export
// [[Rcpp::export]]
SEXP TpsSetTargetValue(int moduleCode, double value) {

  TwRetVal rv = TwTpsSetTargetValue(moduleCode, value);

  return TwRetValString(rv);
}

// TpsGetNbrModules ------------------------------------------------------------
//' Gets the number of controllable modules.
//'
//' \code{TpsGetNbrModules} gets the number of controllable modules.
//'
//' @family TPS functions
//' @export
// [[Rcpp::export]]
SEXP TpsGetNbrModules() {

  IntegerVector value(1);

  TwRetVal rv = TwTpsGetNbrModules(&value[0]);

  if (rv != TwSuccess) {
    return TwRetValString(rv);
  }

  return value;
}

// TpsGetModuleCodes -----------------------------------------------------------
//' Gets the module codes of all controllable TPS modules.
//'
//' \code{TpsGetModuleCodes} gets the module codes of all controllable TPS modules.
//'
//' @family TPS functions
//' @export
// [[Rcpp::export]]
SEXP TpsGetModuleCodes() {

  int nbrModules;

  TwRetVal rv = TwTpsGetNbrModules(&nbrModules);

  if (rv != TwSuccess) {
    return TwRetValString(rv);
  }

  IntegerVector moduleCodeBuffer(nbrModules);

  rv = TwTpsGetModuleCodes(&moduleCodeBuffer[0], nbrModules);

  if (rv != TwSuccess) {
    return TwRetValString(rv);
  }

  return moduleCodeBuffer;
}

// TpsInitialize ---------------------------------------------------------------
//' Initializes TPS.
//'
//' \code{TpsInitialize} initializes the TPS.
//'
//' @family TPS functions
//' @export
// [[Rcpp::export]]
SEXP TpsInitialize() {

  TwRetVal rv = TwTpsInitialize();

  return TwRetValString(rv);
}

// TpsSetAllVoltages -----------------------------------------------------------
//' Sets all voltages.
//'
//' \code{TpsSetAllVoltages} sets all voltages.
//'
//' @family TPS functions
//' @export
// [[Rcpp::export]]
SEXP TpsSetAllVoltages() {

  TwRetVal rv = TwTpsSetAllVoltages();

  return TwRetValString(rv);
}

// TpsShutdown -----------------------------------------------------------------
//' Shuts down TPS.
//'
//' \code{TpsShutdown} shuts down the TPS.
//'
//' @family TPS functions
//' @export
// [[Rcpp::export]]
SEXP TpsShutdown() {

  TwRetVal rv = TwTpsShutdown();

  return TwRetValString(rv);
}

// TpsGetStatus ----------------------------------------------------------------
//' Gets the status of the TPS.
//'
//' \code{TpsGetStatus} gets the status of the TPS.
//'
//' @return List with the status information: connected?, initialized?, shutdown?,
//' ion mode changable?, ion mode supported?, current ion mode?.
//'
//' @family TPS functions
//' @export
// [[Rcpp::export]]
SEXP TpsGetStatus() {

  int status;

  TwRetVal rv = TwTpsGetStatus(&status);
  if (rv != TwSuccess) {
    return TwRetValString(rv);
  }

  bool connected = (status & 0x01); // hex for 0000 0001
  bool initialized = (status & 0x02); // hex for 0000 0010
  bool shutdown = (status & 0x04); // hex for 0000 0100
  bool ionmode_changable = (status & 0x08); // hex for 0000 1000
  bool ionmode_supported = (status & 0x10); // hex for 0001 0000
  int current_ionmode = (status & 0x60) >> 5; // hex for 0110 0000
  StringVector string = StringVector::create("Inconsistent", "Undetermined",
                                             "Positive", "Negative");

  List result;
  result["connected"] = connected;
  result["initialized"] = initialized;
  result["shutdown"] = shutdown;
  result["ionmode_changable"] = ionmode_changable;
  result["ionmode_supported"] = ionmode_supported;
  result["current_ionmode"] = as<std::string>(string[current_ionmode-1]);

  return result;
}

// TpsLoadSetFile --------------------------------------------------------------
//' Loads a TPS set file from disk and sets all values.
//'
//' \code{TpsLoadSetFile} loads a TPS set file from disk and sets all values.
//'
//' @param setfile Path/filename of the set file to load.
//' @section Warning:
//' This does not just load the file (as the function name might suggest), but
//' also immediately sets all values.
//'
//' @family TPS functions
//' @export
// [[Rcpp::export]]
SEXP TpsLoadSetFile(SEXP setFile) {

  char *cFilename = RtoCstring(setFile);

  TwRetVal rv = TwTpsLoadSetFile(cFilename);

  return TwRetValString(rv);
}

// TpsSaveSetFile --------------------------------------------------------------
//' Saves the current TPS settings to a file.
//'
//' \code{TpsSaveSetFile} saves the current TPS settings to a file.
//'
//' @param setfile Path/filename of the set file to save.
//'
//' @family TPS functions
//' @export
// [[Rcpp::export]]
SEXP TpsSaveSetFile(SEXP setFile) {

  char *cFilename = RtoCstring(setFile);

  TwRetVal rv = TwTpsSaveSetFile(cFilename);

  return TwRetValString(rv);
}

// TpsGetActiveFilament --------------------------------------------------------
//' Gets the currently active filament.
//'
//' \code{TpsGetActiveFilament} gets the currently active filament.
//'
//' Note that \code{TpsGetMonitorValue} does not work to query the filament
//' number. Use this function instead.
//'
//' @return Returns 0 for Filament 1, 1 for Filament 2.
//'
//' @family TPS functions
//' @export
// [[Rcpp::export]]
SEXP TpsGetActiveFilament() {

  IntegerVector filament(1);

  TwRetVal rv = TwTpsGetActiveFilament(&filament[0]);

  if (rv != TwSuccess) {
    return TwRetValString(rv);
  }

  return filament;
}

// TpsSetActiveFilament --------------------------------------------------------
//' Sets the active filament.
//'
//' \code{TpsSetActiveFilament} sets the active filament.
//'
//' Note that \code{TpsSetTargetValue} does not work to set the filament
//' number. Use this function instead.
//'
//' @param activeFilament 0 for Filament 1, 1 for Filament 2.
//'
//' @family TPS functions
//' @export
// [[Rcpp::export]]
SEXP TpsSetActiveFilament(int activeFilament) {

  TwRetVal rv = TwTpsSetActiveFilament(activeFilament);

  return TwRetValString(rv);
}

// TpsGetModuleLimits ----------------------------------------------------------
//' Gets the limits for a given TPS module.
//'
//' \code{TpsGetModuleLimits} gets the (ion mode dependent) limits for a given TPS module. Only
//' works for TPS2.
//'
//' @param moduleCode Module code.
//'
//' @family TPS functions
//' @export
// [[Rcpp::export]]
SEXP TpsGetModuleLimits(int moduleCode) {

  NumericVector limit(2);

  TwRetVal rv = TwTpsGetModuleLimits(moduleCode, &limit[0], &limit[1]);
  if (rv != TwSuccess) {
    return TwRetValString(rv);
  }

  return limit;
}

// TpsChangeIonMode ------------------------------------------------------------
//' Changes ion mode and sets target values to 0.
//'
//' \code{TpsChangeIonMode} changes ion mode (and sets target values to 0).
//'
//' Note: this is an undocumented function of TofDaqDll.dll.
//'
//' @param ionMode 0: positive ion mode, 1: negative ion mode
//'
//' @family TPS functions
//' @export
// [[Rcpp::export]]
SEXP TpsChangeIonMode(int ionMode) {

  TwRetVal rv = TwTpsChangeIonMode(ionMode);

  return TwRetValString(rv);
}

#endif

// Not implemented: TwContinueAcquisition --------------------------------------
// Not implemented: TwManualContinueNeeded -------------------------------------
// Not implemented: TwLockBuf --------------------------------------------------
// Not implemented: TwUnLockBuf ------------------------------------------------
// Not implemented: TwIssueDio4Pulse -------------------------------------------
// Not implemented: TwSetDio4State ---------------------------------------------
// Not implemented: TwConfigVarNbrMemories -------------------------------------
// Not implemented: TwSetMassCalib ---------------------------------------------
// Not implemented: TwSetMassCalibEx -------------------------------------------
// Not implemented: TwGetSharedMemory ------------------------------------------
// Not implemented: TwGetMassCalib ---------------------------------------------
// Not implemented: TwGetMassCalibEx -------------------------------------------
