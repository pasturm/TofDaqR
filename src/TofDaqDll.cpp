#include <Rcpp.h>
using namespace Rcpp;
#ifdef _WIN32
#include <TofDaqDll.h>
#include "TofDaqR.h"
#endif

// InitializeDll ---------------------------------------------------------------
//' Initializes the TofDaqDll.dll.
//'
//' \code{InitializeDll} initializes the TofDaqDll.dll. It is usually not necessary
//' to call \code{InitializeDll} explicitly, as it is called automatically by
//' functions that need the DLL to be in an initialized state.
//' @export
// [[Rcpp::export]]
void InitializeDll() {
#ifdef _WIN32
  TwRetVal rv = TwInitializeDll();

  if (rv != TwSuccess) {
    stop(TranslateReturnValue(rv));
  }
#else
  stop("This function is only implemented on Windows.");
#endif
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
#ifdef _WIN32
  TwCleanupDll();
#else
  stop("This function is only implemented on Windows.");
#endif
}

// GetDllVersion ---------------------------------------------------------------
//' Gets the version number of the TofDaq API.
//'
//' \code{GetDllVersion} gets the version number of the TofDaq API.
//' @export
// [[Rcpp::export]]
double GetDllVersion() {
#ifdef _WIN32
  double rv = TwGetDllVersion();
  return rv;
#else
  stop("This function is only implemented on Windows.");
#endif
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
#ifdef _WIN32
  bool rv = TwTofDaqRunning();

  return rv;
#else
  stop("This function is only implemented on Windows.");
#endif
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
#ifdef _WIN32
  bool rv = TwDaqActive();

  return rv;
#else
  stop("This function is only implemented on Windows.");
#endif
}

// StartAcquisition ------------------------------------------------------------
//' Starts an acquisition.
//'
//' \code{StartAcquisition} starts an acquisition.
//' @export
// [[Rcpp::export]]
void StartAcquisition() {
#ifdef _WIN32
  TwRetVal rv = TwStartAcquisition();

  if (rv != TwSuccess) {
    stop(TranslateReturnValue(rv));
  }
#else
  stop("This function is only implemented on Windows.");
#endif
}

// StopAcquisition -------------------------------------------------------------
//' Stops the current acquisition.
//'
//' \code{StopAcquisition} stops the current acquisition.
//' @export
// [[Rcpp::export]]
void StopAcquisition() {
#ifdef _WIN32
  TwRetVal rv = TwStopAcquisition();

  if (rv != TwSuccess) {
    stop(TranslateReturnValue(rv));
  }
#else
  stop("This function is only implemented on Windows.");
#endif
}

// ContinueAcquisition ---------------------------------------------------------
//' Signals to the TofDaq recorder to continue an acquisition.
//'
//' \code{ContinueAcquisition} signals to the TofDaq recorder to continue an
//' acquisition.
//'
//' This is a legacy function that was used with some Acqiris DAQ cards,
//' where every block was armed by software. The feature is enabled by setting
//' the parameter ManualContinueEveryNMemories to a value > 0. All latest DAQ
//' devices operate in a streaming mode in order to achieve 100 \% duty cycle and
//' TofDaq recorder no longer has per block control of the DAQ progress.
//'
//' @export
// [[Rcpp::export]]
void ContinueAcquisition() {
#ifdef _WIN32
  TwRetVal rv = TwContinueAcquisition();

  if (rv != TwSuccess) {
    stop(TranslateReturnValue(rv));
  }
#else
  stop("This function is only implemented on Windows.");
#endif
}

// ManualContinueNeeded --------------------------------------------------------
//' Indicates if the TofDaq recorder expects a continue event.
//'
//' \code{ManualContinueNeeded} indicates if the TofDaq recorder expects a
//' continue event (see \code{\link{ContinueAcquisition}}).
//'
//' This is a legacy function that was used with some Acqiris DAQ cards,
//' where every block was armed by software. The feature is enabled by setting
//' the parameter ManualContinueEveryNMemories to a value > 0. All latest DAQ
//' devices operate in a streaming mode in order to achieve 100 \% duty cycle and
//' TofDaq recorder no longer has per block control of the DAQ progress.
//'
//' @return \code{TRUE} or \code{FALSE}.
//' @export
// [[Rcpp::export]]
bool ManualContinueNeeded() {
#ifdef _WIN32
  bool rv = TwManualContinueNeeded();

  return rv;
#else
  stop("This function is only implemented on Windows.");
#endif
}

// CloseTofDaqRec --------------------------------------------------------------
//' Closes the TofDaq recorder application.
//'
//' \code{CloseTofDaqRec} closes the TofDaq recorder application.
//' @export
// [[Rcpp::export]]
void CloseTofDaqRec() {
#ifdef _WIN32
  TwRetVal rv = TwCloseTofDaqRec();

  if (rv != TwSuccess) {
    stop(TranslateReturnValue(rv));
  }
#else
  stop("This function is only implemented on Windows.");
#endif
}

// IssueDio4Pulse --------------------------------------------------------------
//' Issues a TTL pulse on the digital output line 4.
//'
//' \code{IssueDio4Pulse} issues a TTL pulse on the digital output line 4
//' specified by a delay and a pulse width.
//'
//' Note that in order for this command to work the Dio4Mode parameter must be
//' set to 2 (pulsed) or 3 (manual).
//'
//' @param delay Delay before issuing pulse in ms.
//' @param width Pulse width in ms.
//'
//' @export
// [[Rcpp::export]]
void IssueDio4Pulse(int delay, int width) {
#ifdef _WIN32
  TwRetVal rv = TwIssueDio4Pulse(delay, width);

  if (rv != TwSuccess) {
    stop(TranslateReturnValue(rv));
  }
#else
  stop("This function is only implemented on Windows.");
#endif
}

// SetDio4State ----------------------------------------------------------------
//' Switches the digital output line 4 between states.
//'
//' \code{SetDio4State} switches the digital output line 4 between states.
//'
//' Note that in order for this command to work the Dio4Mode parameter must be
//' set to 2 (pulsed) or 3 (manual).
//'
//' @param state 0: idle state, 1 (or any value other than 0) active state.
//'
//' @export
// [[Rcpp::export]]
void SetDio4State(int state) {
#ifdef _WIN32
  TwRetVal rv = TwSetDio4State(state);

  if (rv != TwSuccess) {
    stop(TranslateReturnValue(rv));
  }
#else
  stop("This function is only implemented on Windows.");
#endif
}

// InitializeDaqDevice ---------------------------------------------------------
//' Initializes the DAQ board.
//'
//' \code{InitializeDaqDevice} initializes the DAQ board (this is also done at
//' startup of TofDaqRec.exe). This can take up to 8 seconds depending on the
//' actual DAQ hardware.
//' @export
// [[Rcpp::export]]
void InitializeDaqDevice() {
#ifdef _WIN32
  TwRetVal rv = TwInitializeDaqDevice();

  if (rv != TwSuccess) {
    stop(TranslateReturnValue(rv));
  }
#else
  stop("This function is only implemented on Windows.");
#endif
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
#ifdef _WIN32
  TwSetTimeout(timeout);
#else
  stop("This function is only implemented on Windows.");
#endif
}

// GetTimeout ------------------------------------------------------------------
//' Gets the timeout.
//'
//' \code{GetTimeout} gets the current timeout value (in ms).
//' @export
// [[Rcpp::export]]
int GetTimeout() {
#ifdef _WIN32
  int rv = TwGetTimeout();

  return rv;
#else
  stop("This function is only implemented on Windows.");
#endif
}

// AutoSetupDaqDevice ----------------------------------------------------------
//' Auto setup routine for the DAQ device.
//'
//' \code{AutoSetupDaqDevice} sets up AP240 or Ndigo5G DAQ device.
//'
//' This function is only functional for AP240 and Ndigo5G hardware. For all
//' other setups it returns success immediately but does not do anything.
//'
//' @export
// [[Rcpp::export]]
void AutoSetupDaqDevice() {
#ifdef _WIN32
  TwRetVal rv = TwAutoSetupDaqDevice();

  if (rv != TwSuccess) {
    stop(TranslateReturnValue(rv));
  }
#else
  stop("This function is only implemented on Windows.");
#endif
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
void OnDemandMassCalibration(int action) {
#ifdef _WIN32
  TwRetVal rv = TwOnDemandMassCalibration(action);

  if (rv != TwSuccess) {
    stop(TranslateReturnValue(rv));
  }
#else
  stop("This function is only implemented on Windows.");
#endif
}

// DioStartDelayActive ---------------------------------------------------------
//' Checks if TofDaq recorder has received a start signal.
//'
//' \code{DioStartDelayActive} checks if TofDaq recorder has received a start
//' signal and is waiting for DioStartDelay to pass.
//'
//' Signals \code{TRUE} when a start signal was received and "DioStartDelay" is
//' active. Goes \code{FALSE} when start delay has expired. In combination with
//' "TwWaitingForDioStartSignal" can be used to indicate to the user what the
//' current experiment state is when digital start signal (and delay) is used.
//'
//' @return \code{TRUE} or \code{FALSE}.
//' @export
// [[Rcpp::export]]
bool DioStartDelayActive() {
#ifdef _WIN32
  bool rv = TwDioStartDelayActive();

  return rv;
#else
  stop("This function is only implemented on Windows.");
#endif
}

// SendDioStartSignal ----------------------------------------------------------
//' Sends a digital start signal.
//'
//' \code{SendDioStartSignal} sends a digital start signal.
//'
//' Software override of digital start signal as configured with DioStart...
//' TofDaq recorder parameters. See also \code{\link{WaitingForDioStartSignal}}
//' for finding out when this function can be called successfully.
//'
//' @return \code{TRUE} or \code{FALSE}.
//' @export
// [[Rcpp::export]]
void SendDioStartSignal() {
#ifdef _WIN32
  TwRetVal rv = TwSendDioStartSignal();

  if (rv != TwSuccess) {
    stop(TranslateReturnValue(rv));
  }
#else
  stop("This function is only implemented on Windows.");
#endif
}

// WaitingForDioStartSignal ----------------------------------------------------
//' Checks if TofDaq recorder is waiting for a start signal.
//'
//' \code{WaitingForDioStartSignal} checks if TofDaq recorder is waiting for a
//' start signal.
//'
//' Allows to query whether TofDaq recorder is currently waiting for a digital
//' start signal (or the SW override signal, see \code{\link{SendDioStartSignal}}). Note
//' that a \code{\link{StartAcquisition}} needs to be issued before this function can return
//' \code{TRUE}. It can take several seconds from the moment \code{\link{StartAcquisition}} is
//' called until this function returns \code{TRUE}.
//'
//' @return \code{TRUE} or \code{FALSE}.
//' @export
// [[Rcpp::export]]
bool WaitingForDioStartSignal() {
#ifdef _WIN32
  bool rv = TwWaitingForDioStartSignal();

  return rv;
#else
  stop("This function is only implemented on Windows.");
#endif
}

// SaturationWarning -----------------------------------------------------------
//' Checks if the signal is saturating the DAQ system.
//'
//' \code{SaturationWarning} checks if the signal is saturating the analog input
//' of the DAQ system.
//'
//' Signals saturation of the input signal (saturation due to analog input
//' limitations, not MCP limit or non-linearity in signal which may happen at
//' significantly lower signal values).
//'
//' @return \code{TRUE} or \code{FALSE}.
//' @export
// [[Rcpp::export]]
bool SaturationWarning() {
#ifdef _WIN32
  bool rv = TwSaturationWarning();

  return rv;
#else
  stop("This function is only implemented on Windows.");
#endif
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
void ShowConfigWindow(int ConfigWindowIndex) {
#ifdef _WIN32
  TwRetVal rv = TwShowConfigWindow(ConfigWindowIndex);

  if (rv != TwSuccess) {
    stop(TranslateReturnValue(rv));
  }
#else
  stop("This function is only implemented on Windows.");
#endif
}

// LoadIniFile -----------------------------------------------------------------
//' Loads a configuration file.
//'
//' \code{LoadIniFile} loads a configuration file (*.ini) from disk.
//'
//' @param IniFile Path/filename of the configuration file. If no path is
//' specified, the TofDaq recorder directory will be used. If \code{IniFile} is
//' an empty string or \code{NULL}, "TwApiTmpIni.ini" will be loaded.
//' @export
// [[Rcpp::export]]
void LoadIniFile(Nullable<Rcpp::String> IniFile = R_NilValue) {
#ifdef _WIN32
  char *cFilename;
  if (IniFile.isNotNull()) {
    std::string str = as<std::string>(IniFile);
    if (!str.empty()) {
      cFilename = StringToChar(str);
    } else {
    cFilename = NULL;
    }
  } else {
    cFilename = NULL;
  }
  TwRetVal rv = TwLoadIniFile(cFilename);

  if (rv != TwSuccess) {
    stop(TranslateReturnValue(rv));
  }
#else
  stop("This function is only implemented on Windows.");
#endif
}

// SaveIniFile -----------------------------------------------------------------
//' Saves the current configuration (*.ini) to disk.
//'
//' \code{SaveIniFile} saves the current configuration (*.ini) to disk.
//'
//' @param IniFile Path/filename of the configuration file. If no path is
//' specified, the file will be saved in the TofDaq recorder directory.
//' If \code{IniFile} is an empty string or \code{NULL}, "TwApiTmpIni.ini"
//' will be used. If a path is specified, existing files cannot be overwritten.
//' @export
// [[Rcpp::export]]
void SaveIniFile(Nullable<Rcpp::String> IniFile = R_NilValue) {
#ifdef _WIN32
  char *cFilename;
  if (IniFile.isNotNull()) {
    std::string str = as<std::string>(IniFile);
    if (!str.empty()) {
      cFilename = StringToChar(str);
    } else {
      cFilename = NULL;
    }
  } else {
    cFilename = NULL;
  }
  TwRetVal rv = TwSaveIniFile(cFilename);

  if (rv != TwSuccess) {
    stop(TranslateReturnValue(rv));
  }
#else
  stop("This function is only implemented on Windows.");
#endif
}

// GetDaqParameter -------------------------------------------------------------
//' Gets a single parameter as a string.
//'
//' \code{GetDaqParameter} gets a single parameter as a string.
//'
//' @param Parameter Parameter name as a string. See
//' \href{https://htmlpreview.github.io/?https://github.com/pasturm/TofDaqR/blob/master/tools/doc/TofDaqDll.htm#parameter_list}{TofDaq API documentation}
//' for a list of all available parameters.
//' @export
// [[Rcpp::export]]
String GetDaqParameter(std::string Parameter) {
#ifdef _WIN32
  char *cParameter = StringToChar(Parameter);

  std::string str = TwGetDaqParameter(cParameter);

  return wrap(str);
#else
  stop("This function is only implemented on Windows.");
#endif
}

// GetDaqParameterInt ----------------------------------------------------------
//' Gets a single integer parameter.
//'
//' \code{GetDaqParameterInt} gets a single integer parameter.
//'
//' @param Parameter Parameter name as a string. See
//' \href{https://htmlpreview.github.io/?https://github.com/pasturm/TofDaqR/blob/master/tools/doc/TofDaqDll.htm#parameter_list}{TofDaq API documentation}
//' for a list of all available parameters.
//' @export
// [[Rcpp::export]]
int GetDaqParameterInt(std::string Parameter) {
#ifdef _WIN32
  char *cParameter = StringToChar(Parameter);

  int result = TwGetDaqParameterInt(cParameter);

  return result;
#else
  stop("This function is only implemented on Windows.");
#endif
}

// GetDaqParameterBool ---------------------------------------------------------
//' Gets a single boolean parameter.
//'
//' \code{GetDaqParameterBool} gets a single boolean parameter.
//'
//' @param Parameter Parameter name as a string. See
//' \href{https://htmlpreview.github.io/?https://github.com/pasturm/TofDaqR/blob/master/tools/doc/TofDaqDll.htm#parameter_list}{TofDaq API documentation}
//' for a list of all available parameters.
//' @export
// [[Rcpp::export]]
bool GetDaqParameterBool(std::string Parameter) {
#ifdef _WIN32
  char *cParameter = StringToChar(Parameter);

  bool result = TwGetDaqParameterBool(cParameter);

  return result;
#else
  stop("This function is only implemented on Windows.");
#endif
}

// GetDaqParameterFloat --------------------------------------------------------
//' Gets a single float parameter.
//'
//' \code{GetDaqParameterFloat} gets a single float parameter.
//'
//' @param Parameter Parameter name as a string. See
//' \href{https://htmlpreview.github.io/?https://github.com/pasturm/TofDaqR/blob/master/tools/doc/TofDaqDll.htm#parameter_list}{TofDaq API documentation}
//' for a list of all available parameters.
//' @export
// [[Rcpp::export]]
double GetDaqParameterFloat(std::string Parameter) {
#ifdef _WIN32
  char *cParameter = StringToChar(Parameter);

  float result = TwGetDaqParameterFloat(cParameter);

  return (double)result;
#else
  stop("This function is only implemented on Windows.");
#endif
}

// GetDaqParameterInt64 --------------------------------------------------------
//' Gets a single int64 parameter as a string.
//'
//' \code{GetDaqParameterInt64} gets a single int64 parameter as a string.
//'
//' The return string can be converted to integer64 using
//' \code{\link[bit64:as.integer64.character]{bit64::as.integer64()}}.
//'
//' @param Parameter Parameter name as a string. See
//' \href{https://htmlpreview.github.io/?https://github.com/pasturm/TofDaqR/blob/master/tools/doc/TofDaqDll.htm#parameter_list}{TofDaq API documentation}
//' for a list of all available parameters.
//' @export
// [[Rcpp::export]]
String GetDaqParameterInt64(std::string Parameter) {
#ifdef _WIN32
  char *cParameter = StringToChar(Parameter);

  __int64 result = TwGetDaqParameterInt64(cParameter);

  std::stringstream ss;
  ss << result;

  return ss.str();
#else
  stop("This function is only implemented on Windows.");
#endif
}

// GetDaqParameterDouble -------------------------------------------------------
//' Gets a single double parameter.
//'
//' \code{GetDaqParameterDouble} gets a single double parameter.
//'
//' @param Parameter Parameter name as a string. See
//' \href{https://htmlpreview.github.io/?https://github.com/pasturm/TofDaqR/blob/master/tools/doc/TofDaqDll.htm#parameter_list}{TofDaq API documentation}
//' for a list of all available parameters.
//' @export
// [[Rcpp::export]]
double GetDaqParameterDouble(std::string Parameter) {
#ifdef _WIN32
  char *cParameter = StringToChar(Parameter);

  double result = TwGetDaqParameterDouble(cParameter);

  return result;
#else
  stop("This function is only implemented on Windows.");
#endif
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
//' \href{https://htmlpreview.github.io/?https://github.com/pasturm/TofDaqR/blob/master/tools/doc/TofDaqDll.htm#parameter_list}{TofDaq API documentation}
//' for a list of all available parameters.
//' @export
// [[Rcpp::export]]
int GetDaqParameterIntRef(std::string Parameter) {
#ifdef _WIN32
  char *cParameter = StringToChar(Parameter);

  int Value;

  TwRetVal rv = TwGetDaqParameterIntRef(cParameter, &Value);
  if (rv != TwSuccess) {
    stop(TranslateReturnValue(rv));
  }

  return Value;
#else
  stop("This function is only implemented on Windows.");
#endif
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
//' \href{https://htmlpreview.github.io/?https://github.com/pasturm/TofDaqR/blob/master/tools/doc/TofDaqDll.htm#parameter_list}{TofDaq API documentation}
//' for a list of all available parameters.
//' @export
// [[Rcpp::export]]
bool GetDaqParameterBoolRef(std::string Parameter) {
#ifdef _WIN32
  char *cParameter = StringToChar(Parameter);

  bool Value;

  TwRetVal rv = TwGetDaqParameterBoolRef(cParameter, &Value);
  if (rv != TwSuccess) {
    stop(TranslateReturnValue(rv));
  }

  return Value;
#else
  stop("This function is only implemented on Windows.");
#endif
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
//' \href{https://htmlpreview.github.io/?https://github.com/pasturm/TofDaqR/blob/master/tools/doc/TofDaqDll.htm#parameter_list}{TofDaq API documentation}
//' for a list of all available parameters.
//' @export
// [[Rcpp::export]]
float GetDaqParameterFloatRef(std::string Parameter) {
#ifdef _WIN32
  char *cParameter = StringToChar(Parameter);

  float Value;

  TwRetVal rv = TwGetDaqParameterFloatRef(cParameter, &Value);
  if (rv != TwSuccess) {
    stop(TranslateReturnValue(rv));
  }

  return Value;
#else
  stop("This function is only implemented on Windows.");
#endif
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
//' \href{https://htmlpreview.github.io/?https://github.com/pasturm/TofDaqR/blob/master/tools/doc/TofDaqDll.htm#parameter_list}{TofDaq API documentation}
//' for a list of all available parameters.
//' @export
// [[Rcpp::export]]
String GetDaqParameterInt64Ref(std::string Parameter) {
#ifdef _WIN32
  char *cParameter = StringToChar(Parameter);

  __int64 Value;

  TwRetVal rv = TwGetDaqParameterInt64Ref(cParameter, &Value);
  if (rv != TwSuccess) {
    stop(TranslateReturnValue(rv));
  }

  std::stringstream ss;
  ss << Value;

  return wrap(ss.str());
#else
  stop("This function is only implemented on Windows.");
#endif
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
//' \href{https://htmlpreview.github.io/?https://github.com/pasturm/TofDaqR/blob/master/tools/doc/TofDaqDll.htm#parameter_list}{TofDaq API documentation}
//' for a list of all available parameters.
//' @export
// [[Rcpp::export]]
double GetDaqParameterDoubleRef(std::string Parameter) {
#ifdef _WIN32
  char *cParameter = StringToChar(Parameter);

  double Value;

  TwRetVal rv = TwGetDaqParameterDoubleRef(cParameter, &Value);
  if (rv != TwSuccess) {
    stop(TranslateReturnValue(rv));
  }

  return Value;
#else
  stop("This function is only implemented on Windows.");
#endif
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
//' \href{https://htmlpreview.github.io/?https://github.com/pasturm/TofDaqR/blob/master/tools/doc/TofDaqDll.htm#parameter_list}{TofDaq API documentation}
//' for a list of all available parameters.
//' @export
// [[Rcpp::export]]
String GetDaqParameterStringRef(std::string Parameter) {
#ifdef _WIN32
  char *cParameter = StringToChar(Parameter);
  char Value[256] = {};

  TwRetVal rv = TwGetDaqParameterStringRef(cParameter, Value);
  if (rv != TwSuccess) {
    stop(TranslateReturnValue(rv));
  }

  std::string str(Value);

  return wrap(str);
#else
  stop("This function is only implemented on Windows.");
#endif
}

// SetDaqParameter -------------------------------------------------------------
//' Sets a single parameter.
//'
//' \code{SetDaqParameter} sets a single parameter.
//'
//' @param Parameter Parameter name as a string. See
//' \href{https://htmlpreview.github.io/?https://github.com/pasturm/TofDaqR/blob/master/tools/doc/TofDaqDll.htm#parameter_list}{TofDaq API documentation}
//' for a list of all available parameters.
//' @param ValueString Value as a string.
//' @export
// [[Rcpp::export]]
void SetDaqParameter(std::string Parameter, std::string ValueString) {
#ifdef _WIN32
  char *cParameter = StringToChar(Parameter);
  char *cValueString = StringToChar(ValueString);

  TwRetVal rv = TwSetDaqParameter(cParameter, cValueString);

  if (rv != TwSuccess) {
    stop(TranslateReturnValue(rv));
  }
#else
  stop("This function is only implemented on Windows.");
#endif
}

// SetDaqParameterInt ----------------------------------------------------------
//' Sets a single parameter with an integer value.
//'
//' \code{SetDaqParameterInt} sets a single parameter with an integer value.
//'
//' @param Parameter Parameter name as a string. See
//' \href{https://htmlpreview.github.io/?https://github.com/pasturm/TofDaqR/blob/master/tools/doc/TofDaqDll.htm#parameter_list}{TofDaq API documentation}
//' for a list of all available parameters.
//' @param Value Integer value.
//' @export
// [[Rcpp::export]]
void SetDaqParameterInt(std::string Parameter, int Value) {
#ifdef _WIN32
  char *cParameter = StringToChar(Parameter);

  TwRetVal rv = TwSetDaqParameterInt(cParameter, Value);

  if (rv != TwSuccess) {
    stop(TranslateReturnValue(rv));
  }
#else
  stop("This function is only implemented on Windows.");
#endif
}

// SetDaqParameterBool ---------------------------------------------------------
//' Sets a single parameter with a boolean value.
//'
//' \code{SetDaqParameterBool} sets a single parameter with a boolean value.
//'
//' @param Parameter Parameter name as a string. See
//' \href{https://htmlpreview.github.io/?https://github.com/pasturm/TofDaqR/blob/master/tools/doc/TofDaqDll.htm#parameter_list}{TofDaq API documentation}
//' for a list of all available parameters.
//' @param Value \code{TRUE} or \code{FALSE}.
//' @export
// [[Rcpp::export]]
void SetDaqParameterBool(std::string Parameter, bool Value) {
#ifdef _WIN32
  char *cParameter = StringToChar(Parameter);

  TwRetVal rv = TwSetDaqParameterBool(cParameter, Value);

  if (rv != TwSuccess) {
    stop(TranslateReturnValue(rv));
  }
#else
  stop("This function is only implemented on Windows.");
#endif
}

// SetDaqParameterFloat --------------------------------------------------------
//' Sets a single parameter with a float value.
//'
//' \code{SetDaqParameterFloat} sets a single parameter with a float value.
//'
//' @param Parameter Parameter name as a string. See
//' \href{https://htmlpreview.github.io/?https://github.com/pasturm/TofDaqR/blob/master/tools/doc/TofDaqDll.htm#parameter_list}{TofDaq API documentation}
//' for a list of all available parameters.
//' @param Value Numeric value.
//' @export
// [[Rcpp::export]]
void SetDaqParameterFloat(std::string Parameter, double Value) {
#ifdef _WIN32
  char *cParameter = StringToChar(Parameter);

  TwRetVal rv = TwSetDaqParameterFloat(cParameter, (float)Value);

  if (rv != TwSuccess) {
    stop(TranslateReturnValue(rv));
  }
#else
  stop("This function is only implemented on Windows.");
#endif
}

// SetDaqParameterInt64 ----------------------------------------------------------
//' Sets a single parameter with an int64 value.
//'
//' \code{SetDaqParameterInt64} sets a single parameter with an int64 value
//' (passed as a string).
//'
//' @param Parameter Parameter name as a string. See
//' \href{https://htmlpreview.github.io/?https://github.com/pasturm/TofDaqR/blob/master/tools/doc/TofDaqDll.htm#parameter_list}{TofDaq API documentation}
//' for a list of all available parameters.
//' @param Value int64 value passed as a string.
//' @export
// [[Rcpp::export]]
void SetDaqParameterInt64(std::string Parameter, std::string Value) {
#ifdef _WIN32
  char *cParameter = StringToChar(Parameter);

  std::stringstream ss(Value);
  __int64 int64value;
  ss >> int64value;

  TwRetVal rv = TwSetDaqParameterInt64(cParameter, int64value);

  if (rv != TwSuccess) {
    stop(TranslateReturnValue(rv));
  }
#else
  stop("This function is only implemented on Windows.");
#endif
}

// SetDaqParameterDouble -------------------------------------------------------
//' Sets a single parameter with a double value.
//'
//' \code{SetDaqParameterDouble} sets a single parameter with a double value.
//'
//' @param Parameter Parameter name as a string. See
//' \href{https://htmlpreview.github.io/?https://github.com/pasturm/TofDaqR/blob/master/tools/doc/TofDaqDll.htm#parameter_list}{TofDaq API documentation}
//' for a list of all available parameters.
//' @param Value Numeric value.
//' @export
// [[Rcpp::export]]
void SetDaqParameterDouble(std::string Parameter, double Value) {
#ifdef _WIN32
  char *cParameter = StringToChar(Parameter);

  TwRetVal rv = TwSetDaqParameterDouble(cParameter, Value);

  if (rv != TwSuccess) {
    stop(TranslateReturnValue(rv));
  }
#else
  stop("This function is only implemented on Windows.");
#endif
}

// ConfigVarNbrMemories --------------------------------------------------------
//' Enables and configures the "variable NbrMemories" feature.
//'
//' \code{ConfigVarNbrMemories} enables and configures the "variable NbrMemories"
//' feature.
//'
//' @param Enable \code{TRUE} to enable or \code{FALSE} to disable "variable NbrMemories"
//' feature.
//' @param StepAtBuf Buf indices for each step.
//' @param NbrMemoriesForStep NbrMemories value for each step.
//' @export
// [[Rcpp::export]]
void ConfigVarNbrMemories(bool Enable, IntegerVector StepAtBuf,
                          IntegerVector NbrMemoriesForStep) {
#ifdef _WIN32
  int NbrSteps = StepAtBuf.size();
  int NbrSteps2 = NbrMemoriesForStep.size();
  if (NbrSteps != NbrSteps2) {
    stop("StepAtBuf and NbrMemoriesForStep must be the same length.");
  }

  TwRetVal rv = TwConfigVarNbrMemories(Enable, NbrSteps, &StepAtBuf[0],
                                       &NbrMemoriesForStep[0]);

  if (rv != TwSuccess) {
    stop(TranslateReturnValue(rv));
  }
#else
  stop("This function is only implemented on Windows.");
#endif
}

// SetMassCalib ----------------------------------------------------------------
//' Configures the mass calibration that will be used for the next acquisition.
//'
//' \code{SetMassCalib} configures the mass calibration that will be used for
//' the next acquisition(s). If \code{nbrParams} is 0, the calibration parameters are
//' determined by the TofDaq recorder based on the mass, tof and weight arrays.
//' If calibration parameters and calibration point information is supplied the
//' calibration parameters define the calibration (no "sanity" check is
//' performed whether the point information yields the same mass calibration
//' parameters).
//'
//' \tabular{cl}{
//' mode \tab Mass calibration function \cr
//' 0 \tab \eqn{i = p_1 \sqrt(m) + p_2} \cr
//' 1 \tab \eqn{i = p_1/\sqrt(m) + p_2} \cr
//' 2 \tab \eqn{i = p_1 m^{p_3} + p_2} \cr
//' 3 \tab \eqn{i = p_1 \sqrt(m) + p_2 + p_3 (m - p_4)^2} \cr
//' 4 \tab \eqn{i = p_1 \sqrt(m) + p_2 + p_3 m^2 + p_4 m + p_5} \cr
//' 5 \tab \eqn{m = p_1 i^2 + p_2 i + p_3}
//' }
//' Note: Modes 3 and 4 are flawed. Don't use them. In mode 3 the fit does not
//' converge well, because of a bug (parameters not correctly initialized).
//' Mode 4 is two sequential fits, first mode 0, then a quadratic fit to the
//' residuals, which is an inferior implementation of mode 3. Mode 1 is for FTMS
//' data.
//'
//' @param mode Mass calibration function to use.
//' @param nbrParams Number of mass calibration parameters.
//' @param p Vector with mass calibration parameters.
//' @param mass Vector with mass of the calibration points.
//' @param tof Vector with TOF sample index of the calibration points.
//' @param weight Vector with weight of the calibration points.
//'
//' @export
// [[Rcpp::export]]
void SetMassCalib(int mode, int nbrParams, NumericVector p, NumericVector mass,
                  NumericVector tof, NumericVector weight) {
#ifdef _WIN32
  int nbrPoints = mass.size();

  TwRetVal rv = TwSetMassCalib(mode, nbrParams, &p[0], nbrPoints, &mass[0],
                               &tof[0], &weight[0]);

  if (rv != TwSuccess) {
    stop(TranslateReturnValue(rv));
  }
#else
  stop("This function is only implemented on Windows.");
#endif
}

// SetMassCalib2 ---------------------------------------------------------------
//' Configures the mass calibration that will be used for the next acquisition.
//'
//' \code{SetMassCalib2} configures the mass calibration that will be used for
//' the next acquisition(s). If \code{nbrParams} is 0, the calibration parameters are
//' determined by the TofDaq recorder based on the mass, tof and weight arrays.
//' If calibration parameters and calibration point information is supplied the
//' calibration parameters define the calibration (no "sanity" check is
//' performed whether the point information yields the same mass calibration
//' parameters).
//'
//' \tabular{cl}{
//' mode \tab Mass calibration function \cr
//' 0 \tab \eqn{i = p_1 \sqrt(m) + p_2} \cr
//' 1 \tab \eqn{i = p_1/\sqrt(m) + p_2} \cr
//' 2 \tab \eqn{i = p_1 m^{p_3} + p_2} \cr
//' 3 \tab \eqn{i = p_1 \sqrt(m) + p_2 + p_3 (m - p_4)^2} \cr
//' 4 \tab \eqn{i = p_1 \sqrt(m) + p_2 + p_3 m^2 + p_4 m + p_5} \cr
//' 5 \tab \eqn{m = p_1 i^2 + p_2 i + p_3}
//' }
//' Note: Modes 3 and 4 are flawed. Don't use them. In mode 3 the fit does not
//' converge well, because of a bug (parameters not correctly initialized).
//' Mode 4 is two sequential fits, first mode 0, then a quadratic fit to the
//' residuals, which is an inferior implementation of mode 3. Mode 1 is for FTMS
//' data.
//'
//' @param mode Mass calibration function to use.
//' @param nbrParams Number of mass calibration parameters.
//' @param p Vector with mass calibration parameters.
//' @param mass Vector with mass of the calibration points.
//' @param tof Vector with TOF sample index of the calibration points.
//' @param weight Vector with weight of the calibration points.
//'
//' @export
// [[Rcpp::export]]
void SetMassCalib2(int mode, int nbrParams, NumericVector p,
                   NumericVector mass, NumericVector tof, NumericVector weight) {
#ifdef _WIN32
  int nbrPoints = mass.size();

  TwRetVal rv = TwSetMassCalib2(mode, nbrParams, &p[0], nbrPoints, &mass[0],
                      &tof[0], &weight[0]);

  if (rv != TwSuccess) {
    stop(TranslateReturnValue(rv));
  }
#else
  stop("This function is only implemented on Windows.");
#endif
}

// SetMassCalibEx --------------------------------------------------------------
//' Configures the mass calibration that will be used for the next acquisition.
//'
//' \code{SetMassCalibEx} configures the mass calibration that will be used for
//' the next acquisition(s). If \code{nbrParams} is 0, the calibration parameters are
//' determined by the TofDaq recorder based on the mass, tof and weight arrays.
//' If calibration parameters and calibration point information is supplied the
//' calibration parameters define the calibration (no "sanity" check is
//' performed whether the point information yields the same mass calibration
//' parameters). Labels to identify compound names/formulas used for
//' calibration have a maximum length of 255 characters.
//'
//' \tabular{cl}{
//' mode \tab Mass calibration function \cr
//' 0 \tab \eqn{i = p_1 \sqrt(m) + p_2} \cr
//' 1 \tab \eqn{i = p_1/\sqrt(m) + p_2} \cr
//' 2 \tab \eqn{i = p_1 m^{p_3} + p_2} \cr
//' 3 \tab \eqn{i = p_1 \sqrt(m) + p_2 + p_3 (m - p_4)^2} \cr
//' 4 \tab \eqn{i = p_1 \sqrt(m) + p_2 + p_3 m^2 + p_4 m + p_5} \cr
//' 5 \tab \eqn{m = p_1 i^2 + p_2 i + p_3}
//' }
//' Note: Modes 3 and 4 are flawed. Don't use them. In mode 3 the fit does not
//' converge well, because of a bug (parameters not correctly initialized).
//' Mode 4 is two sequential fits, first mode 0, then a quadratic fit to the
//' residuals, which is an inferior implementation of mode 3. Mode 1 is for FTMS
//' data.
//'
//' @param mode Mass calibration function to use.
//' @param nbrParams Number of mass calibration parameters.
//' @param p Vector with mass calibration parameters.
//' @param mass Vector with mass of the calibration points.
//' @param tof Vector with TOF sample index of the calibration points.
//' @param weight Vector with weight of the calibration points.
//' @param label Vector with labels of the calibration points.
//'
//' @export
// [[Rcpp::export]]
void SetMassCalibEx(int mode, int nbrParams, NumericVector p,
                    NumericVector mass, NumericVector tof, NumericVector weight,
                    StringVector label) {
#ifdef _WIN32
  int nbrPoints = mass.size();

  if (nbrPoints != label.size()) {
    stop("mass, tof, weight and label must be the same length.");
  }
  char *cLabel = new char[256 * nbrPoints];
  for( int i=0; i < nbrPoints; i++ ) {
    std::string str(label[i]);
    strncpy(&cLabel[i*256], str.c_str(), 256);
  }

  TwRetVal rv = TwSetMassCalibEx(mode, nbrParams, &p[0], nbrPoints, &mass[0],
                      &tof[0], &weight[0], cLabel);
  delete[] cLabel;

  if (rv != TwSuccess) {
    stop(TranslateReturnValue(rv));
  }
#else
  stop("This function is only implemented on Windows.");
#endif
}

// SetMassCalib2Ex -------------------------------------------------------------
//' Configures the mass calibration that will be used for the next acquisition.
//'
//' \code{SetMassCalib2Ex} configures the mass calibration that will be used for
//' the next acquisition(s). If \code{nbrParams} is 0, the calibration parameters are
//' determined by the TofDaq recorder based on the mass, tof and weight arrays.
//' If calibration parameters and calibration point information is supplied the
//' calibration parameters define the calibration (no "sanity" check is
//' performed whether the point information yields the same mass calibration
//' parameters). Labels to identify compound names/formulas used for
//' calibration have a maximum length of 255 characters.
//'
//' \tabular{cl}{
//' mode \tab Mass calibration function \cr
//' 0 \tab \eqn{i = p_1 \sqrt(m) + p_2} \cr
//' 1 \tab \eqn{i = p_1/\sqrt(m) + p_2} \cr
//' 2 \tab \eqn{i = p_1 m^{p_3} + p_2} \cr
//' 3 \tab \eqn{i = p_1 \sqrt(m) + p_2 + p_3 (m - p_4)^2} \cr
//' 4 \tab \eqn{i = p_1 \sqrt(m) + p_2 + p_3 m^2 + p_4 m + p_5} \cr
//' 5 \tab \eqn{m = p_1 i^2 + p_2 i + p_3}
//' }
//' Note: Modes 3 and 4 are flawed. Don't use them. In mode 3 the fit does not
//' converge well, because of a bug (parameters not correctly initialized).
//' Mode 4 is two sequential fits, first mode 0, then a quadratic fit to the
//' residuals, which is an inferior implementation of mode 3. Mode 1 is for FTMS
//' data.
//'
//' @param mode Mass calibration function to use.
//' @param nbrParams Number of mass calibration parameters.
//' @param p Vector with mass calibration parameters.
//' @param mass Vector with mass of the calibration points.
//' @param tof Vector with TOF sample index of the calibration points.
//' @param weight Vector with weight of the calibration points.
//' @param label Vector with labels of the calibration points.
//'
//' @export
// [[Rcpp::export]]
void SetMassCalib2Ex(int mode, int nbrParams, NumericVector p,
                     NumericVector mass, NumericVector tof, NumericVector weight,
                     StringVector label) {
#ifdef _WIN32
  int nbrPoints = mass.size();

  if (nbrPoints != label.size()) {
    stop("mass, tof, weight and label must be the same length.");
  }
  char *cLabel = new char[256 * nbrPoints];
  for( int i=0; i < nbrPoints; i++ ) {
    std::string str(label[i]);
    strncpy(&cLabel[i*256], str.c_str(), 256);
  }

  TwRetVal rv = TwSetMassCalib2Ex(mode, nbrParams, &p[0], nbrPoints, &mass[0],
                                 &tof[0], &weight[0], cLabel);
  delete[] cLabel;

  if (rv != TwSuccess) {
    stop(TranslateReturnValue(rv));
  }
#else
  stop("This function is only implemented on Windows.");
#endif
}

// ConfigureForSingleIonMeasurement --------------------------------------------
//' Configures TofDaq recorder for single ion measurements.
//'
//' \code{ConfigureForSingleIonMeasurement} configures TofDaq recorder for single
//' ion measurements.
//'
//' The function returns \code{nbrBits} and \code{negativeSignal} that need to be passed to
//' the single ion setup function in the tool library (all other parameters can
//' be read directly from TofDaq). Note that it is the user's responsibility to
//' backup current recorder settings before calling this function (and to revert
//' to this sertting after the SI run).
//'
//' @return List with \code{nbrBits} and \code{negativeSignal} to pass to SI config function.
//' @export
// [[Rcpp::export]]
List ConfigureForSingleIonMeasurement() {
#ifdef _WIN32
  int nbrBits;
  bool negativeSignal;

  TwRetVal rv = TwConfigureForSingleIonMeasurement(&nbrBits, &negativeSignal);
  if (rv != TwSuccess) {
    stop(TranslateReturnValue(rv));
  }

  List result;
  result["nbrBits"] = nbrBits;
  result["negativeSignal"] = negativeSignal;
  return result;
#else
  stop("This function is only implemented on Windows.");
#endif
}

// GetDescriptor ---------------------------------------------------------------
//' Gets various information about the active acquisition.
//'
//' \code{GetDescriptor} retrieves the current TSharedMemoryDesc structure.
//' TSharedMemoryDesc contains various static information about the active
//' acquisition that can be retrieved by \code{GetDaqParameter} functions but
//' also information of DAQ progress.
//' See
//' \href{https://htmlpreview.github.io/?https://github.com/pasturm/TofDaqR/blob/master/tools/doc/TofDaqDll.htm}{TofDaq API documentation}
//' for more details.
//'
//' int64 and unsigned int64 parameters are returned as string. They can be
//' converted to integer64 using \code{\link[bit64:as.integer64.character]{bit64::as.integer64()}}.
//'
//' @return A list containing the TSharedMemoryDesc structure
//' @export
// [[Rcpp::export]]
List GetDescriptor() {
#ifdef _WIN32
  //get the current TSharedMemoryDesc structure.
  TSharedMemoryDesc pBufDesc;
  TwRetVal rv = TwGetDescriptor(&pBufDesc);
  if (rv != TwSuccess) {
    stop(TranslateReturnValue(rv));
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
#else
  stop("This function is only implemented on Windows.");
#endif
}

// GetPeakParameters -----------------------------------------------------------
//' Gets parameters for a given peak.
//'
//' \code{GetPeakParameters} gets parameters for a given peak.
//'
//' @param PeakIndex Index of peak (zero-based numbering).
//' @return A list with the peak paramters \emph{label}, \emph{mass}, \emph{loMass} and \emph{hiMass}.
//' @export
// [[Rcpp::export]]
List GetPeakParameters(int PeakIndex) {
#ifdef _WIN32
  TPeakPar PeakPar;
  TwRetVal rv = TwGetPeakParameters(&PeakPar, PeakIndex);

  if (rv != TwSuccess) {
    stop(TranslateReturnValue(rv));
  }

  List result;
  result["label"] = PeakPar.label;
  result["mass"] = PeakPar.mass;
  result["loMass"] = PeakPar.loMass;
  result["hiMass"] = PeakPar.hiMass;

  return result;
#else
  stop("This function is only implemented on Windows.");
#endif
}

// ReleaseSharedMemory ---------------------------------------------------------
//' Manually releases the shared memory acquisition buffers.
//'
//' \code{ReleaseSharedMemory} manually releases the shared memory acquisition
//' buffers. This is needed if \code{\link{KeepSharedMemMapped}} has been set to
//' \code{TRUE}.
//' @export
// [[Rcpp::export]]
void ReleaseSharedMemory() {
#ifdef _WIN32
  TwRetVal rv = TwReleaseSharedMemory();

  if (rv != TwSuccess) {
    stop(TranslateReturnValue(rv));
  }
#else
  stop("This function is only implemented on Windows.");
#endif
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
void WaitForNewData(int timeout, bool WaitForEventReset = true) {
#ifdef _WIN32
  TSharedMemoryDesc pBufDesc;
  TSharedMemoryPointer pShMem;

  TwRetVal rv = TwWaitForNewData(timeout, &pBufDesc, &pShMem, WaitForEventReset);

  if (rv != TwSuccess) {
    stop(TranslateReturnValue(rv));
  }
#else
  stop("This function is only implemented on Windows.");
#endif
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
void WaitForEndOfAcquisition(int timeout) {
#ifdef _WIN32
  TwRetVal rv = TwWaitForEndOfAcquisition(timeout);

  if (rv != TwSuccess) {
    stop(TranslateReturnValue(rv));
  }
#else
  stop("This function is only implemented on Windows.");
#endif
}

// GetMassCalib ----------------------------------------------------------------
//' Returns information about the current mass calibration.
//'
//' \code{GetMassCalib} returns information about the mass calibration currently
//' used in TofDaq recorder.
//'
//' @return List with calibration parameters and calibration points.
//'
//' @export
// [[Rcpp::export]]
List GetMassCalib() {
#ifdef _WIN32
  int mode;
  int nbrParams = 0;
  int nbrPoints = 0;

  TwRetVal rv = TwGetMassCalib(&mode, &nbrParams, NULL, &nbrPoints, NULL, NULL, NULL);

  if (rv != TwValueAdjusted) {
    stop(TranslateReturnValue(rv));
  }

  NumericVector p(nbrParams);
  NumericVector mass(nbrPoints);
  NumericVector tof(nbrPoints);
  NumericVector weight(nbrPoints);

  rv = TwGetMassCalib(&mode, &nbrParams, &p[0], &nbrPoints, &mass[0], &tof[0],
                      &weight[0]);

  if (rv != TwSuccess) {
    stop(TranslateReturnValue(rv));
  }

  List result;
  result["mode"] = mode;
  result["p"] = p;
  result["mass"] = mass;
  result["tof"] = tof;
  result["weight"] = weight;
  return result;
#else
  stop("This function is only implemented on Windows.");
#endif
}

// GetMassCalib2 ---------------------------------------------------------------
//' Returns information about the current mass calibration.
//'
//' \code{GetMassCalib2} returns information about the mass calibration currently
//' used in TofDaq recorder.
//'
//' @return List with calibration parameters and calibration points.
//'
//' @export
// [[Rcpp::export]]
List GetMassCalib2() {
#ifdef _WIN32
  int mode;
  int nbrParams = 0;
  int nbrPoints = 0;

  TwRetVal rv = TwGetMassCalib2(&mode, &nbrParams, NULL, &nbrPoints, NULL, NULL, NULL);

  if (rv != TwValueAdjusted) {
    stop(TranslateReturnValue(rv));
  }

  NumericVector p(nbrParams);
  NumericVector mass(nbrPoints);
  NumericVector tof(nbrPoints);
  NumericVector weight(nbrPoints);

  rv = TwGetMassCalib2(&mode, &nbrParams, &p[0], &nbrPoints, &mass[0], &tof[0],
                      &weight[0]);

  if (rv != TwSuccess) {
    stop(TranslateReturnValue(rv));
  }

  List result;
  result["mode"] = mode;
  result["p"] = p;
  result["mass"] = mass;
  result["tof"] = tof;
  result["weight"] = weight;
  return result;
#else
  stop("This function is only implemented on Windows.");
#endif
}

// GetMassCalibEx --------------------------------------------------------------
//' Returns information about the current mass calibration.
//'
//' \code{GetMassCalibEx} returns information about the mass calibration currently
//' used in TofDaq recorder.
//'
//' This is the same as \code{\link{GetMassCalib}}, but additionally also returns
//' the labels of the calibration points.
//'
//' @return List with calibration parameters, calibration points and labels.
//'
//' @export
// [[Rcpp::export]]
List GetMassCalibEx() {
#ifdef _WIN32
  int mode;
  int nbrParams = 0;
  int nbrPoints = 0;

  TwRetVal rv = TwGetMassCalibEx(&mode, &nbrParams, NULL, &nbrPoints, NULL, NULL,
                                 NULL, NULL);

  if (rv != TwValueAdjusted) {
    stop(TranslateReturnValue(rv));
  }

  NumericVector p(nbrParams);
  NumericVector mass(nbrPoints);
  NumericVector tof(nbrPoints);
  NumericVector weight(nbrPoints);
  char *label = new char[256 * nbrPoints];

  rv = TwGetMassCalibEx(&mode, &nbrParams, &p[0], &nbrPoints, &mass[0], &tof[0],
                      &weight[0], label);

  if (rv != TwSuccess) {
    delete[] label;
    stop(TranslateReturnValue(rv));
  }

  CharacterVector labelArray(nbrPoints);
  std::string str(label, 256 * nbrPoints);
  delete[] label;

  for (int i = 0; i < nbrPoints; ++i) {
    labelArray[i] = str.substr(i*256, 256);
  }

  List result;
  result["mode"] = mode;
  result["p"] = p;
  result["mass"] = mass;
  result["tof"] = tof;
  result["weight"] = weight;
  result["label"] = labelArray;
  return result;
#else
  stop("This function is only implemented on Windows.");
#endif
}

// GetMassCalib2Ex -------------------------------------------------------------
//' Returns information about the current mass calibration.
//'
//' \code{GetMassCalib2Ex} returns information about the mass calibration currently
//' used in TofDaq recorder.
//'
//' This is the same as \code{\link{GetMassCalib2}}, but additionally also returns
//' the labels of the calibration points.
//'
//' @return List with calibration parameters, calibration points and labels.
//'
//' @export
// [[Rcpp::export]]
List GetMassCalib2Ex() {
#ifdef _WIN32
  int mode;
  int nbrParams = 0;
  int nbrPoints = 0;

  TwRetVal rv = TwGetMassCalib2Ex(&mode, &nbrParams, NULL, &nbrPoints, NULL, NULL,
                                 NULL, NULL);

  if (rv != TwValueAdjusted) {
    stop(TranslateReturnValue(rv));
  }

  NumericVector p(nbrParams);
  NumericVector mass(nbrPoints);
  NumericVector tof(nbrPoints);
  NumericVector weight(nbrPoints);
  char *label = new char[256 * nbrPoints];

  rv = TwGetMassCalib2Ex(&mode, &nbrParams, &p[0], &nbrPoints, &mass[0], &tof[0],
                        &weight[0], label);

  if (rv != TwSuccess) {
    delete[] label;
    stop(TranslateReturnValue(rv));
  }

  CharacterVector labelArray(nbrPoints);
  std::string str(label, 256 * nbrPoints);
  delete[] label;

  for (int i = 0; i < nbrPoints; ++i) {
    labelArray[i] = str.substr(i*256, 256);
  }

  List result;
  result["mode"] = mode;
  result["p"] = p;
  result["mass"] = mass;
  result["tof"] = tof;
  result["weight"] = weight;
  result["label"] = labelArray;
  return result;
#else
  stop("This function is only implemented on Windows.");
#endif
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
NumericVector GetSumSpectrumFromShMem(bool Normalize = true) {
#ifdef _WIN32
  TSharedMemoryDesc desc;
  TwRetVal rv = TwGetDescriptor(&desc);
  if (rv != TwSuccess) {
    stop(TranslateReturnValue(rv));
  }

  NumericVector Spectrum(desc.NbrSamples);

  rv = TwGetSumSpectrumFromShMem(&Spectrum[0], Normalize);
  if (rv != TwSuccess) {
    stop(TranslateReturnValue(rv));
  }

  return Spectrum;
#else
  stop("This function is only implemented on Windows.");
#endif
}

// GetSumSpectrumFromShMem2 ----------------------------------------------------
//' Sum spectrum from shared memory.
//'
//' \code{GetSumSpectrumFromShMem2} gets the sum spectrum from shared memory.
//'
//' @param Normalize If \code{FALSE} the spectrum is reported as sum,
//' if \code{TRUE} (default) the spectrum is normalized to counts per extraction.
//' @export
// [[Rcpp::export]]
NumericVector GetSumSpectrumFromShMem2(bool Normalize = true) {
#ifdef _WIN32
  TSharedMemoryDesc desc;
  TwRetVal rv = TwGetDescriptor(&desc);
  if (rv != TwSuccess) {
    stop(TranslateReturnValue(rv));
  }

  NumericVector Spectrum(desc.NbrSamples);

  rv = TwGetSumSpectrumFromShMem2(&Spectrum[0], Normalize);
  if (rv != TwSuccess) {
    stop(TranslateReturnValue(rv));
  }

  return Spectrum;
#else
  stop("This function is only implemented on Windows.");
#endif
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
//' if \code{TRUE} (default) the spectrum is normalized to counts per extraction
//' (ignored and assumed \code{FALSE} if used with \code{SegmentIndex = SegmentEndIndex = -1}).
//' @return A vector containing the mass spectrum or an array containing the
//' block of mass spectra if \code{SegmentIndex = SegmentEndIndex = -1}.
//' @export
// [[Rcpp::export]]
SEXP GetTofSpectrumFromShMem(int SegmentIndex, int SegmentEndIndex,
                             int BufIndex, bool Normalize = true) {
#ifdef _WIN32
  TSharedMemoryDesc desc;
  TwRetVal rv = TwGetDescriptor(&desc);
  if (rv != TwSuccess) {
    stop(TranslateReturnValue(rv));
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
    stop(TranslateReturnValue(rv));
  }

  return wrap(Spectrum);
#else
  stop("This function is only implemented on Windows.");
#endif
}

// GetTofSpectrumFromShMem2 ----------------------------------------------------
//' Single TOF spectrum from shared memory.
//'
//' \code{GetTofSpectrumFromShMem2} reads a single TOF spectrum (possibly
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
//' if \code{TRUE} (default) the spectrum is normalized to counts per extraction
//' (ignored and assumed \code{FALSE} if used with \code{SegmentIndex = SegmentEndIndex = -1}).
//' @return A vector containing the mass spectrum or an array containing the
//' block of mass spectra if \code{SegmentIndex = SegmentEndIndex = -1}.
//' @export
// [[Rcpp::export]]
SEXP GetTofSpectrumFromShMem2(int SegmentIndex, int SegmentEndIndex,
                             int BufIndex, bool Normalize = true) {
#ifdef _WIN32
  TSharedMemoryDesc desc;
  TwRetVal rv = TwGetDescriptor(&desc);
  if (rv != TwSuccess) {
    stop(TranslateReturnValue(rv));
  }

  int specLen;
  if (SegmentIndex == -1 && SegmentEndIndex == -1) {
    specLen = desc.NbrSamples*desc.NbrSegments;
  } else {
    specLen = desc.NbrSamples;
  }

  std::vector<float> Spectrum(specLen);
  rv = TwGetTofSpectrumFromShMem2(&Spectrum[0], SegmentIndex, SegmentEndIndex,
                                 BufIndex, Normalize);

  if (rv != TwSuccess) {
    stop(TranslateReturnValue(rv));
  }

  return wrap(Spectrum);
#else
  stop("This function is only implemented on Windows.");
#endif
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
NumericVector GetSpecXaxisFromShMem(int Type) {
#ifdef _WIN32
  //get descriptor of file
  TSharedMemoryDesc desc;
  TwRetVal rv = TwGetDescriptor(&desc);
  if (rv != TwSuccess) {
    stop(TranslateReturnValue(rv));
  }

  NumericVector SpecAxis(desc.NbrSamples);
  double maxMass = 0.0;
  rv = TwGetSpecXaxisFromShMem(&SpecAxis[0], Type, NULL, maxMass);
  if (rv != TwSuccess) {
    stop(TranslateReturnValue(rv));
  }

  return SpecAxis;
#else
  stop("This function is only implemented on Windows.");
#endif
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
List GetStickSpectrumFromShMem(int SegmentIndex, int SegmentEndIndex,
                               int BufIndex) {
#ifdef _WIN32
  TSharedMemoryDesc desc;
  TwRetVal rv = TwGetDescriptor(&desc);
  if (rv != TwSuccess) {
    stop(TranslateReturnValue(rv));
  }

  std::vector<float> Spectrum(desc.NbrPeaks);
  std::vector<float> Masses(desc.NbrPeaks);

  rv = TwGetStickSpectrumFromShMem(&Spectrum[0], &Masses[0], SegmentIndex,
                                   SegmentEndIndex, BufIndex);
  if (rv != TwSuccess) {
    stop(TranslateReturnValue(rv));
  }

  List result;
  result["Spectrum"] = Spectrum;
  result["Masses"] = Masses;

  return result;
#else
  stop("This function is only implemented on Windows.");
#endif
}

// GetStickSpectrumFromShMem2 --------------------------------------------------
//' Single stick spectrum from shared memory.
//'
//' \code{GetStickSpectrumFromShMem2} reads a single stick spectrum from shared
//' memory.
//'
//' @param SegmentIndex Segment start index of data to fetch.
//' @param SegmentEndIndex Segment end index of data to fetch.
//' @param BufIndex Buf index of data to fetch.
//' @return A list containing the stick spectrum and corresponding masses.
//' @export
// [[Rcpp::export]]
List GetStickSpectrumFromShMem2(int SegmentIndex, int SegmentEndIndex,
                               int BufIndex) {
#ifdef _WIN32
  TSharedMemoryDesc desc;
  TwRetVal rv = TwGetDescriptor(&desc);
  if (rv != TwSuccess) {
    stop(TranslateReturnValue(rv));
  }

  std::vector<float> Spectrum(desc.NbrPeaks);
  std::vector<float> Masses(desc.NbrPeaks);

  rv = TwGetStickSpectrumFromShMem2(&Spectrum[0], &Masses[0], SegmentIndex,
                                   SegmentEndIndex, BufIndex);
  if (rv != TwSuccess) {
    stop(TranslateReturnValue(rv));
  }

  List result;
  result["Spectrum"] = Spectrum;
  result["Masses"] = Masses;

  return result;
#else
  stop("This function is only implemented on Windows.");
#endif
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
#ifdef _WIN32
  TSharedMemoryDesc desc;
  TwRetVal rv = TwGetDescriptor(&desc);
  if (rv != TwSuccess) {
    stop(TranslateReturnValue(rv));
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
    stop(TranslateReturnValue(rv));
  }

  return wrap(SegmentProfile);
#else
  stop("This function is only implemented on Windows.");
#endif
}

// GetSegmentProfileFromShMem2 -------------------------------------------------
//' Segment profile for a given peak and buf index from shared memory.
//'
//' \code{GetSegmentProfileFromShMem2} reads the segment profile for a given
//' peak and buf index from shared memory. Use -1 for \code{PeakIndex} to get
//' segment profiles of all peaks.
//'
//' @param PeakIndex Index of peak to fetch segment profile from. All peaks are
//' read if \code{PeakIndex = -1}.
//' @param BufIndex Buf index of data to fetch.
//' @return A vector containing the segment profile(s).
//' @export
// [[Rcpp::export]]
SEXP GetSegmentProfileFromShMem2(int PeakIndex, int BufIndex) {
#ifdef _WIN32
  TSharedMemoryDesc desc;
  TwRetVal rv = TwGetDescriptor(&desc);
  if (rv != TwSuccess) {
    stop(TranslateReturnValue(rv));
  }

  int profileLen;
  if (PeakIndex == -1) {
    profileLen = desc.NbrSegments*desc.NbrPeaks;
  } else {
    profileLen = desc.NbrSegments;
  }

  std::vector<float> SegmentProfile(profileLen);
  rv = TwGetSegmentProfileFromShMem2(&SegmentProfile[0], PeakIndex, BufIndex);
  if (rv != TwSuccess) {
    stop(TranslateReturnValue(rv));
  }

  return wrap(SegmentProfile);
#else
  stop("This function is only implemented on Windows.");
#endif
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
double GetBufTimeFromShMem(int BufIndex, int WriteIndex) {
#ifdef _WIN32
  double BufTime;

  TwRetVal rv = TwGetBufTimeFromShMem(&BufTime, BufIndex, WriteIndex);
  if (rv != TwSuccess) {
    stop(TranslateReturnValue(rv));
  }

  return BufTime;
#else
  stop("This function is only implemented on Windows.");
#endif
}

// Data storage functions ------------------------------------------------------

// AddLogEntry -----------------------------------------------------------------
//' Adds an entry to the acquisition log.
//'
//' \code{AddLogEntry} adds an entry to the acquisition log.
//'
//' @param LogEntryText Log text (max. 255 characters).
//' @param LogEntryTime Log entry time (number of 100-nanosecond intervals since
//' January 1, 1601 UTC, Windows FILETIME) passed as a string. Set it to "0" for
//' "now".
//'
//' @family Data storage functions
//' @export
// [[Rcpp::export]]
void AddLogEntry(std::string LogEntryText, std::string LogEntryTime) {
#ifdef _WIN32
  char *cLogEntryText = StringToChar(LogEntryText);

  std::stringstream ss(LogEntryTime);
  unsigned __int64 cTime;
  ss >> cTime;

  TwRetVal rv = TwAddLogEntry(cLogEntryText, cTime);

  if (rv != TwSuccess) {
    stop(TranslateReturnValue(rv));
  }
#else
  stop("This function is only implemented on Windows.");
#endif
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
void AddAttributeInt(std::string Object, std::string AttributeName, int Value) {
#ifdef _WIN32
  char *cObject = StringToChar(Object);
  char *cAttributeName = StringToChar(AttributeName);

  TwRetVal rv = TwAddAttributeInt(cObject, cAttributeName, Value);

  if (rv != TwSuccess) {
    stop(TranslateReturnValue(rv));
  }
#else
  stop("This function is only implemented on Windows.");
#endif
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
void AddAttributeDouble(std::string Object, std::string AttributeName,
                        double Value) {
#ifdef _WIN32
  char *cObject = StringToChar(Object);
  char *cAttributeName = StringToChar(AttributeName);

  TwRetVal rv = TwAddAttributeDouble(cObject, cAttributeName, Value);

  if (rv != TwSuccess) {
    stop(TranslateReturnValue(rv));
  }
#else
  stop("This function is only implemented on Windows.");
#endif
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
void AddAttributeString(std::string Object, std::string AttributeName,
                        std::string Value) {
#ifdef _WIN32
  char *cObject = StringToChar(Object);
  char *cAttributeName = StringToChar(AttributeName);
  char *cValue = StringToChar(Value);

  TwRetVal rv = TwAddAttributeString(cObject, cAttributeName, cValue);

  if (rv != TwSuccess) {
    stop(TranslateReturnValue(rv));
  }
#else
  stop("This function is only implemented on Windows.");
#endif
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
//' text description of elements. If \code{ElementDescription} is \code{NULL}
//' the dataset "Info" is not created.
//' @param CompressionLevel ZLIB compression level (0-9) for dataset creation.
//' If the dataset at Location already exists this parameter has no effect.
//'
//' @family Data storage functions
//' @export
// [[Rcpp::export]]
void AddUserData(std::string Location, int NbrElements, NumericVector Data,
                 Nullable<Rcpp::StringVector> ElementDescription = R_NilValue,
                 int CompressionLevel = 0) {
#ifdef _WIN32
  char *cLocation = StringToChar(Location);
  char *cElementDescription = new char[256 * NbrElements];

  if (ElementDescription.isNotNull()) {
    StringVector strvec(ElementDescription); // https://stackoverflow.com/questions/43388698/rcpp-how-can-i-get-the-size-of-a-rcppnullable-numericvector
    for( int i=0; i < NbrElements; i++ ) {
      std::string str(strvec[i]);
      strncpy(&cElementDescription[i*256], str.c_str(), 256);
    }
  } else {
    cElementDescription = NULL;
  }
  TwRetVal rv = TwAddUserData(cLocation, NbrElements, cElementDescription,
                              &Data[0], CompressionLevel);
  delete[] cElementDescription;

  if (rv != TwSuccess) {
    stop(TranslateReturnValue(rv));
  }
#else
  stop("This function is only implemented on Windows.");
#endif
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
//' @param Data Vector of length \code{NbrElements*NbrRows} containing the data to be
//' stored in dataset "Data".
//' @param ElementDescription Vector of length \code{NbrElements} containing the
//' text description of elements. If \code{ElementDescription} is \code{NULL}
//' the dataset "Info" is not created.
//' @param CompressionLevel ZLIB compression level (0-9) for dataset creation.
//' If the dataset at Location already exists this parameter has no effect.
//'
//' @family Data storage functions
//' @export
// [[Rcpp::export]]
void AddUserDataMultiRow(std::string Location, int NbrElements, int NbrRows,
                         NumericVector Data,
                         Nullable<Rcpp::StringVector> ElementDescription = R_NilValue,
                         int CompressionLevel = 0) {
#ifdef _WIN32
  char *cLocation = StringToChar(Location);
  char *cElementDescription = new char[256 * NbrElements];

  if (ElementDescription.isNotNull()) {
    StringVector strvec(ElementDescription); // https://stackoverflow.com/questions/43388698/rcpp-how-can-i-get-the-size-of-a-rcppnullable-numericvector
    for( int i=0; i < NbrElements; i++ ) {
      std::string str(strvec[i]);
      strncpy(&cElementDescription[i*256], str.c_str(), 256);
    }
  } else {
    cElementDescription = NULL;
  }

  TwRetVal rv = TwAddUserDataMultiRow(cLocation, NbrElements, NbrRows,
                                      cElementDescription, &Data[0],
                                      CompressionLevel);
  delete[] cElementDescription;

  if (rv != TwSuccess) {
    stop(TranslateReturnValue(rv));
  }
#else
  stop("This function is only implemented on Windows.");
#endif
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
//' text description of elements. If \code{ElementDescription} is \code{NULL}
//' the dataset "TwInfo" is not created.
//' @param CompressionLevel Compression level used for data storage (0: no
//' compression, 1-9: increasing levels of compression (and CPU load)).
//'
//' @family Data storage functions
//' @export
// [[Rcpp::export]]
void RegisterUserDataBuf(std::string Location, int NbrElements,
                         Nullable<Rcpp::StringVector> ElementDescription = R_NilValue,
                         int CompressionLevel = 0) {
#ifdef _WIN32
  char *cLocation = StringToChar(Location);
  char *cElementDescription = new char[256 * NbrElements];

  if (ElementDescription.isNotNull()) {
    StringVector strvec(ElementDescription); // https://stackoverflow.com/questions/43388698/rcpp-how-can-i-get-the-size-of-a-rcppnullable-numericvector
    for( int i=0; i < NbrElements; i++ ) {
      std::string str(strvec[i]);
      strncpy(&cElementDescription[i*256], str.c_str(), 256);
    }
  } else {
    cElementDescription = NULL;
  }

  TwRetVal rv = TwRegisterUserDataBuf(cLocation, NbrElements,
                                      cElementDescription, CompressionLevel);
  delete[] cElementDescription;

  if (rv != TwSuccess) {
    stop(TranslateReturnValue(rv));
  }
#else
  stop("This function is only implemented on Windows.");
#endif
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
//' text description of elements. If \code{ElementDescription} is \code{NULL}
//' the dataset "TwInfo" is not created.
//' @param CompressionLevel Compression level used for data storage (0: no
//' compression, 1-9: increasing levels of compression (and CPU load)).
//'
//' @family Data storage functions
//' @export
// [[Rcpp::export]]
void RegisterUserDataWrite(std::string Location, int NbrElements,
                           Nullable<Rcpp::StringVector> ElementDescription = R_NilValue,
                           int CompressionLevel = 0) {
#ifdef _WIN32
  char *cLocation = StringToChar(Location);
  char *cElementDescription = new char[256 * NbrElements];

  if (ElementDescription.isNotNull()) {
    StringVector strvec(ElementDescription); // https://stackoverflow.com/questions/43388698/rcpp-how-can-i-get-the-size-of-a-rcppnullable-numericvector
    for( int i=0; i < NbrElements; i++ ) {
      std::string str(strvec[i]);
      strncpy(&cElementDescription[i*256], str.c_str(), 256);
    }
  } else {
    cElementDescription = NULL;
  }

  TwRetVal rv = TwRegisterUserDataWrite(cLocation, NbrElements,
                                        cElementDescription, CompressionLevel);
  delete[] cElementDescription;

  if (rv != TwSuccess) {
    stop(TranslateReturnValue(rv));
  }
#else
  stop("This function is only implemented on Windows.");
#endif
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
//' text description of elements. If \code{ElementDescription} is \code{NULL}
//' the dataset "TwInfo" is not created.
//'
//' @family Data storage functions
//' @export
// [[Rcpp::export]]
void RegisterUserDataNoStore(std::string Location, int NbrElements,
                             Nullable<Rcpp::StringVector> ElementDescription = R_NilValue) {
#ifdef _WIN32
  char *cLocation = StringToChar(Location);
  char *cElementDescription = new char[256 * NbrElements];

  if (ElementDescription.isNotNull()) {
    StringVector strvec(ElementDescription); // https://stackoverflow.com/questions/43388698/rcpp-how-can-i-get-the-size-of-a-rcppnullable-numericvector
    for( int i=0; i < NbrElements; i++ ) {
      std::string str(strvec[i]);
      strncpy(&cElementDescription[i*256], str.c_str(), 256);
    }
  } else {
    cElementDescription = NULL;
  }

  TwRetVal rv = TwRegisterUserDataNoStore(cLocation, NbrElements,
                                          cElementDescription);
  delete[] cElementDescription;

  if (rv != TwSuccess) {
    stop(TranslateReturnValue(rv));
  }
#else
  stop("This function is only implemented on Windows.");
#endif
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
void UnregisterUserData(std::string Location) {
#ifdef _WIN32
  char *cLocation = StringToChar(Location);

  TwRetVal rv = TwUnregisterUserData(cLocation);

  if (rv != TwSuccess) {
    stop(TranslateReturnValue(rv));
  }
#else
  stop("This function is only implemented on Windows.");
#endif
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
void UpdateUserData(std::string Location, int NbrElements, NumericVector Data) {
#ifdef _WIN32
  char *cLocation = StringToChar(Location);

  TwRetVal rv = TwUpdateUserData(cLocation, NbrElements, &Data[0]);

  if (rv != TwSuccess) {
    stop(TranslateReturnValue(rv));
  }
#else
  stop("This function is only implemented on Windows.");
#endif
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
NumericVector ReadRegUserData(std::string Location, int NbrElements) {
#ifdef _WIN32
  char *cLocation = StringToChar(Location);

  NumericVector Data(NbrElements);

  TwRetVal rv = TwReadRegUserData(cLocation, NbrElements, &Data[0]);

  if (rv != TwSuccess) {
    stop(TranslateReturnValue(rv));
  }

  return Data;
#else
  stop("This function is only implemented on Windows.");
#endif

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
int QueryRegUserDataSize(std::string Location) {
#ifdef _WIN32
  char *cLocation = StringToChar(Location);

  int NbrElements;

  TwRetVal rv = TwQueryRegUserDataSize(cLocation, &NbrElements);

  if (rv != TwSuccess) {
    stop(TranslateReturnValue(rv));
  }

  return NbrElements;
#else
  stop("This function is only implemented on Windows.");
#endif
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
List GetRegUserDataSources() {
#ifdef _WIN32
  int arrayLength = 0;

  TwRetVal rv = TwGetRegUserDataSources(&arrayLength, NULL, NULL, NULL);
  if (rv != TwValueAdjusted) {
    stop(TranslateReturnValue(rv));
  }

  char *location = new char[256 * arrayLength];

  IntegerVector nbrElements(arrayLength);
  IntegerVector type(arrayLength);

  rv = TwGetRegUserDataSources(&arrayLength, location, &nbrElements[0], &type[0]);

  if (rv != TwSuccess) {
    delete[] location;
    stop(TranslateReturnValue(rv));
  }

  CharacterVector locationArray(arrayLength);
  std::string str(location, 256 * arrayLength);
  delete[] location;

  for (int i = 0; i < arrayLength; ++i) {
    locationArray[i] = str.substr(i*256, 256);
  }

  List result;
  result["location"] = locationArray;
  result["nbrElements"] = nbrElements;
  result["type"] = type;

  return result;
#else
  stop("This function is only implemented on Windows.");
#endif
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
CharacterVector GetRegUserDataDesc(std::string Location) {
#ifdef _WIN32
  char *cLocation = StringToChar(Location);

  int nbrElements = 0;

  TwRetVal rv = TwGetRegUserDataDesc(cLocation, &nbrElements, NULL);
  if (rv != TwValueAdjusted) {
    stop(TranslateReturnValue(rv));
  }

  char *elementDescription = new char[256 * nbrElements];

  rv = TwGetRegUserDataDesc(cLocation, &nbrElements, elementDescription);
  if (rv != TwSuccess) {
    delete[] elementDescription;
    stop(TranslateReturnValue(rv));
  }

  CharacterVector descriptionArray(nbrElements);
  std::string str(elementDescription, 256 * nbrElements);
  delete[] elementDescription;

  for (int i = 0; i < nbrElements; ++i) {
    descriptionArray[i] = str.substr(i*256, 256);
  }

  return descriptionArray;
#else
  stop("This function is only implemented on Windows.");
#endif
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
void KeepFileOpen(bool keepOpen) {
#ifdef _WIN32
  TwRetVal rv = TwKeepFileOpen(keepOpen);

  if (rv != TwSuccess) {
    stop(TranslateReturnValue(rv));
  }
#else
  stop("This function is only implemented on Windows.");
#endif
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
void TpsConnect() {
#ifdef _WIN32
  TwRetVal rv = TwTpsConnect();

  if (rv != TwSuccess) {
    stop(TranslateReturnValue(rv));
  }
#else
  stop("This function is only implemented on Windows.");
#endif
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
void TpsConnect2(std::string ip, int type) {
#ifdef _WIN32
  char *cFilename = StringToChar(ip);

  TwRetVal rv = TwTpsConnect2(cFilename, type);

  if (rv != TwSuccess) {
    stop(TranslateReturnValue(rv));
  }
#else
  stop("This function is only implemented on Windows.");
#endif
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
void TpsDisconnect() {
#ifdef _WIN32
  TwRetVal rv = TwTpsDisconnect();

  if (rv != TwSuccess) {
    stop(TranslateReturnValue(rv));
  }
#else
  stop("This function is only implemented on Windows.");
#endif
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
double TpsGetMonitorValue(int moduleCode) {
#ifdef _WIN32
  double value;

  TwRetVal rv = TwTpsGetMonitorValue(moduleCode, &value);
  if (rv != TwSuccess) {
    stop(TranslateReturnValue(rv));
  }

  return value;
#else
  stop("This function is only implemented on Windows.");
#endif
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
double TpsGetTargetValue(int moduleCode) {
#ifdef _WIN32
  double value;

  TwRetVal rv = TwTpsGetTargetValue(moduleCode, &value);

  if (rv != TwSuccess) {
    stop(TranslateReturnValue(rv));
  }

  return value;
#else
  stop("This function is only implemented on Windows.");
#endif
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
double TpsGetLastSetValue(int moduleCode) {
#ifdef _WIN32
  double value;

  TwRetVal rv = TwTpsGetLastSetValue(moduleCode, &value);

  if (rv != TwSuccess) {
    stop(TranslateReturnValue(rv));
  }

  return value;
#else
  stop("This function is only implemented on Windows.");
#endif
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
void TpsSetTargetValue(int moduleCode, double value) {
#ifdef _WIN32
  TwRetVal rv = TwTpsSetTargetValue(moduleCode, value);

  if (rv != TwSuccess) {
    stop(TranslateReturnValue(rv));
  }
#else
  stop("This function is only implemented on Windows.");
#endif
}

// TpsGetNbrModules ------------------------------------------------------------
//' Gets the number of controllable modules.
//'
//' \code{TpsGetNbrModules} gets the number of controllable modules.
//'
//' @family TPS functions
//' @export
// [[Rcpp::export]]
int TpsGetNbrModules() {
#ifdef _WIN32
  int value;

  TwRetVal rv = TwTpsGetNbrModules(&value);

  if (rv != TwSuccess) {
    stop(TranslateReturnValue(rv));
  }

  return value;
#else
  stop("This function is only implemented on Windows.");
#endif
}

// TpsGetModuleCodes -----------------------------------------------------------
//' Gets the module codes of all controllable TPS modules.
//'
//' \code{TpsGetModuleCodes} gets the module codes of all controllable TPS modules.
//'
//' @family TPS functions
//' @export
// [[Rcpp::export]]
IntegerVector TpsGetModuleCodes() {
#ifdef _WIN32
  int nbrModules;

  TwRetVal rv = TwTpsGetNbrModules(&nbrModules);

  if (rv != TwSuccess) {
    stop(TranslateReturnValue(rv));
  }

  IntegerVector moduleCodeBuffer(nbrModules);

  rv = TwTpsGetModuleCodes(&moduleCodeBuffer[0], nbrModules);

  if (rv != TwSuccess) {
    stop(TranslateReturnValue(rv));
  }

  return moduleCodeBuffer;
#else
  stop("This function is only implemented on Windows.");
#endif
}

// TpsInitialize ---------------------------------------------------------------
//' Initializes TPS.
//'
//' \code{TpsInitialize} initializes the TPS.
//'
//' @family TPS functions
//' @export
// [[Rcpp::export]]
void TpsInitialize() {
#ifdef _WIN32
  TwRetVal rv = TwTpsInitialize();

  if (rv != TwSuccess) {
    stop(TranslateReturnValue(rv));
  }
#else
  stop("This function is only implemented on Windows.");
#endif
}

// TpsSetAllVoltages -----------------------------------------------------------
//' Sets all voltages.
//'
//' \code{TpsSetAllVoltages} sets all voltages.
//'
//' @family TPS functions
//' @export
// [[Rcpp::export]]
void TpsSetAllVoltages() {
#ifdef _WIN32
  TwRetVal rv = TwTpsSetAllVoltages();

  if (rv != TwSuccess) {
    stop(TranslateReturnValue(rv));
  }
#else
  stop("This function is only implemented on Windows.");
#endif
}

// TpsShutdown -----------------------------------------------------------------
//' Shuts down TPS.
//'
//' \code{TpsShutdown} shuts down the TPS.
//'
//' @family TPS functions
//' @export
// [[Rcpp::export]]
void TpsShutdown() {
#ifdef _WIN32
  TwRetVal rv = TwTpsShutdown();

  if (rv != TwSuccess) {
    stop(TranslateReturnValue(rv));
  }
#else
  stop("This function is only implemented on Windows.");
#endif
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
List TpsGetStatus() {
#ifdef _WIN32
  int status;

  TwRetVal rv = TwTpsGetStatus(&status);
  if (rv != TwSuccess) {
    stop(TranslateReturnValue(rv));
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
#else
  stop("This function is only implemented on Windows.");
#endif
}

// TpsLoadSetFile --------------------------------------------------------------
//' Loads a TPS set file and sets all values.
//'
//' \code{TpsLoadSetFile} loads a TPS set file and sets all values.
//'
//' @param setFile Path/filename of the set file to load.
//' @section Warning:
//' This does not just load the file (as the function name might suggest), but
//' also immediately sets all values.
//'
//' @family TPS functions
//' @export
// [[Rcpp::export]]
void TpsLoadSetFile(std::string setFile) {
#ifdef _WIN32
  char *cFilename = StringToChar(setFile);

  TwRetVal rv = TwTpsLoadSetFile(cFilename);

  if (rv != TwSuccess) {
    stop(TranslateReturnValue(rv));
  }
#else
  stop("This function is only implemented on Windows.");
#endif
}

// TpsLoadSetFile2 -------------------------------------------------------------
//' Loads a TPS set file and sets some values.
//'
//' \code{TpsLoadSetFile2} loads a TPS set file and only sets whitelisted
//' values or all values except blacklisted RC codes.
//'
//' The only 3 supported modes to call this function are:
//' \enumerate{
//'   \item \code{TpsLoadSetFile2(setFile, NULL, NULL)}, this is the same as \code{\link{TpsLoadSetFile}}
//'   \item \code{TpsLoadSetFile2(setFile, blackListArray, NULL)} sets the values from setFile except RC codes in blackListArray
//'   \item \code{TpsLoadSetFile2(setFile, NULL, whiteListArray)} sets only the values from setFile that are also in whiteListArray
//' }
//'
//' @param setFile Path/filename of the set file to load.
//' @param blackListArray RC code array for blacklist.
//' @param whiteListArray RC code array for whitelist .
//'
//' @family TPS functions
//' @export
// [[Rcpp::export]]
void TpsLoadSetFile2(std::string setFile,
                     Nullable<Rcpp::IntegerVector> blackListArray,
                     Nullable<Rcpp::IntegerVector> whiteListArray) {
#ifdef _WIN32
  char *cFilename = StringToChar(setFile);

  // TODO: copied from H5SetMassCalibDynamic -> is IntegerVector blackListArray_(blackListArray); really needed?
  int blackListLength;
  int whiteListLength;
  int *p_blackListArray;
  int *p_whiteListArray;
  if (blackListArray.isNotNull()) {
    IntegerVector blackListArray_(blackListArray);
    blackListLength = blackListArray_.size();
    p_blackListArray = &blackListArray_[0];
  } else {
    blackListLength = 0;
    p_blackListArray = nullptr;
  }
  if (whiteListArray.isNotNull()) {
    IntegerVector whiteListArray_(whiteListArray);
    whiteListLength = whiteListArray_.size();
    p_whiteListArray = &whiteListArray_[0];
  } else {
    whiteListLength = 0;
    p_whiteListArray = nullptr;
  }

  TwRetVal rv = TwTpsLoadSetFile2(cFilename, p_blackListArray, blackListLength,
                                 p_whiteListArray, whiteListLength);

  if (rv != TwSuccess) {
    stop(TranslateReturnValue(rv));
  }
#else
  stop("This function is only implemented on Windows.");
#endif
}

// TpsSaveSetFile --------------------------------------------------------------
//' Saves the current TPS settings to a file.
//'
//' \code{TpsSaveSetFile} saves the current TPS settings to a file.
//'
//' @param setFile Path/filename of the set file to save.
//'
//' @family TPS functions
//' @export
// [[Rcpp::export]]
void TpsSaveSetFile(std::string setFile) {
#ifdef _WIN32
  char *cFilename = StringToChar(setFile);

  TwRetVal rv = TwTpsSaveSetFile(cFilename);

  if (rv != TwSuccess) {
    stop(TranslateReturnValue(rv));
  }
#else
  stop("This function is only implemented on Windows.");
#endif
}

// TpsSaveSetFileRc ------------------------------------------------------------
//' Saves TPS set values with a RC code to a file.
//'
//' \code{TpsSaveSetFileRc} saves the current TPS settings to a file. Only
//' values with assigned RC codes will be saved.
//'
//' Note: set files saved with this function can not be loaded through the TPS
//' web GUI, only \code{\link{TpsLoadSetFile}} and \code{\link{TpsLoadSetFile2}}
//' understand this format.
//'
//' @param setFile Path/filename of the set file to save.
//'
//' @family TPS functions
//' @export
// [[Rcpp::export]]
void TpsSaveSetFileRc(std::string setFile) {
#ifdef _WIN32
  char *cFilename = StringToChar(setFile);

  TwRetVal rv = TwTpsSaveSetFileRc(cFilename);

  if (rv != TwSuccess) {
   stop(TranslateReturnValue(rv));
  }
#else
  stop("This function is only implemented on Windows.");
#endif
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
int TpsGetActiveFilament() {
#ifdef _WIN32
  int filament;

  TwRetVal rv = TwTpsGetActiveFilament(&filament);

  if (rv != TwSuccess) {
    stop(TranslateReturnValue(rv));
  }

  return filament;
#else
  stop("This function is only implemented on Windows.");
#endif
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
void TpsSetActiveFilament(int activeFilament) {
#ifdef _WIN32
  TwRetVal rv = TwTpsSetActiveFilament(activeFilament);

  if (rv != TwSuccess) {
    stop(TranslateReturnValue(rv));
  }
#else
  stop("This function is only implemented on Windows.");
#endif
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
NumericVector TpsGetModuleLimits(int moduleCode) {
#ifdef _WIN32
  NumericVector limit(2);

  TwRetVal rv = TwTpsGetModuleLimits(moduleCode, &limit[0], &limit[1]);
  if (rv != TwSuccess) {
    stop(TranslateReturnValue(rv));
  }

  return limit;
#else
  stop("This function is only implemented on Windows.");
#endif
}

// TpsChangeIonMode ------------------------------------------------------------
//' Changes ion mode and sets target values to 0.
//'
//' \code{TpsChangeIonMode} changes ion mode (and sets target values to 0).
//'
//' Note: This is an undocumented function of TofDaqDll.dll which must be used
//' with care. First shut the TPS down, then wait >30 s before changing the ion
//' mode. Not conforming to this might cause hardware damage.
//'
//' @param ionMode 0: positive ion mode, 1: negative ion mode
//'
//' @family TPS functions
//' @export
// [[Rcpp::export]]
void TpsChangeIonMode(int ionMode) {
#ifdef _WIN32
  TwRetVal rv = TwTpsChangeIonMode(ionMode);

  if (rv != TwSuccess) {
    stop(TranslateReturnValue(rv));
  }
#else
  stop("This function is only implemented on Windows.");
#endif
}

// TpsGetNmtState --------------------------------------------------------------
//' Queries the NMT state of a CANopen node.
//'
//' \code{TpsGetNmtState} queries the NMT (Network Management) state of the
//' CANopen node associated with RC code \code{moduleCode}. Possible returned
//' NMT states are: 0 (0x00, boot up), 4 (0x04, stopped), 5 (0x05, operational)
//' and 127 (0x7f, pre-operational).
//'
//' @param moduleCode Module code.
//'
//' @family TPS functions
//' @export
// [[Rcpp::export]]
int TpsGetNmtState(int moduleCode) {
#ifdef _WIN32
  int nmtState;

  TwRetVal rv = TwTpsGetNmtState(moduleCode, &nmtState);

  if (rv != TwSuccess) {
   stop(TranslateReturnValue(rv));
  }

 return nmtState;
#else
  stop("This function is only implemented on Windows.");
#endif
}

// TpsSetNmtCmd --------------------------------------------------------------
//' Sets the NMT state of a CANopen node.
//'
//' \code{TpsSetNmtCmd} sets the NMT (Network Management) state of the
//' CANopen node associated with RC code \code{moduleCode}. Valid values for
//' nmtState are 1 (0x01, operational), 2 (0x02, stop), 128 (0x80, pre-operational),
//' 129 (0x81, reset node) or 130 (0x82, reset communication).
//'
//' @param moduleCode Module code.
//' @param nmtState New NMT state for node .
//'
//' @family TPS functions
//' @export
// [[Rcpp::export]]
void TpsSetNmtCmd(int moduleCode, int nmtState) {
#ifdef _WIN32

  TwRetVal rv = TwTpsSetNmtCmd(moduleCode, nmtState);

  if (rv != TwSuccess) {
   stop(TranslateReturnValue(rv));
  }
#else
  stop("This function is only implemented on Windows.");
#endif
}

// TpsGetModuleProperties ----------------------------------------------------------------
//' Gets capabilities and label for a given RC code.
//'
//' \code{TpsGetModuleProperties} gets capabilities and label for a given RC code.
//'
//' @param moduleCode Module code.
//'
//' @return List with the properties (hasMonitor, isSettable, isTrigger) and
//' the label associated with \code{moduleCode} (can come from HW or cfg).
//'
//' @family TPS functions
//' @export
// [[Rcpp::export]]
List TpsGetModuleProperties(int moduleCode) {
#ifdef _WIN32
  int properties;
  char label[256] = {};

  TwRetVal rv = TwTpsGetModuleProperties(moduleCode, &properties, label);
  if (rv != TwSuccess) {
    stop(TranslateReturnValue(rv));
  }

  bool hasMonitor = (properties & 0x01); // hex for 0001
  bool isSettable = (properties & 0x02); // hex for 0010
  bool isTrigger = (properties & 0x04); // hex for 0100

  std::string str(label);

  List result;
  result["hasMonitor"] = hasMonitor;
  result["isSettable"] = isSettable;
  result["isTrigger"] = isTrigger;
  result["label"] = wrap(str);

return result;
#else
  stop("This function is only implemented on Windows.");
#endif
}

// Not implemented -------------------------------------------------------------
// TwGetSharedMemory
// TwSetRegUserDataTarget
// TwGetRegUserDataTargetRange
// TwGenerateSegmentProfilesFromEventData
// TwTpsSendPdo
// from API 1.99r1586:
// TwSetMassCalibInShMem
// TwGetMassCalibFromShMem
