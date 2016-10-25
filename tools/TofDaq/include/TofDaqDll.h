#ifndef TofDaqDllH
#define TofDaqDllH

#ifdef TOFDAQDLL_EXPORTS
#define TOFWERK_DAQ_API __declspec(dllexport)
#else
#define TOFWERK_DAQ_API __declspec(dllimport)
#endif

// #ifdef _MSC_VER
#if defined(_WIN32) && (defined(_MSC_VER) || defined(__GNUC__))
#define TwAddAttributeDouble          _TwAddAttributeDouble
#define TwAddAttributeInt             _TwAddAttributeInt
#define TwAddAttributeString          _TwAddAttributeString
#define TwAddLogEntry                 _TwAddLogEntry
#define TwAddUserData                 _TwAddUserData
#define TwAddUserDataMultiRow         _TwAddUserDataMultiRow
#define TwCleanupDll                  _TwCleanupDll
#define TwCloseTofDaqRec              _TwCloseTofDaqRec
#define TwConfigVarNbrMemories        _TwConfigVarNbrMemories
#define TwContinueAcquisition         _TwContinueAcquisition
#define TwDaqActive                   _TwDaqActive
#define TwDioStartDelayActive         _TwDioStartDelayActive
#define TwGetBufTimeFromShMem         _TwGetBufTimeFromShMem
#define TwGetDaqParameter             _TwGetDaqParameter
#define TwGetDaqParameterBool         _TwGetDaqParameterBool
#define TwGetDaqParameterBoolRef      _TwGetDaqParameterBoolRef
#define TwGetDaqParameterDouble       _TwGetDaqParameterDouble
#define TwGetDaqParameterDoubleRef    _TwGetDaqParameterDoubleRef
#define TwGetDaqParameterFloat        _TwGetDaqParameterFloat
#define TwGetDaqParameterFloatRef     _TwGetDaqParameterFloatRef
#define TwGetDaqParameterInt          _TwGetDaqParameterInt
#define TwGetDaqParameterInt64        _TwGetDaqParameterInt64
#define TwGetDaqParameterInt64Ref     _TwGetDaqParameterInt64Ref
#define TwGetDaqParameterIntRef       _TwGetDaqParameterIntRef
#define TwGetDaqParameterStringRef    _TwGetDaqParameterStringRef
#define TwGetDescriptor               _TwGetDescriptor
#define TwGetDllVersion               _TwGetDllVersion
#define TwGetMassCalib                _TwGetMassCalib
#define TwGetMassCalibEx              _TwGetMassCalibEx
#define TwGetMassCalib2               _TwGetMassCalib2
#define TwGetMassCalib2Ex             _TwGetMassCalib2Ex
#define TwGetPeakParameters           _TwGetPeakParameters
#define TwGetRegUserDataDesc          _TwGetRegUserDataDesc
#define TwGetRegUserDataSources       _TwGetRegUserDataSources
#define TwGetSegmentProfileFromShMem  _TwGetSegmentProfileFromShMem
#define TwGetSegmentProfileFromShMem2 _TwGetSegmentProfileFromShMem2
#define TwGetSharedMemory             _TwGetSharedMemory
#define TwGetSpecXaxisFromShMem       _TwGetSpecXaxisFromShMem
#define TwGetStickSpectrumFromShMem   _TwGetStickSpectrumFromShMem
#define TwGetStickSpectrumFromShMem2  _TwGetStickSpectrumFromShMem2
#define TwGetSumSpectrumFromShMem     _TwGetSumSpectrumFromShMem
#define TwGetSumSpectrumFromShMem2    _TwGetSumSpectrumFromShMem2
#define TwGetTimeout                  _TwGetTimeout
#define TwGetTofSpectrumFromShMem     _TwGetTofSpectrumFromShMem
#define TwGetTofSpectrumFromShMem2    _TwGetTofSpectrumFromShMem2
#define TwInitializeDaqDevice         _TwInitializeDaqDevice
#define TwInitializeDll               _TwInitializeDll
#define TwIssueDio4Pulse              _TwIssueDio4Pulse
#define TwLoadIniFile                 _TwLoadIniFile
#define TwLockBuf                     _TwLockBuf
#define TwManualContinueNeeded        _TwManualContinueNeeded
#define TwReadRegUserData             _TwReadRegUserData
#define TwQueryRegUserDataSize        _TwQueryRegUserDataSize
#define TwRegisterUserDataBuf         _TwRegisterUserDataBuf
#define TwRegisterUserDataWrite       _TwRegisterUserDataWrite
#define TwRegisterUserDataNoStore     _TwRegisterUserDataNoStore
#define TwReleaseSharedMemory         _TwReleaseSharedMemory
#define TwSaveIniFile                 _TwSaveIniFile
#define TwSetDaqParameter             _TwSetDaqParameter
#define TwSetDaqParameterBool         _TwSetDaqParameterBool
#define TwSetDaqParameterDouble       _TwSetDaqParameterDouble
#define TwSetDaqParameterFloat        _TwSetDaqParameterFloat
#define TwSetDaqParameterInt          _TwSetDaqParameterInt
#define TwSetDaqParameterInt64        _TwSetDaqParameterInt64
#define TwSetDio4State                _TwSetDio4State
#define TwSetMassCalib                _TwSetMassCalib
#define TwSetMassCalibEx              _TwSetMassCalibEx
#define TwSetMassCalib2               _TwSetMassCalib2
#define TwSetMassCalib2Ex             _TwSetMassCalib2Ex
#define TwSetTimeout                  _TwSetTimeout
#define TwShowAdvTiming               _TwShowAdvTiming
#define TwShowBasicTiming             _TwShowBasicTiming
#define TwShowConfigWindow            _TwShowConfigWindow
#define TwShowDaqWindow               _TwShowDaqWindow
#define TwShowFilesPaths              _TwShowFilesPaths
#define TwShowHdf5Config              _TwShowHdf5Config
#define TwShowMassCalib               _TwShowMassCalib
#define TwStartAcquisition            _TwStartAcquisition
#define TwStopAcquisition             _TwStopAcquisition
#define TwTofDaqRunning               _TwTofDaqRunning
#define TwTpsChangeIonMode			  _TwTpsChangeIonMode
#define TwTpsConnect                  _TwTpsConnect
#define TwTpsConnect2                 _TwTpsConnect2
#define TwTpsDisconnect               _TwTpsDisconnect
#define TwTpsGetActiveFilament        _TwTpsGetActiveFilament
#define TwTpsGetLastSetValue          _TwTpsGetLastSetValue
#define TwTpsGetModuleCodes           _TwTpsGetModuleCodes
#define TwTpsGetModuleLimits          _TwTpsGetModuleLimits
#define TwTpsGetMonitorValue          _TwTpsGetMonitorValue
#define TwTpsGetNbrModules            _TwTpsGetNbrModules
#define TwTpsGetStatus                _TwTpsGetStatus
#define TwTpsGetTargetValue           _TwTpsGetTargetValue
#define TwTpsInitialize               _TwTpsInitialize
#define TwTpsIonModeShutdown          _TwTpsIonModeShutdown
#define TwTpsLoadSetFile              _TwTpsLoadSetFile
#define TwTpsSaveSetFile              _TwTpsSaveSetFile
#define TwTpsSetActiveFilament        _TwTpsSetActiveFilament
#define TwTpsSetAllVoltages           _TwTpsSetAllVoltages
#define TwTpsSetTargetValue           _TwTpsSetTargetValue
#define TwTpsShutdown                 _TwTpsShutdown
#define TwUnLockBuf                   _TwUnLockBuf
#define TwUnregisterUserData          _TwUnregisterUserData
#define TwUpdateUserData              _TwUpdateUserData
#define TwWaitForEndOfAcquisition     _TwWaitForEndOfAcquisition
#define TwWaitForNewData              _TwWaitForNewData
#define TwAutoSetupDaqDevice          _TwAutoSetupDaqDevice
#define TwOnDemandMassCalibration     _TwOnDemandMassCalibration
#define TwKeepFileOpen				  _TwKeepFileOpen
#endif

#include "TofIpcStrucs.h"

//structure containing pointers to the TofDaqRec data in shared memory. Be aware that the
//data format of Raw* pointers is acquisition hardware / acquisition mode dependant.
typedef struct
	{
	double* SumSpectrum;
	double* SumSpectrum2;
	float** TofData;
	float** TofData2;
	float* PeakData;
	float* PeakData2;
	double* Timing;
	unsigned int* RawData32Ch1;
	unsigned short* RawData16Ch1;
	char* RawData8Ch1;
	unsigned int* RawData32Ch2;
	unsigned short* RawData16Ch2;
	char* RawData8Ch2;
	} TSharedMemoryPointer;
////////////////////////////////////////////////////////////////////////////////

#ifdef __cplusplus
extern "C" {
#endif
////////////////////////////////////////////////////////////////////////////////
TOFWERK_DAQ_API TwRetVal TwInitializeDll(void);
//initializes various stuff for dll (done previously in DllMain)
//arguments: none
//possible return values:	TwSuccess				if DLL was initialized
//							TwError					if DLL was already initialized
////////////////////////////////////////////////////////////////////////////////
TOFWERK_DAQ_API void TwCleanupDll(void);
//cleans up the various things produced by TwInitializeDll()
//arguments: none

////////////////////////////////////////////////////////////////////////////////
TOFWERK_DAQ_API bool TwTofDaqRunning(void);
//returns true if TofDaqRec.exe is running, false otherwise
//arguments: none

////////////////////////////////////////////////////////////////////////////////
TOFWERK_DAQ_API bool TwDaqActive(void);
//returns true if TofDaqRec.exe is working on an acquisition, false otherwise
//arguments: none

////////////////////////////////////////////////////////////////////////////////
//TOFWERK_DAQ_API char* TwTranslateReturnValue(TwRetVal ReturnValue);
//translates a TwRetVal into a character string
//arguments: TwRetVal from any TofDaq API function

////////////////////////////////////////////////////////////////////////////////
TOFWERK_DAQ_API TwRetVal TwStartAcquisition(void);
//starts an acquisition with the current settings
//arguments: none
//possible return values:	TwDaqRecNotRunning		if no recording application is found
//							TwAcquisitionActive		if there is already an active acquisition
//							TwTimeout   			if no successfull start was confirmed within TW_TIMEOUT ms
//							TwSuccess   			if Daq started normally
//							TwError					event abandoned or untreated error

////////////////////////////////////////////////////////////////////////////////
TOFWERK_DAQ_API TwRetVal TwStopAcquisition(void);
//stops the current acquisition
//arguments: none
//possible return values:	TwDaqRecNotRunning		if no recording application is found
//							TwNoActiveAcquisition   if there is no active acquisition
//							TwTimeout   			if no successfull stop was confirmed within TW_TIMEOUT ms
//							TwSuccess   			if acquisition was successfully stopped
//							TwError					untreated error

////////////////////////////////////////////////////////////////////////////////
TOFWERK_DAQ_API TwRetVal TwCloseTofDaqRec(void);
//closes the Recorder application
//arguments: none
//possible return values:	TwDaqRecNotRunning		if no recording application is found
//							TwTimeout   			no confirmation of closing within TW_TIMEOUT ms
//							TwSuccess   			if TofDaqRec.exe was successfully closed
//							TwError					untreated error

////////////////////////////////////////////////////////////////////////////////
TOFWERK_DAQ_API TwRetVal TwShowDaqWindow(void);
//shows the configuration dialog for Daq hardware configuration and autosave selection
//arguments: none
//possible return values:	TwSuccess   			Success
//							TwError					Error

////////////////////////////////////////////////////////////////////////////////
TOFWERK_DAQ_API TwRetVal TwShowBasicTiming(void);
//shows the configuration dialog for basic Tof timing
//arguments: none
//possible return values:	TwSuccess  				Success
//							TwError					Error

////////////////////////////////////////////////////////////////////////////////
TOFWERK_DAQ_API TwRetVal TwShowAdvTiming(void);
//shows the configuration dialog for timing mode selection and advanced timing parameters
//arguments: none
//possible return values:	TwSuccess   			Success
//							TwError					Error

////////////////////////////////////////////////////////////////////////////////
TOFWERK_DAQ_API TwRetVal TwShowFilesPaths(void);
//shows the configuration dialog for data path, .ini and mass table save path and viewer application
//arguments: none
//possible return values:	TwSuccess   			Success
//							TwError					Error

////////////////////////////////////////////////////////////////////////////////
TOFWERK_DAQ_API TwRetVal TwShowMassCalib(void);
//shows the configuration dialog for mass calibration/integration settings
//arguments: none
//possible return values:	TwSuccess   			Success
//							TwError					Error

////////////////////////////////////////////////////////////////////////////////
TOFWERK_DAQ_API TwRetVal TwShowHdf5Config(void);
//shows the configuration dialog for HDF5 related settings and HDF user block
//arguments: none
//possible return values:	TwSuccess   			Success
//							TwError					Error
////////////////////////////////////////////////////////////////////////////////
TOFWERK_DAQ_API TwRetVal TwLoadIniFile(char* IniFile);
//loads the configuration file (.ini). Not possible during acquisition.
//arguments: c string with an absolute path and file name
//possible return values:	TwDaqRecNotRunning		if no recording application is found
//							TwAcquisitionActive		if there is an active acquisition
//							TwFileNotFound   		if the specified ini file was not found
//							TwTimeout   			if ini file load was not confirmed within TW_TIMEOUT ms
//							TwSuccess   			if ini file was successfully loaded
////////////////////////////////////////////////////////////////////////////////
TOFWERK_DAQ_API TwRetVal TwSaveIniFile(char* IniFile);
//saves the the current configuration to IniFile. Not possible during acquisition.
//arguments: c string with an absolute path and file name
//possible return values:	TwDaqRecNotRunning		if no recording application is found
//							TwAcquisitionActive		if there is an active acquisition
//							TwFileNotFound   		if the specified directory does not exist
//							TwTimeout   			if ini file save was not confirmed within TW_TIMEOUT ms
//							TwSuccess   			if ini file was successfully saved
//							TwError					untreated error
////////////////////////////////////////////////////////////////////////////////
TOFWERK_DAQ_API TwRetVal TwGetDescriptor(TSharedMemoryDesc* pBufDesc);
//retrieves the current TSharedMemoryDesc structure
//arguments: pointer to user allocated TSharedMemoryDesc structure (see TofIpcStrucs.h)
//possible return values:	TwDaqRecNotRunning		if no recording application is found
//							TwSuccess   			if descriptor structure was successfully transferred
//							TwError					if pBufDesc is a NULL pointer
////////////////////////////////////////////////////////////////////////////////
TOFWERK_DAQ_API TwRetVal TwGetPeakParameters(TPeakPar* PeakPar, int PeakIndex);
//retrieves peak parameters
//arguments: pointer to user allocated TPeakPar structure (see TofIpcStrucs.h), index of peak of interest
//possible return values:	TwDaqRecNotRunning		if no recording application is found
//							TwSuccess   			if PeakPar structure was successfully transferred
//							TwOutOfBounds  			if PeakIndex is not in range 0 to NbrPeaks-1
//							TwError					if peak parameters could not be retrieved
////////////////////////////////////////////////////////////////////////////////
/*
TOFWERK_DAQ_API TwRetVal TwGetLastRecordedSpectrum(int SegmentIndex, int* NbrPoints, float** Spectrum);
//creates a copy of one mass spectrum in the last recorded buf (data stays valid until next call of this function)
//arguments: index of segment to copy, pointer to int (stores the number of data points in the spectrum),
//			 pointer to float pointer where space is allocated to hold the spectrum
//possible return values:	TwDaqRecNotRunning 		if no recording application is found
//							TwNoData   				if no data was found
//							TwOutOfBounds  			if SegmentIndex is not in range 0 to NbrSegments-1
//							TwSuccess   			if spectrum was successfully copied
//							TwError					if NbrPoints or Spectrum is a NULL pointer

////////////////////////////////////////////////////////////////////////////////
TOFWERK_DAQ_API TwRetVal TwGetSumSpectrum(int* NbrSpectra, int* NbrPoints, double** Spectrum);
//creates a copy of the sum spectrum  (data stays valid until next call of this function)
//arguments: pointer to int (number of spectra that were summed), pointer to int (stores the number of data points in the spectrum),
//			 pointer to double pointer where space is allocated to hold the spectrum
//possible return values:	TwDaqRecNotRunning		if no recording application is found
//							TwNoData   				if no data was found
//							TwSuccess   			if spectrum was successfully copied
//							TwError					if NbrSpectra, NbrPoints or Spectrum is a NULL pointer
*/
////////////////////////////////////////////////////////////////////////////////
TOFWERK_DAQ_API TwRetVal TwGetSharedMemory(TSharedMemoryPointer* pShMem, bool KeepSharedMemMapped);
//get pointers to shared memory for direct data access
//arguments: pointer to user allocated SharedMemoryPointer structure, KeepSharedMemMapped indicates that
//TofDaqDll should keep the shared memory mapped (even after the end of the acquisition). If KeepSharedMemMapped
//is true, it should be released manually with the TwReleaseSharedMemory() function. TwReleaseSharedMemory() is called
//automatically by the next call to TwStartAcquisition(). However there may be interferences if more than one instance of
//TofDaqDll communicates with the recording application.
//If KeepSharedMemMapped is false TofDaqDll will automatically unmap the shared memory buffers at the end of the acquisition.
//possible return values:	TwDaqRecNotRunning		if no recording application is found
//							TwNoData   				if no data was found
//							TwSuccess   			if shared memory was successfully mapped
//							TwError					if pShMem is a NULL pointer or memory mapping failed

////////////////////////////////////////////////////////////////////////////////
TOFWERK_DAQ_API TwRetVal TwWaitForNewData(int Timeout, TSharedMemoryDesc* pBufDesc, TSharedMemoryPointer* pShMem, bool WaitForEventReset);
//wait for new data and get updated memory descriptor and shared memory pointer
//arguments: timeout in millseconds, pointer to user allocated TSharedMemoryDesc structure (see TofIpcStrucs.h),
//			 pointer to user allocated SharedMemoryPointer structure (valid only until end of acquisition, see discussion)
// 			 in TwGetSharedMemory())
//			 WaitForEventReset waits until the dataavailable event is reset before returning
//possible return values:	TwDaqRecNotRunning		if no recording application is found
//							TwNoData   				if no data was found
//							TwSuccess   			if shared memory was successfully mapped
//							TwError					if BufferPointer or pBufDesc is a NULL pointer or untreated error
//							TwTimeout				if no new data became available during Timeout ms
////////////////////////////////////////////////////////////////////////////////
TOFWERK_DAQ_API TwRetVal TwSetMassCalib(int mode, int nbrParams, double* p, int nbrPoints, double* mass, double* tof, double* weight);
TOFWERK_DAQ_API TwRetVal TwSetMassCalibEx(int mode, int nbrParams, double* p, int nbrPoints, double* mass, double* tof, double* weight, char* label);
TOFWERK_DAQ_API TwRetVal TwSetMassCalib2(int mode, int nbrParams, double* p, int nbrPoints, double* mass, double* tof, double* weight);
TOFWERK_DAQ_API TwRetVal TwSetMassCalib2Ex(int mode, int nbrParams, double* p, int nbrPoints, double* mass, double* tof, double* weight, char* label);
//change mass calibration for next acquisition
//arguments: mass calibration parameters a and b, mass1, time1, mass2, time2 as in ini files.
//			 parameters a and b determine the new mass calibration, not m1,m2,t1,t2 !
//possible return values:	TwDaqRecNotRunning		if no recording application is found
//							TwAcquisitionActive   	if there is an active acquisition
//							TwSuccess   			if new mass calib parameters were transmitted to Recorder
//							TwTimeout				if no confirmation was received within TW_TIMEOUT ms
////////////////////////////////////////////////////////////////////////////////
TOFWERK_DAQ_API TwRetVal TwAddLogEntry(char* LogEntryText, unsigned __int64 LogEntryTime);
//adds a log entry to the current hdf file
//arguments: c string with entry text (max. 255 chars), UNIX timestamp (or 0 for now)
//possible return values:	TwDaqRecNotRunning		if no recording application is found
//							TwNoActiveAcquisition   if there is no active acquisition
//							TwSuccess   			if the log entry was successfully transmitted
//							TwTimeout				if no confirmation was received within TW_TIMEOUT ms
////////////////////////////////////////////////////////////////////////////////
TOFWERK_DAQ_API TwRetVal TwLockBuf(int TimeOut, int BufToLock);
//lock a given buf of data (TofData and PeakData). Make sure to unlock the buf once processing is
//complete, otherwise the acquisition will stall.
//arguments: timeout in ms, index of buf to lock
//possible return values:	TwDaqRecNotRunning		if no recording application is found
//							TwNoActiveAcquisition   if there is no active acquisition
//						    TwOutOfBounds   		if BufToLock is not in range 0 to NbrBufs-1
//							TwNoData  				if no data has been recorded into the buf to lock
//							TwSuccess   			lock successfully installed
//							TwTimeout				if the lock could not be installed during TimeOut
////////////////////////////////////////////////////////////////////////////////
TOFWERK_DAQ_API TwRetVal TwUnLockBuf(int BufToUnlock);
//unlock a given buf of data (TofData and PeakData).
//arguments: index of buf to unlock
//possible return values:	TwDaqRecNotRunning		if no recording application is found
//							TwNoActiveAcquisition   if there is no active acquisition
//						    TwOutOfBounds   		if BufToLock is not in range 0 to NbrBufs-1
//							TwSuccess   			lock successfully removed
//							TwError   				untreated error
////////////////////////////////////////////////////////////////////////////////
TOFWERK_DAQ_API TwRetVal TwAddAttributeInt(char* Object, char* AttributeName, int Value);
TOFWERK_DAQ_API	TwRetVal TwAddAttributeDouble(char* Object, char* AttributeName, double Value);
TOFWERK_DAQ_API TwRetVal TwAddAttributeString(char* Object, char* AttributeName, char* Value);
//attaches an attribute (int, double or string) into the current HDF5 file
//arguments: Object to attach attribute, key/value pair for attribute
//			 Object max. 255 chars, Name max. 127 chars, StringValue max. 255 chars
//possible return values:	TwDaqRecNotRunning	   	if no recording application is found
//							TwNoActiveAcquisition   if there is no active acquisition
//							TwSuccess 				attribute successfully added
//							TwError   				untreated error
////////////////////////////////////////////////////////////////////////////////
TOFWERK_DAQ_API TwRetVal TwSetDaqParameter(char* Parameter, char* ValueString);
//all parameters can be set by strings, helper functions to work directly with the corresponding datatypes:
TOFWERK_DAQ_API TwRetVal TwSetDaqParameterInt(char* Parameter, int Value);
TOFWERK_DAQ_API TwRetVal TwSetDaqParameterBool(char* Parameter, bool Value);
TOFWERK_DAQ_API TwRetVal TwSetDaqParameterFloat(char* Parameter, float Value);
TOFWERK_DAQ_API TwRetVal TwSetDaqParameterInt64(char* Parameter, __int64 Value);
TOFWERK_DAQ_API TwRetVal TwSetDaqParameterDouble(char* Parameter, double Value);
//sets a single DAQ parameter.
//arguments: parameter name as string (max. 127 chars)
//			 parameter value to set (max. 255 chars)
//possible return values:	TwDaqRecNotRunning	   	if no recording application is found
//							TwAcquisitionActive   	if there is an active acquisition
//							TwSuccess 				parameter set successfully
//							TwValueAdjusted			the supplied value was adjusted
//							TwInvalidParameter		unknown parameter
//							TwInvalidValue 			value not compatible with parameter
////////////////////////////////////////////////////////////////////////////////
TOFWERK_DAQ_API char* TwGetDaqParameter(char* Parameter);
//gets a single DAQ parameter.
//arguments: parameter name as string (max. 127 chars)
//
//possible return values:	"TwDaqRecNotRunning"	if no recording application is found
//							"TwInvalidParameter"	unknown parameter
//							"TwTimeout"				timeout
//							result (as string)		if OK
////////////////////////////////////////////////////////////////////////////////
TOFWERK_DAQ_API int TwGetDaqParameterInt(char* Parameter);
TOFWERK_DAQ_API bool TwGetDaqParameterBool(char* Parameter);
TOFWERK_DAQ_API float TwGetDaqParameterFloat(char* Parameter);
TOFWERK_DAQ_API __int64 TwGetDaqParameterInt64(char* Parameter);
TOFWERK_DAQ_API double TwGetDaqParameterDouble(char* Parameter);
//gets a single DAQ parameter as int, bool, float, int64, double.
//may be convenient, but almost no error checking available
//arguments: parameter name as string (max. 127 chars)
////////////////////////////////////////////////////////////////////////////////
TOFWERK_DAQ_API TwRetVal TwGetDaqParameterIntRef(char* Parameter, int* Value);
TOFWERK_DAQ_API TwRetVal TwGetDaqParameterBoolRef(char* Parameter, bool* Value);
TOFWERK_DAQ_API TwRetVal TwGetDaqParameterFloatRef(char* Parameter, float* Value);
TOFWERK_DAQ_API TwRetVal TwGetDaqParameterInt64Ref(char* Parameter, __int64* Value);
TOFWERK_DAQ_API TwRetVal TwGetDaqParameterDoubleRef(char* Parameter, double* Value);
TOFWERK_DAQ_API TwRetVal TwGetDaqParameterStringRef(char* Parameter, char* Value);
//gets a single DAQ parameter as int, bool, float, int64.
//the string version expects a char buffer of length 256
//standard interface with error checking
//arguments: parameter name as string (max. 127 chars), Value pointer of the appropriate type
//
//possible return values:	TwDaqRecNotRunning	   	if no recording application is found
//                          TwSuccess 				parameter got successfully
//							TwInvalidParameter		unknown parameter
//							TwInvalidValue 			value type not compatible with parameter
//							TwTimeOut				timeout
////////////////////////////////////////////////////////////////////////////////
TOFWERK_DAQ_API TwRetVal TwConfigVarNbrMemories(bool Enable, int NbrSteps, int* StepAtBuf, int* NbrMemoriesForStep);
//configures the variable NbrMemories feature. No validity check performed at setup (invalid StepAtBuf values
//are eliminated at DAQ start)
//arguments: enable feature (independant of config), number of changes in acquisition, array of length NbrSteps
//			 containing Buf indices of steps, array of length NbrSteps with the new NbrMemories value.
//
//possible return values:   TwDaqRecNotRunning	   	if no recording application is found
//                          TwAcquisitionActive   	if there is an active acquisition
//                          TwSuccess 				setup successful
//							TwTimeOut				timeout
//							TwError					error, NULL pointers with NbrSteps > 0, exception when reading values....
////////////////////////////////////////////////////////////////////////////////
TOFWERK_DAQ_API TwRetVal TwShowConfigWindow(int ConfigWindowIndex);
//unified ShowConfigWindow function (TwShowDaqConfig, TwShowBasicTiming etc now use this function)
//arguments: index of config page to show
//
//possible return values:   TwSuccess 				success
//							TwError					error
////////////////////////////////////////////////////////////////////////////////
//TOFWERK_DAQ_API TwRetVal TwGetPressure(double* pressure, int channel);
//read out of one channel from maxiGauge
//arguments: double (by ref.) to hold pressure, channel (1-6)
//
//possible return values:   TwSuccess 				success
//                        	TwDaqRecNotRunning	   	if no recording application is found
//							TwOutOfBounds			if channel is not in range 1-6
//							TwNoData				if no connected maxiGauge was found
//							TwError					error (e.g. NULL pointer)
////////////////////////////////////////////////////////////////////////////////
//TOFWERK_DAQ_API double TwTof2Mass(double TofSample, int MassCalibMode, double* p);
//TOFWERK_DAQ_API double TwMass2Tof(double Mass, int MassCalibMode, double* p);
//helper functions to convert between samples and mass
//arguments: mass/sample, mass calibration parameter a, mass calibration parameter b
//
//possible return values:   mass/sample (TwMass2Tof returns 0.0 for negative masses)
////////////////////////////////////////////////////////////////////////////////
TOFWERK_DAQ_API TwRetVal TwWaitForEndOfAcquisition(int timeout);
//waits for the end of the current acquisition (end of last run for multi-run acquisitions)
//arguments: timeout in ms
//
//possible return values:   TwSuccess 				success (acquisition finished)
//							TwNoActiveAcquisition   if there is no active acquisition
//                        	TwDaqRecNotRunning	   	if no recording application is found
//							TwTimeout				acquisition did not finish during timeout
////////////////////////////////////////////////////////////////////////////////
TOFWERK_DAQ_API TwRetVal TwReleaseSharedMemory(void);
//manually release the shared memory acquisition buffers (only needed in special cases,see TwGetSharedMemory())
//arguments:
//
//possible return values:   TwSuccess 				success
//                        	TwDaqRecNotRunning	   	if no recording application is found
//                        	TwError				   	if unmapping of shared memory failed
////////////////////////////////////////////////////////////////////////////////
TOFWERK_DAQ_API TwRetVal TwTpsConnect(void);
//connect to remote control enabled TPSController software
//arguments:
//
//possible return values:   TwSuccess 				success
//                        	TwError				   	if connect failed (TPS Controller not running or not in RC mode)
////////////////////////////////////////////////////////////////////////////////
TOFWERK_DAQ_API TwRetVal TwTpsConnect2(char* ip, int type);
//connect to new TPS
//arguments:
//
//possible return values:   TwSuccess 				success
//                        	TwError				   	if connect failed (TPS Controller not running or not in RC mode)
////////////////////////////////////////////////////////////////////////////////
TOFWERK_DAQ_API TwRetVal TwTpsDisconnect(void);
//disconnect from remote control enabled TPSController software
//arguments:
//
//possible return values:   TwSuccess 				success
//                        	TwError				   	if disconnect failed (TPS Controller not connected)
////////////////////////////////////////////////////////////////////////////////
TOFWERK_DAQ_API TwRetVal TwTpsGetMonitorValue(int moduleCode, double* value);
//retrieves the last reported monitor value
//arguments: moduleCode, reference to double to store value
//
//possible return values:   TwSuccess 				success
//                        	TwError				   	if not connected
//							TwOutOfBounds			if module does not exist
////////////////////////////////////////////////////////////////////////////////
TOFWERK_DAQ_API TwRetVal TwTpsGetTargetValue(int moduleCode, double* value);
//retrieves the last reported target value
//arguments: moduleCode, reference to double to store value
//
//possible return values:   TwSuccess 				success
//                        	TwError				   	if not connected
//							TwOutOfBounds			if module does not exist
////////////////////////////////////////////////////////////////////////////////
TOFWERK_DAQ_API TwRetVal TwTpsGetLastSetValue(int moduleCode, double* value);
//retrieves the last reported "last set" value
//arguments: moduleCode, reference to double to store value
//
//possible return values:   TwSuccess 				success
//                        	TwError				   	if not connected
//							TwOutOfBounds			if module does not exist
////////////////////////////////////////////////////////////////////////////////
TOFWERK_DAQ_API TwRetVal TwTpsSetTargetValue(int moduleCode, double value);
//retrieves the last reported target value
//arguments: moduleCode, reference to double to store value
//
//possible return values:   TwSuccess 				success
//                        	TwError				   	if not connected
//							TwOutOfBounds			if module does not exist
////////////////////////////////////////////////////////////////////////////////
TOFWERK_DAQ_API TwRetVal TwTpsGetNbrModules(int* nbrModules);
//get nbr of modules
//arguments: reference to int nbrModules
//
//possible return values:   TwSuccess 				success
//                        	TwError				   	if not connected
////////////////////////////////////////////////////////////////////////////////
TOFWERK_DAQ_API TwRetVal TwTpsGetModuleCodes(int* moduleCodeBuffer, int bufferLength);
//get module codes od all tps modules
//arguments: pointer to int buffer of length bufferLength
//
//possible return values:   TwSuccess 				success
//                        	TwError				   	if not connected
////////////////////////////////////////////////////////////////////////////////
TOFWERK_DAQ_API TwRetVal TwTpsInitialize();
//get module codes od all tps modules
//arguments:
//
//possible return values:   TwSuccess 				success
//                        	TwError				   	if not connected
////////////////////////////////////////////////////////////////////////////////
TOFWERK_DAQ_API TwRetVal TwTpsSetAllVoltages();
//get module codes od all tps modules
//arguments:
//
//possible return values:   TwSuccess 				success
//                        	TwError				   	if not connected
////////////////////////////////////////////////////////////////////////////////
TOFWERK_DAQ_API TwRetVal TwTpsShutdown();
//get module codes od all tps modules
//arguments:
//
//possible return values:   TwSuccess 				success
//                        	TwError				   	if not connected
////////////////////////////////////////////////////////////////////////////////
TOFWERK_DAQ_API TwRetVal TwTpsGetStatus(int* status);
//get module codes od all tps modules
//arguments: status by reference
//
//possible return values:   TwSuccess 				success
//                        	TwError				   	if not connected
////////////////////////////////////////////////////////////////////////////////
TOFWERK_DAQ_API double TwGetDllVersion(void);
//returns the DLL version
//arguments: none
//possible return values: 	version number (double)
////////////////////////////////////////////////////////////////////////////////
TOFWERK_DAQ_API TwRetVal TwGetMassCalib(int* mode, int* nbrParams, double* p, int* nbrPoints, double* mass, double* tof, double* weight);
TOFWERK_DAQ_API TwRetVal TwGetMassCalibEx(int* mode, int* nbrParams, double* p, int* nbrPoints, double* mass, double* tof, double* weight, char* label);
TOFWERK_DAQ_API TwRetVal TwGetMassCalib2(int* mode, int* nbrParams, double* p, int* nbrPoints, double* mass, double* tof, double* weight);
TOFWERK_DAQ_API TwRetVal TwGetMassCalib2Ex(int* mode, int* nbrParams, double* p, int* nbrPoints, double* mass, double* tof, double* weight, char* label);
//get the current mass claibration parameters and/or calibration points
//arguments: MCa, MCb pointers to doubles (if either is NULL, mass calibration parameters are not returned)
//           nbrPoints pointer to int holding the number of alements for buffers mass and tof, when the function
//           returns nbrPoints is the actual number of calibration points copied to the buffers mass and tof.
//			 mass, tof pointer to double buffers.
//possible return values:	TwDaqRecNotRunning		if no recording application is found
//							TwSuccess   			if new mass calib parameters were transmitted to Recorder, nbrPoints gives the actual number of points in mass/tof
//							TwValueAdjusted			if there are more calibration points than nbrPoints, only the first nbrPoints mass/tof values are copied to buffers.
//							TwTimeout				if no confirmation was received within TW_TIMEOUT ms
//							TwError					if nbrPoints > 0 and mass or tof == NULL
////////////////////////////////////////////////////////////////////////////////
TOFWERK_DAQ_API TwRetVal TwContinueAcquisition(void);
//signals to TofDaqRecorder to continue the current acquisition
//arguments: none
//possible return values:	TwDaqRecNotRunning		if no recording application is found
//							TwNoActiveAcquisition	if no active acquisition is detected
//							TwSuccess   			if DAQ continues
//							TwTimeout				if no confirmation was received within TW_TIMEOUT ms
//							TwError					if DAQ is not ready to accept a continue event
////////////////////////////////////////////////////////////////////////////////
TOFWERK_DAQ_API bool TwManualContinueNeeded(void);
//indicates whether TofDaqRec expects a continue event (see TwContinueAcquisition)
//arguments: none
//possible return values: 	true					TofDaq needs a continue event
//							false					no action needed
////////////////////////////////////////////////////////////////////////////////
TOFWERK_DAQ_API TwRetVal TwGetSumSpectrumFromShMem(double* Spectrum, bool Normalize);
TOFWERK_DAQ_API TwRetVal TwGetSumSpectrumFromShMem2(double* Spectrum, bool Normalize);
//copies the contents of the sum spectrum buffer to the array Spectrum. If normalize is false the buffer contents are transferred
//as sum, else the values are normalized with the number of tof extractions contributing to the sum.
//arguments: double array of minimal length NbrSamples, normalize flag
//possible return values:	TwDaqRecNotRunning		if no recording application is found
//							TwNoData   				if no data was found
//							TwSuccess   			if spectrum was successfully copied
//							TwError					if Spectrum is a NULL pointer
////////////////////////////////////////////////////////////////////////////////
TOFWERK_DAQ_API TwRetVal TwGetTofSpectrumFromShMem(float* Spectrum, int SegmentIndex, int SegmentEndIndex, int BufIndex, bool Normalize);
TOFWERK_DAQ_API TwRetVal TwGetTofSpectrumFromShMem2(float* Spectrum, int SegmentIndex, int SegmentEndIndex, int BufIndex, bool Normalize);
//copies a single mass spectrum to the array Spectrum. If normalize is false the buffer contents are transferred
//as sum, else the values are normalized with the number of tof extractions contributing to the spectrum.
//arguments: float array of minimal length NbrSamples, segment index, buf index and  normalize flag
//possible return values:	TwDaqRecNotRunning		if no recording application is found
//							TwNoData   				if no data was found
//							TwSuccess   			if spectrum was successfully copied
//							TwError					if Spectrum is a NULL pointer
//							TwOutOfBounds			if segment or buf index is out of range
////////////////////////////////////////////////////////////////////////////////
TOFWERK_DAQ_API TwRetVal TwGetStickSpectrumFromShMem(float* Spectrum, float* Masses, int SegmentIndex, int SegmentEndIndex, int BufIndex);
TOFWERK_DAQ_API TwRetVal TwGetStickSpectrumFromShMem2(float* Spectrum, float* Masses, int SegmentIndex, int SegmentEndIndex, int BufIndex);
//copies a single stick spectrum to the array Spectrum.
//arguments: float arrays of minimal length NbrPeaks for intensities and masses (can be NULL if not needed), segment index and buf index.
//possible return values:	TwDaqRecNotRunning		if no recording application is found
//							TwNoData   				if no data was found
//							TwSuccess   			if spectrum was successfully copied
//							TwError					if Spectrum is a NULL pointer
//							TwOutOfBounds			if segment or buf index is out of range
////////////////////////////////////////////////////////////////////////////////
TOFWERK_DAQ_API TwRetVal TwGetSpecXaxisFromShMem(double* SpecAxis, int Type, char* UnitLabel, double maxMass);
//copies a X axis value for every sample to SpecAxis. available types: 0=sample index, 1=mass, 2=TOF, 3= freq
//arguments: double array of min. length NbrSamples, integer for axis type, char buffer for unit label (max. 63 char)
//possible return values:	TwDaqRecNotRunning		if no recording application is found
//							TwSuccess   			if axis was successfully copied
//							TwError					if SpecAxis is a NULL pointer
//							TwOutOfBounds			if Type is out of valid range
////////////////////////////////////////////////////////////////////////////////
TOFWERK_DAQ_API TwRetVal TwGetSegmentProfileFromShMem(float* SegmentProfile, int PeakIndex, int BufIndex);
TOFWERK_DAQ_API TwRetVal TwGetSegmentProfileFromShMem2(float* SegmentProfile, int PeakIndex, int BufIndex);
//copies the segment profile for a given peak and buf index to the float array.
//arguments: float array with at least NbrPeaks elements, peak index, buf index
//			 pass -1 for peak index in order to get segment profiles for all peaks. SegmentProfile must point
//			 to at least NbrPeaks*NbrSegments floats.
//possible return values:	TwDaqRecNotRunning		if no recording application is found
//							TwNoData   				if no data was found
//							TwSuccess   			if data was successfully copied
//							TwError					if SegmentProfile is a NULL pointer
//							TwOutOfBounds			if PeakIndex or BufIndex is out of range
////////////////////////////////////////////////////////////////////////////////
TOFWERK_DAQ_API TwRetVal TwGetBufTimeFromShMem(double* BufTime, int BufIndex, int WriteIndex);
//copies the time stamp for a given buf and write to BufTime.
//arguments: double (by ref.), buf index, write index
//possible return values:	TwDaqRecNotRunning		if no recording application is found
//							TwNoData   				if no data was found
//							TwSuccess   			if data was successfully copied
//							TwError					if BufTime is a NULL pointer
//							TwOutOfBounds			if PeakIndex or BufIndex is out of range
////////////////////////////////////////////////////////////////////////////////
TOFWERK_DAQ_API TwRetVal TwAddUserData(char* Location, int NbrElements, char* ElementDescription, double* Data, int CompressionLevel);
//creates datasets "Data" and "Info" at group location, if the datasets exist, a new line is added
//arguments: string giving the group where datasets are created, int with number elements (per call to this function),
//char buffer of size 256*NbrElements, pointer to NbrElements doubles with the actual data
//possible return values:	TwDaqRecNotRunning		if no recording application is found
//							TwNoActiveAcquisition   if there is no active acquisition
//							TwSuccess   			if data was successfully written to file
//							TwFileNotFound			if the current acquisition writes no HDF5 file
//							TwError					general error
//							TwTimeout				if adding the data times out
//							TwOutOfBounds			if the dataset at location exists already but has different NbrElements
////////////////////////////////////////////////////////////////////////////////
TOFWERK_DAQ_API TwRetVal TwAddUserDataMultiRow(char* Location, int NbrElements, int NbrRows, char* ElementDescription, double* Data, int CompressionLevel);
//creates datasets "Data" and "Info" at group location, if the datasets exist, a new line is added
//arguments: string giving the group where datasets are created, int with number elements (per call to this function), number of rows to add
//char buffer of size 256*NbrElements, pointer to NbrElements doubles with the actual data
//possible return values:	TwDaqRecNotRunning		if no recording application is found
//							TwNoActiveAcquisition   if there is no active acquisition
//							TwSuccess   			if data was successfully written to file
//							TwFileNotFound			if the current acquisition writes no HDF5 file
//							TwError					general error
//							TwTimeout				if adding the data times out
//							TwOutOfBounds			if the dataset at location exists already but has different NbrElements
////////////////////////////////////////////////////////////////////////////////
TOFWERK_DAQ_API TwRetVal TwRegisterUserDataBuf(char* Location, int NbrElements, char* ElementDescription, int CompressionLevel);
//registers a shared memory region of NbrElements doubles to be written to the data file for every recorded buf.
//Use TwUnregisterUserData and TwUpdateUserData functions to unregister / change values.
//arguments: string giving the group where datasets are created, int with number elements in the shared memory,
//char buffer of size 256*NbrElements, integer compression level (0: no compression, 1-9: zlib compression level)
//possible return values:	TwDaqRecNotRunning		if no recording application is found
//							TwActiveAcquisition     if there is an active acquisition
//							TwSuccess   			if register was successfull
//							TwError					general error
////////////////////////////////////////////////////////////////////////////////
TOFWERK_DAQ_API TwRetVal TwRegisterUserDataWrite(char* Location, int NbrElements, char* ElementDescription, int CompressionLevel);
//registers a shared memory region of NbrElements doubles to be written to the data file for every recorded write.
//Use TwUnregisterUserData and TwUpdateUserData functions to unregister / change values.
//arguments: string giving the group where datasets are created, int with number elements in the shared memory,
//char buffer of size 256*NbrElements, integer compression level (0: no compression, 1-9: zlib compression level)
//possible return values:	TwDaqRecNotRunning		if no recording application is found
//							TwActiveAcquisition     if there is an active acquisition
//							TwSuccess   			if register was successfull
//							TwError					general error
////////////////////////////////////////////////////////////////////////////////
TOFWERK_DAQ_API TwRetVal TwRegisterUserDataNoStore(char* Location, int NbrElements, char* ElementDescription);
//registers a shared memory region of NbrElements doubles to be written to the data file for every recorded write.
//Use TwUnregisterUserData and TwUpdateUserData functions to unregister / change values.
//arguments: string giving the group where datasets are created, int with number elements in the shared memory,
//char buffer of size 256*NbrElements
//possible return values:	TwDaqRecNotRunning		if no recording application is found
//							TwActiveAcquisition     if there is an active acquisition
//							TwSuccess   			if register was successfull
//							TwError					general error
////////////////////////////////////////////////////////////////////////////////
TOFWERK_DAQ_API TwRetVal TwUnregisterUserData(char* Location);
//unregisters a shared memory region set up before with TwRegisterUserDataBuf or TwRegisterUserDataWrite functions.
//arguments: string giving the group that was used to register the shared memory
//possible return values:	TwDaqRecNotRunning		if no recording application is found
//							TwActiveAcquisition     if there is an active acquisition
//							TwSuccess   			if unregister was successfull
//							TwError					general error
////////////////////////////////////////////////////////////////////////////////
TOFWERK_DAQ_API TwRetVal TwUpdateUserData(char* Location, int NbrElements, double* Data);
//updates values in a shared memory region registered with TwRegisterUserDataBuf or TwRegisterUserDataWrite functions.
//arguments: string giving the group that was used to register the shared memory, number elements and a pointer to the actual data (at least NbrElements long).
//possible return values:	TwDaqRecNotRunning		if no recording application is found
//							TwSuccess   			if update was successfull
//							TwError					general error
////////////////////////////////////////////////////////////////////////////////
TOFWERK_DAQ_API TwRetVal TwReadRegUserData(char* Location, int NbrElements, double* Data);
//updates values in a shared memory region registered with TwRegisterUserDataBuf or TwRegisterUserDataWrite functions.
//arguments: string giving the group that was used to register the shared memory, number elements and a pointer to the actual data (at least NbrElements long).
//possible return values:	TwDaqRecNotRunning		if no recording application is found
//							TwSuccess   			if read was successfull
//							TwError					general error
////////////////////////////////////////////////////////////////////////////////
TOFWERK_DAQ_API TwRetVal TwQueryRegUserDataSize(char* Location, int* NbrElements);
////////////////////////////////////////////////////////////////////////////////
TOFWERK_DAQ_API TwRetVal TwIssueDio4Pulse(int delay, int width);
////////////////////////////////////////////////////////////////////////////////
TOFWERK_DAQ_API TwRetVal TwSetDio4State(int state);
////////////////////////////////////////////////////////////////////////////////
TOFWERK_DAQ_API TwRetVal TwGetRegUserDataSources(int* arrayLength, char* location, int* nbrElements, int* type);
//queries currently registered data sources. If location == nbrElements == type == NULL, *arrayLength holds the
//total number of registered data sources on function return.
////////////////////////////////////////////////////////////////////////////////
TOFWERK_DAQ_API TwRetVal TwGetRegUserDataDesc(char* location, int* nbrElements, char* elementDescription);
//queries descriptions of a registered data source.
////////////////////////////////////////////////////////////////////////////////
TOFWERK_DAQ_API bool TwDioStartDelayActive(void);
//returns true if TofDaqRec.exe received start command but a start delay is configured, false otherwise
//arguments: none
////////////////////////////////////////////////////////////////////////////////
TOFWERK_DAQ_API TwRetVal TwInitializeDaqDevice(void);
//(re)initializes the current DAQ board
//arguments: none
////////////////////////////////////////////////////////////////////////////////
TOFWERK_DAQ_API TwRetVal TwTpsGetActiveFilament(int* activeFilament);
////////////////////////////////////////////////////////////////////////////////
TOFWERK_DAQ_API TwRetVal TwTpsSetActiveFilament(int activeFilament);
////////////////////////////////////////////////////////////////////////////////
TOFWERK_DAQ_API void TwSetTimeout(int timeout);
////////////////////////////////////////////////////////////////////////////////
TOFWERK_DAQ_API int TwGetTimeout(void);
////////////////////////////////////////////////////////////////////////////////
TOFWERK_DAQ_API TwRetVal TwTpsLoadSetFile(char* setFile);
////////////////////////////////////////////////////////////////////////////////
TOFWERK_DAQ_API TwRetVal TwTpsSaveSetFile(char* setFile);
////////////////////////////////////////////////////////////////////////////////
TOFWERK_DAQ_API TwRetVal TwTpsIonModeShutdown();
////////////////////////////////////////////////////////////////////////////////
TOFWERK_DAQ_API TwRetVal TwTpsGetModuleLimits(int moduleCode, double* minLimit, double* maxLimit);
////////////////////////////////////////////////////////////////////////////////
TOFWERK_DAQ_API TwRetVal TwAutoSetupDaqDevice(void);
////////////////////////////////////////////////////////////////////////////////
TOFWERK_DAQ_API TwRetVal TwOnDemandMassCalibration(int action);
////////////////////////////////////////////////////////////////////////////////
TOFWERK_DAQ_API TwRetVal TwKeepFileOpen(bool keepOpen);
////////////////////////////////////////////////////////////////////////////////
TOFWERK_DAQ_API TwRetVal TwTpsChangeIonMode(int ionMode);
////////////////////////////////////////////////////////////////////////////////

#ifdef __cplusplus
}
#endif

#endif

