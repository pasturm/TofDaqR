#ifndef TwH5DllH
#define TwH5DllH

// #ifdef _WIN32
#if defined(_WIN32) && !defined(__GNUC__)
#ifdef TWH5DLL_EXPORTS
#define TOFWERK_H5_API __declspec(dllexport)
#else
#define TOFWERK_H5_API __declspec(dllimport)
#endif
#endif
#ifdef __GNUC__
#ifdef TWH5DLL_EXPORTS
#define TOFWERK_H5_API __attribute__((visibility("default")))
#else
#define TOFWERK_H5_API 
#endif
#endif

// #if defined(_WIN32) && defined(_MSC_VER)
#if defined(_WIN32) && (defined(_MSC_VER) || defined(__GNUC__))

#define TwChangePeakDataInit					_TwChangePeakDataInit
#define TwChangePeakDataWrite					_TwChangePeakDataWrite
#define TwChangePeakDataFinalize				_TwChangePeakDataFinalize
#define TwChangePeakTable						_TwChangePeakTable
#define TwChangePeakTable2						_TwChangePeakTable2
#define TwChangePeakTableFromFile				_TwChangePeakTableFromFile
#define TwChangePeakTableFromFile2				_TwChangePeakTableFromFile2
#define TwCloseAll								_TwCloseAll
#define TwCloseH5								_TwCloseH5
#define TwFreeEventListData						_TwFreeEventListData
#define TwFreeEventListData2					_TwFreeEventListData2
#define TwGetAcquisitionLogFromH5				_TwGetAcquisitionLogFromH5
#define TwGetBufProfile2FromH5					_TwGetBufProfile2FromH5
#define TwGetBufProfileFromH5					_TwGetBufProfileFromH5
#define TwGetBufTimeFromH5						_TwGetBufTimeFromH5
#define TwGetBufWriteProfile2FromH5				_TwGetBufWriteProfile2FromH5
#define TwGetBufWriteProfile2FromH5_2			_TwGetBufWriteProfile2FromH5_2
#define TwGetBufWriteProfileFromH5				_TwGetBufWriteProfileFromH5
#define TwGetBufWriteProfileFromH5_2			_TwGetBufWriteProfileFromH5_2
#define TwGetDoubleAttributeFromH5				_TwGetDoubleAttributeFromH5
#define TwGetEventListBlobFromH5				_TwGetEventListBlobFromH5
#define TwGetEventListDataFromH5				_TwGetEventListDataFromH5
#define TwGetEventListSpectrumFromH5			_TwGetEventListSpectrumFromH5
#define TwGetFloatAttributeFromH5				_TwGetFloatAttributeFromH5
#define TwGetH5Descriptor						_TwGetH5Descriptor
#define TwGetInt64AttributeFromH5				_TwGetInt64AttributeFromH5
#define TwGetIntAttributeFromH5					_TwGetIntAttributeFromH5
#define TwGetPeakData							_TwGetPeakData
#define TwGetPeakData2							_TwGetPeakData2
#define TwGetPeakParametersFromH5				_TwGetPeakParametersFromH5
#define TwGetRegUserDataFromH5					_TwGetRegUserDataFromH5
#define TwGetSegmentProfile2FromH5				_TwGetSegmentProfile2FromH5
#define TwGetSegmentProfileFromH5				_TwGetSegmentProfileFromH5
#define TwGetSpecXaxisFromH5					_TwGetSpecXaxisFromH5
#define TwGetStickSpectrum2FromH5				_TwGetStickSpectrum2FromH5
#define TwGetStickSpectrumFromH5				_TwGetStickSpectrumFromH5
#define TwGetStringAttributeFromH5				_TwGetStringAttributeFromH5
#define TwGetSumSpectrum2FromH5					_TwGetSumSpectrum2FromH5
#define TwGetSumSpectrumFromH5					_TwGetSumSpectrumFromH5
#define TwGetTimingData							_TwGetTimingData
#define TwGetTofData							_TwGetTofData
#define TwGetTofData2							_TwGetTofData2
#define TwGetTofSpectrum2FromH5					_TwGetTofSpectrum2FromH5
#define TwGetTofSpectrumFromH5					_TwGetTofSpectrumFromH5
#define TwGetUint64AttributeFromH5				_TwGetUint64AttributeFromH5
#define TwGetUintAttributeFromH5				_TwGetUintAttributeFromH5
#define TwGetUserDataFromH5						_TwGetUserDataFromH5
#define TwGetWriteProfile2FromH5				_TwGetWriteProfile2FromH5
#define TwGetWriteProfileFromH5					_TwGetWriteProfileFromH5
#define TwH5GetMassCalibPar						_TwH5GetMassCalibPar
#define TwH5SetMassCalib						_TwH5SetMassCalib
#define TwH5SetMassCalibEx						_TwH5SetMassCalibEx
#define TwH5SetMassCalib2						_TwH5SetMassCalib2
#define TwH5SetMassCalib2Ex						_TwH5SetMassCalib2Ex
#define TwH5TranslateReturnValue				_TwH5TranslateReturnValue
#define TwMultiPeakFitIntegration				_TwMultiPeakFitIntegration
#define TwReadRawData							_TwReadRawData
#define TwSetDoubleAttributeInH5				_TwSetDoubleAttributeInH5
#define TwSetFloatAttributeInH5					_TwSetFloatAttributeInH5
#define TwSetInt64AttributeInH5					_TwSetInt64AttributeInH5
#define TwSetIntAttributeInH5					_TwSetIntAttributeInH5
#define TwSetStringAttributeInH5				_TwSetStringAttributeInH5
#define TwSetUint64AttributeInH5				_TwSetUint64AttributeInH5
#define TwSetUintAttributeInH5					_TwSetUintAttributeInH5
#define TwH5AddLogEntry         				_TwH5AddLogEntry
#define TwH5AddUserDataMultiRow					_TwH5AddUserDataMultiRow
#define TwH5SetMassCalibDynamic					_TwH5SetMassCalibDynamic
#define TwGenerateSegmentProfilesFromEventList	_TwGenerateSegmentProfilesFromEventList
#define TwGetRegUserDataSourcesFromH5           _TwGetRegUserDataSourcesFromH5
#define TwH5MakePaletteImage					_TwH5MakePaletteImage
#define TwH5MakeTrueColorImage					_TwH5MakeTrueColorImage

#endif

#include "TwH5DescStruct.h"
#include "TofIpcStrucs.h"

#ifdef __cplusplus
extern "C" {
#endif
////////////////////////////////////////////////////////////////////////////////
//TOFWERK_H5_API double TwH5Tof2Mass(double TofSample, int MassCalibMode, double* par);
//TOFWERK_H5_API double TwH5Mass2Tof(double Mass, int MassCalibMode, double* par);
//helper functions to convert between samples and mass
//arguments: mass/sample, mass calibration parameter a, mass calibration parameter b
//
//possible return values:   mass/sample (TwMass2Tof returns 0.0 for negative masses)
////////////////////////////////////////////////////////////////////////////////
TOFWERK_H5_API TwRetVal TwGetH5Descriptor(char* Filename, TwH5Desc* descriptor);
//returns a descriptor structure for the Tofwerk HDF5 file
//arguments: char* with the filename, HDF5 descriptor structure (by ref)
//possible return values:	TwFileNotFound			if the data file was not found (or access denied)
//							TwSuccess   			if data was successfully copied
//							TwError					if descriptor is a NULL pointer
////////////////////////////////////////////////////////////////////////////////
TOFWERK_H5_API TwRetVal TwCloseH5(char* Filename);
//closes an open HDF5 file (alternative to TwCleanUpDll() which closes all open files)
//arguments: char* with the filename, HDF5 descriptor structure (by ref)
//possible return values:	TwFileNotFound			if the data file was not open
//							TwSuccess   			if the file was closed
////////////////////////////////////////////////////////////////////////////////
TOFWERK_H5_API TwRetVal TwCloseAll(void);
//closes all open HDF5 file
//arguments: none
//possible return values:	TwSuccess   			if all files were closed
////////////////////////////////////////////////////////////////////////////////
TOFWERK_H5_API TwRetVal TwGetSumSpectrumFromH5(char* Filename, double* Spectrum, bool Normalize);
TOFWERK_H5_API TwRetVal TwGetSumSpectrum2FromH5(char* Filename, double* Spectrum, bool Normalize);
//gets the sum spectrum for a given HDF5 file. If normalize is false the sum spectrum is read
//as sum, else the values are normalized with the number of tof extractions contributing to the sum.
//arguments: char* with the filename, double array of minimal length NbrSamples, normalize flag
//possible return values:	TwFileNotFound			if the data file was not found (or access denied)
//							TwNoData   				if no data was found
//							TwSuccess   			if spectrum was successfully copied
//							TwError					if Spectrum is a NULL pointer
////////////////////////////////////////////////////////////////////////////////
TOFWERK_H5_API TwRetVal TwGetTofSpectrumFromH5(char* Filename, float* Spectrum, int SegmentIndex, int SegmentEndIndex,
																				 int BufIndex, int BufEndIndex,
																				 int WriteIndex, int WriteEndIndex,
																				 bool BufWriteLinked, bool Normalize);
TOFWERK_H5_API TwRetVal TwGetTofSpectrum2FromH5(char* Filename, float* Spectrum, int SegmentIndex, int SegmentEndIndex,
																				 int BufIndex, int BufEndIndex,
																				 int WriteIndex, int WriteEndIndex,
																				 bool BufWriteLinked, bool Normalize);
//copies a single mass spectrum to the array Spectrum. If normalize is false the contents are transferred
//as sum, else the values are normalized with the number of tof extractions contributing to the spectrum.
//arguments: char* with the filename, float array of minimal length NbrSamples, segment index, segment end index,
//buf index, buf end index, write index, write end index, bufs/writes linked flag and normalize flag
//possible return values:	TwFileNotFound			if the data file was not found (or access denied)
//							TwNoData   				if no data was found
//							TwSuccess   			if spectrum was successfully copied
//							TwError					if Spectrum is a NULL pointer
//							TwOutOfBounds			if any of the indices is out of bounds
////////////////////////////////////////////////////////////////////////////////
TOFWERK_H5_API TwRetVal TwGetStickSpectrumFromH5(char* Filename, float* Spectrum, int SegmentIndex, int SegmentEndIndex,
																				   int BufIndex, int BufEndIndex,
																				   int WriteIndex, int WriteEndIndex,
																				   bool BufWriteLinked, bool Normalize);
TOFWERK_H5_API TwRetVal TwGetStickSpectrum2FromH5(char* Filename, float* Spectrum, int SegmentIndex, int SegmentEndIndex,
																				   int BufIndex, int BufEndIndex,
																				   int WriteIndex, int WriteEndIndex,
																				   bool BufWriteLinked, bool Normalize);
//copies a single stick spectrum to the array Spectrum. If normalize is false the contents are transferred
//as sum, else the values are normalized with the number of tof extractions contributing to the spectrum.
//arguments: char* with the filename, float array of minimal length NbrPeaks, segment index, segment end index,
//buf index, buf end index, write index, write end index bufs/writes linked flag and normalize flag
//possible return values:	TwFileNotFound			if the data file was not found (or access denied)
//							TwNoData   				if no data was found
//							TwSuccess   			if spectrum was successfully copied
//							TwError					if Spectrum is a NULL pointer
//							TwOutOfBounds			if any of the indices is out of bounds
////////////////////////////////////////////////////////////////////////////////
TOFWERK_H5_API TwRetVal TwGetPeakParametersFromH5(char* Filename, TPeakPar* PeakPar, int PeakIndex);
//retrieves peak parameters from HDF5 file
//arguments: pointer to user allocated TPeakPar structure (see TofIpcStrucs.h), index of peak of interest
//possible return values:	TwFileNotFound			if the data file was not found (or access denied)
//							TwSuccess   			if PeakPar structure was successfully transferred
//							TwOutOfBounds  			if PeakIndex is not in range 0 to NbrPeaks-1
//							TwError					if peak parameters could not be retrieved or PeakPar=NULL
////////////////////////////////////////////////////////////////////////////////
TOFWERK_H5_API TwRetVal TwGetBufTimeFromH5(char* Filename, double* BufTime, int BufIndex, int WriteIndex);
//copies the time stamp for a given buf and write to BufTime.
//arguments: char* with the filename, double (by ref.), buf index, write index
//possible return values:	TwFileNotFound			if the data file was not found (or access denied)
//							TwSuccess   			if data was successfully copied
//							TwError					if BufTime is a NULL pointer
//							TwOutOfBounds			if BufIndex or WriteIndex is out of range
////////////////////////////////////////////////////////////////////////////////
TOFWERK_H5_API TwRetVal TwGetSpecXaxisFromH5(char* Filename, double* SpecAxis, int Type, char* UnitLabel, double maxMass, int writeIndex);
//copies a X axis value for every sample to SpecAxis. available types: 0=sample index, 1=mass, 2=TOF, 3= freq
//arguments: char* with the filename,double array of min. length NbrSamples, integer for axis type, char buffer for unit label (max. 63 char)
//possible return values:	TwFileNotFound			if the data file was not found (or access denied)
//							TwSuccess   			if axis was successfully copied
//							TwError					if SpecAxis is a NULL pointer
//							TwOutOfBounds			if Type is out of valid range
////////////////////////////////////////////////////////////////////////////////
TOFWERK_H5_API TwRetVal TwGetSegmentProfileFromH5(char* Filename, float* SegmentProfile, int PeakIndex, int BufStartIndex, int BufEndIndex, int WriteStartIndex, int WriteEndIndex, bool BufWriteLinked);
TOFWERK_H5_API TwRetVal TwGetSegmentProfile2FromH5(char* Filename, float* SegmentProfile, int PeakIndex, int BufStartIndex, int BufEndIndex, int WriteStartIndex, int WriteEndIndex, bool BufWriteLinked);
//copies the segment profile for a given peak, buf and write index to the float array.
//arguments: char* with the filename, float array with at least NbrPeaks elements, peak index, buf index, write index
//			 pass -1 for peak index in order to get segment profiles for all peaks. SegmentProfile must point
//			 to at least NbrPeaks*NbrSegments floats.
//possible return values:	TwFileNotFound			if the data file was not found (or access denied)
//							TwSuccess   			if data was successfully copied
//							TwError					if SegmentProfile is a NULL pointer
//							TwOutOfBounds			if PeakIndex, BufIndex and/or WriteIndex is out of range
////////////////////////////////////////////////////////////////////////////////
TOFWERK_H5_API TwRetVal TwGetBufWriteProfileFromH5(char* Filename, float* Profile, int PeakIndex, int SegmentStartIndex, int SegmentEndIndex);
TOFWERK_H5_API TwRetVal TwGetBufWriteProfileFromH5_2(char* Filename, float* Profile, int PeakIndex, int SegmentStartIndex, int SegmentEndIndex, TwProgressCallback* callback, TwProgressCallback2* callback2, void* cb2UserData);
TOFWERK_H5_API TwRetVal TwGetBufWriteProfile2FromH5(char* Filename, float* Profile, int PeakIndex, int SegmentStartIndex, int SegmentEndIndex);
TOFWERK_H5_API TwRetVal TwGetBufWriteProfile2FromH5_2(char* Filename, float* Profile, int PeakIndex, int SegmentStartIndex, int SegmentEndIndex, TwProgressCallback* callback, TwProgressCallback2* callback2, void* cb2UserData);
//copies the buf/write profile (linked) for a given peak and segment index to the float array.
//arguments: char* with the filename, float array with at least NbrBufs*NbrWrites elements, peak index, buf index, write index
//			 pass -1 for peak index in order to get buf/write profiles for all peaks. Profile must point
//			 to at least NbrPeaks*NbrBufs*NbrWrites floats.
//possible return values:	TwFileNotFound			if the data file was not found (or access denied)
//							TwSuccess   			if data was successfully copied
//							TwError					if Profile is a NULL pointer
//							TwOutOfBounds			if PeakIndexand/or SegmentIndex is out of range
////////////////////////////////////////////////////////////////////////////////
TOFWERK_H5_API TwRetVal TwGetBufProfileFromH5(char* Filename, float* Profile, int PeakIndex, int SegmentStartIndex, int SegmentEndIndex, int WriteStartIndex, int WriteEndIndex);
TOFWERK_H5_API TwRetVal TwGetBufProfile2FromH5(char* Filename, float* Profile, int PeakIndex, int SegmentStartIndex, int SegmentEndIndex, int WriteStartIndex, int WriteEndIndex);
////////////////////////////////////////////////////////////////////////////////
TOFWERK_H5_API TwRetVal TwGetWriteProfileFromH5(char* Filename, float* Profile, int PeakIndex, int SegmentStartIndex, int SegmentEndIndex, int BufStartIndex, int BufEndIndex);
TOFWERK_H5_API TwRetVal TwGetWriteProfile2FromH5(char* Filename, float* Profile, int PeakIndex, int SegmentStartIndex, int SegmentEndIndex, int BufStartIndex, int BufEndIndex);
////////////////////////////////////////////////////////////////////////////////
//lower level functions that allow direct dataset access
//most TwGet...FromH5 functions use these
////////////////////////////////////////////////////////////////////////////////
TOFWERK_H5_API TwRetVal TwGetTofData(char* Filename, int sampleOffset, int sampleCount,
									  int segOffset, int segCount,
									  int bufOffset, int bufCount, int writeOffset, int writeCount,
									  float* dataBuffer);
TOFWERK_H5_API TwRetVal TwGetTofData2(char* Filename, int sampleOffset, int sampleCount,
									  int segOffset, int segCount,
									  int bufOffset, int bufCount, int writeOffset, int writeCount,
									  float* dataBuffer);

TOFWERK_H5_API TwRetVal TwGetTimingData(char* Filename, int bufOffset, int bufCount,
										 int writeOffset, int writeCount, double* dataBuffer);

TOFWERK_H5_API TwRetVal TwGetPeakData(char* Filename, int peakOffset, int peakCount,
									   int segOffset, int segCount,
									   int bufOffset, int bufCount, int writeOffset, int writeCount,
									   float* dataBuffer);
TOFWERK_H5_API TwRetVal TwGetPeakData2(char* Filename, int peakOffset, int peakCount,
									   int segOffset, int segCount,
									   int bufOffset, int bufCount, int writeOffset, int writeCount,
									   float* dataBuffer);
////////////////////////////////////////////////////////////////////////////////
//attribute functions
////////////////////////////////////////////////////////////////////////////////
TOFWERK_H5_API TwRetVal TwGetIntAttributeFromH5(char* Filename, char* location, char* name, int* attribute);
TOFWERK_H5_API TwRetVal TwGetUintAttributeFromH5(char* Filename, char* location, char* name, unsigned int* attribute);
TOFWERK_H5_API TwRetVal TwGetInt64AttributeFromH5(char* Filename, char* location, char* name, int64_t* attribute);
TOFWERK_H5_API TwRetVal TwGetUint64AttributeFromH5(char* Filename, char* location, char* name, uint64_t* attribute);
TOFWERK_H5_API TwRetVal TwGetFloatAttributeFromH5(char* Filename, char* location, char* name, float* attribute);
TOFWERK_H5_API TwRetVal TwGetDoubleAttributeFromH5(char* Filename, char* location, char* name, double* attribute);
TOFWERK_H5_API TwRetVal TwGetStringAttributeFromH5(char* Filename, char* location, char* name, char* attribute);
//read signle attributes from the HDF5 file
//arguments: char* for filename, attribute location and attribute name, attribute by reference of corresponding type
//possible return values:	TwFileNotFound			if the data file was not found (or access denied)
//							TwSuccess   			if data was successfully copied
//							TwError					if attribute is a NULL pointer or attribute not found in file
////////////////////////////////////////////////////////////////////////////////
TOFWERK_H5_API TwRetVal TwSetIntAttributeInH5(char* Filename, char* location, char* name, int attribute);
TOFWERK_H5_API TwRetVal TwSetUintAttributeInH5(char* Filename, char* location, char* name, unsigned int attribute);
TOFWERK_H5_API TwRetVal TwSetInt64AttributeInH5(char* Filename, char* location, char* name, int64_t attribute);
TOFWERK_H5_API TwRetVal TwSetUint64AttributeInH5(char* Filename, char* location, char* name, uint64_t attribute);
TOFWERK_H5_API TwRetVal TwSetFloatAttributeInH5(char* Filename, char* location, char* name, float attribute);
TOFWERK_H5_API TwRetVal TwSetDoubleAttributeInH5(char* Filename, char* location, char* name, double attribute);
TOFWERK_H5_API TwRetVal TwSetStringAttributeInH5(char* Filename, char* location, char* name, char* attribute);
//write single attribute to the HDF5 file
//arguments: char* for filename, attribute location and attribute name, attribute by reference of corresponding type
//possible return values:	TwFileNotFound			if the data file was not found
//							TwSuccess   			if attribute was successfully added
//							TwError					if attribute was not written
////////////////////////////////////////////////////////////////////////////////
TOFWERK_H5_API TwRetVal TwGetRegUserDataFromH5(char* Filename, char* location, int bufIndex, int writeIndex,
											   int* bufLength, double* buffer, char* description);
//reads an entry from a registred data source (created by TwRegisterUserData... functions)
//arguments: char* for filename, char* for location in HDF5 file, bufIndex and writeIndex of interest, bufferLength (by reference),
//           pointer to data buffer (must point to at least bufLength doubles), char* to store descriptions (minimal length is 256*buflength, can be NULL)
//possible return values:    TwFileNotFound			if the data file was not found
//                           TwSuccess   			if data was copied successfully
//							 TwValueAdjusted		if bufLength does not match data dimensions (bufLength will be adjusted to actual data dimensions)
//							 TwOutOfBounds			if buf or write index is out of range
//							 TwError 				if location doesn't hold user data or buffer is NULL
////////////////////////////////////////////////////////////////////////////////
TOFWERK_H5_API TwRetVal TwChangePeakDataInit(char* Filename, TPeakPar* newPeakPar, int nbrNewPeakPar, int options);
////////////////////////////////////////////////////////////////////////////////
TOFWERK_H5_API TwRetVal TwChangePeakDataWrite(char* Filename, int peakOffset, int peakCount, int segOffset, int segCount, int bufOffset, int bufCount, int writeOffset, int writeCount, float* data, float* data2);
////////////////////////////////////////////////////////////////////////////////
TOFWERK_H5_API TwRetVal TwChangePeakDataFinalize(char* Filename);
////////////////////////////////////////////////////////////////////////////////
TOFWERK_H5_API TwRetVal TwChangePeakTable(char* Filename, TPeakPar* newPeakPar, int nbrNewPeakPar, int compressionLevel, TwProgressCallback* callback);
//creates new PeakTable and (re)calculates PeakData accordingly
////////////////////////////////////////////////////////////////////////////////
TOFWERK_H5_API TwRetVal TwChangePeakTableFromFile(char* Filename, char* massTable, int compressionLevel, TwProgressCallback* callback);
//creates new PeakTable and (re)calculates PeakData accordingly
////////////////////////////////////////////////////////////////////////////////
TOFWERK_H5_API TwRetVal TwChangePeakTable2(char* Filename, TPeakPar* newPeakPar, int nbrNewPeakPar, int compressionLevel, TwProgressCallback2* callback, void* cbUserData);
//creates new PeakTable and (re)calculates PeakData accordingly
////////////////////////////////////////////////////////////////////////////////
TOFWERK_H5_API TwRetVal TwChangePeakTableFromFile2(char* Filename, char* massTable, int compressionLevel, TwProgressCallback2* callback, void* cbUserData);
//creates new PeakTable and (re)calculates PeakData accordingly
////////////////////////////////////////////////////////////////////////////////
TOFWERK_H5_API TwRetVal TwReadRawData (char* Filename, int channel, int bufIndex, int writeIndex, int* bufferSize, void* buffer);
//reads raw data from a data file
////////////////////////////////////////////////////////////////////////////////
TOFWERK_H5_API char* TwH5TranslateReturnValue(TwRetVal ReturnValue);
//translates a TwRetVal into a character string
//arguments: TwRetVal from any TofDaq API function
////////////////////////////////////////////////////////////////////////////////
TOFWERK_H5_API TwRetVal TwH5SetMassCalib (char* Filename, int mode, int nbrParams, double* p, int nbrPoints, double* mass, double* tof, double* weight);
TOFWERK_H5_API TwRetVal TwH5SetMassCalibEx (char* Filename, int mode, int nbrParams, double* p, int nbrPoints, double* mass, double* tof, double* weight, char* label);
TOFWERK_H5_API TwRetVal TwH5SetMassCalib2 (char* Filename, int mode, int nbrParams, double* p, int nbrPoints, double* mass, double* tof, double* weight);
TOFWERK_H5_API TwRetVal TwH5SetMassCalib2Ex (char* Filename, int mode, int nbrParams, double* p, int nbrPoints, double* mass, double* tof, double* weight, char* label);
//changes global mass calibration in data file
////////////////////////////////////////////////////////////////////////////////
TOFWERK_H5_API TwRetVal TwGetUserDataFromH5(char* filename, char* location, int* rowIndex, int* nbrElements, double* buffer, char* elementDescription);
//reads user data (stored with TwAddUserData) from the file
////////////////////////////////////////////////////////////////////////////////
TOFWERK_H5_API TwRetVal TwGetAcquisitionLogFromH5(char* filename, int index, int64_t* timestamp, char* logText);
////////////////////////////////////////////////////////////////////////////////
TOFWERK_H5_API TwRetVal TwGetEventListSpectrumFromH5(char* filename, int segmentIndex, int bufIndex, int writeIndex, int* bufferSize, unsigned int* buffer);
////////////////////////////////////////////////////////////////////////////////
TOFWERK_H5_API TwRetVal TwGetEventListDataFromH5(char* filename, int segStartIndex, int segEndIndex, int bufStartIndex, int bufEndIndex,
												 int writeStartIndex, int writeEndIndex, unsigned int*** dataBuffer, unsigned int** dataLength);
////////////////////////////////////////////////////////////////////////////////
TOFWERK_H5_API TwRetVal TwGetEventListBlobFromH5(char* filename, int** blob, int blobLength, unsigned int*** dataBuffer, unsigned int** dataLength, int* totalSpectraInBlobs);
////////////////////////////////////////////////////////////////////////////////
TOFWERK_H5_API TwRetVal TwFreeEventListData();
////////////////////////////////////////////////////////////////////////////////
TOFWERK_H5_API TwRetVal TwFreeEventListData2(char* filename);
////////////////////////////////////////////////////////////////////////////////
TOFWERK_H5_API TwRetVal TwH5GetMassCalibPar(char* Filename, int segmentIndex, int bufIndex, int writeIndex, int* mode, int* nbrParams, double* p);
////////////////////////////////////////////////////////////////////////////////
TOFWERK_H5_API TwRetVal TwMultiPeakFitIntegration(char* Filename, int option, TwProgressCallback2* callback, void* cbUserData);
////////////////////////////////////////////////////////////////////////////////
TOFWERK_H5_API TwRetVal TwH5AddLogEntry(char* filename, char* logEntryText, uint64_t logEntryTime);
////////////////////////////////////////////////////////////////////////////////
TOFWERK_H5_API TwRetVal TwH5AddUserDataMultiRow(char* filename, char* location, int nbrElements, int nbrRows, char* elementDescription, double* data, int compressionLevel);
////////////////////////////////////////////////////////////////////////////////
TOFWERK_H5_API TwRetVal TwH5SetMassCalibDynamic(char* Filename, int segmentIndex, int bufIndex, int writeIndex, int nbrParams, double* par, int nbrStat, double* stat);
////////////////////////////////////////////////////////////////////////////////
TOFWERK_H5_API TwRetVal TwGenerateSegmentProfilesFromEventList(char* Filename, int nbrProfiles, double* startMass, double* endMass, int bufStartIndex, int bufEndIndex, int writeStartIndex, int writeEndIndex, float* data, bool startEndInSamples);
////////////////////////////////////////////////////////////////////////////////
TOFWERK_H5_API TwRetVal TwGetRegUserDataSourcesFromH5(char* Filename, int* nbrSources, char* sourceLocation, int* sourceLength, bool* hasDesc, int* type);
////////////////////////////////////////////////////////////////////////////////
TOFWERK_H5_API TwRetVal TwH5MakePaletteImage(char* filename, char* location, float* data, int width, int height, int palette, unsigned char paletteOffset, bool paletteInvert, float* dataMinMax, float gammaVal);
////////////////////////////////////////////////////////////////////////////////
TOFWERK_H5_API TwRetVal TwH5MakeTrueColorImage(char* filename, char* location, float* data, int width, int height, unsigned char* paletteOffset, bool* paletteInvert, float* dataMinMax, float* gammaVal, unsigned char* specialColours);
////////////////////////////////////////////////////////////////////////////////
#ifdef __cplusplus
}
#endif

#endif

