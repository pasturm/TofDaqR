#ifndef TwToolDllH
#define TwToolDllH

#ifdef _WIN32
#ifdef TWTOOLDLL_EXPORTS
#define TOFWERK_TOOL_API __declspec(dllexport)
#else
#define TOFWERK_TOOL_API __declspec(dllimport)
#endif
#else
#ifdef TWTOOLDLL_EXPORTS
#define TOFWERK_TOOL_API __attribute__((visibility("default")))
#else
#define TOFWERK_TOOL_API 
#endif	
#endif


#if defined(_WIN32) && !defined(__BORLANDC__)

#define TwBruteForceCalibrate          _TwBruteForceCalibrate
#define TwDecomposeMass                _TwDecomposeMass
#define TwEncImsCleanup                _TwEncImsCleanup
#define TwEncImsCorrelateMultiProfiles _TwEncImsCorrelateMultiProfiles
#define TwEncImsCorrelateProfile      _TwEncImsCorrelateProfile
#define TwEvalMultiPeak               _TwEvalMultiPeak
#define TwEvalResolution              _TwEvalResolution
#define TwEvalSinglePeak              _TwEvalSinglePeak
#define TwFitResolution               _TwFitResolution
#define TwFitSinglePeak               _TwFitSinglePeak
#define TwFitSinglePeak2              _TwFitSinglePeak2
#define TwGetComposition              _TwGetComposition
#define TwGetIsotopePattern           _TwGetIsotopePattern
#define TwGetIsotopePattern2          _TwGetIsotopePattern2
#define TwGetMassCalibInfo            _TwGetMassCalibInfo
#define TwGetMoleculeMass             _TwGetMoleculeMass
#define TwMass2Tof                    _TwMass2Tof
#define TwMassCalibrate               _TwMassCalibrate
#define TwMultiPeakFit                _TwMultiPeakFit
#define TwSiCleanup                   _TwSiCleanup
#define TwSiEvalPhd                   _TwSiEvalPhd
#define TwSiFitPhd                    _TwSiFitPhd
#define TwSiFitRateFromPhd            _TwSiFitRateFromPhd
#define TwSiGetHistogram              _TwSiGetHistogram
#define TwSiGetHistogramAmp           _TwSiGetHistogramAmp
#define TwSiGetSumHistogram           _TwSiGetSumHistogram
#define TwSiGetSumHistogramAmp        _TwSiGetSumHistogramAmp
#define TwSiInitializeHistograms      _TwSiInitializeHistograms
#define TwSiProcessSpectrum           _TwSiProcessSpectrum
#define TwSiResetHistograms           _TwSiResetHistograms
#define TwSiSetProcessingOptions      _TwSiSetProcessingOptions
#define TwTof2Mass                    _TwTof2Mass
#define TwTranslateReturnValue        _TwTranslateReturnValue
#define TwNistLibrarySearch        	  _TwNistLibrarySearch
#define TwNistLibraryQueryResult      _TwNistLibraryQueryResult
#define TwFindTpsIp					  _TwFindTpsIp
#define TwGenerateImageData			  _TwGenerateImageData
#define TwImagePaletteOffset2Value    _TwImagePaletteOffset2Value
#define TwImageValue2PaletteOffset    _TwImageValue2PaletteOffset
#define TwImageGetPaletteRGB		  _TwImageGetPaletteRGB

#endif

#include "TofIpcStrucs.h"  //for TwRetVal

#ifdef __cplusplus
extern "C" {
#endif
//---------------------------------------------------------------------------
//peak fitting
TOFWERK_TOOL_API TwRetVal TwFitSinglePeak(int nbrDataPoints, double* yVals, double* xVals, int peakType,
										  double* blOffset, double* blSlope, double* amplitude, double* fwhmLo,
										  double* fwhmHi, double* peakPos, double* mu);
//---------------------------------------------------------------------------
TOFWERK_TOOL_API TwRetVal TwFitSinglePeak2(int nbrDataPoints, double* yVals, double* xVals, int peakType, double* param);
//---------------------------------------------------------------------------
TOFWERK_TOOL_API double TwEvalSinglePeak(double xVal, double* param);
//---------------------------------------------------------------------------
//TOF resolving power
TOFWERK_TOOL_API TwRetVal TwFitResolution(int nbrPoints, double* mass, double* resolution, double* R0, double* m0, double* dm);
//---------------------------------------------------------------------------
TOFWERK_TOOL_API double TwEvalResolution(double R0, double m0, double dm, double mass);
//---------------------------------------------------------------------------
//chemical utility functions
TOFWERK_TOOL_API TwRetVal TwGetMoleculeMass(char* molecule, double* mass);
//---------------------------------------------------------------------------
TOFWERK_TOOL_API TwRetVal TwGetIsotopePattern(char* molecule, double abundanceLimit, int* nbrIsotopes, double* isoMass, double* isoAbundance);
//---------------------------------------------------------------------------
TOFWERK_TOOL_API TwRetVal TwGetIsotopePattern2(char* molecule, double abundanceLimit, int* nbrIsotopes, double* isoMass, double* isoAbundance);
//---------------------------------------------------------------------------
TOFWERK_TOOL_API TwRetVal TwDecomposeMass(double targetMass, double tolerance, int nbrAtoms, double* atomMass, char* atomLabel,
										  int nbrFilters, int* elementIndex1, int* elementIndex2, double* filterMinVal, double* filterMaxVal,
										  int* nbrCompomers);
//---------------------------------------------------------------------------
TOFWERK_TOOL_API TwRetVal TwGetComposition(int index, char* sumFormula, int* sumFormulaLength, double* mass, double* massError);
//---------------------------------------------------------------------------
//mass calibration
TOFWERK_TOOL_API double TwTof2Mass(double tofSample, int massCalibMode, double* p);
//---------------------------------------------------------------------------
TOFWERK_TOOL_API double TwMass2Tof(double mass, int massCalibMode, double* p);
//---------------------------------------------------------------------------
TOFWERK_TOOL_API TwRetVal TwMassCalibrate(int massCalibMode, int nbrPoints, double* mass, double* tof, double* weight,
										  int* nbrParams, double* p, double* legacyA, double* legacyB);
//---------------------------------------------------------------------------
TOFWERK_TOOL_API TwRetVal TwBruteForceCalibrate(int specLength, double* spectrum, int nbrPeaks, bool extendSearch,
												int triggerDelay, int pulseWidth, double sampleInterval, double* a, double* b);
//---------------------------------------------------------------------------
TOFWERK_TOOL_API TwRetVal TwGetMassCalibInfo(int massCalibMode, char* description, int* nbrParams);
//---------------------------------------------------------------------------
//return value 2 text
TOFWERK_TOOL_API char* TwTranslateReturnValue(TwRetVal returnValue);
//---------------------------------------------------------------------------
//encoded IMS functions
TOFWERK_TOOL_API TwRetVal TwEncImsCorrelateProfile(float* profile, int opMode, int* par);
//---------------------------------------------------------------------------
TOFWERK_TOOL_API TwRetVal TwEncImsCorrelateMultiProfiles(void* data, int dataType, int nbrProfiles, int calcMethod, float correlateThreshold, int opMode, int* par);
//---------------------------------------------------------------------------
TOFWERK_TOOL_API TwRetVal TwEncImsCleanup(void);
//---------------------------------------------------------------------------
//single ion functions
TOFWERK_TOOL_API TwRetVal TwSiFitPhd(int nbrPoints, double* intensity, double* counts, double* fwhm, double* a, double* par);
//---------------------------------------------------------------------------
TOFWERK_TOOL_API double TwSiEvalPhd(double* par, double intensity);
//---------------------------------------------------------------------------
TOFWERK_TOOL_API TwRetVal TwSiFitRateFromPhd(int nbrPoints, double* intensity, double* counts, double* siPar, double* rate, double* fitCounts, int* nbrNIonTraces, double** nIonTrace);
//---------------------------------------------------------------------------
TOFWERK_TOOL_API TwRetVal TwSiInitializeHistograms(int nbrHistograms, double* loMass, double* hiMass, int* specType);
//---------------------------------------------------------------------------
TOFWERK_TOOL_API TwRetVal TwSiSetProcessingOptions(char* option, double value, int specType);
//---------------------------------------------------------------------------
TOFWERK_TOOL_API TwRetVal TwSiProcessSpectrum(float* spectrum, int nbrSamples, int specType, float* blFromData, float* thrFromData);
//---------------------------------------------------------------------------
TOFWERK_TOOL_API TwRetVal TwSiGetHistogram(int histogramIndex, float* intensity, unsigned int* counts, unsigned int* arrayLength, unsigned int* spectrumCount, double* meanValue);
//---------------------------------------------------------------------------
TOFWERK_TOOL_API TwRetVal TwSiGetHistogramAmp(int histogramIndex, float* intensity, unsigned int* counts, unsigned int* arrayLength, unsigned int* spectrumCount, double* meanValue);
//---------------------------------------------------------------------------
TOFWERK_TOOL_API TwRetVal TwSiGetSumHistogram(int specType, float* intensity, unsigned int* counts, unsigned int* arrayLength, unsigned int* spectrumCount, double* meanValue, double minMass, double maxMass, double minRate, double maxRate);
//---------------------------------------------------------------------------
TOFWERK_TOOL_API TwRetVal TwSiGetSumHistogramAmp(int specType, float* intensity, unsigned int* counts, unsigned int* arrayLength, unsigned int* spectrumCount, double* meanValue, double minMass, double maxMass, double minRate, double maxRate);
//---------------------------------------------------------------------------
TOFWERK_TOOL_API TwRetVal TwSiResetHistograms(void);
//---------------------------------------------------------------------------
TOFWERK_TOOL_API TwRetVal TwSiCleanup(void);
//---------------------------------------------------------------------------
//multi peak fitting
TOFWERK_TOOL_API TwRetVal TwMultiPeakFit(int nbrDataPoints, double* dataX, double* dataY, int nbrPeaks, double* mass, double* intensity, double* commonPar, int options);
//---------------------------------------------------------------------------
TOFWERK_TOOL_API double TwEvalMultiPeak(double x, int nbrPeaks, double* mass, double* intensity, double* commonPar);
//---------------------------------------------------------------------------
//NIST library lookup
#ifdef _WIN32
TOFWERK_TOOL_API TwRetVal TwNistLibrarySearch(int nbrSticks, double* stickMass, double* stickIntensity, bool interactive, unsigned int* nbrResults);
//---------------------------------------------------------------------------
TOFWERK_TOOL_API TwRetVal TwNistLibraryQueryResult(int index, int propertyIndex, int* valueLen, char* value);
//---------------------------------------------------------------------------
TOFWERK_TOOL_API TwRetVal TwFindTpsIp(char* tpsSerial, int timeout, int* hostStrLen, char* hostStr);
//---------------------------------------------------------------------------
#endif
//---------------------------------------------------------------------------
TOFWERK_TOOL_API TwRetVal TwGenerateImageData(int width, int height, float* dataIn, int palette, unsigned char* paletteOffset, bool* paletteInvert, float* dataMinMax, float* gammaVal, unsigned char* specialColours, unsigned char* imageOut);
//---------------------------------------------------------------------------
TOFWERK_TOOL_API double TwImagePaletteOffset2Value(unsigned char pixelValue, double min, double max, double gamma);
//---------------------------------------------------------------------------
TOFWERK_TOOL_API unsigned char TwImageValue2PaletteOffset(double dataValue, double min, double max, double gamma);
//---------------------------------------------------------------------------
TOFWERK_TOOL_API TwRetVal TwImageGetPaletteRGB(int palette, unsigned char paletteValue, unsigned char offset, bool invert, unsigned char* rgb);
//---------------------------------------------------------------------------

#ifdef __cplusplus
}
#endif

#endif


