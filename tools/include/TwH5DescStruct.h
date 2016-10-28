#ifndef TwH5DescStructH
#define TwH5DescStructH

#ifdef  __BORLANDC__
#include <stdint.h>
#else
#include <inttypes.h>
#endif

//structure containing the descriptor for Tofwerk HDF5 files
typedef struct
	{
	int32_t nbrSamples;
	int32_t nbrPeaks;
	int32_t nbrWaveforms;
	int32_t nbrSegments;
	int32_t nbrBlocks;
	int32_t nbrMemories;
	int32_t nbrBufs;
	int32_t nbrWrites;
	int32_t nbrLogEntries;
	bool secondTof;
	bool hasSumSpectrum;
	bool hasSumSpectrum2;
	bool hasBufTimes;
	bool hasTofData;
	bool hasTofData2;
	bool hasPeakData;
	bool hasPeakData2;
	bool hasTpsData;
	bool hasNbrMemories;
	bool hasPressureData;
	bool hasLogData;
	bool hasMassCalibData;
	bool hasMassCalib2Data;
	bool hasCh1RawData;
	bool hasCh2RawData;
	bool hasRawDataDesc;
	bool hasEventList;
	int16_t segIlf;
	int32_t eventListMaxElementLength;
	int32_t daqMode;
	int32_t acquisitionMode;
	int32_t massCalibMode;
	int32_t massCalibMode2;
	int32_t nbrCalibParams;
	int32_t nbrCalibParams2;
	int32_t nbrCubes; //should be logically between nbrBufs and nbrWrites, but here it does not change the overall descriptor structure size
	double p[16];
	double p2[16];
	double tofPeriod;
	double blockPeriod;
	float sampleInterval;
	float singleIonSignal;
	float singleIonSignal2;
	char dummy[4];
	} TwH5Desc;

#ifdef _WIN32
typedef bool (__cdecl TwProgressCallback)(double progress);
typedef bool (__cdecl TwProgressCallback2)(double progress, void* userData);
#else
typedef bool(TwProgressCallback)(double progress) __attribute__((cdecl));
typedef bool(TwProgressCallback2)(double progress, void* userData) __attribute__((cdecl));
#endif


#endif