#ifndef TofIpcStrucsH
#define TofIpcStrucsH

#ifdef  __BORLANDC__
#include <stdint.h>
#else
#include <inttypes.h>
#endif


////////////////////////////////////////////////////////////////////////////////
//most TofDaqDll functions return one of these values. Use TwTranslateReturnValue
//function to get a char string of the return value.
typedef enum{ TwDaqRecNotRunning,
				TwAcquisitionActive,
				TwNoActiveAcquisition,
				TwFileNotFound,
				TwSuccess,
				TwError,
				TwOutOfBounds,
				TwNoData,
				TwTimeout,
				TwValueAdjusted,
				TwInvalidParameter,
				TwInvalidValue,
				TwAborted
				} TwRetVal;
////////////////////////////////////////////////////////////////////////////////


//descriptor structure containing all information to access data in shared memory
//and fields dedicated to interprocess communication

typedef struct
{
int32_t   NbrSamples;
int32_t	  NbrRawSamples;
int32_t   NbrPeaks;
int32_t   NbrWaveforms;
int32_t   NbrSegments;
int32_t   NbrBlocks;
int32_t   NbrMemories;
int32_t   NbrBufs;
int32_t   NbrWrites;
int32_t   NbrRuns;

int32_t   iWaveform;
int32_t   iSegment;
int32_t	  iBlock;
int32_t   iMemory;
int32_t   iBuf;
int32_t   iWrite;
int32_t   iRun;

int32_t   TotalBufsRecorded;
int32_t   TotalBufsProcessed;
int32_t   TotalBufsWritten;
int32_t   OverallBufsProcessed;
int32_t   TotalNbrMemories;
int32_t	  TotalMemoriesProcessed;

uint32_t RawDataRecordedBuf1;
uint32_t RawDataRecordedBuf2;
uint32_t RawDataLastElementInBuffer1;
uint32_t RawDataLastElementInBuffer2;
uint32_t RawDataProcessedBuf1;
uint32_t RawDataProcessedBuf2;
uint32_t RawDataWrittenBuf1;
uint32_t RawDataWrittenBuf2;

float SampleInterval;
int32_t   TofPeriod;
int32_t NbrCubes; //logically should be positioned between NbrBufs and NbrWrites, however due to data alignment rules, positioning here does not change the size of the descriptor structure
int64_t BlockPeriod;
int64_t BlockPulseDelay;
int64_t BlockDelay;
float SingleIonSignal;
float SingleIonSignal2;
int32_t MassCalibMode;
int32_t MassCalibMode2;
int32_t NbrMassCalibParams;
int32_t NbrMassCalibParams2;
double p[16];
double p2[16];
float R0;
float dm;
float m0;
bool SecondTof;

char  chIniFileName[256];
char  CurrentDataFileName[256];
uint8_t segIlf;
uint16_t iCube;   //logically should be positioned between iBuf and iWrite (and int), however due to data alignment rules, positioning here does not change the size of the descriptor structure
int32_t DaqMode;
int32_t AcquisitionMode;
int32_t CombineMode;
int32_t RecalibFreq;
char AcquisitionLogText[256];  		   //this is also used to pass string attributes
uint64_t AcquisitionLogTime;
uint64_t TimeZero;

uint32_t ExternalLock;        			   //points to NbrBuf bools (when DaqActive!)
uint32_t ProcessingLevel;
//add atributes
int32_t AttributeType; 					   // 0=int, 1=double, 2=string
char AttributeObject[256];
char AttributeName[128];
int32_t AttributeInt;       			   //int attributes are passed here
double AttributeDouble; 			   //double attributes are passed here
//var NbrMemories config
bool EnableVarNbrMemories;
char dummy3[3];
int32_t NbrSteps;
int32_t CurrentStepAtBuf;
int32_t NbrMemoriesForCurrentStep;
} TSharedMemoryDesc;

//Peak parameter structure
typedef struct
{
char label[64];    //peak label
float mass;        //exact mass
float loMass;      //lower integration boundary (in mass units)
float hiMass;      //upper integration boundary (in mass units)
} TPeakPar;

////Pressure parameter structure
//typedef struct
//{
//char channelName[5]; 	//4 chars for name + /0
//char sensorType[10]; 	//9 chars for type + /0
//float pressure;
//int status;
//} TPressureGauge;

#endif
