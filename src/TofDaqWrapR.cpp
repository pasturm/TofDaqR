#include <R.h>
#include <Rdefines.h>
#include <TofIpcStrucs.h>
#include <TwToolDll.h>
#include <TwH5Dll.h>
#ifdef _WIN32
#include <TofDaqDll.h>
#endif
#include <cmath>        // std::abs, pow, sqrt, etc.

#ifdef BUILD_DLL
#define DLL_EXPORT __attribute__((visibility("default")))
#else
#define DLL_EXPORT
#endif


// function which gets a C-string form a R string
char* RtoCstring(SEXP rstring)
{
  // Protection not needed if rstring is an R function argument, otherwise use PROTECT(string = AS_CHARACTER(string));
  int stringlength = strlen(CHAR(STRING_ELT(rstring, 0)));
  char *cstring;
  cstring = R_alloc(stringlength + 1, sizeof(char));  // stringlength + 1 to account for null termination
  strcpy(cstring, CHAR(STRING_ELT(rstring, 0)));  // copy rstring to cstring
  cstring[stringlength] = '\0';  // null terminate for safety
  return(cstring);
}

#ifdef __cplusplus
extern "C" {
#endif

  // TwH5Dll ------------------------------------------------------------------
  DLL_EXPORT SEXP GetH5Descriptor(SEXP filename)
  {
    //get the filename in a useable form
    char *cFilename = RtoCstring(filename);

    //get descriptor of file
    TwH5Desc desc;
    if (TwGetH5Descriptor(cFilename, &desc) != TwSuccess)
    {
      TwCloseH5(cFilename);
      return(R_NilValue);
    }
    TwCloseH5(cFilename);

    int descLength = 40;
    char const* names[40] = { "nbrSamples",
                        "nbrPeaks",
                        "nbrWaveforms",
                        "nbrSegments",
                        "nbrBlocks",
                        "nbrMemories",
                        "nbrBufs",
                        "nbrWrites",
                        "nbrLogEntries",
                        "hasSumSpectrum",
                        "hasSumSpectrum2",
                        "hasBufTimes",
                        "hasTofData",
                        "hasTofData2",
                        "hasPeakData",
                        "hasPeakData2",
                        "hasTpsData",
                        "hasNbrMemories",
                        "hasPressureData",
                        "hasLogData",
                        "hasMassCalibData",
                        "hasMassCalib2Data",
                        "hasCh1RawData",
                        "hasCh2RawData",
                        "hasRawDataDesc",
                        "hasEventList",
                        "eventListMaxElementLength",
                        "daqMode",
                        "acquisitionMode",
                        "massCalibMode",
                        "massCalibMode2",
                        "nbrCalibParams",
                        "nbrCalibParams2",
                        "p",
                        "p2",
                        "tofPeriod",
                        "blockPeriod",
                        "sampleInterval",
                        "singleIonSignal",
                        "singleIonSignal2"
    };

    int *p_rdesc;
    double *pn_rdesc;
    SEXP rdesc, list, list_names;

    // Creating a list with 40 vector elements:
    PROTECT(list = allocVector(VECSXP, descLength));

    // creating descriptor elements
    PROTECT(rdesc = NEW_INTEGER(1));
    p_rdesc = INTEGER_POINTER(rdesc);
    p_rdesc[0] = desc.nbrSamples;
    SET_VECTOR_ELT(list, 0, rdesc);
    UNPROTECT(1);

    PROTECT(rdesc = NEW_INTEGER(1));
    p_rdesc = INTEGER_POINTER(rdesc);
    p_rdesc[0] = desc.nbrPeaks;
    SET_VECTOR_ELT(list, 1, rdesc);
    UNPROTECT(1);

    PROTECT(rdesc = NEW_INTEGER(1));
    p_rdesc = INTEGER_POINTER(rdesc);
    p_rdesc[0] = desc.nbrWaveforms;
    SET_VECTOR_ELT(list, 2, rdesc);
    UNPROTECT(1);

    PROTECT(rdesc = NEW_INTEGER(1));
    p_rdesc = INTEGER_POINTER(rdesc);
    p_rdesc[0] = desc.nbrSegments;
    SET_VECTOR_ELT(list, 3, rdesc);
    UNPROTECT(1);

    PROTECT(rdesc = NEW_INTEGER(1));
    p_rdesc = INTEGER_POINTER(rdesc);
    p_rdesc[0] = desc.nbrBlocks;
    SET_VECTOR_ELT(list, 4, rdesc);
    UNPROTECT(1);

    PROTECT(rdesc = NEW_INTEGER(1));
    p_rdesc = INTEGER_POINTER(rdesc);
    p_rdesc[0] = desc.nbrMemories;
    SET_VECTOR_ELT(list, 5, rdesc);
    UNPROTECT(1);

    PROTECT(rdesc = NEW_INTEGER(1));
    p_rdesc = INTEGER_POINTER(rdesc);
    p_rdesc[0] = desc.nbrBufs;
    SET_VECTOR_ELT(list, 6, rdesc);
    UNPROTECT(1);

    PROTECT(rdesc = NEW_INTEGER(1));
    p_rdesc = INTEGER_POINTER(rdesc);
    p_rdesc[0] = desc.nbrWrites;
    SET_VECTOR_ELT(list, 7, rdesc);
    UNPROTECT(1);

    PROTECT(rdesc = NEW_INTEGER(1));
    p_rdesc = INTEGER_POINTER(rdesc);
    p_rdesc[0] = desc.nbrLogEntries;
    SET_VECTOR_ELT(list, 8, rdesc);
    UNPROTECT(1);

    PROTECT(rdesc = NEW_LOGICAL(1));
    p_rdesc = LOGICAL_POINTER(rdesc);
    p_rdesc[0] = desc.hasSumSpectrum;
    SET_VECTOR_ELT(list, 9, rdesc);
    UNPROTECT(1);

    PROTECT(rdesc = NEW_LOGICAL(1));
    p_rdesc = LOGICAL_POINTER(rdesc);
    p_rdesc[0] = desc.hasSumSpectrum2;
    SET_VECTOR_ELT(list, 10, rdesc);
    UNPROTECT(1);

    PROTECT(rdesc = NEW_LOGICAL(1));
    p_rdesc = LOGICAL_POINTER(rdesc);
    p_rdesc[0] = desc.hasBufTimes;
    SET_VECTOR_ELT(list, 11, rdesc);
    UNPROTECT(1);

    PROTECT(rdesc = NEW_LOGICAL(1));
    p_rdesc = LOGICAL_POINTER(rdesc);
    p_rdesc[0] = desc.hasTofData;
    SET_VECTOR_ELT(list, 12, rdesc);
    UNPROTECT(1);

    PROTECT(rdesc = NEW_LOGICAL(1));
    p_rdesc = LOGICAL_POINTER(rdesc);
    p_rdesc[0] = desc.hasTofData2;
    SET_VECTOR_ELT(list, 13, rdesc);
    UNPROTECT(1);

    PROTECT(rdesc = NEW_LOGICAL(1));
    p_rdesc = LOGICAL_POINTER(rdesc);
    p_rdesc[0] = desc.hasPeakData;
    SET_VECTOR_ELT(list, 14, rdesc);
    UNPROTECT(1);

    PROTECT(rdesc = NEW_LOGICAL(1));
    p_rdesc = LOGICAL_POINTER(rdesc);
    p_rdesc[0] = desc.hasPeakData2;
    SET_VECTOR_ELT(list, 15, rdesc);
    UNPROTECT(1);

    PROTECT(rdesc = NEW_LOGICAL(1));
    p_rdesc = LOGICAL_POINTER(rdesc);
    p_rdesc[0] = desc.hasTpsData;
    SET_VECTOR_ELT(list, 16, rdesc);
    UNPROTECT(1);

    PROTECT(rdesc = NEW_LOGICAL(1));
    p_rdesc = LOGICAL_POINTER(rdesc);
    p_rdesc[0] = desc.hasNbrMemories;
    SET_VECTOR_ELT(list, 17, rdesc);
    UNPROTECT(1);

    PROTECT(rdesc = NEW_LOGICAL(1));
    p_rdesc = LOGICAL_POINTER(rdesc);
    p_rdesc[0] = desc.hasPressureData;
    SET_VECTOR_ELT(list, 18, rdesc);
    UNPROTECT(1);

    PROTECT(rdesc = NEW_LOGICAL(1));
    p_rdesc = LOGICAL_POINTER(rdesc);
    p_rdesc[0] = desc.hasLogData;
    SET_VECTOR_ELT(list, 19, rdesc);
    UNPROTECT(1);

    PROTECT(rdesc = NEW_LOGICAL(1));
    p_rdesc = LOGICAL_POINTER(rdesc);
    p_rdesc[0] = desc.hasMassCalibData;
    SET_VECTOR_ELT(list, 20, rdesc);
    UNPROTECT(1);

    PROTECT(rdesc = NEW_LOGICAL(1));
    p_rdesc = LOGICAL_POINTER(rdesc);
    p_rdesc[0] = desc.hasMassCalib2Data;
    SET_VECTOR_ELT(list, 21, rdesc);
    UNPROTECT(1);

    PROTECT(rdesc = NEW_LOGICAL(1));
    p_rdesc = LOGICAL_POINTER(rdesc);
    p_rdesc[0] = desc.hasCh1RawData;
    SET_VECTOR_ELT(list, 22, rdesc);
    UNPROTECT(1);

    PROTECT(rdesc = NEW_LOGICAL(1));
    p_rdesc = LOGICAL_POINTER(rdesc);
    p_rdesc[0] = desc.hasCh2RawData;
    SET_VECTOR_ELT(list, 23, rdesc);
    UNPROTECT(1);

    PROTECT(rdesc = NEW_LOGICAL(1));
    p_rdesc = LOGICAL_POINTER(rdesc);
    p_rdesc[0] = desc.hasRawDataDesc;
    SET_VECTOR_ELT(list, 24, rdesc);
    UNPROTECT(1);

    PROTECT(rdesc = NEW_LOGICAL(1));
    p_rdesc = LOGICAL_POINTER(rdesc);
    p_rdesc[0] = desc.hasEventList;
    SET_VECTOR_ELT(list, 25, rdesc);
    UNPROTECT(1);

    PROTECT(rdesc = NEW_INTEGER(1));
    p_rdesc = INTEGER_POINTER(rdesc);
    p_rdesc[0] = desc.eventListMaxElementLength;
    SET_VECTOR_ELT(list, 26, rdesc);
    UNPROTECT(1);

    PROTECT(rdesc = NEW_INTEGER(1));
    p_rdesc = INTEGER_POINTER(rdesc);
    p_rdesc[0] = desc.daqMode;
    SET_VECTOR_ELT(list, 27, rdesc);
    UNPROTECT(1);

    PROTECT(rdesc = NEW_INTEGER(1));
    p_rdesc = INTEGER_POINTER(rdesc);
    p_rdesc[0] = desc.acquisitionMode;
    SET_VECTOR_ELT(list, 28, rdesc);
    UNPROTECT(1);

    PROTECT(rdesc = NEW_INTEGER(1));
    p_rdesc = INTEGER_POINTER(rdesc);
    p_rdesc[0] = desc.massCalibMode;
    SET_VECTOR_ELT(list, 29, rdesc);
    UNPROTECT(1);

    PROTECT(rdesc = NEW_INTEGER(1));
    p_rdesc = INTEGER_POINTER(rdesc);
    p_rdesc[0] = desc.massCalibMode2;
    SET_VECTOR_ELT(list, 30, rdesc);
    UNPROTECT(1);

    PROTECT(rdesc = NEW_NUMERIC(1));
    pn_rdesc = NUMERIC_POINTER(rdesc);
    pn_rdesc[0] = desc.nbrCalibParams;
    SET_VECTOR_ELT(list, 31, rdesc);
    UNPROTECT(1);

    PROTECT(rdesc = NEW_NUMERIC(1));
    pn_rdesc = NUMERIC_POINTER(rdesc);
    pn_rdesc[0] = desc.nbrCalibParams2;
    SET_VECTOR_ELT(list, 32, rdesc);
    UNPROTECT(1);

    PROTECT(rdesc = NEW_NUMERIC(16));
    pn_rdesc = NUMERIC_POINTER(rdesc);
    for (int i = 0; i < 16; ++i)
    {
      if (desc.p[i] < 1e-200)
      {
        pn_rdesc[i] = 0;
      }
      else
      {
        pn_rdesc[i] = desc.p[i];
      }
    }
    SET_VECTOR_ELT(list, 33, rdesc);
    UNPROTECT(1);

    PROTECT(rdesc = NEW_NUMERIC(16));
    pn_rdesc = NUMERIC_POINTER(rdesc);
    for (int i = 0; i < 16; ++i)
    {
      if (desc.p2[i] < 1e-200)
      {
        pn_rdesc[i] = 0;
      }
      else
      {
        pn_rdesc[i] = desc.p2[i];
      }
    }
    SET_VECTOR_ELT(list, 34, rdesc);
    UNPROTECT(1);

    PROTECT(rdesc = NEW_NUMERIC(1));
    pn_rdesc = NUMERIC_POINTER(rdesc);
    pn_rdesc[0] = desc.tofPeriod;
    SET_VECTOR_ELT(list, 35, rdesc);
    UNPROTECT(1);

    PROTECT(rdesc = NEW_NUMERIC(1));
    pn_rdesc = NUMERIC_POINTER(rdesc);
    pn_rdesc[0] = desc.blockPeriod;
    SET_VECTOR_ELT(list, 36, rdesc);
    UNPROTECT(1);

    PROTECT(rdesc = NEW_NUMERIC(1));
    pn_rdesc = NUMERIC_POINTER(rdesc);
    pn_rdesc[0] = desc.sampleInterval;
    SET_VECTOR_ELT(list, 37, rdesc);
    UNPROTECT(1);

    PROTECT(rdesc = NEW_NUMERIC(1));
    pn_rdesc = NUMERIC_POINTER(rdesc);
    pn_rdesc[0] = desc.singleIonSignal;
    SET_VECTOR_ELT(list, 38, rdesc);
    UNPROTECT(1);

    PROTECT(rdesc = NEW_NUMERIC(1));
    pn_rdesc = NUMERIC_POINTER(rdesc);
    if (desc.singleIonSignal2 < 1e-30)
    {
      pn_rdesc[0] = 0;
    }
    else
    {
      pn_rdesc[0] = desc.singleIonSignal2;
    }
    SET_VECTOR_ELT(list, 39, rdesc);
    UNPROTECT(1);

    // Creating a character string vector
    // of the "names" attribute of the
    // objects in out list:
    PROTECT(list_names = allocVector(STRSXP, descLength));
    int i;
    for (i = 0; i < descLength; i++)
      SET_STRING_ELT(list_names, i, mkChar(names[i]));

    // and attaching the vector names:
    setAttrib(list, R_NamesSymbol, list_names);
    UNPROTECT(2);
    return list;
  }

  DLL_EXPORT SEXP CloseH5(SEXP filename)
  {
    //get the filename in a useable form
    char *cFilename = RtoCstring(filename);

    TwRetVal rv = TwCloseH5(cFilename);
    SEXP rValue;
    PROTECT(rValue = NEW_INTEGER(1));
    int *pValue = INTEGER_POINTER(rValue);
    pValue[0] = rv;
    UNPROTECT(1);

    return rValue;
  }

  DLL_EXPORT SEXP CloseAll()
  {
    TwRetVal rv = TwCloseAll();
    SEXP rValue;
    PROTECT(rValue = NEW_INTEGER(1));
    int *pValue = INTEGER_POINTER(rValue);
    pValue[0] = rv;
    UNPROTECT(1);

    return rValue;
  }

  DLL_EXPORT SEXP GetSumSpectrumFromH5(SEXP filename, SEXP normalize)
  {
    bool norm = LOGICAL_VALUE(normalize);

    //get the filename in a useable form
    char *cFilename = RtoCstring(filename);

    //get descriptor of file
    TwH5Desc desc;
    if (TwGetH5Descriptor(cFilename, &desc) != TwSuccess)
    {
      TwCloseH5(cFilename);
      return(R_NilValue);
    }

    //prepare return array
    SEXP sumSpec;
    double *sumSpecBuf;
    int sumSpecLen = desc.nbrSamples;
    PROTECT(sumSpec = NEW_NUMERIC(sumSpecLen));
    sumSpecBuf = NUMERIC_POINTER(sumSpec);
    //read the actual values directly into sumSpecBuf
    if (TwGetSumSpectrumFromH5(cFilename, sumSpecBuf, norm) != TwSuccess)
    {
      TwCloseH5(cFilename);
      UNPROTECT(1);
      return(R_NilValue);
    }
    TwCloseH5(cFilename);
    UNPROTECT(1);
    return sumSpec;
  }

  DLL_EXPORT SEXP GetTofSpectrumFromH5(SEXP filename, SEXP segStart, SEXP segEnd, SEXP bufStart, SEXP bufEnd, SEXP writeStart, SEXP writeEnd, SEXP linked, SEXP normalize)
  {
    bool norm = LOGICAL_VALUE(normalize);
    bool linkBufWrite = LOGICAL_VALUE(linked);
    int segStartIndex = INTEGER_VALUE(segStart);
    int segEndIndex = INTEGER_VALUE(segEnd);
    int bufStartIndex = INTEGER_VALUE(bufStart);
    int bufEndIndex = INTEGER_VALUE(bufEnd);
    int writeStartIndex = INTEGER_VALUE(writeStart);
    int writeEndIndex = INTEGER_VALUE(writeEnd);

    //get the filename in a useable form
    char *cFilename = RtoCstring(filename);

    //get descriptor of file
    TwH5Desc desc;
    if (TwGetH5Descriptor(cFilename, &desc) != TwSuccess)
    {
      TwCloseH5(cFilename);
      return(R_NilValue);
    }
    //prepare return array
    SEXP tofSpec;
    double *tofSpecBuf;
    int specLen = desc.nbrSamples;
    PROTECT(tofSpec = NEW_NUMERIC(specLen));
    tofSpecBuf = NUMERIC_POINTER(tofSpec);
    //allocate temp float buffer
    float *temp = new float[specLen];
    if (TwGetTofSpectrumFromH5(cFilename, temp, segStartIndex, segEndIndex, bufStartIndex, bufEndIndex,
                               writeStartIndex, writeEndIndex, linkBufWrite, norm) != TwSuccess)
    {
      delete[] temp;
      temp = nullptr;
      TwCloseH5(cFilename);
      UNPROTECT(1);
      return(R_NilValue);
    }
    //copy to result
    for (int i = 0; i < specLen; ++i)
    {
      tofSpecBuf[i] = temp[i];
    }
    delete[] temp;
    temp = nullptr;

    TwCloseH5(cFilename);
    UNPROTECT(1);
    return tofSpec;
  }

  DLL_EXPORT SEXP GetStickSpectrumFromH5(SEXP filename, SEXP segStart, SEXP segEnd, SEXP bufStart, SEXP bufEnd, SEXP writeStart, SEXP writeEnd, SEXP linked, SEXP normalize)
  {
    bool norm = LOGICAL_VALUE(normalize);
    bool linkBufWrite = LOGICAL_VALUE(linked);
    int segStartIndex = INTEGER_VALUE(segStart);
    int segEndIndex = INTEGER_VALUE(segEnd);
    int bufStartIndex = INTEGER_VALUE(bufStart);
    int bufEndIndex = INTEGER_VALUE(bufEnd);
    int writeStartIndex = INTEGER_VALUE(writeStart);
    int writeEndIndex = INTEGER_VALUE(writeEnd);

    //get the filename in a useable form
    char *cFilename = RtoCstring(filename);

    //get descriptor of file
    TwH5Desc desc;
    if (TwGetH5Descriptor(cFilename, &desc) != TwSuccess)
    {
      TwCloseH5(cFilename);
      return(R_NilValue);
    }
    //prepare return array
    SEXP stickSpec;
    double *stickSpecBuf;
    int specLen = desc.nbrPeaks;
    PROTECT(stickSpec = NEW_NUMERIC(specLen));
    stickSpecBuf = NUMERIC_POINTER(stickSpec);
    //allocate temp float buffer
    float *temp = new float[specLen];
    if (TwGetStickSpectrumFromH5(cFilename, temp, segStartIndex, segEndIndex, bufStartIndex, bufEndIndex,
                                 writeStartIndex, writeEndIndex, linkBufWrite, norm) != TwSuccess)
    {
      TwCloseH5(cFilename);
      UNPROTECT(1);
      return(R_NilValue);
    }
    //copy to result
    for (int i = 0; i < specLen; ++i)
    {
      stickSpecBuf[i] = temp[i];
    }
    delete[] temp;
    temp = nullptr;

    TwCloseH5(cFilename);
    UNPROTECT(1);
    return stickSpec;
  }

  DLL_EXPORT SEXP GetPeakParametersFromH5(SEXP filename, SEXP peakIndex)
  {
    //get the filename in a useable form
    char *cFilename = RtoCstring(filename);

    int index = INTEGER_VALUE(peakIndex);

    //get descriptor to know nbrPeaks
    TwH5Desc desc;
    if (TwGetH5Descriptor(cFilename, &desc) != TwSuccess || index < -1 || index >= desc.nbrPeaks)
    {
      return(R_NilValue);
    }

    //allocating R resources needed to return list

    //label string
    SEXP label;
    if (index != -1)
    {
      PROTECT(label = allocVector(STRSXP, 1));
    }
    else
    {
      PROTECT(label = allocVector(STRSXP, desc.nbrPeaks - 1));
    }

    //mass
    SEXP mass;
    double *p_mass;
    if (index != -1)
    {
      PROTECT(mass = NEW_NUMERIC(1));
    }
    else
    {
      PROTECT(mass = NEW_NUMERIC(desc.nbrPeaks-1));
    }
    p_mass = NUMERIC_POINTER(mass);

    //loMass
    SEXP loMass;
    double *p_loMass;
    if (index != -1)
    {
      PROTECT(loMass = NEW_NUMERIC(1));
    }
    else
    {
      PROTECT(loMass = NEW_NUMERIC(desc.nbrPeaks-1));
    }
    p_loMass = NUMERIC_POINTER(loMass);

    //hiMass
    SEXP hiMass;
    double *p_hiMass;
    if (index != -1)
    {
      PROTECT(hiMass = NEW_NUMERIC(1));
    }
    else
    {
      PROTECT(hiMass = NEW_NUMERIC(desc.nbrPeaks-1));
    }
    p_hiMass = NUMERIC_POINTER(hiMass);

    //fill in actual values
    TPeakPar param;
    if (index != -1)
    {
      TwRetVal rv = TwGetPeakParametersFromH5(cFilename, &param, index);
      if (rv == TwSuccess)
      {
        SET_STRING_ELT(label, 0, mkChar(param.label));
        p_mass[0] = param.mass;
        p_loMass[0] = param.loMass;
        p_hiMass[0] = param.hiMass;
      }
      else
      {
        SET_STRING_ELT(label, 0, mkChar("error"));
        p_mass[0] = 0.0;
        p_loMass[0] = 0.0;
        p_hiMass[0] = 0.0;
      }
    }
    else
    {
      for (int j=0; j<desc.nbrPeaks-1; ++j)
      {
        TwGetPeakParametersFromH5(cFilename, &param, j);
        SET_STRING_ELT(label, j, mkChar(param.label));
        p_mass[j] = param.mass;
        p_loMass[j] = param.loMass;
        p_hiMass[j] = param.hiMass;
      }
    }

    //names for peakPar elements
    SEXP listNames;
    char const* names[4] = { "label", "mass", "loMass", "hiMass" };
    PROTECT(listNames = allocVector(STRSXP, 4));
    for (int j = 0; j < 4; ++j)
    {
      SET_STRING_ELT(listNames, j, mkChar(names[j]));
    }
    //list
    SEXP peakPar;
    PROTECT(peakPar = allocVector(VECSXP, 4));
    SET_VECTOR_ELT(peakPar, 0, label);
    SET_VECTOR_ELT(peakPar, 1, mass);
    SET_VECTOR_ELT(peakPar, 2, loMass);
    SET_VECTOR_ELT(peakPar, 3, hiMass);
    setAttrib(peakPar, R_NamesSymbol, listNames);

    TwCloseH5(cFilename);
    UNPROTECT(6);
    return peakPar;
  }

  DLL_EXPORT SEXP GetBufTimeFromH5(SEXP filename, SEXP bufIndex, SEXP writeIndex)
  {
    int cBufIndex = INTEGER_VALUE(bufIndex);
    int cWriteIndex = INTEGER_VALUE(writeIndex);

    //get the filename in a useable form
    char *cFilename = RtoCstring(filename);

    //prepare return array
    SEXP bufTime;
    double *cBufTime;
    PROTECT(bufTime = NEW_NUMERIC(1));
    cBufTime = NUMERIC_POINTER(bufTime);
    //read the actual values directly into cBufTime
    if (TwGetBufTimeFromH5(cFilename, cBufTime, cBufIndex, cWriteIndex) != TwSuccess)
    {
      TwCloseH5(cFilename);
      UNPROTECT(1);
      return(R_NilValue);
    }
    TwCloseH5(cFilename);
    UNPROTECT(1);
    return bufTime;
  }

  DLL_EXPORT SEXP GetSpecXaxisFromH5(SEXP filename, SEXP type, SEXP writeIndex)
  {
    //get the filename in a useable form
    char *cFilename = RtoCstring(filename);

    int index = INTEGER_VALUE(writeIndex);
    int xtype = INTEGER_VALUE(type);

    //get descriptor to know spec length
    TwH5Desc desc;
    if (TwGetH5Descriptor(cFilename, &desc) != TwSuccess)
    {
      TwCloseH5(cFilename);
      return(R_NilValue);
    }

    //prepare return array
    SEXP massAxis;
    double *massAxisBuf;
    int specLen = desc.nbrSamples;
    PROTECT(massAxis = NEW_NUMERIC(specLen));
    massAxisBuf = NUMERIC_POINTER(massAxis);

    TwRetVal rv = TwGetSpecXaxisFromH5(cFilename, massAxisBuf, xtype, NULL, 0.0, index);
    if (rv != TwSuccess)
    {
      TwCloseH5(cFilename);
      UNPROTECT(1);
      return(R_NilValue);
    }
    TwCloseH5(cFilename);
    UNPROTECT(1);
    return massAxis;
  }

  DLL_EXPORT SEXP GetSegmentProfileFromH5(SEXP filename, SEXP peakIndex, SEXP bufStart, SEXP bufEnd, SEXP writeStart, SEXP writeEnd, SEXP linked)
  {
    bool linkBufWrite = LOGICAL_VALUE(linked);
    int index = INTEGER_VALUE(peakIndex);
    int bufStartIndex = INTEGER_VALUE(bufStart);
    int bufEndIndex = INTEGER_VALUE(bufEnd);
    int writeStartIndex = INTEGER_VALUE(writeStart);
    int writeEndIndex = INTEGER_VALUE(writeEnd);

    //get the filename in a useable form
    char *cFilename = RtoCstring(filename);

    //get descriptor of file
    TwH5Desc desc;
    if (TwGetH5Descriptor(cFilename, &desc) != TwSuccess)
    {
      TwCloseH5(cFilename);
      return(R_NilValue);
    }
    //prepare return array
    SEXP segProfile;
    double *segProfileBuf;
    if (index != -1)
    {
      int specLen = desc.nbrSegments;
      PROTECT(segProfile = NEW_NUMERIC(specLen));
      segProfileBuf = NUMERIC_POINTER(segProfile);
      //allocate temp float buffer
      float *temp = new float[specLen];
      if (TwGetSegmentProfileFromH5(cFilename, temp, index, bufStartIndex, bufEndIndex,
                                    writeStartIndex, writeEndIndex, linkBufWrite) != TwSuccess)
      {
        TwCloseH5(cFilename);
        UNPROTECT(1);
        return(R_NilValue);
      }
      //copy to result
      for (int i = 0; i < specLen; ++i)
      {
        segProfileBuf[i] = temp[i];
      }
      delete[] temp;
      temp = nullptr;
    }
    else
    {
      int specLen = desc.nbrSegments*desc.nbrPeaks;
      PROTECT(segProfile = NEW_NUMERIC(specLen));
      segProfileBuf = NUMERIC_POINTER(segProfile);
      //allocate temp float buffer
      float *temp = new float[specLen];
      if (TwGetSegmentProfileFromH5(cFilename, temp, index, bufStartIndex, bufEndIndex,
                                    writeStartIndex, writeEndIndex, linkBufWrite) != TwSuccess)
      {
        TwCloseH5(cFilename);
        UNPROTECT(1);
        return(R_NilValue);
      }
      //copy to result
      for (int i = 0; i < specLen; ++i)
      {
        segProfileBuf[i] = temp[i];
      }
      delete[] temp;
      temp = nullptr;
    }

    TwCloseH5(cFilename);
    UNPROTECT(1);
    return segProfile;
  }

  // Not implemented: TwGetBufWriteProfileFromH5
  // Not implemented: TwGetBufWriteProfileFromH5_2

  DLL_EXPORT SEXP GetRegUserDataSourcesFromH5(SEXP filename)
  {
    //get the filename in a useable form
    char *cFilename = RtoCstring(filename);

    //get reg user data length
    int nbrSources = 0;
    if (TwGetRegUserDataSourcesFromH5(cFilename, &nbrSources, NULL, NULL, NULL, NULL) != TwValueAdjusted)
    {
      TwCloseH5(cFilename);
      return(R_NilValue);
    }

    //prepare return array
    SEXP sourceLocation, sourceLength, hasDesc, type, list, list_names;

    char *p_sourceLocation = new char[256 * nbrSources];
    memset(p_sourceLocation, 0, 256 * nbrSources);
    int *p_sourceLength;
    bool *p_hasDesc = new bool[nbrSources];  // needed because LOGICAL_POINTER(x) is int and not bool
    int *p_type;

    PROTECT(sourceLength = NEW_INTEGER(nbrSources));
    p_sourceLength = INTEGER_POINTER(sourceLength);
    UNPROTECT(1);

    PROTECT(type = NEW_INTEGER(nbrSources));
    p_type = INTEGER_POINTER(type);
    UNPROTECT(1);

    //TwRetVal rv = TwGetRegUserDataSourcesFromH5(cFilename, &nbrSources, p_sourceLocation, p_sourceLength, p_hasDesc, p_type);
    TwRetVal rv = TwGetRegUserDataSourcesFromH5(cFilename, &nbrSources, p_sourceLocation, p_sourceLength, p_hasDesc, p_type);
    if (rv != TwSuccess)
    {
      delete[] p_sourceLocation;
      p_sourceLocation = nullptr;
      delete[] p_hasDesc;
      p_hasDesc = nullptr;
      return(R_NilValue);
    }

    PROTECT(sourceLocation = allocVector(STRSXP, nbrSources));
    for (int i = 0; i < nbrSources; ++i)
    {
      SET_STRING_ELT(sourceLocation, i, mkChar(&p_sourceLocation[i * 256]));
    }
    UNPROTECT(1);

    int *p_hasDesc2;
    PROTECT(hasDesc = NEW_LOGICAL(nbrSources));
    p_hasDesc2 = LOGICAL_POINTER(hasDesc);
    for (int i = 0; i < nbrSources; ++i)
    {
      p_hasDesc2[i] = p_hasDesc[i];
    }
    UNPROTECT(1);

    delete[] p_sourceLocation;
    p_sourceLocation = nullptr;
    delete[] p_hasDesc;
    p_hasDesc = nullptr;

    PROTECT(list = allocVector(VECSXP, 4));
    SET_VECTOR_ELT(list, 0, sourceLocation);
    SET_VECTOR_ELT(list, 1, sourceLength);
    SET_VECTOR_ELT(list, 2, hasDesc);
    SET_VECTOR_ELT(list, 3, type);
    UNPROTECT(1);

    // list names
    char const* names[4] = { "sourceLocation", "sourceLength", "hasDesc", "type" };
    PROTECT(list_names = allocVector(STRSXP, 4));
    int i;
    for (i = 0; i < 4; i++)
      SET_STRING_ELT(list_names, i, mkChar(names[i]));
    setAttrib(list, R_NamesSymbol, list_names);
    UNPROTECT(1);

    return list;
  }

  DLL_EXPORT SEXP GetRegUserDataFromH5(SEXP filename, SEXP location, SEXP bufIndex, SEXP writeIndex, SEXP readDescription)
  {
    //get the filename in a useable form
    char *cFilename = RtoCstring(filename);

    //get the data location in a useable form
    char *cLocation = RtoCstring(location);

    int bIndex = INTEGER_VALUE(bufIndex);
    int wIndex = INTEGER_VALUE(writeIndex);
    bool readDesc = LOGICAL_VALUE(readDescription);

    //get reg user data length
    int regDataSize = 0;
    if (TwGetRegUserDataFromH5(cFilename, cLocation, bIndex, wIndex, &regDataSize, NULL, NULL) != TwValueAdjusted)
    {
      TwCloseH5(cFilename);
      return(R_NilValue);
    }

    //prepare return array
    SEXP regUserData, list, list_names;
    double *regUserDataBuf;

    PROTECT(regUserData = NEW_NUMERIC(regDataSize));
    regUserDataBuf = NUMERIC_POINTER(regUserData);

    if (readDesc)
    {
      char *description = new char[256 * regDataSize];
      memset(description, 0, 256 * regDataSize);

      TwRetVal rv = TwGetRegUserDataFromH5(cFilename, cLocation, bIndex, wIndex, &regDataSize, regUserDataBuf, description);
      TwCloseH5(cFilename);
      if (rv != TwSuccess)
      {
        UNPROTECT(1);
        delete[] description;
        description = nullptr;
        return(R_NilValue);
      }

      SEXP rNames;
      PROTECT(rNames = allocVector(STRSXP, regDataSize));
      for (int i = 0; i < regDataSize; ++i)
      {
        SET_STRING_ELT(rNames, i, mkChar(&description[i * 256]));
      }
      UNPROTECT(1);

      delete[] description;
      description = nullptr;

      PROTECT(list = allocVector(VECSXP, 2));
      SET_VECTOR_ELT(list, 0, regUserData);
      SET_VECTOR_ELT(list, 1, rNames);
      UNPROTECT(1);

      // list names
      char const* names[2] = { "data", "description" };
      PROTECT(list_names = allocVector(STRSXP, 2));
      int i;
      for (i = 0; i < 2; i++)
        SET_STRING_ELT(list_names, i, mkChar(names[i]));
      setAttrib(list, R_NamesSymbol, list_names);

      UNPROTECT(2);

      return list;
    }
    else
    {
      TwRetVal rv = TwGetRegUserDataFromH5(cFilename, cLocation, bIndex, wIndex, &regDataSize, regUserDataBuf, NULL);
      TwCloseH5(cFilename);
      if (rv != TwSuccess)
      {
        UNPROTECT(1);
        return(R_NilValue);
      }

      PROTECT(list = allocVector(VECSXP, 2));
      SET_VECTOR_ELT(list, 0, regUserData);
      SET_VECTOR_ELT(list, 1, R_NilValue);
      UNPROTECT(1);

      // list names
      char const* names[2] = { "data", "description" };
      PROTECT(list_names = allocVector(STRSXP, 2));
      int i;
      for (i = 0; i < 2; i++)
        SET_STRING_ELT(list_names, i, mkChar(names[i]));
      setAttrib(list, R_NamesSymbol, list_names);

      UNPROTECT(2);

      return list;
    }
  }

  DLL_EXPORT SEXP GetTofData(SEXP filename, SEXP sampleOffset, SEXP sampleCount, SEXP segOffset, SEXP segCount, SEXP bufOffset, SEXP bufCount, SEXP writeOffset, SEXP writeCount)
  {
    //get the filename in a useable form
    char *cFilename = RtoCstring(filename);

    int sampleOff = INTEGER_VALUE(sampleOffset);
    int sampleC = INTEGER_VALUE(sampleCount);
    int segOff= INTEGER_VALUE(segOffset);
    int segC = INTEGER_VALUE(segCount);
    int bufOff = INTEGER_VALUE(bufOffset);
    int bufC = INTEGER_VALUE(bufCount);
    int writeOff = INTEGER_VALUE(writeOffset);
    int writeC = INTEGER_VALUE(writeCount);

    int nx = sampleC*segC*bufC*writeC;

    //prepare return array
    SEXP data;
    double *dataBuffer;
    PROTECT(data = NEW_NUMERIC(nx));
    dataBuffer = NUMERIC_POINTER(data);
    //allocate temp float buffer
    float *temp = new float[nx];

    TwRetVal rv = TwGetTofData(cFilename, sampleOff, sampleC, segOff, segC, bufOff, bufC, writeOff, writeC, temp);
    if (rv != TwSuccess)
    {
      TwCloseH5(cFilename);
      UNPROTECT(1);
      return(R_NilValue);
    }
    //copy to result
    for (int i = 0; i < nx; ++i)
    {
      dataBuffer[i] = temp[i];
    }
    delete[] temp;
    temp = nullptr;

    TwCloseH5(cFilename);
    UNPROTECT(1);

    return data;
  }

  DLL_EXPORT SEXP GetTofData2(SEXP filename, SEXP sampleOffset, SEXP sampleCount, SEXP segOffset, SEXP segCount, SEXP bufOffset, SEXP bufCount, SEXP writeOffset, SEXP writeCount)
  {
    //get the filename in a useable form
    char *cFilename = RtoCstring(filename);

    int sampleOff = INTEGER_VALUE(sampleOffset);
    int sampleC = INTEGER_VALUE(sampleCount);
    int segOff = INTEGER_VALUE(segOffset);
    int segC = INTEGER_VALUE(segCount);
    int bufOff = INTEGER_VALUE(bufOffset);
    int bufC = INTEGER_VALUE(bufCount);
    int writeOff = INTEGER_VALUE(writeOffset);
    int writeC = INTEGER_VALUE(writeCount);

    int nx = sampleC*segC*bufC*writeC;

    //prepare return array
    SEXP data;
    double *dataBuffer;
    PROTECT(data = NEW_NUMERIC(nx));
    dataBuffer = NUMERIC_POINTER(data);
    //allocate temp float buffer
    float *temp = new float[nx];

    TwRetVal rv = TwGetTofData2(cFilename, sampleOff, sampleC, segOff, segC, bufOff, bufC, writeOff, writeC, temp);
    if (rv != TwSuccess)
    {
      TwCloseH5(cFilename);
      UNPROTECT(1);
      return(R_NilValue);
    }
    //copy to result
    for (int i = 0; i < nx; ++i)
    {
      dataBuffer[i] = temp[i];
    }
    delete[] temp;
    temp = nullptr;

    TwCloseH5(cFilename);
    UNPROTECT(1);

    return data;
  }

  DLL_EXPORT SEXP GetPeakData(SEXP filename, SEXP peakOffset, SEXP peakCount, SEXP segOffset, SEXP segCount, SEXP bufOffset, SEXP bufCount, SEXP writeOffset, SEXP writeCount)
  {
    //get the filename in a useable form
    char *cFilename = RtoCstring(filename);

    int peakOff = INTEGER_VALUE(peakOffset);
    int peakC = INTEGER_VALUE(peakCount);
    int segOff = INTEGER_VALUE(segOffset);
    int segC = INTEGER_VALUE(segCount);
    int bufOff = INTEGER_VALUE(bufOffset);
    int bufC = INTEGER_VALUE(bufCount);
    int writeOff = INTEGER_VALUE(writeOffset);
    int writeC = INTEGER_VALUE(writeCount);

    int nx = peakC*segC*bufC*writeC;

    //prepare return array
    SEXP data;
    double *dataBuffer;
    PROTECT(data = NEW_NUMERIC(nx));
    dataBuffer = NUMERIC_POINTER(data);
    //allocate temp float buffer
    float *temp = new float[nx];

    TwRetVal rv = TwGetPeakData(cFilename, peakOff, peakC, segOff, segC, bufOff, bufC, writeOff, writeC, temp);
    if (rv != TwSuccess)
    {
      TwCloseH5(cFilename);
      UNPROTECT(1);
      return(R_NilValue);
    }
    //copy to result
    for (int i = 0; i < nx; ++i)
    {
      dataBuffer[i] = temp[i];
    }
    delete[] temp;
    temp = nullptr;

    TwCloseH5(cFilename);
    UNPROTECT(1);

    return data;
  }

  DLL_EXPORT SEXP GetPeakData2(SEXP filename, SEXP peakOffset, SEXP peakCount, SEXP segOffset, SEXP segCount, SEXP bufOffset, SEXP bufCount, SEXP writeOffset, SEXP writeCount)
  {
    //get the filename in a useable form
    char *cFilename = RtoCstring(filename);

    int peakOff = INTEGER_VALUE(peakOffset);
    int peakC = INTEGER_VALUE(peakCount);
    int segOff = INTEGER_VALUE(segOffset);
    int segC = INTEGER_VALUE(segCount);
    int bufOff = INTEGER_VALUE(bufOffset);
    int bufC = INTEGER_VALUE(bufCount);
    int writeOff = INTEGER_VALUE(writeOffset);
    int writeC = INTEGER_VALUE(writeCount);

    int nx = peakC*segC*bufC*writeC;

    //prepare return array
    SEXP data;
    double *dataBuffer;
    PROTECT(data = NEW_NUMERIC(nx));
    dataBuffer = NUMERIC_POINTER(data);
    //allocate temp float buffer
    float *temp = new float[nx];

    TwRetVal rv = TwGetPeakData2(cFilename, peakOff, peakC, segOff, segC, bufOff, bufC, writeOff, writeC, temp);
    if (rv != TwSuccess)
    {
      TwCloseH5(cFilename);
      UNPROTECT(1);
      return(R_NilValue);
    }
    //copy to result
    for (int i = 0; i < nx; ++i)
    {
      dataBuffer[i] = temp[i];
    }
    delete[] temp;
    temp = nullptr;

    TwCloseH5(cFilename);
    UNPROTECT(1);

    return data;
  }

  DLL_EXPORT SEXP GetTimingData(SEXP filename, SEXP bufOffset, SEXP bufCount, SEXP writeOffset, SEXP writeCount)
  {
    //get the filename in a useable form
    char *cFilename = RtoCstring(filename);

    int bufOff = INTEGER_VALUE(bufOffset);
    int bufC = INTEGER_VALUE(bufCount);
    int writeOff = INTEGER_VALUE(writeOffset);
    int writeC = INTEGER_VALUE(writeCount);

    int nx = bufC*writeC;

    //prepare return array
    SEXP data;
    double *dataBuffer;
    PROTECT(data = NEW_NUMERIC(nx));
    dataBuffer = NUMERIC_POINTER(data);

    TwRetVal rv = TwGetTimingData(cFilename,  bufOff, bufC, writeOff, writeC, dataBuffer);
    if (rv != TwSuccess)
    {
      TwCloseH5(cFilename);
      UNPROTECT(1);
      return(R_NilValue);
    }

    TwCloseH5(cFilename);
    UNPROTECT(1);

    return data;
  }

  DLL_EXPORT SEXP GetIntAttributeFromH5(SEXP filename, SEXP location, SEXP name)
  {
    //get the filename in a useable form
    char *cFilename = RtoCstring(filename);

    //get the data location in a useable form
    char *cLocation = RtoCstring(location);

    //get the attribute in a useable form
    char *cName = RtoCstring(name);

    int value;
    TwRetVal rv = TwGetIntAttributeFromH5(cFilename, cLocation, cName, &value);
    if (rv != TwSuccess)
    {
      TwCloseH5(cFilename);
      return(R_NilValue);
    }

    SEXP rValue;
    PROTECT(rValue = NEW_INTEGER(1));
    int *pValue = INTEGER_POINTER(rValue);
    pValue[0] = value;
    UNPROTECT(1);

    TwCloseH5(cFilename);

    return rValue;
  }

  DLL_EXPORT SEXP GetUintAttributeFromH5(SEXP filename, SEXP location, SEXP name)
  {
    //get the filename in a useable form
    char *cFilename = RtoCstring(filename);

    //get the data location in a useable form
    char *cLocation = RtoCstring(location);

    //get the attribute in a useable form
    char *cName = RtoCstring(name);

    unsigned int value;
    TwRetVal rv = TwGetUintAttributeFromH5(cFilename, cLocation, cName, &value);
    if (rv != TwSuccess)
    {
      TwCloseH5(cFilename);
      return(R_NilValue);
    }
    // read unsigned int as numeric
    SEXP rValue;
    PROTECT(rValue = NEW_NUMERIC(1));
    double *pValue = NUMERIC_POINTER(rValue);
    pValue[0] = value;
    UNPROTECT(1);

    TwCloseH5(cFilename);

    return rValue;
  }

  DLL_EXPORT SEXP GetInt64AttributeFromH5(SEXP filename, SEXP location, SEXP name)
  {
    //get the filename in a useable form
    char *cFilename = RtoCstring(filename);

    //get the data location in a useable form
    char *cLocation = RtoCstring(location);

    //get the attribute in a useable form
    char *cName = RtoCstring(name);

    int64_t value;
    TwRetVal rv = TwGetInt64AttributeFromH5(cFilename, cLocation, cName, &value);
    if (rv != TwSuccess)
    {
      TwCloseH5(cFilename);
      return(R_NilValue);
    }
    // read int64 as string
    char *pint64char = new char[64];
    memset(pint64char, 0, 64);

    sprintf(pint64char, "%lld", value);
    SEXP rValue;
    PROTECT(rValue = NEW_CHARACTER(1));
    SET_STRING_ELT(rValue, 0, mkChar(pint64char));
    UNPROTECT(1);

    delete[] pint64char;
    pint64char = nullptr;

    TwCloseH5(cFilename);

    return rValue;
  }

  DLL_EXPORT SEXP GetUint64AttributeFromH5(SEXP filename, SEXP location, SEXP name)
  {
    //get the filename in a useable form
    char *cFilename = RtoCstring(filename);

    //get the data location in a useable form
    char *cLocation = RtoCstring(location);

    //get the attribute in a useable form
    char *cName = RtoCstring(name);

    uint64_t value;
    TwRetVal rv = TwGetUint64AttributeFromH5(cFilename, cLocation, cName, &value);
    if (rv != TwSuccess)
    {
      TwCloseH5(cFilename);
      return(R_NilValue);
    }
    // read int64 as string
    char *pint64char = new char[64];
    memset(pint64char, 0, 64);

    sprintf(pint64char, "%llud", value);
    SEXP rValue;
    PROTECT(rValue = NEW_CHARACTER(1));
    SET_STRING_ELT(rValue, 0, mkChar(pint64char));
    UNPROTECT(1);

    delete[] pint64char;
    pint64char = nullptr;

    TwCloseH5(cFilename);

    return rValue;
  }

  DLL_EXPORT SEXP GetFloatAttributeFromH5(SEXP filename, SEXP location, SEXP name)
  {
    //get the filename in a useable form
    char *cFilename = RtoCstring(filename);

    //get the data location in a useable form
    char *cLocation = RtoCstring(location);

    //get the attribute in a useable form
    char *cName = RtoCstring(name);

    float value;
    TwRetVal rv = TwGetFloatAttributeFromH5(cFilename, cLocation, cName, &value);
    if (rv != TwSuccess)
    {
      TwCloseH5(cFilename);
      return(R_NilValue);
    }
    SEXP rValue;
    PROTECT(rValue = NEW_NUMERIC(1));
    double *pValue = NUMERIC_POINTER(rValue);
    pValue[0] = value;
    UNPROTECT(1);

    TwCloseH5(cFilename);

    return rValue;
  }

  DLL_EXPORT SEXP GetDoubleAttributeFromH5(SEXP filename, SEXP location, SEXP name)
  {
    //get the filename in a useable form
    char *cFilename = RtoCstring(filename);

    //get the data location in a useable form
    char *cLocation = RtoCstring(location);

    //get the attribute in a useable form
    char *cName = RtoCstring(name);

    double value;
    TwRetVal rv = TwGetDoubleAttributeFromH5(cFilename, cLocation, cName, &value);
    if (rv != TwSuccess)
    {
      TwCloseH5(cFilename);
      return(R_NilValue);
    }
    SEXP rValue;
    PROTECT(rValue = NEW_NUMERIC(1));
    double *pValue = NUMERIC_POINTER(rValue);
    pValue[0] = value;
    UNPROTECT(1);

    TwCloseH5(cFilename);

    return rValue;
  }

  DLL_EXPORT SEXP GetStringAttributeFromH5(SEXP filename, SEXP location, SEXP name)
  {
    //get the filename in a useable form
    char *cFilename = RtoCstring(filename);

    //get the data location in a useable form
    char *cLocation = RtoCstring(location);

    //get the attribute in a useable form
    char *cName = RtoCstring(name);

    char *buffer = new char[256];
    memset(buffer, 0, 256);

    TwRetVal rv = TwGetStringAttributeFromH5(cFilename, cLocation, cName, buffer);
    if (rv != TwSuccess)
    {
      TwCloseH5(cFilename);
      delete[] buffer;
      buffer = nullptr;
      return(R_NilValue);
    }

    SEXP rValue;
    PROTECT(rValue = allocVector(STRSXP, 1));
    SET_STRING_ELT(rValue, 0, mkChar(buffer));
    UNPROTECT(1);

    TwCloseH5(cFilename);
    delete[] buffer;
    buffer = nullptr;

    return rValue;
  }

  DLL_EXPORT SEXP SetIntAttributeInH5(SEXP filename, SEXP location, SEXP name, SEXP attribute)
  {
    //get the filename in a useable form
    char *cFilename = RtoCstring(filename);

    //get the data location in a useable form
    char *cLocation = RtoCstring(location);

    //get the attribute in a useable form
    char *cName = RtoCstring(name);

    int value = INTEGER_VALUE(attribute);
    TwRetVal rv = TwSetIntAttributeInH5(cFilename, cLocation, cName, value);
    TwCloseH5(cFilename);

    SEXP rValue;
    PROTECT(rValue = NEW_INTEGER(1));
    int *pValue = INTEGER_POINTER(rValue);
    pValue[0] = rv;
    UNPROTECT(1);

    return rValue;
  }

  // Not implemented: TwSetUintAttributeInH5
  // Not implemented: TwSetInt64AttributeInH5
  // Not implemented: TwSetUInt64AttributeInH5
  // Not implemented: TwSetFloatAttributeInH5

  DLL_EXPORT SEXP SetDoubleAttributeInH5(SEXP filename, SEXP location, SEXP name, SEXP attribute)
  {
    //get the filename in a useable form
    char *cFilename = RtoCstring(filename);

    //get the data location in a useable form
    char *cLocation = RtoCstring(location);

    //get the attribute in a useable form
    char *cName = RtoCstring(name);

    double value = NUMERIC_VALUE(attribute);
    TwRetVal rv = TwSetDoubleAttributeInH5(cFilename, cLocation, cName, value);
    TwCloseH5(cFilename);

    SEXP rValue;
    PROTECT(rValue = NEW_INTEGER(1));
    int *pValue = INTEGER_POINTER(rValue);
    pValue[0] = rv;
    UNPROTECT(1);

    return rValue;
  }

  DLL_EXPORT SEXP SetStringAttributeInH5(SEXP filename, SEXP location, SEXP name, SEXP attribute)
  {
    //get the filename in a useable form
    char *cFilename = RtoCstring(filename);

    //get the data location in a useable form
    char *cLocation = RtoCstring(location);

    //get the name in a useable form
    char *cName = RtoCstring(name);

    //get the string in a useable form
    char *cvalue = RtoCstring(attribute);

    TwRetVal rv = TwSetStringAttributeInH5(cFilename, cLocation, cName, cvalue);
    TwCloseH5(cFilename);

    SEXP rValue;
    PROTECT(rValue = NEW_INTEGER(1));
    int *pValue = INTEGER_POINTER(rValue);
    pValue[0] = rv;
    UNPROTECT(1);

    return rValue;
  }

  // Not implemented: TwChangePeakTable
  // Not implemented: TwChangePeakTable2
  // Not implemented: TwChangePeakFromFile
  // Not implemented: TwChangePeakFromFile2
  // Not implemented: TwProgressCallback
  // Not implemented: TwProgressCallback2
  // Not implemented: TwReadRawData
  // Not implemented: TwH5SetMassCalib
  // Not implemented: TwH5SetMassCalibEx

  DLL_EXPORT SEXP GetUserDataFromH5(SEXP filename, SEXP location, SEXP rowIndex)
  {
    //get the filename in a useable form
    char *cFilename = RtoCstring(filename);

    //get the data location in a useable form
    char *cLocation = RtoCstring(location);

    int rIndex = INTEGER_VALUE(rowIndex);

    //get buffer data length
    int bufferSize = 0;
    int tmp = 0;
    if (TwGetUserDataFromH5(cFilename, cLocation, &tmp, &bufferSize, NULL, NULL) != TwValueAdjusted)
    {
      TwCloseH5(cFilename);
      return(R_NilValue);
    }

    //prepare return array
    SEXP userData, list, list_names;
    double *userDataBuf;

    char *description = new char[256 * bufferSize];
    memset(description, 0, 256 * bufferSize);

    PROTECT(userData = NEW_NUMERIC(bufferSize));
    userDataBuf = NUMERIC_POINTER(userData);

    TwRetVal rv = TwGetUserDataFromH5(cFilename, cLocation, &rIndex, &bufferSize, userDataBuf, description);
    TwCloseH5(cFilename);
    if (rv != TwSuccess)
    {
      UNPROTECT(1);
      delete[] description;
      return(R_NilValue);
    }

    SEXP rNames;
    PROTECT(rNames = allocVector(STRSXP, bufferSize));
    for (int i = 0; i < bufferSize; ++i)
    {
      SET_STRING_ELT(rNames, i, mkChar(&description[i * 256]));
    }
    UNPROTECT(1);

    delete[] description;
    description = nullptr;

    PROTECT(list = allocVector(VECSXP, 2));
    SET_VECTOR_ELT(list, 0, userData);
    SET_VECTOR_ELT(list, 1, rNames);
    UNPROTECT(1);

    // list names
    char const* names[2] = { "data", "description" };
    PROTECT(list_names = allocVector(STRSXP, 2));
    int i;
    for (i = 0; i < 2; i++)
      SET_STRING_ELT(list_names, i, mkChar(names[i]));
    setAttrib(list, R_NamesSymbol, list_names);

    UNPROTECT(2);

    return list;
  }

  DLL_EXPORT SEXP GetAcquisitionLogFromH5(SEXP filename, SEXP index)
  {
    //get the filename in a useable form
    char *cFilename = RtoCstring(filename);

    int cindex = INTEGER_VALUE(index);

    //prepare return array
    SEXP logtxt, time, list, list_names;
    int64_t timestamp;

    char *logText = new char[256];
    memset(logText, 0, 256);

    TwRetVal rv = TwGetAcquisitionLogFromH5(cFilename, cindex, &timestamp, logText);
    TwCloseH5(cFilename);
    if (rv != TwSuccess)
    {
      TwCloseH5(cFilename);
      delete[] logText;
      logText = nullptr;
      return(R_NilValue);
    }

    PROTECT(logtxt = allocVector(STRSXP, 1));
    SET_STRING_ELT(logtxt, 0, mkChar(logText));
    UNPROTECT(1);

    delete[] logText;
    logText = nullptr;

    // read int64 as string
    char *pint64char = new char[64];
    memset(pint64char, 0, 64);

    sprintf(pint64char, "%lld", timestamp);
    PROTECT(time = NEW_CHARACTER(1));
    SET_STRING_ELT(time, 0, mkChar(pint64char));
    UNPROTECT(1);

    delete[] pint64char;
    pint64char = nullptr;

    PROTECT(list = allocVector(VECSXP, 2));
    SET_VECTOR_ELT(list, 0, time);
    SET_VECTOR_ELT(list, 1, logtxt);
    UNPROTECT(1);

    // list names
    char const* names[2] = { "timestamp", "logText" };
    PROTECT(list_names = allocVector(STRSXP, 2));
    int i;
    for (i = 0; i < 2; i++)
      SET_STRING_ELT(list_names, i, mkChar(names[i]));
    setAttrib(list, R_NamesSymbol, list_names);

    UNPROTECT(1);

    return list;
  }

  DLL_EXPORT SEXP GetEventListSpectrumFromH5(SEXP filename, SEXP segmentIndex, SEXP bufIndex, SEXP writeIndex)
  {
    int cSegmentIndex = INTEGER_VALUE(segmentIndex);
    int cBufIndex = INTEGER_VALUE(bufIndex);
    int cWriteIndex = INTEGER_VALUE(writeIndex);

    //get the filename in a useable form
    char *cFilename = RtoCstring(filename);

    //get bufferSize
    int bufferSize = 0;
    if (TwGetEventListSpectrumFromH5(cFilename, cSegmentIndex, cBufIndex, cWriteIndex, &bufferSize, NULL) != TwValueAdjusted)
    {
      TwCloseH5(cFilename);
      return(R_NilValue);
    }

    //prepare return array
    SEXP eventList;
    double *dbuffer;
    PROTECT(eventList = NEW_NUMERIC(bufferSize));
    dbuffer = NUMERIC_POINTER(eventList);

    unsigned int *buffer = new unsigned int[bufferSize];;
    TwRetVal rv = TwGetEventListSpectrumFromH5(cFilename, cSegmentIndex, cBufIndex, cWriteIndex, &bufferSize, buffer);
    if (rv != TwSuccess)
    {
      TwCloseH5(cFilename);
      UNPROTECT(1);
      return(R_NilValue);
    }

    //copy to result
    for (int i = 0; i < bufferSize; ++i)
    {
      dbuffer[i] = buffer[i];
    }
    delete[] buffer;
    buffer = nullptr;

    TwCloseH5(cFilename);
    UNPROTECT(1);
    return eventList;
  }

  // Not implemented: TwGetEventListDataFromH5
  // Not implemented: TwGetEventListBlobFromH5
  // Not implemented: TwFreeEventListData
  // Not implemented: TwFreeEventListData2
  // Not implemented: TwH5GetMassCalibPar
  // Not implemented: TwMultiPeakFitIntegration
  // Not implemented: TwH5AddLogEntry
  // Not implemented: TwH5AddUserDataMultiRow
  // Not implemented: TwH5SetMassCalibDynamic
  // Not implemented: TwGenerateSegmentProfilesFromEventList

  // TwToolDll ----------------------------------------------------------------

  // Not implemented: TwTranslateReturnValue

  DLL_EXPORT SEXP FitSinglePeak(SEXP yVals, SEXP xVals, SEXP peakType, SEXP blOffset, SEXP blSlope, SEXP amplitude, SEXP fwhmLo, SEXP fwhmHi, SEXP peakPos, SEXP mu)
  {
    int nx = length(xVals);

    double *p_xVals = NUMERIC_POINTER(xVals);
    double *p_yVals = NUMERIC_POINTER(yVals);

    int ptype = INTEGER_VALUE(peakType);
    double bloff = NUMERIC_VALUE(blOffset);
    double blsl = NUMERIC_VALUE(blSlope);
    double amp = NUMERIC_VALUE(amplitude);
    double fwlo = NUMERIC_VALUE(fwhmLo);
    double fwhi = NUMERIC_VALUE(fwhmHi);
    double ppos = NUMERIC_VALUE(peakPos);
    double m = NUMERIC_VALUE(mu);

    TwRetVal rv = TwFitSinglePeak(nx, p_yVals, p_xVals, ptype, &bloff, &blsl, &amp, &fwlo, &fwhi, &ppos, &m);

    if (rv != TwSuccess)
    {
      return(R_NilValue);
    }

    // Creating a list as return value:
    SEXP list, list_names;
    PROTECT(list = allocVector(VECSXP, 7));

    // creating list elements
    SEXP rValue;
    double *pValue;

    PROTECT(rValue = NEW_NUMERIC(1));
    pValue = NUMERIC_POINTER(rValue);
    pValue[0] = bloff;
    SET_VECTOR_ELT(list, 0, rValue);
    UNPROTECT(1);

    PROTECT(rValue = NEW_NUMERIC(1));
    pValue = NUMERIC_POINTER(rValue);
    pValue[0] = blsl;
    SET_VECTOR_ELT(list, 1, rValue);
    UNPROTECT(1);

    PROTECT(rValue = NEW_NUMERIC(1));
    pValue = NUMERIC_POINTER(rValue);
    pValue[0] = amp;
    SET_VECTOR_ELT(list, 2, rValue);
    UNPROTECT(1);

    PROTECT(rValue = NEW_NUMERIC(1));
    pValue = NUMERIC_POINTER(rValue);
    pValue[0] = fwlo;
    SET_VECTOR_ELT(list, 3, rValue);
    UNPROTECT(1);

    PROTECT(rValue = NEW_NUMERIC(1));
    pValue = NUMERIC_POINTER(rValue);
    pValue[0] = fwhi;
    SET_VECTOR_ELT(list, 4, rValue);
    UNPROTECT(1);

    PROTECT(rValue = NEW_NUMERIC(1));
    pValue = NUMERIC_POINTER(rValue);
    pValue[0] = ppos;
    SET_VECTOR_ELT(list, 5, rValue);
    UNPROTECT(1);

    PROTECT(rValue = NEW_NUMERIC(1));
    pValue = NUMERIC_POINTER(rValue);
    pValue[0] = m;
    SET_VECTOR_ELT(list, 6, rValue);
    UNPROTECT(1);

    // list names
    char const* names[7] = { "blOffset", "blSlope", "amplitude", "fwhmLo", "fwhmHi", "peakPos", "mu" };
    PROTECT(list_names = allocVector(STRSXP, 7));
    int i;
    for (i = 0; i < 7; i++)
      SET_STRING_ELT(list_names, i, mkChar(names[i]));
    setAttrib(list, R_NamesSymbol, list_names);

    UNPROTECT(2);

    return list;
  }

  DLL_EXPORT SEXP FitSinglePeak2(SEXP yVals, SEXP xVals, SEXP peakType, SEXP param)
  {
    int nx = length(xVals);

    double *p_xVals = NUMERIC_POINTER(xVals);
    double *p_yVals = NUMERIC_POINTER(yVals);

    int ptype = INTEGER_VALUE(peakType);

    double *p_param = NUMERIC_POINTER(param);

    TwRetVal rv = TwFitSinglePeak2(nx, p_yVals, p_xVals, ptype, p_param);

    if (rv != TwSuccess)
    {
      return(R_NilValue);
    }

    // Creating a list as return value:
    SEXP list, list_names;
    PROTECT(list = allocVector(VECSXP, 7));

    // creating list elements
    SEXP rValue;
    double *pValue;

    for (int j = 0; j < 7; ++j)
    {
      PROTECT(rValue = NEW_NUMERIC(1));
      pValue = NUMERIC_POINTER(rValue);
      pValue[0] = p_param[j];
      SET_VECTOR_ELT(list, j, rValue);
      UNPROTECT(1);
    }

    // list names
    char const* names[7] = { "blOffset", "blSlope", "amplitude", "fwhmLo", "fwhmHi", "peakPos", "mu" };
    PROTECT(list_names = allocVector(STRSXP, 7));
    int i;
    for (i = 0; i < 7; i++)
      SET_STRING_ELT(list_names, i, mkChar(names[i]));
    setAttrib(list, R_NamesSymbol, list_names);

    UNPROTECT(2);

    return list;
  }

  DLL_EXPORT SEXP EvalSinglePeak(SEXP xVals, SEXP blOffset, SEXP blSlope, SEXP amplitude, SEXP fwhmLo, SEXP fwhmHi, SEXP peakPos, SEXP mu)
  {
    int nx = length(xVals);

    double *p_xVals = NUMERIC_POINTER(xVals);

    double bloff = NUMERIC_VALUE(blOffset);
    double blsl = NUMERIC_VALUE(blSlope);
    double amp = NUMERIC_VALUE(amplitude);
    double fwlo = NUMERIC_VALUE(fwhmLo);
    double fwhi = NUMERIC_VALUE(fwhmHi);
    double ppos = NUMERIC_VALUE(peakPos);
    double m = NUMERIC_VALUE(mu);

    double *param = new double[7];
    param[0] = bloff;
    param[1] = blsl;
    param[2] = amp;
    param[3] = fwlo;
    param[4] = fwhi;
    param[5] = ppos;
    param[6] = m;

    //prepare return array
    SEXP yValsFit;
    double *p_yValsFit;
    PROTECT(yValsFit = NEW_NUMERIC(nx));
    p_yValsFit = NUMERIC_POINTER(yValsFit);

    for (int j = 0; j<nx; ++j)
    {
      p_yValsFit[j] = TwEvalSinglePeak(p_xVals[j], param);
    }

    UNPROTECT(1);

    return yValsFit;
  }

  // Not implemented: TwMultiPeakFit
  // Not implemented: TwEvalMultiPeak
  // Not implemented: TwFitResolution
  // Not implemented: TwEvalResolution

  DLL_EXPORT SEXP GetMoleculeMass(SEXP molecule)
  {
    //get the molecule formula in a useable form
    char *cMolecule = RtoCstring(molecule);

    SEXP mass;
    double *p_mass;
    PROTECT(mass = NEW_NUMERIC(1));
    p_mass = NUMERIC_POINTER(mass);

    TwRetVal rv = TwGetMoleculeMass(cMolecule, p_mass);
    if (rv != TwSuccess)
    {
      UNPROTECT(1);
      return(R_NilValue);
    }

    UNPROTECT(1);
    return mass;
  }

  DLL_EXPORT SEXP GetIsotopePattern(SEXP molecule, SEXP abundanceLimit)
  {
    //get the molecule formula in a useable form
    char *cMolecule = RtoCstring(molecule);

    double cAbundanceLimit = NUMERIC_VALUE(abundanceLimit);
    int NbrIsotopes = 0;

    // get number of isotopes
    TwRetVal rv = TwGetIsotopePattern(cMolecule, cAbundanceLimit, &NbrIsotopes, NULL, NULL);
    if (rv != TwValueAdjusted)
    {
      return(R_NilValue);
    }

    // prepare return arrays
    SEXP isoMass, isoAbundance, list, list_names, rValue;
    int listLength = 3;
    char const* names[3] = { "isoMass", "isoAbundance", "TwRetVal" };
    PROTECT(list = allocVector(VECSXP, listLength));

    PROTECT(isoMass = NEW_NUMERIC(NbrIsotopes));
    double *isoMassBuf = NUMERIC_POINTER(isoMass);
    PROTECT(isoAbundance = NEW_NUMERIC(NbrIsotopes));
    double *isoAbundanceBuf = NUMERIC_POINTER(isoAbundance);
    PROTECT(rValue = NEW_INTEGER(1));
    int *pValue = INTEGER_POINTER(rValue);

    //allocate double buffers
    double *cisoMass = new double[NbrIsotopes];
    double *cisoAbundance = new double[NbrIsotopes];

    rv = TwGetIsotopePattern(cMolecule, cAbundanceLimit, &NbrIsotopes, cisoMass, cisoAbundance);
    if (rv != TwSuccess)
    {
      UNPROTECT(4);
      return(R_NilValue);
    }

    //copy to result
    for (int i = 0; i < NbrIsotopes; ++i)
    {
      isoMassBuf[i] = cisoMass[i];
      isoAbundanceBuf[i] = cisoAbundance[i];
    }
    delete[] cisoMass;
    cisoMass = nullptr;
    delete[] cisoAbundance;
    cisoAbundance = nullptr;

    SET_VECTOR_ELT(list, 0, isoMass);
    SET_VECTOR_ELT(list, 1, isoAbundance);
    pValue[0] = rv;
    SET_VECTOR_ELT(list, 2, rValue);

    // Create character string vector of the "names" attribute
    PROTECT(list_names = allocVector(STRSXP, listLength));
    int i;
    for (i = 0; i < listLength; i++)
      SET_STRING_ELT(list_names, i, mkChar(names[i]));
    // attach vector names
    setAttrib(list, R_NamesSymbol, list_names);
    UNPROTECT(5);

    return list;

  }

  DLL_EXPORT SEXP GetIsotopePattern2(SEXP molecule, SEXP abundanceLimit)
  {
    //get the molecule formula in a useable form
    char *cMolecule = RtoCstring(molecule);

    double cAbundanceLimit = NUMERIC_VALUE(abundanceLimit);
    int NbrIsotopes = 0;

    // get number of isotopes
    TwRetVal rv = TwGetIsotopePattern2(cMolecule, cAbundanceLimit, &NbrIsotopes, NULL, NULL);
    if (rv != TwValueAdjusted)
    {
      return(R_NilValue);
    }

    // prepare return arrays
    SEXP isoMass, isoAbundance, list, list_names, rValue;
    int listLength = 3;
    char const* names[3] = { "isoMass", "isoAbundance", "TwRetVal" };
    PROTECT(list = allocVector(VECSXP, listLength));

    PROTECT(isoMass = NEW_NUMERIC(NbrIsotopes));
    double *isoMassBuf = NUMERIC_POINTER(isoMass);
    PROTECT(isoAbundance = NEW_NUMERIC(NbrIsotopes));
    double *isoAbundanceBuf = NUMERIC_POINTER(isoAbundance);
    PROTECT(rValue = NEW_INTEGER(1));
    int *pValue = INTEGER_POINTER(rValue);

    //allocate double buffers
    double *cisoMass = new double[NbrIsotopes];
    double *cisoAbundance = new double[NbrIsotopes];

    rv = TwGetIsotopePattern2(cMolecule, cAbundanceLimit, &NbrIsotopes, cisoMass, cisoAbundance);
    if (rv != TwSuccess)
    {
      UNPROTECT(4);
      return(R_NilValue);
    }

    //copy to result
    for (int i = 0; i < NbrIsotopes; ++i)
    {
      isoMassBuf[i] = cisoMass[i];
      isoAbundanceBuf[i] = cisoAbundance[i];
    }
    delete[] cisoMass;
    cisoMass = nullptr;
    delete[] cisoAbundance;
    cisoAbundance = nullptr;

    SET_VECTOR_ELT(list, 0, isoMass);
    SET_VECTOR_ELT(list, 1, isoAbundance);
    pValue[0] = rv;
    SET_VECTOR_ELT(list, 2, rValue);

    // Create character string vector of the "names" attribute
    PROTECT(list_names = allocVector(STRSXP, listLength));
    int i;
    for (i = 0; i < listLength; i++)
      SET_STRING_ELT(list_names, i, mkChar(names[i]));
    // attach vector names
    setAttrib(list, R_NamesSymbol, list_names);
    UNPROTECT(5);

    return list;

  }

  // Not implemented: TwDecomposeMass
  // Not implemented: TwGetComposition
  // Not implemented: TwNistLibrarySearch
  // Not implemented: TwNistLibraryQueryResult

  DLL_EXPORT SEXP Tof2Mass(SEXP tofSample, SEXP massCalibMode, SEXP p)
  {
    //get input parameters
    double *tof = NUMERIC_POINTER(tofSample);

    int nbrSamples;
    nbrSamples = length(tofSample);

    int mode = INTEGER_VALUE(massCalibMode);
    double *par = NUMERIC_POINTER(p);

    //prepare return value
    SEXP massR;
    double *massC;
    PROTECT(massR = NEW_NUMERIC(nbrSamples));
    massC = NUMERIC_POINTER(massR);

    for (int i = 0; i < nbrSamples; ++i)
    {
      massC[i] = TwTof2Mass(tof[i], mode, par);
    }
    UNPROTECT(1);
    return massR;
  }

  DLL_EXPORT SEXP Mass2Tof(SEXP mass, SEXP massCalibMode, SEXP p)
  {
    //get input parameters
    double *massC = NUMERIC_POINTER(mass);

    int nbrSamples = length(mass);

    int mode = INTEGER_VALUE(massCalibMode);
    double *par = NUMERIC_POINTER(p);

    //prepare return value
    SEXP tofR;
    double *tofC;
    PROTECT(tofR = NEW_NUMERIC(nbrSamples));
    tofC = NUMERIC_POINTER(tofR);

    for (int i = 0; i < nbrSamples; ++i)
    {
      tofC[i] = TwMass2Tof(massC[i], mode, par);
    }
    UNPROTECT(1);
    return tofR;
  }

  DLL_EXPORT SEXP MassCalibrate(SEXP massCalibMode, SEXP mass, SEXP tof, SEXP weight)
  {
    //get input parameters
    int mode = INTEGER_VALUE(massCalibMode);

    int nx = length(mass);

    double *p_mass = NUMERIC_POINTER(mass);
    double *p_tof = NUMERIC_POINTER(tof);
    double *p_weight = NUMERIC_POINTER(weight);

    int nbrParams;

    // get the number of parameters
    char *description = new char[64];
    memset(description, 0, 64);
    if (TwGetMassCalibInfo(mode, description, &nbrParams) != TwSuccess)
    {
      delete[] description;
      description = nullptr;
      return(R_NilValue);
    }
    delete[] description;
    description = nullptr;

    SEXP param;
    double *p_param;
    PROTECT(param = NEW_NUMERIC(nbrParams));
    p_param = NUMERIC_POINTER(param);

    TwMassCalibrate(mode, nx, p_mass, p_tof, p_weight, &nbrParams, p_param, NULL, NULL);

    UNPROTECT(1);
    return(param);
  }

  // Not implemented: TwBruteForceCalibrate
  // Not implemented: TwGetMassCalibInfo
  // Not implemented: TwEncImsCorrelateProfile
  // Not implemented (removed in 1.98): TwEncImsSharpenProfile
  // Not implemented (removed in 1.98): TwEncImsDenoiseProfile
  // Not implemented: TwEncImsCorrelateMultiProfiles
  // Not implemented: TwEncImsCleanup

  DLL_EXPORT SEXP SiInitializeHistograms(SEXP loMass, SEXP hiMass, SEXP specType)
  {
    int nbrHist = length(loMass);
    double *p_loMass = NUMERIC_POINTER(loMass);
    double *p_hiMass = NUMERIC_POINTER(hiMass);
    int *p_specType = INTEGER_POINTER(specType);

    TwRetVal rv = TwSiInitializeHistograms(nbrHist, p_loMass, p_hiMass, p_specType);

    SEXP rValue;
    PROTECT(rValue = NEW_INTEGER(1));
    int *pValue = INTEGER_POINTER(rValue);
    pValue[0] = rv;
    UNPROTECT(1);

    return rValue;
  }

  DLL_EXPORT SEXP SiSetProcessingOptions(SEXP option, SEXP value, SEXP specType)
  {
    double cvalue = NUMERIC_VALUE(value);
    int spect = INTEGER_VALUE(specType);

    //get the option string in a useable form
    char *coption = RtoCstring(option);

    TwRetVal rv = TwSiSetProcessingOptions(coption, cvalue, spect);

    SEXP rValue;
    PROTECT(rValue = NEW_INTEGER(1));
    int *pValue = INTEGER_POINTER(rValue);
    pValue[0] = rv;
    UNPROTECT(1);

    return rValue;
  }

  DLL_EXPORT SEXP SiProcessSpectrum(SEXP spectrum, SEXP specType)
  {
    double *p_spectrum = NUMERIC_POINTER(spectrum);

    int nbrSamples = length(spectrum);
    int spect = INTEGER_VALUE(specType);

    //allocate float buffer
    float *p_spectrumf = new float[nbrSamples];
    //copy p_spectrum to p_spectrumf
    for (int i = 0; i < nbrSamples; ++i)
    {
      p_spectrumf[i] = p_spectrum[i];
    }

    //prepare return array
    SEXP blFromData, thrFromData, list, list_names;
    float cblFromData;
    float cthrFromData;

    if (TwSiProcessSpectrum(p_spectrumf, nbrSamples, spect, &cblFromData, &cthrFromData) != TwSuccess)
    {
      delete[] p_spectrumf;
      p_spectrumf = nullptr;
      return(R_NilValue);
    }

    delete[] p_spectrumf;
    p_spectrumf = nullptr;

    PROTECT(blFromData = NEW_NUMERIC(1));
    double *p_blFromData = NUMERIC_POINTER(blFromData);
    p_blFromData[0] = cblFromData;
    UNPROTECT(1);
    PROTECT(thrFromData = NEW_NUMERIC(1));
    double *p_thrFromData = NUMERIC_POINTER(thrFromData);
    p_thrFromData[0] = cthrFromData;
    UNPROTECT(1);

    // Creating a list with 2 vector elements:
    PROTECT(list = allocVector(VECSXP, 2));
    SET_VECTOR_ELT(list, 0, blFromData);
    SET_VECTOR_ELT(list, 1, thrFromData);

    // Creating a character string vector of the "names" attribute of the
    // objects in out list:
    char const* names[2] = { "baseline", "threshold" };
    PROTECT(list_names = allocVector(STRSXP, 2));
    int i;
    for (i = 0; i < 2; i++)
    {
      SET_STRING_ELT(list_names, i, mkChar(names[i]));
    }

    // and attaching the vector names:
    setAttrib(list, R_NamesSymbol, list_names);
    UNPROTECT(2);
    return list;

  }

  DLL_EXPORT SEXP SiGetHistogram(SEXP histogramIndex)
  {
    int histidx = INTEGER_VALUE(histogramIndex);

    //prepare return array
    SEXP intensity, counts, spectrumCount, meanValue, list, list_names;
    unsigned int arrayLength = -1;  // initialize to maximum value: 2^32-1
    unsigned int speccount;
    double meanval;

    //get length of array
    if (TwSiGetHistogram(histidx, NULL, NULL, &arrayLength, &speccount, &meanval) != TwValueAdjusted)
    {
      return(R_NilValue);
    }

    double *intens;
    double *cnts;  // get unsigned int as double

    PROTECT(intensity = NEW_NUMERIC(arrayLength));
    intens = NUMERIC_POINTER(intensity);
    PROTECT(counts = NEW_NUMERIC(arrayLength));
    cnts = NUMERIC_POINTER(counts);

    //allocate temp float and unsigned int buffers
    float *temp1 = new float[arrayLength];
    unsigned int *temp2 = new unsigned int[arrayLength];

    TwRetVal rv = TwSiGetHistogram(histidx, temp1, temp2, &arrayLength, &speccount, &meanval);
    //rv = TwSiGetHistogram(histidx, temp1, temp2, &arrayLength, &speccount, &meanval);
    if (rv != TwSuccess)
    {
      UNPROTECT(2);
      delete[] temp1;
      temp1 = nullptr;
      delete[] temp2;
      temp2 = nullptr;
      return(R_NilValue);
    }
    //copy to result
    for (unsigned int i = 0; i < arrayLength; ++i)
    {
      intens[i] = temp1[i];
      cnts[i] = temp2[i];
    }

    delete[] temp1;
    temp1 = nullptr;
    delete[] temp2;
    temp2 = nullptr;

    PROTECT(meanValue = NEW_NUMERIC(1));
    double *pValue = NUMERIC_POINTER(meanValue);
    pValue[0] = meanval;
    UNPROTECT(1);

    // read unsigned int as numeric
    PROTECT(spectrumCount = NEW_NUMERIC(1));
    pValue = NUMERIC_POINTER(spectrumCount);
    pValue[0] = speccount;
    UNPROTECT(1);

    PROTECT(list = allocVector(VECSXP, 4));
    SET_VECTOR_ELT(list, 0, intensity);
    SET_VECTOR_ELT(list, 1, counts);
    SET_VECTOR_ELT(list, 2, spectrumCount);
    SET_VECTOR_ELT(list, 3, meanValue);
    UNPROTECT(1);

    // list names
    char const* names[4] = { "intensity", "counts", "spectrumCount", "meanValue" };
    PROTECT(list_names = allocVector(STRSXP, 4));
    int i;
    for (i = 0; i < 4; i++)
      SET_STRING_ELT(list_names, i, mkChar(names[i]));
    setAttrib(list, R_NamesSymbol, list_names);

    UNPROTECT(3);

    return list;
  }

  DLL_EXPORT SEXP SiGetSumHistogram(SEXP specType, SEXP minMass, SEXP maxMass, SEXP minRate, SEXP maxRate)
  {
    int spect = INTEGER_VALUE(specType);
    double minm = NUMERIC_VALUE(minMass);
    double maxm = NUMERIC_VALUE(maxMass);
    double minr = NUMERIC_VALUE(minRate);
    double maxr = NUMERIC_VALUE(maxRate);

    //prepare return array
    SEXP intensity, counts, spectrumCount, meanValue, list, list_names;
    unsigned int arrayLength = -1;  // initialize to maximum value: 2^32-1
    unsigned int speccount;
    double meanval;

    //get length of array
    if (TwSiGetSumHistogram(spect, NULL, NULL, &arrayLength, &speccount, &meanval, minm, maxm, minr, maxr) != TwValueAdjusted)
    {
      return(R_NilValue);
    }

    double *intens;
    double *cnts;  // get unsigned int as double

    PROTECT(intensity = NEW_NUMERIC(arrayLength));
    intens = NUMERIC_POINTER(intensity);
    PROTECT(counts = NEW_NUMERIC(arrayLength));
    cnts = NUMERIC_POINTER(counts);

    //allocate temp float and unsigned int buffers
    float *temp1 = new float[arrayLength];
    unsigned int *temp2 = new unsigned int[arrayLength];

    TwRetVal rv = TwSiGetSumHistogram(spect, temp1, temp2, &arrayLength, &speccount, &meanval, minm, maxm, minr, maxr);
    if (rv != TwSuccess)
    {
      UNPROTECT(2);
      delete[] temp1;
      temp1 = nullptr;
      delete[] temp2;
      temp2 = nullptr;
      return(R_NilValue);
    }
    //copy to result
    for (unsigned int i = 0; i < arrayLength; ++i)
    {
      intens[i] = temp1[i];
      cnts[i] = temp2[i];
    }

    delete[] temp1;
    temp1 = nullptr;
    delete[] temp2;
    temp2 = nullptr;

    PROTECT(meanValue = NEW_NUMERIC(1));
    double *pValue = NUMERIC_POINTER(meanValue);
    pValue[0] = meanval;
    UNPROTECT(1);

    // read unsigned int as numeric
    PROTECT(spectrumCount = NEW_NUMERIC(1));
    pValue = NUMERIC_POINTER(spectrumCount);
    pValue[0] = speccount;
    UNPROTECT(1);

    PROTECT(list = allocVector(VECSXP, 4));
    SET_VECTOR_ELT(list, 0, intensity);
    SET_VECTOR_ELT(list, 1, counts);
    SET_VECTOR_ELT(list, 2, spectrumCount);
    SET_VECTOR_ELT(list, 3, meanValue);
    UNPROTECT(1);

    // list names
    char const* names[4] = { "intensity", "counts", "spectrumCount", "meanValue" };
    PROTECT(list_names = allocVector(STRSXP, 4));
    int i;
    for (i = 0; i < 4; i++)
      SET_STRING_ELT(list_names, i, mkChar(names[i]));
    setAttrib(list, R_NamesSymbol, list_names);

    UNPROTECT(3);

    return list;
  }

  DLL_EXPORT SEXP SiResetHistograms()
  {
    TwRetVal rv = TwSiResetHistograms();
    SEXP rValue;
    PROTECT(rValue = NEW_INTEGER(1));
    int *pValue = INTEGER_POINTER(rValue);
    pValue[0] = rv;
    UNPROTECT(1);

    return rValue;
  }

  DLL_EXPORT SEXP SiCleanup()
  {
    TwRetVal rv = TwSiCleanup();
    SEXP rValue;
    PROTECT(rValue = NEW_INTEGER(1));
    int *pValue = INTEGER_POINTER(rValue);
    pValue[0] = rv;
    UNPROTECT(1);

    return rValue;
  }

  DLL_EXPORT SEXP SiFitPhd(SEXP intensity, SEXP counts)
  {
    double *intens = NUMERIC_POINTER(intensity);
    double *cnts = NUMERIC_POINTER(counts);

    int nx = length(intensity);

    //prepare return array
    SEXP fwhm, a, par, list, list_names;

    PROTECT(fwhm = NEW_NUMERIC(1));
    double *p_fwhm = NUMERIC_POINTER(fwhm);
    PROTECT(a = NEW_NUMERIC(1));
    double *p_a = NUMERIC_POINTER(a);
    PROTECT(par = NEW_NUMERIC(4));
    double *p_par = NUMERIC_POINTER(par);

    TwRetVal rv = TwSiFitPhd(nx, intens, cnts, p_fwhm, p_a, p_par);

    if (rv != TwSuccess)
    {
      UNPROTECT(3);
      return(R_NilValue);
    }

    // Creating a list as return value:
    PROTECT(list = allocVector(VECSXP, 3));
    SET_VECTOR_ELT(list, 0, fwhm);
    SET_VECTOR_ELT(list, 1, a);
    SET_VECTOR_ELT(list, 2, par);
    UNPROTECT(4);

    // list names
    char const* names[3] = { "fwhm", "a", "par" };
    PROTECT(list_names = allocVector(STRSXP, 3));
    int i;
    for (i = 0; i < 3; i++)
      SET_STRING_ELT(list_names, i, mkChar(names[i]));
    setAttrib(list, R_NamesSymbol, list_names);
    UNPROTECT(1);

    return list;
  }

  DLL_EXPORT SEXP SiEvalPhd(SEXP par, SEXP intensity)
  {
    int nx = length(intensity);

    double *cpar = NUMERIC_POINTER(par);
    double *intens = NUMERIC_POINTER(intensity);

    //prepare return array
    SEXP yValsFit;
    double *p_yValsFit;
    PROTECT(yValsFit = NEW_NUMERIC(nx));
    p_yValsFit = NUMERIC_POINTER(yValsFit);

    for (int j = 0; j<nx; ++j)
    {
      p_yValsFit[j] = TwSiEvalPhd(cpar, intens[j]);
    }

    UNPROTECT(1);

    return yValsFit;
  }

  DLL_EXPORT SEXP SiFitRateFromPhd(SEXP intensity, SEXP counts, SEXP siPar)
  {
    double *intens = NUMERIC_POINTER(intensity);
    double *cnts = NUMERIC_POINTER(counts);
    double *cpar = NUMERIC_POINTER(siPar);

    int nx = length(intensity);

    //prepare return array

    SEXP rate, fitCounts;
    double *p_rate;
    double *p_fitCounts;
    PROTECT(rate = NEW_NUMERIC(1));
    p_rate = NUMERIC_POINTER(rate);
    PROTECT(fitCounts = NEW_NUMERIC(nx));
    p_fitCounts = NUMERIC_POINTER(fitCounts);

    TwRetVal rv = TwSiFitRateFromPhd(nx, intens, cnts, cpar, p_rate, p_fitCounts, 0, NULL);

    if (rv != TwSuccess)
    {
      UNPROTECT(2);
      return(R_NilValue);
    }

    // Creating a list as return value:
    SEXP list, list_names;
    PROTECT(list = allocVector(VECSXP, 2));
    SET_VECTOR_ELT(list, 0, rate);
    SET_VECTOR_ELT(list, 1, fitCounts);
    UNPROTECT(1);

    // list names
    char const* names[2] = { "rate", "fitCounts" };
    PROTECT(list_names = allocVector(STRSXP, 2));
    int i;
    for (i = 0; i < 2; i++)
      SET_STRING_ELT(list_names, i, mkChar(names[i]));
    setAttrib(list, R_NamesSymbol, list_names);

    UNPROTECT(3);

    return list;
  }

#ifdef _WIN32
  DLL_EXPORT SEXP FindTpsIp(SEXP TpsSerial, SEXP timeout)
  {
    //get the TpsSerial in a useable form
    char *cTpsSerial = RtoCstring(TpsSerial);

    int tout = INTEGER_VALUE(timeout);
    int hostStrLen = 15;

    char *buffer = new char[15];
    memset(buffer, 0, 15);

    TwRetVal rv = TwFindTpsIp(cTpsSerial, tout, &hostStrLen, buffer);

    int listLength = 2;
    char const* names[2] = { "TwRetVal", "IP" };

    SEXP rValue, list, list_names;

    PROTECT(list = allocVector(VECSXP, listLength));

    // TwRetVal return value
    PROTECT(rValue = NEW_INTEGER(1));
    int *pValue = INTEGER_POINTER(rValue);
    pValue[0] = rv;
    SET_VECTOR_ELT(list, 0, rValue);
    UNPROTECT(1);

    PROTECT(rValue = allocVector(STRSXP, 1));
    if (rv != TwSuccess)
      SET_VECTOR_ELT(list, 1, R_NilValue);
    else
    {
      SET_STRING_ELT(rValue, 0, mkChar(buffer));
      SET_VECTOR_ELT(list, 1, rValue);
    }
    UNPROTECT(1);

    delete[] buffer;
    buffer = nullptr;

    // Create character string vector of the "names" attribute
    PROTECT(list_names = allocVector(STRSXP, listLength));
    int i;
    for (i = 0; i < listLength; i++)
      SET_STRING_ELT(list_names, i, mkChar(names[i]));
    // attach vector names
    setAttrib(list, R_NamesSymbol, list_names);
    UNPROTECT(2);

    return list;
  }

  // TofDaqDll ----------------------------------------------------------------

  // Control functions --------------------------------------------------------

  DLL_EXPORT SEXP InitializeDll()
  {
    TwRetVal rv = TwInitializeDll();
    SEXP rValue;
    PROTECT(rValue = NEW_INTEGER(1));
    int *pValue = INTEGER_POINTER(rValue);
    pValue[0] = rv;
    UNPROTECT(1);

    return rValue;
  }

  DLL_EXPORT SEXP CleanupDll()
  {
    void TwCleanupDll(void);
    return(R_NilValue);
  }

  DLL_EXPORT SEXP GetDllVersion()
  {
    double rv = TwGetDllVersion();
    SEXP rValue;
    PROTECT(rValue = NEW_NUMERIC(1));
    double *pValue = NUMERIC_POINTER(rValue);
    pValue[0] = rv;
    UNPROTECT(1);

    return rValue;
  }

  DLL_EXPORT SEXP TofDaqRunning()
  {
    bool rv = TwTofDaqRunning();
    SEXP rValue;
    PROTECT(rValue = NEW_LOGICAL(1));
    int *pValue = LOGICAL_POINTER(rValue);
    pValue[0] = rv;
    UNPROTECT(1);

    /*SEXP rValue;
    PROTECT(rValue = NEW_LOGICAL(1));
    int *pValue = LOGICAL_POINTER(rValue);
    pValue[0] = TwTofDaqRunning();
    UNPROTECT(1);*/

    /*int rv = TwTofDaqRunning();
    SEXP rValue;
    PROTECT(rValue = NEW_INTEGER(1));
    int *pValue = INTEGER_POINTER(rValue);
    pValue[0] = rv;
    UNPROTECT(1);*/

    return rValue;
  }

  DLL_EXPORT SEXP DaqActive()
  {
    bool rv = TwDaqActive();
    SEXP rValue;
    PROTECT(rValue = NEW_LOGICAL(1));
    int *pValue = LOGICAL_POINTER(rValue);
    pValue[0] = rv;
    UNPROTECT(1);

    return rValue;
  }

  DLL_EXPORT SEXP StartAcquisition()
  {
    TwRetVal rv = TwStartAcquisition();
    SEXP rValue;
    PROTECT(rValue = NEW_INTEGER(1));
    int *pValue = INTEGER_POINTER(rValue);
    pValue[0] = rv;
    UNPROTECT(1);

    return rValue;
  }

  DLL_EXPORT SEXP StopAcquisition()
  {
    TwRetVal rv = TwStopAcquisition();
    SEXP rValue;
    PROTECT(rValue = NEW_INTEGER(1));
    int *pValue = INTEGER_POINTER(rValue);
    pValue[0] = rv;
    UNPROTECT(1);

    return rValue;
  }

  // Not implemented: TwContinueAcquisition
  // Not implemented: TwManualContinueNeeded

  DLL_EXPORT SEXP CloseTofDaqRec()
  {
    TwRetVal rv = TwCloseTofDaqRec();
    SEXP rValue;
    PROTECT(rValue = NEW_INTEGER(1));
    int *pValue = INTEGER_POINTER(rValue);
    pValue[0] = rv;
    UNPROTECT(1);

    return rValue;
  }

  // Not implemented: TwLockBuf
  // Not implemented: TwUnLockBuf
  // Not implemented: TwIssueDio4Pulse
  // Not implemented: TwSetDio4State

  DLL_EXPORT SEXP InitializeDaqDevice()
  {
    TwRetVal rv = TwInitializeDaqDevice();
    SEXP rValue;
    PROTECT(rValue = NEW_INTEGER(1));
    int *pValue = INTEGER_POINTER(rValue);
    pValue[0] = rv;
    UNPROTECT(1);

    return rValue;
  }

  DLL_EXPORT SEXP SetTimeout(SEXP timeout)
  {
    //get input parameters
    int tout = INTEGER_VALUE(timeout);

    TwSetTimeout(tout);
    return(R_NilValue);
  }

  DLL_EXPORT SEXP GetTimeout()
  {
    int rv = TwGetTimeout();
    SEXP rValue;
    PROTECT(rValue = NEW_INTEGER(1));
    int *pValue = INTEGER_POINTER(rValue);
    pValue[0] = rv;
    UNPROTECT(1);

    return rValue;
  }

  DLL_EXPORT SEXP AutoSetupDaqDevice()
  {
    TwRetVal rv = TwAutoSetupDaqDevice();
    SEXP rValue;
    PROTECT(rValue = NEW_INTEGER(1));
    int *pValue = INTEGER_POINTER(rValue);
    pValue[0] = rv;
    UNPROTECT(1);

    return rValue;
  }

  DLL_EXPORT SEXP OnDemandMassCalibration(SEXP action)
  {
    //get input parameters
    int act = INTEGER_VALUE(action);

    TwRetVal rv = TwOnDemandMassCalibration(act);
    SEXP rValue;
    PROTECT(rValue = NEW_INTEGER(1));
    int *pValue = INTEGER_POINTER(rValue);
    pValue[0] = rv;
    UNPROTECT(1);

    return rValue;
  }

  // Configuration functions --------------------------------------------------
  DLL_EXPORT SEXP ShowConfigWindow(SEXP ConfigWindowIndex)
  {
    //get input parameters
    int index = INTEGER_VALUE(ConfigWindowIndex);

    TwRetVal rv = TwShowConfigWindow(index);
    SEXP rValue;
    PROTECT(rValue = NEW_INTEGER(1));
    int *pValue = INTEGER_POINTER(rValue);
    pValue[0] = rv;
    UNPROTECT(1);

    return rValue;
  }

  DLL_EXPORT SEXP LoadIniFile(SEXP IniFile)
  {
    if (strlen(CHAR(STRING_ELT(IniFile, 0))) == 0)
    {
      TwRetVal rv = TwLoadIniFile(NULL);
      SEXP rValue;
      PROTECT(rValue = NEW_INTEGER(1));
      int *pValue = INTEGER_POINTER(rValue);
      pValue[0] = rv;
      UNPROTECT(1);
      return rValue;
    }
    else
    {
      //get the filename in a useable form
      char *cFilename = RtoCstring(IniFile);

      TwRetVal rv = TwLoadIniFile(cFilename);
      SEXP rValue;
      PROTECT(rValue = NEW_INTEGER(1));
      int *pValue = INTEGER_POINTER(rValue);
      pValue[0] = rv;
      UNPROTECT(1);
      return rValue;
    }
  }

  DLL_EXPORT SEXP SaveIniFile(SEXP IniFile)
  {
    if (strlen(CHAR(STRING_ELT(IniFile, 0))) == 0)
    {
      TwRetVal rv = TwSaveIniFile(NULL);
      SEXP rValue;
      PROTECT(rValue = NEW_INTEGER(1));
      int *pValue = INTEGER_POINTER(rValue);
      pValue[0] = rv;
      UNPROTECT(1);
      return rValue;
    }
    else
    {
      //get the filename in a useable form
      char *cFilename = RtoCstring(IniFile);

      TwRetVal rv = TwSaveIniFile(cFilename);
      SEXP rValue;
      PROTECT(rValue = NEW_INTEGER(1));
      int *pValue = INTEGER_POINTER(rValue);
      pValue[0] = rv;
      UNPROTECT(1);
      return rValue;
    }
  }

  DLL_EXPORT SEXP GetDaqParameter(SEXP Parameter)
  {
    //get the Parameter in a useable form
    char *cParameter = RtoCstring(Parameter);

    char *buffer = new char[256];
    memset(buffer, 0, 256);

    buffer = TwGetDaqParameter(cParameter);

    SEXP rValue;
    PROTECT(rValue = allocVector(STRSXP, 1));
    SET_STRING_ELT(rValue, 0, mkChar(buffer));
    UNPROTECT(1);

    delete[] buffer;
    buffer = nullptr;

    return rValue;
  }

  DLL_EXPORT SEXP SetDaqParameter(SEXP Parameter, SEXP ValueString)
  {
    //get the Parameter in a useable form
    char *cParameter = RtoCstring(Parameter);

    //get the ValueString in a useable form
    char *cValueString = RtoCstring(ValueString);

    TwRetVal rv = TwSetDaqParameter(cParameter, cValueString);
    SEXP rValue;
    PROTECT(rValue = NEW_INTEGER(1));
    int *pValue = INTEGER_POINTER(rValue);
    pValue[0] = rv;
    UNPROTECT(1);

    return rValue;
  }

  // Not implemented: TwConfigVarNbrMemories
  // Not implemented: TwSetMassCalib
  // Not implemented: TwSetMassCalibEx

  // Data access functions ----------------------------------------------------
  DLL_EXPORT SEXP GetDescriptor()
  {
    //get the current TSharedMemoryDesc structure.
    TSharedMemoryDesc desc;
    if (TwGetDescriptor(&desc) != TwSuccess)
    {
      return(R_NilValue);
    }

    int descLength = 68;
    char const* names[68] = { "NbrSamples",
                        "NbrRawSamples",
                        "NbrPeaks",
                        "NbrWaveforms",
                        "NbrSegments",
                        "NbrBlocks",
                        "NbrMemories",
                        "NbrBufs",
                        "NbrWrites",
                        "NbrRuns",
                        "iWaveform",
                        "iSegment",
                        "iBlock",
                        "iMemory",
                        "iBuf",
                        "iWrite",
                        "iRun",
                        "TotalBufsRecorded",
                        "TotalBufsProcessed",
                        "TotalBufsWritten",
                        "OverallBufsProcessed",
                        "TotalNbrMemories",
                        "TotalMemoriesProcessed",
                        "RawDataRecordedBuf1",
                        "RawDataRecordedBuf2",
                        "RawDataLastElementInBuffer1",
                        "RawDataLastElementInBuffer2",
                        "RawDataProcessedBuf1",
                        "RawDataProcessedBuf2",
                        "RawDataWrittenBuf1",
                        "RawDataWrittenBuf2",
                        "SampleInterval",
                        "TofPeriod",
                        "BlockPeriod",
                        "BlockPulseDelay",
                        "BlockDelay",
                        "SingleIonSignal",
                        "SingleIonSignal2",
                        "MassCalibMode",
                        "MassCalibMode2",
                        "NbrMassCalibParams",
                        "NbrMassCalibParams2",
                        "p",
                        "p2",
                        "R0",
                        "dm",
                        "m0",
                        "SecondTof",
                        "chIniFileName",
                        "CurrentDataFileName",
                        "DaqMode",
                        "AcquisitionMode",
                        "CombineMode",
                        "RecalibFreq",
                        "AcquisitionLogText",
                        "AcquisitionLogTime",
                        "TimeZero",
                        "ExternalLock",
                        "ProcessingLevel",
                        "AttributeType",
                        "AttributeObject",
                        "AttributeName",
                        "AttributeInt",
                        "AttributeDouble",
                        "EnableVarNbrMemories",
                        "NbrSteps",
                        "CurrentStepAtBuf",
                        "NbrMemoriesForCurrentStep"
    };

    int *p_rdesc;
    double *pn_rdesc;
    char *pchar256_rdesc = new char[256];
    memset(pchar256_rdesc, 0, 256);
    char *pchar128_rdesc = new char[128];
    memset(pchar128_rdesc, 0, 128);
    char *pint64char = new char[64];
    memset(pint64char, 0, 64);
    SEXP rdesc, list, list_names;

    // Creating a list with 68 vector elements:
    PROTECT(list = allocVector(VECSXP, descLength));

    // creating descriptor elements
    PROTECT(rdesc = NEW_INTEGER(1));
    p_rdesc = INTEGER_POINTER(rdesc);
    p_rdesc[0] = desc.NbrSamples;
    SET_VECTOR_ELT(list, 0, rdesc);
    UNPROTECT(1);

    PROTECT(rdesc = NEW_INTEGER(1));
    p_rdesc = INTEGER_POINTER(rdesc);
    p_rdesc[0] = desc.NbrRawSamples;
    SET_VECTOR_ELT(list, 1, rdesc);
    UNPROTECT(1);

    PROTECT(rdesc = NEW_INTEGER(1));
    p_rdesc = INTEGER_POINTER(rdesc);
    p_rdesc[0] = desc.NbrPeaks;
    SET_VECTOR_ELT(list, 2, rdesc);
    UNPROTECT(1);

    PROTECT(rdesc = NEW_INTEGER(1));
    p_rdesc = INTEGER_POINTER(rdesc);
    p_rdesc[0] = desc.NbrWaveforms;
    SET_VECTOR_ELT(list, 3, rdesc);
    UNPROTECT(1);

    PROTECT(rdesc = NEW_INTEGER(1));
    p_rdesc = INTEGER_POINTER(rdesc);
    p_rdesc[0] = desc.NbrSegments;
    SET_VECTOR_ELT(list, 4, rdesc);
    UNPROTECT(1);

    PROTECT(rdesc = NEW_INTEGER(1));
    p_rdesc = INTEGER_POINTER(rdesc);
    p_rdesc[0] = desc.NbrBlocks;
    SET_VECTOR_ELT(list, 5, rdesc);
    UNPROTECT(1);

    PROTECT(rdesc = NEW_INTEGER(1));
    p_rdesc = INTEGER_POINTER(rdesc);
    p_rdesc[0] = desc.NbrMemories;
    SET_VECTOR_ELT(list, 6, rdesc);
    UNPROTECT(1);

    PROTECT(rdesc = NEW_INTEGER(1));
    p_rdesc = INTEGER_POINTER(rdesc);
    p_rdesc[0] = desc.NbrBufs;
    SET_VECTOR_ELT(list, 7, rdesc);
    UNPROTECT(1);

    PROTECT(rdesc = NEW_INTEGER(1));
    p_rdesc = INTEGER_POINTER(rdesc);
    p_rdesc[0] = desc.NbrWrites;
    SET_VECTOR_ELT(list, 8, rdesc);
    UNPROTECT(1);

    PROTECT(rdesc = NEW_INTEGER(1));
    p_rdesc = INTEGER_POINTER(rdesc);
    p_rdesc[0] = desc.NbrRuns;
    SET_VECTOR_ELT(list, 9, rdesc);
    UNPROTECT(1);

    PROTECT(rdesc = NEW_INTEGER(1));
    p_rdesc = INTEGER_POINTER(rdesc);
    p_rdesc[0] = desc.iWaveform;
    SET_VECTOR_ELT(list, 10, rdesc);
    UNPROTECT(1);

    PROTECT(rdesc = NEW_INTEGER(1));
    p_rdesc = INTEGER_POINTER(rdesc);
    p_rdesc[0] = desc.iSegment;
    SET_VECTOR_ELT(list, 11, rdesc);
    UNPROTECT(1);

    PROTECT(rdesc = NEW_INTEGER(1));
    p_rdesc = INTEGER_POINTER(rdesc);
    p_rdesc[0] = desc.iBlock;
    SET_VECTOR_ELT(list, 12, rdesc);
    UNPROTECT(1);

    PROTECT(rdesc = NEW_INTEGER(1));
    p_rdesc = INTEGER_POINTER(rdesc);
    p_rdesc[0] = desc.iMemory;
    SET_VECTOR_ELT(list, 13, rdesc);
    UNPROTECT(1);

    PROTECT(rdesc = NEW_INTEGER(1));
    p_rdesc = INTEGER_POINTER(rdesc);
    p_rdesc[0] = desc.iBuf;
    SET_VECTOR_ELT(list, 14, rdesc);
    UNPROTECT(1);

    PROTECT(rdesc = NEW_INTEGER(1));
    p_rdesc = INTEGER_POINTER(rdesc);
    p_rdesc[0] = desc.iWrite;
    SET_VECTOR_ELT(list, 15, rdesc);
    UNPROTECT(1);

    PROTECT(rdesc = NEW_INTEGER(1));
    p_rdesc = INTEGER_POINTER(rdesc);
    p_rdesc[0] = desc.iRun;
    SET_VECTOR_ELT(list, 16, rdesc);
    UNPROTECT(1);

    PROTECT(rdesc = NEW_INTEGER(1));
    p_rdesc = INTEGER_POINTER(rdesc);
    p_rdesc[0] = desc.TotalBufsRecorded;
    SET_VECTOR_ELT(list, 17, rdesc);
    UNPROTECT(1);

    PROTECT(rdesc = NEW_INTEGER(1));
    p_rdesc = INTEGER_POINTER(rdesc);
    p_rdesc[0] = desc.TotalBufsProcessed;
    SET_VECTOR_ELT(list, 18, rdesc);
    UNPROTECT(1);

    PROTECT(rdesc = NEW_INTEGER(1));
    p_rdesc = INTEGER_POINTER(rdesc);
    p_rdesc[0] = desc.TotalBufsWritten;
    SET_VECTOR_ELT(list, 19, rdesc);
    UNPROTECT(1);

    PROTECT(rdesc = NEW_INTEGER(1));
    p_rdesc = INTEGER_POINTER(rdesc);
    p_rdesc[0] = desc.OverallBufsProcessed;
    SET_VECTOR_ELT(list, 20, rdesc);
    UNPROTECT(1);

    PROTECT(rdesc = NEW_INTEGER(1));
    p_rdesc = INTEGER_POINTER(rdesc);
    p_rdesc[0] = desc.TotalNbrMemories;
    SET_VECTOR_ELT(list, 21, rdesc);
    UNPROTECT(1);

    PROTECT(rdesc = NEW_INTEGER(1));
    p_rdesc = INTEGER_POINTER(rdesc);
    p_rdesc[0] = desc.TotalMemoriesProcessed;
    SET_VECTOR_ELT(list, 22, rdesc);
    UNPROTECT(1);

    // read unsigned int as numeric
    PROTECT(rdesc = NEW_NUMERIC(1));
    pn_rdesc = NUMERIC_POINTER(rdesc);
    pn_rdesc[0] = desc.RawDataRecordedBuf1;
    SET_VECTOR_ELT(list, 23, rdesc);
    UNPROTECT(1);

    // read unsigned int as numeric
    PROTECT(rdesc = NEW_NUMERIC(1));
    pn_rdesc = NUMERIC_POINTER(rdesc);
    pn_rdesc[0] = desc.RawDataRecordedBuf2;
    SET_VECTOR_ELT(list, 24, rdesc);
    UNPROTECT(1);

    // read unsigned int as numeric
    PROTECT(rdesc = NEW_NUMERIC(1));
    pn_rdesc = NUMERIC_POINTER(rdesc);
    pn_rdesc[0] = desc.RawDataLastElementInBuffer1;
    SET_VECTOR_ELT(list, 25, rdesc);
    UNPROTECT(1);

    // read unsigned int as numeric
    PROTECT(rdesc = NEW_NUMERIC(1));
    pn_rdesc = NUMERIC_POINTER(rdesc);
    pn_rdesc[0] = desc.RawDataRecordedBuf2;
    SET_VECTOR_ELT(list, 26, rdesc);
    UNPROTECT(1);

    // read unsigned int as numeric
    PROTECT(rdesc = NEW_NUMERIC(1));
    pn_rdesc = NUMERIC_POINTER(rdesc);
    pn_rdesc[0] = desc.RawDataProcessedBuf1;
    SET_VECTOR_ELT(list, 27, rdesc);
    UNPROTECT(1);

    // read unsigned int as numeric
    PROTECT(rdesc = NEW_NUMERIC(1));
    pn_rdesc = NUMERIC_POINTER(rdesc);
    pn_rdesc[0] = desc.RawDataProcessedBuf2;
    SET_VECTOR_ELT(list, 28, rdesc);
    UNPROTECT(1);

    // read unsigned int as numeric
    PROTECT(rdesc = NEW_NUMERIC(1));
    pn_rdesc = NUMERIC_POINTER(rdesc);
    pn_rdesc[0] = desc.RawDataWrittenBuf1;
    SET_VECTOR_ELT(list, 29, rdesc);
    UNPROTECT(1);

    // read unsigned int as numeric
    PROTECT(rdesc = NEW_NUMERIC(1));
    pn_rdesc = NUMERIC_POINTER(rdesc);
    pn_rdesc[0] = desc.RawDataWrittenBuf2;
    SET_VECTOR_ELT(list, 30, rdesc);
    UNPROTECT(1);

    PROTECT(rdesc = NEW_NUMERIC(1));
    pn_rdesc = NUMERIC_POINTER(rdesc);
    pn_rdesc[0] = desc.SampleInterval;
    SET_VECTOR_ELT(list, 31, rdesc);
    UNPROTECT(1);

    PROTECT(rdesc = NEW_INTEGER(1));
    p_rdesc = INTEGER_POINTER(rdesc);
    p_rdesc[0] = desc.TofPeriod;
    SET_VECTOR_ELT(list, 32, rdesc);
    UNPROTECT(1);

    //// read int64 as numeric
    //PROTECT(rdesc = NEW_NUMERIC(1));
    //pn_rdesc = NUMERIC_POINTER(rdesc);
    //pn_rdesc[0] = desc.BlockPeriod;
    //SET_VECTOR_ELT(list, 33, rdesc);
    //UNPROTECT(1);

    // read unsigned int64 as string
    PROTECT(rdesc = NEW_CHARACTER(1));
    sprintf(pint64char, "%lld", desc.BlockPeriod);
    SET_STRING_ELT(rdesc, 0, mkChar(pint64char));
    SET_VECTOR_ELT(list, 33, rdesc);
    UNPROTECT(1);

    // read unsigned int64 as string
    PROTECT(rdesc = NEW_CHARACTER(1));
    sprintf(pint64char, "%lld", desc.BlockPulseDelay);
    SET_STRING_ELT(rdesc, 0, mkChar(pint64char));
    SET_VECTOR_ELT(list, 34, rdesc);
    UNPROTECT(1);

    // read unsigned int64 as string
    PROTECT(rdesc = NEW_CHARACTER(1));
    sprintf(pint64char, "%lld", desc.BlockDelay);
    SET_STRING_ELT(rdesc, 0, mkChar(pint64char));
    SET_VECTOR_ELT(list, 35, rdesc);
    UNPROTECT(1);

    PROTECT(rdesc = NEW_NUMERIC(1));
    pn_rdesc = NUMERIC_POINTER(rdesc);
    pn_rdesc[0] = desc.SingleIonSignal;
    SET_VECTOR_ELT(list, 36, rdesc);
    UNPROTECT(1);

    PROTECT(rdesc = NEW_NUMERIC(1));
    pn_rdesc = NUMERIC_POINTER(rdesc);
    pn_rdesc[0] = desc.SingleIonSignal2;
    SET_VECTOR_ELT(list, 37, rdesc);
    UNPROTECT(1);

    PROTECT(rdesc = NEW_INTEGER(1));
    p_rdesc = INTEGER_POINTER(rdesc);
    p_rdesc[0] = desc.MassCalibMode;
    SET_VECTOR_ELT(list, 38, rdesc);
    UNPROTECT(1);

    PROTECT(rdesc = NEW_INTEGER(1));
    p_rdesc = INTEGER_POINTER(rdesc);
    p_rdesc[0] = desc.MassCalibMode2;
    SET_VECTOR_ELT(list, 39, rdesc);
    UNPROTECT(1);

    PROTECT(rdesc = NEW_INTEGER(1));
    p_rdesc = INTEGER_POINTER(rdesc);
    p_rdesc[0] = desc.NbrMassCalibParams;
    SET_VECTOR_ELT(list, 40, rdesc);
    UNPROTECT(1);

    PROTECT(rdesc = NEW_INTEGER(1));
    p_rdesc = INTEGER_POINTER(rdesc);
    p_rdesc[0] = desc.NbrMassCalibParams2;
    SET_VECTOR_ELT(list, 41, rdesc);
    UNPROTECT(1);

    PROTECT(rdesc = NEW_NUMERIC(16));
    pn_rdesc = NUMERIC_POINTER(rdesc);
    for (int i = 0; i < 16; ++i)
    {
      if (desc.p[i] < 1e-200)
      {
        pn_rdesc[i] = 0;
      }
      else
      {
        pn_rdesc[i] = desc.p[i];
      }
    }
    SET_VECTOR_ELT(list, 42, rdesc);
    UNPROTECT(1);

    PROTECT(rdesc = NEW_NUMERIC(16));
    pn_rdesc = NUMERIC_POINTER(rdesc);
    for (int i = 0; i < 16; ++i)
    {
      if (desc.p2[i] < 1e-200)
      {
        pn_rdesc[i] = 0;
      }
      else
      {
        pn_rdesc[i] = desc.p2[i];
      }
    }
    SET_VECTOR_ELT(list, 43, rdesc);
    UNPROTECT(1);

    PROTECT(rdesc = NEW_NUMERIC(1));
    pn_rdesc = NUMERIC_POINTER(rdesc);
    pn_rdesc[0] = desc.R0;
    SET_VECTOR_ELT(list, 44, rdesc);
    UNPROTECT(1);

    PROTECT(rdesc = NEW_NUMERIC(1));
    pn_rdesc = NUMERIC_POINTER(rdesc);
    pn_rdesc[0] = desc.dm;
    SET_VECTOR_ELT(list, 45, rdesc);
    UNPROTECT(1);

    PROTECT(rdesc = NEW_NUMERIC(1));
    pn_rdesc = NUMERIC_POINTER(rdesc);
    pn_rdesc[0] = desc.m0;
    SET_VECTOR_ELT(list, 46, rdesc);
    UNPROTECT(1);

    PROTECT(rdesc = NEW_LOGICAL(1));
    p_rdesc = LOGICAL_POINTER(rdesc);
    p_rdesc[0] = desc.SecondTof;
    SET_VECTOR_ELT(list, 47, rdesc);
    UNPROTECT(1);

    PROTECT(rdesc = NEW_CHARACTER(1));
    pchar256_rdesc = desc.chIniFileName;
    SET_STRING_ELT(rdesc, 0, mkChar(pchar256_rdesc));
    SET_VECTOR_ELT(list, 48, rdesc);
    UNPROTECT(1);

    PROTECT(rdesc = NEW_CHARACTER(1));
    pchar256_rdesc = desc.CurrentDataFileName;
    SET_STRING_ELT(rdesc, 0, mkChar(pchar256_rdesc));
    SET_VECTOR_ELT(list, 49, rdesc);
    UNPROTECT(1);

    PROTECT(rdesc = NEW_INTEGER(1));
    p_rdesc = INTEGER_POINTER(rdesc);
    p_rdesc[0] = desc.DaqMode;
    SET_VECTOR_ELT(list, 50, rdesc);
    UNPROTECT(1);

    PROTECT(rdesc = NEW_INTEGER(1));
    p_rdesc = INTEGER_POINTER(rdesc);
    p_rdesc[0] = desc.AcquisitionMode;
    SET_VECTOR_ELT(list, 51, rdesc);
    UNPROTECT(1);

    PROTECT(rdesc = NEW_INTEGER(1));
    p_rdesc = INTEGER_POINTER(rdesc);
    p_rdesc[0] = desc.CombineMode;
    SET_VECTOR_ELT(list, 52, rdesc);
    UNPROTECT(1);

    PROTECT(rdesc = NEW_INTEGER(1));
    p_rdesc = INTEGER_POINTER(rdesc);
    p_rdesc[0] = desc.RecalibFreq;
    SET_VECTOR_ELT(list, 53, rdesc);
    UNPROTECT(1);

    PROTECT(rdesc = NEW_CHARACTER(1));
    pchar256_rdesc = desc.AcquisitionLogText;
    SET_STRING_ELT(rdesc, 0, mkChar(pchar256_rdesc));
    SET_VECTOR_ELT(list, 54, rdesc);
    UNPROTECT(1);

    // read unsigned int64 as string
    PROTECT(rdesc = NEW_CHARACTER(1));
    sprintf(pint64char, "%llu", desc.AcquisitionLogTime);
    SET_STRING_ELT(rdesc, 0, mkChar(pint64char));
    SET_VECTOR_ELT(list, 55, rdesc);
    UNPROTECT(1);

    // read unsigned int64 as string
    PROTECT(rdesc = NEW_CHARACTER(1));
    sprintf(pint64char, "%llu", desc.TimeZero);
    SET_STRING_ELT(rdesc, 0, mkChar(pint64char));
    SET_VECTOR_ELT(list, 56, rdesc);
    UNPROTECT(1);

    PROTECT(rdesc = NEW_LOGICAL(1));
    p_rdesc = LOGICAL_POINTER(rdesc);
    p_rdesc[0] = desc.ExternalLock;
    SET_VECTOR_ELT(list, 57, rdesc);
    UNPROTECT(1);

    // read unsigned int as numeric
    PROTECT(rdesc = NEW_NUMERIC(1));
    pn_rdesc = NUMERIC_POINTER(rdesc);
    pn_rdesc[0] = desc.ProcessingLevel;
    SET_VECTOR_ELT(list, 58, rdesc);
    UNPROTECT(1);

    PROTECT(rdesc = NEW_INTEGER(1));
    p_rdesc = INTEGER_POINTER(rdesc);
    p_rdesc[0] = desc.AttributeType;
    SET_VECTOR_ELT(list, 59, rdesc);
    UNPROTECT(1);

    PROTECT(rdesc = NEW_CHARACTER(1));
    pchar256_rdesc = desc.AttributeObject;
    SET_STRING_ELT(rdesc, 0, mkChar(pchar256_rdesc));
    SET_VECTOR_ELT(list, 60, rdesc);
    UNPROTECT(1);

    PROTECT(rdesc = NEW_CHARACTER(1));
    pchar128_rdesc = desc.AttributeName;
    SET_STRING_ELT(rdesc, 0, mkChar(pchar128_rdesc));
    SET_VECTOR_ELT(list, 61, rdesc);
    UNPROTECT(1);

    PROTECT(rdesc = NEW_INTEGER(1));
    p_rdesc = INTEGER_POINTER(rdesc);
    p_rdesc[0] = desc.AttributeInt;
    SET_VECTOR_ELT(list, 62, rdesc);
    UNPROTECT(1);

    PROTECT(rdesc = NEW_NUMERIC(1));
    pn_rdesc = NUMERIC_POINTER(rdesc);
    pn_rdesc[0] = desc.AttributeDouble;
    SET_VECTOR_ELT(list, 63, rdesc);
    UNPROTECT(1);

    PROTECT(rdesc = NEW_LOGICAL(1));
    p_rdesc = LOGICAL_POINTER(rdesc);
    p_rdesc[0] = desc.EnableVarNbrMemories;
    SET_VECTOR_ELT(list, 64, rdesc);
    UNPROTECT(1);

    PROTECT(rdesc = NEW_INTEGER(1));
    p_rdesc = INTEGER_POINTER(rdesc);
    p_rdesc[0] = desc.NbrSteps;
    SET_VECTOR_ELT(list, 65, rdesc);
    UNPROTECT(1);

    PROTECT(rdesc = NEW_INTEGER(1));
    p_rdesc = INTEGER_POINTER(rdesc);
    p_rdesc[0] = desc.CurrentStepAtBuf;
    SET_VECTOR_ELT(list, 66, rdesc);
    UNPROTECT(1);

    PROTECT(rdesc = NEW_INTEGER(1));
    p_rdesc = INTEGER_POINTER(rdesc);
    p_rdesc[0] = desc.NbrMemoriesForCurrentStep;
    SET_VECTOR_ELT(list, 67, rdesc);
    UNPROTECT(1);

    delete[] pchar256_rdesc;
    pchar256_rdesc = nullptr;
    delete[] pchar128_rdesc;
    pchar128_rdesc = nullptr;
    delete[] pint64char;
    pint64char = nullptr;

    // Creating a character string vector
    // of the "names" attribute of the
    // objects in out list:
    PROTECT(list_names = allocVector(STRSXP, descLength));
    int i;
    for (i = 0; i < descLength; i++)
      SET_STRING_ELT(list_names, i, mkChar(names[i]));

    // and attaching the vector names:
    setAttrib(list, R_NamesSymbol, list_names);
    UNPROTECT(2);
    return list;
  }

  DLL_EXPORT SEXP GetPeakParameters(SEXP PeakIndex)
  {
    int peakidx = INTEGER_VALUE(PeakIndex);

    int listLength = 5;
    char const* names[5] = { "TwRetVal", "label", "mass", "loMass", "hiMass" };

    SEXP rValue, list, list_names;

    TPeakPar PeakPar;
    TwRetVal rv = TwGetPeakParameters(&PeakPar, peakidx);

    PROTECT(list = allocVector(VECSXP, listLength));

    // TwRetVal return value
    PROTECT(rValue = NEW_INTEGER(1));
    int *pValue = INTEGER_POINTER(rValue);
    pValue[0] = rv;
    SET_VECTOR_ELT(list, 0, rValue);
    UNPROTECT(1);

    // label
    char *pchar64_rdesc = new char[64];
    memset(pchar64_rdesc, 0, 64);

    PROTECT(rValue = NEW_CHARACTER(1));
    if (rv != TwSuccess)
      SET_VECTOR_ELT(list, 1, R_NilValue);
    else
    {
      pchar64_rdesc = PeakPar.label;
      SET_STRING_ELT(rValue, 0, mkChar(pchar64_rdesc));
      SET_VECTOR_ELT(list, 1, rValue);
    }
    UNPROTECT(1);

    delete[] pchar64_rdesc;
    pchar64_rdesc = nullptr;

    // mass
    PROTECT(rValue = NEW_NUMERIC(1));
    if (rv != TwSuccess)
      SET_VECTOR_ELT(list, 1, R_NilValue);
    else
    {
      double *pValue = NUMERIC_POINTER(rValue);
      pValue[0] = PeakPar.mass;
      SET_VECTOR_ELT(list, 2, rValue);
    }
    UNPROTECT(1);

    // loMass
    PROTECT(rValue = NEW_NUMERIC(1));
    if (rv != TwSuccess)
      SET_VECTOR_ELT(list, 1, R_NilValue);
    else
    {
      double *pValue = NUMERIC_POINTER(rValue);
      pValue[0] = PeakPar.loMass;
      SET_VECTOR_ELT(list, 3, rValue);
    }
    UNPROTECT(1);

    // hiMass
    PROTECT(rValue = NEW_NUMERIC(1));
    if (rv != TwSuccess)
      SET_VECTOR_ELT(list, 1, R_NilValue);
    else
    {
      double *pValue = NUMERIC_POINTER(rValue);
      pValue[0] = PeakPar.hiMass;
      SET_VECTOR_ELT(list, 4, rValue);
    }
    UNPROTECT(1);

    // Create character string vector of the "names" attribute
    PROTECT(list_names = allocVector(STRSXP, listLength));
    int i;
    for (i = 0; i < listLength; i++)
      SET_STRING_ELT(list_names, i, mkChar(names[i]));
    // attach vector names
    setAttrib(list, R_NamesSymbol, list_names);
    UNPROTECT(2);

    return list;
  }

  // Not implemented: TwGetSharedMemory

  DLL_EXPORT SEXP ReleaseSharedMemory()
  {
    TwRetVal rv = TwReleaseSharedMemory();
    SEXP rValue;
    PROTECT(rValue = NEW_INTEGER(1));
    int *pValue = INTEGER_POINTER(rValue);
    pValue[0] = rv;
    UNPROTECT(1);

    return rValue;
  }

  DLL_EXPORT SEXP WaitForNewData(SEXP timeout, SEXP WaitForEventReset)
  {
    int tout = INTEGER_VALUE(timeout);
    TSharedMemoryDesc tsmd;
    TSharedMemoryPointer tsmp;
    bool wfer = LOGICAL_VALUE(WaitForEventReset);

    TwRetVal rv = TwWaitForNewData(tout, &tsmd, &tsmp, wfer);
    SEXP rValue;
    PROTECT(rValue = NEW_INTEGER(1));
    int *pValue = INTEGER_POINTER(rValue);
    pValue[0] = rv;
    UNPROTECT(1);

    return rValue;
  }

  DLL_EXPORT SEXP WaitForEndOfAcquisition(SEXP timeout)
  {
    int tout = INTEGER_VALUE(timeout);

    TwRetVal rv = TwWaitForEndOfAcquisition(tout);
    SEXP rValue;
    PROTECT(rValue = NEW_INTEGER(1));
    int *pValue = INTEGER_POINTER(rValue);
    pValue[0] = rv;
    UNPROTECT(1);

    return rValue;
  }

  // Not implemented: TwGetMassCalib
  // Not implemented: TwGetMassCalibEx

  DLL_EXPORT SEXP GetSumSpectrumFromShMem(SEXP Normalize)
  {
    bool norm = LOGICAL_VALUE(Normalize);

    //get descriptor of file
    TSharedMemoryDesc desc;
    if (TwGetDescriptor(&desc) != TwSuccess)
    {
      return(R_NilValue);
    }

    //prepare return array
    SEXP sumSpec;
    double *sumSpecBuf;
    int sumSpecLen = desc.NbrSamples;
    PROTECT(sumSpec = NEW_NUMERIC(sumSpecLen));
    sumSpecBuf = NUMERIC_POINTER(sumSpec);
    //read the actual values directly into sumSpecBuf
    if (TwGetSumSpectrumFromShMem(sumSpecBuf, norm) != TwSuccess)
    {
      UNPROTECT(1);
      return(R_NilValue);
    }
    UNPROTECT(1);
    return sumSpec;
  }

  DLL_EXPORT SEXP GetTofSpectrumFromShMem(SEXP SegmentIndex, SEXP SegmentEndIndex, SEXP BufIndex, SEXP Normalize)
  {
    bool norm = LOGICAL_VALUE(Normalize);
    int segidx = INTEGER_VALUE(SegmentIndex);
    int segendidx = INTEGER_VALUE(SegmentEndIndex);
    int bufidx = INTEGER_VALUE(BufIndex);

    //get descriptor of file
    TSharedMemoryDesc desc;
    if (TwGetDescriptor(&desc) != TwSuccess)
    {
      return(R_NilValue);
    }

    //prepare return array
    SEXP tofSpec;
    double *tofSpecBuf;
    int specLen;
    if (segidx == -1 && segendidx == -1)
    {
      specLen = desc.NbrSamples*desc.NbrSegments;
    }
    else
    {
      specLen = desc.NbrSamples;
    }
    PROTECT(tofSpec = NEW_NUMERIC(specLen));
    tofSpecBuf = NUMERIC_POINTER(tofSpec);
    //allocate temp float buffer
    float *temp = new float[specLen];
    TwRetVal rv = TwGetTofSpectrumFromShMem(temp, segidx, segendidx, bufidx, norm);
    //copy to result
    for (int i = 0; i < specLen; ++i)
    {
      tofSpecBuf[i] = temp[i];
    }
    delete[] temp;
    temp = nullptr;
    UNPROTECT(1);
    if (rv != TwSuccess) { return(R_NilValue); }
    return tofSpec;
  }

  DLL_EXPORT SEXP GetSpecXaxisFromShMem(SEXP Type)
  {
    int typ = INTEGER_VALUE(Type);
    double maxMass = 0.0;

    //get descriptor of file
    TSharedMemoryDesc desc;
    if (TwGetDescriptor(&desc) != TwSuccess)
    {
      return(R_NilValue);
    }

    //prepare return array
    SEXP SpecAxis;
    double *SpecAxisBuf;
    int SpecAxisLen = desc.NbrSamples;
    PROTECT(SpecAxis = NEW_NUMERIC(SpecAxisLen));
    SpecAxisBuf = NUMERIC_POINTER(SpecAxis);
    //read the actual values directly into SpecAxisBuf
    if (TwGetSpecXaxisFromShMem(SpecAxisBuf, typ, NULL, maxMass) != TwSuccess)
    {
      UNPROTECT(1);
      return(R_NilValue);
    }
    UNPROTECT(1);
    return SpecAxis;
  }

  DLL_EXPORT SEXP GetStickSpectrumFromShMem(SEXP SegmentIndex, SEXP SegmentEndIndex, SEXP BufIndex)
  {
    int segidx = INTEGER_VALUE(SegmentIndex);
    int segendidx = INTEGER_VALUE(SegmentEndIndex);
    int bufidx = INTEGER_VALUE(BufIndex);

    //get descriptor of file
    TSharedMemoryDesc desc;
    if (TwGetDescriptor(&desc) != TwSuccess)
    {
      return(R_NilValue);
    }

    //prepare return array
    SEXP stickSpec, Masses, list, list_names;
    double *stickSpecBuf;
    double *MassesBuf;
    int specLen = desc.NbrPeaks;
    PROTECT(stickSpec = NEW_NUMERIC(specLen));
    stickSpecBuf = NUMERIC_POINTER(stickSpec);
    PROTECT(Masses = NEW_NUMERIC(specLen));
    MassesBuf = NUMERIC_POINTER(Masses);
    //allocate temp float buffer
    float *temp1 = new float[specLen];
    float *temp2 = new float[specLen];
    if (TwGetStickSpectrumFromShMem(temp1, temp2, segidx, segendidx, bufidx) != TwSuccess)
    {
      UNPROTECT(2);
      return(R_NilValue);
    }

    //copy to result
    for (int i = 0; i < specLen; ++i)
    {
      stickSpecBuf[i] = temp1[i];
      MassesBuf[i] = temp2[i];
    }
    delete[] temp1;
    temp1 = nullptr;
    delete[] temp2;
    temp2 = nullptr;

    // Creating a list with 2 vector elements:
    PROTECT(list = allocVector(VECSXP, 2));
    SET_VECTOR_ELT(list, 0, stickSpec);
    SET_VECTOR_ELT(list, 1, Masses);

    // Creating a character string vector of the "names" attribute of the
    // objects in out list:
    char const* names[2] = { "Spectrum", "Masses" };
    PROTECT(list_names = allocVector(STRSXP, 2));
    int i;
    for (i = 0; i < 2; i++)
    {
      SET_STRING_ELT(list_names, i, mkChar(names[i]));
    }

    // and attaching the vector names:
    setAttrib(list, R_NamesSymbol, list_names);
    UNPROTECT(4);
    return list;
  }

  DLL_EXPORT SEXP GetSegmentProfileFromShMem(SEXP PeakIndex, SEXP BufIndex)
  {
    int peakidx = INTEGER_VALUE(PeakIndex);
    int bufidx = INTEGER_VALUE(BufIndex);

    //get descriptor of file
    TSharedMemoryDesc desc;
    if (TwGetDescriptor(&desc) != TwSuccess)
    {
      return(R_NilValue);
    }

    //prepare return array
    SEXP SegmentProfile;
    double *SegmentProfileBuf;
    int profileLen;
    if (peakidx == -1)
    {
      profileLen = desc.NbrSegments*desc.NbrPeaks;
    }
    else
    {
      profileLen = desc.NbrSegments;
    }
    PROTECT(SegmentProfile = NEW_NUMERIC(profileLen));
    SegmentProfileBuf = NUMERIC_POINTER(SegmentProfile);
    //allocate temp float buffer
    float *temp = new float[profileLen];
    TwRetVal rv = TwGetSegmentProfileFromShMem(temp, peakidx, bufidx);
    //copy to result
    for (int i = 0; i < profileLen; ++i)
    {
      SegmentProfileBuf[i] = temp[i];
    }
    delete[] temp;
    temp = nullptr;
    UNPROTECT(1);
    if (rv != TwSuccess) { return(R_NilValue); }
    return SegmentProfile;
  }

  DLL_EXPORT SEXP GetBufTimeFromShMem(SEXP BufIndex, SEXP WriteIndex)
  {
    int bufidx = INTEGER_VALUE(BufIndex);
    int writeidx = INTEGER_VALUE(WriteIndex);

    //prepare return value
    SEXP BufTime;
    double *BufTimeBuf;
    PROTECT(BufTime = NEW_NUMERIC(1));
    BufTimeBuf = NUMERIC_POINTER(BufTime);

    //read the actual values directly into SpecAxisBuf
    if (TwGetBufTimeFromShMem(BufTimeBuf, bufidx, writeidx) != TwSuccess)
    {
      UNPROTECT(1);
      return(R_NilValue);
    }
    UNPROTECT(1);
    return BufTime;
  }

  // Data stroage functions ---------------------------------------------------
  DLL_EXPORT SEXP AddLogEntry(SEXP LogEntryText)
  {
    //get the LogEntryText in a useable form
    char *cLogEntryText = RtoCstring(LogEntryText);

    unsigned __int64 cTime = 0;  // only "now" is supported due to difficulty with casting int64 values

    TwRetVal rv = TwAddLogEntry(cLogEntryText, cTime);
    SEXP rValue;
    PROTECT(rValue = NEW_INTEGER(1));
    int *pValue = INTEGER_POINTER(rValue);
    pValue[0] = rv;
    UNPROTECT(1);

    return rValue;
  }

  DLL_EXPORT SEXP AddAttributeInt(SEXP Object, SEXP AttributeName, SEXP Value)
  {
    //get the Object in a useable form
    char *cObject = RtoCstring(Object);

    //get the AttributeName in a useable form
    char *cAttributeName = RtoCstring(AttributeName);

    int cValue = INTEGER_VALUE(Value);

    TwRetVal rv = TwAddAttributeInt(cObject, cAttributeName, cValue);
    SEXP rValue;
    PROTECT(rValue = NEW_INTEGER(1));
    int *pValue = INTEGER_POINTER(rValue);
    pValue[0] = rv;
    UNPROTECT(1);

    return rValue;
  }

  DLL_EXPORT SEXP AddAttributeDouble(SEXP Object, SEXP AttributeName, SEXP Value)
  {
    //get the Object in a useable form
    char *cObject = RtoCstring(Object);

    //get the AttributeName in a useable form
    char *cAttributeName = RtoCstring(AttributeName);

    double cValue = NUMERIC_VALUE(Value);

    TwRetVal rv = TwAddAttributeDouble(cObject, cAttributeName, cValue);
    SEXP rValue;
    PROTECT(rValue = NEW_INTEGER(1));
    int *pValue = INTEGER_POINTER(rValue);
    pValue[0] = rv;
    UNPROTECT(1);

    return rValue;
  }

  DLL_EXPORT SEXP AddAttributeString(SEXP Object, SEXP AttributeName, SEXP Value)
  {
    //get the Object in a useable form
    char *cObject = RtoCstring(Object);

    //get the AttributeName in a useable form
    char *cAttributeName = RtoCstring(AttributeName);

    //get the Value in a useable form
    char *cValue = RtoCstring(Value);

    TwRetVal rv = TwAddAttributeString(cObject, cAttributeName, cValue);
    SEXP rValue;
    PROTECT(rValue = NEW_INTEGER(1));
    int *pValue = INTEGER_POINTER(rValue);
    pValue[0] = rv;
    UNPROTECT(1);

    return rValue;
  }

  DLL_EXPORT SEXP AddUserData(SEXP Location, SEXP NbrElements, SEXP ElementDescription, SEXP Data, SEXP CompressionLevel)
  {
    char *cLocation = RtoCstring(Location);

    int cNbrElements = INTEGER_VALUE(NbrElements);

    char *cElementDescription;
    if (strlen(CHAR(STRING_ELT(ElementDescription, 0))) == 0)
    {
      cElementDescription = NULL;
    }
    else
    {
      cElementDescription = RtoCstring(ElementDescription);
    }

    double *pData = NUMERIC_POINTER(Data);

    int cCompressionLevel = INTEGER_VALUE(CompressionLevel);

    TwRetVal rv = TwAddUserData(cLocation, cNbrElements, cElementDescription, pData, cCompressionLevel);
    SEXP rValue;
    PROTECT(rValue = NEW_INTEGER(1));
    int *pValue = INTEGER_POINTER(rValue);
    pValue[0] = rv;
    UNPROTECT(1);

    return rValue;
  }

  DLL_EXPORT SEXP AddUserDataMultiRow(SEXP Location, SEXP NbrElements, SEXP NbrRows, SEXP ElementDescription, SEXP Data, SEXP CompressionLevel)
  {
    char *cLocation = RtoCstring(Location);

    int cNbrElements = INTEGER_VALUE(NbrElements);

    int cNbrRows = INTEGER_VALUE(NbrRows);

    char *cElementDescription;
    if (strlen(CHAR(STRING_ELT(ElementDescription, 0))) == 0)
    {
      cElementDescription = NULL;
    }
    else
    {
      cElementDescription = RtoCstring(ElementDescription);
    }

    double *pData = NUMERIC_POINTER(Data);

    int cCompressionLevel = INTEGER_VALUE(CompressionLevel);

    TwRetVal rv = TwAddUserDataMultiRow(cLocation, cNbrElements, cNbrRows, cElementDescription, pData, cCompressionLevel);
    SEXP rValue;
    PROTECT(rValue = NEW_INTEGER(1));
    int *pValue = INTEGER_POINTER(rValue);
    pValue[0] = rv;
    UNPROTECT(1);

    return rValue;
  }

  DLL_EXPORT SEXP RegisterUserDataBuf(SEXP Location, SEXP NbrElements, SEXP ElementDescription, SEXP CompressionLevel)
  {
    char *cLocation = RtoCstring(Location);

    int cNbrElements = INTEGER_VALUE(NbrElements);

    char *cElementDescription;
    if (strlen(CHAR(STRING_ELT(ElementDescription, 0))) == 0)
    {
      cElementDescription = NULL;
    }
    else
    {
      cElementDescription = RtoCstring(ElementDescription);
    }

    int cCompressionLevel = INTEGER_VALUE(CompressionLevel);

    TwRetVal rv = TwRegisterUserDataBuf(cLocation, cNbrElements, cElementDescription, cCompressionLevel);
    SEXP rValue;
    PROTECT(rValue = NEW_INTEGER(1));
    int *pValue = INTEGER_POINTER(rValue);
    pValue[0] = rv;
    UNPROTECT(1);

    return rValue;
  }

  DLL_EXPORT SEXP RegisterUserDataWrite(SEXP Location, SEXP NbrElements, SEXP ElementDescription, SEXP CompressionLevel)
  {
    char *cLocation = RtoCstring(Location);

    int cNbrElements = INTEGER_VALUE(NbrElements);

    char *cElementDescription;
    if (strlen(CHAR(STRING_ELT(ElementDescription, 0))) == 0)
    {
      cElementDescription = NULL;
    }
    else
    {
      cElementDescription = RtoCstring(ElementDescription);
    }

    int cCompressionLevel = INTEGER_VALUE(CompressionLevel);

    TwRetVal rv = TwRegisterUserDataWrite(cLocation, cNbrElements, cElementDescription, cCompressionLevel);
    SEXP rValue;
    PROTECT(rValue = NEW_INTEGER(1));
    int *pValue = INTEGER_POINTER(rValue);
    pValue[0] = rv;
    UNPROTECT(1);

    return rValue;
  }

  DLL_EXPORT SEXP RegisterUserDataNoStore(SEXP Location, SEXP NbrElements, SEXP ElementDescription)
  {
    char *cLocation = RtoCstring(Location);

    int cNbrElements = INTEGER_VALUE(NbrElements);

    char *cElementDescription;
    if (strlen(CHAR(STRING_ELT(ElementDescription, 0))) == 0)
    {
      cElementDescription = NULL;
    }
    else
    {
      cElementDescription = RtoCstring(ElementDescription);
    }

    TwRetVal rv = TwRegisterUserDataNoStore(cLocation, cNbrElements, cElementDescription);
    SEXP rValue;
    PROTECT(rValue = NEW_INTEGER(1));
    int *pValue = INTEGER_POINTER(rValue);
    pValue[0] = rv;
    UNPROTECT(1);

    return rValue;
  }

  DLL_EXPORT SEXP UnregisterUserData(SEXP Location)
  {
    char *cLocation = RtoCstring(Location);

    TwRetVal rv = TwUnregisterUserData(cLocation);
    SEXP rValue;
    PROTECT(rValue = NEW_INTEGER(1));
    int *pValue = INTEGER_POINTER(rValue);
    pValue[0] = rv;
    UNPROTECT(1);

    return rValue;
  }

  DLL_EXPORT SEXP UpdateUserData(SEXP Location, SEXP NbrElements, SEXP Data)
  {
    char *cLocation = RtoCstring(Location);

    int cNbrElements = INTEGER_VALUE(NbrElements);

    double *pData = NUMERIC_POINTER(Data);

    TwRetVal rv = TwUpdateUserData(cLocation, cNbrElements, pData);
    SEXP rValue;
    PROTECT(rValue = NEW_INTEGER(1));
    int *pValue = INTEGER_POINTER(rValue);
    pValue[0] = rv;
    UNPROTECT(1);

    return rValue;
  }

  DLL_EXPORT SEXP ReadRegUserData(SEXP Location, SEXP NbrElements)
  {
    char *cLocation = RtoCstring(Location);

    int cNbrElements = INTEGER_VALUE(NbrElements);

    //prepare return array
    SEXP Data;
    double *pData;
    PROTECT(Data = NEW_NUMERIC(cNbrElements));
    pData = NUMERIC_POINTER(Data);
    UNPROTECT(1);

    TwRetVal rv = TwReadRegUserData(cLocation, cNbrElements, pData);
    if (rv != TwSuccess)
    {
      return(R_NilValue);
    }

    return Data;

  }

  DLL_EXPORT SEXP QueryRegUserDataSize(SEXP Location)
  {
    char *cLocation = RtoCstring(Location);

    //prepare return array
    SEXP NbrElements;
    int *pNbrElements;
    PROTECT(NbrElements = NEW_INTEGER(1));
    pNbrElements = INTEGER_POINTER(NbrElements);
    UNPROTECT(1);

    TwRetVal rv = TwQueryRegUserDataSize(cLocation, pNbrElements);
    if (rv != TwSuccess)
    {
      return(R_NilValue);
    }

    return NbrElements;

  }

  DLL_EXPORT SEXP GetRegUserDataSources()
  {
    int arrayLength = 0;

    if (TwGetRegUserDataSources(&arrayLength, NULL, NULL, NULL) != TwValueAdjusted)
    {
      return(R_NilValue);
    }

    //prepare return array
    SEXP location, nbrElements, type, list, list_names;

    char *p_location = new char[256 * arrayLength];
    memset(p_location, 0, 256 * arrayLength);
    int *p_nbrElements,  *p_type;

    PROTECT(nbrElements = NEW_INTEGER(arrayLength));
    p_nbrElements = INTEGER_POINTER(nbrElements);
    UNPROTECT(1);

    PROTECT(type = NEW_INTEGER(arrayLength));
    p_type = INTEGER_POINTER(type);
    UNPROTECT(1);

    TwRetVal rv = TwGetRegUserDataSources(&arrayLength, p_location, p_nbrElements, p_type);
    if (rv != TwSuccess)
    {
      delete[] p_location;
      p_location = nullptr;
      return(R_NilValue);
    }

    PROTECT(location = allocVector(STRSXP, arrayLength));
    for (int i = 0; i < arrayLength; ++i)
    {
      SET_STRING_ELT(location, i, mkChar(&p_location[i * 256]));
    }
    UNPROTECT(1);

    delete[] p_location;
    p_location = nullptr;

    PROTECT(list = allocVector(VECSXP, 3));
    SET_VECTOR_ELT(list, 0, location);
    SET_VECTOR_ELT(list, 1, nbrElements);
    SET_VECTOR_ELT(list, 2, type);
    UNPROTECT(1);

    // list names
    char const* names[3] = { "location", "nbrElements", "type" };
    PROTECT(list_names = allocVector(STRSXP, 3));
    int i;
    for (i = 0; i < 3; i++)
      SET_STRING_ELT(list_names, i, mkChar(names[i]));
    setAttrib(list, R_NamesSymbol, list_names);
    UNPROTECT(1);

    return list;
  }

  DLL_EXPORT SEXP GetRegUserDataDesc(SEXP Location)
  {
    char *cLocation = RtoCstring(Location);

    int nbrElements = 0;

    if (TwGetRegUserDataDesc(cLocation, &nbrElements, NULL) != TwValueAdjusted)
    {
      return(R_NilValue);
    }

    //prepare return array
    SEXP ElementDescription;

    char *cElementDescription = new char[256 * nbrElements];
    memset(cElementDescription, 0, 256 * nbrElements);

    TwRetVal rv = TwGetRegUserDataDesc(cLocation, &nbrElements, cElementDescription);
    if (rv != TwSuccess)
    {
      delete[] cElementDescription;
      cElementDescription = nullptr;
      return(R_NilValue);
    }

    PROTECT(ElementDescription = allocVector(STRSXP, nbrElements));
    for (int i = 0; i < nbrElements; ++i)
    {
      SET_STRING_ELT(ElementDescription, i, mkChar(&cElementDescription[i * 256]));
    }
    UNPROTECT(1);

    delete[] cElementDescription;
    cElementDescription = nullptr;

    return ElementDescription;
  }

  DLL_EXPORT SEXP KeepFileOpen(SEXP keepOpen)
  {
    bool ko = LOGICAL_VALUE(keepOpen);

    TwRetVal rv = TwKeepFileOpen(ko);
    SEXP rValue;
    PROTECT(rValue = NEW_INTEGER(1));
    int *pValue = INTEGER_POINTER(rValue);
    pValue[0] = rv;
    UNPROTECT(1);

    return rValue;

  }

  // TPS control functions ----------------------------------------------------
  DLL_EXPORT SEXP TpsConnect()
  {
    TwRetVal rv = TwTpsConnect();
    SEXP rValue;
    PROTECT(rValue = NEW_INTEGER(1));
    int *pValue = INTEGER_POINTER(rValue);
    pValue[0] = rv;
    UNPROTECT(1);

    return rValue;
  }

  DLL_EXPORT SEXP TpsConnect2(SEXP ip, SEXP type)
  {
    int xtype = INTEGER_VALUE(type);

    //get the ip in a useable form
    if (xtype == 1)
    {
      char *cFilename = RtoCstring(ip);

      TwRetVal rv = TwTpsConnect2(cFilename, xtype);
      SEXP rValue;
      PROTECT(rValue = NEW_INTEGER(1));
      int *pValue = INTEGER_POINTER(rValue);
      pValue[0] = rv;
      UNPROTECT(1);

      return rValue;
    }
    else if (xtype == 0)
    {
      TwRetVal rv = TwTpsConnect();
      SEXP rValue;
      PROTECT(rValue = NEW_INTEGER(1));
      int *pValue = INTEGER_POINTER(rValue);
      pValue[0] = rv;
      UNPROTECT(1);

      return rValue;
    }
    else
    { return(R_NilValue); }
  }

  DLL_EXPORT SEXP TpsDisconnect()
  {
    TwRetVal rv = TwTpsDisconnect();
    SEXP rValue;
    PROTECT(rValue = NEW_INTEGER(1));
    int *pValue = INTEGER_POINTER(rValue);
    pValue[0] = rv;
    UNPROTECT(1);

    return rValue;
  }

  DLL_EXPORT SEXP TpsGetMonitorValue(SEXP moduleCode)
  {
    int module = INTEGER_VALUE(moduleCode);

    double value;
    TwTpsGetMonitorValue(module, &value);

    SEXP rValue;
    PROTECT(rValue = NEW_NUMERIC(1));
    double *pValue = NUMERIC_POINTER(rValue);
    pValue[0] = value;
    UNPROTECT(1);

    return rValue;
  }

  DLL_EXPORT SEXP TpsGetTargetValue(SEXP moduleCode)
  {
    int module = INTEGER_VALUE(moduleCode);

    double value;
    TwTpsGetTargetValue(module, &value);

    SEXP rValue;
    PROTECT(rValue = NEW_NUMERIC(1));
    double *pValue = NUMERIC_POINTER(rValue);
    pValue[0] = value;
    UNPROTECT(1);

    return rValue;
  }

  DLL_EXPORT SEXP TpsGetLastSetValue(SEXP moduleCode)
  {
    int module = INTEGER_VALUE(moduleCode);

    double value;
    TwRetVal rv = TwTpsGetLastSetValue(module, &value);
    if (rv != TwSuccess)
    {
      return(R_NilValue);
    }

    SEXP rValue;
    PROTECT(rValue = NEW_NUMERIC(1));
    double *pValue = NUMERIC_POINTER(rValue);
    pValue[0] = value;
    UNPROTECT(1);

    return rValue;
  }

  DLL_EXPORT SEXP TpsSetTargetValue(SEXP moduleCode, SEXP value)
  {
    int module = INTEGER_VALUE(moduleCode);
    double cvalue = NUMERIC_VALUE(value);
    TwRetVal rv = TwTpsSetTargetValue(module, cvalue);

    SEXP rValue;
    PROTECT(rValue = NEW_INTEGER(1));
    int *pValue = INTEGER_POINTER(rValue);
    pValue[0] = rv;
    UNPROTECT(1);

    return rValue;
  }

  DLL_EXPORT SEXP TpsGetNbrModules()
  {
    int nbrModules;
    TwRetVal rv = TwTpsGetNbrModules(&nbrModules);
    if (rv != TwSuccess)
    {
      return(R_NilValue);
    }
    SEXP rValue;
    PROTECT(rValue = NEW_INTEGER(1));
    int *pValue = INTEGER_POINTER(rValue);
    pValue[0] = nbrModules;
    UNPROTECT(1);

    return rValue;
  }

  DLL_EXPORT SEXP TpsGetModuleCodes()
  {
    //get number of modules
    int nbrModules;
    TwRetVal rv1 = TwTpsGetNbrModules(&nbrModules);
    if (rv1 != TwSuccess)
    {
      return(R_NilValue);
    }

    //prepare return array
    SEXP moduleCodes;
    int *moduleCodeBuffer;
    int bufferLength = nbrModules;
    PROTECT(moduleCodes = NEW_INTEGER(bufferLength));
    moduleCodeBuffer = INTEGER_POINTER(moduleCodes);
    TwRetVal rv2 = TwTpsGetModuleCodes(moduleCodeBuffer, bufferLength);
    if (rv2 != TwSuccess)
    {
      return(R_NilValue);
    }
    UNPROTECT(1);
    return moduleCodes;
  }

  DLL_EXPORT SEXP TpsInitialize()
  {
    TwRetVal rv = TwTpsInitialize();
    SEXP rValue;
    PROTECT(rValue = NEW_INTEGER(1));
    int *pValue = INTEGER_POINTER(rValue);
    pValue[0] = rv;
    UNPROTECT(1);

    return rValue;
  }

  DLL_EXPORT SEXP TpsSetAllVoltages()
  {
    TwRetVal rv = TwTpsSetAllVoltages();
    SEXP rValue;
    PROTECT(rValue = NEW_INTEGER(1));
    int *pValue = INTEGER_POINTER(rValue);
    pValue[0] = rv;
    UNPROTECT(1);

    return rValue;
  }

  DLL_EXPORT SEXP TpsShutdown()
  {
    TwRetVal rv = TwTpsShutdown();
    SEXP rValue;
    PROTECT(rValue = NEW_INTEGER(1));
    int *pValue = INTEGER_POINTER(rValue);
    pValue[0] = rv;
    UNPROTECT(1);

    return rValue;
  }

  DLL_EXPORT SEXP TpsGetStatus()
  {
    int status;
    TwRetVal rv = TwTpsGetStatus(&status);
    if (rv != TwSuccess)
    {
      return(R_NilValue);
    }
    SEXP rValue;
    PROTECT(rValue = NEW_INTEGER(1));
    int *pValue = INTEGER_POINTER(rValue);
    pValue[0] = status;
    UNPROTECT(1);

    return rValue;
  }

  DLL_EXPORT SEXP TpsLoadSetFile(SEXP setFile)
  {
    //get the filename in a useable form
    char *cFilename = RtoCstring(setFile);

    TwRetVal rv = TwTpsLoadSetFile(cFilename);
    SEXP rValue;
    PROTECT(rValue = NEW_INTEGER(1));
    int *pValue = INTEGER_POINTER(rValue);
    pValue[0] = rv;
    UNPROTECT(1);

    return rValue;
  }

  DLL_EXPORT SEXP TpsSaveSetFile(SEXP setFile)
  {
    //get the filename in a useable form
    char *cFilename = RtoCstring(setFile);

    TwRetVal rv = TwTpsSaveSetFile(cFilename);
    SEXP rValue;
    PROTECT(rValue = NEW_INTEGER(1));
    int *pValue = INTEGER_POINTER(rValue);
    pValue[0] = rv;
    UNPROTECT(1);

    return rValue;
  }

  DLL_EXPORT SEXP TpsGetActiveFilament()
  {
    int filament;
    TwTpsGetActiveFilament(&filament);

    SEXP rValue;
    PROTECT(rValue = NEW_INTEGER(1));
    int *pValue = INTEGER_POINTER(rValue);
    pValue[0] = filament;
    UNPROTECT(1);

    return rValue;
  }

  DLL_EXPORT SEXP TpsSetActiveFilament(SEXP activeFilament)
  {
    int filament = INTEGER_VALUE(activeFilament);

    TwRetVal rv = TwTpsSetActiveFilament(filament);

    SEXP rValue;
    PROTECT(rValue = NEW_INTEGER(1));
    int *pValue = INTEGER_POINTER(rValue);
    pValue[0] = rv;
    UNPROTECT(1);

    return rValue;
  }

  DLL_EXPORT SEXP TpsGetModuleLimits(SEXP moduleCode)
  {
    int module = INTEGER_VALUE(moduleCode);

    double minLimit;
    double maxLimit;
    TwRetVal rv = TwTpsGetModuleLimits(module, &minLimit, &maxLimit);
    if (rv != TwSuccess)
    {
      return(R_NilValue);
    }

    SEXP rValue;
    PROTECT(rValue = NEW_NUMERIC(2));
    double *pValue = NUMERIC_POINTER(rValue);
    pValue[0] = minLimit;
    pValue[1] = maxLimit;
    UNPROTECT(1);

    return rValue;
  }

  DLL_EXPORT SEXP TpsChangeIonMode(SEXP ionMode)
  {
    int mode = INTEGER_VALUE(ionMode);

    TwRetVal rv = TwTpsChangeIonMode(mode);
    SEXP rValue;
    PROTECT(rValue = NEW_INTEGER(1));
    int *pValue = INTEGER_POINTER(rValue);
    pValue[0] = rv;
    UNPROTECT(1);

    return rValue;
  }

  // Additional functions -----------------------------------------------------
  DLL_EXPORT SEXP KeepSharedMemMapped(SEXP keepMapped)
  {
    // Note: This is part of TwGetSharedMemory, but here extracted as an
    // idependent function.

    TSharedMemoryPointer pShMem;
    bool keep = LOGICAL_VALUE(keepMapped);

    TwRetVal rv = TwGetSharedMemory(&pShMem, keep);
    SEXP rValue;
    PROTECT(rValue = NEW_INTEGER(1));
    int *pValue = INTEGER_POINTER(rValue);
    pValue[0] = rv;
    UNPROTECT(1);

    return rValue;
  }

  DLL_EXPORT SEXP SiProcessSpectrumFromShMem(SEXP specType, SEXP BufIndex)
  {
    int bufidx = INTEGER_VALUE(BufIndex);
    int spect = INTEGER_VALUE(specType);

    //get descriptor of file
    TSharedMemoryDesc desc;
    if (TwGetDescriptor(&desc) != TwSuccess)
    {
      return(R_NilValue);
    }

    int specLen = desc.NbrSamples;

    //allocate temp float buffer
    float *spectrum = new float[specLen];
    TwGetTofSpectrumFromShMem(spectrum, 0, 0, bufidx, false);

    //prepare return array
    SEXP blFromData, thrFromData, list, list_names;
    float cblFromData;
    float cthrFromData;

    if (TwSiProcessSpectrum(spectrum, specLen, spect, &cblFromData, &cthrFromData) != TwSuccess)
    {
      delete[] spectrum;
      spectrum = nullptr;
      return(R_NilValue);
    }

    delete[] spectrum;
    spectrum = nullptr;

    PROTECT(blFromData = NEW_NUMERIC(1));
    double *p_blFromData = NUMERIC_POINTER(blFromData);
    p_blFromData[0] = cblFromData;
    UNPROTECT(1);
    PROTECT(thrFromData = NEW_NUMERIC(1));
    double *p_thrFromData = NUMERIC_POINTER(thrFromData);
    p_thrFromData[0] = cthrFromData;
    UNPROTECT(1);

    // Creating a list with 2 vector elements:
    PROTECT(list = allocVector(VECSXP, 2));
    SET_VECTOR_ELT(list, 0, blFromData);
    SET_VECTOR_ELT(list, 1, thrFromData);

    // Creating a character string vector of the "names" attribute of the
    // objects in out list:
    char const* names[2] = { "baseline", "threshold" };
    PROTECT(list_names = allocVector(STRSXP, 2));
    int i;
    for (i = 0; i < 2; i++)
    {
      SET_STRING_ELT(list_names, i, mkChar(names[i]));
    }

    // and attaching the vector names:
    setAttrib(list, R_NamesSymbol, list_names);
    UNPROTECT(2);
    return list;

  }
#endif

  DLL_EXPORT SEXP tof(SEXP toftype, SEXP drift, SEXP pulse, SEXP mass, SEXP x, SEXP v)
  {
    const double amu = 1.660538921e-27;  // atomic mass unit(kg)
    const double e = 1.60217657e-19;  // elementary charge(C)

    //get the Parameter in a useable form
    char *cToftype = RtoCstring(toftype);

    double cdrift = NUMERIC_VALUE(drift);
    double cpulse = NUMERIC_VALUE(pulse);
    double cmass = NUMERIC_VALUE(mass);

    double *cx = NUMERIC_POINTER(x);
    double *cv = NUMERIC_POINTER(v);

    int nbr;
    nbr = length(x);
    if (nbr != length(v))
    {
      return(R_NilValue);
    }

    double Vdrift = -std::abs(cdrift);
    double Vpush = std::abs(cpulse);
    double Vpull = -Vpush;

    double d1 = 0.004 + 0.0035;  // pull to push distance
    //Note: here d1 covers region from push plate to pull grid and d2 = 0
    // -> tof calculation also works for k<0, i.e. with ions which are beyond reference grid

    double u1 = Vpush - Vpull;
    double u3 = Vpull - Vdrift;
    double d3, d4, d5, d6, d7;
    double u5, u6;
    double u7 = 0;  // initialized to prevent -Wmaybe-uninitialized warning.

    if (strcmp(cToftype, "HTOF") == 0)
    {
      d3 = 0.0065;
      d4 = 0.507 + 0.511;
      d5 = 0.0165;
      d6 = 0.0713 - d5;
    }
    else if (strcmp(cToftype, "HTOF-W") == 0)
    {
      d3 = 0.0065;
      d4 = 0.507 + 0.511 + 2 * 0.5125;
      d5 = 0.0165;
      d6 = 0.0713 - d5;
      d7 = 0.0085;
      u7 = 2 * Vpush;
    }
    else if (strcmp(cToftype, "LTOF") == 0)
    {
      d3 = 0.014;
      d4 = 1.0465 + 1.052;
      d5 = 0.0385;
      d6 = 0.1593 - d5;
    }
    else if (strcmp(cToftype, "CTOF") == 0)
    {
      d3 = 0.006;
      d4 = 2 * 0.1435;
      d5 = 0.017;
      d6 = 0.0165;
    }
    else
    {
      return(R_NilValue);
    }

    // calculate u5 and u6
    double x0 = d1 - 0.001;
    double k0 = x0 / d1;
    double p0 = k0 + u3 / u1;
    double a, b;
    if (strcmp(cToftype, "HTOF-W") == 0)
    {
      a = d1 / u1*pow(k0, -0.5) + d3 / u3*(pow(p0, -0.5) - pow(k0, -0.5)) - d4 / u1 / 2 * pow(p0, -1.5) + 2 * d7 / u7*pow(p0, -0.5);
      b = d1 / u1 / 2 * pow(k0, -1.5) + d3 / u3 / 2 * (pow(p0, -1.5) - pow(k0, -1.5)) - d4 / u1 * 3 / 4 * pow(p0, -2.5) + d7 / u7*pow(p0, -1.5);
      u5 = (a - 2 * p0*b + 4 * d5 / u1*pow(p0, -1.5)) / (-2 * b / u1);
      u6 = (-4 * d6*pow(p0 - u5 / u1, -0.5)) / (a + 4 * d5 / u5*(pow(p0, -0.5) - pow(p0 - u5 / u1, -0.5)));
    }
    else
    {
      a = d1 / u1*pow(k0, -0.5) + d3 / u3*(pow(p0, -0.5) - pow(k0, -0.5)) - d4 / u1 / 2 * pow(p0, -1.5);
      b = d1 / u1 / 2 * pow(k0, -1.5) + d3 / u3 / 2 * (pow(p0, -1.5) - pow(k0, -1.5)) - d4 / u1 * 3 / 4 * pow(p0, -2.5);
      u5 = (a - 2 * p0*b + 2 * d5 / u1*pow(p0, -1.5)) / (-2 * b / u1);
      u6 = (-2 * d6*pow(p0 - u5 / u1, -0.5)) / (a + 2 * d5 / u5*(pow(p0, -0.5) - pow(p0 - u5 / u1, -0.5)));
    }

    //allocate xi and k buffers
    double *xi = new double[nbr];
    double *k = new double[nbr];

    for (int i = 0; i < nbr; ++i)
    {
      xi[i] = cv[i] * sqrt((cmass*amu) / (2 * e*u1));
      k[i] = (x0 - cx[i]) / d1;
    }

    //prepare return value
    SEXP timeOfFlight;
    double *ctimeOfFlight;
    PROTECT(timeOfFlight = NEW_NUMERIC(nbr));
    ctimeOfFlight = NUMERIC_POINTER(timeOfFlight);

    // calculate time-of-flight
    if (strcmp(cToftype, "HTOF-W") == 0)
    {
      for (int i = 0; i < nbr; ++i)
      {
        ctimeOfFlight[i] = sqrt((cmass*amu) / (2 * e)) *
          (2 * d1 / sqrt(u1)*(sqrt(xi[i] * xi[i] + k[i]) - xi[i]) +
          2 * d3 / u3*(sqrt(u1*xi[i] * xi[i] + k[i] * u1 + u3) - sqrt(u1*xi[i] * xi[i] + k[i] * u1)) +
          d4 / sqrt(u1*xi[i] * xi[i] + k[i] * u1 + u3) +
          2 * 4 * d5 / u5*(sqrt(u1*xi[i] * xi[i] + k[i] * u1 + u3) - sqrt(u1*xi[i] * xi[i] + k[i] * u1 + u3 - u5)) +
          2 * 4 * d6 / u6*sqrt(u1*xi[i] * xi[i] + k[i] * u1 + u3 - u5) +
          4 * d7 / u7*sqrt(u1*xi[i] * xi[i] + k[i] * u1 + u3));
      }
    }
    else
    {
      for (int i = 0; i < nbr; ++i)
      {
        ctimeOfFlight[i] = sqrt((cmass*amu) / (2 * e)) *
          (2 * d1 / sqrt(u1)*(sqrt(xi[i] * xi[i] + k[i]) - xi[i]) +
          2 * d3 / u3*(sqrt(u1*xi[i] * xi[i] + k[i] *u1 + u3) - sqrt(u1*xi[i] * xi[i] + k[i] *u1)) +
          d4 / sqrt(u1*xi[i] * xi[i] + k[i] *u1 + u3) +
          4 * d5 / u5*(sqrt(u1*xi[i] * xi[i] + k[i] *u1 + u3) - sqrt(u1*xi[i] * xi[i] + k[i] *u1 + u3 - u5)) +
          4 * d6 / u6*sqrt(u1*xi[i] * xi[i] + k[i] *u1 + u3 - u5));
      }
    }

    delete[] xi;
    xi = nullptr;
    delete[] k;
    k = nullptr;

    UNPROTECT(1);

    return timeOfFlight;
  }

  DLL_EXPORT SEXP EventList2TofSpec(SEXP events, SEXP clockPeriod, SEXP sampleInterval, SEXP nbrSamples)
  {

    unsigned int n = length(events);
    double *pevents = NUMERIC_POINTER(events);
    double cclockPeriod = NUMERIC_VALUE(clockPeriod);
    double csampleInterval = NUMERIC_VALUE(sampleInterval);
    int cnbrSamples = INTEGER_VALUE(nbrSamples);

    //prepare return array
    SEXP spectrum;
    double *cspectrum;
    PROTECT(spectrum = NEW_NUMERIC(cnbrSamples));
    cspectrum = NUMERIC_POINTER(spectrum);

    for (unsigned int i = 0; i < n; ++i) {
      const unsigned int timestamp = (unsigned int)pevents[i] & 0xFFFFFF;
      //convert timestamp to a sample index
      const unsigned int sampleIndex = (unsigned int)(timestamp / (csampleInterval / cclockPeriod) + 0.5);
      const unsigned int dataLength = (unsigned int)pevents[i] >> 24;
      if (dataLength == 0) { //TDC data -> each event is 1 count
        cspectrum[sampleIndex] += 1.0f;
      }
      else { //ADC data (timestamp is time of first sample in packet)
        for (unsigned int j = 0; j < dataLength; ++j) {
          ++i;
          unsigned int mask = (unsigned int)pevents[i];
          float* adcData = (float*)& mask;
          cspectrum[sampleIndex + j] += *adcData;
        }
      }
    }

    UNPROTECT(1);

    return spectrum;
  }

  DLL_EXPORT SEXP DecodeEventList(SEXP events, SEXP clockPeriod, SEXP sampleInterval)
  {

    unsigned int n = length(events);
    double *pevents = NUMERIC_POINTER(events);
    double cclockPeriod = NUMERIC_VALUE(clockPeriod);
    double csampleInterval = NUMERIC_VALUE(sampleInterval);

    //prepare return list
    int listLength = 2;
    char const* names[2] = { "sampleindex", "value" };

    SEXP sampleindex, value, list, list_names;

    PROTECT(list = allocVector(VECSXP, listLength));

    // temporary arrays
    double *temp1 = new double[n];
    double *temp2 = new double[n];

    unsigned int k = 0;

    for (unsigned int i = 0; i < n; ++i) {
      const unsigned int timestamp = (unsigned int)pevents[i] & 0xFFFFFF;
      const unsigned int dataLength = (unsigned int)pevents[i] >> 24;
      if (dataLength == 0) {  // TDC data -> each event is 1 count
        // convert timestamp to a sample index
        temp1[k] = (int)(timestamp * cclockPeriod / csampleInterval + 0.5);
        temp2[k] = 1;
        k++;
      }
      else { // ADC data or "Ndigo TDC" data (timestamp is time of first sample in packet)
        for (unsigned int j = 0; j < dataLength; ++j) {
          ++i;
          const unsigned int mask = (unsigned int)pevents[i];
          float* adcData = (float*)& mask;
          temp1[k] = (int)(timestamp * cclockPeriod / csampleInterval + 0.5) + j;
          temp2[k] = *adcData;
          k++;
        }
      }
    }

    PROTECT(sampleindex = NEW_NUMERIC(k));
    PROTECT(value = NEW_NUMERIC(k));
    double *csampleindex = NUMERIC_POINTER(sampleindex);
    double *cvalue = NUMERIC_POINTER(value);

    for (unsigned int i = 0; i < k; i++) {
      csampleindex[i] = temp1[i];
      cvalue[i] = temp2[i];
    }

    delete[] temp1;
    temp1 = nullptr;
    delete[] temp2;
    temp2 = nullptr;

    SET_VECTOR_ELT(list, 0, sampleindex);
    SET_VECTOR_ELT(list, 1, value);

    // Create character string vector of the "names" attribute
    PROTECT(list_names = allocVector(STRSXP, listLength));
    int i;
    for (i = 0; i < listLength; i++)
      SET_STRING_ELT(list_names, i, mkChar(names[i]));
    // attach vector names
    setAttrib(list, R_NamesSymbol, list_names);
    UNPROTECT(4);

    return list;
  }

#ifdef __cplusplus
}
#endif
