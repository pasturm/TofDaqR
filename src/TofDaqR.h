#ifndef TOFDAQR_H
#define TOFDAQR_H

// convert std::string to char*
char* StringToChar(std::string str);

// convert TwRetVal to std::string
std::string TranslateReturnValue(TwRetVal rv);

#endif
