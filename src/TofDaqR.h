#ifndef TOFDAQR_H
#define TOFDAQR_H

// convert R string to C string
char* RtoCstring(SEXP rstring);

// convert TwRetVal to String
String TwRetValString(TwRetVal rv);

#endif
