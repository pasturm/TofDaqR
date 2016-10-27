# TofDaqR
R Interface to TOFWERK TofDaq API.

The TofDaqR package provides an interface to the TOFWERK TofDaq API (TofDaqDll, 
TwH5Dll, TwToolDll for Windows, and libtwh5, libtwtool for macOS and Linux). 
TofDaqR functions can be used to communicate with the TofDaq recorder 
application, to access data in realtime and to control the TOF Power Supply. 
Further, there are functions to access and modify recorded data stored in HDF5 
files and additional "tool" functions, e.g. for peak fitting.

## Version history
See [NEWS](https://github.com/pasturm/TofDaqR/blob/master/NEWS).

## Installation
```
install.packages("devtools")
devtools::install_github("pasturm/TofDaqR")
```

Notes:

* On Windows, [Rtools](https://cran.r-project.org/bin/windows/Rtools/) is required to compile and install TofDaqR.

* On Mac, [Xcode](https://developer.apple.com/xcode/) is required to compile and install TofDaqR.

* On Linux, `sudo apt-get install r-base-dev` (or similar) is required to compile and install TofDaqR.
