# TofDaqR
R Interface to the TOFWERK TofDaq API.

The TofDaqR package provides a R interface to the TOFWERK TofDaq API (TofDaqDll, 
TwH5Dll, TwToolDll for Windows, and libtwh5, libtwtool for macOS and Linux). 
TofDaqR functions can be used to communicate with the TofDaq recorder 
application, to access data in realtime and to control the TOF Power Supply. 
Further, there are post-processing functions to access and modify recorded data 
stored in HDF5 files and additional "tool" functions, e.g. for peak fitting.

## Installation
### Binary
For Windows (64-bit) a binary package is available at [releases](https://github.com/pasturm/TofDaqR/releases). Download and install 
the zip file.

### Source
```
install.packages("devtools")
devtools::install_github("pasturm/TofDaqR", args = "--clean")
```

Notes:

* On Windows, [Rtools](https://cran.r-project.org/bin/windows/Rtools/) is required to install TofDaqR from source.

* On Mac, [Xcode](https://developer.apple.com/xcode/) is required to install TofDaqR.

* On Linux, `sudo apt-get install r-base-dev` (or similar) is required to install TofDaqR.

## Version history
See [NEWS](https://github.com/pasturm/TofDaqR/blob/master/NEWS).

