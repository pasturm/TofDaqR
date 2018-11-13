# TofDaqR
R Interface to the [TOFWERK TofDaq API](https://www.tofwerk.com/software/tofdaq/).

The TofDaqR package provides a R interface to the TOFWERK TofDaq API, which consists of libraries for communication with the TofDAQ recorder application, data file access and general (time-of-flight) mass spectrometry related utility functions. 

Features:

* Acquisition setup and control
* Real-time data access
* Adding (structured) custom data to the data file alongside the TOF data
* Control of the TOF Power Supply
* Quick access to TOFWERK HDF5 data files without the need to study the file format details
* Add additional data to existing data files (e.g. post-processing results)
* Peak fitting functions
* Mass calibration
* Chemistry functions (molecular mass and isotope pattern calculation)
* Single ion analysis functions

## Installation
To install the latest binary version on Windows, run:
```
install.packages("https://github.com/pasturm/TofDaqR/releases/download/v0.3.8/TofDaqR_0.3.8.zip", repos = NULL)
```
On macOS run:
```
install.packages("https://github.com/pasturm/TofDaqR/releases/download/v0.3.8/TofDaqR_0.3.8.904.tgz", repos = NULL)
```

### Source installation
To install the latest development version from source (on Windows, macOS and Linux), run:
```
install.packages("devtools")
devtools::install_github("pasturm/TofDaqR", clean = TRUE)
```
Installing from source requires [Rtools](https://cran.r-project.org/bin/windows/Rtools/) on Windows,  [Xcode](https://developer.apple.com/xcode/) on macOS and `sudo apt-get install r-base-dev` (or similar) on Linux. 

### Release notes
See [NEWS](https://github.com/pasturm/TofDaqR/blob/master/NEWS.md) for release notes.
