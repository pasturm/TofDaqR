# TofDaqR
R Interface to the [TOFWERK TofDaq API](https://www.tofwerk.com/software/tofdaq/).

The TofDaqR package provides a R interface to the TOFWERK TofDaq API, which consists of libraries for communication with the TofDAQ recorder application, data file access and general (time-of-flight) mass spectrometry related utility functions. 

Features:

* Acquisition setup and control
* Real-time data access
* Control of the TOF Power Supply
* Quick access to TOFWERK HDF5 data files
* Add additional data to existing data files
* Peak fitting functions
* Mass calibration
* Chemistry functions (molecular mass and isotope pattern calculation)
* Single ion analysis functions

## Installation
To install the latest binary version on Windows, run
```
install.packages("https://github.com/pasturm/TofDaqR/releases/download/v0.3.6/TofDaqR_0.3.6.zip", repos = NULL)
```

To install the latest development version from source (on Windows, macOS or Linux), run
```
install.packages("devtools")
devtools::install_github("pasturm/TofDaqR", args = "--clean")
```
Installing from source requires [Rtools](https://cran.r-project.org/bin/windows/Rtools/) on Windows,  [Clang](https://cran.r-project.org/bin/macosx/) on macOS and `sudo apt-get install r-base-dev` (or similar) on Linux. 

See [NEWS](https://github.com/pasturm/TofDaqR/blob/master/NEWS) for release notes.
