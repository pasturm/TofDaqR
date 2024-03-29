# TofDaqR
<!-- badges: start -->
[![Project Status: Active – The project has reached a stable, usable state and is being actively developed.](https://www.repostatus.org/badges/latest/active.svg)](https://www.repostatus.org/#active)
<!-- badges: end -->


## R Interface to the TOFWERK TofDaq API
The TofDaqR package provides a R interface to the [TOFWERK TofDaq API](https://www.tofwerk.com/software/tofdaq/), which consists of libraries for communication with the TofDAQ recorder application, data file access and general (time-of-flight) mass spectrometry related utility functions. 

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

On Windows:
```r
install.packages("https://github.com/pasturm/TofDaqR/releases/download/v0.3.11/TofDaqR_0.3.11.zip", repos = NULL)
```

Alternatively, the latest development version can be installed from source on Windows, macOS (with Intel-based processors) and Linux:
```r
if (!require("remotes")) { install.packages("remotes") }
remotes::install_github("pasturm/TofDaqR", clean = TRUE)
```
Installing from source requires [Rtools](https://cran.r-project.org/bin/windows/Rtools/) on Windows, [Xcode Command Line Tools](https://developer.apple.com/xcode/resources/) (`xcode-select --install`) on macOS and `sudo apt-get install r-base-dev` (or similar) on Linux. 

Please note: the TofDaq API libraries are not (yet) available for ARM-based processors.

## Documentation
```r
help(package = "TofDaqR")  # Index of help pages
help("API-documentation", package = "TofDaqR")  # Links to the TOFWERK TofDaq API documentation
```

## Release notes
See the [NEWS file](https://github.com/pasturm/TofDaqR/blob/master/NEWS.md) for latest release notes.
