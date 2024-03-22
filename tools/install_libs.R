# Script to install the TofDaq API libraries (invoked by Makevars).

path = "https://github.com/pasturm/TofDaqR/blob/master/TofDaqAPI/"
TofDaqAPI = "TofDaq_1.99r1369_API_20210930"

# download TofDaq API files ----------------------------------------------------
if (dir.exists(file.path("../..", TofDaqAPI))) {  # debug
  print("*** install_libs.R: using local TofDaqAPI")
  tmpdir = file.path("../..", TofDaqAPI)
} else {
  tmpzip = tempfile()
  tmpdir = file.path(tempdir(), TofDaqAPI)
  dir.create(tmpdir)
  download.file(paste0(path, TofDaqAPI, ".zip"), tmpzip)
  unzip(tmpzip, exdir = tmpdir)
}

# copy include to tools/ -------------------------------------------------------
rv = file.copy(file.path(tmpdir, "include"), "../tools/", recursive = TRUE)

# add ChangeIonMode (which does not exist by default) to TofDaqDll.h
TofDaqDll.h = readLines(file("../tools/include/TofDaqDll.h"))
n = which(TofDaqDll.h=="TOFWERK_DAQ_API bool TwSaturationWarning(void);")
TofDaqDll.h[n]  = "TOFWERK_DAQ_API bool TwSaturationWarning(void);
////////////////////////////////////////////////////////////////////////////////
TOFWERK_DAQ_API TwRetVal TwTpsChangeIonMode(int ionMode);
////////////////////////////////////////////////////////////////////////////////"
writeLines(TofDaqDll.h, con = "../tools/include/TofDaqDll.h" )

# copy the appropriate dynamic libraries to R_PACKAGE_DIR/libs -----------------
R_PACKAGE_DIR = commandArgs(TRUE)  # passed by Makevars
R_ARCH = Sys.getenv("R_ARCH")
libarch = if (nzchar(R_ARCH)) paste0("libs", R_ARCH) else "libs"
dest = file.path(R_PACKAGE_DIR, libarch)

if (Sys.info()["sysname"] == "Windows") {
  if (R_ARCH == "/x64") {
    files = Sys.glob(file.path(tmpdir,  "bin/windows_x64", "*.dll"))
  } else if (R_ARCH == "/i386") {
    files = Sys.glob(file.path(tmpdir, "bin/windows_x86", "*.dll"))
  }
} else {
  if (Sys.info()["sysname"] == "Linux" & Sys.info()["machine"] == "x86_64") {
    files = Sys.glob(file.path(tmpdir, "bin/linux_x86_64", "*.so"))
  } else if (Sys.info()["sysname"] == "Darwin" & Sys.info()["machine"] == "x86_64") {
    files = Sys.glob(file.path(tmpdir, "bin/macos_universal", "*.dylib"))
  } else if (Sys.info()["sysname"] == "Darwin" & Sys.info()["machine"] == "arm64") {
    stop("ARM-based processors not supported.")
  } else {
    stop("OS not supported.")
  }
}

dir.create(dest, recursive = TRUE, showWarnings = FALSE)
rv = file.copy(files, dest, overwrite = TRUE)

# remove temporary files -------------------------------------------------------
if (!dir.exists(file.path("../..", TofDaqAPI))) {
  unlink(c(tmpzip, tmpdir), recursive = TRUE)
}
