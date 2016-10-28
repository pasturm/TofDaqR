.onLoad <- function(lib, pkg) {
  ver <- as.character(read.dcf(file.path(lib, pkg, "DESCRIPTION"), "Version"))
  packageStartupMessage("TofDaqR ", ver, " loaded.")
}

.onUnload <- function(libpath) {
  CleanupDll()
  library.dynam.unload("TofDaqR", libpath)
}
