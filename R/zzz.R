.onLoad <- function(lib, pkg) {
  ver <- as.character(read.dcf(file.path(lib, pkg, "DESCRIPTION"), "Version"))
  cat("TofDaqR", ver, "loaded.\n")
}

.onUnload <- function(libpath) {
  CleanupDll()
  library.dynam.unload("TofDaqWrapR", libpath)
}
