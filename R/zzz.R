# .onLoad <-function(lib,pkg){
#   ver <- read.dcf(file.path(lib, pkg, "DESCRIPTION"),
#                   "Version")
#   ver <- as.character(ver)
#   library.dynam("TofDaqWrapR", pkg, lib)
#   cat("mypackage", ver,file.path(lib,pkg), "loaded\n")
# }

.onLoad <- function(lib, pkg) {
  library.dynam('TofDaqWrapR', pkg, lib)
}

.onUnload <- function(libpath) {
  CleanupDll()
  # library.dynam.unload('TofDaqWrapR', libpath)  # this often crashes R
  library.dynam.unload('TofDaqR', libpath)
}
