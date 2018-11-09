# copy the relevant dynamic libraries to ${R_PACKAGE_DIR}/libs

R_PACKAGE_DIR = commandArgs(TRUE)  # passed by Makevars
lib_path = "../tools/bin"

if (Sys.info()["sysname"] == "Windows") {
  dest = file.path(R_PACKAGE_DIR, paste0("libs", Sys.getenv("R_ARCH")))
  if (Sys.getenv("R_ARCH") == "/x64") {
    files = Sys.glob(file.path(lib_path, "windows_x64", "*.dll"))
  } else if (Sys.getenv("R_ARCH") == "/i386") {
    files = Sys.glob(file.path(lib_path, "windows_x86", "*.dll"))
  }
} else {
  dest = file.path(R_PACKAGE_DIR, "libs")
  if (Sys.info()["sysname"] == "Linux" & Sys.info()["machine"] == "x86_64") {
    files = Sys.glob(file.path(lib_path, "linux_x86_64", "*.so"))
  } else if (Sys.info()["sysname"] == "Darwin" & Sys.info()["machine"] == "x86_64") {
    files = Sys.glob(file.path(lib_path, "macos_universal", "*.dylib"))
  } else {
    print("OS not supported.")
  }
}

dir.create(dest, recursive = TRUE, showWarnings = FALSE)
file.copy(files, dest, overwrite = TRUE)
