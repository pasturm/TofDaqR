files <- Sys.glob(c(file.path(R_PACKAGE_SOURCE, paste0("src", R_ARCH), "*.dll"),
                    file.path(R_PACKAGE_SOURCE, "src", "*.dll")))
dest <- file.path(R_PACKAGE_DIR, paste0("libs", R_ARCH))
dir.create(dest, recursive = TRUE, showWarnings = FALSE)
file.copy(files, dest, overwrite = TRUE)
if(file.exists("TofDaqR.dll"))
  file.copy("TofDaqR.dll", dest, overwrite = TRUE)
