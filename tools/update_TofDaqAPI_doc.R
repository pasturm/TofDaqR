# Script to update the documentation files of the TofDaq API in tools/doc.
# This needs to be run manually once if the TofDaq API version is updated.

path = "./TofDaqAPI/"
TofDaqAPI = "TofDaq_1.99r1369_API_20210930"

# unzip ------------------------------------------------------------------------
tmpdir = tempdir()
if (!dir.exists(tmpdir)) { dir.create(tmpdir) }
unzip(paste0(path, TofDaqAPI, ".zip"), exdir = tmpdir)

# remove unnecessary jquery.php, js.php and prompt.js from doc/ ----------------
unlink(file.path(tmpdir, list.files(tmpdir, "*.php", recursive = TRUE)))
unlink(file.path(tmpdir, list.files(tmpdir, "*.js", recursive = TRUE)))

# copy doc/ to /tools ----------------------------------------------------------
# (assuming the working directory is set to the package source directory)
file.copy(file.path(tmpdir, "doc"), "./tools", recursive = TRUE)

# remove temporary files -------------------------------------------------------
unlink(file.path(tmpdir, TofDaqAPI), recursive = TRUE)
