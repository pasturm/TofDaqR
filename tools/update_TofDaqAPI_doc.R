# Script to update the documentation files of the TofDaq API in tools/doc.
# This needs to be run manually once if the TofDaq API version is updated.

TofDaqAPI = "TofDaq_1.99r759_API_20180912"

# download TofDaq API files from https://soft.tofwerk.com/ and unzip -----------
tmpzip = tempfile()
tmpdir = file.path(tempdir(), TofDaqAPI)
dir.create(tmpdir)
download.file(paste0("https://soft.tofwerk.com/", TofDaqAPI, ".zip"), tmpzip)
unzip(tmpzip, exdir = tmpdir)

# copy the doc folder to tools/doc ---------------------------------------------

# first remove unnecessary jquery.php and js.php from doc/
unlink(file.path(tmpdir, list.files(tmpdir, "*.php", recursive = TRUE)))

# copy doc/ to /tools
# (assuming the working directory is set to the package source directory)
file.copy(file.path(tmpdir, "doc"), "./tools", recursive = TRUE)

# remove temporary files -------------------------------------------------------
unlink(c(tmpzip, tmpdir), recursive = TRUE)
