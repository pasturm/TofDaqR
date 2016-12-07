# Tps1rc ----------------------------------------------------------------
#' TPS1 RC codes and names.
#'
#' Data frame of TPS1 codes (as used in tps1rc.cfg to communicate with the TPS
#' through the API) and corresponding names.
#' @export
Tps1rc <- c(
  1, "RBP",
  2, "RG",
  3, "LENS",
  4, "HM",
  5, "PA",
  6, "MCP",
  9, "DRIFT",
  10, "U-low",
  11, "U-high",
  12, "U+low",
  13, "U+high",
  14, "L2",
  15, "DEFL",
  16, "DEFLFL",
  17, "IONEX",
  18, "L1",
  19, "EP",
  20, "IONCH",
  21, "FIL",
  22, "FILEM",
  23, "IFIL",
  24, "ILIM",
  27, "FILNUM",
  28, "SKIMMER",
  49, "Q1F",
  50, "Q1B",
  51, "Q2F",
  52, "Q2B",
  53, "LSK",
  55, "IMSLENS",
  56, "NOZZLE",
  57, "Q1EP",
  58, "IMSEND",
  62, "RF1",
  63, "RFAMP1",
  64, "RF2",
  65, "RFAMP2",
  66, "LFFREQ1",
  67, "LFAMP1",
  68, "LFFREQ2",
  69, "LFAMP2",
  70, "LFFREQ3",
  71, "LFAMP3",
  72, "LFFREQ4",
  73, "LFAMP4",
  83, "Grid+",
  84, "Grid-",
  85, "NEEDLE",
  86, "IMSHV",
  87, "IMSGRIDSTATE",
  88, "IMS_TR_FREQ",
  90, "IMSGateHV",
  116, "SKIM2",
  117, "REFERENCE",
  121, "CORONA",
  124, "IMR",
  162, "QUADA_COIL_SET",
  163, "QUADA_RFDAC_SET",
  164, "QUADB_COIL_SET",
  165, "QUADB_RFDAC_SET",
  200, "TOF_EXTR2",
  201, "TOF_EXTR1",
  202, "TOF_REF",
  203, "TOF_PULSE",
  600, "IONMODE",
  601, "INTERLOCK",
  602, "HVSUPPLY",
  603, "HVPOS",
  604, "HVNEG",
  605, "SPARE1",
  606, "SPARE2",
  607, "SPARE3",
  608, "MODE",
  609, "EXT24V",
  1016,	"TPSBODYTEMP",
  1017,	"TPSAIRTEMP",
  1018,	"HVATEMP",
  1019, "HVATEMP"
)
Tps1rc <- as.data.frame(matrix(Tps1rc, ncol = 2, byrow = TRUE), stringsAsFactors = FALSE)
names(Tps1rc) <- c("codes", "names")
Tps1rc$codes <-  as.integer(Tps1rc$codes)

# SaveMassTableToFile ----------------------------------------------------------
#' Saves the current peak parameters to a file.
#'
#' \code{SaveMassTableToFile} saves the current peak parameters to a file.
#'
#' Note: this function is not part of the TofDaq API, but is included in the
#' package for convenience.
#'
#' @param filename Path/filename. If not specified, "TmpMassTable.txt" in the
#' current working directory will be used.
#' @export
SaveMassTableToFile <- function(filename = "TmpMassTable.txt") {
  if (file.exists(filename)) {
    file.remove(filename)
    file.create(filename)
  }
  desc <- .Call("GetDescriptor")
  write.table(desc$NbrPeaks, file = filename, quote = FALSE, row.names = FALSE,
              col.names = FALSE)
  for (i in 1:desc$NbrPeaks) {
    rv <- .Call("GetPeakParameters", i-1)
    rv$TwRetVal <- NULL
    write.table(t(unlist(rv)), file = filename, append = TRUE,
                quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)
  }
}
