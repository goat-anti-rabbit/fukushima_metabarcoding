# Project settings shared by all scripts.
# Run the scripts from R, with the project folder as working directory, e.g. source("src/01_input_reads.R")

set.seed(20261007)   # fixed seed for anything random

PATH <- list(
  data    = "data",
  db      = "db",                                                                     # reference databases, see db/README.md
  reads   = function(marker, type) file.path("data", paste0("reads_", marker), type), # type: delivered / input / ...
  objects = "results/objects",
  figures = "results/figures",
  tables  = "results/tables",
  logs    = "logs"
)
