pkg_vec <- list(
  'dplyr',
  'ggnetwork',
  'ggplot2',
  'lwgeom',
  'magrittr',
  'qs',
  'readxl',
  'reshape2',
  'rlang',
  'rprojroot',
  'sf',
  'sfnetworks',
  'stringr',
  'tarchetypes',
  'targets',
  'terra',
  'tibble',
  'tidygraph',
  'tidyterra',
  'data.table')
lapply(pkg_vec, function(p) {
  suppressWarnings(
    suppressPackageStartupMessages(
      library(p, character.only = TRUE, warn.conflicts = FALSE, quietly = FALSE)
    ))
})