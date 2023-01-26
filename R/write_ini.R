#' Write ini file
#'
#' @name write_ini
#' @rdname ini
#' @keywords internal
#' @export
#' @return Writes ini file to directory
NULL

write_ini <- function(BATCH,
                      OUT,
                      HABITAT,
                      LOCATION.FILE,
                      CONNECTION = 'connect_four_neighbors_only = FALSE',
                      CURRENT.MAP = "write_cur_maps = 0",
                      CUMONLY = 'write_cum_cur_map_only = True',
                      VOLTAGE = "write_volt_maps = False",
                      MAP = "write_cum_cur_map_only = False",
                      GRNDFILE,
                      VARSOURCE,
                      SOURCEFILE,
                      VARSOURCEFILE,
                      PARALLELIZE = FALSE,
                      CORES = NULL,
                      solver = 'cholmod',
                      precision = FALSE,
                      is_resistance = TRUE,
                      PAIRS_TO_INCLUDE = "included_pairs_file = (Browse for a file with pairs to include or exclude)",
                      POLYGON_FILE,
                      USE_POLYGONS = FALSE,
                      PAIRS = "use_included_pairs = False",
                      SCENARIO = "scenario = pairwise",
                      FOCALNODE = "set_focal_node_currents_to_zero = False",
                      MAXMAP = "write_max_cur_maps = False",
                      silent = NULL){
  if(PARALLELIZE == TRUE && !is.null(CORES)) {
    PARALLELIZE <- "parallelize = True"
    CORES <- paste0("max_parallel = ", CORES)
  } else {
    PARALLELIZE = "parallelize = False"
    CORES <- "max_parallel = 0"
  }

  if(!is.null(solver)) {
    solver <- "solver = cholmod"
  } else {
    solver <- "solver = cg+amg"
  }

  if(isTRUE(silent)) {
    LOG <- "log_level = critical"
  } else {
    LOG <- "log_level = INFO"
  }

  if(isTRUE(precision == 'single')) {
    precision <- "precision = single"
  } else {
    precision <- "precision = None"
  }

  if(isTRUE(is_resistance)) {
    RESISTANCE <- "connect_using_avg_resistances = False"
    HABITAT_RES <- "habitat_map_is_resistances = True"
  } else {
    RESISTANCE <- "connect_using_avg_resistances = False"
    HABITAT_RES <- "habitat_map_is_resistances = False"
  }

  sink(BATCH)
  cat("[Options for advanced mode]")
  cat("\n")
  cat("ground_file_is_resistances = True")
  cat("\n")
  cat("remove_src_or_gnd = rmvsrc")
  cat("\n")
  cat(GRNDFILE)
  cat("\n")
  cat("use_unit_currents = False")
  cat("\n")
  cat(SOURCEFILE)
  cat("\n")
  cat("use_direct_grounds = False")
  cat("\n")
  cat("\n")
  cat("[Mask file]")
  cat("\n")
  cat("mask_file = (Browse for a raster mask file)")
  cat("\n")
  cat("use_mask = False")
  cat("\n")
  cat("\n")
  cat("[Calculation options]")
  cat("\n")
  cat("low_memory_mode = False")
  cat("\n")
  cat(PARALLELIZE)
  cat("\n")
  cat(CORES)
  cat("\n")
  cat(solver)
  cat("\n")
  cat("print_timings = False")
  cat("\n")
  cat("preemptive_memory_release = False")
  cat("\n")
  cat("print_rusages = False")
  cat("\n")
  cat("\n")
  cat("[Short circuit regions (aka polygons)]")
  cat("\n")
  cat(POLYGON_FILE)
  cat("\n")
  cat(paste0('use_polygons = ', USE_POLYGONS))
  # cat("polygon_file = (Browse for a short-circuit region file)")
  # cat("\n")
  # cat("use_polygons = False")
  cat("\n")
  cat("\n")
  cat("[Options for one-to-all and all-to-one modes]")
  cat("\n")
  cat(VARSOURCE)
  cat("\n")
  cat(VARSOURCEFILE)
  cat("\n")
  cat("\n")
  cat("[Output options]")
  cat("\n")
  cat("set_null_currents_to_nodata = True")
  cat("\n")
  cat(FOCALNODE)
  cat("\n")
  cat("set_null_voltages_to_nodata = True")
  cat("\n")
  cat("compress_grids = False")
  cat("\n")
  cat(VOLTAGE)
  cat("\n")
  cat(CURRENT.MAP)
  cat("\n")
  cat(OUT)
  cat("\n")
  cat(MAP)
  cat("\n")
  cat("log_transform_maps = False")
  cat("\n")
  cat(MAXMAP)
  cat("\n")
  cat("\n")
  cat("[Options for reclassification of habitat data]")
  cat("\n")
  cat("reclass_file = (Browse for file with reclassification data)")
  cat("\n")
  cat("use_reclass_table = False")
  cat("\n")
  cat("\n")
  cat("[Logging Options]")
  cat("\n")
  cat(LOG)
  cat("\n")
  cat("log_file = None")
  cat("\n")
  cat("profiler_log_file = None")
  cat("\n")
  cat("screenprint_log = False")
  cat("\n")
  cat("\n")
  cat("[Options for pairwise and one-to-all and all-to-one modes]")
  cat("\n")
  cat(PAIRS_TO_INCLUDE)
  cat("\n")
  cat(PAIRS)
  cat("\n")
  cat(LOCATION.FILE)
  cat("\n")
  cat("\n")
  cat("[Connection scheme for raster habitat data]")
  cat("\n")
  cat(RESISTANCE)
  cat("\n")
  cat(CONNECTION)
  cat("\n")
  cat("\n")
  cat("[Habitat raster or graph]")
  cat("\n")
  cat(HABITAT_RES)
  cat("\n")
  cat(HABITAT)
  cat("\n")
  cat("\n")
  cat("[Circuitscape mode]")
  cat("\n")
  cat("data_type = raster")
  cat("\n")
  cat(SCENARIO)
  cat("\n")
  cat(precision)
  sink()
}
