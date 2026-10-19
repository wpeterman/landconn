#' Write ini file
#'
#' Internal helper that writes a Circuitscape `.ini` configuration file. Called by
#' [run_cs()]; not intended to be called directly.
#'
#' @name write_ini
#' @rdname ini
#' @keywords internal
#' @return Writes an `.ini` file to disk and returns `NULL` invisibly.
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
  # Write directly to the file; never change the caller's output sink.
  lines <- c("[Options for advanced mode]", "ground_file_is_resistances = True",
             "remove_src_or_gnd = rmvsrc", GRNDFILE, "use_unit_currents = False", SOURCEFILE,
             "use_direct_grounds = False", "", "[Calculation options]",
             paste("parallelize =", PARALLELIZE),
             paste("solver =", if (is.null(solver)) "cg+amg" else "cholmod"),
             "", "[Short circuit regions (aka polygons)]", POLYGON_FILE,
             paste("use_polygons =", USE_POLYGONS), "",
             "[Options for one-to-all and all-to-one modes]", VARSOURCE, VARSOURCEFILE,
             "", "[Output options]", "set_null_currents_to_nodata = True", FOCALNODE,
             "set_null_voltages_to_nodata = True", "compress_grids = False", VOLTAGE,
             CURRENT.MAP, OUT, MAP, "log_transform_maps = False", MAXMAP,
             "", "[Logging Options]", paste("log_level =", if (isTRUE(silent)) "critical" else "INFO"),
             paste("suppress_messages =", isTRUE(silent)), "log_file = None",
             "", "[Options for pairwise and one-to-all and all-to-one modes]",
             PAIRS_TO_INCLUDE, PAIRS, LOCATION.FILE,
             "", "[Connection scheme for raster habitat data]",
             "connect_using_avg_resistances = False", CONNECTION,
             "", "[Habitat raster or graph]", paste("habitat_map_is_resistances =", is_resistance), HABITAT,
             "", "[Circuitscape mode]", "data_type = raster", SCENARIO,
             paste("precision =", if (isTRUE(precision) || identical(precision, "single")) "single" else "double"),
             "use_64bit_indexing = True")
  lines <- sub(" = TRUE$", " = True", lines)
  lines <- sub(" = FALSE$", " = False", lines)
  writeLines(lines, BATCH)
  invisible(NULL)
}
