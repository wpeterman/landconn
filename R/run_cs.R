#' Run Circuitscape from R
#'
#' Prepare focal-node and habitat grids, run the Julia version of Circuitscape,
#' and return resistance values, voltages, or current maps.
#'
#' @param JULIA_HOME Character path to the Julia binary directory. An already
#'   running JuliaConnectoR session is reused; restart R to change Julia binaries.
#' @param r A single-layer `SpatRaster` or `RasterLayer` with square cells.
#'   Values must be finite and positive. Use `NA` for excluded habitat cells.
#' @param input_locs Spatial points or polygons as a `SpatVector`, `sf`, or
#'   supported `sp` object. Points must occupy distinct valid habitat cells.
#'   Polygon cells with one ID form a core. May be `NULL` when `regions` is supplied.
#' @param regions Optional polygons (`sf` or `SpatVector`) defining short-circuit
#'   regions. When supplied, these replace `input_locs`: one valid raster cell is
#'   sampled within each region ID for its focal node. Set a seed before calling
#'   for reproducible sampling. Default is `NULL`.
#' @param field Optional name of a numeric ID column. IDs must be positive integers
#'   no greater than `2^53 - 1`. Without a field, IDs are the feature row numbers.
#'   Supplied IDs are preserved. Invalid fields or missing cores produce an error.
#' @param CurrentMap Logical, whether to write current maps (default `TRUE`).
#'   Requesting `output = "raster"` enables this automatically.
#' @param cumulative_map_only Logical, default `TRUE`: write the cumulative map
#'   without individual pair maps. `FALSE` also reads individual current maps into
#'   the returned raster. Advanced mode has a single current map.
#' @param VoltageMaps Logical, write voltage maps and include them in raster
#'   output (default `FALSE`).
#' @param scenario One of `"pairwise"` (default), `"one-to-all"`, `"all-to-one"`,
#'   or `"advanced"`. Return structures differ by scenario; see Value.
#' @param Neighbor_Connect Either 8 (default, including diagonal neighbors) or 4
#'   (cardinal neighbors only).
#' @param output_dir Existing directory for solver results. Defaults to the R
#'   session temporary directory. Result maps are retained there.
#' @param output_name Single filename stem, without directory separators.
#'   Defaults to the habitat layer name. Use a new stem to keep runs separate.
#' @param output `"raster"` (default) for current/voltage maps, or `"matrix"`
#'   for the scenario's numeric solver result.
#' @param pairs_to_include Optional vector of 0/1 values for pairwise mode, one
#'   per pair in the column-wise lower triangle of a matrix ordered by increasing
#'   node ID. At least one pair must be retained. Other scenarios reject this option.
#' @param variable_source Optional two-column numeric table: node ID, current
#'   strength in amperes. Values must be nonnegative, with at least one positive
#'   source. In one-to-all/all-to-one, unlisted IDs retain Circuitscape's default
#'   unit strength. In advanced mode, only listed IDs inject current and a
#'   `variable_ground` table is required. Not used in pairwise mode.
#' @param variable_ground Advanced-mode two-column numeric table: node ID,
#'   resistance to ground in the same resistance units as the modeled circuit.
#'   Zero gives a direct ground. Source and ground IDs must be disjoint.
#' @param parallel Logical, use Julia threads (default `FALSE`). Threads must be
#'   configured before Julia starts; see Details.
#' @param cores Optional positive integer thread count, used only with
#'   `parallel = TRUE`. It sets `JULIA_NUM_THREADS` before a new Julia session and
#'   must equal the thread count of an existing session. A mismatch is an error.
#' @param cholmod Logical, use the direct CHOLMOD solver (`TRUE`, default) or
#'   the iterative `cg+amg` solver (`FALSE`).
#' @param precision Logical: `FALSE` (default) for double precision or `TRUE`
#'   for single precision. The strings `"double"` and `"single"` are also accepted.
#' @param is_resistance Logical, whether `r` represents resistance (`TRUE`,
#'   default) or conductance (`FALSE`).
#' @param focal_node_current_zero Logical, set active focal-node currents to
#'   zero in maps (default `FALSE`); resistance values are unaffected.
#' @param max_map Logical, also save the maximum current map (default `FALSE`).
#'   This extra file is not included in the returned raster.
#' @param remove_files Logical, remove temporary input grids/configuration on
#'   success and failure (default `TRUE`). Solver output maps are retained.
#' @param silent Logical, request suppressed solver messages (default `TRUE`).
#'   Errors retain the original Julia condition text.
#'
#' @return With `output = "matrix"`: pairwise mode returns a square resistance
#'   matrix with row and column names equal to node IDs. One-to-all/all-to-one
#'   return a two-column matrix (`node_id`, `resistance`), one row per focal ID.
#'   Advanced mode returns the complete voltage matrix in raster row/column order,
#'   without stripping any row or column. Missing/excluded cells follow the
#'   solver's nodata convention.
#'
#'   With `output = "raster"`: a `SpatRaster` for SpatRaster input, or a legacy
#'   Raster object for RasterLayer input. The first layer is cumulative current
#'   (the single current grid in advanced mode). Requested individual currents
#'   follow, then voltages. Layers carry the habitat CRS and file-based names.
#'   Current maps express modeled current under the specified sources and grounds;
#'   they are not calibrated movement probabilities.
#'
#' @details Focal and short-circuit grids use `-9999` nodata and double precision,
#'   preserving IDs and preventing empty cells from becoming positive unsigned
#'   values. The Julia ASCII reader is temporarily set to double precision to
#'   prevent large IDs from being merged during import, then its prior setting
#'   is restored. Every core must cover valid habitat. Coordinate systems must agree
#'   when both inputs declare one; undeclared coordinates are assumed aligned.
#'
#'   In advanced mode, one source or ground is attached to the first represented
#'   cell of each listed core. Polygon cells of that ID are already collapsed into
#'   one circuit node, so the requested source is a total current per core, not
#'   per cell. A point ID spanning several cells is ambiguous in advanced mode and
#'   is rejected. Ground values are resistances, including when `r` is conductance.
#'
#'   For parallel work, start a fresh R session and set `cores`, or set
#'   `JULIA_NUM_THREADS` before any Julia call. The actual thread count is checked;
#'   it cannot be changed in a running Julia process. The R environment variables
#'   are restored after the call, but the Julia process persists. Verified locally
#'   with Julia 1.12.6 and Circuitscape 5.17.1. Compatibility with other versions
#'   requires running the optional integration suite.
#'
#'   See `vignette("circuitscape-workflows", package = "landconn")` for polygon,
#'   source/ground, output and parallel examples, and
#'   <https://docs.circuitscape.org/Circuitscape.jl/latest/options/> for solver options.
#' @seealso [omni_cs()], [julia_packages()]
#' @examples
#' \dontrun{
#' library(terra)
#' r <- rast(nrows = 8, ncols = 8, xmin = 0, xmax = 8, ymin = 0, ymax = 8,
#'           crs = "EPSG:3857")
#' values(r) <- 1
#' pts <- vect(cbind(c(1.5, 6.5), c(1.5, 6.5)), crs = crs(r))
#' pts$core_id <- c(11L, 44L)
#' cores <- buffer(pts, width = 0.6)
#' distances <- run_cs(Sys.getenv("JULIA_BINDIR"), r, cores, field = "core_id",
#'                     CurrentMap = FALSE, output = "matrix")
#' distances # resistance matrix named by core ID
#' }
#' @export
#' @author Bill Peterman <peterman.73@@osu.edu>
#'
#' @importFrom stats dist
#' @importFrom raster reclassify rasterToPoints writeRaster raster coordinates stack crs rasterize
#' @importFrom sp coordinates
#' @importFrom terra geomtype rast isFALSE isTRUE crs crs<- rasterize writeRaster crop ext global extend vect
#' @importFrom methods as
#' @importFrom JuliaConnectoR juliaCall juliaImport juliaEval juliaSetupOk
#' @importFrom sf as_Spatial st_as_sf st_geometry_type st_sample
#' @importFrom utils write.table
#'
run_cs <- function(JULIA_HOME = NULL,
                   r,
                   input_locs = NULL,
                   regions = NULL,
                   field = NULL,
                   CurrentMap = TRUE,
                   cumulative_map_only = TRUE,
                   VoltageMaps = FALSE,
                   scenario = c('pairwise', 'advanced', 'one-to-all', 'all-to-one'),
                   Neighbor_Connect = 8,
                   output_dir = NULL,
                   output_name = NULL,
                   output = "raster",
                   pairs_to_include = NULL,
                   variable_source = NULL,
                   variable_ground = NULL,
                   parallel = FALSE,
                   cores = NULL,
                   cholmod = TRUE,
                   precision = FALSE,
                   is_resistance = TRUE,
                   focal_node_current_zero = FALSE,
                   max_map = FALSE,
                   remove_files = TRUE,
                   silent = TRUE) {

  scenario <- match.arg(scenario)
  output <- match.arg(output, c("raster", "matrix"))
  input_raster <- inherits(r, "Raster")
  if (input_raster) r <- terra::rast(r)
  if (!inherits(r, "SpatRaster") || terra::nlyr(r) != 1L) {
    stop("`r` must be a single-layer RasterLayer or SpatRaster.")
  }
  rv <- terra::values(r, mat = FALSE)
  if (!any(is.finite(rv)) || any(!is.na(rv) & (!is.finite(rv) | rv <= 0))) {
    stop("Habitat values must be finite and positive, or NA for excluded cells.")
  }
  if (!isTRUE(all.equal(terra::res(r)[1], terra::res(r)[2]))) {
    stop("Circuitscape ASCII grids require square raster cells.")
  }
  if (is.null(JULIA_HOME) || !is.character(JULIA_HOME) || length(JULIA_HOME) != 1L) {
    stop("`JULIA_HOME` must specify the Julia binary directory.")
  }
  if (is.null(input_locs) && is.null(regions)) stop("Supply `input_locs` or `regions`.")
  for (nm in c("CurrentMap", "cumulative_map_only", "VoltageMaps", "parallel", "cholmod",
               "is_resistance", "focal_node_current_zero", "max_map", "remove_files", "silent")) {
    val <- get(nm)
    if (!is.logical(val) || length(val) != 1L || is.na(val)) stop("`", nm, "` must be TRUE or FALSE.")
  }
  if (length(Neighbor_Connect) != 1L || is.na(Neighbor_Connect) || !Neighbor_Connect %in% c(4, 8)) {
    stop("`Neighbor_Connect` must be 4 or 8.")
  }
  precision <- if (identical(precision, TRUE)) "single" else if (identical(precision, FALSE)) "double" else
    match.arg(precision, c("double", "single"))
  prepared <- .cs_prepare(r, input_locs, regions, field)
  ids <- prepared$ids
  if (scenario != "advanced" && length(ids) < 2L) stop("At least two distinct focal IDs are required.")
  if (!is.null(variable_ground) && scenario != "advanced") stop("`variable_ground` is only used in advanced mode.")
  if (!is.null(variable_source) && scenario == "pairwise") stop("Variable source strengths are not used in pairwise mode.")
  if (scenario == "advanced" && (is.null(variable_source) || is.null(variable_ground))) {
    stop("Advanced mode requires both `variable_source` and `variable_ground`.")
  }
  if (!is.null(variable_source)) variable_source <- .cs_table(variable_source, ids, "variable_source")
  if (!is.null(variable_ground)) variable_ground <- .cs_table(variable_ground, ids, "variable_ground", ground = TRUE)
  if (scenario == "advanced" && any(variable_source[, 1] %in% variable_ground[, 1])) {
    stop("Advanced sources and grounds must have disjoint IDs.")
  }
  if (!is.null(pairs_to_include)) {
    if (scenario != "pairwise") stop("`pairs_to_include` is supported only in pairwise mode.")
    if (!is.numeric(pairs_to_include) && !is.logical(pairs_to_include)) stop("`pairs_to_include` must contain 0 or 1.")
    if (length(pairs_to_include) != choose(length(ids), 2) || anyNA(pairs_to_include) ||
        any(!pairs_to_include %in% c(0, 1)) || !any(pairs_to_include == 1)) {
      stop("`pairs_to_include` needs one 0/1 per lower-triangle pair, including at least one retained pair.")
    }
  }
  if (is.null(output_dir)) output_dir <- tempdir()
  if (!dir.exists(output_dir)) stop("`output_dir` must be an existing directory.")
  output_dir <- normalizePath(output_dir, winslash = "/", mustWork = TRUE)
  if (is.null(output_name)) output_name <- names(r)[1]
  if (length(output_name) != 1L || is.na(output_name) || !nzchar(output_name) ||
      grepl("[/\\\\\r\n]", output_name)) stop("`output_name` must be a single filename stem.")
  if (output == "raster") CurrentMap <- TRUE

  # Keep input cleanup separate from user-requested result maps. A failure at any
  # later stage still releases all temporary input grids and their sidecars.
  work <- tempfile("landconn_cs_")
  dir.create(work)
  if (remove_files) on.exit(unlink(work, recursive = TRUE), add = TRUE)
  else message("Circuitscape input files retained in ", normalizePath(work, winslash = "/"))
  files <- file.path(work, c("habitat.asc", "focal.asc", "regions.asc", "source.asc", "ground.asc", "strengths.txt", "pairs.txt", "run.ini"))
  names(files) <- c("habitat", "focal", "regions", "source", "ground", "strengths", "pairs", "ini")
  files <- chartr("\\", "/", files)
  write_grid <- function(x, path) terra::writeRaster(x, path, datatype = "FLT8S", NAflag = -9999, overwrite = TRUE)
  write_grid(r, files["habitat"])
  write_grid(prepared$focal, files["focal"])
  use_poly <- !is.null(prepared$regions)
  if (use_poly) write_grid(prepared$regions, files["regions"])
  source_line <- "source_file = None"
  ground_line <- "ground_file = None"
  strength_line <- "variable_source_file = None"
  if (scenario == "advanced") {
    write_grid(.cs_advanced_grid(variable_source, prepared$focal, prepared$regions), files["source"])
    write_grid(.cs_advanced_grid(variable_ground, prepared$focal, prepared$regions), files["ground"])
    source_line <- paste("source_file =", files["source"])
    ground_line <- paste("ground_file =", files["ground"])
  } else if (!is.null(variable_source)) {
    utils::write.table(variable_source, files["strengths"], row.names = FALSE, col.names = FALSE, quote = FALSE)
    strength_line <- paste("variable_source_file =", files["strengths"])
  }
  pair_line <- "included_pairs_file = None"
  if (!is.null(pairs_to_include)) {
    pair_indices <- which(lower.tri(matrix(0, length(ids), length(ids))), arr.ind = TRUE)
    pairs <- cbind(ids[pair_indices[, 1]], ids[pair_indices[, 2]])
    writeLines("mode include", files["pairs"])
    utils::write.table(pairs[pairs_to_include == 1, , drop = FALSE], files["pairs"],
                       append = TRUE, row.names = FALSE, col.names = FALSE, quote = FALSE)
    pair_line <- paste("included_pairs_file =", files["pairs"])
  }

  # Julia threads are fixed when the Julia process starts. Set its startup
  # environment before launching, and restore the R caller's environment afterward.
  old_env <- Sys.getenv(c("JULIA_BINDIR", "JULIA_NUM_THREADS"), unset = NA_character_)
  on.exit({
    for (nm in names(old_env)) {
      if (is.na(old_env[[nm]])) Sys.unsetenv(nm) else do.call(Sys.setenv, stats::setNames(list(old_env[[nm]]), nm))
    }
  }, add = TRUE)
  if (!is.null(cores)) {
    if (!parallel) stop("Supply `cores` only with `parallel = TRUE`.")
    if (!is.numeric(cores) || length(cores) != 1L || !is.finite(cores) || cores < 1 || cores != floor(cores)) {
      stop("`cores` must be a positive integer.")
    }
    Sys.setenv(JULIA_NUM_THREADS = as.character(cores))
  }
  if (!juliaSetupOk()) Sys.setenv(JULIA_BINDIR = JULIA_HOME)
  if (!juliaSetupOk()) stop("Check the Julia binary directory in `JULIA_HOME`.")
  result <- tryCatch({
    juliaEval("using Circuitscape")
    if (parallel) {
      threads <- as.integer(juliaEval("Threads.nthreads()"))
      if (!is.null(cores) && threads != cores) stop("The running Julia session has ", threads,
          " threads, not ", cores, ". Restart R and set JULIA_NUM_THREADS before starting Julia.")
      if (threads < 2L) stop("Parallel execution requires Julia to start with at least two threads.")
    }
    write_ini(BATCH = files["ini"], OUT = paste("output_file =", file.path(output_dir, output_name)),
              HABITAT = paste("habitat_file =", files["habitat"]),
              LOCATION.FILE = paste("point_file =", files["focal"]),
              CONNECTION = paste("connect_four_neighbors_only =", Neighbor_Connect == 4),
              GRNDFILE = ground_line, VARSOURCE = paste("use_variable_source_strengths =", scenario != "advanced" && !is.null(variable_source)),
              SOURCEFILE = source_line, VARSOURCEFILE = strength_line,
              MAP = paste("write_cum_cur_map_only =", cumulative_map_only),
              VOLTAGE = paste("write_volt_maps =", VoltageMaps), CURRENT.MAP = paste("write_cur_maps =", CurrentMap),
              PAIRS_TO_INCLUDE = pair_line, PAIRS = paste("use_included_pairs =", !is.null(pairs_to_include)),
              PARALLELIZE = parallel, solver = if (cholmod) "cholmod" else NULL, precision = precision,
              POLYGON_FILE = paste("polygon_file =", if (use_poly) files["regions"] else "None"),
              USE_POLYGONS = use_poly, SCENARIO = paste("scenario =", scenario),
              FOCALNODE = paste("set_focal_node_currents_to_zero =", focal_node_current_zero),
              MAXMAP = paste("write_max_cur_maps =", max_map), silent = silent, is_resistance = is_resistance)
    .cs_result(.cs_compute(files["ini"]), scenario)
  }, error = function(e) stop("Circuitscape failed: ", conditionMessage(e), call. = FALSE))
  if (output == "matrix") return(result)

  # Select only ASCII grids, never projection sidecars. Prefix matching uses
  # literal stems so punctuation in an output name cannot act as a regex.
  main_file <- file.path(output_dir, paste0(output_name, if (scenario == "advanced") "_curmap.asc" else "_cum_curmap.asc"))
  if (!file.exists(main_file)) stop("Circuitscape did not write the requested cumulative current map.")
  paths <- main_file
  available <- list.files(output_dir, full.names = TRUE)
  stems <- basename(available)
  ascii <- endsWith(tolower(stems), ".asc")
  if (!cumulative_map_only) paths <- c(paths, available[ascii & startsWith(stems, paste0(output_name, "_curmap_"))])
  if (VoltageMaps) paths <- c(paths, available[ascii & (startsWith(stems, paste0(output_name, "_voltmap_")) |
                                                      stems == paste0(output_name, "_voltmap.asc"))])
  ans <- terra::rast(paths)
  names(ans) <- sub("\\.asc$", "", basename(paths))
  names(ans)[1] <- output_name
  terra::crs(ans) <- terra::crs(r)
  if (input_raster) return(methods::as(ans, "Raster"))
  ans
}
