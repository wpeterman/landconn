# Circuitscape preparation helpers. All grids share the habitat geometry and use
# signed-capable double precision when serialized, including their nodata cells.
.cs_vector <- function(x, label) {
  if (!inherits(x, "SpatVector")) {
    x <- tryCatch(terra::vect(x), error = function(e) {
      stop("`", label, "` must contain spatial points or polygons.", call. = FALSE)
    })
  }
  if (!nrow(x) || !terra::geomtype(x) %in% c("points", "polygons")) {
    stop("`", label, "` must contain nonempty points or polygons.", call. = FALSE)
  }
  x
}

.cs_ids <- function(x, field) {
  if (is.null(field)) return(seq_len(nrow(x)))
  if (!is.character(field) || length(field) != 1L || !field %in% names(x)) {
    stop("`field` must name an existing numeric ID column.", call. = FALSE)
  }
  ids <- x[[field]][[1]]
  if (!is.numeric(ids) || any(!is.finite(ids) | ids <= 0 | ids != floor(ids)) ||
      any(ids > 2^53 - 1)) {
    stop("`field` IDs must be finite positive integers no greater than 2^53 - 1.", call. = FALSE)
  }
  ids
}

.cs_prepare <- function(r, input_locs, regions, field) {
  x <- .cs_vector(if (is.null(regions)) input_locs else regions,
                  if (is.null(regions)) "input_locs" else "regions")
  if (!is.null(regions) && terra::geomtype(x) != "polygons") {
    stop("`regions` must contain polygons.", call. = FALSE)
  }
  if (nzchar(terra::crs(r)) && nzchar(terra::crs(x)) && !terra::same.crs(r, x)) {
    stop("Locations and habitat must use the same coordinate reference system.", call. = FALSE)
  }
  ids <- .cs_ids(x, field)
  x$.landconn_id <- ids
  if (terra::geomtype(x) == "points") {
    cells <- terra::cellFromXY(r, terra::crds(x))
    if (anyNA(cells) || anyDuplicated(cells)) {
      stop("Each point must occupy a distinct cell inside the habitat extent.", call. = FALSE)
    }
  }
  grid <- terra::rasterize(x, r, field = ".landconn_id", background = NA)
  # A core outside the habitat or entirely inside its missing mask cannot be used.
  grid <- terra::mask(grid, r)
  vals <- terra::values(grid, mat = FALSE)
  if (!setequal(unique(ids), unique(vals[!is.na(vals)]))) {
    stop("Every ID must cover at least one valid habitat cell; check small, masked, or overlapping cores.", call. = FALSE)
  }
  poly <- if (terra::geomtype(x) == "polygons") grid else NULL
  if (!is.null(regions)) {
    # Sample raster cells rather than continuous points: each draw is guaranteed
    # to lie in its own represented region and in a valid habitat cell.
    focal <- rep(NA_real_, terra::ncell(r))
    for (id in unique(ids)) {
      available <- which(!is.na(vals) & vals == id)
      cell <- available[sample.int(length(available), 1L)]
      focal[cell] <- id
    }
    grid <- terra::setValues(terra::rast(r), focal)
  }
  list(focal = grid, regions = poly, ids = sort(unique(ids)))
}

.cs_table <- function(x, ids, name, ground = FALSE) {
  x <- as.matrix(x)
  if (!is.numeric(x) || ncol(x) != 2L || !nrow(x) || any(!is.finite(x)) ||
      anyDuplicated(x[, 1]) || any(!x[, 1] %in% ids) ||
      any(x[, 2] < 0) || (!ground && !any(x[, 2] > 0))) {
    stop("`", name, "` must be a two-column numeric table of unique existing IDs and nonnegative values",
         if (!ground) ", including a positive source strength." else ".", call. = FALSE)
  }
  x
}

.cs_advanced_grid <- function(table, focal, polygons) {
  vals <- terra::values(focal, mat = FALSE)
  out <- rep(NA_real_, length(vals))
  for (i in seq_len(nrow(table))) {
    cells <- which(!is.na(vals) & vals == table[i, 1])
    if (length(cells) > 1L && is.null(polygons)) {
      stop("Advanced mode needs one point cell per ID, or polygon cores collapsed into short-circuit regions.", call. = FALSE)
    }
    # One injection/ground per collapsed core preserves the requested total
    # current or ground resistance regardless of the number of polygon cells.
    out[cells[1]] <- table[i, 2]
  }
  terra::setValues(terra::rast(focal), out)
}

.cs_result <- function(out, scenario) {
  if (!is.matrix(out) || !is.numeric(out)) stop("Circuitscape returned an unexpected result.", call. = FALSE)
  if (scenario == "pairwise") {
    if (nrow(out) < 3L || ncol(out) != nrow(out)) {
      stop("Circuitscape returned an invalid pairwise matrix.", call. = FALSE)
    }
    ids <- as.character(out[-1, 1])
    ans <- out[-1, -1, drop = FALSE]
    dimnames(ans) <- list(ids, as.character(out[1, -1]))
    return(ans)
  }
  if (scenario %in% c("one-to-all", "all-to-one")) {
    if (ncol(out) != 2L) stop("Circuitscape returned an invalid per-node table.", call. = FALSE)
    colnames(out) <- c("node_id", "resistance")
  }
  # Advanced mode returns the complete voltage grid, without an ID header.
  out
}

.cs_compute <- function(config) {
  # Use GDAL's configuration API: updating Julia ENV alone is insufficient on
  # Windows when GDAL and Julia use different C runtimes. Native readers in newer
  # solver versions need no GDAL override.
  old <- juliaEval('isdefined(Circuitscape, :ArchGDAL) ? Circuitscape.ArchGDAL.getconfigoption("AAIGRID_DATATYPE") : nothing')
  if (!is.null(old)) {
    juliaEval('Circuitscape.ArchGDAL.setconfigoption("AAIGRID_DATATYPE", "Float64")')
    on.exit({
      restore <- if (!nzchar(old)) 'Circuitscape.ArchGDAL.clearconfigoption("AAIGRID_DATATYPE")' else
        paste0('Circuitscape.ArchGDAL.setconfigoption("AAIGRID_DATATYPE", ', encodeString(old, quote = '\"'), ')')
      tryCatch(juliaEval(restore), error = function(e) {
        warning("Could not restore Julia's ASCII reader setting: ", conditionMessage(e), call. = FALSE)
      })
    }, add = TRUE)
  }
  juliaCall("compute", config)
}
