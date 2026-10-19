# Inspect the actual ASCII grids at the Julia boundary without requiring Julia.
cs_grid_fixture <- function() {
  r <- terra::rast(nrows = 6, ncols = 6, xmin = 0, xmax = 6,
                   ymin = 0, ymax = 6, crs = "EPSG:3857")
  terra::values(r) <- 1
  patches <- terra::as.polygons(terra::rast(nrows = 2, ncols = 2,
                                           xmin = 1, xmax = 5,
                                           ymin = 1, ymax = 5,
                                           crs = terra::crs(r)))
  patches$id <- c(11L, 22L, 33L, 44L)
  list(r = r, patches = patches)
}

test_that("run_cs writes nodata outside focal nodes and short-circuit regions", {
  fixture <- cs_grid_fixture()
  ids <- fixture$patches$id
  seen <- NULL
  local_mocked_bindings(
    juliaSetupOk = function() TRUE,
    juliaEval = function(...) NULL,
    juliaCall = function(fun, config) {
      ini <- readLines(config, warn = FALSE)
      read_grid <- function(key) {
        path <- sub(paste0("^", key, " = "), "",
                    grep(paste0("^", key, " = "), ini, value = TRUE))
        lines <- readLines(path)
        # Read raw tokens: a GIS reader can hide an incorrect nodata encoding.
        nodata <- as.numeric(strsplit(trimws(lines[6]), "[[:space:]]+")[[1]][2])
        cells <- as.numeric(strsplit(paste(trimws(lines[-seq_len(6)]),
                                           collapse = " "), "[[:space:]]+")[[1]])
        list(nodata = nodata, cells = cells)
      }
      seen <<- list(focal = read_grid("point_file"))
      if ("use_polygons = True" %in% ini) {
        seen$regions <<- read_grid("polygon_file")
      }
      # Circuitscape's matrix includes a leading row and column of node IDs.
      distances <- matrix(1, length(ids), length(ids))
      diag(distances) <- 0
      rbind(c(0, ids), cbind(ids, distances))
    }
  )

  for (format in c("SpatVector", "sf", "automatic", "large_ids", "points")) {
    patches <- fixture$patches
    ids <- patches$id
    field <- "id"
    if (format == "large_ids") {
      # Double precision avoids narrowing valid IDs to signed 32-bit integers.
      patches$id <- ids <- c(2147483648, 2147483650, 2147483652, 2147483654)
    }
    locs <- patches
    if (format == "sf") locs <- sf::st_as_sf(patches)
    if (format == "automatic") {
      field <- NULL
      ids <- seq_len(nrow(patches))
    }
    if (format == "points") {
      locs <- terra::centroids(patches)
      ids <- seq_len(nrow(patches))
      field <- NULL
    }
    seen <- NULL
    result <- run_cs(JULIA_HOME = "mock", r = fixture$r, input_locs = locs,
                     field = field, CurrentMap = FALSE,
                     output = "matrix")
    expect_equal(unname(result), 1 - diag(4))
    for (grid in seen) {
      expect_equal(grid$nodata, -9999, info = format)
      expect_equal(sort(unique(grid$cells[grid$cells != grid$nodata])),
                   sort(ids), info = format)
      expect_true(any(grid$cells == grid$nodata), info = format)
    }
    seen <- lapply(seen, function(grid) {
      replace(grid$cells, grid$cells == grid$nodata, NA_real_)
    })
    if (format %in% c("SpatVector", "sf", "automatic", "large_ids")) {
      if (inherits(locs, "sf")) locs <- terra::vect(locs)
      expected <- terra::values(terra::rasterize(locs, fixture$r,
                                                field = if (is.null(field)) seq_len(4) else field),
                                 mat = FALSE)
      expect_equal(seen$focal, expected)
      expect_equal(seen$regions, expected)
    } else {
      expect_equal(sum(!is.na(seen$focal)), 4)
      expect_null(seen$regions)
    }
  }
})
