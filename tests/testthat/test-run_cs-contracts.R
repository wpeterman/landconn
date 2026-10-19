cs_contract_fixture <- function() {
  r <- terra::rast(nrows = 8, ncols = 8, xmin = 0, xmax = 8,
                   ymin = 0, ymax = 8, crs = "EPSG:3857")
  terra::values(r) <- 1
  p <- terra::as.polygons(terra::rast(nrows = 2, ncols = 2, xmin = 2, xmax = 6,
                                     ymin = 2, ymax = 6, crs = terra::crs(r)))
  p$id <- c(44L, 11L, 33L, 22L)
  list(r = r, p = p, points = terra::centroids(p))
}

test_that("all point classes preserve unsorted nonconsecutive IDs", {
  f <- cs_contract_fixture()
  for (x in list(f$points, sf::st_as_sf(f$points), methods::as(f$points, "Spatial"))) {
    z <- .cs_prepare(f$r, x, NULL, "id")
    expect_equal(z$ids, c(11, 22, 33, 44))
    values <- terra::values(z$focal, mat = FALSE)
    expect_equal(values[terra::cellFromXY(f$r, terra::crds(f$points))], f$p$id)
  }
})

test_that("region sampling keeps IDs and avoids nodata cells", {
  f <- cs_contract_fixture()
  f$r[19] <- NA
  for (x in list(f$p, sf::st_as_sf(f$p))) {
    z <- .cs_prepare(f$r, NULL, x, "id")
    values <- terra::values(z$focal, mat = FALSE)
    expect_equal(sort(values[!is.na(values)]), c(11, 22, 33, 44))
    expect_equal(sum(!is.na(values)), 4)
    expect_true(is.na(values[19]))
  }
})

test_that("invalid spatial IDs and collisions fail before Julia", {
  f <- cs_contract_fixture()
  expect_error(.cs_prepare(f$r, f$points, NULL, "missing"), "existing")
  f$points$id[1] <- 0
  expect_error(.cs_prepare(f$r, f$points, NULL, "id"), "positive integers")
  dup <- f$points[c(1, 1), ]
  expect_error(.cs_prepare(f$r, dup, NULL, NULL), "distinct cell")
  outside <- terra::shift(f$points, dx = 100)
  expect_error(.cs_prepare(f$r, outside, NULL, NULL), "inside")
})

test_that("solver results retain all scenario-specific rows and IDs", {
  m <- rbind(c(0, 11, 44), c(11, 0, 3), c(44, 3, 0))
  expect_equal(.cs_result(m, "pairwise"), matrix(c(0, 3, 3, 0), 2,
                                               dimnames = list(c("11", "44"), c("11", "44"))))
  nodes <- cbind(c(11, 22, 33, 44), 1:4)
  for (mode in c("one-to-all", "all-to-one")) {
    z <- .cs_result(nodes, mode)
    expect_equal(unname(z), nodes)
    expect_named(as.data.frame(z), c("node_id", "resistance"))
  }
  expect_identical(.cs_result(matrix(1:6 + 0, 2), "advanced"), matrix(1:6 + 0, 2))
})

test_that("sources and grounds preserve core totals and validate tables", {
  f <- cs_contract_fixture()
  p <- .cs_prepare(f$r, f$p, NULL, "id")
  t <- .cs_table(cbind(c(11, 22), c(2, 3)), p$ids, "variable_source")
  x <- .cs_advanced_grid(t, p$focal, p$regions)
  expect_equal(sum(terra::values(x), na.rm = TRUE), 5)
  expect_equal(sum(!is.na(terra::values(x))), 2)
  expect_error(.cs_table(cbind(99, 1), p$ids, "variable_source"), "existing IDs")
  expect_error(.cs_table(matrix(1:6, 2), p$ids, "variable_ground", TRUE), "two-column")
})

test_that("run_cs honors option files and leaves caller sinks and environment intact", {
  f <- cs_contract_fixture()
  captured <- NULL
  local_mocked_bindings(
    juliaSetupOk = function() TRUE,
    juliaEval = function(code) if (code == "Threads.nthreads()") 2L else NULL,
    juliaCall = function(fun, path) {
      ini <- readLines(path)
      value <- function(key) sub(paste0("^", key, " = "), "", grep(paste0("^", key, " = "), ini, value = TRUE))
      captured <<- list(ini = ini, path = path, pairs = readLines(value("included_pairs_file")))
      rbind(c(0, 11, 22, 33, 44), cbind(c(11, 22, 33, 44), 1 - diag(4)))
    })
  old <- Sys.getenv(c("JULIA_BINDIR", "JULIA_NUM_THREADS"), unset = NA_character_)
  log <- tempfile()
  sink(log)
  depth <- sink.number()
  tryCatch({
    z <- run_cs("mock", f$r, f$points, field = "id", output = "matrix",
                CurrentMap = FALSE, precision = TRUE, pairs_to_include = c(1, 0, 0, 0, 0, 0),
                parallel = TRUE, cores = 2)
    expect_equal(sink.number(), depth)
  }, finally = sink())
  expect_identical(Sys.getenv(names(old), unset = NA_character_), old)
  expect_true("precision = single" %in% captured$ini)
  expect_true("use_included_pairs = True" %in% captured$ini)
  expect_equal(captured$pairs[1], "mode include")
  expect_match(captured$pairs[2], "22.*11")
  expect_false(file.exists(captured$path))
  expect_equal(rownames(z), as.character(c(11, 22, 33, 44)))
})

test_that("advanced configuration uses separate source and ground grids", {
  f <- cs_contract_fixture()
  captured <- NULL
  local_mocked_bindings(juliaSetupOk = function() TRUE, juliaEval = function(...) NULL,
    juliaCall = function(fun, path) {
      captured <<- readLines(path)
      matrix(0, 8, 8)
    })
  z <- run_cs("mock", f$r, f$p, field = "id", scenario = "advanced", output = "matrix",
               CurrentMap = FALSE, variable_source = cbind(11, 2), variable_ground = cbind(44, 0))
  expect_equal(dim(z), c(8, 8))
  expect_match(grep("^source_file", captured, value = TRUE), "source.asc")
  expect_match(grep("^ground_file =", captured, value = TRUE), "ground.asc")
  expect_true("use_variable_source_strengths = False" %in% captured)
})

test_that("original solver errors survive and temporary inputs are removed", {
  f <- cs_contract_fixture()
  path <- NULL
  local_mocked_bindings(juliaSetupOk = function() TRUE, juliaEval = function(...) NULL,
    juliaCall = function(fun, config) { path <<- config; stop("solver sentinel") })
  expect_error(run_cs("mock", f$r, f$points, output = "matrix"), "solver sentinel")
  expect_false(file.exists(path))
})

test_that("the Julia ASCII reader setting is restored even after solver failure", {
  setting <- "Float32"
  local_mocked_bindings(
    juliaEval = function(code) {
      if (grepl("getconfigoption", code, fixed = TRUE)) return(setting)
      setting <<- if (grepl("Float64", code, fixed = TRUE)) "Float64" else "Float32"
      invisible(NULL)
    },
    juliaCall = function(...) {
      expect_identical(setting, "Float64")
      stop("reader failure sentinel")
    })
  expect_error(.cs_compute("config.ini"), "reader failure sentinel")
  expect_identical(setting, "Float32")
})

test_that("additional output loading ignores projection sidecars", {
  f <- cs_contract_fixture()
  out <- tempfile(); dir.create(out)
  on.exit(unlink(out, recursive = TRUE))
  local_mocked_bindings(juliaSetupOk = function() TRUE, juliaEval = function(...) NULL,
    juliaCall = function(fun, path) {
      for (suffix in c("_cum_curmap", "_curmap_11_22", "_voltmap_11_22")) {
        terra::writeRaster(f$r, file.path(out, paste0("test.name", suffix, ".asc")), overwrite = TRUE)
      }
      rbind(c(0, 11, 22, 33, 44), cbind(c(11, 22, 33, 44), 1 - diag(4)))
    })
  z <- run_cs("mock", f$r, f$points, output_dir = out, output_name = "test.name",
               cumulative_map_only = FALSE, VoltageMaps = TRUE)
  expect_s4_class(z, "SpatRaster")
  expect_equal(terra::nlyr(z), 3)
  expect_true(terra::same.crs(z, f$r))
  legacy <- run_cs("mock", raster::raster(f$r), f$points, output_dir = out,
                    output_name = "test.name", cumulative_map_only = FALSE, VoltageMaps = TRUE)
  expect_true(inherits(legacy, "Raster"))
  expect_equal(raster::nlayers(legacy), 3)
})

test_that("Julia setup installs only Circuitscape and restores its environment", {
  calls <- character()
  local_mocked_bindings(juliaEval = function(code, ...) { calls <<- c(calls, code); NULL },
                        .package = "JuliaConnectoR")
  old <- Sys.getenv("JULIA_BINDIR", unset = NA_character_)
  julia_packages(tempdir(), install_libraries = TRUE)
  expect_identical(Sys.getenv("JULIA_BINDIR", unset = NA_character_), old)
  expect_equal(length(calls), 2)
  expect_match(calls[1], 'Pkg.add("Circuitscape")', fixed = TRUE)
  expect_identical(calls[2], "using Circuitscape")
  expect_false(any(grepl("ConScape|Plots", calls)))
})

test_that("omni_cs preserves holes and creates exactly the requested unique nodes", {
  f <- cs_contract_fixture()
  f$r[20] <- NA
  local_mocked_bindings(run_cs = function(r, input_locs, ...) {
    expect_true(is.na(terra::values(terra::crop(r, f$r))[20]))
    expect_equal(nrow(input_locs), 19)
    expect_equal(length(unique(terra::cellFromXY(r, terra::crds(input_locs)))), 19)
    r
  })
  z <- omni_cs("mock", f$r, n_nodes = 19)
  expect_true(is.na(terra::values(z)[20]))
  expect_equal(as.vector(terra::ext(z)), as.vector(terra::ext(f$r)))
})
