test_that("Circuitscape scenarios, maps, and parallel solves work with real Julia", {
  skip_if(Sys.getenv("LANDCONN_RUN_JULIA") != "true", "Set LANDCONN_RUN_JULIA=true for optional Julia integration.")
  home <- Sys.getenv("JULIA_BINDIR")
  skip_if(!nzchar(home), "Set JULIA_BINDIR to the Julia bin directory.")
  r <- terra::rast(nrows = 8, ncols = 8, xmin = 0, xmax = 8, ymin = 0, ymax = 8,
                   crs = "EPSG:3857")
  terra::values(r) <- 1
  pts <- terra::vect(cbind(c(1.5, 6.5, 6.5, 1.5), c(1.5, 1.5, 6.5, 6.5)), crs = terra::crs(r))
  pts$id <- c(44L, 11L, 33L, 22L)
  polygons <- terra::buffer(pts, width = .6)
  out <- tempfile(); dir.create(out)
  on.exit(unlink(out, recursive = TRUE))
  reader_before <- JuliaConnectoR::juliaEval('using Circuitscape; Circuitscape.ArchGDAL.getconfigoption("AAIGRID_DATATYPE")')
  run <- function(...) run_cs(JULIA_HOME = home, r = r, input_locs = polygons,
                               field = "id", output_dir = out, ...)
  serial <- run(output = "matrix", CurrentMap = FALSE, output_name = "serial")
  expect_equal(dim(serial), c(4, 4))
  expect_equal(rownames(serial), c("11", "22", "33", "44"))
  expect_equal(diag(serial), stats::setNames(rep(0, 4), rownames(serial)))
  parallel <- run(output = "matrix", CurrentMap = FALSE, output_name = "parallel", parallel = TRUE, cores = 2)
  expect_equal(parallel, serial, tolerance = 1e-8)
  expect_equal(as.integer(JuliaConnectoR::juliaEval("Threads.nthreads()")), 2L)
  selected <- run(output = "matrix", CurrentMap = FALSE, output_name = "selected",
                    pairs_to_include = c(1, 0, 0, 0, 0, 0))
  expect_equal(selected, serial[c("11", "22"), c("11", "22")], tolerance = 1e-8)
  maps <- run(output_name = "maps", cumulative_map_only = FALSE, VoltageMaps = TRUE)
  expect_s4_class(maps, "SpatRaster")
  expect_equal(terra::nlyr(maps), 13)
  expect_true(terra::same.crs(maps, r))
  for (scenario in c("one-to-all", "all-to-one")) {
    z <- run(scenario = scenario, output = "matrix", CurrentMap = FALSE,
               output_name = scenario, variable_source = cbind(c(11, 22, 33, 44), 1))
    expect_equal(dim(z), c(4, 2))
    expect_equal(z[, "node_id"], c(11, 22, 33, 44))
  }
  v <- run(scenario = "advanced", output = "matrix", CurrentMap = FALSE,
             output_name = "advanced", variable_source = cbind(11, 1), variable_ground = cbind(44, 0))
  expect_equal(dim(v), c(8, 8))
  expect_true(all(is.finite(v)))
  current <- run(scenario = "advanced", output_name = "advanced_map", VoltageMaps = TRUE,
                   variable_source = cbind(11, 1), variable_ground = cbind(44, 0))
  expect_equal(terra::nlyr(current), 2)
  large <- polygons
  large$id <- as.numeric(large$id) + 2147483648
  large_result <- run_cs(home, r, large, field = "id", CurrentMap = FALSE,
                          output = "matrix", output_dir = out, output_name = "large_ids")
  expect_equal(rownames(large_result), as.character(c(11, 22, 33, 44) + 2147483648))
  expect_equal(unname(large_result), unname(serial), tolerance = 1e-8)
  expect_identical(JuliaConnectoR::juliaEval('using Circuitscape; Circuitscape.ArchGDAL.getconfigoption("AAIGRID_DATATYPE")'), reader_before)
  r[28] <- NA
  omni <- omni_cs(home, r, n_nodes = 12, output_dir = out, output_name = "omni")
  expect_true(is.na(terra::values(omni)[28]))
})
