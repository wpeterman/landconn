############ RUN JULIA ############

#' Prepare data and run Julia
#'
#' Create an omnidirectional current flow map
#'
#' @param JULIA_HOME Path to the folder containing the Julia binary (See Details).
#' @param r Accepts a RasterLayer object or SpatRaster.
#' @param n_nodes How many points to surround landscape. (Default = 20)
#' @param output_dir Directory where CS results will be written. Only specify if Circuitscape current map outputs are requested. It is critical that there are NO SPACES in the directory, as this may cause the function to fail.
#' @param output_name Name to be applied to output files. Defaults to using the raster layer names.
#' @param cholmod (Logical; Default = TRUE). Should the cholmod solver be used?
#' @param precision (Logical; Default = FALSE). Should experimental single precision method be used? See details.
#' @param is_resistance Default = TRUE. Is the landscape represented as a resistance (TRUE) or conductance (FALSE) surface?
#' @param remove_files Remove temporary files. Default is TRUE.
#' @param silent Printing of output from Circuitscape will be suppressed unless an error occurs (Default = TRUE)
#'
#' @return Function will return an omnidirectional current flow SpatRaster.
#'
#' @details There is extensive documentation for Circuitscape here: https://docs.circuitscape.org/Circuitscape.jl/latest/ . This function follows guidelines outlined by Koen et al. 2014 (https://besjournals.onlinelibrary.wiley.com/doi/full/10.1111/2041-210X.12197). Briefly, this function will extend a resistance surface and fill that spaces with random positive values. Then, points are distributed around this landscape and an 'all-to-one' Circuitscape analysis is run. Current surfaces are then trimmed to their original extent.
#' @export
#' @examples
#' f <- system.file("data/resist.tif", package = "landconn")
#' r <- terra::rast(f)
#' jl_home <- "C:/Users/peterman.73/AppData/Local/Programs/Julia-1.10.5/bin/"
#' omni_r <- omni_cs(JULIA_HOME = jl_home,
#'                   r,
#'                   n_nodes = 20,
#'                   output_dir = tempdir(),
#'                   output_name = 'test',
#'                   cholmod = TRUE,
#'                   precision = FALSE,
#'                   is_resistance = TRUE,
#'                   remove_files = TRUE,
#'                   silent = TRUE)
#' plot(omni_r)
#'
#' @usage
#' omni_cs(JULIA_HOME = NULL,
#'         r,
#'         n_nodes = 20,
#'         output_dir = NULL,
#'         output_name = NULL,
#'         cholmod = TRUE,
#'         precision = FALSE,
#'         is_resistance = TRUE,
#'         remove_files = TRUE,
#'         silent = TRUE)

omni_cs <- function(JULIA_HOME = NULL,
                    r,
                    n_nodes = 20,
                    output_dir = NULL,
                    output_name = NULL,
                    cholmod = TRUE,
                    precision = FALSE,
                    is_resistance = TRUE,
                    remove_files = TRUE,
                    silent = TRUE) {

  if(class(r) == 'RasterLayer'){
    r <- terra::rast(r)
  }

  d <- dim(r)[-3]
  # e <- ext(r)
  # e_x <- e *1.25
  d <- floor(d*0.25)
  r_x <- terra::extend(r, d, fill = NA)
  n_ <- global(is.na(r_x), fun='sum')[[1]]
  r_x[is.na(r_x)] <- rnorm(n = n_, 5, 1)

  e <- ext(r_x) * 0.995
  x_ <- seq(e[1], e[2], length.out = ceiling(n_nodes / 4))
  y_ <- seq(e[3], e[4], length.out = ceiling(n_nodes / 4))
  xy1 <- expand.grid(x_, c(e[3:4]))
  xy <- rbind(expand.grid(c(e[1:2]), y_), xy1)
  pts <- vect(xy, geom = c('Var1', 'Var2'))
  # plot(r_x)
  # plot(pts, add = T, col = 'red')
  xy_sf <- sf::st_as_sf(pts)

  suppressMessages(omni_out <- run_cs(JULIA_HOME = JULIA_HOME,
                                      rast = r_x,
                                      input_locs = xy_sf,
                                      output_dir,
                                      output_name,
                                      cholmod,
                                      precision,
                                      is_resistance,
                                      remove_files,
                                      silent,
                                      regions = NULL,
                                      field = NULL,
                                      CurrentMap = TRUE,
                                      cumulative_map_only = TRUE,
                                      VoltageMaps = FALSE,
                                      scenario = 'all-to-one',
                                      Neighbor_Connect = 8,
                                      output = "raster",
                                      pairs_to_include = NULL,
                                      variable_source = NULL,
                                      variable_ground = NULL,
                                      parallel = FALSE,
                                      cores = NULL,
                                      focal_node_current_zero = FALSE,
                                      max_map = FALSE))

  omni_rast <- terra::crop(omni_out, r)
  return(omni_rast)
}
