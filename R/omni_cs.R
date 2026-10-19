############ RUN JULIA ############

#' Create an omnidirectional current map
#'
#' Create an omnidirectional current flow map
#'
#' @param JULIA_HOME Path to the folder containing the Julia binary (See Details).
#' @param r Accepts a RasterLayer object or SpatRaster.
#' @param n_nodes Integer number of distinct perimeter cells to use as focal nodes
#'   (default 20, minimum 2). Must not exceed the number of perimeter cells.
#' @param output_dir Existing result directory, defaulting to `tempdir()`.
#' @param output_name Name to be applied to output files. Defaults to using the raster layer names.
#' @param cholmod (Logical; Default = TRUE). Should the cholmod solver be used?
#' @param precision Logical, `FALSE` for double precision (default), `TRUE` for
#'   single precision. Also accepts `"double"` or `"single"`.
#' @param is_resistance Default = TRUE. Is the landscape represented as a resistance (TRUE) or conductance (FALSE) surface?
#' @param parallel Logical, use Julia threads (default `FALSE`). See [run_cs()].
#' @param cores Optional thread count for `parallel = TRUE`. A running Julia
#'   session must already have that many threads; see [run_cs()].
#' @param remove_files Remove temporary files. Default is TRUE.
#' @param silent Printing of output from Circuitscape will be suppressed unless an error occurs (Default = TRUE)
#' @param ... Not currently implemented
#'
#' @return Function will return an omnidirectional current flow SpatRaster.
#'
#' @importFrom stats rnorm
#' @importFrom terra global ext vect crop extend rast
#'
#' @details This is an adaptation of the perimeter-node approach of Koen et al.
#'   (2014), not a reproduction of their random-node, pairwise analysis. This
#'   implementation uses exactly `n_nodes` cells spaced around the raster perimeter
#'   and all-to-one solves. Each focal node is grounded in turn while the other
#'   nodes supply current, and currents are accumulated across solves.
#'
#'   The raster is extended by 25% of its row and column counts on each side
#'   (rounded down, with at least one added row/column). Only new buffer cells
#'   receive positive `abs(rnorm(mean = 5, sd = 1))` values; missing cells inside
#'   the original extent remain excluded. The output is cropped and masked to the
#'   original raster. Set a seed before calling to reproduce the random buffer.
#'   Repeat with different seeds to assess buffer sensitivity. High current is
#'   conditional on the specified circuit and does not establish movement rates
#'   or the effect of habitat removal.
#'
#'   See `vignette("circuitscape-workflows", package = "landconn")` and
#'   \doi{10.1111/2041-210X.12197} for the original perimeter-node study.
#' @seealso [run_cs()], [julia_packages()]
#' @export
#' @examples
#' \dontrun{
#' library(terra)
#' library(landconn)
#' f <- system.file("extdata/resist.tif", package = "landconn")
#' r <- rast(f)
#'
#' ## Set this to the Julia 'bin' folder on your own machine
#' jl_home <- "C:/Users/<user>/AppData/Local/Programs/Julia-1.10.5/bin/"
#'
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
#' }

omni_cs <- function(JULIA_HOME = NULL,
                    r,
                    n_nodes = 20,
                    output_dir = NULL,
                    output_name = NULL,
                    cholmod = TRUE,
                    precision = FALSE,
                    is_resistance = TRUE,
                    parallel = FALSE,
                    cores = NULL,
                    remove_files = TRUE,
                    silent = TRUE,
                    ...) {

  if (inherits(r, "Raster")) r <- terra::rast(r)
  if (!inherits(r, "SpatRaster") || terra::nlyr(r) != 1L) stop("`r` must be a single-layer raster.")
  if (!is.numeric(n_nodes) || length(n_nodes) != 1L || !is.finite(n_nodes) ||
      n_nodes < 2 || n_nodes != floor(n_nodes)) stop("`n_nodes` must be an integer of at least 2.")
  padding <- pmax(1L, floor(dim(r)[1:2] * 0.25))
  r_x <- terra::extend(r, padding, fill = NA)
  # The extent indicator is 1 even at original nodata cells. Only cells added by
  # extending the domain receive simulated resistance/conductance values.
  inside <- terra::extend(terra::setValues(terra::rast(r), 1), r_x)
  outside <- is.na(terra::values(inside, mat = FALSE))
  r_x[outside] <- pmax(abs(stats::rnorm(sum(outside), 5, 1)), .Machine$double.eps)
  nr <- nrow(r_x); nc <- ncol(r_x)
  ring <- c(seq_len(nc), seq.int(2L, nr) * nc,
            (nr - 1L) * nc + rev(seq_len(nc - 1L)),
            if (nr > 2L) (rev(seq.int(2L, nr - 1L)) - 1L) * nc + 1L else integer())
  if (n_nodes > length(ring)) stop("`n_nodes` exceeds the number of distinct perimeter cells.")
  cells <- ring[floor(seq(0, length(ring), length.out = n_nodes + 1L)[- (n_nodes + 1L)]) + 1L]
  pts <- terra::vect(terra::xyFromCell(r_x, cells), crs = terra::crs(r))
  omni_out <- run_cs(JULIA_HOME = JULIA_HOME, r = r_x, input_locs = pts,
                     output_dir = output_dir, output_name = output_name,
                     cholmod = cholmod, precision = precision, is_resistance = is_resistance,
                     parallel = parallel, cores = cores, remove_files = remove_files,
                     silent = silent, scenario = "all-to-one", output = "raster")
  ans <- terra::crop(omni_out, r)
  terra::mask(ans, r)
}
