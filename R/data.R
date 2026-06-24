#------------------------------------------------------------------------------
# Sample datasets

#' Sample resistance raster
#'
#' A small example raster representing resistance to movement across a landscape,
#' where higher cell values are harder to move through. It is shipped as a GeoTIFF in
#' the package's `extdata` folder and is used in the examples and the getting-started
#' vignette. Load it with [terra::rast()] on the path returned by [system.file()].
#'
#' @format A single-layer GeoTIFF (`.tif`) read as a `SpatRaster` by [terra::rast()].
#'
#' @examples
#' f <- system.file("extdata/resist.tif", package = "landconn")
#' r <- terra::rast(f)
#' terra::plot(r)
#'
#' @name resist.tif
#'
