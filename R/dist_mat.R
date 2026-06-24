# Distance matrix -----------------------------------------

#' Calculate pairwise distance matrix
#'
#' Build a square matrix of straight-line (Euclidean) distances between a set of
#' locations, ready to pass to [ifc()].
#'
#' @param coords Provide a two-column array, matrix, or data.frame with xy coordinates of locations. The first column is x (easting/longitude) and the second is y (northing/latitude), in whatever distance units you want the output to be in.
#' @return An object of class `land_dist`: a square symmetric matrix whose entry
#'   \[i, j\] is the distance between site i and site j, with zeros on the diagonal.
#'   The matrix is in the same units as `coords`. It also carries the `"matrix"` class,
#'   so it behaves like an ordinary matrix and can be passed straight to [ifc()] or
#'   [lower()]. See [print.land_dist()], [summary.land_dist()], and [plot.land_dist()]
#'   for the methods that summarize and display it.
#'
#' @examples
#' set.seed(123)
#'
#' ## Create random coordinates
#' loc <- data.frame(x = runif(10,0,10),
#'                   y = runif(10,0,10))
#'
#' ## Create distance matrix
#' d_mat <- dist_mat(coords = loc)
#'
#' ## Inspect it with the class methods
#' summary(d_mat)
#'
#' @usage dist_mat(coords)
#'
#' @seealso [lower()] to extract the lower triangle, [ifc()] to compute connectivity.
#'
#' @export
#' @author Bill Peterman <peterman.73@@osu.edu>

dist_mat <- function(coords) {

  if(dim(coords)[2] != 2){
    stop("`coords` must be two columns representing xy coordinates")
  }

  d <- as.matrix(dist(coords))

  ## Return a classed matrix. Keeping `"matrix"` in the class vector means the object
  ## still works anywhere an ordinary matrix does (e.g. `ifc()`, `lower()`).
  structure(d, class = c("land_dist", "matrix"))

}
