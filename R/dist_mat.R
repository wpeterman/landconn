# Distance matrix -----------------------------------------

#' Calculate pairwise distance matrix
#'
#' @param coords Provide a two-column array, matrix, or data.frame with xy coordinates of locations.
#' @return A square distance matrix
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
#' @usage dist_mat(coords)
#'
#' @export
#' @author Bill Peterman <peterman.73@@osu.edu>

dist_mat <- function(coords) {

  if(dim(coords)[2] != 2){
    stop("`coords` must be two columns representing xy coordinates")
  }

  d <- as.matrix(dist(coords))

  return(d)

}
