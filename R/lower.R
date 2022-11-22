# Lower matrix -----------------------------------------
#'
#' Extract the lower half of a distance matrix
#'
#' @param dist_mat Square distance matrix
#' @return Vector of values from lower half of matrix
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
#' d_vec <- lower(d_mat)
#'
#' @usage lower(dist_mat)
#'
#' @export
#' @author Bill Peterman <peterman.73@@osu.edu>

lower <- function(dist_mat) {

  if(dim(dist_mat)[1] != dim(dist_mat)[1]){
    stop(cat("The number of rows does not equal the number of columns. \n`dist_mat` must be a square matrix"))
  }

  d <- dist_mat[lower.tri(dist_mat)]

  return(d)

}
