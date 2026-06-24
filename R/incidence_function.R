
# Incidence Function Connectivity -----------------------------------------

#' Calculate connectivity using incidence functions
#'
#' @param alpha The average dispersal distance of your organism. Should be in the same units as the distances in the `dist_mat`. See Details
#' @param dist_mat A square pairwise distance matrix
#' @param model Numeric value 1-4 corresponding to how connectivity will be calculated. See Details
#' @param occ_sites A 0/1 vector equal in length to the number of sites in the `dist_mat`. If provided, only sites indicated with a 1 will contribute to connectivity. If not provided, all sites contribute to connectivity.
#' @param patch_area If provided, can be used to calculate connectivity as a function of the area of surrounding patches (Model 2), focal patches (Model 3), or both (Model 4). Provide a numeric vector equal in length to the number of sites in the `dist_mat`, in whatever area units are meaningful to you. See Details
#' @param scale Logical (Default = FALSE). If TRUE, connectivity estimates will be scaled to have a maximum of 1.
#' @param area_transform How patch area enters Models 2-4. `"none"` (Default) uses `patch_area` directly, matching the equations in Details. `"log"` uses `log(patch_area)` instead, which dampens the influence of very large patches. See Details
#'
#' @return An object of class `ifc`: a named numeric vector of relative connectivity
#'   values, one per site, in the same order as the rows of `dist_mat`. Higher values
#'   indicate sites that are better connected to the rest of the network. The model
#'   number, `alpha`, `area_transform`, and whether the result was scaled are stored as
#'   attributes and shown by `print()` and `summary()`. Because the object is a numeric
#'   vector underneath, it can be used directly in arithmetic, correlation, or plotting.
#'   See [print.ifc()], [summary.ifc()], and [plot.ifc()] for the methods that help you
#'   read it.
#'
#' @details `alpha` is used in conjunction with `dist_mat` in a negative exponential equation to calculate pairwise probability of connectivity:
#' \deqn{C ~ exp(-(1 / alpha) * dist_mat)}.
#'
#' ##### CONNECTIVITY MODELS #####
#' model 1 --> Distance only \cr
#' Calculate connectivity as a function of distance.\cr
#' Hypothesis: Connectivity of patches declines as they become more isolated from other patches
#' \deqn{C_i = \sum_{j \ne i }p_j \text{exp} (-\alpha d_{ij})}
#'
#' model 2 --> Contributing area \cr
#' Calculate connectivity, weighting by the area of contributing patches.\cr
#' Hypothesis: Larger patches contribute more individuals
#' \deqn{C_i = \sum_{j \ne i }p_j \text{exp} (-\alpha d_{ij})A_j}
#'
#' model 3 --> Focal area \cr
#' Calculate connectivity, weighting by the area of the focal patches.\cr
#' Hypothesis: Larger patches have greater attraction
#' \deqn{C_i = A_i\sum_{j \ne i }p_j \text{exp} (-\alpha d_{ij})}
#'
#' model 4 --> Contributing & Focal area \cr
#' Hypothesis: Larger patches contribute more individuals and larger patches have greater attraction \cr
#' \deqn{C_i = A_i\sum_{j \ne i }p_j \text{exp} (-\alpha d_{ij})A_j}
#'
#' \eqn{p_j} (`occ_sites`) is the occupancy of site j and modifies which patches contribute to the connectivity of a patch
#'
#' In the equations above, \eqn{A_j} and \eqn{A_i} are patch areas entered directly, which is the default (`area_transform = "none"`). Setting `area_transform = "log"` substitutes \eqn{\log(A_j)} and \eqn{\log(A_i)} in Models 2-4. Note that `log()` of an area less than 1 is negative, so use `"log"` only when your area values are comfortably greater than 1. Model 1 does not use patch area and is unaffected by this argument.
#'
#' @examples
#' set.seed(123)
#'
#' ## Create random coordinates
#' loc <- data.frame(x = runif(10,0,10),
#'                   y = runif(10,0,10))
#'
#' ## Create distance matrix
#' dist_mat <- as.matrix(dist(loc))
#'
#' ## Set alpha and indicate occupied sites
#' alpha <- 3
#' occ_sites <- c(0,0,1,1,0,1,0,0,1,1)
#'
#' ## Patch Areas
#' p_a <- runif(10, 5, 50)
#'
#'
#' ## Calculate connectivity
#' ## Model 1
#' (c1 <- ifc(alpha,
#'            dist_mat,
#'            model = 1))
#'
#' ## Model 1, scaled
#' (c1s <- ifc(alpha,
#'             dist_mat,
#'             model = 1,
#'             scale = TRUE))
#'
#' ## Model 2
#' (c2 <- ifc(alpha,
#'            dist_mat,
#'            model = 2,
#'            patch_area = p_a))
#'
#' ## Model 3
#' (c3 <- ifc(alpha,
#'            dist_mat,
#'            model = 3,
#'            patch_area = p_a))
#'
#' ## Model 4
#' (c4 <- ifc(alpha,
#'            dist_mat,
#'            model = 4,
#'            patch_area = p_a))
#'
#' ## Model 4, Occupied sites only
#' (c4b <- ifc(alpha,
#'             dist_mat,
#'             model = 4,
#'             occ_sites = occ_sites,
#'             patch_area = p_a))
#'
#' ## Model 2 with log-transformed area (pre-0.3.0 behavior)
#' (c2log <- ifc(alpha,
#'               dist_mat,
#'               model = 2,
#'               patch_area = p_a,
#'               area_transform = "log"))
#'
#' ## The result is an `ifc` object with print, summary, and plot methods
#' summary(c4)
#' plot(c4)
#'
#' @usage
#' ifc(alpha,
#'     dist_mat,
#'     model,
#'     occ_sites = NULL,
#'     patch_area = NULL,
#'     scale = FALSE,
#'     area_transform = c("none", "log"))
#'
#' @seealso [dist_mat()] to build the input distance matrix, [print.ifc()],
#'   [summary.ifc()], and [plot.ifc()] for the methods.
#'
#' @export
#' @author Bill Peterman <peterman.73@@osu.edu>

ifc <- function(alpha,
                dist_mat,
                model,
                occ_sites = NULL,
                patch_area = NULL,
                scale = FALSE,
                area_transform = c("none", "log")) {

  area_transform <- match.arg(area_transform)

  ## Checks
  if(inherits(dist_mat, "matrix")){
    if(dim(dist_mat)[1] != dim(dist_mat)[2]){
      stop("`dist_mat` must be a square distance matrix")
    }
  } else {
    stop("`dist_mat` must be a square distance matrix")
  }


  if(missing(alpha) || !is.numeric(alpha)){
    stop('A numeric value for `alpha` must be provided')
  }

  if(missing(model) || !is.numeric(model)){
    stop('A numeric value for `model` must be provided')
  }

  if(!(model >= 1 & model <= 4)){
    stop('`model` must be a numeric value 1-4. See Details.')
  }

  if(model > 1 && is.null(patch_area)){
    stop('`patch_area` must be provided for models 2, 3, and 4. See Details.')
  }

  if(!is.null(occ_sites)){
    if(length(occ_sites) != dim(dist_mat)[1]){
      stop("occ_sites must be a vector of equal length as sites in the dist_mat")
    }
  }

  if(!is.null(patch_area)){
    if(length(patch_area) != dim(dist_mat)[1]){
      stop("patch_area must be a vector of equal length as sites in the dist_mat")
    }
  }


  # Connectivity probability ------------------------------------------------

  c_prob <-
    exp(-(1 / alpha) * dist_mat)  ## calculate connectivity probability
  c_prob <- as.matrix(c_prob)
  diag(c_prob) <- 0

  if(is.null(occ_sites)){
    occ_sites <- rep(1, dim(dist_mat)[1])
  }

  c_prob <- sweep(c_prob, 2, occ_sites, "*")

  ## Area weight entering Models 2-4. `"none"` uses area directly (the documented
  ## equations); `"log"` dampens large patches but goes negative for areas < 1.
  if(!is.null(patch_area)){
    w <- if(area_transform == "log") log(patch_area) else patch_area
  }



  # ** Method 1 -------------------------------------------------------------
  if(model == 1){

    connectivity <- rowSums(c_prob)

  }

  # ** Method 2 -------------------------------------------------------------
  if(model == 2){

    connectivity <- rowSums(sweep(c_prob, 2, w, "*"))

  }

  # ** Method 3 -------------------------------------------------------------

  if(model == 3){

    connectivity <- w * rowSums(c_prob)

  }

  # ** Method 4 -------------------------------------------------------------

  if(model == 4){

    connectivity <- rowSums(w * sweep(c_prob, 2, w, "*"))

  }

  if(isTRUE(scale)){
    connectivity <- connectivity / max(connectivity)
  }

  ## Return a classed `ifc` object: a numeric vector carrying its model context.
  if(is.null(names(connectivity))){
    names(connectivity) <- seq_along(connectivity)
  }

  structure(connectivity,
            class = "ifc",
            model = model,
            alpha = alpha,
            area_transform = area_transform,
            scaled = isTRUE(scale))
}
