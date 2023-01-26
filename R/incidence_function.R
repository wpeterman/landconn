
# Incidence Function Connectivity -----------------------------------------

#' Calculate connectivity using incidence functions
#'
#' @param alpha The average dispersal distance of your organism. Should be in the same units as the distances in the `dist_mat`. See Details
#' @param dist_mat A square pairwise distance matrix
#' @param model Numeric value 1-4 corresponding to how connectivity will be calculated. See Details
#' @param occ_sites A 0/1 vector equal in length to the number of sites in the `dist_mat`. If provided, only sites indicated with a 1 will contribute to connectivity. If not provided, all sites contribute to connectivity.
#' @param patch_area If provided, can be used calculate connectivity as a function of surrounding patched (Model 2), focal patches (Model 3), or both (Model 4). See details
#' @param scale Logical (Default = FALSE). If TRUE, connectivity estimates will be scaled to have a maximum of 1.
#'
#' @return A vector of relative connectivity values for each site
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
#' @usage ifc(alpha,
#'            dist_mat,
#'            model,
#'            occ_sites = NULL,
#'            patch_area = NULL,
#'            scale = FALSE)
#'
#' @export
#' @author Bill Peterman <peterman.73@@osu.edu>

ifc <- function(alpha,
                dist_mat,
                model,
                occ_sites = NULL,
                patch_area = NULL,
                scale = FALSE) {

  ## Checks
  if(any(class(dist_mat) == "matrix")){
    if(dim(dist_mat)[1] != dim(dist_mat)[2]){
      stop("`dist_mat` must be a square distance matrix")
    }
  } else {
    stop("`dist_mat` must be a square distance matrix")
  }


  if(!exists('alpha')){
    stop('A numeric value for `alpha` must be provided')
  }

  if(!exists('model') | !is.numeric(model)){
    stop('A numeric value for `model` must be provided')
  }

  if(!(model >= 1 & model <= 4)){
    stop('`model` must be a numeric value 1-4. See Details.')
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



  # ** Method 1 -------------------------------------------------------------
  if(model == 1){

    connectivity <- rowSums(c_prob)

  }

  # ** Method 2 -------------------------------------------------------------
  if(model == 2){

    connectivity <- rowSums(sweep(c_prob, 2, log(patch_area), "*"))

  }

  # ** Method 3 -------------------------------------------------------------

  if(model == 3){

    connectivity <- log(patch_area) * rowSums(c_prob)

  }

  # ** Method 4 -------------------------------------------------------------

  if(model == 4){

    connectivity <- rowSums(log(patch_area) * sweep(c_prob, 2, log(patch_area), "*"))

  }

  if(isTRUE(scale)){
    return(connectivity / max(connectivity))
  } else {
    connectivity
  }
}


# EXAMPLE -----------------------------------------------------------------
# set.seed(123)
# loc <- data.frame(x = runif(10,0,10),
#                   y = runif(10,0,10))
# dist_mat <- as.matrix(dist(loc))
# alpha <- 3
# occ_sites <- c(0,0,1,1,0,1,0,0,1,1)
#
# (c1 <- ifc(alpha,
#            dist_mat,
#            model = 5))
# (c1s <- ifc(alpha,
#             dist_mat,
#             model = 1,
#             scale = T))
#
# plot(c1 ~ c1s)
# cor(c1, c1s)
#
# (c1b <- ifc(alpha,
#             dist_mat,
#             model = 1,
#             occ_sites,
#             scale = T))
#
# ## View probability
# d <- seq(1, max(dist_mat), length.out = 25)
# prob <- exp(-(1 / alpha) * d)
# plot(prob ~ d, type = 'l')
