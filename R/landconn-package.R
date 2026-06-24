#' landconn: Landscape Connectivity
#'
#' Tools to calculate and explore landscape connectivity, developed for teaching and
#' research. The package covers two complementary approaches: patch-based incidence
#' function connectivity, and resistance-surface connectivity through 'Circuitscape'.
#'
#' @section Main functions:
#' \describe{
#'   \item{[dist_mat()]}{Build a pairwise (Euclidean) distance matrix from site
#'     coordinates. Returns a `land_dist` object.}
#'   \item{[ifc()]}{Calculate incidence function connectivity from a distance matrix,
#'     with four model formulations. Returns an `ifc` object with `print`, `summary`,
#'     and `plot` methods.}
#'   \item{[lower()]}{Extract the lower triangle of a distance matrix as a vector.}
#'   \item{[run_cs()]}{Run 'Circuitscape' from R (via 'JuliaConnectoR') to produce
#'     current maps or resistance distances across a resistance surface.}
#'   \item{[omni_cs()]}{Create an omnidirectional current flow map following the
#'     approach of Koen et al. (2014).}
#'   \item{[julia_packages()]}{Start Julia from R and load the libraries the
#'     'Circuitscape' tools depend on.}
#' }
#'
#' @section Getting started:
#' See `vignette("getting-started", package = "landconn")` for a worked walkthrough.
#'
#' @keywords internal
"_PACKAGE"
