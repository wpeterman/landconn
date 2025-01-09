#' Launch the Julia installation from R and load the necessary libraries in Julia
#'
#' This function starts a Julia session from R and imports the `ConScape` and `Circuitscape` library
#' to be used from R. It assumes that Julia is already installed in the system
#' and the path to its executables is given as the `julia_path` argument.
#' The first time the function is ran, it is best to set `install_libraries = TRUE`
#' to install the libraries in Julia.
#'
#' @param julia_path `[character]` \cr The directory for the Julia bin, e.g.
#' "C:/Programs/Julia-1.9.3/bin".
#' @param install_libraries `[logical]` \cr If `FALSE`, Julia will be
#' launched and the required libraries will be loaded without installing them.
#' If `TRUE` (default), the libraries will be (re-)installed in Julia.
#'
#' @return NULL.
#' @export
#'
julia_packages <- function(julia_path, install_libraries = TRUE) {

  Sys.setenv(JULIA_BINDIR = julia_path)

  if (install_libraries){
    Pkg <- JuliaConnectoR::juliaImport("Pkg")
    JuliaConnectoR::juliaEval("Pkg.add(\"ConScape\")")
    JuliaConnectoR::juliaEval("Pkg.add(\"SparseArrays\")")
    JuliaConnectoR::juliaEval("Pkg.add(\"Plots\")")
    JuliaConnectoR::juliaEval("Pkg.add(\"Circuitscape\")")
    # JuliaConnectoR::juliaEval("Pkg.add(\ name=\"Circuitscape\"\ ,version=\"5.13.3\")")

  }
  SA <- JuliaConnectoR::juliaImport("SparseArrays")
  CS <- JuliaConnectoR::juliaImport("ConScape")
  CIRC <- JuliaConnectoR::juliaImport("Circuitscape")
}
#' @importFrom JuliaConnectoR juliaImport juliaEval juliaImport
