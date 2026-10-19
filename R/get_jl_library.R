#' Launch the Julia installation from R and load the necessary libraries in Julia
#'
#' This function starts a Julia session from R and loads the `Circuitscape` library
#' to be used from R. It assumes that Julia is already installed in the system
#' and the path to its executables is given as the `julia_path` argument.
#' Set `install_libraries = TRUE` when Circuitscape needs to be installed in the
#' active Julia project. Later calls can load it without installation.
#'
#' @param julia_path `[character]` \cr The directory for the Julia bin, e.g.
#' "C:/path/to/Julia/bin".
#' @param install_libraries `[logical]` \cr If `FALSE`, Julia will be
#' launched and the required libraries will be loaded without installing them.
#' If `TRUE` (default), the libraries will be (re-)installed in Julia.
#'
#' @return Returns `NULL` invisibly after starting Julia and loading (optionally
#'   installing) Circuitscape. No ConScape or plotting libraries are installed.
#'   The caller's JULIA_BINDIR value is restored; the Julia process persists.
#'
#' @details Install only when needed, since installation changes the active Julia
#'   package environment and requires network access. For repeatable projects,
#'   manage versions in a Julia project environment. The integration tests currently
#'   pin Circuitscape 5.17.1 and use Julia 1.12.
#' @seealso [run_cs()], [omni_cs()]
#' @examples
#' \dontrun{
#' ## Point to the Julia 'bin' folder on your own machine.
#' ## Use install_libraries = TRUE the first time, FALSE thereafter.
#' julia_packages(julia_path = "C:/path/to/Julia/bin",
#'                install_libraries = TRUE)
#' }
#' @export
#'
julia_packages <- function(julia_path, install_libraries = TRUE) {

  if (!is.character(julia_path) || length(julia_path) != 1L || !dir.exists(julia_path)) {
    stop("`julia_path` must be an existing Julia binary directory.")
  }
  if (!is.logical(install_libraries) || length(install_libraries) != 1L || is.na(install_libraries)) {
    stop("`install_libraries` must be TRUE or FALSE.")
  }
  old <- Sys.getenv("JULIA_BINDIR", unset = NA_character_)
  on.exit(if (is.na(old)) Sys.unsetenv("JULIA_BINDIR") else Sys.setenv(JULIA_BINDIR = old))
  Sys.setenv(JULIA_BINDIR = julia_path)
  if (install_libraries) JuliaConnectoR::juliaEval('using Pkg; Pkg.add("Circuitscape")')
  JuliaConnectoR::juliaEval("using Circuitscape")
  invisible(NULL)
}
#' @importFrom JuliaConnectoR juliaImport juliaEval juliaImport
