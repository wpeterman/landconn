# landconn 0.3.0
* `ifc()` now returns an `ifc` object with `print()`, `summary()`, and `plot()` methods.
* `dist_mat()` now returns a `land_dist` object with `print()`, `summary()`, and `plot()` methods. It still behaves as an ordinary matrix.
* **Behavior change:** `ifc()` Models 2-4 now use patch area directly by default, matching the documented incidence function equations. The previous log-area behavior is available with the new `area_transform = "log"` argument.
* `ifc()` now errors clearly when `model` or `patch_area` is missing or invalid (the previous checks never triggered).
* Fixed `lower()` so it correctly rejects non-square matrices.
* Circuitscape examples (`run_cs()`, `omni_cs()`, `julia_packages()`) are now wrapped in `\dontrun{}` so they are not executed by `R CMD check`.
* Added a getting-started vignette and package-level help (`?landconn`).
* Added a `testthat` test suite.
* Moved the sample raster to `inst/extdata/resist.tif`; load it with `system.file("extdata/resist.tif", package = "landconn")`.
* Removed the stale `inst/UNICOR` files (UNICOR support was dropped in 0.2.1) and the bundled `inst/Conefor_command_line` binaries.
* DESCRIPTION cleanup: `Authors@R`, revised Description, dependencies split into Imports/Suggests.
* `write_ini()` is now internal (no longer exported).

# landconn 0.2.3
* Updated functions to use `JuliaConnectoR`

# landconn 0.2.2
* Added `omni_cs` function to calculated omnidirectional current flows across resistance surfaces.

# landconn 0.2.1
* Removed reticulate and Unicor functions

# landconn 0.2.0
* Added JuliaConnectoR function

# landconn 0.1.4
* Updated code and added conefor files

# landconn 0.1.3
* Updated code and imports of necessary functions

# landconn 0.1.2

* Update `ifc` function: Simplified to only requiring one patch area vector; corrected error in connectivity matrix calculation

# landconn 0.1.1

* Added function to calculate pairwise distance matrices


# landconn 0.1.0

* Added a `NEWS.md` file to track changes to the package.
