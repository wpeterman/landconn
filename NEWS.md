# landconn 0.5.0
* `ifc_optim()` can now compare several `ifc()` models at once: pass a vector to `model` (for example `1:4`), or both area transforms, and it fits the best `alpha` for each and returns an `ifc_optim_set` with a ranked model selection table.
* Added `ifc_modsel()` to build a model selection table (AIC, AICc, BIC, delta, Akaike weights) from `ifc_optim` fits, with `K` counting the GLM parameters plus the optimized `alpha`.
* Added information criteria methods for `ifc_optim`: `logLik()`, `nobs()`, `AIC()`, `BIC()`, and a new generic `AICc()` (with a default method usable on `glm`/`lm`). All count the `alpha` parameter, so an `ifc_optim` fit compares correctly against an ordinary `glm`.

# landconn 0.4.0
* Added `ifc_optim()` to estimate the incidence scale parameter `alpha` by optimizing the fit of a GLM whose predictor is the connectivity from `ifc()`. Reports a profile likelihood confidence interval and an optional parametric bootstrap interval, and returns an `ifc_optim` object with `print()`, `summary()`, and `plot()` (the alpha profile) methods.

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
