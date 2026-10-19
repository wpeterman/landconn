# landconn review remediation

This batch addresses the 19 findings in the October 19, 2026 review of 0.6.1.
The work is grouped by behavior; one patch release and commit covers the batch.

## Spatial inputs and Circuitscape contracts

| Finding | Resolution | Regression evidence |
|---|---|---|
| F01 | Normalize supported spatial classes once, preserve explicit IDs, validate fields and cells, and name pairwise matrix dimensions by solver IDs. Force double precision at the Julia ASCII-read boundary and restore its prior setting. | Point-class, reader-restoration and ordering checks in `test-run_cs-contracts.R`; live named matrices and IDs above the signed 32-bit range. |
| F02 | Sample represented, valid raster cells per region ID, retaining the region/node relationship. | Region sampling with an internal missing cell and sf/SpatVector inputs. |
| F03 | Return pairwise, per-node and advanced results according to their actual shapes, without stripping scenario-specific rows. | Synthetic shape regressions and live tests of all four scenarios. |
| F04 | Separate advanced source/ground grids from per-node strength tables. One injection/ground per collapsed core preserves the requested total and ground resistance. | Config/grid checks plus live polygon strengths and advanced voltage/current runs. |
| F05 | Fill only newly added buffer cells, retain original masks through solves and return maps, and construct unique perimeter nodes. | Mask/node-count regressions and live omnidirectional map with a hole. |
| F10 | Read only ASCII raster output files, exclude sidecars, and assemble consistent raster classes and CRS. | Fake sidecar regression, legacy Raster output test, and live current/voltage layers. |
| F11 | Honor lower-triangle pair selection and logical/string precision; serialize booleans in the solver's accepted form. Configure/check real Julia thread counts. | Configuration checks; serial/two-thread agreement; live pair selection. |
| F12 | Write configuration without sinks, preserve original solver errors, clean temporary inputs on failure, and restore R environment variables. | Sink, environment and failing-solver regressions. |

## Statistical fitting and uncertainty

| Finding | Resolution | Regression evidence |
|---|---|---|
| F06 | Refit GLM/NB calls while preserving formula, weights, offsets, subsets and controls; retain unmarked update behavior. | Independent weighted/offset/subset GLM reference; weighted NB template. |
| F07 | Separate network size from likelihood sample size, align simulated responses to fitted rows, and preserve missingness/subsets. | Missing-response nobs/BIC and bootstrap tests; supported unmarked response-matrix bootstraps. |
| F08 | Reject non-converged, separated or non-finite candidates; retain reasons and report rejected-fit counts. | Deliberately insufficient-iteration GLM and engine diagnostics checks. |
| F09 | Reject transformed/structured GLM/NB response expressions when bootstrap is requested, rather than silently transforming simulated values twice. | Early, actionable transformed-response bootstrap rejection. Ordinary optimization remains supported. |
| F14 | Validate distances, coordinates, alpha, model codes, occupancy, areas and optimization controls. Leave all-zero scaled connectivity at zero. | Scientific-domain and zero-network regressions. |
| F15 | Support single-season unmarked occu/pcount explicitly and require connectivity in the state formula, not detection. | Unused-formula rejection, missingness restoration, actual occu/pcount bootstrap refits. |
| F16 | Check response/row/weight/likelihood comparability, reject entirely undefined rankings, retain unranked rows explicitly, and report disjoint/truncated profile support. | Incompatible-response and undefined-AICc errors; interval metadata tests. |

## Mathematics, documentation and maintenance

| Finding | Resolution | Validation |
|---|---|---|
| F13 | Use `exp(-distance / alpha)` consistently in all four equations and distinguish association scale from independently measured dispersal. | Regenerated help and independent valid-input formula regressions. |
| F17 | Align the reversed heatmap rows with site labels and display the distance range/color interpretation. | Executed plotting examples and visual inspection of rendered vignette figures. |
| F18 | Document all returned fields, scenario-specific outputs, supported engines, masks, node counts, stochastic buffer, limitations, and unknown dataset provenance. Add the raster help topic and two self-contained guides. Qualify the Koen et al. adaptation. | Documentation regeneration, three built vignettes, runnable examples, spelling/URL checks and figure inspection. |
| F19 | Add R-release/platform and R-devel CI, optional Windows/Linux Julia integration with a pinned solver version, a technical word list, and minimal Julia installation. | Local R suites and source-archive checks, optional real-solver tests, workflow syntax review and setup mocks. Remote CI runs remain pending a push. |

## Intentional contract corrections

- Invalid IDs, disappearing cores, colliding point cells, invalid scientific inputs,
  unsupported options and unsupported model classes now fail explicitly.
- Focal IDs are preserved. Pairwise returns have ID dimension names; one-to-all
  and all-to-one return complete two-column tables. Advanced mode returns a full
  voltage grid. SpatRaster/legacy Raster return classes follow the input.
- Advanced tables contain per-core total source current and resistance to ground,
  with disjoint source/ground IDs. Only one point cell per ID is allowed unless
  polygon cells are collapsed into a short-circuit region.
- `omni_cs()` now preserves original nodata and uses exactly the requested number
  of unique perimeter cells. This intentionally corrects the previous circuit.
- Bootstrap response transformations are explicitly restricted. Rejected fits,
  failed replicates, disjoint profile support and truncated limits are visible.
- Incompatible model comparisons and undefined best-model rankings are rejected.
  Older optimization objects without comparison metadata must be refit.
- Julia installation loads only Circuitscape. Thread counts are checked against
  the actual process, and R environment variables/output redirection are restored.

## Validation record

- Regenerated help with `devtools::document()`. The installed R regression suite
  passed 257 assertions with no failures or warnings. Its optional Julia test was
  skipped there and run separately against the installed release archive.
- Real Julia 1.12.6 / Circuitscape 5.17.1 integration: 20 assertions passed, with
  no failures, warnings or skips. Includes polygon IDs above 2^31, four scenarios,
  selected pairs, serial/two-thread agreement, output maps and internal nodata.
- All four documentation URLs passed; spelling was clean. Both CI workflow files
  parsed successfully. Three vignettes built, HTML anchors were checked, and seven
  embedded figures were visually inspected. Full browser-page rendering was not
  independently verified.
- Release archive SHA-256:
  `71801925644d53b71d34161d36999ccc57fe4ee94be082b47a94c6135ace439c`.
  Its 57 source, help and test files match the checkout (normalizing line endings);
  DESCRIPTION differs only by the standard build formatting and generated fields.
- `devtools::check(args = "--no-manual", document = FALSE)` on R 4.6.1:
  **0 errors, 0 warnings, 0 notes**. Examples, longer `donttest` examples,
  installed tests and all three vignette rebuilds passed.
- The exact same archive checked with `--no-manual --as-cran` on R 4.5.3:
  **0 errors, 0 warnings, 2 notes**. Both were administrative/environment notes:
  CRAN incoming classified it as a new submission, and the network time check
  could not verify current time. All 257 assertions and three vignette rebuilds
  passed. The R 4.6.1 devtools check uses its standard disabled incoming checks;
  the R 4.5.3 direct archive check left those enabled.

Local checks do not establish Linux/macOS, R-devel, or hosted CI success; those
jobs are configured for execution when the changes are pushed. Julia compatibility
is verified for the explicitly tested version, not for arbitrary future releases.
