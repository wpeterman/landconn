
# landconn

<!-- badges: start -->
<!-- badges: end -->

This package contains functions to calculate landscape connectivity. These functions have been developed for personal teaching and research purposes and are subject to update/change as needed for these purposes. Documentation and examples may be limited.

## Installation

You can install the development version of `landconn` like so:

``` r
devtools::install_github('wpeterman/landconn')
```

## Example: incidence function connectivity

``` r
library(landconn)

set.seed(123)
## Create random coordinates
loc <- data.frame(x = runif(10, 0, 10),
                  y = runif(10, 0, 10))

## Create a pairwise distance matrix (a `land_dist` object)
d_mat <- dist_mat(loc)

## Set alpha (average dispersal distance) and indicate occupied sites
alpha <- 3
occ_sites <- c(0, 0, 1, 1, 0, 1, 0, 0, 1, 1)

## Patch areas
p_a <- runif(10, 5, 50)

## Model 1: distance only
c1 <- ifc(alpha, d_mat, model = 1)

## Model 1, scaled to a maximum of 1
c1s <- ifc(alpha, d_mat, model = 1, scale = TRUE)

## Model 2: weight by contributing patch area
c2 <- ifc(alpha, d_mat, model = 2, patch_area = p_a)

## Model 3: weight by focal patch area
c3 <- ifc(alpha, d_mat, model = 3, patch_area = p_a)

## Model 4: weight by both
c4 <- ifc(alpha, d_mat, model = 4, patch_area = p_a)

## Model 4, occupied sites only
c4b <- ifc(alpha, d_mat, model = 4, occ_sites = occ_sites, patch_area = p_a)

## The result is an `ifc` object with print, summary, and plot methods
c4
summary(c4)
plot(c4)
```

As of version 0.3.0, Models 2 to 4 use patch area directly (the standard incidence
function formulation). To reproduce the earlier log-area behavior, pass
`area_transform = "log"`.

## Circuitscape

`landconn` can also run [Circuitscape](https://docs.circuitscape.org/Circuitscape.jl/latest/)
from R through `JuliaConnectoR` (functions `run_cs()` and `omni_cs()`). These require a
local Julia installation and are demonstrated in
`vignette("getting-started", package = "landconn")`.

