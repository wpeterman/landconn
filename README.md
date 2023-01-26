
# landconn

<!-- badges: start -->
<!-- badges: end -->

This package contains functions to calculate landscape connectivity. These functions have been developed for personal teaching and research purposes and are subject to update/change as needed for these purposes. Documentation and examples may be limited.

## Installation

You can install the development version of `landconn` like so:

``` r
devtools::install_github('wpeterman/landconn')
```

## Example – Incidence function connectivity

This example is copied from the `ifc` example
``` r
library(landconn)

set.seed(123)
## Create random coordinates
loc <- data.frame(x = runif(10,0,10),
                  y = runif(10,0,10))
## Create distance matrix
dist_mat <- as.matrix(dist(loc))

## Set alpha and indicate occupied sites
alpha <- 3
occ_sites <- c(0,0,1,1,0,1,0,0,1,1)

## Patch Areas
p_a <- runif(10, 5, 50)

## Calculate connectivity
## Model 1
(c1 <- ifc(alpha,
           dist_mat,
           model = 1))
           
## Model 1, scaled
(c1s <- ifc(alpha,
            dist_mat,
            model = 1,
            scale = T))

## Model 2
(c2 <- ifc(alpha,
           dist_mat,
           model = 2,
           contrib_area = p_a))

## Model 3
(c3 <- ifc(alpha,
           dist_mat,
           model = 3,
           focal_area = p_a))

## Model 4
(c4 <- ifc(alpha,
           dist_mat,
           model = 4,
           contrib_area = p_a,
           focal_area = p_a))

## Model 4, Occupied sites only
(c4b <- ifc(alpha,
            dist_mat,
            model = 4,
            occ_sites = occ_sites,
            contrib_area = p_a,
            focal_area = p_a))
```

