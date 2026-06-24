# S3 methods for ifc and land_dist objects ------------------------------------

# ---- ifc methods ----

#' Methods for `ifc` connectivity objects
#'
#' Display, summarize, and plot the output of [ifc()]. An `ifc` object is a named
#' numeric vector of relative connectivity values (one per site) that also records
#' which model produced it.
#'
#' @param x,object An object of class `ifc`, as returned by [ifc()].
#' @param sort Logical (Default = TRUE). If TRUE, `plot()` orders sites from most to
#'   least connected, which makes the ranking easy to read. If FALSE, sites are shown
#'   in their original order.
#' @param ... Additional arguments passed to the underlying [graphics::barplot()]
#'   (for `plot()`) or ignored (for `print()` and `summary()`).
#'
#' @return
#'   `print()` returns its input invisibly. `summary()` invisibly returns a list with
#'   the model context and the distribution of connectivity values, and prints a
#'   readable summary as a side effect. `plot()` draws a barplot and returns the bar
#'   midpoints invisibly.
#'
#' @details
#'   Connectivity values are *relative*: they are only meaningful in comparison with
#'   one another within the same call, not as absolute quantities. A site with twice
#'   the value of another is roughly twice as connected to the rest of the network
#'   under the chosen model. If the object was produced with `scale = TRUE`, values run
#'   from 0 to 1, where 1 is the most connected site.
#'
#' @examples
#' set.seed(123)
#' loc <- data.frame(x = runif(10, 0, 10), y = runif(10, 0, 10))
#' d <- dist_mat(loc)
#' conn <- ifc(alpha = 3, dist_mat = d, model = 1)
#'
#' conn            # print method
#' summary(conn)   # distribution and most/least connected sites
#' plot(conn)      # sorted barplot
#'
#' @name ifc-methods
#' @seealso [ifc()]
NULL

#' @rdname ifc-methods
#' @exportS3Method print ifc
print.ifc <- function(x, ...) {
  model <- attr(x, "model")
  cat("Incidence function connectivity (`ifc`)\n")
  cat("  Model:          ", model, "\n", sep = "")
  cat("  alpha:          ", attr(x, "alpha"), "\n", sep = "")
  if(model > 1){
    cat("  area_transform: ", attr(x, "area_transform"), "\n", sep = "")
  }
  cat("  scaled:         ", attr(x, "scaled"), "\n", sep = "")
  cat("  sites:          ", length(x), "\n\n", sep = "")
  cat("Relative connectivity per site:\n")
  print(stats::setNames(as.numeric(x), names(x)))
  invisible(x)
}

#' @rdname ifc-methods
#' @exportS3Method summary ifc
summary.ifc <- function(object, ...) {
  v <- stats::setNames(as.numeric(object), names(object))
  ord <- order(v, decreasing = TRUE)
  n_show <- min(3L, length(v))

  out <- list(
    model = attr(object, "model"),
    alpha = attr(object, "alpha"),
    area_transform = attr(object, "area_transform"),
    scaled = attr(object, "scaled"),
    n_sites = length(v),
    distribution = summary(v),
    most_connected = v[utils::head(ord, n_show)],
    least_connected = v[utils::tail(ord, n_show)]
  )

  cat("Incidence function connectivity (`ifc`) summary\n")
  cat("  Model ", out$model, " | alpha = ", out$alpha,
      if(out$model > 1) paste0(" | area_transform = ", out$area_transform) else "",
      " | scaled = ", out$scaled, "\n", sep = "")
  cat("  ", out$n_sites, " sites\n\n", sep = "")
  cat("Distribution of relative connectivity:\n")
  print(out$distribution)
  cat("\nMost connected site(s):\n")
  print(out$most_connected)
  cat("\nLeast connected site(s):\n")
  print(out$least_connected)

  invisible(out)
}

#' @rdname ifc-methods
#' @exportS3Method plot ifc
plot.ifc <- function(x, sort = TRUE, ...) {
  v <- stats::setNames(as.numeric(x), names(x))
  if(isTRUE(sort)){
    v <- v[order(v, decreasing = TRUE)]
  }
  ylab <- if(isTRUE(attr(x, "scaled"))){
    "Relative connectivity (scaled 0-1)"
  } else {
    "Relative connectivity"
  }
  bp <- graphics::barplot(v,
                          ylab = ylab,
                          xlab = "Site",
                          main = paste0("ifc model ", attr(x, "model")),
                          ...)
  invisible(bp)
}


# ---- land_dist methods ----

#' Methods for `land_dist` distance matrices
#'
#' Display, summarize, and plot the square distance matrix returned by [dist_mat()].
#'
#' @param x,object An object of class `land_dist`, as returned by [dist_mat()].
#' @param ... Additional arguments passed to [graphics::image()] (for `plot()`) or
#'   ignored (for `print()` and `summary()`).
#'
#' @return
#'   `print()` returns its input invisibly. `summary()` invisibly returns the
#'   distribution of the unique pairwise distances (the lower triangle) and prints it.
#'   `plot()` draws a heatmap of the matrix and returns `NULL` invisibly.
#'
#' @details
#'   The diagonal of a `land_dist` matrix is zero (each site is zero distance from
#'   itself). `summary()` reports only the lower triangle, which holds each pairwise
#'   distance once, so the statistics describe the spread of distances between distinct
#'   sites.
#'
#' @examples
#' set.seed(123)
#' loc <- data.frame(x = runif(10, 0, 10), y = runif(10, 0, 10))
#' d <- dist_mat(loc)
#'
#' d              # print method
#' summary(d)     # spread of pairwise distances
#' plot(d)        # heatmap
#'
#' @name land_dist-methods
#' @seealso [dist_mat()], [lower()]
NULL

#' @rdname land_dist-methods
#' @exportS3Method print land_dist
print.land_dist <- function(x, ...) {
  n <- nrow(x)
  cat("Pairwise distance matrix (`land_dist`)\n")
  cat("  ", n, " sites (", n, " x ", n, ")\n\n", sep = "")
  if(n <= 10){
    print(unclass(round(x, 3)))
  } else {
    cat("Showing top-left 6 x 6 corner:\n")
    print(round(unclass(x)[1:6, 1:6], 3))
    cat("...\n")
  }
  invisible(x)
}

#' @rdname land_dist-methods
#' @exportS3Method summary land_dist
summary.land_dist <- function(object, ...) {
  d <- lower(object)
  cat("Pairwise distance matrix (`land_dist`) summary\n")
  cat("  ", nrow(object), " sites, ", length(d), " unique pairwise distances\n\n",
      sep = "")
  cat("Distribution of pairwise distances:\n")
  out <- summary(d)
  print(out)
  invisible(out)
}

#' @rdname land_dist-methods
#' @exportS3Method plot land_dist
plot.land_dist <- function(x, ...) {
  m <- unclass(x)
  n <- nrow(m)
  ## image() plots with the first row at the bottom; flip so site 1 reads top-left.
  graphics::image(1:n, 1:n, t(m)[, n:1],
                  xlab = "Site", ylab = "Site",
                  main = "Pairwise distances",
                  ...)
  invisible(NULL)
}
