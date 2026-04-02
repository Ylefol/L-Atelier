###############################################################################
########### Grid Utility Functions ###########
###############################################################################

#' Create Coarse Parameter Grid
#'
#' @description Generates a coarse grid of parameter combinations for initial
#' exploration in cross-validation or grid search. Creates all combinations
#' of the specified values across n parameters.
#'
#' @param n_params Integer. Number of parameters to create grid for.
#' @param values Numeric vector. Values to use for each parameter
#'   (default = c(0.3, 0.5, 0.7)).
#'
#' @return Data frame with n_params columns (param1, param2, ..., paramN)
#'   containing all combinations of the specified values. Each row represents
#'   one parameter combination. Total rows = length(values)^n_params.
#'
#' @details
#' This function is useful for initial broad exploration of parameter space
#' before refining the search with a fine grid around promising regions.
#'
#' For example, with 3 parameters and 3 values each, this creates 27 (3^3)
#' combinations. This coarse grid helps identify the general region of
#' optimal parameters without excessive computation.
#'
#' @export
#'
#' @examples
#' \dontrun{
#' # Create coarse grid for 3 tau parameters
#' coarse_grid <- DEMETER_create_coarse_grid(n_params = 3)
#' # Returns 27 combinations of (0.3, 0.5, 0.7)
#'
#' # Custom coarse values
#' coarse_grid <- DEMETER_create_coarse_grid(
#'   n_params = 2,
#'   values = c(0.1, 0.5, 0.9)
#' )
#'
#' }
DEMETER_create_coarse_grid <- function(n_params, values = c(0.3, 0.5, 0.7)) {

  # Input validation
  if (!is.numeric(n_params) || n_params < 1) {
    stop("n_params must be a positive integer")
  }

  if (!is.numeric(values) || length(values) == 0) {
    stop("values must be a numeric vector with at least one element")
  }

  # Create named list for expand.grid
  grid_args <- vector("list", n_params)
  names(grid_args) <- paste0("param", 1:n_params)

  # Fill each parameter with the same values
  for (i in 1:n_params) {
    grid_args[[i]] <- values
  }

  # Generate grid
  grid <- do.call(expand.grid, grid_args)

  cat(sprintf("[DEMETER] Created coarse grid: %d combinations across %d parameters\n",
              nrow(grid), n_params))

  return(grid)
}


#' Create Fine Parameter Grid Around Specific Ranges
#'
#' @description Generates a fine grid of parameter combinations within
#' specified ranges. Used for refined search after coarse grid identifies
#' promising parameter regions.
#'
#' @param ranges List of numeric vectors. Each element is a 2-element vector
#'   c(min, max) specifying the range for one parameter.
#'   Length of list determines number of parameters.
#' @param n_points Integer. Number of points to generate within each range
#'   (default = 5). Points are evenly spaced using seq().
#'
#' @return Data frame with length(ranges) columns (param1, param2, ..., paramN)
#'   containing all combinations of evenly-spaced points within the specified
#'   ranges. Total rows = n_points^length(ranges).
#'
#' @details
#' This function enables focused exploration of promising parameter regions
#' identified by coarse grid search. Each parameter can have a different
#' range, allowing asymmetric refinement.
#'
#' Points within each range are generated using seq(min, max, length.out = n_points),
#' ensuring even spacing including the endpoints.
#'
#' @export
#'
#' @examples
#' \dontrun{
#' # Fine grid around promising coarse values
#' fine_grid <- DEMETER_create_fine_grid(
#'   ranges = list(
#'     c(0.4, 0.6),  # param1: refine around 0.5
#'     c(0.6, 0.8),  # param2: refine around 0.7
#'     c(0.4, 0.6)   # param3: refine around 0.5
#'   ),
#'   n_points = 5
#' )
#' # Returns 125 (5^3) combinations
#'
#' # Different ranges per parameter
#' fine_grid <- DEMETER_create_fine_grid(
#'   ranges = list(c(0.1, 0.3), c(0.7, 0.9)),
#'   n_points = 10
#' )
#'
#' }
DEMETER_create_fine_grid <- function(ranges, n_points = 5) {

  # Input validation
  if (!is.list(ranges) || length(ranges) == 0) {
    stop("ranges must be a non-empty list")
  }

  # Validate each range
  for (i in 1:length(ranges)) {
    if (!is.numeric(ranges[[i]]) || length(ranges[[i]]) != 2) {
      stop(sprintf("ranges[[%d]] must be a numeric vector of length 2", i))
    }
    if (ranges[[i]][1] >= ranges[[i]][2]) {
      stop(sprintf("ranges[[%d]]: min must be less than max", i))
    }
  }

  if (!is.numeric(n_points) || n_points < 2) {
    stop("n_points must be an integer >= 2")
  }

  # Create named list for expand.grid
  n_params <- length(ranges)
  grid_args <- vector("list", n_params)
  names(grid_args) <- paste0("param", 1:n_params)

  # Generate evenly-spaced points for each parameter
  for (i in 1:n_params) {
    grid_args[[i]] <- seq(ranges[[i]][1], ranges[[i]][2], length.out = n_points)
  }

  # Generate grid
  grid <- do.call(expand.grid, grid_args)

  cat(sprintf("[DEMETER] Created fine grid: %d combinations across %d parameters\n",
              nrow(grid), n_params))

  return(grid)
}


#' Create Ranges Around Best Parameter Values
#'
#' @description Helper function to automatically generate ranges for fine grid
#' search centered around best parameter values from coarse grid search.
#'
#' @param best_values Numeric vector. Best parameter values from coarse search.
#' @param width Numeric. Total width of range around each value (default = 0.2).
#'   Range will be `value - width/2` to `value + width/2`.
#' @param bounds Numeric vector of length 2. Lower and upper bounds to clip
#'   ranges (default = c(0, 1)). Useful for parameters with valid ranges
#'   (e.g., tau in 0 to 1).
#'
#' @return List of numeric vectors. Each element is c(min, max) for one
#'   parameter, suitable for input to DEMETER_create_fine_grid().
#'
#' @details
#' This convenience function automates the creation of fine grid ranges by:
#' \itemize{
#'   \item Centering each range around the corresponding best value
#'   \item Ensuring symmetric exploration (±width/2)
#'   \item Clipping to valid parameter bounds
#' }
#'
#' This is particularly useful in two-stage grid search:
#' \enumerate{
#'   \item Run coarse grid CV to find best parameters
#'   \item Use this function to create ranges around best values
#'   \item Run fine grid CV for refined tuning
#' }
#'
#' @export
#'
#' @examples
#' \dontrun{
#' # After coarse CV, best tau values were (0.5, 0.7, 0.5)
#' ranges <- DEMETER_create_ranges_around_values(
#'   best_values = c(0.5, 0.7, 0.5),
#'   width = 0.2,
#'   bounds = c(0, 1)
#' )
#' # Returns: list(c(0.4, 0.6), c(0.6, 0.8), c(0.4, 0.6))
#'
#' # Then create fine grid
#' fine_grid <- DEMETER_create_fine_grid(ranges, n_points = 5)
#'
#' # Larger search width
#' ranges <- DEMETER_create_ranges_around_values(
#'   best_values = c(0.3, 0.8),
#'   width = 0.4,
#'   bounds = c(0, 1)
#' )
#' # Returns: list(c(0.1, 0.5), c(0.6, 1.0))
#'
#' }
DEMETER_create_ranges_around_values <- function(best_values,
                                                 width = 0.2,
                                                 bounds = c(0, 1)) {

  # Input validation
  if (!is.numeric(best_values) || length(best_values) == 0) {
    stop("best_values must be a numeric vector with at least one element")
  }

  if (!is.numeric(width) || width <= 0) {
    stop("width must be a positive number")
  }

  if (!is.numeric(bounds) || length(bounds) != 2) {
    stop("bounds must be a numeric vector of length 2")
  }

  if (bounds[1] >= bounds[2]) {
    stop("bounds[1] must be less than bounds[2]")
  }

  # Create ranges
  ranges <- vector("list", length(best_values))

  for (i in 1:length(best_values)) {
    # Calculate range centered on best value
    min_val <- best_values[i] - width / 2
    max_val <- best_values[i] + width / 2

    # Clip to bounds
    min_val <- max(min_val, bounds[1])
    max_val <- min(max_val, bounds[2])

    # Ensure min < max (in case best_value is at boundary)
    if (min_val >= max_val) {
      warning(sprintf("Parameter %d: best_value (%.2f) is at boundary. Range set to [%.2f, %.2f]",
                      i, best_values[i], bounds[1], bounds[2]))
      min_val <- bounds[1]
      max_val <- bounds[2]
    }

    ranges[[i]] <- c(min_val, max_val)
  }

  cat(sprintf("[DEMETER] Created %d ranges around best values\n", length(ranges)))
  for (i in 1:length(ranges)) {
    cat(sprintf("  Param %d: [%.3f, %.3f] (centered on %.3f)\n",
                i, ranges[[i]][1], ranges[[i]][2], best_values[i]))
  }

  return(ranges)
}
