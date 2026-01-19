library(ggplot2)

###############################################################################
########### Color Utilities ###########
###############################################################################

#' Generate Color Gradient
#'
#' @description Generates a gradient of colors from a dark base color to lighter
#' shades. Useful for plotting multiple samples within the same group where
#' visual distinction is needed but with a consistent color theme.
#'
#' @param base_color Character string. The base (darkest) color as a hex code
#'   (e.g., "#08519c") or R color name (e.g., "darkblue").
#' @param n Integer. Number of colors to generate (default = 5).
#' @param direction Character string. Direction of gradient:
#'   "to_light" (default) = dark to light,
#'   "to_dark" = light to dark.
#'
#' @return Character vector of n hex color codes.
#'
#' @details
#' The function converts the base color to HSL (Hue, Saturation, Lightness)
#' color space, then generates n colors by interpolating the lightness value
#' while preserving hue and saturation.
#'
#' For circos plots with multiple samples per group:
#' - Use "to_light" direction
#' - Plot lightest color first, darkest last (dark on top)
#' - This creates clean visuals where similar data is "hidden" under darker lines
#'
#' @export
#'
#' @examples
#' # Generate 6 shades of blue
#' blues <- AETHER_generate_gradient("#08519c", n = 6)
#'
#' # Generate 4 shades of red, light to dark
#' reds <- AETHER_generate_gradient("#d62728", n = 4, direction = "to_dark")
#'
#' # Preview colors
#' barplot(rep(1, 6), col = blues, border = NA)
#'
AETHER_generate_gradient <- function(base_color, n = 5, direction = "to_light") {

  if (n < 1) {
    stop("n must be at least 1")
  }

  if (n == 1) {
    return(base_color)
  }

  if (!direction %in% c("to_light", "to_dark")) {
    stop("direction must be 'to_light' or 'to_dark'")
  }

  # Convert base color to RGB
  rgb_vals <- col2rgb(base_color)[, 1] / 255

  # Convert RGB to HSL
  r <- rgb_vals[1]
  g <- rgb_vals[2]
  b <- rgb_vals[3]

  max_val <- max(r, g, b)
  min_val <- min(r, g, b)
  l <- (max_val + min_val) / 2

  if (max_val == min_val) {
    h <- 0
    s <- 0
  } else {
    d <- max_val - min_val
    s <- if (l > 0.5) d / (2 - max_val - min_val) else d / (max_val + min_val)

    h <- if (max_val == r) {
      ((g - b) / d + (if (g < b) 6 else 0)) / 6
    } else if (max_val == g) {
      ((b - r) / d + 2) / 6
    } else {
      ((r - g) / d + 4) / 6
    }
  }

  # Generate lightness values
  # For "to_light": start at base lightness, go toward 0.95 (very light)
  # For "to_dark": start at 0.95, go toward base lightness
  if (direction == "to_light") {
    l_start <- l
    l_end <- min(0.95, l + (1 - l) * 0.85)  # Don't go pure white
    lightness_vals <- seq(l_start, l_end, length.out = n)
  } else {
    l_start <- min(0.95, l + (1 - l) * 0.85)
    l_end <- l
    lightness_vals <- seq(l_start, l_end, length.out = n)
  }

  # Convert HSL back to RGB for each lightness value
  hsl_to_rgb <- function(h, s, l) {
    if (s == 0) {
      return(c(l, l, l))
    }

    hue_to_rgb <- function(p, q, t) {
      if (t < 0) t <- t + 1
      if (t > 1) t <- t - 1
      if (t < 1/6) return(p + (q - p) * 6 * t)
      if (t < 1/2) return(q)
      if (t < 2/3) return(p + (q - p) * (2/3 - t) * 6)
      return(p)
    }

    q <- if (l < 0.5) l * (1 + s) else l + s - l * s
    p <- 2 * l - q

    r <- hue_to_rgb(p, q, h + 1/3)
    g <- hue_to_rgb(p, q, h)
    b <- hue_to_rgb(p, q, h - 1/3)

    return(c(r, g, b))
  }

  # Generate colors
  colors <- sapply(lightness_vals, function(l_val) {
    rgb_out <- hsl_to_rgb(h, s, l_val)
    rgb(rgb_out[1], rgb_out[2], rgb_out[3])
  })

  return(colors)
}


#' Default Color Palettes
#'
#' @description Returns default color palettes for common omics types and
#' experimental groups. These can be used as base colors for gradient generation.
#'
#' @param type Character string. Type of palette:
#'   "omics" = colors by omics type (ATACseq, CHIPseq, RNAseq),
#'   "groups" = colors by common group names (WT, KO, Control, Treatment).
#'
#' @return Named character vector of hex color codes.
#'
#' @export
#'
#' @examples
#' # Get omics palette
#' omics_colors <- AETHER_default_palette("omics")
#' # Returns: ATACseq="#2ca02c", CHIPseq="#1f77b4", RNAseq="#ff7f0e"
#'
#' # Get group palette
#' group_colors <- AETHER_default_palette("groups")
#'
AETHER_default_palette <- function(type = "omics") {

  if (type == "omics") {
    return(c(
      ATACseq = "#2ca02c",   # Green
      CHIPseq = "#1f77b4",   # Blue
      RNAseq = "#ff7f0e",    # Orange
      CUTnTAG = "#9467bd",   # Purple
      CUTnRUN = "#8c564b",   # Brown
      Methylation = "#e377c2" # Pink
    ))
  } else if (type == "groups") {
    return(c(
      WT = "#1f77b4",        # Blue
      KO = "#d62728",        # Red
      Control = "#7f7f7f",   # Gray
      Treatment = "#9467bd", # Purple
      Treated = "#9467bd",   # Purple (alias)
      Untreated = "#7f7f7f", # Gray (alias)
      Healthy = "#2ca02c",   # Green
      Disease = "#d62728"    # Red
    ))
  } else {
    stop("type must be 'omics' or 'groups'")
  }
}
