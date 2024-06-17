#' @include fctBio.r
#' @importFrom magrittr %>%
#' @import data.table
#' @importFrom fastmatch %fin%

NULL

#' ggplot2 theme
#'
#' Simple theme
#' @export
THEME_NEXOMIS <- ggplot2::theme_classic() + ggplot2::theme( # nolint
  panel.background = ggplot2::element_rect(fill = "#FAFAFA"),
  axis.text.x = ggplot2::element_text(
    angle = 45, size = 10, vjust = 1, hjust = 1),
  axis.text.y = ggplot2::element_text(
    angle = 0, size = 10, vjust = 1, hjust = 1),
  strip.text = ggplot2::element_text(size = 12, color = "black"),
  legend.title = ggplot2::element_text(size = 12, color = "black"),
  legend.text = ggplot2::element_text(size = 10, color = "black"),
  plot.title = ggplot2::element_text(hjust = 0.5))

DOWN_COLOR <- "#313695" # nolint
UP_COLOR <- "#A50026" # nolint
MID_COLOR <- "#FFFFBF" # nolint
NA_COLOR <- "#474747" # nolint

UP_TO_DOWN_GRADIENT_COLORS <- c( # nolint
  UP_COLOR, "#D73027", "#F46D43",
  "#FDAE61", "#FEE090", MID_COLOR,
  "#E0F3F8", "#ABD9E9", "#74ADD1",
  "#4575B4", DOWN_COLOR)
