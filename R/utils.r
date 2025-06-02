#' @include fctBio.r
#' @include globals.r
NULL

#' fctBio ggplot2 theme
#'
#' Simple theme
#' @export
THEME_NEXOMIS <- ggplot2::theme_classic() +
  ggplot2::theme(
    panel.background = ggplot2::element_rect(fill = "#FAFAFA"),
    axis.text.x = ggplot2::element_text(
      angle = 45.0, size = 10.0, vjust = 1.0, hjust = 1.0
    ),
    axis.text.y = ggplot2::element_text(
      angle = 0.0, size = 10.0, vjust = 1.0, hjust = 1.0
    ),
    strip.text = ggplot2::element_text(size = 12.0, color = "black"),
    legend.title = ggplot2::element_text(size = 12.0, color = "black"),
    legend.text = ggplot2::element_text(size = 10.0, color = "black"),
    plot.title = ggplot2::element_text(hjust = 0.5)
  )

DOWN_COLOR <- "#313695"
UP_COLOR <- "#A50026"
MID_COLOR <- "#FFFFBF"
NA_COLOR <- "#474747"

UP_TO_DOWN_GRADIENT_COLORS <- c(
  UP_COLOR, "#D73027", "#F46D43",
  "#FDAE61", "#FEE090", MID_COLOR,
  "#E0F3F8", "#ABD9E9", "#74ADD1",
  "#4575B4", DOWN_COLOR
)
