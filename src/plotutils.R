library(scales)
library(latex2exp)
library(ggplot2)
library(reshape2)
library(ggpubr)
library(ggExtra)
library(grid)
library(gridExtra)
library(ggpubr)
library(gtable)
library(ggnewscale)
library(knitr)
ggplot2::theme_set(ggplot2::theme_bw())

crop <- function(filename) {
  knitr::plot_crop(filename)
}

default_dims <- function(dims = NULL) {
  if (is.null(dims)) {
    dims <- list()
  }
  if (is.null(dims$width)) {
    dims$width = 5
  }
  if (is.null(dims$height)) {
    dims$height = 5
  }
  return(dims)
}

save_plot <- function(p, filename, dims = NULL) {

  dims <- default_dims(dims)

  ggsave(plot = p, file = filename,
         device = "pdf",
         width = dims$width, height = dims$height,
         units = "in",
         dpi = 600)
  crop(filename)
}

add_if_not_null <- function(p, cond, value) {
  if (!is.null(cond)) {
    p <- p + value
  }
  return(p)
}