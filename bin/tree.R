#!/usr/local/bin/Rscript --vanilla

library(ggtree)
library(ape)

branch_scale_factor <- 0.05

strip_extension <- function(x) sub("\\.[^.]+$", "", x)
read_cluster_tips <- function(path, tree_tips) {
  lines <- trimws(readLines(path, warn = FALSE))
  lines <- lines[nzchar(lines)]

  clusters <- lapply(lines, function(line) {
    tips <- regmatches(line, gregexpr("[^[:space:]]+\\.fna", line, perl = TRUE))[[1]]
    unique(strip_extension(tips))
  })

  duplicated_tips <- unique(unlist(clusters)[duplicated(unlist(clusters))])
  if (length(duplicated_tips) > 0) {
    stop(
      "Some tips were assigned to multiple clusters in ",
      path,
      ": ",
      paste(duplicated_tips, collapse = ", ")
    )
  }

  clusters <- lapply(clusters, intersect, y = tree_tips)
  Filter(function(tips) length(tips) > 1, clusters)
}

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) {
  stop("Usage: tree.R <tree_file> <cluster_file> [output_png]")
}

tree_file <- args[1]
annotation_file <- args[2]
output_file <- if (length(args) >= 3) args[3] else "tree.png"

tree <- read.tree(tree_file)
tree$tip.label <- strip_extension(tree$tip.label)

cluster_tips <- read_cluster_tips(annotation_file, tree$tip.label)
cluster_marker_nodes <- if (length(cluster_tips) > 0) {
  vapply(cluster_tips, getMRCA, integer(1), phy = tree)
} else {
  integer(0)
}

plot_tree <- tree
tree_depth <- NULL
display_xmax <- NULL
scale_bar_width <- NULL
scale_bar_label <- NULL

if (!is.null(tree$edge.length)) {
  tree_depth <- max(node.depth.edgelength(tree))
  plot_tree$edge.length <- tree$edge.length * branch_scale_factor
  display_xmax <- max(node.depth.edgelength(plot_tree)) + (tree_depth * 0.15)

  pretty_widths <- pretty(c(0, tree_depth), n = 4)
  pretty_widths <- pretty_widths[pretty_widths > 0]

  if (length(pretty_widths) > 0) {
    scale_bar_width <- pretty_widths[1] * branch_scale_factor
    scale_bar_label <- format(signif(pretty_widths[1], 3), scientific = FALSE)
  }
}

tip_label_offset <- if (is.null(tree_depth)) 0 else tree_depth * 0.01

p <- ggtree(plot_tree) +
  geom_tiplab(
    aes(label = label),
    align = TRUE,
    linetype = "dotted",
    linesize = 0.2,
    offset = tip_label_offset,
    size = 3.3
  )

if (length(cluster_marker_nodes) > 0) {
  cluster_markers <- subset(p$data, node %in% cluster_marker_nodes, c("x", "y"))
  p <- p + ggplot2::geom_point(
    data = cluster_markers,
    mapping = ggplot2::aes(x = x, y = y),
    shape = 17,
    size = 3.9,
    color = "red3",
    inherit.aes = FALSE
  )
}

if (!is.null(display_xmax)) {
  p <- p + ggplot2::coord_cartesian(xlim = c(0, display_xmax), clip = "off")
}

if (!is.null(tree$node.label) && any(nzchar(tree$node.label))) {
  p <- p + geom_text2(
    aes(subset = !isTip & !is.na(label) & nzchar(label), label = label),
    hjust = -0.25,
    size = 3
  )
}

if (!is.null(scale_bar_width) && !is.null(scale_bar_label)) {
  scale_y <- -0.7

  p <- p +
    geom_treescale(
      x = 0,
      y = scale_y,
      width = scale_bar_width,
      label = scale_bar_label,
      offset = scale_bar_width * 0.2,
      color = "black"
    ) +
    ggplot2::annotate(
      "point",
      x = scale_bar_width * 2.6,
      y = scale_y,
      shape = 17,
      size = 3.9,
      color = "red3"
    ) +
    ggplot2::annotate(
      "text",
      x = scale_bar_width * 3.0,
      y = scale_y,
      label = "SNP Cluster",
      hjust = 0,
      vjust = 0.5,
      size = 3
    )
}

plot_width <- max(6, min(12, 3 + max(nchar(tree$tip.label)) * 0.09))
plot_height <- max(3.5, min(10, 1.4 + Ntip(tree) * 0.55))
right_margin <- 10 + max(nchar(tree$tip.label)) * 1.2

p <- p + ggplot2::theme(
  plot.margin = ggplot2::margin(5.5, right_margin, 30, 5.5),
  legend.position = "none"
)

ggsave(output_file, plot = p, width = plot_width, height = plot_height, units = "in", dpi = 300, limitsize = FALSE)
