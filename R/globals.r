#' @importFrom fastmatch %fin%
#' @importFrom stats na.omit
#' @importFrom stats p.adjust
#' @importFrom stats phyper
#' @importFrom stats setNames
#' @import data.table
#' @import utils

utils::globalVariables(
  c(
    "genes", "pval", "batch", "ann_name", ".", "ann_name", "cluster",
    "extra_x_term", "group", "intra_x_term", "lfc_abs_lim", "min_signif",
    "name", "parent_Ngenes", "parent_genes", "parents", "qval_bh",
    "qval_bonferroni", "temp", "term", "type", "uniprot", "x_cluster"
  )
)
