#' @importFrom fastmatch %fin%
#' @importFrom magrittr %>%
#' @importFrom data.table :=
#' @importFrom data.table .SD
#' @importFrom rlang .data

utils::globalVariables(
  c(
    "genes", "pval", "batch", "ann_name", ".", "ann_name", "cluster",
    "extra_x_term", "group", "intra_x_term", "lfc_abs_lim", "min_signif",
    "name", "parent_Ngenes", "parent_genes", "parents", "qval_bh",
    "qval_bonferroni", "temp", "term", "type", "uniprot", "x_cluster"
  )
)
