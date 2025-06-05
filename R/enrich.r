#' @include fctBio.r
#' @include utils.r
#' @include globals.r
NULL

#' Simple Gene Set Enrichment
#'
#' This function performs enrichment analysis using topology-based annotations.
#' It applies the parent-child method for over-representation analysis.
#' Multi-testing correction is performed using BH and Bonferroni
#' methods for terms with a minimal p-value below 1e-4 by default. P-values
#' are calculated using the hypergeometric test. This function is intended
#' for internal use due to its reliance on integer IDs for genes and terms;
#' later a new implementation may allow direct use of this function.
#'
#' @param gene_input Input gene list converted as an INTEGER vector.
#' @param gene_universe Optional, Universe gene list as an INTEGER vector.
#' @param ann Modified annotation table as returned from \link{load_ann_space}
#' but with
#' @param lim_pmin Minimum value for pmin (lowest possible p-value).
#'   This limit is used for filtering before p-value correction.
#' @param classic Boolean, if TRUE performs analysis relative to the universe
#'   (vs. parent-child method).
#' @param log_level Logging level (see logging package). Default is WARN.
#' @param hard_pmin_filter Boolean, if TRUE, terms with pmin > lim_pmin are
#'   removed from outputs.
#' @export
enrich <- function(
  gene_input,
  ann,
  gene_universe = NULL,
  lim_pmin = 0.05,
  classic = FALSE,
  log_level = "WARN",
  hard_pmin_filter = TRUE
) {

  logging::basicConfig(level = log_level)

  logging::loginfo("Loading and filtering database based on type")

  assertthat::assert_that(
    is.integer(gene_input),
    all(! is.na(gene_input)),
    data.table::is.data.table(ann),
    all(purrr::map_lgl(ann$genes, is.integer)),
    all(purrr::map_lgl(ann$genes, function(x) all(! is.na(x)))),
    all(purrr::map_lgl(ann$parents, is.integer)),
    all(purrr::map_lgl(ann$parents, function(x) all(! is.na(x)))),
    is.integer(ann$term),
    all(! is.na(ann$term))
  )

  if (! is.null(gene_universe)) {
    assertthat::assert_that(
      is.integer(gene_universe),
      all(! is.na(gene_universe))
    )
  }

  lim_corr <- lim_pmin

  logging::logdebug("test Loading lists")

  allgenes <- stats::na.omit(unique(unlist(ann$genes)))

  if (is.null(gene_universe)) {
    universe <- allgenes
  } else {
    universe <- unique(gene_universe)
  }

  gi <- unique(gene_input)

  n_gi <- c(length(gi))
  n_universe <- c(length(universe))

  logging::loginfo("Filtering gene of interests based on universe")

  if (! is.null(gene_universe)) {
    universe <- intersect(universe, allgenes)
  }

  gi <- intersect(gi, universe)

  n <- length(universe) # Total number of genes considered
  k <- length(gi) # Number of gene of interest

  n_gi <- c(n_gi, k)
  n_universe <- c(n_universe, n)

  data_numbers <- data.frame(gi = n_gi, n_universe = n_universe)

  row.names(data_numbers) <- c("Input", "GI in Universe")

  logging::loginfo("Filtering database based on universe")

  gi <- sort(gi)
  universe <- sort(universe)

  ann[, genes := lapply(
    ann$genes,
    function(x) x[x %fin% universe]
  )]

  logging::loginfo(
    "Filtering database intesection with gene of interest based on universe"
  )

  ann$GI <- lapply(ann$genes, function(x) x[x %fin% gi])

  get_count_from_ids_dict <- c(n, k)
  names(get_count_from_ids_dict) <- c("genes", "GI")

  # Warning parent-child intersection Method !

  logging::loginfo("Denombring intersection at term level")

  ann$NGI <- vapply(ann$GI, length, USE.NAMES = FALSE, FUN.VALUE = 1L)
  ann$Ngenes <- vapply(ann$genes, length, USE.NAMES = FALSE,
    FUN.VALUE = 1L
  )

  logging::loginfo("Denombring intersection at parent level")

  if (classic) {
    ann$parent_NGI <- k
    ann$parent_Ngenes <- n
  } else {

    ann$parent_genes <- lapply(
      ann$parent_genes,
      function(x) {
        if (length(x) > 0L) {
          x[x %fin% universe]
        } else {
          universe
        }
      }
    )

    ann$parent_GI <- lapply(
      ann$parent_genes,
      function(x) x[x %fin% gi]
    )

    ann$parent_NGI <- vapply(ann$parent_GI, length, USE.NAMES = FALSE,
      FUN.VALUE = 1L
    )
    ann$parent_Ngenes <- vapply(ann$parent_genes, length, USE.NAMES = FALSE,
      FUN.VALUE = 1L
    )

    ann[parent_Ngenes == 0, parent_Ngenes := n] #nolint

    ann$parent_GI <- NULL
    ann$parent_genes <- NULL

  }

  get_pval_fisher <- function(n_genes, n_gi, parent_n_genes, parent_n_gi) {
    stats::phyper(
      n_gi - 1L,
      parent_n_gi,
      parent_n_genes - parent_n_gi,
      n_genes,
      lower.tail = FALSE,
      log.p = FALSE
    )
  }

  get_pmin <- function(n_genes, parent_n_genes, parent_n_gi) {
    get_pval_fisher(n_genes, min(n_genes, parent_n_gi),
      parent_n_genes, parent_n_gi
    )
  }

  logging::loginfo("Computing minimum probabilities and filtering out")

  ann[, pmin := vapply(
    seq_len(nrow(ann)),
    function(i) {
      get_pmin(
        ann[["Ngenes"]][i],
        ann[["parent_Ngenes"]][i],
        ann[["parent_NGI"]][i]
      )
    },
    FUN.VALUE = 1.0,
    USE.NAMES = FALSE
  )]

  logging::logdebug("Computing probabilities")

  ann[, pval := vapply(
    seq_len(nrow(ann)),
    function(i) {
      get_pval_fisher(
        ann[["Ngenes"]][i],
        ann[["NGI"]][i],
        ann[["parent_Ngenes"]][i],
        ann[["parent_NGI"]][i]
      )
    },
    FUN.VALUE = 1.0,
    USE.NAMES = FALSE
  )]

  ann[, `:=`(qval_bonferroni = as.numeric(NA), qval_bh = as.numeric(NA))]

  logging::logdebug("BH pvalue correction")

  ann[pmin < lim_corr, qval_bh := stats::p.adjust(pval, method = "BH")] # nolint
  logging::logdebug("Bonferoni pvalue correction")

  ann[pmin < lim_corr, qval_bonferroni := stats::p.adjust(pval, method = "bonferroni")] # nolint # nolint

  logging::logdebug("Sorting & output")

  # With data.table, you can sort by a column using setorder()
  data.table::setorder(ann, qval_bh) # nolint

  if (hard_pmin_filter) {
    return(invisible(
      list(enrich = ann[pmin < lim_corr], mapping = data_numbers,
        universe = universe
      )
    ))
  } else {
    return(invisible(
      list(enrich = ann, mapping = data_numbers, universe = universe)
    ))
  }

}
