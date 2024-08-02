#' @include utils.r

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
#' @param regex Regex pattern to remove from gene name from input lists
#'   (useful for deleting accession version info which should not be present
#'   on database files).
#' @param type_keep Character vector of term types to retain (all types
#'   are kept by default); e.g., Biological_process for gene ontology.
#' @param classic Boolean, if TRUE performs analysis relative to the universe
#'   (vs. parent-child method).
#' @param log_level Logging level (see logging package). Default is WARN.
#' @export
enrich <- function(gene_input, ann, gene_universe = NULL, lim_pmin = 0.05,
  classic = FALSE, log_level = "WARN", hard_pmin_filter = TRUE) {

  logging::basicConfig(level = log_level)

  logging::loginfo("Loading and filtering database based on type")

  assertthat::assert_that(is.integer(gene_input))
  assertthat::assert_that(all(! is.na(gene_input)))
  if (! is.null(gene_universe)) {
    assertthat::assert_that(is.integer(gene_universe))
    assertthat::assert_that(all(! is.na(gene_universe)))
  }
  assertthat::assert_that(data.table::is.data.table(ann))
  assertthat::assert_that(
    all(purrr::map_lgl(ann$genes, is.integer)))
  assertthat::assert_that(
    all(purrr::map_lgl(ann$genes, function(x) all(! is.na(x)))))
  assertthat::assert_that(
    all(purrr::map_lgl(ann$parents, is.integer)))
  assertthat::assert_that(
    all(purrr::map_lgl(ann$parents, function(x) all(! is.na(x)))))
  assertthat::assert_that(is.integer(ann$term))
  assertthat::assert_that(all(! is.na(ann$term)))

  lim_corr <- lim_pmin

  logging::logdebug("test Loading lists")

  allgenes <- na.omit(unique(unlist(ann$genes)))

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

  ann[, genes := lapply( #nolint
    ann$genes,
    function(x) x[x %fin% universe]
  )]

  logging::loginfo(
    "Filtering database intesection with gene of interest based on universe")

  ann$GI <- lapply(ann$genes, function(x) x[x %fin% gi])

  get_count_from_ids_dict <- c(n, k)
  names(get_count_from_ids_dict) <- c("genes", "GI")

  # Warning parent-child intersection Method !

  logging::loginfo("Denombring intersection at term level")

  ann$NGI <- vapply(ann$GI, length, USE.NAMES = FALSE, FUN.VALUE = integer(1))
  ann$Ngenes <- vapply(ann$genes, length, USE.NAMES = FALSE,
    FUN.VALUE = integer(1))

  logging::loginfo("Denombring intersection at parent level")

  if (classic) {
    ann$parent_NGI <- k
    ann$parent_Ngenes <- n
  } else {

    ann$parent_genes <- lapply(
      ann$parent_genes,
      function(x) {
        if (length(x) > 0) {
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
      FUN.VALUE = integer(1))
    ann$parent_Ngenes <- vapply(ann$parent_genes, length, USE.NAMES = FALSE,
      FUN.VALUE = integer(1))

    ann[parent_Ngenes == 0, parent_Ngenes := n] #nolint

    ann$parent_GI <- NULL
    ann$parent_genes <- NULL

  }

  get_pval_fisher <- function(n_genes, n_gi, parent_n_genes, parent_n_gi) {
    phyper(
      n_gi - 1,
      parent_n_gi,
      parent_n_genes - parent_n_gi,
      n_genes,
      lower.tail = FALSE,
      log.p = FALSE
    )
  }

  get_pmin <- function(n_genes, parent_n_genes, parent_n_gi) {
    get_pval_fisher(n_genes, min(n_genes, parent_n_gi),
      parent_n_genes, parent_n_gi)
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
    FUN.VALUE = numeric(1),
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
    FUN.VALUE = numeric(1),
    USE.NAMES = FALSE
  )]

  ann[, `:=`(qval_bonferroni = as.numeric(NA), qval_bh = as.numeric(NA))]

  logging::logdebug("BH pvalue correction")

  ann[pmin < lim_corr, qval_bh := p.adjust(pval, method = "BH")] # nolint
  logging::logdebug("Bonferoni pvalue correction")

  ann[pmin < lim_corr, qval_bonferroni := p.adjust(pval, method = "bonferroni")] # nolint # nolint

  logging::logdebug("Sorting & output")

  # With data.table, you can sort by a column using setorder()
  data.table::setorder(ann, qval_bh) # nolint

  if (hard_pmin_filter) {
    return(invisible(
    list(enrich = ann[pmin < lim_corr], mapping = data_numbers,
      universe = universe)
  ))
  } else {
    return(invisible(
    list(enrich = ann, mapping = data_numbers, universe = universe)
  ))
  }

}

#' Nested Enrichment Analysis Class
#'
#' @description
#' The `NestedEnrich` class automates the process of performing enrichment analysis
#' on nested data frames where each row represents a distinct set of conditions
#' under which gene lists are analyzed. This class is particularly useful for
#' complex experimental designs where gene data is segmented into multiple batches,
#' groups, or types, and where each segment might require separate enrichment analysis
#' against different annotation backgrounds.
#'
#' @param in_batch Vector of batch codes; used to select or report only gene lists
#' or enrichment results for specified batches.
#' @param in_group Vector of group codes; used to select or report only gene lists
#' or enrichment results for specified groups.
#' @param in_type Vector of type codes; used to select or report only gene lists
#' or enrichment results for specified types.
#' @param in_ann_name Vector of annotation name codes; used to select or report
#' only enrichment results for specified annotations.
#' @param p_type Type of p-value to consider (e.g., raw p-value, adjusted p-value).
#' @param max_cluster Maximum number of clusters to display or analyze.
#' @param max_term_per_cluster Maximum number of terms per cluster to display or analyze.
#' @export

NestedEnrich <- R6::R6Class("NestedEnrich", # nolint
  public <- list(
    #' @description 
    #' Initialize `NestedEnrich` object.
    #' 
    #' The constructor method for `NestedEnrich` sets up the enrichment analysis
    #' environment by loading and preparing data. It handles the initialization
    #' of data structures required for nested enrichment analysis based on
    #' specified parameters.
    #' 
    #' @param nested_df A nested tibble containing gene sets organized by batch,
    #' group, and type. Each row should correspond to a unique analysis context,
    #' containing a nested dataframe of gene identifiers.
    #'   - `batch`: First level of grouping, specifies the batch context of the gene list.
    #'   - `group`: Second level of grouping, specifies the group context.
    #'   - `type`: Third level of grouping, specifies the type context.
    #'   - `data`: Column name containing the nested dataframe with gene identifiers;
    #'     must match the `data_col` argument.
    #'   - `data_univ`: Optional; column name for the nested dataframe with the 
    #'     "universe" of genes against which enrichment is assessed; 
    #      must match the `data_univ_col` argument.
    #' @param ann_space Annotation space data loaded from `load_ann_space`, defining
    #' the annotation context for enrichment analysis.
    #' @param batch_col Optional; specifies the column name for batch if different from 'batch'.
    #' @param group_col Optional; specifies the column name for group if different from 'group'.
    #' @param type_col Optional; specifies the column name for type.
    #' @param batch_labels Optional; specifies labels for batches if renaming is required.
    #' @param group_labels Optional; specifies labels for groups if renaming is required.
    #' @param data_col Name of the column in `nested_df` containing gene data.
    #' @param data_univ_col Optional; name of the column containing the universe of genes.
    #' @param data_nested_id Column identifier in the nested data frames for gene identifiers.
    #' @param regex Regular expression to clean gene identifiers (e.g., to remove version numbers).
    #' @param lim_pmin Minimum p-value threshold; gene sets below this threshold are not retained
    #' unless `hard_pmin_filter` is set to FALSE.
    #' @param log_level Logging level as defined by the logging package.
    #' @param do_cluster_analysis Boolean; if FALSE, skips the automatic clustering analysis.
    #' @param hard_pmin_filter Boolean; if TRUE, removes gene sets with p-values below `lim_pmin`.
    #' @param classic Boolean; if TRUE, performs analysis relative to the universe of genes.
    #' @param ... Additional arguments passed to the `enrich` function.
    #' @return A new `NestedEnrich` object initialized with the provided data and settings.
    initialize = function(nested_df, ann_space, batch_col = NULL,
      group_col = NULL, type_col = NULL, batch_labels = NULL, 
      group_labels = NULL, data_col = "data", data_univ_col = NULL,
      data_nested_id = "uniprot", regex = "-.*", lim_pmin = 0.05,
      classic = FALSE, log_level = "WARN", do_cluster_analysis = TRUE,
      hard_pmin_filter = TRUE, ...) {

      logging::basicConfig(level = log_level)
      logging::loginfo("Initialize cluster")
      logging::loginfo("Parse inputs")
      private$data_nested_id <- data_nested_id
      nested_df <- tibble::as_tibble(nested_df)
      if (! is.null(batch_col)) {
        nested_df <- dplyr::rename_with(
          nested_df,
          function(x) ("batch"),
          tidyselect::all_of(batch_col)
        )
      }
      if (! is.null(group_col)) {
        nested_df <- dplyr::rename_with(
          nested_df,
          function(x) ("group"),
          tidyselect::all_of(group_col)
        )
      }
      if (! is.null(type_col)) {
        nested_df <- dplyr::rename_with(
          nested_df,
          function(x) ("type"),
          tidyselect::all_of(type_col)
        )
      }
      if (! ("batch" %in% names(nested_df) & "group" %in% names(nested_df))) {
        logging::logerror("batch and/or group undefined")
        stop()
      }
      private$batch_labels <- unique(nested_df$batch)
      names(private$batch_labels) <- private$batch_labels

      if (! is.null(batch_labels)) {
        private$batch_labels[names(batch_labels)] <- batch_labels
      }

      private$group_labels <- unique(nested_df$group)
      names(private$group_labels) <- private$group_labels

      if (is.null(group_labels)) {
        private$group_labels[names(group_labels)] <- group_labels
      }
      nested_dt <- data.table::data.table(
        batch = as.factor(nested_df$batch),
        group = as.factor(nested_df$group),
        type = as.factor(nested_df$type),
        gene_inputs = lapply(
          nested_df[[data_col]],
          function(df) {
            df[[data_nested_id]]
          }
        ),
        gene_universe =
          if (is.null(data_univ_col)) {
            NULL
          } else {
            lapply(
              nested_df[[data_univ_col]],
              function(df) {
                df[[data_nested_id]]
              }
            )
          }
      )
      logging::loginfo("Copy ann_space")
      private$raw_ann <- data.table::copy(ann_space)
      logging::loginfo("Parse genes with regex")
      if (regex != "") {
        nested_dt$gene_inputs <- lapply(nested_dt$gene_inputs, function(x) {
          unique(gsub(regex, "", x, perl = TRUE))
        })
        if (!is.null(data_univ_col)) {
          nested_dt$gene_universe <-
            lapply(nested_dt$gene_universe, function(x) {
              unique(gsub(regex, "", x, perl = TRUE))
        })
        }
      }

      logging::loginfo("Build gene IDs")
      all_genes_raw_ann <- unique(unlist(private$raw_ann[["genes"]]))
      all_genes_inputs <- unique(unlist(nested_dt$gene_inputs))
      if (! is.null(data_univ_col)) {
        all_genes_inputs <- unique(c(all_genes_inputs,
          unique(unlist(nested_dt$gene_universe))
        ))
      }
      # get gene names to ID
      private$gene_ids <-
        sort(na.omit(unique(c(all_genes_raw_ann, all_genes_inputs))))


      logging::loginfo("match gene ids to integer")
      private$raw_ann$genes <- lapply(
        private$raw_ann$genes,
        function(x) {
          fastmatch::fmatch(na.omit(x), private$gene_ids)
        }
      )
      if (! classic) {
        private$raw_ann$parent_genes <- lapply(
          private$raw_ann$parent_genes,
          function(x) {
            fastmatch::fmatch(na.omit(x), private$gene_ids)
          }
        )
      }
      nested_dt$gene_inputs <- lapply(nested_dt$gene_inputs, function(x) {
        fastmatch::fmatch(na.omit(x), private$gene_ids)
      })
      if (! is.null(data_univ_col)) {
        nested_dt$gene_universe <- lapply(nested_dt$gene_universe, function(x) {
          fastmatch::fmatch(na.omit(x), private$gene_ids)
        })
      }

      logging::loginfo("term features to index")
      # must not be sorted to be able to retrieve data
      private$term_ids <- private$raw_ann$term
      # must not be sorted to be able to retrieve data
      private$term_names <- private$raw_ann$name
      private$term_types <- private$raw_ann$type
      private$raw_ann$parents <- lapply(private$raw_ann$parents, function(x) {
          as.integer(
            na.omit(fastmatch::fmatch(x, private$term_ids))
          )
        }
      )
      private$raw_ann$term <- fastmatch::fmatch(
        private$raw_ann$term, private$term_ids)

      logging::loginfo("filter raw_ann")
      private$raw_ann <-
        private$raw_ann[, c("term", "genes", "parent_genes", "parents",
          "ann_name"), with = FALSE]

      logging::loginfo("prepare combinations")
      combinations <- expand.grid(
        i = seq_len(nrow(nested_dt)),
        ann_name = levels(ann_space$ann_name)
      )
      combinations_list <- split(combinations, seq_len(nrow(combinations)))
      raw_ann <- private$raw_ann

      logging::loginfo("Enrichment loop")

      fx <- function(i, ann_name_input) {
        nested_dt_row <- as.data.table(nested_dt[i, , drop = FALSE])
        nested_dt_row$ann_name <- ann_name_input
        ann <- raw_ann[ann_name == ann_name_input, ]
        ann[, ann_name := NULL]
        gene_input <- nested_dt_row$gene_inputs[[1]]
        if (is.null(data_univ_col)) {
          gene_universe <- NULL
        } else {
          gene_universe <- nested_dt_row$gene_universe[[1]]
        }
        res <- enrich(gene_input, ann, gene_universe = gene_universe,
          classic = classic, lim_pmin = lim_pmin, log_level = log_level,
          hard_pmin_filter = hard_pmin_filter)
        nested_dt_row$mapping <- list((res$mapping))
        nested_dt_row$enrich <- list(res$enrich)
        nested_dt_row$universe <- list(res$universe)
        nested_dt_row
      }

      private$results <- rbindlist(
        lapply(
          combinations_list,
          function(x) (fx(x[["i"]],x[["ann_name"]]))
        ), fill = TRUE
      )
      logging::loginfo("Enrichment terminated.")

      if (do_cluster_analysis) {
        logging::loginfo("Start clustering analysis")
        self$filter_and_set_significant_results()
        logging::loginfo("Finished clustering analysis")
      } else {
        self$filter_and_set_significant_results(
          0.05, p_type = "qval_bonferroni", min_signif_term_for_clust = Inf)
      }
    },

    #' @description
    #' set batch labels
    #' @param batch_labels batch_lables (named vector ordered)
    set_batch_labels = function(batch_labels) {
      private$batch_labels <- batch_labels
    },

    #' @description
    #' set group labels
    #' @param group_labels group_labels (named vector ordered)
    set_group_labels = function(group_labels) {
      private$group_labels <- group_labels
    },

    #' @description
    #' Return results
    #' @param verbose Replace integer with IDs and add names for terms
    #' @return Nested data frame with results
    get_results = function(verbose = FALSE) {
      results <- data.table::copy(private$results)
      if (verbose) {
        results[, gene_inputs := lapply(gene_inputs, function(x) {
          private$gene_ids[x]
        })]
        results[, universe := lapply(universe, function(x) {
          private$gene_ids[x]
        })]
        results[, enrich := lapply(enrich, function(dt) {
          new_dt <- data.table::copy(dt)
          new_dt[, genes := lapply(genes, function(x) {
            private$gene_ids[x]
          })]
          new_dt[, GI := lapply(GI, function(x) {
            private$gene_ids[x]
          })]
          new_dt[, name := sapply(term, function(x) {
            private$term_names[x]
          }, USE.NAMES = FALSE)]
          new_dt[, type := sapply(term, function(x) {
            private$term_types[x]
          }, USE.NAMES = FALSE)]
          new_dt[, term := sapply(term, function(x) {
            private$term_ids[x]
          }, USE.NAMES = FALSE)]
          new_dt[, parents := lapply(parents, function(x) {
            private$term_ids[x]
          })]
          data.table::setcolorder(new_dt,
            c("term", "name", "type", setdiff(names(new_dt), c("term", "name", "type"))))
          return(new_dt)
        })]
      }
      return(results)
    },

    #' @description
    #' Return result table after filtering
    #' @param verbose Replace integer with IDs and add names for terms
    #' @return Nested data frame with results
    filter_and_get_results = function(
      in_batch, in_group, in_type, in_ann_name, verbose = FALSE) {
      self$get_results(verbose)[
        batch %fin% in_batch &
        group %fin% in_group &
        type %fin% in_type &
        ann_name %fin% in_ann_name
      ]
    },

    #' @description
    #' Return results unnested
    #' @param verbose Replace integer with IDs and add names for terms
    #' @return Unnested data frame with results
    unnest_and_get_results = function(verbose = FALSE) {
      results <- self$get_results(verbose)
      data.table::rbindlist(lapply(
        seq_len(nrow(results)),
        function(i) {
          result <- data.table::copy(results[["enrich"]][[i]])
          result$batch <- results[["batch"]][i]
          result$group <- results[["group"]][i]
          result$type <- results[["type"]][i]
          result$ann_name <- results[["ann_name"]][i]
          result
        }
      ))
    },

    #' @description
    #' filter results based on pvalue and set significant results
    #' @param p_max_enrich max p-value to allow
    #' @param p_type type of pvalue to use:
    #' - "pval" : raw p-value
    #' - "qval_bh" : q-value using Benjamini Hochberg
    #' - "qval_bonferroni": q-value using Bonferoni
    #' @param build_and_set_i_matrix build a set term-gene incidence matrix
    #' @param build_and_set_p_matrix build a set term-sample pvalue matrix
    #' @param build_and_set_hclust built and set hclust with default parameters
    #' @param min_signif_term_for_clust minimun number of term to start clusters
    filter_and_set_significant_results = function(p_max_enrich = 0.05,
      p_type = "qval_bh", build_and_set_i_matrix = TRUE,
      build_and_set_p_matrix = FALSE, build_and_set_hclust = TRUE,
      min_signif_term_for_clust = 10) {

      if (! p_type %in% c("pval", "qval_bh", "qval_bonferroni")) {
        logging::logerror("`p_type` must be pval, pval_bh or pval_bonferroni")
        stop()
      }

      # reset
      private$i_matrix <- NULL
      private$p_matrix <- NULL
      private$hclust <- NULL
      private$clusters <- NULL

      private$significant_results <- self$unnest_and_get_results()

      private$significant_results[,
        significant := fifelse(
          is.na(.SD[[p_type]]) | .SD[[p_type]] > p_max_enrich, FALSE, TRUE)]

      private$significant_terms <- unique(
        private$significant_results$term[
          private$significant_results$significant]
      )

      private$significant_results <- private$significant_results[
        term %fin% private$significant_terms
      ]

      if (length(private$significant_terms) <= min_signif_term_for_clust) {
        if (is.infinite(min_signif_term_for_clust)) {
          logging::loginfo("Clustering skipped")
        } else {
          logging::logwarn("Cluster not done (signif term < threshold)")
        }
        if (length(private$significant_terms) > 0) {
          private$clusters <- rep(1L, length(private$significant_terms))
          names(private$clusters) <- private$term_ids[private$significant_terms]
        } else {
          logging::logwarn("No significant term with used threshold")
        }
      } else {
        if (build_and_set_i_matrix) {
          self$build_and_set_i_matrix()
        }
        if (build_and_set_p_matrix) {
          self$build_and_set_p_matrix()
        }
        if (build_and_set_hclust) {
          if (! build_and_set_i_matrix) {
            logging::logwarn(paste(
              "build_and_set_hclust with default parameters does require the",
              "build_and_set_i_matrix option to be TRUE"
            ))
          } else {
            self$build_and_set_hclust()
          }
        }
      }
    },

    #' @description
    #' Get the significant terms
    get_significant_terms = function() {
      private$term_ids[private$significant_terms]
    },

    #' @description
    #' get significant results
    #' @return tibble with significant results only
    get_significant_results = function() {
      results <- data.table::copy(private$significant_results)
      results[, `:=`(
        term = private$term_ids[.SD[["term"]]],
        name = private$term_names[.SD[["term"]]],
        ann_type = private$term_types[.SD[["term"]]],
        genes = lapply(
          .SD[["genes"]],
          function(x) private$gene_ids[x]
        ),
        GI = lapply(
          .SD[["GI"]],
          function(x) private$gene_ids[x]
        )
      )]
      results
    },

    #' @description
    #' Build and set the incidence matrix from the terms present in significant
    #' results.
    #'
    #' The raw annotation database is taken into account to build up this matrix
    #' with the union of all universes as base.
    #'
    #' This function requires the filter_and_set_significant_results method to
    #' be ran before
    build_and_set_i_matrix = function() {

      if (is.null(private$significant_results)) {
        logging::logerror(paste(
          "This function requires to run first",
          "the filter_and_set_significant_results method"))
        stop()
      }

      logging::logdebug("start incidence matrix")

      filtered_raw_ann <- private$raw_ann[term %fin% private$significant_terms]

      filtered_raw_ann[, `:=`(
        term = private$term_ids[.SD[["term"]]],
        genes = lapply(
          .SD[["genes"]],
          function(x) private$gene_ids[x]
        )
      )]

      allgenes <- unique(unlist(filtered_raw_ann$genes))

      private$i_matrix <- do.call(rbind, lapply(
        filtered_raw_ann$term,
        function(in_term) {
          matrix(
            as.integer(allgenes %fin%
              filtered_raw_ann[term == in_term, genes[[1]]]),
            nrow = 1,
            dimnames = list(in_term, NULL)
          )
        }
      ))

      colnames(private$i_matrix) <- allgenes
    },

    #' @description
    #' Build and set the pval matrix from the terms present in significant
    #' results.
    #'
    #' This function requires the filter_and_set_significant_results method to
    #' be ran before
    #' @param p_type type of pvalue to use:
    #' - "pval" : raw p-value
    #' - "qval_bh" : q-value using Benjamini Hochberg
    #' - "qval_bonferroni": q-value using Bonferoni
    build_and_set_p_matrix = function(p_type = "pval") {
      if (is.null(private$significant_results)) {
        logging::logerror(paste(
          "This function requires to run first",
          "the filter_and_set_significant_results method"))
        stop()
      }

      dt <- self$get_significant_results(
        )[, .(term, batch, group, type, get(p_type))]
      setnames(dt, names(dt)[5], p_type)
      dt[, col := paste(batch, group, type, sep = "_")]
      dt_wide <- dcast(dt, term ~ col, value.var = p_type)

      rn <- dt_wide$term
      dt_wide[, term := NULL]
      private$p_matrix <- as.matrix(dt_wide)
      rownames(private$p_matrix) <- rn
    },

    #' @description
    #' get incidence matrix
    get_i_matrix = function() {
      private$i_matrix
    },

    #' @description
    #' get incidence matrix
    get_p_matrix = function() {
      private$p_matrix
    },

    #' @description
    #' Make hierarchical classification using pvclust.
    #'
    #' IMPORTANT: Cluster ar automatically assigned to 8 clusters if not using
    #' pvclust or using a p-value threshold of 0.8 in case of using pvclust.
    #' Methods are available to modify these defaults.
    #' @param method_dist pvclust param preset (based on philentropy::distance
    #' function)
    #' @param method_hclust pvclust param preset
    #' @param using_pvclust whether to use pvclust toi build the cluster
    #' @param matrix_type type of matrix to use for classification
    #' - "incidence" : Term gene incidence
    #' - "pval" : Term List p-value
    #' - "both" : both list concatenated
    #' @param dim_reduce number of dimension to reduce the data (0 or False
    #' means no reduction). If < 1, then the tol prcomp args will be used.
    #' else if > 1 the rank. arg will be used to define the number of PCs to
    #' keep.
    #' @param prcomp_args prcomp args that are passed if dim_reduce is not
    #' 0/False
    #' @param ... additional parameters for pvclust or hclust
    #' @return NULL
    build_and_set_hclust = function(method_dist = "cosine", #nolint
      method_hclust = "ward.D2", using_pvclust = FALSE, dim_reduce = 6,
      matrix_type = "incidence", prcomp_args = list() ,...) { #nolint
      private$using_pvclust <- using_pvclust

      # build matrix
      if (matrix_type == "incidence") {
        mat <- private$i_matrix
      } else if (matrix_type == "pval") {
        mat <- private$p_matrix
      } else if (matrix_type == "both") {
        mat <- cbind(private$i_matrix, private$p_matrix)
      } else {
        logging::logerror("Incorrect value for matrix_type")
        stop()
      }

      if (dim_reduce) {
        prcomp_args$x <- mat
        if (dim_reduce < 1) {
          prcomp_args$tol <- dim_reduce
        } else {
          prcomp_args$rank. <- as.integer(dim_reduce)
        }
        pca <- do.call(prcomp, prcomp_args)
        mat <- pca$x
      }

      if (using_pvclust) {
        private$hclust <- pvclust::pvclust(
          t(mat),
          method.dist = function(x) {
            philentropy::distance(
              t(x), method = method_dist,
              as.dist.obj = TRUE,
              mute.message = TRUE,
              use.row.names = TRUE
            )
          },
          method.hclust = method_hclust,
          ...
        )
      } else {
        tdist <- philentropy::distance(
          mat,
          method = method_dist,
          as.dist.obj = TRUE,
          mute.message = TRUE,
          use.row.names = TRUE
        )
        private$hclust <- hclust(tdist, method = method_hclust)
      }
      if (using_pvclust) {
        self$cut_hclust_and_set_clusters(0.8)
      } else {
        self$cut_hclust_and_set_clusters(8)
      }
    },

    #' @description
    #' get hierarchical cluster
    get_hclust = function() {
      private$hclust
    },

    #' @description
    #' get clusters
    get_clusters = function() {
      private$clusters
    },

    #' @description
    #' plot hierachical clustering of terms
    #' @param cluster_rect Whether to add or not cluster rectangle
    #' @param rect_linetype linetype for cluster rectangle
    #' @param rect_linewidth size for cluster rectangle
    #' @param ggdendro_args list of arguments passed to ggdendrogram
    #' @param ... passed to get_label_dict
    #' @return graphical output from ggdendro::ggdendrogram
    plot_hclust = function(cluster_rect = TRUE, rect_linetype = "solid",
      rect_linewidth = 1, ggdendro_args = list(), ...) {
      if (private$using_pvclust) {
        hc <- private$hclust$hclust
      } else {
        hc <- private$hclust
      }
      if (cluster_rect) {
        dendro <- as.dendrogram(hc)
      }
      term2label <- self$get_label_dict(...)
      hc$labels <- term2label[hc$labels]
      defaults <- modifyList(list(rotate = TRUE), ggdendro_args)
      args <- modifyList(defaults, list(data = hc, ...))
      dend <- do.call(ggdendro::ggdendrogram, args)

      if (length(unique(private$clusters[duplicated(private$clusters)])) == 0) {
        cluster_rect <- FALSE
      }

      if (cluster_rect) {
        max_height <- max(dendextend::get_nodes_attr(dendro, "height"))
        rect_data <- purrr::map_dfr(
          unique(private$clusters[duplicated(private$clusters)]),
          function(cluster_id) {
            cluster_members <-
              names(private$clusters)[private$clusters == cluster_id]
            cluster_leaves <- match(cluster_members, labels(dendro))
            height <- max(
              dendextend::get_nodes_attr(
                dendextend::prune(
                  dendro,
                  labels(dendro)[!labels(dendro) %in% cluster_members]
                ),
                "height"
              )
            )
            data.frame(
              xmin = min(cluster_leaves) - 0.3,
              xmax = max(cluster_leaves) + 0.3,
              ymin = - max_height * 0.01,
              ymax = height - (max_height * 0.01),
              cluster = cluster_id
            )
          }
        )

        rect_data <- dplyr::arrange(rect_data, .data$xmin)
        rect_data$color <- "a"
        rect_data$color[seq(2, nrow(rect_data), 2)] <- "b"

        # Plot the dendrogram
        dend <- dend +
          ggplot2::geom_rect(
            data = rect_data,
            ggplot2::aes(
              xmin = .data$xmin, xmax = .data$xmax,
              ymin = .data$ymin, ymax = .data$ymax,
              color = .data$color),
            linetype = rect_linetype,
            linewidth = rect_linewidth,
            fill = NA
          )  +
        ggplot2::guides(color = "none") +
        ggplot2::geom_label(data = rect_data, 
          mapping = ggplot2::aes((.data$xmax + .data$xmin)/2, 
            .data$ymax, label = .data$cluster),
        fill="#FFFFFFD0")
      }
      return(dend)
    },

    #' @description
    #' Cut hclust and set cluster
    #' cut hclust return a named integer vector. Each term is linked with its
    #' cluster id. Finally, the orphan terms are given an additional id.
    #' @param value value to cut the tree
    #'   - p-values for pvclust (see pvclust::pvpick)
    #'   - number of groups for hclust
    #' @param min_size for pvclust only, the minimum size of clusters
    #' @param max_size for pvclust only, the minimum size of clusters
    #' @param rm_redundancy_method Method
    #' - "largest" : select the largest group for a term to be in only one
    #' cluster
    #' - "smallest" : select the smallest group for a term to be in only one
    #' cluster
    #' @param max_size_only_to_parents whether max_size should apply only to
    #' cluster that are parent (include other clusters)
    #' @param min_size_only_to_children whether min_size should apply only to
    #' cluster that are child (are included by other clusters)
    #' @param filter_max_size_first whether to filter based on max_size before
    #' min_size. The order might have an impact on clusters kept.
    #' @param pvclust_pv probability value (pv) for pvclust
    #' (see pvclust::pvpick)
    #' @return named character as detailed in description
    cut_hclust_and_set_clusters = function(value, min_size = 5, max_size = 8,
      max_size_only_to_parents = TRUE, min_size_only_to_children = TRUE,
      filter_max_size_first = TRUE, rm_redundancy_method = "largest",
      pvclust_pv = "au") {
      if (value < 0) {
        logging::logdebug("`value` cannot be < 0")
        stop()
      }
      if (! rm_redundancy_method %in% c("largest", "smallest")) {
        logging::logdebug("`pvclust_unique_method` wrong value")
        stop()
      }
      if (! pvclust_pv %in% c("si", "au", "bp")) {
        logging::logdebug("`pvclust_pv` wrong value")
        stop()
      }
      if (private$using_pvclust) {
        if (value > 1) {
          logging::logdebug("`value` cannot be > 1 with pvclust")
          stop()
        }

        # get edge below pvals
        hclust_pvpick <-
          pvclust::pvpick(private$hclust, alpha = value, max.only = FALSE,
            pv = pvclust_pv)

        # for each cluster count the number of terms
        # note: hclust_pvpick$clusters is a list
        tab_cluster <- tibble::tibble(
          clusters = hclust_pvpick$clusters,
          n = purrr::map_int(
            hclust_pvpick$clusters,
            length
          )
        )

        get_n_parents <- function(tab_cluster) {
          purrr::map_int(
            tab_cluster$clusters,
            function(x) {
              sum(purrr::map_lgl(
                tab_cluster$clusters,
                function(y) {
                  all(x %in% y)
                }
              )) - 1L
            }
          )
        }


        get_n_children <- function(tab_cluster) {
          purrr::map_int(
            tab_cluster$clusters,
            function(x) {
            sum(purrr::map_lgl(
                tab_cluster$clusters,
                function(y) {
                all(y %in% x)
                }
            )) - 1L
            }
          )
        }

        filter_min_size <- function(tab_cluster, min_size_only_to_children) {
          tab_cluster$n_parents <- get_n_parents(tab_cluster)
          if (min_size_only_to_children) {
            results <- dplyr::filter(tab_cluster,
              .data$n >= .env$min_size &
              .data$n_parents > 0
            )
          } else {
            results <- dplyr::filter(tab_cluster,
              .data$n >= .env$min_size
            )
          }
          return(results)
        }

        filter_max_size <- function(tab_cluster, max_size_only_to_parents) {
          tab_cluster$n_children <- get_n_children(tab_cluster)
          if (max_size_only_to_parents) {
            results <- dplyr::filter(tab_cluster,
              .data$n <= .env$max_size &
              .data$n_children > 0
            )
          } else {
            results <- dplyr::filter(tab_cluster,
              .data$n <= .env$max_size
            )
          }
          return(results)
        }

        if (filter_max_size_first) {
          tab_cluster <- filter_max_size(tab_cluster, max_size_only_to_parents)
          tab_cluster <- filter_min_size(tab_cluster, min_size_only_to_children)
        } else {
          tab_cluster <- filter_min_size(tab_cluster, max_size_only_to_parents)
          tab_cluster <- filter_max_size(tab_cluster, min_size_only_to_children)
        }

        tab_cluster$n_parents <- get_n_parents(tab_cluster)
        tab_cluster$n_children <- get_n_children(tab_cluster)

        if (rm_redundancy_method == "largest") {
          tab_cluster <- dplyr::filter(tab_cluster, .data$n_parents == 0)
        } else {
          tab_cluster <- dplyr::filter(tab_cluster, .data$n_children == 0)
        }

        cls <- unlist(purrr::map(
          seq_len(length(tab_cluster$clusters)),
          function(x) {
            tmp <- rep(x, length(tab_cluster$clusters[[x]]))
            names(tmp) <- tab_cluster$clusters[[x]]
            tmp
          }

        ))
        cl <- length(tab_cluster$clusters) + 1
        orphans <-
          row.names(
            private$i_matrix
          )[! row.names(private$i_matrix) %in% names(cls)]
        cls_bis <- cl:(cl + length(orphans) - 1)
        names(cls_bis) <- orphans
        private$clusters <- c(cls, cls_bis)
      } else {
        value <- as.integer(value)
        if (value < 1) {
          logging::logdebug("`value` cannot be < 1 without pvclust")
          stop()
        }
        private$clusters <- cutree(private$hclust, k = value)
      }
      if (length(private$clusters) == 0) {
        logging::logwarn("No cluster was found with this parametrization")
      }
    },

    #' @description
    #' annotate and get significant results with cluster
    annotate_clusters_and_get_significant_results = function() {
      if (is.null(private$clusters)) {
        logging::logdebug("clusters not set yet")
        stop()
      }
      results <- self$get_significant_results()
      results[, cluster := private$clusters[.SD[["term"]]]]
      results
    },

    #' @description
    #' is the clustering built with pvclust
    #' @return boolean or NULL if no clustering
    is_using_pvclust = function() {
      private$using_pvclust
    },


    #' @description
    #' Prepare a translate dictionary with the different type of possible ids.
    #' @param to type of label to use:
    #'  - "id"
    #'  - "reduced_label"
    #'  - "reduced_label_with_code"
    #'  - "reduced_label_with_id"
    #'  - "full_name"
    #'  - "empty"
    #' @param append_cluster_id append clutser id to labels, default is false
    #' @param append_sep separator between ids
    #' @param reduced_label_max_size maximum length for term label
    #' number of gene per term is used for ordering.
    #' @return translte diuctionary
    get_label_dict = function(to = "reduced_label_with_code",
      append_cluster_id = FALSE,
      append_sep = "#",
      reduced_label_max_size = 40) {
      results <- setNames(switch(to,
        id = private$term_ids,
        reduced_label = paste(
          stringr::str_sub(private$term_names, end = reduced_label_max_size),
          "...", sep = ""),
        reduced_label_with_code = paste(
          stringr::str_sub(private$term_names, end = reduced_label_max_size),
          "...",
          private$raw_ann$ann_name,
          sep = ":"
        ),
        reduced_label_with_id = paste(
          stringr::str_sub(private$term_names, end = reduced_label_max_size),
          "...",
          private$term_ids,
          sep = ":"
        ),
        full_name = private$term_names,
        empty = rep("", length(private$term_ids))
      ), private$term_ids)
      if (append_cluster_id) {
        results = paste0(results, append_sep, private$clusters[private$term_ids])
        names(results) = private$term_ids
      }
      results
    },

    #' @description
    #' prepare results for plot
    #' @param max_cluster maximum number of cluster to plot (default is
    #' no maximum)
    #' @param max_term_per_cluster maximum number of term per cluster (default
    #' is no maximum)
    #' @param ordered_by_pval whether to order terms using p-value. If False the
    #' number of gene per term is used for ordering.
    #' @param keep_only_signif if False all enrichment will be kept if term is
    #' significant in at least one place. In graph value not significant will
    #' appears differently (cross in place of a cricle)
    #' @param in_batch vector of batch code to keep
    #' @param in_group vector of group code to keep
    #' @param in_type vector of type code to keep
    #' @param ... passed to get_label_dict
    #' @return ggplot2 object
    prepare_results_for_plot = function(max_cluster = NULL,
      max_term_per_cluster = NULL, ordered_by_pval = TRUE,
      keep_only_signif = TRUE, in_batch = NULL, in_group = NULL, in_type = NULL,
      ...) {
      dt <- self$annotate_clusters_and_get_significant_results()

      dt <- dt[order(-pval)]

      if (keep_only_signif) {
        dt <- dt[dt$significant]
      }

      term2label <- self$get_label_dict(...)

      cl_qvals <- dt[, .(min_p = min(pval)), by = cluster][order(min_p)]

      cl_sorted <- cl_qvals$cluster

      dt[, cluster := factor(cluster, levels = cl_sorted)]


      if (! is.null(private$batch_labels)) {
          dt[, batch := factor(private$batch_labels[as.character(batch)],
            levels = private$batch_labels)]
      }

      if (! is.null(private$group_labels)) {
          dt[, group := factor(private$group_labels[as.character(group)],
            levels = private$group_labels)]
      }

      if (! is.null(max_cluster)) {
        if (length(cl_sorted) > max_cluster) {
          cl_top <- cl_sorted[1: max_cluster]
        } else {
          cl_top <- cl_sorted
        }
        dt <- dt[cluster %fin% cl_top]
      }

      data_sum <- dt[,
        .(min_p = min(pval, na.rm = TRUE)), by = .(term, cluster, name)]
      data_sum <- data_sum[order(min_p)]

      terms_order_pval <- data_sum$term

      if (! is.null(max_term_per_cluster)) {
        data_sum <- data_sum[, head(.SD, max_term_per_cluster), by = cluster,
          .SDcols = c("min_p", "term")]
        top_terms <- data_sum$term
        dt <- dt[term %fin% top_terms]
      }

      if (! is.null(in_batch)) {
        dt <- dt[as.character(batch) %fin% in_batch]
      }
      if (! is.null(in_group)) {
        dt <- dt[as.character(group) %fin% in_group]
      }
      if (! is.null(in_type)) {
        dt <- dt[as.character(type) %fin% in_type]
      }

      # Summarise Ngenes, parent_Ngenes by term and cluster, then arrange by
      # Ngenes
      df_size <- unique(dt[, .(term, cluster, Ngenes, parent_Ngenes)])
      df_size <- df_size[,
        .(Ngenes = max(Ngenes), parent_Ngenes = max(parent_Ngenes)),
        by = .(term, cluster)]
      df_size <- df_size[order(-Ngenes)]
      term_order_per_ngenes <- as.character(df_size$term)

      df_size <- melt(df_size, id.vars = c("term", "cluster"))

      if (ordered_by_pval) {
        term_order <- terms_order_pval
      } else {
        term_order <- term_order_per_ngenes
      }

      dt[, `:=`(
        tlabel = term2label[term],
        term = factor(term, levels = rev(term_order))
      )]
      df_size[, `:=`(
        tlabel = term2label[term],
        term = factor(term, levels = rev(term_order))
      )]

      list(
        pval = dt[order(term)],
        size = df_size[order(term)],
        term2label = term2label
      )
    },

    #' @description
    #' plot_enrich with pval
    #' @param ... arguments passed to `prepare_results_for_plot`
    #' @return ggplot2 object
    plot_with_pval = function(p_type = "qval_bonferroni", ...) {

      if (! p_type %in% c("pval", "qval_bh", "qval_bonferroni")) {
        logging::logerror("`p_type` must be pval, pval_bh or pval_bonferroni")
        stop()
      }

      list2color <- c(
        up = UP_COLOR,
        upregulated = UP_COLOR,
        Up = UP_COLOR,
        Upregulated = UP_COLOR,
        UP = UP_COLOR,
        UPREGULATED = UP_COLOR,
        down = DOWN_COLOR,
        downregulated = DOWN_COLOR,
        Down = DOWN_COLOR,
        Downregulated = DOWN_COLOR,
        DOWN = DOWN_COLOR,
        DOWNREGULATED = DOWN_COLOR,
        regulated = NA_COLOR,
        all = NA_COLOR,
        deregulated = NA_COLOR,
        deg = NA_COLOR,
        Regulated = NA_COLOR,
        All = NA_COLOR,
        Deregulated = NA_COLOR,
        Deg = NA_COLOR,
        REGULATED = NA_COLOR,
        ALL = NA_COLOR,
        DEREGULATED = NA_COLOR,
        DEG = NA_COLOR
      )

      plist <- self$prepare_results_for_plot(...)
      df <- plist$pval
      actual_list2color <- list2color[levels(df$type)]
      actual_list2color[is.na(actual_list2color)] <- "#000000"

      ggplot2::ggplot(
          df,
          ggplot2::aes(
              .data$type,
              .data$term,
              size = -log10(.data[[p_type]]),
              fill = .data$type,
              color = .data$type,
              shape = as.character(.data$significant)
          )
      ) +
      ggplot2::geom_point() +
      ggplot2::scale_y_discrete(labels = plist$term2label) +
      ggh4x::facet_nested(
        formula("cluster~batch+group"),
        space = "free_y",
        scales = "free_y") +
      THEME_NEXOMIS +
      ggplot2::theme(axis.text.x = ggplot2::element_blank()) +
      ggplot2::labs(x = NULL, y = NULL, size = "score") +
      ggplot2::theme(
        legend.position = "bottom", legend.justification = "right") +
      ggplot2::scale_fill_manual(values = actual_list2color) +
      ggplot2::scale_color_manual(values = actual_list2color) +
      ggplot2::scale_shape_manual(values = c("TRUE" = 19, "FALSE" = 4)) +
      ggplot2::guides(
        fill = ggplot2::guide_legend(nrow = 2, byrow = TRUE),
        size = ggplot2::guide_legend(nrow = 2, byrow = TRUE),
        shape = "none"
      )
    },

    #' @description
    #' plot_enrich with size
    #' @param ... arguments passed to `prepare_results_for_plot`
    #' @return ggplot2 object
    plot_with_size = function(...) {

      plist <- self$prepare_results_for_plot(...)
      df_size <- plist$size

      log10_minor_break <- function(...) {
        function(x) {
          minx <- floor(min(log10(x), na.rm = TRUE)) - 1
          maxx <- ceiling(max(log10(x), na.rm = TRUE)) + 1
          n_major      <- maxx - minx + 1
          major_breaks <- seq(minx, maxx, by = 1)
          minor_breaks <-
            rep(log10(seq(1, 9, by = 1)), times = n_major) +
            rep(major_breaks, each = 9)
          return(10^(minor_breaks))
        }
      }

      t <- c("1", "", "", "", "", "", "", "", "",
          "10", "", "", "", "", "", "", "", "",
          "100", "", "", "", "", "", "", "", "",
          "1e3", "", "", "", "", "", "", "", "",
          "1e4")

      names(t) <- log10_minor_break()(c(10, 1000))[1:37]

      ggplot2::ggplot(
        df_size,
        ggplot2::aes(.data$term, .data$value, alpha = .data$variable)) +
        ggh4x::facet_nested(
          formula("cluster ~ ."), space = "free", scales = "free") +
        ggplot2::coord_flip() +
        ggplot2::scale_y_log10(breaks = as.numeric(names(t)), labels = t) +
        ggplot2::scale_x_discrete(labels = plist$term2label) +
        ggplot2::geom_bar(stat = "identity", position = "identity") +
        ggplot2::labs(y = NULL, x = NULL, title = NULL) +
        ggplot2::guides(
          size = "none", shape = "none", fill = "none", color = "none",
          alpha = ggplot2::guide_legend(title = "", nrow = 2, byrow = TRUE)
        ) +
        THEME_NEXOMIS +
        ggplot2::scale_alpha_manual(
          values = c(parent_Ngenes = 0.35, Ngenes = 1),
          labels = c(parent_Ngenes = "supersets", Ngenes = "set")
        ) +
        ggplot2::theme(legend.position = "bottom",
          legend.justification = "right"
        )

    },

    #' @description
    #' Combination of enrichment plots with p-value per lists and term sizes
    #' @param rel_widths relative width of the combined plots
    #' @param ... arguments passed to `prepare_results_for_plot`
    plot_combined = function(rel_widths = c(6, 1), ...) {
      cowplot::plot_grid(
        self$plot_with_pval(...) +
          ggplot2::theme(strip.text.y = ggplot2::element_blank()),
        self$plot_with_size(...) +
          ggplot2::theme(axis.text.y = ggplot2::element_blank()),
    nrow = 1, align = "h", rel_widths = rel_widths, axis = "bt")
    },

    #' @description
    #' count gene per cluster
    #' @param id2name dictionary to translate gene ids
    #' @param new_name_label name for the column of genes after translation
    #' @param verbose replace column names with meaningful ones.
    #' - gene -> Gene
    #' - cluster -< Cluster
    #' - intra_x_term -> Intra
    #' - extra_x_term -> Extra
    #' - x_cluster -> Cluster
    #' - x_batch -> Batch
    #' - x_group -> Group
    #' @return count table wih the following names:
    #' - gene: [character] gene 
    #' - cluster: [integer] cluster (set of terms) identitified by its number
    #' - intra_x_term: [integer] Number of occurence for the gene in the terms
    #'   intra-cluster
    #' - extra_x_term: [integer] Number of occurence for the gene in the terms globally
    #' - x_cluster:  [integer] Number of occurence for the gene in clusters
    #' - x_batch: Number of occurence for the gene in batches
    #' - x_group: Number of occurence for the gene in groups
    count_gene_per_cluster = function(id2name = NULL,
      new_name_label = "Gene Name", verbose = FALSE) {

      dt_gene_x_cluster_x_term <- data.table::rbindlist(lapply(
        unique(private$clusters),
        function(x) {
          terms_in_cluster <- names(private$clusters)[private$clusters == x]
          indices <- fastmatch::fmatch(terms_in_cluster, private$term_ids)
          dt <- private$raw_ann[term %fin% indices]
          result <- data.table::rbindlist(lapply(
            seq_len(nrow(dt)),
            function(i) {
              dt_row <- dt[i, ]
              data.table::data.table(
                term = dt_row$term,
                gene = dt_row$genes[[1]]
              )
            }
          ))
          result[, cluster := x]
          return(result)
        }
      ))

      dt_gene_x_design <- data.table::rbindlist(lapply(
        unique(private$results$batch),
        function(batch_in) {
          data.table::rbindlist(lapply(
            unique(private$results$group),
            function(group_in) {
              dt <- private$results[batch == batch_in & group == group_in]
              data.table::data.table(
                batch = batch_in,
                group = group_in,
                gene = unique(unlist(dt$gene_input))
              )
            }
          ))
        }
      ))

      dt_intra_x_term <- dt_gene_x_cluster_x_term[,
        .(intra_x_term = data.table::uniqueN(term)), by = .(gene, cluster)]
      dt_extra_x_term <- dt_gene_x_cluster_x_term[,
        .(extra_x_term = data.table::uniqueN(term)), by = .(gene)]
      dt_x_cluster <- dt_gene_x_cluster_x_term[,
        .(x_cluster = data.table::uniqueN(cluster)), by = .(gene)]
      dt_x_batch <- dt_gene_x_design[, .(x_batch = data.table::uniqueN(batch)),
        by = .(gene)]
      dt_x_group <- dt_gene_x_design[, .(x_group = data.table::uniqueN(
        interaction(batch, group, drop = TRUE))), by = .(gene)]
      summary_dt <- data.table::merge.data.table(dt_intra_x_term,
        dt_extra_x_term, by = "gene")
      summary_dt <- data.table::merge.data.table(summary_dt, dt_x_cluster, by = "gene")
      summary_dt <- data.table::merge.data.table(summary_dt, dt_x_batch, by = "gene")
      summary_dt <- data.table::merge.data.table(summary_dt, dt_x_group, by = "gene")

      summary_dt[, gene := private$gene_ids[gene]]
      
      if (! is.null(id2name)) {
        # How to make it to the 3rd column position ?
        summary_dt[, (new_name_label) := as.character(id2name[gene])]
      }
      data.table::setorder(summary_dt,
        cluster, - intra_x_term, - extra_x_term, gene)

      if (verbose) {
        old_names <- c("gene", "cluster", "intra_x_term", "extra_x_term",
          "x_cluster", "x_batch", "x_group")
        new_names <- c("Gene", "Cluster", "Intra",
          "Extra", "Clusters", "Batches", "Groups")
        if (! is.null(id2name)) {
          old_names <- c(old_names, new_name_label)
          new_names <- c(new_names, new_name_label)
        }
        data.table::setnames(summary_dt, old_names, new_names)
      }
      if (! is.null(id2name)) {
        current_cols <- names(summary_dt)
        desired_order <-
          c(current_cols[1], new_name_label, current_cols[2:(length(current_cols)-1)])
        data.table::setcolorder(summary_dt, desired_order)
      }

      return(summary_dt)
    },

    #' Build Cluster Summary Table
    #'
    #' Constructs a summary table of the enrichment clusters.
    #'
    #' @description
    #' This method processes the results from previous enrichment analyses to compile
    #' a detailed summary for each cluster. It focuses on aggregating minimum q-values
    #' (adjusted p-values) for each cluster to assess the overall significance, alongside
    #' other relevant term details.
    #'
    #' @param verbose Boolean, if TRUE, column names in the returned table are more descriptive.
    #' The renaming adjusts column headers to be more human-readable and suitable for reports
    #' or presentations. Specific changes include:
    #'   - "cluster" becomes "Cluster"
    #'   - "term" becomes "Term"
    #'   - "name" becomes "Name"
    #'   - "ann_name" becomes "Database"
    #'   - "Ngenes" becomes "# of Genes"
    #'   - "min_qval" becomes "Min q-value"
    #' @return A tibble table containing a summary of the clusters with the following columns:
    #'   - `cluster`: Cluster identifier (numeric or character based on input data).
    #'   - `term`: Unique identifier for each term within a cluster.
    #'   - `name`: Descriptive name of the term.
    #'   - `ann_name`: Source database from which the term was derived.
    #'   - `Ngenes`: Count of genes associated with each term.
    #'   - `min_qval`: Minimum q-value observed within each cluster, reflecting the lowest
    #'     adjusted p-value computed for terms within the cluster.
    #'   - Additional columns representing each combination of batch, group, and type, formatted
    #'     as 'batch|group|type', each containing the q-value for the term within that context.
    #'
    build_cluster_summary_table = function(verbose = FALSE) {

      select_res <- c(
        "cluster", "term", "name", "ann_name", "Ngenes", "batch", "group",
        "type", "qval_bonferroni"
      )

      scl_table <- self$prepare_results_for_plot()$pval[, ..select_res]
      if (verbose) {
        scl_table[, category := paste(
          private$batch_labels[batch],
          private$group_labels[group],
          type, sep = "|")]
      } else {
        scl_table[, category := paste(batch, group, type, sep = "|")]
      }
      
      scl_table[, c("batch", "group", "type") := NULL]
      
      cl_sum <- scl_table[, .(min_q_val_cl = min(qval_bonferroni)), by = cluster]
      cl2min <- cl_sum$min_q_val_cl
      setattr(cl2min, "names", as.character(cl_sum$cluster))
      scl_table[["min_q_val_cl"]] <- cl2min[as.character(scl_table[["cluster"]])]
      setattr(scl_table[["min_q_val_cl"]], "names", NULL)

      cl_sum <- scl_table[, .(min_q_val_term = min(qval_bonferroni)), by = term]
      cl2min <- cl_sum$min_q_val_term
      setattr(cl2min, "names", as.character(cl_sum$term))
      scl_table[["min_q_val_term"]] <- cl2min[as.character(scl_table[["term"]])]
      setattr(scl_table[["min_q_val_term"]], "names", NULL)

      wide_scl_table <- dcast(scl_table, ... ~ category,
        value.var = "qval_bonferroni")

      data.table::setorder(wide_scl_table, min_q_val_cl, min_q_val_cat)

      setattr(wide_scl_table$cluster, "names", NULL)
      if (verbose) {
        setnames(wide_scl_table, old = c("cluster", "term", "name", "ann_name",
          "Ngenes", "min_q_val_cl", "min_q_val_term"),
          new = c("Cluster", "Term", "Name", "Database", "# of Genes", "Cluster min q-value", "Term min q-value"))
      }

      return(wide_scl_table)
    },

    #' @description
    #' write the enrichment data into xlsx files
    #' @param output_folder output folder (default is .)
    #' @param id2name dictionary to translate gene ids
    #' @param new_name_label name for the column of genes after translation
    #' @param write_cluster_summary whether to write or not the summary
    #' of all cluster
    #' @return None
    write_xlsx = function(output_folder = ".", id2name = NULL,
      new_name_label = "Tested Gene Name", write_cluster_summary = TRUE) {
      new_names <- c(
        "Term ID",
        "Term name",
        "Cluster",
        "P value",
        "Q value (Bonferroni)",
        "Q value (Benj. Hoch.)",
        "Term type",
        "Term parents IDs",
        "Tested Gene ID",
        "# of tested genes",
        "# of genes in term",
        "Lowest p value possible"
      )
      if (! is.null(id2name)) {
        new_names <- c(new_names, new_name_label)
      }
      if (!dir.exists(output_folder)) {
        dir.create(output_folder)
      }
      d_input <- dplyr::distinct(
        self$get_results(verbose = TRUE)[,
          c("batch", "group", "type", "ann_name")])
      for (batch in unique(d_input$batch)) {
        d_input_batch <- dplyr::filter(d_input, .data$batch == .env$batch)
        for (group in unique(d_input_batch$group)) {
          v_tables <- list()
          d_input_batch_group <- dplyr::filter(d_input_batch,
            .data$group == .env$group)
          for (i in seq_along(d_input_batch_group)) {
            type <- d_input_batch_group$type[i]
            ann_name <- d_input_batch_group$ann_name[i]
            sheet_name <- paste(type, ann_name, sep = "_")
            v_tables[[sheet_name]] <-
              dplyr::filter(self$get_results(verbose = TRUE),
                .data$batch == .env$batch &
                .data$group == .env$group &
                .data$type == .env$type &
                .data$ann_name == .env$ann_name)$enrich[[1]]
            v_tables[[sheet_name]][["cluster"]] <-
              private$clusters[v_tables[[sheet_name]][["term"]]]
            v_tables[[sheet_name]] <- v_tables[[sheet_name]][,
                c("term",
                "name",
                "cluster",
                "pval",
                "qval_bonferroni",
                "qval_bh",
                "type",
                "parents",
                "GI",
                "NGI",
                "Ngenes",
                "pmin")]
            v_tables[[sheet_name]][["parents"]] <- purrr::map_chr(
              v_tables[[sheet_name]][["parents"]],
              function(x) paste(x, collapse = ";"))
            if (! is.null(id2name)) {
              v_tables[[sheet_name]][["tGI"]] <- purrr::map_chr(
                v_tables[[sheet_name]][["GI"]],
                function(x) paste(na.omit(id2name[x]), collapse = ";"))
            }
            v_tables[[sheet_name]]
            v_tables[[sheet_name]][["GI"]] <- purrr::map_chr(
              v_tables[[sheet_name]][["GI"]],
              function(x) paste(x, collapse = ";"))
            names(v_tables[[sheet_name]]) <- new_names
            class(v_tables[[sheet_name]][["P value"]]) <- "scientific"
            class(v_tables[[sheet_name]][["Q value (Bonferroni)"]]) <-
              "scientific"
            class(v_tables[[sheet_name]][["Q value (Benj. Hoch.)"]]) <-
              "scientific"
            class(v_tables[[sheet_name]][["Lowest p value possible"]]) <-
              "scientific"
          }
          options(openxlsx.numFmt = "#,#0.00")
          openxlsx::write.xlsx(v_tables,
            paste(
              output_folder,
              paste(batch, group, "fctBio.xlsx", sep = "_"),
              sep = "/"
            ),
            asTable = TRUE, tableStyle = "TableStyleLight1")
        }
      }
      if (write_cluster_summary) {
        scl_table <- self$build_cluster_summary_table(verbose = TRUE)
        for (qname in names(scl_table)[6:ncol(scl_table)]) {
            class(scl_table[[qname]]) <- "scientific"
        }
        options(openxlsx.numFmt = "#,#0.00")
        openxlsx::write.xlsx(
          list(
            Summary = scl_table,
            Occurences = self$count_gene_per_cluster(id2name = id2name,
              new_name_label = new_name_label, verbose = TRUE)
          ),
          paste(
              output_folder,
              "cluster_summary_fctBio.xlsx",
              sep = "/"
            ),
          asTable = TRUE, tableStyle = "TableStyleLight1"
        )
      }
    }
  ),
  private <- list(
    results = NULL,
    data_nested_id = "uniprot",
    significant_results = NULL,
    significant_terms = NULL,
    clusters = NULL,
    ann_sources = NULL,
    ann_types = NULL,
    batch_labels = NULL,
    group_labels = NULL,
    i_matrix = NULL,
    p_matrix = NULL,
    hclust = NULL,
    using_pvclust = NULL,
    raw_ann = NULL,
    gene_ids = NULL,
    term_ids = NULL,
    term_names = NULL,
    term_types = NULL
  )
)
