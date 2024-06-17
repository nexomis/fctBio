#' @include utils.r

NULL

#' Query expression tags
#' @description
#'
#' Function to build a query to find a list of genes.
#'
#' @param nested_df nested tibble with
#' - gene sets given in a tibble in `data_col` variable and in `data_nested_id`
#' nested variable
#' - eventually enrichment data in the enrich column
#' from the filter_and_get_nested_results method of NestedEnrich object
#' @param filter_name filter to apply
#' - batch : select genes if present in batch(s)
#' - batch_label : select genes if present in batch label(s)
#' - group : select genes if present in group(s)
#' - group_label : select genes if present in group label(s)
#' - term_id : select genes if part of a term based on the `enrich_gene_col`.
#' - term_name : select genes if part of a term based on the `enrich_gene_col`.
#' - cluster : select genes if part of a cluster (significant term).
#' - intra_occurences : select genes based on their minimal occurences in terms
#' intra cluster
#' - extra_occurences : select genes based on their minimal occurences in terms
#' - cluster_occurences : select genes based on their minimal occurences in
#' cluster
#' @param value value or vector of value to be filtered on based on filter_name
#' @param genes_per_clusters
#' from the count_gene_per_cluster method of NestedEnrich object with columns:
#' - "ID"
#' - "Cluster"
#' - "# in enriched term intra"
#' - "# in batch"
#' - "# in group"
#' - "# in enriched term"
#' - "# in cluster"
#' @param batch_labels batch labels
#' @param group_labels group labels
#' @param enrich_gene_col coloum to select genes (default is genes)
#' - GI for gene of intermterest
#' - genes for all genes in terms
#' @param data_col see `nested_df`
#' @param data_univ_col same as `data_col` but for the "universe"
#' @param data_nested_id see `nested_df`
#' @param enrich_col column for the nested enrichment data. Default is "enrich".
#' Is set to null if not in nested_df
#' @export
query_genes <- function(nested_df, filter_name, value,
  genes_per_clusters = NULL,
  batch_labels = NULL, group_labels = NULL,
  data_col = "data", data_nested_id = "uniprot",
  enrich_gene_col = "genes"
) {

  if (! all(value %in% query_genes_choices(nested_df, filter_name,
  genes_per_clusters, batch_labels, group_labels))) {
    warning("value not a valid choice")
  }

  if (filter_name == "batch_label" && is.null(batch_labels)) {
    warning("filtering batch_label without label dict")
    batch_labels <- query_genes_choices(nested_df, filter_name)
    names(batch_labels) <- batch_labels
  }

  if (filter_name == "group_label" && is.null(group_labels)) {
    warning("filtering group_label without label dict")
    group_labels <- query_genes_choices(nested_df, filter_name)
    names(group_labels) <- group_labels
  }

  switch(filter_name,
    batch = {
      unique(unlist(purrr::map(
        dplyr::filter(nested_df, .data$batch %in% value)[[data_col]],
        function(df) {
          df[[data_nested_id]]
        }
      )))
    },
    batch_label = {
      unique(unlist(purrr::map(
        dplyr::filter(nested_df,
          as.character(batch_labels[.data$batch]) %in% value)[[data_col]],
        function(df) {
          df[[data_nested_id]]
        }
      )))
    },
    group = {
      unique(unlist(purrr::map(
        dplyr::filter(nested_df, .data$group %in% value)[[data_col]],
        function(df) {
          df[[data_nested_id]]
        }
      )))
    },
    group_label = {
      unique(unlist(purrr::map(
        dplyr::filter(nested_df,
          as.character(group_labels[.data$group]) %in% value)[[data_col]],
        function(df) {
          df[[data_nested_id]]
        }
      )))
    },
    term_id = {
      if ("enrich" %in% names(nested_df)) {
        choices <- unique(unlist(purrr::map(
          nested_df[["enrich"]],
          function(df) {
            unique(unlist(
              dplyr::filter(df, .data$term %in% .env$value)[[enrich_gene_col]]
            ))
          }
        )))
      } else {
        choices <- character()
      }
      choices
    },
    term_name = {
      if ("enrich" %in% names(nested_df)) {
        choices <- unique(unlist(purrr::map(
          nested_df[["enrich"]],
          function(df) {
            unique(unlist(
              dplyr::filter(df, .data$name %in% .env$value)[[enrich_gene_col]]
            ))
          }
        )))
      } else {
        choices <- character()
      }
      choices
    },
    cluster = {
      if (is.null(genes_per_clusters)) {
        choices <- character()
      } else {
        choices <- unique(dplyr::filter(genes_per_clusters,
          .data[["Cluster"]] %in% value)[["ID"]])
      }
      choices
    },
    intra_occurences = {
      if (is.null(genes_per_clusters)) {
        choices <- integer()
      } else {
        choices <- unique(dplyr::filter(genes_per_clusters,
          .data[["# in enriched term intra"]] >= value)[["ID"]])
      }
      choices
    },
    extra_occurences = {
      if (is.null(genes_per_clusters)) {
        choices <- integer()
      } else {
        choices <- unique(dplyr::filter(genes_per_clusters,
          .data[["# in enriched term"]] >= value)[["ID"]])
      }
      choices
    },
    cluster_occurences = {
      if (is.null(genes_per_clusters)) {
        choices <- integer()
      } else {
        choices <- unique(dplyr::filter(genes_per_clusters,
          .data[["# in cluster"]] >= value)[["ID"]])
      }
      choices
    }
  )
}

#' Get possible choices for a query
#' @description
#' get choices to query genes (see query_genes function)
#' @inheritParams query_genes
#' @export
query_genes_choices <- function(nested_df, filter_name,
  genes_per_clusters = NULL,
  batch_labels = NULL, group_labels = NULL
) {
  switch(filter_name,
    batch = {
      unique(nested_df$batch)
    },
    batch_label = {
      if (is.null(batch_labels)) {
        choices <- unique(nested_df$batch)
      } else {
        choices <- as.character(batch_labels[unique(nested_df$batch)])
      }
      choices
    },
    group = {
      unique(nested_df$group)
    },
    group_label = {
      if (is.null(group_labels)) {
        choices <- unique(nested_df$group)
      } else {
        choices <- as.character(group_labels[unique(nested_df$group)])
      }
      choices
    },
    term_id = {
      if ("enrich" %in% names(nested_df)) {
        choices <- sort(unique(unlist(purrr::map(
          nested_df[["enrich"]],
          function(df) {
            df[["term"]]
          }
        ))))
      } else {
        choices <- character()
      }
      choices
    },
    term_name = {
      if ("enrich" %in% names(nested_df)) {
        choices <- sort(unique(unlist(purrr::map(
          nested_df$enrich,
          function(df) {
            df$name
          }
        ))))
      } else {
        choices <- character()
      }
      choices
    },
    cluster = {
      if (is.null(genes_per_clusters)) {
        choices <- character()
      } else {
        choices <- unique(genes_per_clusters[["Cluster"]])
      }
      choices
    },
    intra_occurences = {
      if (is.null(genes_per_clusters)) {
        choices <- integer()
      } else {
        choices <- seq(0,max(genes_per_clusters[["# in enriched term intra"]]))
      }
      choices
    },
    extra_occurences = {
      if (is.null(genes_per_clusters)) {
        choices <- integer()
      } else {
        choices <- seq(0, max(genes_per_clusters[["# in enriched term"]]))
      }
      choices
    },
    cluster_occurences = {
      if (is.null(genes_per_clusters)) {
        choices <- integer()
      } else {
        choices <- seq(0, max(genes_per_clusters[["# in cluster"]]))
      }
      choices
    }
  )
}
