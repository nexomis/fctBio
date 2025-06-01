#' @include utils.r

NULL

#' Load Annotation Space
#'
#' This function loads gene annotations from specified sources, such as
#' databases or files, and processes them into a structured format suitable for
#' further analysis.

#' Annotations typically define the notion of "terms" which can represent
#' biological pathways, communities of genes in coexpression networks, or other
#' relevant biological groupings. These terms might be organized hierarchically
#' using directed acyclic graphs (DAGs) in systems like Gene Ontology or
#' Reactome, allowing for the representation of complex biological
#' relationships.
#'
#' @param ann_sources Named character vector with keys as database codes and
#' values as paths to tabular files containing annotations.
#' Each file must include the following columns:
#'   - `term`: Short identifier of the term.
#'   - `genes`: Semicolon-separated list of gene identifiers associated with
#'     each term.
#'   - `type`: Character string defining the category of the term, e.g.,
#'     Biological Process, Molecular Function, or Cellular Component in the
#'     context of Gene Ontology.
#'   - `name`: Descriptive name or long identifier of the term.
#'   - `parents`: Comma-separated list of parent terms' identifiers, reflecting
#'     the hierarchical structure in databases using DAGs to organize terms.
#' @param ann_types Named character vector with keys as database codes and
#' values specifying
#' the types of terms to retain from the annotations. The default behavior is to
#' keep all types. This parameter filters terms based on their `type` field.
#' @param compute_parent_genes Logical, indicating whether to compute a
#' collective list of genes for each term by incorporating genes from its parent
#' terms. Defaults to TRUE.
#' @param use_parents_union Logical, indicating how to handle gene lists from
#' parent terms:
#'   - If TRUE, computes the union of genes across all parent terms.
#'   - If FALSE, computes the intersection of genes across all parent terms.
#'
#' @return A data.table object structured as follows:
#'   - `term`: Short identifier of the term.
#'   - `name`: Long name or descriptive identifier of the term.
#'   - `genes`: List of gene identifiers associated with each term.
#'   - `parents`: List of parent terms' identifiers.
#'   - `parent_genes`: Computed list of genes derived from parent terms, present
#'     only if
#'     `compute_parent_genes` is TRUE.
#'   - `type`: Category of the term, consistent with the `type` field in the
#'     input files.
#'   - `ann_name`: Database code from which each term was sourced, aligning with
#'     keys in `ann_sources`.
#' @export
load_ann_space <- function(
  ann_sources, ann_types, compute_parent_genes = TRUE,
  use_parents_union = TRUE
) {
  data.table::rbindlist(
    lapply(
      names(ann_sources),
      function(ann_name_input) {

        ann_source <- ann_sources[ann_name_input]
        if (ann_name_input %in% names(ann_types)) {
          type_keep <- ann_types[ann_name_input]
        } else {
          type_keep <- NULL
        }

        ann <- data.table::fread(file = ann_source, sep = "\t", header = TRUE,
          stringsAsFactors = FALSE
        )

        if (!is.null(type_keep)) {
          filtering_vector <- vapply(
            ann$type,
            function(x) {
              any(
                stringr::str_split_1(x, ";") %fin% type_keep
              )
            },
            FUN.VALUE = logical(1L),
            USE.NAMES = FALSE
          )
          ann <- ann[filtering_vector]
        }

        ann[, genes := lapply(
          .SD[["genes"]],
          function(x) (stringr::str_split_1(x, ";"))
        )]

        ann[,
          ann_name := factor(ann_name_input, levels = names(ann_sources))
        ]
        if (! "type" %in% names(ann)) {
          ann[, type := factor("")]
        }

        if (! compute_parent_genes) {
          rcols <- c("term", "name", "genes", "parents", "type", "ann_name")
        } else {

          ann[, parents := lapply(
            .SD[["parents"]],
            function(x) {
              if (x == "") {
                character()
              } else {
                stringr::str_split_1(x, "[;,]{1}")
              }
            }
          )]

          named_genes <- setNames(ann$genes, ann$term)
          ann[, temp := lapply(parents, function(x) named_genes[x])]
          if (use_parents_union) {
            ann[, parent_genes :=
              lapply(temp, function(x) {
                Reduce(
                  function(a, b) {
                    collapse::funique(c(a, b))
                  },
                  x
                )
              })
            ]
          } else {
            ann[, parent_genes :=
              lapply(temp, function(x) {
                Reduce(
                  function(a, b) {
                    a <- c(a, b)
                    a[collapse::fduplicated(a)]
                  },
                  x
                )
              })
            ]
          }
          ann[, temp := NULL]
          rcols <- c("term", "name", "genes", "parents", "parent_genes",
            "type", "ann_name"
          )
        }
        ann[, rcols]
      }
    )
  )
}
