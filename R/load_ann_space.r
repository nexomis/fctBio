#' @include utils.r

NULL

#' Load annotation space
#' @param ann_sources named character vector with database code and database
#' path (see \link{enrich}) ann_source arg
#' @param ann_types named character vector with database code as keys and
#' database types to keep (see \link{enrich}) type_keep arg. Default is to
#' keep all.
#' @param compute_parent_genes whether to compute or not parent genes
#' @return a data.table
#' @export
load_ann_space <- function(ann_sources, ann_types, compute_parent_genes = TRUE,
use_parents_union = TRUE) {
  rbindlist(
    lapply(
      names(ann_sources),
      function(ann_name_input) {

        ann_source <- ann_sources[ann_name_input]
        if (ann_name_input %in% names(ann_types)) {
          type_keep <- ann_types[ann_name_input]
        } else {
          type_keep <- NULL
        }

        ann <- fread(file = ann_source, sep = "\t", header = TRUE,
          stringsAsFactors = FALSE)

        if (!is.null(type_keep)) {
          filtering_vector <- vapply(
            ann$type,
            function(x) {
              any(
                stringr::str_split_1(x, ";") %fin% type_keep
              )
            },
            FUN.VALUE = logical(1),
            USE.NAMES = FALSE
          )
          ann <- ann[filtering_vector]
        }

        ann[, genes := lapply(
          .SD[["genes"]],
          function(x) (stringr::str_split_1(x, ";"))
        )]

        ann[,
          ann_name := factor(ann_name_input, levels = names(ann_sources))]
        if (! "type" %in% names(ann)) {
          ann[, type := factor("")]
        }

        if (! compute_parent_genes) {
          return(
            ann[, c(
              "term", "name", "genes", "parents", "type", "ann_name")]
          )
        } else {

          ann[, parents := lapply(
            .SD[["parents"]],
            function(x) {
              if (x == "") {
                return(character())
              } else {
                return(stringr::str_split_1(x, "[;,]{1}"))
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
          return(
            ann[, c("term", "name", "genes", "parents", "parent_genes",
              "type", "ann_name")]
          )
        }
      }
    )
  )
}