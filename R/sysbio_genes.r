#' @include aaa.r
#' @include utils.r

NULL

#' Extend Gene List Using TRRUST Database
#'
#' @description
#' `extend_trrust` extends a list of input genes by adding linked transcription
#' factors (TFs) and/or target genes based on the TRRUST database. This function
#' allows for flexible gene network expansion by including regulatory relationships
#' from the TRRUST (Transcriptional Regulatory Relationships Unraveled by
#' Sentence-based Text mining) database.
#'
#' The function can selectively add transcription factors that regulate the input
#' genes and/or target genes that are regulated by the input genes. It also
#' provides options to filter relationships based on literature support and to
#' retain input genes that have no regulatory relationships in the database.
#'
#' @param trrust_db A character string specifying the path to the TRRUST database
#' in TSV format. Defaults to the package's internal TRRUST data file.
#' @param input_genes A character vector of input gene symbols to extend.
#' @param add_tf Logical; whether to add transcription factors that regulate
#' the input genes. When `TRUE`, includes TFs that have the input genes as targets,
#' while preserving any input genes that are already TFs in the database.
#' Default is `TRUE`.
#' @param add_target Logical; whether to add target genes that are regulated
#' by the input genes. When `TRUE`, includes targets of input genes that act as TFs,
#' while preserving any input genes that are already targets in the database.
#' Default is `TRUE`.
#' @param keep_input_without_hit Logical; whether to retain input genes that
#' have no regulatory relationships in the database (neither as TFs nor as targets).
#' Default is `TRUE`.
#' @param ignore_single_pmid Logical; whether to exclude regulatory relationships
#' that are supported by only one PubMed ID reference. When `TRUE`, only
#' relationships with multiple literature references are considered. Default is `FALSE`.
#'
#' @return A character vector of unique gene symbols representing the extended
#' gene list. The returned genes include the original input genes (if
#' `keep_input_without_hit` is `TRUE`) plus any additional TFs and/or targets
#' based on the specified parameters.
#'
#' @export
#' @examples
#' # Basic usage with default parameters
#' extended_genes <- extend_trrust(
#'   input_genes = c("BAX", "MYC", "TP53")
#' )
#'
#' # Add only transcription factors, exclude single PMID relationships
#' tf_extended <- extend_trrust(
#'   input_genes = c("BAX", "MYC", "TP53"),
#'   add_tf = TRUE,
#'   add_target = FALSE,
#'   ignore_single_pmid = TRUE
#' )
#'
#' # Add only targets, don't keep input genes without hits
#' target_extended <- extend_trrust(
#'   input_genes = c("BAX", "MYC", "TP53"),
#'   add_tf = FALSE,
#'   add_target = TRUE,
#'   keep_input_without_hit = FALSE
#' )
extend_trrust <- function(
  input_genes,
  trrust_db = system.file(
    "extdata/human_network_annotations/trrust/data.tab", package = "fctBio"
  ),
  add_tf = TRUE,
  add_target = TRUE,
  keep_input_without_hit = TRUE,
  ignore_single_pmid = FALSE
) {

  # Input validation
  if (missing(input_genes) || length(input_genes) == 0L) {
    stop("input_genes must be provided and non-empty", call. = FALSE)
  }

  if (!file.exists(trrust_db)) {
    stop("TRRUST database file not found: ", trrust_db, call. = FALSE)
  }

  # Load TRRUST database
  trrust_data <- read.delim(trrust_db, header = TRUE, stringsAsFactors = FALSE)

  # Filter by PMID count if requested
  if (ignore_single_pmid) {
    trrust_data$pmid_count <- vapply(
      strsplit(as.character(trrust_data$ref_PMID), ";"),
      length,
      integer(1L)
    )
    trrust_data_filtered <- trrust_data[trrust_data$pmid_count > 1L, ]
    trrust_data_filtered$pmid_count <- NULL
  } else {
    trrust_data_filtered <- trrust_data
  }

  # Initialize results with input genes if requested
  results <- if (keep_input_without_hit) input_genes else character(0L)

  # Add transcription factors regulating the input genes
  if (add_tf) {
    # TFs that regulate input genes (input genes as targets)
    tf_regulating_input <- unique(
      trrust_data_filtered[trrust_data_filtered$target %in% input_genes, "tf"]
    )

    # Input genes that are TFs in the database
    tf_in_input <- unique(
      trrust_data[trrust_data$tf %in% input_genes, "tf"]
    )

    results <- unique(c(results, tf_regulating_input, tf_in_input))
  }

  # Add target genes regulated by the input genes
  if (add_target) {
    # Targets regulated by input genes (input genes as TFs)
    target_under_input <- unique(
      trrust_data_filtered[trrust_data_filtered$tf %in% input_genes, "target"]
    )

    # Input genes that are targets in the database
    target_in_input <- unique(
      trrust_data[trrust_data$target %in% input_genes, "target"]
    )

    results <- unique(c(results, target_under_input, target_in_input))
  }

  return(unique(results))
}
