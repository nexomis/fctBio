# TODO : default comportement more basic (same as web interface).
#        define this function in 2 distinct step by 2 function.
#        First in order to create 5 vectors (saved in 5 collumn of one dataframe in order to facilate junction and comparison) :
#          3 about input (tf, target, na) and 2 about db request (tf, target).
#        And a second function in order to manae output result from user option settings.
#
#        Finaaly, in future, a last function must be create in order to summarize all db request function sysbio in a single.
#        In this order, it is more easy if variable name is universal.
#        PS: detailled output isn't usefull



#' Extend a list of genes by adding linked transcription factors (TFs) and/or target genes.
#'
#' @param trrust_db Path to the TRRUST database in tsv format.
#' @param input_genes A character vector of input genes.
#' @param add_tf Logical, whether to add transcription factors regulating the input genes (keeping 'tf' already present in input genes).
#' @param add_target Logical, whether to add target genes regulated by the input genes (keeping 'target' already present in input genes).
#' @param keep_input_without_hit Logical, whether to keep all input genes (including 'tf', 'target' and unannotated gene).
#' @param ignore_single_pmid Logical, whether to ignore relations supported by only one PMID reference.
#' @param detailled_output [not implemented] To return in addition of vector of gene names, a dataframe with relation of each genes pairs.
#' @return A character vector of extended genes.
#' @examples
#' extended_genes <- extend_trrust(
#'   trrust_db = "path/to/trrust_rawdata.human.tsv",
#'   input_genes = c("BAX", "MYC", "TP53"),
#'   add_tf = TRUE,
#'   add_target = TRUE,
#'   keep_input_without_hit = TRUE,
#'   ignore_single_pmid = FALSE
#' )
#' @export
extend_trrust <- function(trrust_db = system.file("extdata/human_network_annotations/trrust/data.tab", package = "fctBio"),
                          input_genes,
                          add_tf = TRUE,
                          add_target = TRUE,
                          keep_input_without_hit = TRUE,
                          ignore_single_pmid = FALSE,
                          detailled_output = FALSE) {
  
  # Load TRRUST database from tsv file and if specified, exclude relation supported by only 1 PMID reference
  trrust_data <- read.delim(trrust_db, header = TRUE)
  if (ignore_single_pmid) {
    trrust_data$PMID_count <- sapply(strsplit(as.character(trrust_data$ref_PMID), ";"), length)
    trrust_data_filtered <- subset(trrust_data, PMID_count > 1)
    trrust_data_filtered$PMID_count <- NULL
  }
  else {
    trrust_data_filtered <- trrust_data
  }
  
  # Force keeping of original input genes
  results <- if (keep_input_without_hit) input_genes else character(0)
  
  # Add transcription factors regulating the input genes
  if (add_tf) {
    tf_regulating_input <- unique(trrust_data_filtered[trrust_data_filtered$target %in% input_genes, ]$tf)
    tf_in_input <- unique(trrust_data[trrust_data$tf %in% input_genes, ]$tf)
    results <- unique(append(results, c(tf_regulating_input, tf_in_input)))
  }
  
  # Add target genes regulated by the input genes
  if (add_target) {
    target_under_input <- unique(trrust_data_filtered[trrust_data_filtered$tf %in% input_genes, ]$target)
    target_in_input <- unique(trrust_data[trrust_data$target %in% input_genes, ]$target)
    results <- unique(append(results, c(target_under_input, target_in_input)))
  }
  
  return(unique(results))
}






#### DRAFT: ressource to implement 'detailled_output' option for 'extend_network_trrust' function
# # import DB from tab file
# trrust_data <- read.delim("C:/Users/abdel/Downloads/trrust_rawdata.human.tsv", header = FALSE,
#                           col.names = c("TranscriptionFactor", "TargetGene", "RegulationMode", "PMID"))
# 
# # exemple input genes
# input_genes <- c("BAX", "MYC", "TP53")
# 
# ### Format 1
# relations_as_targets <- trrust_data[trrust_data$TargetGene %in% input_genes, ]
# relations_as_targets_df1 <- data.frame(Input_gn = relations_as_targets$TargetGene,
#                                        Direction = "by",
#                                        Pairs_gn = relations_as_targets$TranscriptionFactor,
#                                        Mode = relations_as_targets$RegulationMode,
#                                        Ref_PMID = relations_as_targets$PMID,
#                                        stringsAsFactors = FALSE)
# 
# # add input genes without hit
# no_relation_genes <- setdiff(input_genes, unique(relations_as_targets$TargetGene))
# no_relation_df1 <- data.frame(Input_gn = no_relation_genes, Direction = "NA", Pairs_gn = "NA", Mode = "NA", Ref_PMID = "NA", stringsAsFactors = FALSE)
# dataframe_1 <- rbind(relations_as_targets_df1, no_relation_df1)
# 
# ### Format 2
# relations_as_tfs <- trrust_data[trrust_data$TranscriptionFactor %in% input_genes, ]
# relations_as_tfs_df2 <- data.frame(Input_gn = relations_as_tfs$TranscriptionFactor,
#                                    Type = "TF",
#                                    TF_gn = relations_as_tfs$TranscriptionFactor,
#                                    Target_gn = relations_as_tfs$TargetGene,
#                                    RegulationMode = relations_as_tfs$RegulationMode,
#                                    Ref_PMID = relations_as_tfs$PMID,
#                                    stringsAsFactors = FALSE)
# 
# # add input genes without hit
# no_relation_tfs <- setdiff(input_genes, unique(relations_as_tfs$TranscriptionFactor))
# no_relation_df2 <- data.frame(Input_gn = no_relation_tfs, Type = "NA", TF_gn = "NA", Target_gn = "NA", RegulationMode = "NA", Ref_PMID = "NA", stringsAsFactors = FALSE)
# dataframe_2 <- rbind(relations_as_tfs_df2, no_relation_df2)
# 
# 
