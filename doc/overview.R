## ----setup, include = FALSE---------------------------------------------------
knitr::opts_chunk$set(
  cache = TRUE,
  echo = TRUE,
  out.width = "90%",
  fig.width = 6.5, 
  fig.height = 5.5
)

## ----eval=FALSE---------------------------------------------------------------
#  

## ----echo=TRUE, class.source = "fold-hide"------------------------------------

base <- "go"
file_base <- system.file("extdata",
                         paste("human_functional_annotations",
                               base, "data.tab", sep = "/"),
                         package = "fctBio"
)
tab_base <- read.delim(file_base)

# response to type I interferon
# antiviral innate immune response
l1 <- stringr::str_split(paste(tab_base[tab_base$term %in%
                                          c("GO:0034340", "GO:0140374"), "genes"], collapse = ";"), ";")
l2 <- stringr::str_split(paste(tab_base[tab_base$term %in%
                                          c("GO:0002347"), "genes"], collapse = ";"), ";")
l3 <- stringr::str_split(paste(tab_base[tab_base$term %in%
                                          c("GO:0070914"), "genes"], collapse = ";"), ";")
l4 <- stringr::str_split(paste(tab_base[c(20, 30, 40),
                                        "genes"], collapse = ";"), ";")
l5 <- stringr::str_split(paste(tab_base[tab_base$term %in%
                                          c("GO:0042246"), "genes"], collapse = ";"), ";")
l6 <- stringr::str_split(paste(tab_base[c(25, 35, 45),
                                        "genes"], collapse = ";"), ";")

l1 <- unique(unlist(l1))
l2 <- unique(unlist(l2))
l3 <- unique(unlist(l3))
l4 <- unique(unlist(l4))
l5 <- unique(unlist(l5))
l6 <- unique(unlist(l6))

l1 <- l1[1:floor(0.9 * length(l1))]
l2 <- l2[1:floor(0.9 * length(l2))]
l3 <- l3[1:floor(0.9 * length(l3))]
l4 <- l4[1:floor(0.3 * length(l4))]
l5 <- l5[1:floor(0.9 * length(l5))]
l6 <- l6[1:floor(0.3 * length(l6))]

data <- tibble::tibble(
  batch = c("b1", "b1", "b2", "b2", "b1", "b1", "b2", "b2"),
  group = c("g1", "g2", "g1", "g2", "g1", "g2", "g1", "g2"),
  type = c("upregulated", "upregulated", "upregulated", "upregulated",
           "downregulated", "downregulated", "downregulated", "downregulated"),
  data = list(
    tibble::tibble(
      uniprot = unique(c(l1, l2, l4)),
      symbol = unique(c(l1, l2, l4))
    ),
    tibble::tibble(
      uniprot = unique(c(l1, l3, l4)),
      symbol = unique(c(l1, l3, l4))
    ),
    tibble::tibble(
      uniprot = unique(c(l2, l4)),
      symbol = unique(c(l2, l4))
    ),
    tibble::tibble(
      uniprot = unique(c(l3, l4)),
      symbol = unique(c(l3, l4))
    ),
    tibble::tibble(
      uniprot = unique(c(l5, l6)),
      symbol = unique(c(l5, l6))
    ),
    tibble::tibble(
      uniprot = unique(c(l5, l6)),
      symbol = unique(c(l5, l6))
    ),
    tibble::tibble(
      uniprot = unique(c(l5, l6)),
      symbol = unique(c(l5, l6))
    ),
    tibble::tibble(
      uniprot = unique(c(l5, l6)),
      symbol = unique(c(l5, l6))
    )
  )
)


## ----echo=TRUE----------------------------------------------------------------

ann_sources <- c(
  GO = system.file(
    "extdata",
    paste("human_functional_annotations",
      "go", "data.tab", sep = "/"),
    package = "fctBio"
  ),
  RC = system.file(
    "extdata",
    paste("human_functional_annotations",
      "reactome", "data.tab", sep = "/"),
    package = "fctBio"
  )
)

ann_types <- c(
  GO = c("biological_process")
)

ann_space <- fctBio::load_ann_space(ann_sources, ann_types)

## ----echo=TRUE----------------------------------------------------------------
knitr::kable(head(ann_space,n = 5)[,c("term","name","parents","type","ann_name")])

## -----------------------------------------------------------------------------
str(head(ann_space,n = 2)[,c("genes","parent_genes")])

## -----------------------------------------------------------------------------
knitr::kable(head(
  dplyr::filter(ann_space, .data$ann_name == "RC")
  ,n = 5)[,c("term","name","parents","type","ann_name")])

## -----------------------------------------------------------------------------
test_enrich <- fctBio::NestedEnrich$new(data, ann_space)

## -----------------------------------------------------------------------------
test_enrich$filter_and_set_significant_results(
  0.05, p_type = "qval_bonferroni", min_signif_term_for_clust = Inf)

## -----------------------------------------------------------------------------
test_enrich$build_and_set_i_matrix()
test_enrich$build_and_set_p_matrix()


## -----------------------------------------------------------------------------
head_imatrix <- data.frame(test_enrich$get_i_matrix()[1:20,1:30])

head_imatrix$Row <- rownames(head_imatrix)

# Transform to long format
head_imatrix <- tidyr::pivot_longer(head_imatrix, 
                        cols = -Row, 
                        names_to = "Column", 
                        values_to = "Value")

# Plot using ggplot2
ggplot2::ggplot(head_imatrix, ggplot2::aes(x = Column, y = Row, fill = as.factor(Value))) +
  ggplot2::geom_tile(color = "white") + 
  ggplot2::scale_fill_manual(values = c("0" = "white", "1" = "black")) +
  ggplot2::theme_minimal() +
  ggplot2::guides(fill="none") +
  ggplot2::labs(x = "Columns", y = "Rows", fill = "Value") +
  ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))

## -----------------------------------------------------------------------------
head_pmatrix <- data.frame(test_enrich$get_p_matrix()[1:20, ])

head_pmatrix$Row <- rownames(head_pmatrix)

# Transform to long format
head_pmatrix <- tidyr::pivot_longer(head_pmatrix, 
                        cols = -Row, 
                        names_to = "Column", 
                        values_to = "Value")
head_pmatrix$Value[is.na(head_pmatrix$Value)] <- 1
# Plot using ggplot2
ggplot2::ggplot(head_pmatrix, ggplot2::aes(x = Column, y = Row, fill = -log10(Value))) +
  ggplot2::geom_tile(color = "white") + 
  ggplot2::scale_fill_gradient(high = "black", low = "white") +
  ggplot2::theme_minimal() +
  ggplot2::guides(fill=FALSE) +
  ggplot2::labs(x = "Columns", y = "Rows", fill = "Value") +
  ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))



## -----------------------------------------------------------------------------
test_enrich$build_and_set_hclust(method_dist = "jaccard")
test_enrich$plot_hclust()

## -----------------------------------------------------------------------------
test_enrich$plot_combined(max_term_per_cluster = 4)

## -----------------------------------------------------------------------------
test_enrich$build_and_set_hclust(method_dist = "jaccard", matrix_type = "pval")
test_enrich$plot_hclust()

## -----------------------------------------------------------------------------
test_enrich$plot_combined(max_term_per_cluster = 4)

## -----------------------------------------------------------------------------
test_enrich$build_and_set_hclust(method_dist = "jaccard",using_pvclust = T, parallel = 12L)

## -----------------------------------------------------------------------------
hc <- test_enrich$get_hclust()
plot(hc)

## -----------------------------------------------------------------------------
test_enrich$cut_hclust_and_set_clusters(
  0.9,
  min_size = 2,
  max_size = 8,
  rm_redundancy_method = "largest"
  )
test_enrich$plot_hclust()

## -----------------------------------------------------------------------------
sum_table <- head(test_enrich$build_cluster_summary_table())
names(sum_table) <- stringr::str_replace_all(names(sum_table), "\\|", ":")
kableExtra::scroll_box(knitr::kable(sum_table, digits = 10), width = "100%")


## -----------------------------------------------------------------------------
ct_table <- head(test_enrich$count_gene_per_cluster())
knitr::kable(ct_table)

