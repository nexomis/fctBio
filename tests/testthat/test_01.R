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

test_enrich <- fctBio::NestedEnrich$new(data, ann_space)

test_enrich$filter_and_set_significant_results(
  0.05, p_type = "qval_bonferroni", min_signif_term_for_clust = Inf)

test_that("plot without clustering", {
  vdiffr::expect_doppelganger("default",
    test_enrich$plot_combined())
  vdiffr::expect_doppelganger("with max",
    test_enrich$plot_combined(
      max_cluster = 5,
      max_term_per_cluster = 10
  ))
  vdiffr::expect_doppelganger("size ordered",
    test_enrich$plot_combined(
      ordered_by_pval = FALSE
  ))
  vdiffr::expect_doppelganger("keep unsignif",
    test_enrich$plot_combined(
      keep_only_signif = FALSE
  ))
  vdiffr::expect_doppelganger("in_ filters",
    test_enrich$plot_combined(
      in_batch = c("b1"),
      in_group = c("g1"),
      in_type = c("upregulated", "downregulated")
  ))
  vdiffr::expect_doppelganger("labels id",
    test_enrich$plot_combined(
      to = "id"
  ))
  vdiffr::expect_doppelganger("reduced_label size 20",
    test_enrich$plot_combined(
      to = "reduced_label",
      reduced_label_max_size = 20
  ))
  vdiffr::expect_doppelganger("reduced_label_with_id size 5",
    test_enrich$plot_combined(
      to = "reduced_label_with_id",
      reduced_label_max_size = 5
  ))
  vdiffr::expect_doppelganger("labels full_name",
    test_enrich$plot_combined(
      to = "full_name"
  ))
  test_enrich$set_batch_labels(c(b2 = "batch 2", b1 = "batch 1"))
  test_enrich$set_group_labels(c(g2 = "group 2", g1 = "group 1"))
  vdiffr::expect_doppelganger("batch group labelled",
    test_enrich$plot_combined())
})

library(magrittr)

results_raw <- test_enrich$get_results()
results <- test_enrich$get_results(verbose = TRUE)

testthat::test_that("test enrich structure", {
  testthat::expect_identical(
    c("batch", "group", "type", "gene_inputs", "ann_name", "mapping", "enrich",
      "universe"),
    names(results)
  )
  testthat::expect_identical(
    c("batch", "group", "type", "gene_inputs", "ann_name", "mapping", "enrich",
      "universe"),
    names(results_raw)
  )
  testthat::expect_identical(
    c("factor", "factor", "factor", "list", "factor", "list", "list",
      "list"),
    sapply(names(results), function(n) class(results[[n]]), USE.NAMES = FALSE)
  )
  testthat::expect_identical(
    c("factor", "factor", "factor", "list", "factor", "list", "list",
      "list"),
    sapply(names(results_raw), function(n) class(results_raw[[n]]),
      USE.NAMES = FALSE)
  )
  testthat::expect_identical(
    "character",
    class(results[["gene_inputs"]][[1]])
  )
  testthat::expect_identical(
    "character",
    class(results[["universe"]][[1]])
  )
  testthat::expect_identical(
    "integer",
    class(results_raw[["gene_inputs"]][[1]])
  )
  testthat::expect_identical(
    "integer",
    class(results_raw[["universe"]][[1]])
  )
  testthat::expect_identical(
    c("term", "name", "type", "genes", "parents", "GI", "NGI", "Ngenes", "parent_NGI",
      "parent_Ngenes", "pmin", "pval", "qval_bonferroni", "qval_bh"),
    names(results[["enrich"]][[1]])
  )
  testthat::expect_identical(
    c("term", "genes", "parents", "GI", "NGI", "Ngenes", "parent_NGI",
      "parent_Ngenes", "pmin", "pval", "qval_bonferroni", "qval_bh"),
    names(results_raw[["enrich"]][[1]])
  )
  testthat::expect_identical(
    c("character", "character", "character", "list", "list", "list", "integer", "integer", "integer",
      "integer", "numeric", "numeric", "numeric", "numeric"),
    sapply(
      names(results[["enrich"]][[1]]),
      function(n) class(results[["enrich"]][[1]][[n]]),
      USE.NAMES = FALSE)
  )
  testthat::expect_identical(
    c("integer", "list", "list", "list", "integer", "integer", "integer",
      "integer", "numeric", "numeric", "numeric", "numeric"),
    sapply(
      names(results_raw[["enrich"]][[1]]),
      function(n) class(results_raw[["enrich"]][[1]][[n]]),
      USE.NAMES = FALSE)
  )
  testthat::expect_identical(
    "character",
    class(results[["enrich"]][[1]][["genes"]][[1]])
  )
  testthat::expect_identical(
    "character",
    class(results[["enrich"]][[1]][["parents"]][[1]])
  )
  testthat::expect_identical(
    "character",
    class(results[["enrich"]][[1]][["GI"]][[1]])
  )
    testthat::expect_identical(
    "integer",
    class(results_raw[["enrich"]][[1]][["genes"]][[1]])
  )
  testthat::expect_identical(
    "integer",
    class(results_raw[["enrich"]][[1]][["parents"]][[1]])
  )
  testthat::expect_identical(
    "integer",
    class(results_raw[["enrich"]][[1]][["GI"]][[1]])
  )
})

known_terms <- list(
  l1 = c("GO:0034340", "GO:0140374"),
  l2 = c("GO:0002347"),
  l4 = tab_base[c(20, 30, 40), "term"], # LOW
  l3 = c("GO:0070914"),
  l5 = c("GO:0042246"),
  l6 = tab_base[c(25, 35, 45), "term"]  # LOW
)

b1g1up <- test_enrich$filter_and_get_results("b1", "g1", "upregulated", "GO"
  )$enrich[[1]]

b1g1up_signif <- b1g1up %>%
  dplyr::filter(
    .data$term %in% unlist(known_terms[c("l1", "l2")])
  )

b1g2up <- test_enrich$filter_and_get_results("b1", "g2", "upregulated", "GO"
  )$enrich[[1]]

b1g2up_signif <- b1g2up %>%
  dplyr::filter(
    .data$term %in% unlist(known_terms[c("l1", "l3")])
  )

b2g1up <- test_enrich$filter_and_get_results("b2", "g1", "upregulated", "GO"
  )$enrich[[1]]

b2g1up_signif <- b2g1up %>%
  dplyr::filter(
    .data$term %in% unlist(known_terms[c("l2")])
  )

b2g2up <- test_enrich$filter_and_get_results("b2", "g2", "upregulated", "GO"
  )$enrich[[1]]

b2g2up_signif <- b2g2up %>%
  dplyr::filter(
    .data$term %in% unlist(known_terms[c("l3")])
  )

b1g1down <- test_enrich$filter_and_get_results("b1", "g1", "downregulated", "GO"
  )$enrich[[1]]

b1g1down_signif <- b1g1down %>%
  dplyr::filter(
    .data$term %in% unlist(known_terms[c("l5")])
  )

b1g1up_rc <- test_enrich$filter_and_get_results("b1", "g1", "upregulated", "RC"
  )$enrich[[1]]

testthat::test_that("test enrich expectations", {
  testthat::expect_gt(
    nrow(b1g1up),
    100
  )
  testthat::expect_gt(
    nrow(b1g2up),
    100
  )
  testthat::expect_gt(
    nrow(b2g1up),
    100
  )
  testthat::expect_gt(
    nrow(b2g2up),
    100
  )
  testthat::expect_gt(
    nrow(b1g1down),
    100
  )
  testthat::expect_gt(
    nrow(b1g1up_rc),
    20
  )
  testthat::expect_true(
    all(is.integer(b1g1up$NGI))
  )
  testthat::expect_true(
    all(is.integer(b1g1up$Ngenes))
  )
  testthat::expect_true(
    all(is.integer(b1g1up$parent_NGI))
  )
  testthat::expect_true(
    all(is.integer(b1g1up$parent_Ngenes))
  )
  testthat::expect_true(
    all(is.numeric(b1g1up$pmin))
  )
  testthat::expect_true(
    all(is.numeric(b1g1up$pval))
  )
  testthat::expect_true(
    all(! is.na(b1g1up$term))
  )
  testthat::expect_true(
    all(! is.na(b1g1up$type))
  )
  testthat::expect_true(
    all(! is.na(b1g1up$pmin))
  )
  testthat::expect_true(
    all(! is.na(b1g1up$name))
  )
  testthat::expect_true(
    all(! is.na(b1g1up$pval))
  )
})

testthat::test_that("test enrich pvals", {
  testthat::expect_lt(
    sum(b1g1up_signif$pval),
    0.01
  )
  testthat::expect_lt(
    sum(b1g2up_signif$pval),
    0.01
  )
  testthat::expect_lt(
    sum(b2g1up_signif$pval),
    0.01
  )
  testthat::expect_lt(
    sum(b2g2up_signif$pval),
    0.01
  )
  testthat::expect_lt(
    sum(b1g1down_signif$pval),
    0.01
  )
})

signif_results <- test_enrich$get_significant_results()

test_enrich$filter_and_set_significant_results(
  0.05, p_type = "qval_bonferroni", build_and_set_i_matrix = FALSE)
signif_results_clone <- test_enrich$get_significant_results()
imat_null <- test_enrich$get_i_matrix()
hclust_null <- test_enrich$get_hclust()
clusters_null <- test_enrich$get_clusters()

test_enrich$filter_and_set_significant_results(0.05)
imat <- test_enrich$get_i_matrix()
p_1 <- sum(imat) / (nrow(imat) * ncol(imat))

testthat::test_that("test incidence matrix", {

  testthat::expect_identical(
    signif_results,
    signif_results_clone
  )

  testthat::expect_null(imat_null)
  testthat::expect_null(hclust_null)
  testthat::expect_null(clusters_null)

  testthat::expect_true(
    all(imat %in% c(0, 1))
  )
  testthat::expect_true(
    all(matrixStats::rowSums2(imat) > 0)
  )

  testthat::expect_false(
    is.null(rownames(imat))
  )

  testthat::expect_false(
    is.null(colnames(imat))
  )

  testthat::expect_gt(
    ncol(imat),
    12000
  )

  testthat::expect_gt(
    nrow(imat),
    100
  )

  testthat::expect_gt(
    p_1,
    0.001
  )

  testthat::expect_lt(
    p_1,
    0.15
  )

  vdiffr::expect_doppelganger("default cluster",
    test_enrich$plot_hclust()
  )

  test_enrich$cut_hclust_and_set_clusters(6)

  vdiffr::expect_doppelganger("hclust cut 6",
    test_enrich$plot_hclust()
  )

  vdiffr::expect_doppelganger("plot combined hclust cut 6",
    test_enrich$plot_combined()
  )

})

cgene_count = test_enrich$count_gene_per_cluster()
genes <- cgene_count$gene
gene2Test <- paste0(genes, "test")
names(gene2Test) <- genes
cgene_count_verbose = test_enrich$count_gene_per_cluster(id2name <- gene2Test,
  verbose = TRUE)


expected_cols <- c("gene", "cluster", "intra_x_term", "extra_x_term",
  "x_cluster", "x_batch", "x_group")
expected_cols_verbose <- c("Gene", "Gene Name", "Cluster", "Intra", "Extra", "Clusters",
  "Batches", "Groups")

testthat::test_that("test cgenecount matrix", {

  testthat::expect_identical(
    names(cgene_count),
    expected_cols
  )
  testthat::expect_identical(
    names(cgene_count_verbose),
    expected_cols_verbose
  )

  testthat::expect_identical(
    sapply(names(cgene_count), function(n) class(cgene_count[[n]]),
      USE.NAMES = FALSE),
    c("character","integer","integer","integer","integer","integer","integer")
  )

  expect_true(all(grepl("test$", cgene_count_verbose[["Gene Name"]])))

})

sumt <- test_enrich$build_cluster_summary_table()
sumt_verbose <- test_enrich$build_cluster_summary_table(verbose = TRUE)

temp_dir <- tempdir()

test_enrich$write_xlsx(temp_dir)
file_path <- file.path(temp_dir, "cluster_summary_fctBio.xlsx")
test_that("XLSX file is created", {
  expect_true(file.exists(file_path))
})


# all values below should be character vector of size > 0
ch_batch <- fctBio::query_genes_choices(results,"batch",cgene_count)
ch_batch_label <- fctBio::query_genes_choices(results,"batch_label",cgene_count, batch_labels = c(b2 = "batch 2", b1 = "batch 1"))
ch_group <- fctBio::query_genes_choices(results,"group",cgene_count)
ch_group_label <- fctBio::query_genes_choices(results,"group_label",cgene_count, group_labels = c(g2 = "group 2", g1 = "group 1"))
ch_term_id <- fctBio::query_genes_choices(results,"term_id",cgene_count)
ch_term_name <- fctBio::query_genes_choices(results,"term_name",cgene_count)
ch_cluster <- fctBio::query_genes_choices(results,"cluster",cgene_count)

# all values below should be integer vector of size > 0
ch_intra_occurences <- fctBio::query_genes_choices(results,"intra_occurences",cgene_count)
ch_extra_occurences <- fctBio::query_genes_choices(results,"extra_occurences",cgene_count)
ch_cluster_occurences <- fctBio::query_genes_choices(results,"cluster_occurences",cgene_count)

# all values below should be character vector of size > 0
q_batch <- fctBio::query_genes(results,"batch", "b1",cgene_count)
q_batch_label <- fctBio::query_genes(results,"batch_label", "batch 1",cgene_count, batch_labels = c(b2 = "batch 2", b1 = "batch 1"))
q_group <- fctBio::query_genes(results,"group", "g1",cgene_count)
q_group_label <- fctBio::query_genes(results,"group_label", "group 1",cgene_count, group_labels = c(g2 = "group 2", g1 = "group 1"))
q_term_id <- fctBio::query_genes(results,"term_id", "GO:0140374", cgene_count)
q_cluster <- fctBio::query_genes(results,"cluster", "1", cgene_count)
q_intra_occurences <- fctBio::query_genes(results,"intra_occurences", 1, cgene_count)
q_extra_occurences <- fctBio::query_genes(results,"extra_occurences", 1, cgene_count)
q_cluster_occurences <- fctBio::query_genes(results,"cluster_occurences", 1, cgene_count)

test_that("Check value type and size", {
  # Check for character vectors of size > 0
  expect_true(is.factor(ch_batch) && length(ch_batch) > 0)
  expect_true(is.character(ch_batch_label) && length(ch_batch_label) > 0)
  expect_true(is.factor(ch_group) && length(ch_group) > 0)
  expect_true(is.character(ch_group_label) && length(ch_group_label) > 0)
  expect_true(is.character(ch_term_id) && length(ch_term_id) > 0)
  expect_true(is.character(ch_term_name) && length(ch_term_name) > 0)
  expect_true(is.integer(ch_cluster) && length(ch_cluster) > 0)
  
  # Check for integer vectors of size > 0
  expect_true(is.integer(ch_intra_occurences) && length(ch_intra_occurences) > 0)
  expect_true(is.integer(ch_extra_occurences) && length(ch_extra_occurences) > 0)
  expect_true(is.integer(ch_cluster_occurences) && length(ch_cluster_occurences) > 0)

  # Check for character vectors with specific values and size > 0
  expect_true(is.character(q_batch) && length(q_batch) > 0)
  expect_true(is.character(q_batch_label) && length(q_batch_label) > 0)
  expect_true(is.character(q_group) && length(q_group) > 0)
  expect_true(is.character(q_group_label) && length(q_group_label) > 0)
  expect_true(is.character(q_term_id) && length(q_term_id) > 0)
  expect_true(is.character(q_cluster) && length(q_cluster) > 0)
  expect_true(is.character(q_intra_occurences) && length(q_intra_occurences) > 0)
  expect_true(is.character(q_extra_occurences) && length(q_extra_occurences) > 0)
  expect_true(is.character(q_cluster_occurences) && length(q_cluster_occurences) > 0)
})

## test pval matrix
test_enrich$build_and_set_p_matrix()
test_enrich$build_and_set_hclust(method_dist = "jaccard", matrix_type = "pval")

