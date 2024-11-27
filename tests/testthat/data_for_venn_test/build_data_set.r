# script for building data (input and reference) for unitary test of 'venn_genes' functions


main_test_dir="tests/testthat/data_for_venn_test/"
dir.create(main_test_dir)

deg_table_base <- data.table(
  batch = c("Batch1", "Batch1", "Batch2", "Batch2", "Batch2", "Batch3", "Batch3", "Batch3", "Batch3"),
  group = c("GrpA", "GrpB", "GrpA", "GrpC", "GrpD", "GrpD", "GrpE", "GrpF", "GrpG"),
  lfc_abs_lim = 0.5849625,
  min_signif = c(0.05, 0.05, 0.01, 0.05, 0.05, 0.05, 0.05, 0.01, 0.01),
  type = c("deregulated", "deregulated", "downregulated", "downregulated", "downregulated", "downregulated", "deregulated", "downregulated", "downregulated"),
  data = list(
    data.table(symbol = paste0("GENE", seq(1, 10)),
               uniprot = paste0("PROT", seq(1, 10))),
    data.table(symbol = paste0("GENE", seq(3, 12)),
               uniprot = paste0("PROT", seq(3, 12))),
    data.table(symbol = paste0("GENE", seq(1, 10)),
               uniprot = paste0("PROT", seq(1, 10))),
    data.table(symbol = paste0("GENE", seq(4, 9)),
               uniprot = paste0("PROT", seq(4, 9))),
    data.table(symbol = paste0("GENE", seq(9, 15)),
               uniprot = paste0("PROT", seq(9, 15))),
    data.table(symbol = paste0("GENE", seq(6, 20)),
               uniprot = paste0("PROT", seq(6, 20))),
    data.table(symbol = paste0("GENE", seq(40, 50)),
               uniprot = paste0("PROT", seq(40, 50))),
    data.table(symbol = paste0("GENE", seq(4, 15)),
               uniprot = paste0("PROT", seq(4, 15))),
    data.table(symbol = paste0("GENE", seq(7, 15)),
               uniprot = paste0("PROT", seq(7, 15)))
  )
)

## 01
test_id="01"
dir.create(paste0(main_test_dir, "/", test_id))

deg_table <- deg_table_base
to_keep <- NULL
ext_groups <- NULL

exec_args <- list(deg_table_in = deg_table,
  min_to_displayed_threshold = 2,
  to_keep = to_keep,
  ext_groups = ext_groups,
  mode_ext_groups = "union",
  inter_batch = NULL,
  mode_inter_batch = "union",
  shape = "ellipse",
  path_dir_to_save_plot =  NULL,
  dict_uniprot_to_symbol = setNames(paste0("GENE", 0:100), paste0("PROT", 0:100)))
                   
res <- do.call(compute_euler_plot, exec_args)

saveRDS(deg_table,
        file = paste0(main_test_dir, "/", test_id, "/deg_table_input.rds"))
saveRDS(to_keep,
        file = paste0(main_test_dir, "/", test_id, "/to_keep_input.rds"))
saveRDS(ext_groups,
        file = paste0(main_test_dir, "/", test_id, "/ext_groups_input.rds"))

saveRDS(exec_args,
        file = paste0(main_test_dir, "/", test_id, "/execution_args.rds"))

saveRDS(res$dt_region_content,
        file = paste0(main_test_dir, "/", test_id, "/result_dt_output.rds"))
saveRDS(res$euler_plots,
        file = paste0(main_test_dir, "/", test_id, "/result_plot_output.rds"))

## 02
test_id="02"
dir.create(paste0(main_test_dir, "/", test_id))

deg_table <- deg_table_base
to_keep <- data.table(
  batch = c("Batch_XXXX", "Batch1", "Batch3"),
  group = c("Grp_XXXX", "GrpB", "GrpE")
)
ext_groups <- data.table(
  batch = c("Batch1", "Batch1", "Batch3", "Batch3", "Batch3"),
  group = c("GrpA", "GrpB", "GrpD", "GrpF", "GrpG")
)

exec_args <- list(deg_table_in = deg_table,
  min_to_displayed_threshold = 0.01,
  to_keep = to_keep,
  ext_groups = ext_groups,
  mode_ext_groups = "union",
  inter_batch = NULL,
  mode_inter_batch = "union",
  shape = "ellipse",
  path_dir_to_save_plot =  NULL,
  dict_uniprot_to_symbol = setNames(paste0("GENE", 0:100), paste0("PROT", 0:100)))

res <- do.call(compute_euler_plot, exec_args)

saveRDS(deg_table,
        file = paste0(main_test_dir, "/", test_id, "/deg_table_input.rds"))
saveRDS(to_keep,
        file = paste0(main_test_dir, "/", test_id, "/to_keep_input.rds"))
saveRDS(ext_groups,
        file = paste0(main_test_dir, "/", test_id, "/ext_groups_input.rds"))

saveRDS(exec_args,
        file = paste0(main_test_dir, "/", test_id, "/execution_args.rds"))

saveRDS(res$dt_region_content,
        file = paste0(main_test_dir, "/", test_id, "/result_dt_output.rds"))
saveRDS(res$euler_plots,
        file = paste0(main_test_dir, "/", test_id, "/result_plot_output.rds"))

## 03
test_id="03"
dir.create(paste0(main_test_dir, "/", test_id))

deg_table <- deg_table_base
to_keep <- NULL
ext_groups <- c("PROT3", "PROT6", "PROT9", "PROT12", "PROT15", "PROT18", "PROT21", "PROT41")

exec_args <- list(deg_table_in = deg_table,
                  min_to_displayed_threshold = 0,
                  to_keep = to_keep,
                  ext_groups = ext_groups,
                  mode_ext_groups = "union",
                  inter_batch = NULL,
                  mode_inter_batch = "intersect",
                  shape = "ellipse",
                  path_dir_to_save_plot =  NULL)

res <- do.call(compute_euler_plot, exec_args)

saveRDS(deg_table,
        file = paste0(main_test_dir, "/", test_id, "/deg_table_input.rds"))
saveRDS(to_keep,
        file = paste0(main_test_dir, "/", test_id, "/to_keep_input.rds"))
saveRDS(ext_groups,
        file = paste0(main_test_dir, "/", test_id, "/ext_groups_input.rds"))

saveRDS(exec_args,
        file = paste0(main_test_dir, "/", test_id, "/execution_args.rds"))

saveRDS(res$dt_region_content,
        file = paste0(main_test_dir, "/", test_id, "/result_dt_output.rds"))
saveRDS(res$euler_plots,
        file = paste0(main_test_dir, "/", test_id, "/result_plot_output.rds"))

## 04
test_id="04"
dir.create(paste0(main_test_dir, "/", test_id))

deg_table <- deg_table_base
to_keep <- NULL
ext_groups <- NULL

exec_args <- list(deg_table_in = deg_table,
                  min_to_displayed_threshold = 0,
                  to_keep = to_keep,
                  ext_groups = ext_groups,
                  mode_ext_groups = "union",
                  inter_batch = list("b_1_2_uni" = list("Batch1", "Batch2"),
                                     "b_1_2_inter" = list("Batch1", "Batch2"),
                                     "b_1_2_uni" = list("Batch1", "Batch3")),
                  mode_inter_batch = list("b_1_2_uni" = "union",
                                          "b_1_2_inter" = "intersect",
                                          "b_1_2_uni" = "union"),
                  shape = "ellipse",
                  path_dir_to_save_plot =  NULL)

res <- do.call(compute_euler_plot, exec_args)

saveRDS(deg_table,
        file = paste0(main_test_dir, "/", test_id, "/deg_table_input.rds"))
saveRDS(to_keep,
        file = paste0(main_test_dir, "/", test_id, "/to_keep_input.rds"))
saveRDS(ext_groups,
        file = paste0(main_test_dir, "/", test_id, "/ext_groups_input.rds"))

saveRDS(exec_args,
        file = paste0(main_test_dir, "/", test_id, "/execution_args.rds"))

saveRDS(res$dt_region_content,
        file = paste0(main_test_dir, "/", test_id, "/result_dt_output.rds"))
saveRDS(res$euler_plots,
        file = paste0(main_test_dir, "/", test_id, "/result_plot_output.rds"))

