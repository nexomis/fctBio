#' @include utils.r
#' @include globals.r

NULL

#' Draw Single Euler Plot
#'
#' @description
#' Creates a simple Euler diagram for visualizing overlaps between sets.
#' This function generates an Euler plot with customizable display thresholds
#' and styling options, with optional saving to file.
#'
#' @param sets A named list of character vectors representing the sets to compare.
#' @param min_to_displayed_threshold Numeric; minimal values to display on the plot.
#' If less than 1, interpreted as a proportion of the universe (total unique elements
#' across all sets). If 1 or greater, interpreted as an absolute count. Default is 0.01.
#' @param shape Character; shape of the Euler diagram. Default is "ellipse".
#' @param title Character; title for the plot. Default is NULL.
#' @param save_path_file Character; file path to save the plot. If NULL, plot is
#' not saved. Default is NULL.
#'
#' @return An Euler plot object that can be displayed or further customized.
#'
#' @examples
#' sets <- list(
#'   "Set A" = c("gene1", "gene2", "gene3"),
#'   "Set B" = c("gene2", "gene3", "gene4"),
#'   "Set C" = c("gene3", "gene4", "gene5")
#' )
#' plot <- draw_single_euler_plot(sets, title = "Gene Overlap")
#'
#' @export
draw_single_euler_plot <- function(
  sets,
  min_to_displayed_threshold = 0.01,
  shape = "ellipse",
  title = NULL,
  save_path_file = NULL
) {

  if (!is.list(sets) || is.null(names(sets))) {
    stop("sets must be a named list", call. = FALSE)
  }

  euler_data <- eulerr::euler(sets, shape = shape)

  if (min_to_displayed_threshold < 1L) {
    min_to_displayed <- round(
      length(unique(unlist(sets))) * min_to_displayed_threshold
    )
  } else {
    min_to_displayed <- min_to_displayed_threshold
  }

  qts <- euler_data$original.values
  qts[qts < min_to_displayed] <- ""

  euler_diagram <- plot(
    euler_data,
    quantities = qts,
    lty = 1L:6L,
    edges = "black",
    main = paste0(title, " (min.val: ", min_to_displayed, ")")
  )

  if (!is.null(save_path_file)) {
    ggplot2::ggsave(
      filename = save_path_file,
      plot = euler_diagram,
      width = 11.69,
      height = 8.27
    )
  }

  return(euler_diagram)
}

#' Extract Venn Diagram Zone Content
#'
#' @description
#' Extracts the specific content of each zone in a Venn/Euler diagram,
#' identifying which elements belong exclusively to specific combinations
#' of sets. This function can process single sets or lists of sets and
#' optionally save results to an Excel file.
#'
#' @param sets_list A named list of sets or a list of named lists of sets.
#' Each list will be processed independently to extract zone-specific content.
#' @param sort_group_for_region_naming Logical; whether to alphabetically sort
#' group names when creating region names. If TRUE, sorts group names
#' alphabetically. If FALSE, preserves the original order. Default is TRUE.
#' @param path_xlsx_file_to_save_table Character; file path to save results
#' in an Excel file. If provided, creates one sheet per set. Default is NULL.
#'
#' @return A named list of results for each input set. Each result contains
#' the specific content of each zone as character vectors.
#'
#' @examples
#' set123 <- list(
#'   "Group1" = c("gene1", "gene2", "gene3"),
#'   "Group2" = c("gene2", "gene4"),
#'   "Group3" = c("gene1", "gene4", "gene5")
#' )
#' setABC <- list(
#'   "GroupA" = c("geneA", "geneB", "geneC"),
#'   "GroupB" = c("geneB", "geneC"),
#'   "GroupC" = c("geneC", "geneD", "geneE")
#' )
#' results <- extract_venn_zones_content(
#'   list("Set123" = set123, "SetABC" = setABC)
#' )
#'
#' @export
extract_venn_zones_content <- function(
  sets_list,
  sort_group_for_region_naming = TRUE,
  path_xlsx_file_to_save_table = NULL
) {

  # Convert single set to list format for consistent processing
  if (!is.list(sets_list[[1L]])) {
    sets_list <- list("specific_intersect" = sets_list)
  }

  all_results <- list()

  # Process each set
  for (set_name in names(sets_list)) {
    sets <- sets_list[[set_name]]

    if (!is.list(sets) || is.null(names(sets))) {
      warning("Skipping set '", set_name, "': must be a named list", call. = FALSE)
      next
    }

    # Generate region names for all possible combinations
    if (sort_group_for_region_naming) {
      region_names <- unlist(lapply(seq_along(names(sets)), function(x) {
        utils::combn(sort(names(sets)), x, paste, collapse = "&")
      }))
    } else {
      region_names <- unlist(lapply(seq_along(names(sets)), function(x) {
        utils::combn(names(sets), x, paste, collapse = "&")
      }))
    }

    regions_content <- list()

    # Extract specific content for each region
    for (region in region_names) {
      incl_sets <- unlist(strsplit(region, "&"))
      excl_sets <- setdiff(names(sets), incl_sets)

      # Find intersection of included sets
      region_content <- as.character(na.omit(Reduce(intersect, sets[incl_sets])))

      # Remove elements present in excluded sets
      if (length(excl_sets) > 0L) {
        region_content <- setdiff(region_content, unlist(sets[excl_sets]))
      }

      regions_content[[region]] <- region_content
    }

    all_results[[set_name]] <- regions_content
  }

  # Save results to Excel file if requested
  if (!is.null(path_xlsx_file_to_save_table)) {
    if (!requireNamespace("openxlsx", quietly = TRUE)) {
      warning("openxlsx package required for Excel export", call. = FALSE)
    } else {
      wb <- openxlsx::createWorkbook()
      for (set_name in names(all_results)) {
        openxlsx::addWorksheet(wb, set_name)

        # Convert results to table format
        max_length <- max(sapply(all_results[[set_name]], length))
        region_dt <- as.data.frame(
          do.call(cbind, lapply(all_results[[set_name]], function(x) {
            if (length(x) == 0L) return(rep(NA, max_length))
            c(x, rep(NA, max_length - length(x)))
          }))
        )
        openxlsx::writeData(wb, set_name, region_dt, rowNames = FALSE)
      }
      openxlsx::saveWorkbook(wb, path_xlsx_file_to_save_table, overwrite = TRUE)
    }
  }

  return(all_results)
}

#' Get External Group UniProt IDs
#'
#' @description
#' Extracts UniProt IDs from a differential expression table based on
#' specified batch and group combinations. Supports different modes
#' for combining multiple groups.
#'
#' @param deg_dt A data.table with differential expression results containing
#' columns: batch, group, type, and uniprot. The 'type' column must not contain
#' NA values.
#' @param ext_groups_dt A data.table or data.frame with columns 'batch' and 'group'
#' specifying which groups to extract.
#' @param mode_ext_groups Character; mode for combining groups. Must be one of:
#' "individual" (separate list for each group), "union" (combine all groups),
#' or "intersect" (intersection of all groups).
#' @param deg_dt_name Character; name of the deg_dt for warning messages.
#' Default is NULL.
#'
#' @return A named list of character vectors containing UniProt IDs.
#'
#' @examples
#' \dontrun{
#' deg_data <- data.table(
#'   batch = c("B1", "B1", "B2"),
#'   group = c("G1", "G2", "G1"),
#'   type = c("up", "down", "up"),
#'   uniprot = list(c("P1", "P2"), c("P2", "P3"), c("P1", "P4"))
#' )
#' ext_groups <- data.table(batch = "B1", group = c("G1", "G2"))
#' result <- get_ext_group(deg_data, ext_groups, "union")
#' }
#'
#' @keywords internal
get_ext_group <- function(
  deg_dt,
  ext_groups_dt,
  mode_ext_groups,
  deg_dt_name = NULL
) {

  if (!mode_ext_groups %in% c("individual", "union", "intersect")) {
    stop("mode_ext_groups must be 'individual', 'union', or 'intersect'", call. = FALSE)
  }

  deg_dt_ext <- deg_dt[ext_groups_dt, on = .(batch, group)]

  # Check for missing values
  missing_rows <- deg_dt_ext[is.na(type)]
  if (nrow(missing_rows) > 0L) {
    warning(
      "Some 'ext_groups' values were not found in 'deg_table' for analysis of ",
      deg_dt_name %||% "unnamed dataset", call. = FALSE
    )
    print(missing_rows[, .(batch, group)])
    deg_dt_ext <- deg_dt_ext[!is.na(type)]
  }

  switch(mode_ext_groups,
    individual = {
      ext_sets <- lapply(deg_dt_ext$uniprot, function(uniprot_ids) {
        unique(na.omit(unlist(uniprot_ids)))
      })
      names(ext_sets) <- paste("EXT", deg_dt_ext$batch, deg_dt_ext$group, sep = ":")
      ext_sets
    },
    union = {
      ext_sets <- list(unique(na.omit(unlist(deg_dt_ext$uniprot))))
      names(ext_sets) <- "EXT:union"
      ext_sets
    },
    intersect = {
      ext_sets <- list(as.character(na.omit(Reduce(intersect, deg_dt_ext$uniprot))))
      names(ext_sets) <- "EXT:intersect"
      ext_sets
    }
  )
}

#' Compute and Draw Euler Plot from DEG Table
#'
#' @description
#' Computes and draws Euler plots based on significant differential expression
#' gene (DEG) tables. The analysis is performed at the UniProt level, with
#' NA values ignored. The function supports both intra-batch and inter-batch
#' comparisons with various customization options.
#'
#' @param deg_table_in A data.table with significant DEG results. Must contain
#' columns: batch, group, type, lfc_abs_lim, min_signif, and data (containing
#' uniprot IDs). The 'type' column should not contain NA values. Group and
#' batch names must not contain '&' characters.
#' @param min_to_displayed_threshold Numeric; minimal values to display on plots.
#' If less than 1, interpreted as proportion of universe. Default is 0.01.
#' @param to_keep A data.table or data.frame with columns 'batch' and 'group'
#' to filter the analysis. If NULL, all data is kept. Default is NULL.
#' @param ext_groups External groups to include in the analysis. Can be:
#' a list/vector of UniProt IDs, or a data.table/data.frame with 'batch' and
#' 'group' columns. Default is NULL.
#' @param mode_ext_groups Character; mode for external groups when ext_groups
#' is a table. Must be "individual", "union", or "intersect". Default is "union".
#' @param inter_batch A named list specifying batches to compare for inter-batch
#' analysis. If NULL and mode_inter_batch is not NULL, includes all batches.
#' Default is NULL.
#' @param mode_inter_batch A named list specifying the mode ("union" or "intersect")
#' for each inter-batch comparison. If NULL, no inter-batch comparison is performed.
#' Default is NULL.
#' @param shape Character; shape for Euler diagrams. Default is "ellipse".
#' @param sort_group_for_region_naming Logical; whether to sort group names
#' alphabetically in region names. Default is TRUE.
#' @param path_dir_to_save_plot Character; directory path to save individual
#' plots as PDF files. If NULL, plots are not saved. Default is NULL.
#' @param dict_uniprot_to_symbol Named character vector mapping UniProt IDs
#' to gene symbols. If NULL and 'symbol' column exists in input data, it will
#' be reconstructed from the input table. Default is NULL.
#'
#' @return A list with two elements:
#' \describe{
#'   \item{euler_plots}{A list of Euler plot objects}
#'   \item{dt_region_content}{A data.table with specific content of each region
#'   in the same format as the input deg_table_in}
#' }
#'
#' @details
#' The comparison is performed using UniProt IDs. The function handles both
#' intra-batch comparisons (comparing groups within the same batch) and
#' inter-batch comparisons (comparing batches across different conditions).
#'
#' @examples
#' \dontrun{
#' # Basic usage
#' result <- compute_euler_plot(deg_table)
#'
#' # With filtering and external groups
#' result <- compute_euler_plot(
#'   deg_table_in = my_deg_table,
#'   to_keep = data.table(batch = "B1", group = c("G1", "G2")),
#'   ext_groups = c("P12345", "P67890"),
#'   path_dir_to_save_plot = "plots/"
#' )
#' }
#'
#' @export
compute_euler_plot <- function(
  deg_table_in,
  min_to_displayed_threshold = 0.01,
  to_keep = NULL,
  ext_groups = NULL,
  mode_ext_groups = "union",
  inter_batch = NULL,
  mode_inter_batch = NULL,
  shape = "ellipse",
  sort_group_for_region_naming = TRUE,
  path_dir_to_save_plot = NULL,
  dict_uniprot_to_symbol = NULL
) {

  # Input validation
  if (!data.table::is.data.table(deg_table_in)) {
    stop("deg_table_in must be a data.table", call. = FALSE)
  }

  required_cols <- c("batch", "group", "type", "lfc_abs_lim", "min_signif", "data")
  missing_cols <- setdiff(required_cols, names(deg_table_in))
  if (length(missing_cols) > 0L) {
    stop(
      "Missing required columns: ",
      paste(missing_cols, collapse = ", "),
      call. = FALSE
    )
  }

  # Prevent modification of input table
  deg_table <- data.table::copy(deg_table_in)

  # Remove rows with NA type
  deg_table <- deg_table[!is.na(type)]

  # Build symbol dictionary if needed
  with_symbol <- all(sapply(deg_table$data, function(dt) "symbol" %in% names(dt)))
  if (with_symbol && is.null(dict_uniprot_to_symbol)) {
    dict_uniprot_to_symbol <- unique(data.table::rbindlist(deg_table$data))
    data.table::setorder(dict_uniprot_to_symbol)  # Ensure reproducibility
    dict_uniprot_to_symbol <- stats::setNames(
      dict_uniprot_to_symbol$symbol,
      dict_uniprot_to_symbol$uniprot
    )
  }

  # Check for duplicate UniProt-symbol mappings
  if (!is.null(dict_uniprot_to_symbol)) {
    dup_uniprots <- dict_uniprot_to_symbol[duplicated(names(dict_uniprot_to_symbol))]
    if (length(dup_uniprots) > 0L) {
      warning(
        "Some UniProt IDs are associated with different symbols. ",
        "Using first association only.", call. = FALSE
      )
    }
  }

  # Extract UniProt IDs
  deg_table[, uniprot := lapply(data, function(x) x$uniprot)]
  deg_table[, data := NULL]

  # Handle external groups
  deg_table_unfiltered <- NULL
  ext_sets <- NULL

  if (!is.null(ext_groups)) {
    if (data.table::is.data.table(ext_groups) || is.data.frame(ext_groups)) {
      deg_table_unfiltered <- data.table::copy(deg_table)
    } else if (is.list(ext_groups) || is.vector(ext_groups)) {
      ext_sets <- list(unique(unlist(ext_groups)))
      names(ext_sets) <- "EXT_LIST"
    } else {
      stop("ext_groups must be a data.table, data.frame, list, or vector", call. = FALSE)
    }
  }

  # Filter batches/groups
  if (!is.null(to_keep)) {
    deg_table <- deg_table[to_keep, on = .(batch, group)]
    missing_rows <- deg_table[is.na(type)]
    if (nrow(missing_rows) > 0L) {
      warning("Some 'to_keep' values were not found in 'deg_table'", call. = FALSE)
      print(missing_rows[, .(batch, group)])
      deg_table <- deg_table[!is.na(type)]
    }
  }

  # Create output directory if needed
  if (!is.null(path_dir_to_save_plot) && !dir.exists(path_dir_to_save_plot)) {
    dir.create(path_dir_to_save_plot, recursive = TRUE)
  }

  euler_plots <- list()
  list_new_rows <- list()

  # Intra-batch analysis
  list_deg_table <-
    split(deg_table, by = c("batch", "lfc_abs_lim", "min_signif", "type"))

  for (deg_dt_name in names(list_deg_table)) {
    deg_dt_i <- list_deg_table[[deg_dt_name]]

    lfc_abs_lim_i <- deg_dt_i$lfc_abs_lim[1L]
    min_signif_i <- deg_dt_i$min_signif[1L]
    type_i <- deg_dt_i$type[1L]

    # Create sets for this analysis
    sets <- lapply(deg_dt_i$uniprot, function(uniprot_ids) {
      as.character(na.omit(unlist(uniprot_ids)))
    })
    names(sets) <- deg_dt_i$group

    # Add external sets if needed
    if (!is.null(deg_table_unfiltered)) {
      deg_table_unfiltered_i <- deg_table_unfiltered[
        lfc_abs_lim == lfc_abs_lim_i
        & min_signif == min_signif_i
        & type == type_i
      ]
      ext_sets <- get_ext_group(
        deg_dt = deg_table_unfiltered_i,
        ext_groups_dt = ext_groups,
        mode_ext_groups = mode_ext_groups,
        deg_dt_name = deg_dt_name
      )
    }
    sets <- c(sets, ext_sets)

    # Create Euler plot
    path_file <- if (!is.null(path_dir_to_save_plot)) {
      paste0(path_dir_to_save_plot, "/", deg_dt_name, ".pdf")
    } else {
      NULL
    }

    euler_diagram <- draw_single_euler_plot(
      sets = sets,
      min_to_displayed_threshold = min_to_displayed_threshold,
      title = paste0(deg_dt_name, " | "),
      shape = shape,
      save_path_file = path_file
    )

    euler_plots[[deg_dt_name]] <- euler_diagram

    # Extract region content
    regions_content <- extract_venn_zones_content(
      sets,
      sort_group_for_region_naming = sort_group_for_region_naming
    )

    # Format results
    batch_i <- deg_dt_i$batch[1L]
    for (region_name in names(regions_content$specific_intersect)) {
      region_content <- regions_content$specific_intersect[[region_name]]

      if (with_symbol && !is.null(dict_uniprot_to_symbol)) {
        new_row <- list(
          batch = paste0("_dcmp.", batch_i),
          group = paste0("_sp.", region_name),
          lfc_abs_lim = lfc_abs_lim_i,
          min_signif = min_signif_i,
          type = type_i,
          data = list(data.table::data.table(
            symbol = dict_uniprot_to_symbol[region_content],
            uniprot = region_content
          ))
        )
      } else {
        new_row <- list(
          batch = paste0("_dcmp.", batch_i),
          group = paste0("_sp.", region_name),
          lfc_abs_lim = lfc_abs_lim_i,
          min_signif = min_signif_i,
          type = type_i,
          data = list(data.table::data.table(uniprot = region_content))
        )
      }
      list_new_rows <- append(
        list_new_rows,
        list(new_row)
      )
    }
  }

  # Inter-batch analysis
  if (!is.null(mode_inter_batch)) {
    inter_batch_is_initially_null <- is.null(inter_batch)
    list_deg_table <- split(deg_table, by = c("lfc_abs_lim", "min_signif", "type"))

    for (deg_dt_name in names(list_deg_table)) {
      deg_dt_i <- list_deg_table[[deg_dt_name]]

      lfc_abs_lim_i <- deg_dt_i$lfc_abs_lim[1L]
      min_signif_i <- deg_dt_i$min_signif[1L]
      type_i <- deg_dt_i$type[1L]

      # Add external sets if needed
      if (!is.null(deg_table_unfiltered)) {
        deg_table_unfiltered_i <- deg_table_unfiltered[
          lfc_abs_lim == lfc_abs_lim_i
          & min_signif == min_signif_i
          & type == type_i
        ]
        ext_sets <- get_ext_group(
          deg_dt = deg_table_unfiltered_i,
          ext_groups_dt = ext_groups,
          mode_ext_groups = mode_ext_groups,
          deg_dt_name = deg_dt_name
        )
      }

      # Set up inter-batch comparisons
      if (inter_batch_is_initially_null) {
        inter_batch <- list(rep(
          unique(deg_dt_i$batch),
          length(mode_inter_batch)
        ))
        names(inter_batch) <- names(mode_inter_batch)
      }

      # Process each inter-batch comparison
      for (inter_b_comp_names in names(mode_inter_batch)) {
        deg_dt_i_batch_selected <- deg_dt_i[batch %in% inter_batch[[inter_b_comp_names]]]

        # Flatten UniProt lists by batch
        if (mode_inter_batch[[inter_b_comp_names]] == "union") {
          deg_dt_i_batch_select_flat <- data.table::rbindlist(lapply(
            unique(deg_dt_i_batch_selected$batch), function(x) {
              dt <- deg_dt_i_batch_selected[batch == x]
              data.table::data.table(
                batch = x,
                uniprot = list(as.character(na.omit(unique(unlist(dt$uniprot)))))
              )
            }
          ))
        } else if (mode_inter_batch[[inter_b_comp_names]] == "intersect") {
          deg_dt_i_batch_select_flat <- data.table::rbindlist(lapply(
            unique(deg_dt_i_batch_selected$batch), function(x) {
              dt <- deg_dt_i_batch_selected[batch == x]
              data.table::data.table(
                batch = x,
                uniprot = list(as.character(na.omit(Reduce(intersect, dt$uniprot))))
              )
            }
          ))
        } else {
          stop(
            "Unknown value for mode_inter_batch: '",
            mode_inter_batch[[inter_b_comp_names]], "'", call. = FALSE
          )
        }

        sets <- deg_dt_i_batch_select_flat$uniprot
        names(sets) <- unique(deg_dt_i_batch_selected$batch)
        sets <- c(sets, ext_sets)

        # Create Euler plot
        path_file <- if (!is.null(path_dir_to_save_plot)) {
          paste0(path_dir_to_save_plot, "/Inter-batch_", inter_b_comp_names, ".pdf")
        } else {
          NULL
        }

        euler_diagram <- draw_single_euler_plot(
          sets = sets,
          min_to_displayed_threshold = min_to_displayed_threshold,
          title = paste0("Inter-batch | ", inter_b_comp_names, " | ", deg_dt_name),
          save_path_file = path_file
        )

        euler_plots[[
          paste0("Inter-batch_", inter_b_comp_names, "_", deg_dt_name)
        ]] <- euler_diagram

        # Extract region content
        regions_content <- extract_venn_zones_content(
          sets,
          sort_group_for_region_naming = sort_group_for_region_naming
        )

        # Format results
        batch_i <- inter_b_comp_names
        for (region_name in names(regions_content$specific_intersect)) {
          region_content <- regions_content$specific_intersect[[region_name]]

          if (with_symbol && !is.null(dict_uniprot_to_symbol)) {
            new_row <- list(
              batch = paste0("_dcmp.", batch_i),
              group = paste0("_sp.", region_name),
              lfc_abs_lim = lfc_abs_lim_i,
              min_signif = min_signif_i,
              type = type_i,
              data = list(data.table::data.table(
                symbol = dict_uniprot_to_symbol[region_content],
                uniprot = region_content
              ))
            )
          } else {
            new_row <- list(
              batch = paste0("_dcmp.", batch_i),
              group = paste0("_sp.", region_name),
              lfc_abs_lim = lfc_abs_lim_i,
              min_signif = min_signif_i,
              type = type_i,
              data = list(data.table::data.table(uniprot = region_content))
            )
          }
          list_new_rows <- append(
            list_new_rows,
            list(new_row)
          )
        }
      }
    }
  }

  dt_region_content <- data.table::rbindlist(
    list_new_rows,
    use.names = TRUE,
    fill = FALSE
  )

  return(list(
    euler_plots = euler_plots,
    dt_region_content = dt_region_content
  ))
}
