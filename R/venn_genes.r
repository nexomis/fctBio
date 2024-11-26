
## Define function

# Draw simple and single euler plot about a named list of sets
# 
# @param min_to_displayed_threshold Minimal values displayed. If <1: in proportion of universe (sum of all sets lengths (with redondancy)).
draw_single_euler_plot <- function (sets,
                                    min_to_displayed_threshold = 0.01,
                                    shape = "ellipse",
                                    title = NULL,
                                    save_path_file = NULL) {
  
  euler_data <- eulerr::euler(sets, shape = shape)
  
  if (min_to_displayed_threshold < 1) {
    min_to_displayed <- round(length(unique(unlist(sets))) * min_to_displayed_threshold)
  } else {
    min_to_displayed = min_to_displayed_threshold
  }
  qts <- euler_data$original.values
  qts[qts < min_to_displayed] <- ""
  
  euler_diagram <- plot(euler_data,
                        quantities = qts,
                        lty = 1:6,
                        edges = "black",
                        main = paste0(title, "(min.val: ", min_to_displayed, ")"))
  
  if (!is.null(save_path_file)) {
    ggsave(
      filename = save_path_file,
      plot = euler_diagram,
      width = 11.69,
      height = 8.27)
  }
  
  return(euler_diagram)
}





# Extract specific content of each zone of a Venn/Euler diagram
#
# @param sets_list A named list of sets or a list of named lists of sets. Each list will be processed independently.
# @param path_xlsx_file_to_save_table File path to save results in an Excel file (default: NULL). If not NULL, creates one sheet per set.
# 
# @return A named list of results for each input set or list of sets. Each result contains specific content of each zone.
#
# @examples
# set123 <- list("111" = c("12345", "23456", "34567"),
#                "222" = c("23456", "45678"),
#                "333" = c("12345", "45678", "56789"))
# setABC <- list("AAA" = c("ABCDE", "BCDEF", "DEFGH"),
#                "BBB" = c("BCDEF", "DEFGH"),
#                "CCC" = c("DEFGH", "EFGHI", "FGHIJ"))
# extract_venn_zones_content(list("123" = set123,
#                                 "ABC" = setABC),
#                            "test.xlsx")
#
extract_venn_zones_content <- function(sets_list, path_xlsx_file_to_save_table = NULL) {
  # if only one set as input, transfromr to list of set to simplify function
  if (!is.list(sets_list[[1]])) {
    sets_list <- list("specific_intersect" = sets_list)
  }
  
  all_results <- list()
  # for each global list
  for (set_name in names(sets_list)) {
    sets <- sets_list[[set_name]]
    
    # name all region
    region_names <- unlist(lapply(1:length(names(sets)), function(x) {
      combn(names(sets), x, paste, collapse = "&")
    }))
    
    regions_content <- list()
    
    # specific content of each region
    for (region in region_names) {
      incl_sets <- unlist(strsplit(region, "&"))
      excl_sets <- setdiff(names(sets), incl_sets)
      
      # first, add all elements included in incl_sets. then, remove all elements included in excl_sets
      region_i_contents <- as.character(na.omit(Reduce(intersect, sets[incl_sets])))
      if (length(excl_sets) > 0) {
        region_i_contents <- setdiff(region_i_contents, unlist(sets[excl_sets]))
      }
      
      regions_content[[region]] <- region_i_contents
    }
    
    all_results[[set_name]] <- regions_content
  }
  
  # save result on excel file
  if (!is.null(path_xlsx_file_to_save_table)) {
    wb <- createWorkbook()
    for (set_name in names(all_results)) {
      addWorksheet(wb, set_name)
      # convert result in table format
      max_length <- max(sapply(all_results[[set_name]], length))
      region_dt <- as.data.table(
        do.call(cbind, lapply(all_results[[set_name]], function(x) {
          if (length(x) == 0) return(rep(NA, max_length))
          return(c(x, rep(NA, max_length - length(x))))
        })))
      writeData(wb, set_name, region_dt, rowNames = FALSE)
    }
    saveWorkbook(wb, path_xlsx_file_to_save_table, overwrite = TRUE)
  }
  
  return(all_results)
}





#' Compute and draw euleur plot based on significant DEG table
#'
#' Computed at uniprot level ! 'NA' value has been ignored. (quid of multiple uniprot_id for a single gene_symbol ?!?!?!)
#'
#' @param deg_table_in Recquired. dt with significant results: results of `deseq_pred$cross_args_and_generate_lists()`. (Note: differente, pvalue, log_FC and deregulation type has been managed separatly !). Note: the “type” column should not contain “NA”, otherwise false warnings that some 'to_keep' or 'ext_group' groups/batches will be incorrectly issued.
#' @param to_keep (default: NULL) dt or df with 2 column: 'batch' and 'group'. If 'NULL', no filter, keep all data. Impact inter AND intra batch comparison but no impact 'ext_group'.
#' @param min_to_displayed_threshold (default: '0.01') Minimal values displayed. If <1: in proportion of universe (sum of all sets lengths (with redondancy) - including potentials ext_groups).
#' @param ext_groups (default: NULL) list of uniprot ids used as ext_group (list or vector). If table (dt or df with 2 column: 'batch' and 'group') determine the gene list from 'deg_table_in'. If 'NULL', no external set. Used for inter AND intra batche comparison. For now: ignore other parameters ("type', 'lfc', 'pval').
#' @param mode_ext_groups (default: 'union') Used only if 'ext_groups' is table. Must be 'individual', 'union' or 'intersect'.
#' @param inter_batch (default: NULL) list of named list with batch to compare together. If NULL, and 'mode_inter_batch' not NULL, include all batche. Performed independantly for each 'type', 'pval' and 'lfc'.
#' @param mode_inter_batch (default: NULL) list of named list (same name as 'inter_batch') specifying mode of selection of deg to each batch for each inter-batch comparison ('union' or 'intersect'). If 'NULL', no Inter-batch comparison.
#' @param shape (default: 'ellipse') euler diagram shape.
#' @param path_dir_to_save_plot (default: NULL) directory path where save individual euler diagram on separate pdf. If 'NULL', no saving.
#' @param dict_uniprot_to_symbol (default: NULL) used to reconstitude column 'symbol' of table in 'data' of each result if 'symbol' column is already present in 'deg_table_in' input (if 'symbol column not in input table, this parameters will be ignored). If NULL, reconstructed from input table. 
#' 
#' @return  
#' This function return as first element a list of eullerr plots object ('$euler_plot') and as second element a data.table with specific content of each region on format of 'deg_table_in' ('$dt_region_content').
#' 
#' @details 
#' Comparison is performed using uniprot ids.  
#'
#' TODO:
#'    - in option management of up/down. Although this is already possible in the current state of development by transforming deg_table upstream, it would be a good idea to add the following 2 options)
#'       - summarize analyse: merge up and down on single analyse (after prefixing genes by their status "up_"/"down_")
#'       - very interesting analyse: venn with up and down in same plot on disctinct sets to highlight swtched genes !
#'    - ext_groups:
#'       - add possibility to put mutliple specific gene list as vector (with mode: union, individual, intersect)
#'       - add option to use ext group at inter AND/OR intra batche comparison (actually one sigle group used for both)
#'    - plot: 
#'       - put all venn diagram on same pdf using output/append on new page option !
#'       - manage tittle parameters
#'    - reduce time running: don't perform independantly euleurr plot ! use already computed intersection to put stats to eulleur diagram
#' 
#' @examples
#' Example usage
#' 
#' compute_euler_plot(deg_table = data.table(deseq_pred$cross_args_and_generate_lists(
#'                                             cross_type = c("deregulated"),
#'                                             cross_lfc_abs_lim = c(LFC_threshold),
#'                                             cross_min_signif = c(pvaladj_threshold)
#'                                          )
#'                   )
#'                   
#' @export
compute_euler_plot <- function(deg_table_in,
                               min_to_displayed_threshold = 0.01,
                               to_keep = NULL,
                               ext_groups = NULL,
                               mode_ext_groups = "union",
                               inter_batch = NULL,
                               mode_inter_batch = NULL,
                               shape = "ellipse",
                               path_dir_to_save_plot = NULL,
                               dict_uniprot_to_symbol = NULL) {
  
  # prevent modification of input tables
  deg_table <- copy(deg_table_in)
  
  # build 'dict_uniprot_to_symbol' column 'symbol' is present on all sub_dt of input and if dict_uniprot_to_symbol is null
  with_symbol <- all(sapply(deg_table_base$data, function(dt) "symbol" %in% names(dt)))
  if (with_symbol & is.null(dict_uniprot_to_symbol)) {
    dict_uniprot_to_symbol <- unique(rbindlist(deg_table$data))
    # sort only to ensure reproducibility in the case of nonunique uniprot/symbol combinations
    setorder(dict_uniprot_to_symbol)
    dict_uniprot_to_symbol <- setNames(dict_uniprot_to_symbol$symbol,
                                       dict_uniprot_to_symbol$uniprot)
  }
  # warning if nonunique uniprot/symbol combinations
  dup_uniprots <- dict_uniprot_to_symbol[duplicated(names(dict_uniprot_to_symbol))]
  if (length(dup_uniprots) > 0) {
    warning("WARN: some uniprot ids are associated with different symbol ids: '",
            dup_uniprots,
            "'. Association of symbol id will be based only on the first.")
  }
  
  # extract uniprot ids
  deg_table[, uniprot := lapply(data, function(x) {x$uniprot})]
  deg_table[, data := NULL]
  
  # create exterrnal sets (TODO include pval fc and type on ext_group design ?)
  if (!is.null(ext_groups)){
    if (is.data.table(ext_groups) || is.data.frame(ext_groups)) {
        ext_deg_table <- deg_table[ext_groups, on = .(batch, group)]
        # check missing values
        missing_rows <- ext_deg_table[is.na(type)]
        if (nrow(missing_rows) > 0) {
        warning("WARN: Some 'ext_groups' values were not found in 'deg_table' (not included in 'ext_groups').")
        print(missing_rows[, .(batch, group)])
        }
        if (mode_ext_groups == "individual") {
        ext_sets <- lapply(ext_deg_table$uniprot, function(uniprot_ids) {
            unique(na.omit(unlist(uniprot_ids)))
        })
        names(ext_sets) <- paste("EXT", ext_deg_table$batch, ext_deg_table$group, sep=":")
        } else if (mode_ext_groups == "union") {
        ext_sets <- list(unique(na.omit(unlist(ext_deg_table$uniprot))))
        names(ext_sets) <- "EXT:union"
        } else if (mode_ext_groups == "intersect") {
        ext_sets <- list(as.character(na.omit(Reduce(intersect, ext_deg_table$uniprot))))
        names(ext_sets) <- "EXT:intersect"
        } else {
        stop("Error: unknown value for 'mode_ext_groups': '", mode_ext_groups, "'")
        }
    } else if (is.list(ext_groups) || is.vector(ext_groups)) {
        ext_sets <- list(unique(unlist(ext_groups)))
        names(ext_sets) <- "EXT_LIST"
    } else {
        stop("'ext_groups' must be a data.table, data.frame, list, or vector.")
    }
  } else {
    ext_sets = NULL
  }
  
  # filter batchs/groups
  if (!is.null(to_keep)) {
    deg_table <- deg_table[to_keep, on = .(batch, group)]
    # check missing values
    missing_rows <- deg_table[is.na(type)]
    if (nrow(missing_rows) > 0) {
      warning("WARN: Some 'to_keep' values were not found in 'deg_table' (not included in 'to_keep').")
      print(missing_rows[, .(batch, group)])
    }
  }
  
  if (!is.null(path_dir_to_save_plot) && !dir.exists(path_dir_to_save_plot)) {
    dir.create(path_dir_to_save_plot)
  }
  
  euler_plots <- list()
  list_new_rows_dt_region_content <- list()

  ### Intra-batch
  list_deg_table <- split(deg_table, by = c("batch", "lfc_abs_lim", "min_signif", "type"))
  
  for (deg_dt_name in names(list_deg_table)) {
    deg_dt_i <- list_deg_table[[deg_dt_name]]
    sets <- lapply(deg_dt_i$uniprot, function(uniprot_ids) {
      as.character(na.omit(unlist(uniprot_ids)))
    })
    names(sets) <- deg_dt_i$group 
    
    sets <- c(sets, ext_sets)
    
    path_file = if(!is.null(path_dir_to_save_plot)) {
      paste0(path_dir_to_save_plot, "/", deg_dt_name, ".pdf")
      } else NULL
    euler_diagram <- draw_single_euler_plot(sets,
                                            min_to_displayed_threshold,
                                            title = paste0(deg_dt_name, " | "),
                                            shape = shape,
                                            save_path_file = path_file)

    euler_plots[[deg_dt_name]] <- euler_diagram
    
    regions_content <- extract_venn_zones_content(sets)
    
    # formatting of set list results from “xxx” to initial “deg_table” format (temporarily saved on a list)
    lfc_abs_lim = deg_dt_i$lfc_abs_lim[1]
    min_signif = deg_dt_i$min_signif[1]
    type = deg_dt_i$type[1]
    batch = deg_dt_i$batch[1]  # better than "deg_dt_name" (e.g: 'Batch1.0.5849625.0.05.deregulated')
    for (region_name in names(regions_content$specific_intersect)) {
      region_content <- regions_content$specific_intersect[[region_name]]
      if(with_symbol) {
        new_row <- list(batch = paste0("_dcmp.", batch),
                        group = paste0("_sp.", region_name),
                        lfc_abs_lim = lfc_abs_lim,
                        min_signif = min_signif,
                        type = type,
                        data = list(data.table(dict_uniprot_to_symbol[region_content],
                                               uniprot = region_content))
        )
      } else {
        new_row <- list(batch = paste0("_dcmp.", batch),
                        group = paste0("_sp.", region_name),
                        lfc_abs_lim = lfc_abs_lim,
                        min_signif = min_signif,
                        type = type,
                        data = list(data.table(uniprot = region_content))
        )
      }
      list_new_rows_dt_region_content <- append(list_new_rows_dt_region_content,
                                                list(new_row))
    }
  }
  
  ### Inter-batch
  if (is.null(inter_batch)) {
    inter_batch_is_initially_null = TRUE
  } else {
    inter_batch_is_initially_null = FALSE
  }
  list_deg_table <- split(deg_table, by = c("lfc_abs_lim", "min_signif", "type"))
  
  if (!is.null(mode_inter_batch)) {
    # for each threshold in deg_table
    for (deg_dt_name in names(list_deg_table)) {
      deg_dt_i <- list_deg_table[[deg_dt_name]]
      # if not list batch to include on comparison, include all available batches
      if (inter_batch_is_initially_null) {
        inter_batch <- list(rep(unique(deg_dt_i$batch),
                                n=length(mode_inter_batch)))
        names(inter_batch) <- names(mode_inter_batch)
      }
      # for each recquested inter-batch
      for (inter_b_comp_names in names(mode_inter_batch)) {
        deg_dt_i_batch_selected <- deg_dt_i[batch %in% inter_batch[[inter_b_comp_names]] ]
        # flatten uniprot list by batch (union or interesect)
        if (mode_inter_batch[[inter_b_comp_names]] == "union") {
          deg_dt_i_batch_selected_flattened <- rbindlist(lapply(
            unique(deg_dt_i_batch_selected$batch), function(x) {
              dt <- deg_dt_i_batch_selected[batch == x]
              data.table(
                batch = x,
                uniprot = list(as.character(na.omit(unique(unlist(dt$uniprot)))))
            )
          }))
        } else if (mode_inter_batch[[inter_b_comp_names]] == "intersect") {
          deg_dt_i_batch_selected_flattened <- rbindlist(lapply(
            unique(deg_dt_i_batch_selected$batch), function(x) {
              dt <- deg_dt_i_batch_selected[batch == x]
              data.table(
                batch = x,
                uniprot = list(as.character(na.omit(Reduce(intersect, dt$uniprot))))
            )
          }))
        } else {
          stop("Error: unknown value for specific element in 'mode_inter_batch': '", mode_inter_batch[[inter_b_comp_names]], "'")
        }
        sets = deg_dt_i_batch_selected_flattened$uniprot
        names(sets) <- unique(deg_dt_i_batch_selected$batch)
        sets <- c(sets, ext_sets)
        
        # euler_plot
        path_file = if(!is.null(path_dir_to_save_plot)) {
          paste0(path_dir_to_save_plot, "/Intra-batch_", inter_b_comp_names,".pdf")
        } else NULL
        euler_diagram <- draw_single_euler_plot(sets, min_to_displayed_threshold,
                                                title = paste0("Intra-batch | ", inter_b_comp_names, " | ", deg_dt_name),
                                                save_path_file = path_file)
        euler_plots[[paste0("Intra-batch_", inter_b_comp_names, "_", deg_dt_name)]] <- euler_diagram
        # content by region
        regions_content <- extract_venn_zones_content(sets)
        # formatting of set list results from “xxx” to initial “deg_table” format
        lfc_abs_lim = deg_dt_i$lfc_abs_lim[1]
        min_signif = deg_dt_i$min_signif[1]
        type = deg_dt_i$type[1]
        batch = inter_b_comp_names  # better than "paste0(deg_dt_name, ":", mode_inter_batch[[inter_b_comp_names]])" (e.g: 'Batch1.0.5849625.0.05.deregulated:union')
        for (region_name in names(regions_content$specific_intersect)) {
          region_content <- regions_content$specific_intersect[[region_name]]
          if(with_symbol) {
            new_row <- list(batch = paste0("_dcmp.", batch),
                            group = paste0("_sp.", region_name),
                            lfc_abs_lim = lfc_abs_lim,
                            min_signif = min_signif,
                            type = type,
                            data = list(data.table(symbol = dict_uniprot_to_symbol[region_content],
                                                   uniprot = region_content))
            )
          } else {
            new_row <- list(batch = paste0("_dcmp.", batch),
                            group = paste0("_sp.", region_name),
                            lfc_abs_lim = lfc_abs_lim,
                            min_signif = min_signif,
                            type = type,
                            data = list(data.table(uniprot = region_content))
            )
          }
          list_new_rows_dt_region_content <- append(list_new_rows_dt_region_content,
                                                    list(new_row))
        }                      
      }      
    }
    
  }

  dt_region_content <- rbindlist(list_new_rows_dt_region_content,
                                 use.names = TRUE,
                                 fill = FALSE)
  
  return(list(euler_plots = euler_plots,
              dt_region_content = dt_region_content))
}



