#' @import data.table
#' @import ggplot2
#' @importFrom eulerr euler
#' @importFrom openxlsx createWorkbook addWorksheet writeData saveWorkbook


## Define function

#' Draw simple and single euler plot about a named list of sets
#' @param min_to_displayed_threshold Minimal values displayed. If <1: in proportion of universe (sum of all sets lengths (with redondancy)).
draw_single_euler_plot <- function (sets,
                                    min_to_displayed_threshold = 0.01,
                                    shape = "ellipse",
                                    title = NULL,
                                    save_path_file = NULL) {
  
  euler_data <- euler(sets, shape = shape)
  
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





#' Extract specific content of each zone of a Venn/Euler diagram
#'
#' @param sets_list A named list of sets or a list of named lists of sets. Each list will be processed independently.
#' @param path_xlsx_file_to_save_table File path to save results in an Excel file (default: NULL). If not NULL, creates one sheet per set.
#' 
#' @return A named list of results for each input set or list of sets. Each result contains specific content of each zone.
#'
#' @examples
#' set123 <- list("111" = c("12345", "23456", "34567"),
#'                "222" = c("23456", "45678"),
#'                "333" = c("12345", "45678", "56789"))
#' setABC <- list("AAA" = c("ABCDE", "BCDEF", "DEFGH"),
#'                "BBB" = c("BCDEF", "DEFGH"),
#'                "CCC" = c("DEFGH", "EFGHI", "FGHIJ"))
#' extract_venn_zones_content(list("123" = set123,
#'                                 "ABC" = setABC),
#'                            "test.xlsx")
#'
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
#'
#' @param deg_table Recquired. dt with significant results: results of `deseq_pred$cross_args_and_generate_lists()`.
#' @param to_keep (default: NULL) dt or df with 2 column: 'batch' and 'group'. If 'NULL', no filter, keep all data.
#' @param min_to_displayed_threshold (default: '0.01') Minimal values displayed. If <1: in proportion of universe (sum of all sets lengths (with redondancy) - including potentials ext_groups).
#' @param mode_ext_groups (default: 'union') Must be 'individual', 'union' or 'intersect'.
#' @param ext_groups (default: NULL) dt or df with 2 column: 'batch' and 'group'. If 'NULL', no external set.
#' @param mode_inter_batch (default: NULL) 'union' or 'intersect'. If 'NULL', no Inter-batch euler diagram. Perform inter-batch euler diagram at level of batches: all internal groups are summarized following choosen method.
#' @param shape (default: 'ellipse') euler diagram shape.
#' @param path_dir_to_save_plot (default: NULL) directory path where save individual euler diagram on separate pdf. If 'NULL', no saving.
#' 
#' @return This function return a list of plots.
#' 
#' @details Comparison is performed using uniprot ids.
#' @examples
#' Example usage
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
                               mode_ext_groups = "union",
                               ext_groups = NULL,
                               mode_inter_batch = "union",
                               shape = "ellipse",
                               path_dir_to_save_plot = NULL,
                               path_xlsx_file_to_save_table = NULL) {
  
  # prevent modification of input tables
  deg_table <- copy(deg_table_in)
  
  # extract uniprot ids
  deg_table[, uniprot := lapply(data, function(x) {x$uniprot})]
  deg_table[, data := NULL]
  
  # create exterrnal sets
  if (!is.null(ext_groups)){
    ext_deg_table <- deg_table[ext_groups, on = .(batch, group)]
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
  } else {
    ext_sets = NULL
  }
  
  # filter batchs/groups
  if (!is.null(to_keep)) {
    deg_table <- deg_table[to_keep, on = .(batch, group)]
  }
  
  if (!dir.exists(path_dir_to_save_plot)) {
    dir.create(path_dir_to_save_plot)
  }
  
  sets_list <- list()
  euler_plots <- list()
  
  ### Inter-batch
  list_batches <- unique(deg_table$batch)
  for (batch_i in list_batches) {
    dt_batch <- deg_table[batch == batch_i]
    
    sets <- lapply(dt_batch$uniprot, function(uniprot_ids) {
      as.character(na.omit(unlist(uniprot_ids)))
    })
    names(sets) <- dt_batch$group 
    
    sets <- c(sets, ext_sets)
    
    path_file = if(!is.null(path_dir_to_save_plot)) {
      paste0(path_dir_to_save_plot, "/", batch_i, ".pdf")
      } else NULL
    euler_diagram <- draw_single_euler_plot(sets,
                                            min_to_displayed_threshold,
                                            title = paste0(batch_i, " | "),
                                            shape = shape,
                                            save_path_file = path_file)
  
    euler_plots[[batch_i]] <- euler_diagram
    sets_list[[batch_i]] <- sets
  }
  
  ### Intra-batch: global
  if (!is.null(mode_inter_batch)) {
    if (mode_inter_batch == "union") {
      batch_deg_table <- rbindlist(lapply(unique(deg_table$batch), function(x) {
        dt <- deg_table[batch == x]
        data.table(
          batch = x,
          uniprot = list(as.character(na.omit(unique(unlist(dt$uniprot)))))
        )
      }))
    } else if (mode_inter_batch == "intersect") {
      batch_deg_table <- rbindlist(lapply(unique(deg_table$batch), function(x) {
        dt <- deg_table[batch == x]
        data.table(
          batch = x,
          uniprot = list(as.character(na.omit(Reduce(intersect, dt$uniprot))))
        )
      }))
    } else {
      stop("Error: unknown value for 'mode_inter_batch': '", mode_inter_batch, "'")
    }
    sets = batch_deg_table$uniprot
    names(sets) <- batch_deg_table$batch
    sets <- c(sets, ext_sets)
    
    path_file = if(!is.null(path_dir_to_save_plot)) {
      paste0(path_dir_to_save_plot, "/Intra-batch.pdf")
      } else NULL
    euler_diagram <- draw_single_euler_plot(sets, min_to_displayed_threshold,
                                     title = paste0("Intra-batch", " | "),
                                     save_path_file = path_file)
    
    euler_plots[["Intra-batch"]] <- euler_diagram
    sets_list[["Intra-batch"]] <- sets
    
  }
  
  ### save specific content of each region on xlsx file
  dt_all_regions_content <- extract_venn_zones_content(sets_list,
                                                       path_xlsx_file_to_save_table)
  
  return(list(euler_plots = euler_plots,
              dt_all_regions_content = dt_all_regions_content))
}

