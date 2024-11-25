######## extend_network_trrust ########

# TRRUST database path
file_trrust <- system.file("extdata", 
                           "human_network_annotations/trrust/data.tab", 
                           package = "fctBio")

# Tests definition
test_that("extend_network_trrust function tests", {
  
  # Base case
  genes <- c("BAX", "MYC", "TP53", "XXX_NEXOMIS_XXX")
  extended_genes <- extend_network_trrust(trrust_db = file_trrust,
                                         input_genes = genes,
                                         add_tf = TRUE,
                                         add_target = TRUE,
                                         keep_input_without_hit = TRUE,
                                         ignore_single_pmid = FALSE)
  results <- c("BAX", "MYC", "TP53", "XXX_NEXOMIS_XXX", "AATF", "ABL1", "AHR", "APC", "AR", "ATF3", "ATM", "BCL6", "BDP1", "BIN1", "BRCA1", "CEBPA", "CEBPE", "CHD8", "CNBP", "CREBBP", "CTCF", "CTNNB1", "DLX5", "DMAP1", "DNMT1", "E2F1", "E2F4", "E2F5", "EGR1", "ELF4", "ELL", "ENO1", "EP300", "ERCC2", "ESR1", "ETS1", "ETS2", "EZH2", "FOS", "FOXA1", "FOXM1", "FUBP1", "GFI1", "HDAC1", "HDAC11", "HDAC2", "HDAC3", "HDAC9", "HIC1", "HIF1A", "HIPK2", "HMGA1", "HNF4A", "HOXA1", "HOXA10", "IFI16", "IKBKB", "ING1", "ING4", "IRF1", "JUN", "KLF4", "L3MBTL1", "LEF1", "MAX", "MDM4", "MLLT10", "MLLT3", "MTA1", "MXI1", "MYB", "MYBL2", "MYCN", "NF1", "NFKB1", "NR3C1", "NUPR1", "PAX2", "PAX5", "PAX8", "PER1", "PGR", "PML", "PPARG", "PRDM1", "PTTG1", "RB1", "RBL1", "RELA", "RFX1", "RUNX3", "RUVBL1", "SIRT1", "SMAD3", "SMAD4", "SMAD7", "SNIP1", "SOX6", "SP1", "SRSF1", "STAT1", "STAT3", "STAT4", "TBL1X", "TBP", "TCF3", "TCF4", "TFDP1", "TLX1", "TP63", "TP73", "VHL", "WDR5", "WT1", "WWTR1", "YEATS4", "YY1", "ZBTB2", "ZBTB7A", "ZNF300", "ZNF382", "ASS1", "AXIN2", "BAG2", "BCL2", "BMI1", "BRD7", "CAD", "CBFB", "CCNA2", "CCNB1", "CCND1", "CCND2", "CCNE1", "PROM1", "CD33", "CD38", "CDC25C", "CDC34", "CDCA7", "CDH3", "CDK2", "CDK4", "CDK6", "CDKN1A", "CFLAR", "CHEK1", "CHEK2", "CXCR4", "CYP3A4", "DDX18", "DKK1", "EBNA1BP2", "FASLG", "FMR1", "FTH1P18", "FUT3", "GATA1", "GATA4", "HLA-A", "HLA-B", "HLA-C", "HLA-E", "HLA-F", "HLA-G", "HNRNPA1", "HNRNPA2B1", "HSPA4", "ID2", "IER3", "IFNA1", "IGF2BP1", "IGK", "IL6", "IPO7", "JUNB", "KRAS", "LDHA", "MAPK1", "MME", "MST1", "MTDH", "NDRG1", "NFE2", "NOL7", "ODC1", "PCNA", "PGK1", "PMAIP1", "PRDX3", "PRODH", "PTBP1", "RCC1", "SART3", "SFRP1", "SHMT1", "ST3GAL1", "ST3GAL3", "ST3GAL4", "SURF1", "SURF2", "TERT", "TFAM", "TFAP4", "TFRC", "TNFRSF10B", "TOP1MT", "UBE2C", "VEGFA", "WRN", "XPO1", "ABCB1", "AFP", "AKT1", "ALKBH2", "ALOX5", "ANKRD1", "ANTXR1", "APAF1", "APCS", "AQP3", "ADGRB1", "BBC3", "BCL2L1", "BCL3", "BDKRB2", "BGLAP", "BIRC5", "BRCA2", "BTG2", "BTG3", "BTS1", "CABLES1", "CARM1", "CASP1", "CASP3", "CAV1", "CCNA1", "CCNB2", "CD82", "CDK1", "CDKN1B", "CHUK", "CKM", "CRYAB", "CTSD", "CYFIP2", "DDB1", "DDB2", "DUSP1", "E2F3", "E2F7", "EGFR", "EPHA2", "FAS", "FDXR", "FGF2", "FLT1", "FOXO3", "FRMD5", "FUS", "GADD45A", "GC", "GDF15", "GLI2", "GML", "GPNMB", "GSTP1", "HRAS", "HSF1", "HSP90AB1", "ID1", "IGF1R", "IGFBP3", "KLF2", "MAD1L1", "MAP4", "MCL1", "MCM7", "MDM2", "MET", "MGMT", "MKI67", "MMP1", "MMP2", "MUC2", "NKX3-1", "NLRC4", "NME1", "OGG1", "PDGFRB", "PERP", "PIAS3", "PLAGL1", "PLK1", "PLK3", "POLD1", "PPARGC1A", "PTPA", "PRC1", "PRG3", "PRKAB1", "PTEN", "PTPN13", "RAD51B", "RASSF1", "RECQL4", "REEP5", "RRM2B", "S100B", "SFN", "SIAH1", "SIVA1", "SLC2A1", "SLC6A6", "SSTR2", "TBXAS1", "TCEAL1", "TCF7L2", "THBS1", "TNFRSF10A", "TP53I3", "TRIM22", "TYMS", "UHRF1", "UHRF1BP1", "VCAN", "XPC")
  expect_identical(extended_genes, results)
  
  # Add only TFs and keeping all input genes
  extended_genes <- extend_network_trrust(trrust_db = file_trrust,
                                                  input_genes = genes,
                                                  add_tf = TRUE,
                                                  add_target = FALSE,
                                                  keep_input_without_hit = TRUE,
                                                  ignore_single_pmid = FALSE)
  results <- c("BAX", "MYC", "TP53", "XXX_NEXOMIS_XXX", "AATF", "ABL1", "AHR", "APC", "AR", "ATF3", "ATM", "BCL6", "BDP1", "BIN1", "BRCA1", "CEBPA", "CEBPE", "CHD8", "CNBP", "CREBBP", "CTCF", "CTNNB1", "DLX5", "DMAP1", "DNMT1", "E2F1", "E2F4", "E2F5", "EGR1", "ELF4", "ELL", "ENO1", "EP300", "ERCC2", "ESR1", "ETS1", "ETS2", "EZH2", "FOS", "FOXA1", "FOXM1", "FUBP1", "GFI1", "HDAC1", "HDAC11", "HDAC2", "HDAC3", "HDAC9", "HIC1", "HIF1A", "HIPK2", "HMGA1", "HNF4A", "HOXA1", "HOXA10", "IFI16", "IKBKB", "ING1", "ING4", "IRF1", "JUN", "KLF4", "L3MBTL1", "LEF1", "MAX", "MDM4", "MLLT10", "MLLT3", "MTA1", "MXI1", "MYB", "MYBL2", "MYCN", "NF1", "NFKB1", "NR3C1", "NUPR1", "PAX2", "PAX5", "PAX8", "PER1", "PGR", "PML", "PPARG", "PRDM1", "PTTG1", "RB1", "RBL1", "RELA", "RFX1", "RUNX3", "RUVBL1", "SIRT1", "SMAD3", "SMAD4", "SMAD7", "SNIP1", "SOX6", "SP1", "SRSF1", "STAT1", "STAT3", "STAT4", "TBL1X", "TBP", "TCF3", "TCF4", "TFDP1", "TLX1", "TP63", "TP73", "VHL", "WDR5", "WT1", "WWTR1", "YEATS4", "YY1", "ZBTB2", "ZBTB7A", "ZNF300", "ZNF382")
  expect_identical(extended_genes, results)
  
  # Add only Targets and keeping all input genes
  extended_genes <- extend_network_trrust(trrust_db = file_trrust,
                                                      input_genes = genes,
                                                      add_tf = FALSE,
                                                      add_target = TRUE,
                                                      keep_input_without_hit = TRUE,
                                                      ignore_single_pmid = FALSE)
  results <- c("BAX", "MYC", "TP53", "XXX_NEXOMIS_XXX", "ASS1", "ATF3", "AXIN2", "BAG2", "BCL2", "BMI1", "BRD7", "CAD", "CBFB", "CCNA2", "CCNB1", "CCND1", "CCND2", "CCNE1", "PROM1", "CD33", "CD38", "CDC25C", "CDC34", "CDCA7", "CDH3", "CDK2", "CDK4", "CDK6", "CDKN1A", "CEBPA", "CFLAR", "CHEK1", "CHEK2", "CREBBP", "CXCR4", "CYP3A4", "DDX18", "DKK1", "EBNA1BP2", "EP300", "FASLG", "FMR1", "FOXM1", "FTH1P18", "FUT3", "GATA1", "GATA4", "HDAC2", "HLA-A", "HLA-B", "HLA-C", "HLA-E", "HLA-F", "HLA-G", "HNRNPA1", "HNRNPA2B1", "HSPA4", "ID2", "IER3", "IFNA1", "IGF2BP1", "IGK", "IL6", "IPO7", "JUNB", "KRAS", "LDHA", "MAPK1", "MME", "MST1", "MTDH", "NDRG1", "NFE2", "NOL7", "ODC1", "PCNA", "PGK1", "PMAIP1", "PRDX3", "PRODH", "PTBP1", "RCC1", "SART3", "SFRP1", "SHMT1", "SRSF1", "ST3GAL1", "ST3GAL3", "ST3GAL4", "SURF1", "SURF2", "TERT", "TFAM", "TFAP4", "TFRC", "TLX1", "TNFRSF10B", "TOP1MT", "TP73", "UBE2C", "VEGFA", "WRN", "XPO1", "ABCB1", "AFP", "AKT1", "ALKBH2", "ALOX5", "ANKRD1", "ANTXR1", "APAF1", "APCS", "AQP3", "ADGRB1", "BBC3", "BCL2L1", "BCL3", "BDKRB2", "BGLAP", "BIRC5", "BRCA1", "BRCA2", "BTG2", "BTG3", "BTS1", "CABLES1", "CARM1", "CASP1", "CASP3", "CAV1", "CCNA1", "CCNB2", "CD82", "CDK1", "CDKN1B", "CHUK", "CKM", "CRYAB", "CTNNB1", "CTSD", "CYFIP2", "DDB1", "DDB2", "DNMT1", "DUSP1", "E2F1", "E2F3", "E2F7", "EGFR", "EGR1", "ELF4", "EPHA2", "ESR1", "EZH2", "FAS", "FDXR", "FGF2", "FLT1", "FOXO3", "FRMD5", "FUS", "GADD45A", "GC", "GDF15", "GLI2", "GML", "GPNMB", "GSTP1", "HIPK2", "HNF4A", "HRAS", "HSF1", "HSP90AB1", "ID1", "IFI16", "IGF1R", "IGFBP3", "ING1", "KLF2", "MAD1L1", "MAP4", "MCL1", "MCM7", "MDM2", "MDM4", "MET", "MGMT", "MKI67", "MMP1", "MMP2", "MUC2", "NFKB1", "NKX3-1", "NLRC4", "NME1", "OGG1", "PAX5", "PDGFRB", "PERP", "PIAS3", "PLAGL1", "PLK1", "PLK3", "PML", "POLD1", "PPARGC1A", "PTPA", "PRC1", "PRG3", "PRKAB1", "PTEN", "PTPN13", "PTTG1", "RAD51B", "RASSF1", "RB1", "RECQL4", "REEP5", "RELA", "RRM2B", "S100B", "SFN", "SIAH1", "SIVA1", "SLC2A1", "SLC6A6", "SMAD3", "SSTR2", "STAT1", "STAT3", "TBXAS1", "TCEAL1", "TCF7L2", "THBS1", "TNFRSF10A", "TP53I3", "TRIM22", "TYMS", "UHRF1", "UHRF1BP1", "VCAN", "XPC")
  expect_identical(extended_genes, results)
  
  # Add only TFs (without force to keep all input genes)
  extended_genes <- extend_network_trrust(trrust_db = file_trrust,
                                                 input_genes = genes,
                                                 add_tf = TRUE,
                                                 add_target = FALSE,
                                                 keep_input_without_hit = FALSE,
                                                 ignore_single_pmid = FALSE)
  results <- c("AATF", "ABL1", "AHR", "APC", "AR", "ATF3", "ATM", "BCL6", "BDP1", "BIN1", "BRCA1", "CEBPA", "CEBPE", "CHD8", "CNBP", "CREBBP", "CTCF", "CTNNB1", "DLX5", "DMAP1", "DNMT1", "E2F1", "E2F4", "E2F5", "EGR1", "ELF4", "ELL", "ENO1", "EP300", "ERCC2", "ESR1", "ETS1", "ETS2", "EZH2", "FOS", "FOXA1", "FOXM1", "FUBP1", "GFI1", "HDAC1", "HDAC11", "HDAC2", "HDAC3", "HDAC9", "HIC1", "HIF1A", "HIPK2", "HMGA1", "HNF4A", "HOXA1", "HOXA10", "IFI16", "IKBKB", "ING1", "ING4", "IRF1", "JUN", "KLF4", "L3MBTL1", "LEF1", "MAX", "MDM4", "MLLT10", "MLLT3", "MTA1", "MXI1", "MYB", "MYBL2", "MYC", "MYCN", "NF1", "NFKB1", "NR3C1", "NUPR1", "PAX2", "PAX5", "PAX8", "PER1", "PGR", "PML", "PPARG", "PRDM1", "PTTG1", "RB1", "RBL1", "RELA", "RFX1", "RUNX3", "RUVBL1", "SIRT1", "SMAD3", "SMAD4", "SMAD7", "SNIP1", "SOX6", "SP1", "SRSF1", "STAT1", "STAT3", "STAT4", "TBL1X", "TBP", "TCF3", "TCF4", "TFDP1", "TLX1", "TP53", "TP63", "TP73", "VHL", "WDR5", "WT1", "WWTR1", "YEATS4", "YY1", "ZBTB2", "ZBTB7A", "ZNF300", "ZNF382")
  expect_identical(extended_genes, results)

  # Base case but ingoring relation supported by one single PMID ref
  extended_genes <- extend_network_trrust(trrust_db = file_trrust,
                                                            input_genes = genes,
                                                            add_tf = TRUE,
                                                            add_target = TRUE,
                                                            keep_input_without_hit = TRUE,
                                                            ignore_single_pmid = TRUE)
  results <- c("BAX", "MYC", "TP53", "XXX_NEXOMIS_XXX", "AHR", "BRCA1", "CTCF", "CTNNB1", "E2F1", "ENO1", "ESR1", "FUBP1", "HDAC3", "ING1", "MDM4", "MYB", "MYCN", "NF1", "NFKB1", "PPARG", "PRDM1", "RB1", "RBL1", "RELA", "RFX1", "SIRT1", "SP1", "STAT3", "TBP", "TCF4", "YY1", "CAD", "CCND1", "CDK4", "CDKN1A", "FASLG", "ODC1", "TERT", "APAF1", "BBC3", "BCL2", "BRCA2", "DDB2", "DNMT1", "DUSP1", "EZH2", "FAS", "GADD45A", "IGF1R", "MDM2", "MGMT", "MMP2", "RRM2B", "THBS1", "TNFRSF10B", "TP53I3", "TP73", "UHRF1", "XPC")
  expect_identical(extended_genes, results)

})

######## venn_genes ########

# Function to test 'compute_euler_plot' following several context using fodler containing input obect (deg_table, ext_groups, to_keep and uniprot_to_symbol), output object and arguments used to generate it.
test_compute_euler_plot_from_dir <- function(test_dir) {
  # load input object and arguments
  deg_table_test <- readRDS(file.path(test_dir, "deg_table_input.rds"))
  ext_groups <- readRDS(file.path(test_dir, "ext_groups_input.rds"))
  execution_args <- readRDS(file.path(test_dir, "execution_args.rds"))
  
  # load reference results
  ref_dt_region_content <- readRDS(file.path(test_dir, "result_dt_output.rds"))
  ref_euler_plot <- readRDS(file.path(test_dir, "result_plot_output.rds"))
  
  # execution of 'compute_euler_plot()'
  res <- do.call(compute_euler_plot, execution_args))
  res_ <- res$dt_region_content

  ##### Test of 'dt_region_content' resuls #####
  test_that(paste("Region content structure ('dt_region_content') for ", test_dir), {
    # check structure
    expect_identical(names(res$dt_region_content), names(ref_dt_region_content),
                info = paste("Column names are different on '$dt_region_content' of ", test_dir))
    expect_identical(sapply(res$dt_region_content, class), sapply(ref_dt_region_content, class),
                info = paste("Column type are different on '$dt_region_content' of ", test_dir))
    
    # check meta data of computed analysis
    cols_no_data <- setdiff(names(res$dt_region_content), "data")
    expect_identical(res$dt_region_content[, ..cols_no_data], ref_dt_region_content[, ..cols_no_data],
                info = paste("List of comparisons (based on 'dt_region_content' rows) performed seems different in ", test_dir))

    # check results of computed analyse
    expect__identical(res$dt_region_content$data, ref_dt_region_content$data,
                info = paste("Results (column data) of content on each reagion are different in ", test_dir))
  })

  ##### Test of 'euler_plot' results #####
  test_that(paste("Plot structure ('euler_plot') for ", test_dir), {
    # check structure
    expect_identical(names(res$euler_plot), names(ref_euler_plot),
                info = paste("Plots names are different on '$euler_plot' of ", test_dir))
    expect_identical(sapply(res$euler_plot, class), sapply(ref_euler_plot, class),
                info = paste("Plots type are different on '$euler_plot' of ", test_dir))
    
    # check results of computed analyse
    expect_identical(res$euler_plot$original.values, ref_euler_plot$original.values,
                info = paste("Value in 'original.values' are different for ", test_dir))
    expect_identical(res$euler_plot$fitted.values, ref_euler_plot$fitted.values,
                info = paste("Value in 'fitted.values' are different for ", test_dir))
    expect_identical(res$euler_plot$quantities$labels, ref_euler_plot$quantities$labels,
                info = paste("Value in 'quantities$labels' are different for ", test_dir))
  })
  
  ##### Complete object comparison #####
  test_that(paste("Complete (substance and format) results object comparison for ", test_dir), {
    # dt_region_content
    if (!identical(res$dt_region_content, ref_dt_region_content)) {
      warning(paste("Object results from 'dt_region_content' show differences in ", test_dir, 
                    ". This could concern only the form of the data, to be verified."))
    }
    
    # euler_plot
    if (!identical(res$euler_plot, ref_euler_plot)) {
      warning(paste("LObject results from 'euler_plot' show differences in ", test_dir, 
                    ". This could concern only the form of the data, to be verified."))
    }
  })
  ##### Test of results files (.pdf and .xlsx ?!? #####
}

# Run tests for all subfolders in the 'test_compute_euler_resources' directory
main_test_dir = "tests/testthat/data_for_venn_test/"
test_dirs <- list.dirs(main_test_dir, full.names = TRUE, recursive = FALSE)
for (test_dir in test_dirs) {
  cat("\n", "Testing folder:", test_dir, "\n")
  test_compute_euler_plot_from_dir(test_dir)
}
