#' Differential Abundance Analysis Functions
#' This file contains essential functions for differential abundance analysis in microbiome data
#' 
#' @author Orr Tobaly & Nikol Elyashov
#' @date April 2025

# Required packages 
# library(ANCOMBC)
# library(phyloseq)
# library(dplyr)
# library(ggplot2)

#' Run ANCOM-BC2 Analysis
#'
#' Core function to run ANCOM-BC2 differential abundance analysis
#'
#' @param ps_object A phyloseq object containing the microbiome data
#' @param tax_level Taxonomic level for analysis (default: "Genus")
#' @param formula_str Fixed effects formula as a string
#' @param group_var The grouping variable for comparison
#' @param prv_cut Prevalence cutoff for filtering taxa (default: 0.025)
#' @param seed Random seed for reproducibility (default: 123)
#' @return ANCOM-BC2 results object
#' @export
run_ancombc2 <- function(ps_object, 
                         tax_level = "Genus",
                         formula_str,
                         group_var,
                         prv_cut = 0.025,
                         seed = 123) {
  
  # Set seed for reproducibility
  set.seed(seed)
  
  # Run ANCOM-BC2
  ancombc2_result <- ancombc2(
    data = ps_object,
    tax_level = tax_level,
    fix_formula = formula_str,
    rand_formula = NULL,
    p_adj_method = "holm",
    pseudo_sens = TRUE,
    prv_cut = prv_cut,
    lib_cut = 1000,
    s0_perc = 0.1,
    group = group_var,
    struc_zero = TRUE,
    neg_lb = TRUE,
    alpha = 0.05,
    n_cl = 2,
    global = FALSE,
    pairwise = TRUE
  )
  
  return(ancombc2_result)
}

#' Visualize ANCOM-BC2 Results
#'
#' Creates basic visualizations of significant taxa from ANCOM-BC2 results
#'
#' @param ancombc2_result ANCOM-BC2 results object
#' @param title Title for the visualization
#' @return A list containing plots and significant taxa information
#' @export
visualize_ancombc2_results <- function(ancombc2_result, title = "ANCOM-BC2 Results") {
  # Extract primary results
  res_primary <- ancombc2_result$res
  
  # Create empty list for plots and significant taxa
  plot_list <- list()
  sig_taxa_list <- list()
  
  # Get the column names related to log fold changes and differential expression
  lfc_cols <- grep("^lfc_", colnames(res_primary), value = TRUE)
  diff_cols <- grep("^diff_", colnames(res_primary), value = TRUE)
  
  # Process each comparison
  for (i in seq_along(lfc_cols)) {
    # Get corresponding column names
    lfc_col <- lfc_cols[i]
    diff_col <- diff_cols[i]
    
    # Skip if not in the data
    if (!lfc_col %in% colnames(res_primary) || !diff_col %in% colnames(res_primary)) {
      next
    }
    
    # Extract the comparison name from the column
    comparison <- gsub("^lfc_", "", lfc_col)
    
    # Create data frame for plotting
    plot_data <- data.frame(
      taxon = res_primary$taxon,
      lfc = res_primary[[lfc_col]],
      diff = res_primary[[diff_col]]
    ) %>%
      # Filter for significant taxa
      filter(diff == 1) %>%
      # Sort by absolute log fold change
      arrange(desc(abs(lfc)))
    
    # Skip if no significant taxa
    if (nrow(plot_data) == 0) {
      next
    }
    
    # Create bar plot of significant taxa
    # Limit to top 20 for visibility
    top_n_taxa <- min(nrow(plot_data), 20)
    bar_data <- plot_data[1:top_n_taxa, ]
    
    # Create bar plot
    bar_plot <- ggplot(bar_data, aes(x = reorder(taxon, lfc), y = lfc, 
                                     fill = lfc > 0)) +
      geom_bar(stat = "identity") +
      scale_fill_manual(values = c("#619CFF", "#F8766D"), 
                        name = "Direction", 
                        labels = c("Decreased", "Increased")) +
      labs(
        title = paste0(title, " - ", comparison),
        x = "Taxon",
        y = "Log Fold Change"
      ) +
      coord_flip() +
      theme_minimal()
    
    # Add plot to list
    plot_list[[comparison]] <- bar_plot
    
    # Add significant taxa to list
    sig_taxa_list[[comparison]] <- plot_data
  }
  
  return(list(plots = plot_list, sig_taxa = sig_taxa_list))
}

#' Extract Significant Taxa from ANCOM-BC2 Results
#'
#' Helper function to extract names of significant taxa from ANCOM-BC2 results
#'
#' @param ancombc2_result ANCOM-BC2 results object
#' @param comparison_col The column name for the specific comparison of interest
#' @param sensitivity_filter Whether to apply sensitivity score filtering (default: TRUE)
#' @return A character vector of significant taxa names
#' @export
extract_significant_taxa <- function(ancombc2_result, comparison_col, sensitivity_filter = TRUE) {
  # Extract primary results
  res_primary <- ancombc2_result$res
  
  # Get the column names for differential abundance and sensitivity scores
  diff_col <- paste0("diff_", comparison_col)
  ss_col <- paste0("passed_ss_", comparison_col)
  
  # Check if columns exist
  if(!diff_col %in% colnames(res_primary)) {
    stop(paste("Column", diff_col, "not found in results"))
  }
  
  # Filter for significant taxa
  if(sensitivity_filter && ss_col %in% colnames(res_primary)) {
    sig_taxa <- res_primary$taxon[res_primary[[diff_col]] == 1 & res_primary[[ss_col]] == TRUE]
  } else {
    sig_taxa <- res_primary$taxon[res_primary[[diff_col]] == 1]
  }
  
  return(sig_taxa)
}