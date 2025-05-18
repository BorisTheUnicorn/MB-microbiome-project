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

#' Analyze and Summarize ANCOM-BC2 Model Results
#'
#' This function takes a single ANCOM-BC2 model output, extracts results
#' for a specified primary variable, summarizes significant findings, and
#' generates a bar plot of top differentially abundant taxa.
#'
#' @param ancom_model A list object returned by `ancombc2`.
#' @param model_name A character string to name the model in outputs (e.g., "Model 1: Carriage Group + Batch").
#' @param primary_variable A character string specifying the main variable of interest
#'        (e.g., "carriage_group"). The function will look for LFC and diff columns
#'        related to this variable.
#' @param phyloseq_object A phyloseq object that was used to run the ANCOM-BC2 model.
#'        This is used to determine the levels of the primary_variable for constructing
#'        column names.
#' @param reference_level A character string specifying the reference level of the
#'        `primary_variable`. If NULL (default), the function attempts to infer it or
#'        will require careful naming of `primary_variable_comparison_level`.
#' @param primary_variable_comparison_level A character string specifying the specific
#'        comparison level of the `primary_variable` to analyze (e.g., "Acquisition_carrier"
#'        if "Non_carrier" is the reference for "carriage_group"). If NULL, the function
#'        will try to determine this based on the `primary_variable` and `reference_level`
#'        or use the first non-reference level found.
#' @param alpha_significance The significance level (default: 0.05).
#' @param filter_on_sensitivity Boolean, whether to filter significant taxa by `passed_ss_` columns (default: TRUE).
#' @param top_n_taxa_plot Integer, number of top taxa to display in the bar plot (default: 20).
#' @return A list containing:
#'         - `model_name`: The provided name of the model.
#'         - `model_formula`: The formula used in the ANCOM-BC2 model.
#'         - `taxonomic_level`: The taxonomic level analyzed.
#'         - `primary_comparison_details`: Description of the comparison made.
#'         - `summary_stats`: Data frame with summary statistics.
#'         - `significant_taxa_df`: Data frame of significant taxa for the primary variable.
#'         - `plot_significant_taxa`: ggplot object of the bar plot, or NULL if no significant taxa.
#'         - `error_message`: NULL if successful, or an error message string.
#' @export
#' @examples
#' # Assuming 'ps_baseline' is a phyloseq object and 'model_m1' is an ANCOM-BC2 result
#' # analysis_m1 <- analyze_ancom_model_results(
#' #   ancom_model = model_m1,
#' #   model_name = "Model 1: carriage_group + batch",
#' #   primary_variable = "carriage_group",
#' #   phyloseq_object = ps_baseline, # Crucial for determining levels
#' #   reference_level = "Non_carrier", # Explicitly state the reference
#' #   primary_variable_comparison_level = "Acquisition_carrier",
#' #   alpha_significance = 0.05,
#' #   filter_on_sensitivity = TRUE
#' # )
#' # if(!is.null(analysis_m1$plot_significant_taxa)) print(analysis_m1$plot_significant_taxa)
#' # if(nrow(analysis_m1$significant_taxa_df) > 0) DT::datatable(analysis_m1$significant_taxa_df)

analyze_ancom_model_results <- function(ancom_model,
                                        model_name,
                                        primary_variable,
                                        phyloseq_object, # Added to help determine levels
                                        reference_level = NULL, # Added for clarity
                                        primary_variable_comparison_level = NULL, # Added for explicit control
                                        alpha_significance = 0.05,
                                        filter_on_sensitivity = TRUE,
                                        top_n_taxa_plot = 20) {
  
  # Initialize return list structure
  return_list <- list(
    model_name = model_name,
    model_formula = NULL,
    taxonomic_level = NULL,
    primary_comparison_details = NULL,
    summary_stats = data.frame(
      Description = character(),
      Value = character(),
      stringsAsFactors = FALSE
    ),
    significant_taxa_df = data.frame(),
    plot_significant_taxa = NULL,
    error_message = NULL
  )
  
  # --- 1. Basic Checks ---
  if (is.null(ancom_model) || !is.list(ancom_model)) {
    return_list$error_message <- "ancom_model is NULL or not a list."
    return(return_list)
  }
  if (!is.null(ancom_model$error)) {
    return_list$error_message <- paste("ANCOM-BC2 model generation failed:", ancom_model$error)
    return(return_list)
  }
  if (is.null(ancom_model$res) || !is.data.frame(ancom_model$res)) {
    return_list$error_message <- "ANCOM-BC2 result does not contain a valid 'res' data frame."
    return(return_list)
  }
  if (is.null(phyloseq_object) || !inherits(phyloseq_object, "phyloseq")) {
    return_list$error_message <- "A valid phyloseq_object is required to determine variable levels."
    return(return_list)
  }
  
  
  res_df <- ancom_model$res
  return_list$model_formula <- ifelse(!is.null(ancom_model$formula), ancom_model$formula, "Formula not available")
  return_list$taxonomic_level <- ifelse(!is.null(ancom_model$tax_level), ancom_model$tax_level, "Tax level not available")
  
  # --- 2. Identify LFC, diff, and sensitivity score columns for the primary_variable ---
  
  # Determine the comparison level if not explicitly provided
  if (is.null(primary_variable_comparison_level)) {
    if (!primary_variable %in% names(sample_data(phyloseq_object))) {
      return_list$error_message <- paste0("primary_variable '", primary_variable, "' not found in phyloseq_object sample_data.")
      return(return_list)
    }
    var_levels <- levels(as.factor(get_variable(phyloseq_object, primary_variable)))
    if (is.null(reference_level)) {
      # Try to infer reference (often the first level alphabetically or by factor order)
      reference_level <- var_levels[1]
      message(paste0("Inferring '", reference_level, "' as the reference level for '", primary_variable, "'. Explicitly set reference_level if this is incorrect."))
    }
    comparison_levels <- setdiff(var_levels, reference_level)
    if (length(comparison_levels) == 0) {
      return_list$error_message <- paste0("No comparison levels found for '", primary_variable, "' after excluding reference '", reference_level, "'.")
      return(return_list)
    }
    # Default to the first comparison level if multiple exist and none specified
    primary_variable_comparison_level <- comparison_levels[1]
    if (length(comparison_levels) > 1) {
      message(paste0("Multiple comparison levels exist for '", primary_variable, "'. Analyzing '", primary_variable_comparison_level, "' vs '", reference_level, "'. Specify 'primary_variable_comparison_level' for others."))
    }
  } else {
    if (is.null(reference_level)) {
      # If comparison level is set, but reference is not, try to infer reference
      var_levels <- levels(as.factor(get_variable(phyloseq_object, primary_variable)))
      possible_references <- setdiff(var_levels, primary_variable_comparison_level)
      if(length(possible_references) > 0) {
        reference_level <- possible_references[1] # Take the first possible one
        message(paste0("Inferring '", reference_level, "' as the reference level for '", primary_variable, "' based on comparison level '", primary_variable_comparison_level, "'."))
      } else {
        return_list$error_message <- paste0("Cannot infer reference level when primary_variable_comparison_level is '", primary_variable_comparison_level, "'.")
        return(return_list)
      }
    }
  }
  
  
  # Construct column names based on ANCOM-BC2 naming convention
  # Example: primary_variable = "carriage_group", primary_variable_comparison_level = "Acquisition_carrier"
  # Expected columns: lfc_carriage_groupAcquisition_carrier, diff_carriage_groupAcquisition_carrier
  lfc_col_name <- paste0("lfc_", primary_variable, primary_variable_comparison_level)
  diff_col_name <- paste0("diff_", primary_variable, primary_variable_comparison_level)
  ss_col_name <- paste0("passed_ss_", primary_variable, primary_variable_comparison_level)
  
  # Check if these columns exist
  if (!lfc_col_name %in% names(res_df)) {
    return_list$error_message <- paste0("LFC column '", lfc_col_name, "' not found in model results. Check primary_variable, reference_level, and primary_variable_comparison_level.")
    # Try to list available lfc columns for debugging
    avail_lfc_cols <- grep("^lfc_", names(res_df), value = TRUE)
    if(length(avail_lfc_cols) > 0) message(paste("Available LFC columns:", paste(avail_lfc_cols, collapse=", ")))
    return(return_list)
  }
  if (!diff_col_name %in% names(res_df)) {
    return_list$error_message <- paste0("Differential abundance column '", diff_col_name, "' not found.")
    return(return_list)
  }
  
  return_list$primary_comparison_details <- paste0(primary_variable, ": ", primary_variable_comparison_level, " vs ", reference_level)
  
  # --- 3. Summary Statistics ---
  total_taxa_analyzed <- nrow(res_df)
  
  # Initial significant filter based on diff_col_name
  significant_mask <- res_df[[diff_col_name]] == TRUE & !is.na(res_df[[diff_col_name]])
  
  num_significant_initial <- sum(significant_mask)
  
  # Further filter by sensitivity score if required
  if (filter_on_sensitivity) {
    if (ss_col_name %in% names(res_df)) {
      significant_mask <- significant_mask & (res_df[[ss_col_name]] == TRUE & !is.na(res_df[[ss_col_name]]))
      num_significant_final <- sum(significant_mask)
      sensitivity_note <- " (passed sensitivity score)"
    } else {
      message(paste0("Sensitivity score column '", ss_col_name, "' not found. Reporting based on diff column only."))
      num_significant_final <- num_significant_initial
      sensitivity_note <- " (sensitivity score column not found)"
    }
  } else {
    num_significant_final <- num_significant_initial
    sensitivity_note <- " (sensitivity score not filtered)"
  }
  
  return_list$summary_stats <- rbind(return_list$summary_stats,
                                     data.frame(Description = "Total taxa analyzed by ANCOM-BC2", Value = as.character(total_taxa_analyzed)),
                                     data.frame(Description = paste0("Significant taxa for ", return_list$primary_comparison_details, sensitivity_note), Value = as.character(num_significant_final)))
  
  
  # --- 4. Data Frame of Significant Taxa ---
  if (num_significant_final > 0) {
    sig_taxa_data <- res_df[significant_mask, ]
    
    # Select and rename columns for the output data frame
    # Need to get se, W, p_val, q_val columns corresponding to the lfc_col_name
    # These are typically named by removing "lfc_" and adding the prefix.
    base_comp_name <- gsub("lfc_", "", lfc_col_name)
    
    se_col_name <- paste0("se_", base_comp_name)
    w_col_name <- paste0("W_", base_comp_name)
    pval_col_name <- paste0("p_", base_comp_name)
    qval_col_name <- paste0("q_", base_comp_name)
    # diff_col_name is already defined
    # ss_col_name is already defined
    
    cols_to_select <- c("taxon", lfc_col_name)
    
    # Add other columns if they exist
    if(se_col_name %in% names(sig_taxa_data)) cols_to_select <- c(cols_to_select, se_col_name)
    if(w_col_name %in% names(sig_taxa_data)) cols_to_select <- c(cols_to_select, w_col_name)
    if(pval_col_name %in% names(sig_taxa_data)) cols_to_select <- c(cols_to_select, pval_col_name)
    if(qval_col_name %in% names(sig_taxa_data)) cols_to_select <- c(cols_to_select, qval_col_name)
    cols_to_select <- c(cols_to_select, diff_col_name) # This should always be TRUE for this subset
    if(ss_col_name %in% names(sig_taxa_data)) cols_to_select <- c(cols_to_select, ss_col_name)
    
    # Ensure no duplicate column names (taxon might be the only one guaranteed initially)
    cols_to_select <- unique(cols_to_select)
    
    processed_sig_taxa_df <- sig_taxa_data[, cols_to_select, drop = FALSE]
    
    # Rename columns for consistency
    new_colnames <- c("Taxon", "LFC")
    if(se_col_name %in% names(sig_taxa_data)) new_colnames <- c(new_colnames, "SE")
    if(w_col_name %in% names(sig_taxa_data)) new_colnames <- c(new_colnames, "W_statistic")
    if(pval_col_name %in% names(sig_taxa_data)) new_colnames <- c(new_colnames, "P_value")
    if(qval_col_name %in% names(sig_taxa_data)) new_colnames <- c(new_colnames, "Q_value")
    new_colnames <- c(new_colnames, "Differentially_Abundant") # From diff_col_name
    if(ss_col_name %in% names(sig_taxa_data)) new_colnames <- c(new_colnames, "Passed_Sensitivity_Score")
    
    # Adjust if some columns were not found
    current_selected_count = ncol(processed_sig_taxa_df)
    if(length(new_colnames) > current_selected_count) {
      new_colnames <- new_colnames[1:current_selected_count]
    }
    
    colnames(processed_sig_taxa_df) <- new_colnames[1:ncol(processed_sig_taxa_df)] # Make sure length matches
    
    return_list$significant_taxa_df <- processed_sig_taxa_df %>%
      dplyr::arrange(desc(abs(LFC)))
    
    # --- 5. Visualization ---
    plot_df_viz <- return_list$significant_taxa_df %>%
      dplyr::slice_head(n = top_n_taxa_plot) %>%
      dplyr::mutate(Direction = ifelse(LFC > 0, "Increased", "Decreased"),
                    Taxon = factor(Taxon, levels = Taxon[order(LFC)])) # Order by LFC for plot
    
    if (nrow(plot_df_viz) > 0) {
      p <- ggplot2::ggplot(plot_df_viz, ggplot2::aes(x = Taxon, y = LFC, fill = Direction)) +
        ggplot2::geom_bar(stat = "identity", color = "black") +
        ggplot2::scale_fill_manual(values = c("Increased" = "#F8766D", "Decreased" = "#619CFF")) +
        ggplot2::coord_flip() +
        ggplot2::labs(
          title = paste0("Top ", nrow(plot_df_viz), " Significant Taxa"),
          subtitle = paste0(model_name, "\nComparison: ", return_list$primary_comparison_details),
          x = "Taxon",
          y = "Log Fold Change (LFC)"
        ) +
        ggplot2::theme_minimal(base_size = 10) +
        ggplot2::theme(
          plot.title = ggplot2::element_text(hjust = 0.5, face = "bold"),
          plot.subtitle = ggplot2::element_text(hjust = 0.5),
          legend.position = "bottom"
        )
      return_list$plot_significant_taxa <- p
    }
  } else {
    message(paste0("No significant taxa found for '", return_list$primary_comparison_details, "' in model '", model_name, "' after filtering."))
  }
  
  return(return_list)
}
