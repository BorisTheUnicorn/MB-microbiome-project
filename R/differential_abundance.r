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
                         rand_formula = NULL,
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
    rand_formula = rand_formula,
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

#' Prepare Phyloseq Object for ANCOM-BC2 Analysis with Input Validation
#'
#' Subset a phyloseq object to a specific time point, filter on two comparison groups,
#' remove a specified batch, prune zero-sum samples, and re-level the grouping factor.
#'
#' @param ps_object A phyloseq object containing the microbiome data.
#' @param group_levels Character vector of length 2 specifying the comparison levels for the grouping variable.
#' @param group_var Character string naming the grouping variable in sample_data (default: "carriage_group").
#' @param time_var Character string naming the time variable in sample_data (default: "time").
#' @param time_value Value of the time_var to subset on (default: "0").
#' @param batch_var Character string naming the batch variable in sample_data (default: "batch").
#' @param drop_batch Value of the batch_var to exclude from the analysis (default: "4").
#' @return A phyloseq object filtered and relabeled according to the provided parameters.
#' @export
prepare_phyloseq <- function(ps_object,
                             group_levels,
                             group_var     = "carriage_group",
                             time_var      = "time",
                             time_value    = "0",
                             batch_var     = "batch",
                             drop_batch    = "4") {
  # Validate phyloseq object
  if (!requireNamespace("phyloseq", quietly = TRUE) || !inherits(ps_object, "phyloseq")) {
    stop("ps_object must be a valid phyloseq object.")
  }
  
  # Validate group_levels
  if (!is.character(group_levels) || length(group_levels) != 2) {
    stop("group_levels must be a character vector of length 2 (e.g., c('A','B')).")
  }
  
  # Extract sample data for checks
  sd_df <- as.data.frame(phyloseq::sample_data(ps_object), stringsAsFactors = FALSE)
  
  # Check required columns
  required_vars <- c(group_var, time_var, batch_var)
  missing_vars <- setdiff(required_vars, colnames(sd_df))
  if (length(missing_vars) > 0) {
    stop(paste("Missing required sample_data columns:", paste(missing_vars, collapse = ", ")))  
  }
  
  # Check time_value exists
  if (!time_value %in% sd_df[[time_var]]) {
    stop(sprintf("time_value '%s' not found in column '%s'.", time_value, time_var))
  }
  
  # Warn if drop_batch not present
  if (!drop_batch %in% sd_df[[batch_var]]) {
    warning(sprintf("drop_batch '%s' not found in column '%s'. No batch filtering applied.", drop_batch, batch_var))
  }
  
  # Subset samples and prune zeros
  keep_samples <- sd_df[[time_var]] == time_value &
    sd_df[[group_var]] %in% group_levels &
    sd_df[[batch_var]] != drop_batch
  ps_sub <- phyloseq::prune_samples(keep_samples, ps_object)
  ps_sub <- phyloseq::prune_taxa(phyloseq::taxa_sums(ps_sub) > 0, ps_sub)
  
  # Ensure samples remain
  if (phyloseq::nsamples(ps_sub) == 0) {
    stop("No samples remain after subsetting. Check group_levels, time_value, and batch filters.")
  }
  
  # Re-level grouping factor
  sd_sub <- as.data.frame(phyloseq::sample_data(ps_sub), stringsAsFactors = FALSE)
  sd_sub[[group_var]] <- factor(sd_sub[[group_var]], levels = group_levels)
  phyloseq::sample_data(ps_sub) <- sd_sub
  
  return(ps_sub)
}


#' Run Multiple ANCOM-BC2 Models with Input Validation
#'
#' Loop over a named list of formula strings to run ANCOM-BC2 on a prepared phyloseq object.
#'
#' @param ps_object A phyloseq object that has been prepared (e.g., via prepare_phyloseq()).
#' @param formulas Named list of formula strings. Names will be used in the returned list.
#' @param group_var Character string naming the grouping variable (default: "carriage_group").
#' @param tax_level Taxonomic level for analysis (default: "Genus"). Must exist in phyloseq::rank_names(ps_object).
#' @param prv_cut Prevalence cutoff (a numeric fraction between 0 and 1) for filtering taxa (default: 0.025) ([rdrr.io](https://rdrr.io/github/FrederickHuangLin/ANCOMBC/man/ancombc2.html?utm_source=chatgpt.com)).
#' @param seed Single numeric seed for reproducibility (default: 123).
#' @return A named list of ANCOM-BC2 result objects.
#' @export
run_ancombc2_models <- function(ps_object,
                                formulas,
                                rand_formula = NULL,
                                group_var = "carriage_group",
                                tax_level = "Genus",
                                prv_cut   = 0.025,
                                seed      = 123) {
  # Validate ANCOM wrapper exists
  if (!exists("run_ancombc2", mode = "function")) {
    stop("Function 'run_ancombc2' must be defined and sourced before calling run_ancombc2_models.")
  }
  
  # Validate phyloseq object
  if (!requireNamespace("phyloseq", quietly = TRUE) || !inherits(ps_object, "phyloseq")) {
    stop("ps_object must be a valid phyloseq object.")
  }
  
  # Validate formulas input
  if (!is.list(formulas) || is.null(names(formulas))) {
    stop("'formulas' must be a named list of non-empty formula strings.")
  }
  for (nm in names(formulas)) {
    f <- formulas[[nm]]
    if (!is.character(f) || length(f) != 1 || nchar(f) == 0) {
      stop(sprintf("Formula '%s' must be a non-empty string.", nm))
    }
  }
  
  # Check grouping variable in sample_data
  sd_df <- as.data.frame(phyloseq::sample_data(ps_object), stringsAsFactors = FALSE)
  if (!group_var %in% colnames(sd_df)) {
    stop(sprintf("group_var '%s' not found in sample_data(ps_object).", group_var))
  }
  
  # Validate taxonomic level
  available_ranks <- phyloseq::rank_names(ps_object)
  if (!tax_level %in% available_ranks) {
    stop(sprintf("tax_level '%s' not found. Available ranks: %s",
                 tax_level, paste(available_ranks, collapse = ", ")))
  }
  
  # Validate prevalence cutoff
  if (!is.numeric(prv_cut) || prv_cut <= 0 || prv_cut >= 1) {
    stop("prv_cut must be a numeric fraction strictly between 0 and 1.")
  }
  
  # Validate seed
  if (!is.numeric(seed) || length(seed) != 1) {
    stop("seed must be a single numeric value.")
  }
  
  # Run models
  results <- vector("list", length(formulas))
  names(results) <- names(formulas)
  for (nm in names(formulas)) {
    results[[nm]] <- run_ancombc2(
      ps_object   = ps_object,
      tax_level   = tax_level,
      formula_str = formulas[[nm]],
      rand_formula = rand_formula,
      group_var   = group_var,
      prv_cut     = prv_cut,
      seed        = seed
    )
  }
  
  return(results)
}

######################## CLAUDE GENERAL ANALYSIS FUNCTION  ######################## 

#' Analyze ANCOM-BC2 Results with Summary Tables and Visualizations
#'
#' This function provides a comprehensive analysis of ANCOM-BC2 model results,
#' including summary tables of significant taxa and comparative visualizations.
#'
#' @param ancom_model Output from ancombc2() function
#' @param model_name Character string to identify the model in outputs
#' @param alpha Significance threshold (default: 0.05)
#' @param min_abs_lfc Minimum absolute log fold change to consider (default: 0)
#' @param save_plots Logical, whether to save plots to file (default: FALSE)
#' @param plot_dir Directory to save plots if save_plots = TRUE
#' @return A list containing summary tables and plots
#' @export

analyze_ancom <- function(ancom_model, 
                          model_name = "ANCOM-BC2 Model",
                          alpha = 0.05,
                          min_abs_lfc = 0,
                          save_plots = FALSE,
                          plot_dir = NULL) {
  
  # Check if results exist
  if (is.null(ancom_model$res) || !is.data.frame(ancom_model$res)) {
    stop("No valid results found in ancom_model$res")
  }
  
  res_df <- ancom_model$res
  
  # Identify all contrasts by finding lfc_ columns
  lfc_cols <- grep("^lfc_", colnames(res_df), value = TRUE)
  contrasts <- gsub("^lfc_", "", lfc_cols)
  
  # Initialize storage for results
  all_sig_taxa <- list()
  summary_stats <- data.frame()
  
  # Process each contrast
  for (i in seq_along(contrasts)) {
    contrast <- contrasts[i]
    lfc_col <- paste0("lfc_", contrast)
    diff_col <- paste0("diff_", contrast)
    ss_col <- paste0("passed_ss_", contrast)
    q_col <- paste0("q_", contrast)
    se_col <- paste0("se_", contrast)
    
    # Skip if essential columns don't exist
    if (!all(c(lfc_col, diff_col) %in% colnames(res_df))) next
    
    # Extract significant taxa
    sig_mask <- res_df[[diff_col]] == TRUE & !is.na(res_df[[diff_col]])
    
    # Apply sensitivity filter if column exists
    if (ss_col %in% colnames(res_df)) {
      sig_mask <- sig_mask & res_df[[ss_col]] == TRUE & !is.na(res_df[[ss_col]])
    }
    
    # Apply minimum LFC filter
    sig_mask <- sig_mask & abs(res_df[[lfc_col]]) >= min_abs_lfc
    
    # Get significant taxa for this contrast
    sig_taxa <- res_df[sig_mask, ]
    
    if (nrow(sig_taxa) > 0) {
      # Create summary for this contrast
      contrast_summary <- data.frame(
        contrast = contrast,
        taxon = sig_taxa$taxon,
        lfc = sig_taxa[[lfc_col]],
        stringsAsFactors = FALSE
      )
      
      # Add other columns if they exist
      if (se_col %in% colnames(sig_taxa)) {
        contrast_summary$se <- sig_taxa[[se_col]]
      }
      if (q_col %in% colnames(sig_taxa)) {
        contrast_summary$q_value <- sig_taxa[[q_col]]
      }
      
      all_sig_taxa[[contrast]] <- contrast_summary
      
      # Summary statistics
      summary_stats <- rbind(summary_stats, data.frame(
        contrast = contrast,
        n_significant = nrow(sig_taxa),
        n_increased = sum(sig_taxa[[lfc_col]] > 0),
        n_decreased = sum(sig_taxa[[lfc_col]] < 0),
        mean_abs_lfc = mean(abs(sig_taxa[[lfc_col]])),
        max_abs_lfc = max(abs(sig_taxa[[lfc_col]]))
      ))
    }
  }
  
  # Combine all significant taxa into one data frame
  if (length(all_sig_taxa) > 0) {
    all_sig_df <- do.call(rbind, all_sig_taxa)
    rownames(all_sig_df) <- NULL
  } else {
    message("No significant taxa found in any contrast.")
    return(list(summary_stats = NULL, significant_taxa = NULL, plots = list()))
  }
  
  # Print summary
  cat("\n=== ANCOM-BC2 Analysis Summary for:", model_name, "===\n")
  cat("\nTotal contrasts analyzed:", length(contrasts), "\n")
  cat("Contrasts with significant taxa:", nrow(summary_stats), "\n")
  cat("\nSummary by contrast:\n")
  print(kable(summary_stats, digits = 3) %>% kable_styling())
  
  cat("\nTop 10 taxa by absolute LFC across all contrasts:\n")
  top_taxa <- all_sig_df %>%
    arrange(desc(abs(lfc))) %>%
    head(10)
  print(kable(top_taxa, digits = 3) %>% kable_styling())
  
  # Create visualizations
  plots <- list()
  
  # 1. Heatmap of significant taxa across contrasts
  if (length(unique(all_sig_df$contrast)) > 1) {
    # Prepare data for heatmap
    heatmap_data <- all_sig_df %>%
      select(taxon, contrast, lfc) %>%
      tidyr::pivot_wider(names_from = contrast, values_from = lfc, values_fill = 0)
    
    # Convert to matrix for heatmap
    heatmap_matrix <- as.matrix(heatmap_data[, -1])
    rownames(heatmap_matrix) <- heatmap_data$taxon
    
    # Create heatmap using ggplot2
    heatmap_long <- all_sig_df %>%
      select(taxon, contrast, lfc)
    
    # Order taxa by hierarchical clustering if there are multiple contrasts
    if (ncol(heatmap_matrix) > 1 && nrow(heatmap_matrix) > 1) {
      taxa_order <- hclust(dist(heatmap_matrix))$order
      heatmap_long$taxon <- factor(heatmap_long$taxon, 
                                   levels = rownames(heatmap_matrix)[taxa_order])
    }
    
    p_heatmap <- ggplot(heatmap_long, aes(x = contrast, y = taxon, fill = lfc)) +
      geom_tile(color = "white", size = 0.5) +
      scale_fill_gradient2(low = "blue", mid = "white", high = "red", 
                           midpoint = 0, name = "Log Fold\nChange") +
      theme_minimal() +
      theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
            axis.text.y = element_text(size = 8),
            plot.title = element_text(hjust = 0.5, face = "bold")) +
      labs(title = paste(model_name, "- Significant Taxa Heatmap"),
           x = "Contrast", y = "Taxon")
    
    plots$heatmap <- p_heatmap
    print(p_heatmap)
  }
  
  # 2. Dot plot with confidence intervals (if SE available)
  if ("se" %in% colnames(all_sig_df) && all(!is.na(all_sig_df$se))) {
    # Calculate confidence intervals
    all_sig_df$lfc_lower <- all_sig_df$lfc - 1.96 * all_sig_df$se
    all_sig_df$lfc_upper <- all_sig_df$lfc + 1.96 * all_sig_df$se
    
    # Select top taxa by absolute LFC for clarity
    top_taxa_for_plot <- all_sig_df %>%
      group_by(taxon) %>%
      summarise(max_abs_lfc = max(abs(lfc))) %>%
      arrange(desc(max_abs_lfc)) %>%
      head(20) %>%
      pull(taxon)
    
    plot_data <- all_sig_df %>%
      filter(taxon %in% top_taxa_for_plot) %>%
      mutate(taxon = factor(taxon, levels = rev(top_taxa_for_plot)))
    
    p_dotplot <- ggplot(plot_data, aes(x = lfc, y = taxon, color = contrast)) +
      geom_vline(xintercept = 0, linetype = "dashed", alpha = 0.5) +
      geom_errorbarh(aes(xmin = lfc_lower, xmax = lfc_upper), 
                     height = 0.2, alpha = 0.7) +
      geom_point(size = 3) +
      theme_bw() +
      theme(legend.position = "bottom",
            plot.title = element_text(hjust = 0.5, face = "bold")) +
      labs(title = paste(model_name, "- Effect Sizes with 95% CI"),
           x = "Log Fold Change", y = "Taxon", color = "Contrast")
    
    plots$dotplot <- p_dotplot
    print(p_dotplot)
  }
  
  # 3. Volcano-style plot for each contrast (FIXED VERSION)
  for (contrast in unique(all_sig_df$contrast)) {
    contrast_data <- res_df
    lfc_col <- paste0("lfc_", contrast)
    q_col <- paste0("q_", contrast)
    diff_col <- paste0("diff_", contrast)
    ss_col <- paste0("passed_ss_", contrast)
    
    if (all(c(lfc_col, q_col) %in% colnames(contrast_data))) {
      # Pre-calculate -log10(q-value) to avoid issues with aes_string
      contrast_data$neg_log10_q <- -log10(contrast_data[[q_col]])
      
      # Handle infinite values (when q-value = 0)
      max_finite <- max(contrast_data$neg_log10_q[is.finite(contrast_data$neg_log10_q)], na.rm = TRUE)
      contrast_data$neg_log10_q[!is.finite(contrast_data$neg_log10_q)] <- max_finite * 1.1
      
      # Create significance categories
      contrast_data$significance <- "Not Significant"
      
      if (diff_col %in% colnames(contrast_data)) {
        contrast_data$significance[contrast_data[[diff_col]] == TRUE] <- "Significant"
      }
      
      if (ss_col %in% colnames(contrast_data) && diff_col %in% colnames(contrast_data)) {
        contrast_data$significance[contrast_data[[diff_col]] == TRUE & 
                                     contrast_data[[ss_col]] == TRUE] <- "Sig. & Pass SS"
      }
      
      # Create volcano plot using modern tidy evaluation
      p_volcano <- ggplot(contrast_data, 
                          aes(x = !!sym(lfc_col), 
                              y = neg_log10_q,
                              color = significance)) +
        geom_point(alpha = 0.6, size = 2) +
        geom_vline(xintercept = c(-min_abs_lfc, min_abs_lfc), 
                   linetype = "dashed", alpha = 0.5) +
        geom_hline(yintercept = -log10(alpha), 
                   linetype = "dashed", alpha = 0.5) +
        scale_color_manual(values = c("Not Significant" = "grey60",
                                      "Significant" = "orange",
                                      "Sig. & Pass SS" = "red")) +
        theme_bw() +
        theme(legend.position = "bottom",
              plot.title = element_text(hjust = 0.5, face = "bold")) +
        labs(title = paste(model_name, "-", contrast, "Volcano Plot"),
             x = "Log Fold Change",
             y = expression(-log[10](q-value)),
             color = "Status")
      
      # Add labels for top taxa
      top_taxa_volcano <- contrast_data %>%
        filter(significance == "Sig. & Pass SS") %>%
        arrange(desc(abs(.data[[lfc_col]]))) %>%
        head(10)
      
      if (nrow(top_taxa_volcano) > 0) {
        p_volcano <- p_volcano +
          ggrepel::geom_text_repel(data = top_taxa_volcano,
                                   aes(label = taxon),
                                   size = 3,
                                   box.padding = 0.3,
                                   max.overlaps = 15)
      }
      
      plots[[paste0("volcano_", contrast)]] <- p_volcano
      print(p_volcano)
    }
  }
  
  # Save plots if requested
  if (save_plots && !is.null(plot_dir)) {
    if (!dir.exists(plot_dir)) dir.create(plot_dir, recursive = TRUE)
    
    for (plot_name in names(plots)) {
      filename <- file.path(plot_dir, paste0(model_name, "_", plot_name, ".png"))
      ggsave(filename, plots[[plot_name]], width = 10, height = 8, dpi = 300)
    }
    message("Plots saved to:", plot_dir)
  }
  
  # Return results
  return(list(
    summary_stats = summary_stats,
    significant_taxa = all_sig_df,
    plots = plots,
    full_results = res_df
  ))
}