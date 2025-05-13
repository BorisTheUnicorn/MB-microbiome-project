#' Visualization Functions for Microbiome Analysis
#' This file contains functions for creating various visualizations for
#' microbiome data, particularly for the NM carriage study.
#' 
#' @author Orr Tobaly & Nikol Elyashov
#' @date April 2025

# Required packages --------------------------------------------------------
# library(ggplot2)
# library(patchwork)
# library(viridis)
# library(microbiome)
# library(phyloseq)
# library(dplyr)
# library(tidyr)
# library(ggpubr)

#' Create Beta Diversity Ordination Plot with Enhanced Aesthetics
#'
#' Creates a beta diversity ordination plot using microeco's built-in functions,
#' allowing separate variables for color and shape, and control over point/ellipse aesthetics.
#'
#' @param dataset A microtable object.
#' @param plot_color_var Factor name (from sample_table) to use for coloring points/ellipses.
#' @param plot_shape_var Factor name (from sample_table) to use for point shapes (optional).
#'                         If NULL, shape may default to plot_color_var or a single shape.
#' @param distance_metric Beta diversity measure (default: "bray").
#' @param method Ordination method (default: "PCoA").
#' @param plot_type Plot type (default: c("point", "ellipse")). Can be a vector.
#' @param title Plot title (optional). If NULL, a title will be auto-generated.
#' @param point_size Size of the points in the plot (default: 3). This is passed to microeco's point_size.
#' @param ellipse_alpha Alpha (transparency) for ellipses (default: 0.1). This is passed to microeco's ellipse_chull_alpha.
#' @param color_values Custom color palette for plot_color_var (optional).
#' @param shape_values Custom shape palette for plot_shape_var (optional).
#' @return A ggplot object with the ordination plot.
#' @export
plot_pcoa <- function(dataset,
                      plot_color_var,
                      plot_shape_var = NULL,
                      distance_metric = "bray",
                      method = "PCoA",
                      plot_type = c("point", "ellipse"),
                      title = NULL,
                      point_size = 3,        # User-facing parameter for the wrapper
                      ellipse_alpha = 0.1,   # User-facing parameter for the wrapper
                      color_values = NULL,
                      shape_values = NULL) {
  
  # Ensure the microeco package is available
  # if (!requireNamespace("microeco", quietly = TRUE)) {
  #   stop("Please install the 'microeco' package to use this function.")
  # }
  # if (!requireNamespace("ggplot2", quietly = TRUE)) {
  #   stop("Please install the 'ggplot2' package.")
  # }
  
  # Calculate beta diversity if not already present for the chosen metric
  if (is.null(dataset$beta_diversity[[distance_metric]])) {
    dataset$cal_betadiv(method = distance_metric)
  }
  
  beta_obj <- microeco::trans_beta$new(dataset = dataset,
                                       measure = distance_metric,
                                       group = plot_color_var) # Group primarily for ellipses
  
  beta_obj$cal_ordination(method = method)
  
  if(is.null(title)) {
    title_parts <- paste0("Beta Diversity (", toupper(method), " on ", distance_metric, ") by ", plot_color_var)
    if (!is.null(plot_shape_var) && plot_color_var != plot_shape_var) {
      title_parts <- paste0(title_parts, " and Shape by ", plot_shape_var)
    }
    title <- title_parts
  }
  
  current_plot_shape <- if (is.null(plot_shape_var)) plot_color_var else plot_shape_var
  
  # Create the plot using microeco's expected argument names
  beta_plot <- beta_obj$plot_ordination(
    plot_color = plot_color_var,
    plot_shape = current_plot_shape,
    plot_type = plot_type,
    point_size = point_size,             # Corrected: Pass to microeco's 'point_size'
    ellipse_chull_alpha = ellipse_alpha  # Corrected: Pass to microeco's 'ellipse_chull_alpha'
  ) +
    labs(title = title) +
    theme_bw()
  
  # Apply custom color palette if provided
  if (!is.null(color_values)) {
    beta_plot <- beta_plot + ggplot2::scale_color_manual(values = color_values)
    if (any(c("ellipse", "chull") %in% plot_type)) {
      beta_plot <- beta_plot + ggplot2::scale_fill_manual(values = color_values)
    }
  }
  
  # Apply custom shape palette if provided and shape is mapped to a variable
  if (!is.null(shape_values) && !is.null(plot_shape_var)) {
    beta_plot <- beta_plot + ggplot2::scale_shape_manual(values = shape_values)
  }
  
  return(beta_plot)
}

#' Compare Pre and Post Correction Plots
#'
#' Creates side-by-side comparison plots of data before and after batch correction
#'
#' @param dataset_original Original microtable object
#' @param dataset_corrected Batch-corrected microtable object
#' @param factor_names Vector of factor names to create ordination plots for
#' @param plot_titles Vector of plot titles (optional)
#' @return A list of comparison plots
#' @export
compare_correction_plots <- function(dataset_original, dataset_corrected, 
                                     factor_names, plot_titles = NULL) {
  # Initialize list for plots
  comparison_plots <- list()
  
  # If plot titles not provided, use factor names
  if(is.null(plot_titles)) {
    plot_titles <- paste("Batch Effect by", factor_names)
  }
  
  # Make sure we have the right number of titles
  if(length(plot_titles) != length(factor_names)) {
    plot_titles <- rep(plot_titles, length.out = length(factor_names))
  }
  
  # Create comparison plots for each factor
  for(i in seq_along(factor_names)) {
    factor_name <- factor_names[i]
    
    # Create pre-correction ordination
    pre_plot <- plot_pcoa(
      dataset = dataset_original, 
      factor_name = factor_name,
      title = paste("Pre-correction -", plot_titles[i])
    )
    
    # Create post-correction ordination
    post_plot <- plot_pcoa(
      dataset = dataset_corrected, 
      factor_name = factor_name,
      title = paste("Post-correction -", plot_titles[i])
    )
    
    # Combine plots using patchwork
    comparison_plots[[factor_name]] <- pre_plot + post_plot
  }
  
  return(comparison_plots)
}

#' Plot Alpha Diversity by Group
#'
#' Creates boxplots of alpha diversity indices grouped by a specified variable
#'
#' @param dataset A microtable object
#' @param group_var The grouping variable name
#' @param measures Vector of alpha diversity indices to plot (default: c("Shannon", "Simpson", "Observed"))
#' @param title Plot title (optional)
#' @return A ggplot object with alpha diversity plots
#' @export
plot_alpha_diversity <- function(dataset, group_var, 
                                 measures = c("Shannon", "Simpson", "Observed"),
                                 title = NULL) {
  # Make sure alpha diversity has been calculated
  if(is.null(dataset$alpha_diversity)) {
    dataset$cal_alphadiv(measures = measures)
  }
  
  # Create alpha_div object
  alpha_div_obj <- trans_alpha$new(dataset = dataset, group = group_var)
  
  # Calculate differential test (not mandatory but adds statistical significance)
  alpha_div_obj$cal_diff(method = "wilcox")
  
  # Generate boxplot - using measure parameter here, not during initialization
  alpha_plot <- alpha_div_obj$plot_alpha(
    measure = measures[1],  # Start with first measure
    plot_type = "ggboxplot",
    add_sig = TRUE
  )
  
  # Add title if provided
  if (!is.null(title)) {
    alpha_plot <- alpha_plot + labs(title = title)
  }
  
  return(alpha_plot)
}

#' Create a Stacked Bar Plot of Taxa Composition
#'
#' Creates a stacked bar plot showing the relative abundance of the top taxa.
#' Can handle single or multiple faceting variables.
#'
#' @param dataset A microtable object.
#' @param tax_level Taxonomic level to aggregate to (default: "Genus").
#' @param facet_vars A character string for a single faceting variable OR 
#'                   a character vector for multiple faceting variables (e.g., c("carriage_group", "batch")).
#'                   These should be column names in dataset$sample_table.
#' @param n_taxa Number of top taxa to show individually (default: 15).
#' @param title Plot title (optional).
#' @return A ggplot object with compositional bar plot.
#' @export
plot_taxa_composition <- function(dataset, tax_level = "Genus", facet_vars = NULL,
                                  n_taxa = 15, title = NULL) {
  # Create a clone to avoid modifying original dataset
  dataset_agg <- dataset$clone()
  
  # Ensure faceting variables are factors for ordered plotting if desired by user beforehand
  if (!is.null(facet_vars)) {
    for (var in facet_vars) {
      if (var %in% colnames(dataset_agg$sample_table) && !is.factor(dataset_agg$sample_table[[var]])) {
        # message(paste0("Converting '", var, "' to factor for faceting.")) # Optional message
        dataset_agg$sample_table[[var]] <- as.factor(dataset_agg$sample_table[[var]])
      }
    }
  }
  
  # Merge taxa to the desired taxonomic level if needed
  if (tax_level != "OTU") {
    dataset_agg$merge_taxa(taxa = tax_level)
  }
  
  # Create a trans_abund object
  # The trans_abund$new() function itself does not take grouping/faceting variables for the plot
  # It prepares the abundance data.
  trans_abund_obj <- trans_abund$new(
    dataset = dataset_agg, 
    taxrank = tax_level, 
    ntaxa = n_taxa
  )
  
  # Generate the barplot
  # The `facet` argument in `plot_bar` will receive the character vector directly
  taxa_plot <- trans_abund_obj$plot_bar(
    facet = facet_vars, # Pass the single string or character vector here
    others_color = "grey70", 
    xtext_keep = FALSE 
    # Consider adding other plot_bar arguments if needed, e.g., for ordering within facets if microeco supports it
  )
  
  # Add title if provided
  if (!is.null(title)) {
    taxa_plot <- taxa_plot + labs(title = title)
  }
  
  return(taxa_plot)
}

#' Create a Longitudinal PCoA Plot
#'
#' Creates a PCoA plot showing trajectories of samples across time points
#'
#' @param dataset A microtable object
#' @param subject_var Variable indicating subject ID (default: "pid")
#' @param time_var Variable indicating time point (default: "time")
#' @param group_var Variable for grouping/coloring subjects
#' @param title Plot title (optional)
#' @return A ggplot object with longitudinal PCoA plot
#' @export
plot_longitudinal_pcoa <- function(dataset, subject_var = "pid", time_var = "time",
                                   group_var, title = NULL) {
  # Make sure beta diversity has been calculated
  if(is.null(dataset$beta_diversity)) {
    dataset$cal_betadiv(method = "bray")
  }
  
  # Create beta diversity object
  beta_obj <- trans_beta$new(dataset = dataset, 
                             measure = "bray", 
                             group = group_var)
  
  # Run ordination
  beta_obj$cal_ordination(method = "PCoA")
  
  # Extract ordination data - handle column names correctly
  scores_df <- as.data.frame(beta_obj$res_ordination$scores)
  
  # Get the actual column names (they may not be PC1, PC2)
  pc_cols <- colnames(scores_df)[1:2]
  
  # Add sample IDs as column
  scores_df$SampleID <- rownames(scores_df)
  
  # Add metadata
  metadata <- dataset$sample_table
  metadata_rownames <- rownames(metadata)
  if(length(metadata_rownames) == 0 || any(metadata_rownames == "")) {
    # If rownames aren't set properly in metadata, use first column
    rownames(metadata) <- metadata[,1]
  }
  
  # Join metadata to scores
  scores_df[[subject_var]] <- metadata[scores_df$SampleID, subject_var]
  scores_df[[time_var]] <- metadata[scores_df$SampleID, time_var]
  scores_df[[group_var]] <- metadata[scores_df$SampleID, group_var]
  
  # Convert time to factor if it's not already
  scores_df[[time_var]] <- as.factor(scores_df[[time_var]])
  
  # Extract variance explained
  eig <- beta_obj$res_ordination$model$values$Eigenvalues
  var_explained <- round(100 * eig / sum(eig), 1)
  
  # Create the plot using the actual column names
  p <- ggplot(scores_df, aes_string(x = pc_cols[1], y = pc_cols[2], 
                                    group = subject_var, 
                                    color = group_var)) +
    geom_line(alpha = 0.5) +
    geom_point(aes_string(shape = time_var), size = 3) +
    labs(
      x = paste0(pc_cols[1], " (", var_explained[1], "%)"),
      y = paste0(pc_cols[2], " (", var_explained[2], "%)"),
      title = title,
      color = group_var,
      shape = time_var
    ) +
    theme_bw() +
    theme(
      legend.position = "right",
      plot.title = element_text(hjust = 0.5)
    )
  
  return(p)
}
