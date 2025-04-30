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

#' Create Beta Diversity Ordination
#'
#' Creates beta diversity ordination plot using microeco's built-in functions
#'
#' @param dataset A microtable object
#' @param factor_name Factor to color/group by in the plot
#' @param distance_metric Beta diversity measure (default: "bray")
#' @param method Ordination method (default: "PCoA")
#' @param plot_type Plot type (default: c("point", "ellipse"))
#' @param title Plot title (optional)
#' @param color_values Color values for plotting (optional)
#' @return A ggplot object with the ordination plot
#' @export
plot_pcoa <- function(dataset, 
                      factor_name, 
                      distance_metric = "bray", 
                      method = "PCoA",
                      plot_type = c("point", "ellipse"),
                      title = NULL,
                      color_values = NULL) {
  
  # Create the trans_beta object
  beta_obj <- trans_beta$new(dataset = dataset, 
                             measure = distance_metric, 
                             group = factor_name)
  
  # Run the ordination
  beta_obj$cal_ordination(method = method)
  
  # Generate the plot title if not provided
  if(is.null(title)) {
    title <- paste("Beta Diversity by", factor_name)
  }
  
  # Create the plot
  beta_plot <- beta_obj$plot_ordination(
    plot_type = plot_type,
    plot_color = factor_name,
    plot_shape = factor_name
  ) + 
    labs(title = title)
  
  # Return the plot
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
#' Creates a stacked bar plot showing the relative abundance of the top taxa
#'
#' @param dataset A microtable object
#' @param tax_level Taxonomic level to aggregate to (default: "Genus")
#' @param group_var The grouping variable name
#' @param n_taxa Number of top taxa to show individually (default: 15)
#' @param title Plot title (optional)
#' @return A ggplot object with compositional bar plot
#' @export
plot_taxa_composition <- function(dataset, tax_level = "Genus", group_var,
                                  n_taxa = 15, title = NULL) {
  # Create a clone to avoid modifying original dataset
  dataset_agg <- dataset$clone()
  
  # Merge taxa to the desired taxonomic level if needed
  if(tax_level != "OTU") {
    dataset_agg$merge_taxa(taxa = tax_level)
  }
  
  # Create a trans_abund object
  trans_abund_obj <- trans_abund$new(
    dataset = dataset_agg, 
    taxrank = tax_level, 
    ntaxa = n_taxa
  )
  
  # Generate the barplot
  taxa_plot <- trans_abund_obj$plot_bar(
    facet = group_var,
    others_color = "grey70", 
    xtext_keep = FALSE
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
