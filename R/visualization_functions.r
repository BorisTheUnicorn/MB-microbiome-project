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

#' Create Beta Diversity Ordination Plot with Enhanced Aesthetics and Optional Faceting
#'
#' Creates a beta diversity ordination plot. Uses microeco's built-in functions for non-faceted plots,
#' or ggplot2 directly for faceted plots. Allows separate variables for color and shape.
#'
#' @param dataset A microtable object.
#' @param plot_color_var Factor name (from sample_table) to use for coloring points/ellipses.
#' @param plot_shape_var Factor name (from sample_table) to use for point shapes (optional).
#' @param facet_var Factor name (from sample_table) to use for faceting the plot (optional).
#' @param distance_metric Beta diversity measure (default: "bray").
#' @param method Ordination method (default: "PCoA").
#' @param plot_type Plot type (default: c("point", "ellipse")). Can be a vector.
#' @param title Plot title (optional). If NULL, a title will be auto-generated.
#' @param point_size Size of the points in the plot (default: 3).
#' @param ellipse_alpha Alpha (transparency) for ellipses (default: 0.1).
#' @param color_values Custom color palette for plot_color_var (optional).
#' @param shape_values Custom shape palette for plot_shape_var (optional).
#' @param facet_scales Argument passed to facet_wrap's scales (e.g., "free", "fixed", "free_x", "free_y"). Default "free".
#' @param facet_nrow Number of rows for facet_wrap. Default NULL (ggplot2 decides).
#' @param facet_ncol Number of columns for facet_wrap. Default NULL (ggplot2 decides).
#' @param axes_choice Numeric vector of length 2 indicating which ordination axes to plot (default: c(1, 2)).
#' @return A ggplot object with the ordination plot.
#' @export
plot_pcoa <- function(dataset,
                      plot_color_var,
                      plot_shape_var = NULL,
                      facet_var = NULL,
                      distance_metric = "bray",
                      method = "PCoA",
                      plot_type = c("point", "ellipse"),
                      title = NULL,
                      point_size = 3,
                      ellipse_alpha = 0.1,
                      color_values = NULL,
                      shape_values = NULL,
                      facet_scales = "free",
                      facet_nrow = NULL,
                      facet_ncol = NULL,
                      axes_choice = c(1, 2)) { # New parameter for choosing axes
  
  if (!requireNamespace("microeco", quietly = TRUE)) stop("Please install the 'microeco' package.")
  if (!requireNamespace("ggplot2", quietly = TRUE)) stop("Please install the 'ggplot2' package.")
  if (!requireNamespace("dplyr", quietly = TRUE)) stop("Please install the 'dplyr' package.")
  
  # Ensure chosen axes are valid
  if(length(axes_choice) != 2 || !is.numeric(axes_choice) || any(axes_choice < 1)){
    stop("axes_choice must be a numeric vector of two positive integers, e.g., c(1, 2).")
  }
  
  # Calculate beta diversity if not already present
  # microeco stores measures in dataset$beta_diversity as a list
  if (is.null(dataset$beta_diversity[[distance_metric]])) {
    message(paste0("Calculating ", distance_metric, " beta diversity..."))
    dataset$cal_betadiv(method = distance_metric)
  }
  
  beta_obj <- microeco::trans_beta$new(dataset = dataset,
                                       measure = distance_metric,
                                       group = plot_color_var) # Group for microeco's default ellipse logic
  
  beta_obj$cal_ordination(method = method, ncomp = max(axes_choice)) # Ensure enough components are calculated
  
  # Auto-generate title if not provided
  ordination_method_name <- toupper(method)
  if(is.null(title)) {
    title_parts <- paste0(ordination_method_name, " on ", distance_metric, " by ", plot_color_var)
    if (!is.null(plot_shape_var) && plot_color_var != plot_shape_var) {
      title_parts <- paste0(title_parts, " (Shape: ", plot_shape_var, ")")
    }
    if (!is.null(facet_var)) {
      title_parts <- paste0(title_parts, ", Faceted by ", facet_var)
    }
    title <- title_parts
  }
  
  # If not faceting, use microeco's plotting (generally simpler for basic cases)
  if (is.null(facet_var)) {
    current_plot_shape <- if (is.null(plot_shape_var)) plot_color_var else plot_shape_var
    
    beta_plot <- beta_obj$plot_ordination(
      plot_color = plot_color_var,
      plot_shape = current_plot_shape,
      plot_type = plot_type,
      point_size = point_size,
      ellipse_chull_alpha = ellipse_alpha,
      choices = axes_choice # Pass axes choice to microeco
    ) +
      labs(title = title) +
      theme_bw()
    
    if (!is.null(color_values)) {
      beta_plot <- beta_plot + ggplot2::scale_color_manual(values = color_values)
      if (any(c("ellipse", "chull") %in% plot_type)) {
        beta_plot <- beta_plot + ggplot2::scale_fill_manual(values = color_values)
      }
    }
    if (!is.null(shape_values) && !is.null(plot_shape_var)) {
      beta_plot <- beta_plot + ggplot2::scale_shape_manual(values = shape_values)
    }
    
  } else { # If faceting, construct plot with ggplot2 directly
    
    ordination_scores_all <- as.data.frame(beta_obj$res_ordination$scores)
    if(ncol(ordination_scores_all) < max(axes_choice)){
      stop(paste0("Requested axes (", paste(axes_choice, collapse=", "), ") exceed available ordination components (", ncol(ordination_scores_all), "). Rerun cal_ordination with higher ncomp."))
    }
    ordination_scores <- ordination_scores_all[, axes_choice, drop = FALSE]
    colnames(ordination_scores) <- c("AxisX", "AxisY") # Generic names for plotting
    
    ordination_scores$SampleID <- rownames(ordination_scores_all) # Use original rownames as SampleID
    
    sample_info <- dataset$sample_table
    # Ensure SampleID column exists in sample_info for robust merging
    if (!"SampleID" %in% colnames(sample_info)) {
      if (all(rownames(sample_info) %in% ordination_scores$SampleID)) {
        sample_info$SampleID <- rownames(sample_info)
      } else if (all(dataset$sample_table[[1]] %in% ordination_scores$SampleID)) {
        sample_info$SampleID <- as.character(dataset$sample_table[[1]])
        if(any(duplicated(sample_info$SampleID))) stop("First column of sample_table used as SampleID contains duplicates.")
        message("Using the first column of sample_table as SampleID for merging.")
      } else {
        stop("Cannot find a reliable 'SampleID' column in sample_table, and rownames do not match. Please ensure a 'SampleID' column that matches ordination score rownames exists.")
      }
    }
    
    # Convert SampleID columns to character for robust join
    ordination_scores$SampleID <- as.character(ordination_scores$SampleID)
    sample_info$SampleID <- as.character(sample_info$SampleID)
    
    plot_data <- dplyr::left_join(ordination_scores, sample_info, by = "SampleID")
    
    if(nrow(plot_data) != nrow(ordination_scores)){
      warning("Not all samples in ordination scores were matched in metadata. Check SampleID consistency.")
    }
    # Check for NAs in critical plotting variables after merge
    critical_vars_for_check <- c(plot_color_var, plot_shape_var, facet_var)
    critical_vars_for_check <- critical_vars_for_check[!sapply(critical_vars_for_check, is.null)]
    for(v in critical_vars_for_check){
      if(any(is.na(plot_data[[v]]))) {
        warning(paste("NAs introduced or present in plotting variable '", v, "' after metadata join. This can cause issues with plotting or faceting.", sep=""))
      }
    }
    
    
    # Extract explained variance for axis labels
    eig_values <- NULL
    # Try to get eigenvalues from standard PCoA (ape::pcoa via microeco)
    if(!is.null(beta_obj$res_ordination$values) && !is.null(beta_obj$res_ordination$values$Eigenvalues)){
      eig_values <- beta_obj$res_ordination$values$Eigenvalues
    } else if (!is.null(beta_obj$res_ordination$eig)) { # For vegan ordinations (rda, cca)
      eig_values <- beta_obj$res_ordination$eig
    } else if (!is.null(beta_obj$res_ordination$sdev)) { # For PCA (prcomp via microeco)
      eig_values <- beta_obj$res_ordination$sdev^2
    }
    
    
    axis_label_x <- colnames(ordination_scores_all)[axes_choice[1]]
    axis_label_y <- colnames(ordination_scores_all)[axes_choice[2]]
    
    if(!is.null(eig_values) && length(eig_values) >= max(axes_choice)){
      explained_variance <- eig_values / sum(eig_values) * 100
      axis_label_x <- paste0(colnames(ordination_scores_all)[axes_choice[1]], " [", round(explained_variance[axes_choice[1]], 1), "%]")
      axis_label_y <- paste0(colnames(ordination_scores_all)[axes_choice[2]], " [", round(explained_variance[axes_choice[2]], 1), "%]")
    }
    
    # Define shape aesthetic: if plot_shape_var is NULL, don't map shape to a variable explicitly
    # geom_point will use a default shape.
    shape_aes_string <- if (!is.null(plot_shape_var)) plot_shape_var else shQuote("default_shape_val")
    if(is.null(plot_shape_var)) plot_data$default_shape_val <- "Sample" # dummy column for single shape
    
    beta_plot <- ggplot2::ggplot(plot_data, ggplot2::aes(x = AxisX, y = AxisY)) +
      ggplot2::geom_point(ggplot2::aes_string(color = plot_color_var, shape = shape_aes_string), size = point_size, alpha = 0.7) +
      ggplot2::xlab(axis_label_x) +
      ggplot2::ylab(axis_label_y) +
      ggplot2::labs(title = title) +
      theme_bw() +
      ggplot2::facet_wrap(stats::as.formula(paste("~", facet_var)), scales = facet_scales, nrow = facet_nrow, ncol = facet_ncol)
    
    if (any(c("ellipse") %in% plot_type) && !is.null(plot_color_var)) {
      beta_plot <- beta_plot + ggplot2::stat_ellipse(ggplot2::aes_string(fill = plot_color_var, color = plot_color_var), 
                                                     geom = "polygon", alpha = ellipse_alpha, type = "t")
    }
    if (any(c("chull") %in% plot_type) && !is.null(plot_color_var)) {
      # Ensure columns for chull are correctly referenced, matching 'AxisX' and 'AxisY'
      chull_data <- plot_data %>%
        dplyr::group_by(!!dplyr::sym(facet_var), !!dplyr::sym(plot_color_var)) %>%
        dplyr::filter(n() >= 3) %>% # chull needs at least 3 points
        dplyr::slice(chull(AxisX, AxisY)) %>% # Use generic names here
        dplyr::ungroup()
      if(nrow(chull_data) > 0){
        beta_plot <- beta_plot + ggplot2::geom_polygon(data = chull_data, 
                                                       ggplot2::aes_string(fill = plot_color_var, color = plot_color_var), 
                                                       alpha = ellipse_alpha, show.legend = FALSE)
      } else {
        message("Not enough points to draw convex hulls for some groups after faceting.")
      }
    }
    
    if (!is.null(color_values)) {
      beta_plot <- beta_plot + ggplot2::scale_color_manual(values = color_values)
      if (any(c("ellipse", "chull") %in% plot_type)) {
        beta_plot <- beta_plot + ggplot2::scale_fill_manual(values = color_values)
      }
    }
    # Only add scale_shape_manual if plot_shape_var was actually provided and used for mapping
    if (!is.null(shape_values) && !is.null(plot_shape_var)) {
      beta_plot <- beta_plot + ggplot2::scale_shape_manual(values = shape_values)
    } else if (is.null(plot_shape_var) && !is.null(shape_values) && length(shape_values)==1){
      # If no shape var, but shape_values is given (e.g. a single shape number), apply it
      beta_plot <- beta_plot + ggplot2::scale_shape_manual(values = c("Sample" = shape_values[1])) + guides(shape = "none")
    } else if (is.null(plot_shape_var)){
      beta_plot <- beta_plot + guides(shape = "none") # Hide dummy shape legend
    }
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

#' Plot Alpha Diversity by Group with Optional Faceting
#'
#' Creates boxplots of alpha diversity indices grouped by a specified variable,
#' with an option to facet the plot by another variable.
#'
#' @param dataset A microtable object. Alpha diversity should ideally be
#'   pre-calculated on this specific dataset subset. If not, it will be calculated.
#' @param group_var The grouping variable name from sample_table for comparisons on the x-axis.
#' @param measures A character vector of alpha diversity indices to plot.
#'   Currently, the function plots the first measure if multiple are provided,
#'   or you can loop outside the function to create plots for each measure.
#'   (This behavior is kept from the original for now).
#' @param facet_var Optional. A character string naming a column in sample_table
#'   to use for faceting the plot (e.g., "batch").
#' @param title Plot title (optional).
#' @param ggtheme ggplot theme to apply, e.g. theme_bw(). Default NULL uses microeco's default.
#' @param ... Additional arguments passed to `trans_alpha$plot_alpha()`.
#' @return A ggplot object.
#' @export
plot_alpha_diversity <- function(dataset,
                                 group_var,
                                 measures = c("Shannon", "Simpson", "Observed"),
                                 facet_var = NULL,
                                 title = NULL,
                                 ggtheme = NULL,
                                 ...) {
  if (!inherits(dataset, "microtable")) stop("Input 'dataset' must be a microtable object.")
  if (is.null(dataset$sample_table) || !group_var %in% names(dataset$sample_table)) {
    stop(paste0("group_var '", group_var, "' not found in dataset$sample_table."))
  }
  if (!is.null(facet_var) && !facet_var %in% names(dataset$sample_table)) {
    stop(paste0("facet_var '", facet_var, "' not found in dataset$sample_table."))
  }
  
  # Make sure alpha diversity has been calculated on the input dataset
  # It's best if this is done on the specific subset being plotted
  if (is.null(dataset$alpha_diversity)) {
    message(paste0("Alpha diversity not found in the provided dataset. Calculating for: ", paste(measures, collapse=", ")))
    dataset$cal_alphadiv(measures = measures)
    if(is.null(dataset$alpha_diversity)){
      stop("Failed to calculate alpha diversity. Ensure dataset is valid.")
    }
  }
  
  # Check if the specific measure to plot is available
  plot_measure <- measures[1] # Plotting the first measure by default
  if (!plot_measure %in% colnames(dataset$alpha_diversity)) {
    # If the specific measure is not there, try to calculate it
    message(paste0("Measure '", plot_measure, "' not found in dataset$alpha_diversity. Attempting to calculate it."))
    all_available_measures <- colnames(dataset$alpha_diversity)[!colnames(dataset$alpha_diversity) %in% rownames(dataset$sample_table)]
    dataset$cal_alphadiv(measures = unique(c(plot_measure, all_available_measures)))
    if (!plot_measure %in% colnames(dataset$alpha_diversity)) {
      stop(paste0("Measure '", plot_measure, "' could not be calculated or found."))
    }
  }
  
  
  # Create alpha_div object from the microeco package
  # 'group' in trans_alpha$new is for the primary x-axis grouping
  alpha_div_obj <- microeco::trans_alpha$new(dataset = dataset, group = group_var)
  
  # Calculate differential test (adds statistical significance markers like letters or p-values)
  # This will perform tests between levels of 'group_var'
  # If facet_var is used, these stats are for overall group_var differences, not within facets.
  # For stats within facets, one would typically filter data per facet level first, then test.
  # However, for visualization, comparing group_var across facets is common.
  
  # Only attempt cal_diff if there are at least two groups to compare in group_var
  if (length(unique(dataset$sample_table[[group_var]])) > 1) {
    alpha_div_obj$cal_diff(method = "wilcox") # Example: Wilcoxon test
  } else {
    message(paste0("Only one level found for group_var '", group_var, "'. Skipping differential test calculation for alpha diversity plot."))
  }
  
  
  # Generate boxplot using the first measure in the 'measures' vector
  alpha_p <- alpha_div_obj$plot_alpha(
    measure = plot_measure,
    plot_type = "ggboxplot", # default boxplot
    add_sig = TRUE,
    ... # Pass additional arguments like color_values, add_sig_label etc.
  )
  
  # Add title
  if (!is.null(title)) {
    alpha_p <- alpha_p + labs(title = title)
  } else {
    alpha_p <- alpha_p + labs(title = paste(plot_measure, "Diversity by", group_var))
  }
  
  # Apply custom theme if provided
  if (!is.null(ggtheme)) {
    alpha_p <- alpha_p + ggtheme
  }
  
  # Add faceting if facet_var is specified
  if (!is.null(facet_var)) {
    # Check if the facet_var actually results in multiple facets
    if (length(unique(dataset$sample_table[[facet_var]])) > 1) {
      alpha_p <- alpha_p + facet_wrap(as.formula(paste("~", facet_var)), scales = "free_y") +
        theme(strip.background = element_rect(fill="grey90"), # Example theming for facets
              strip.text = element_text(size = 10))
    } else {
      message(paste0("facet_var '", facet_var, "' has only one level in the provided dataset. Faceting will not be applied or will result in a single panel."))
    }
  }
  
  return(alpha_p)
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

#' Create a Longitudinal PCoA Plot Faceted by a Variable
#'
#' Creates a PCoA plot showing trajectories of samples across time points,
#' faceted by a specified variable (e.g., batch).
#'
#' @param dataset A microtable object.
#' @param subject_var Variable indicating subject ID (default: "pid").
#' @param time_var Variable indicating time point (default: "time").
#' @param group_var Variable for grouping/coloring subjects.
#' @param facet_var Variable in sample_table to facet by (e.g., "batch").
#' @param title Plot title (optional).
#' @param distance_metric Beta diversity metric (default: "bray").
#' @param ordination_method Ordination method (default: "PCoA").
#' @param axes_choice Numeric vector of length 2 indicating which ordination axes to plot (default: c(1, 2)).
#' @param point_size Size of the points (default: 3).
#' @param line_alpha Alpha for trajectory lines (default: 0.5).
#' @param facet_scales Scales argument for facet_wrap (default: "free").
#' @param facet_nrow Number of rows for facet_wrap (default: NULL).
#' @param facet_ncol Number of columns for facet_wrap (default: NULL).
#' @return A ggplot object with faceted longitudinal PCoA plot.
#' @export
plot_longitudinal_pcoa_faceted <- function(dataset,
                                           subject_var = "pid",
                                           time_var = "time",
                                           group_var,
                                           facet_var, # New parameter for faceting
                                           title = NULL,
                                           distance_metric = "bray",
                                           ordination_method = "PCoA",
                                           axes_choice = c(1, 2),
                                           point_size = 3,
                                           line_alpha = 0.5,
                                           facet_scales = "free",
                                           facet_nrow = NULL,
                                           facet_ncol = NULL) {
  
  if (!requireNamespace("microeco", quietly = TRUE)) stop("Please install the 'microeco' package.")
  if (!requireNamespace("ggplot2", quietly = TRUE)) stop("Please install the 'ggplot2' package.")
  if (!requireNamespace("dplyr", quietly = TRUE)) stop("Please install the 'dplyr' package.")
  
  # Ensure the dataset has the required columns
  if (!all(c(subject_var, time_var, group_var, facet_var) %in% colnames(dataset$sample_table))) {
    stop(paste("One or more specified variables (subject_var, time_var, group_var, facet_var) not found in sample_table. Available:", paste(colnames(dataset$sample_table), collapse=", ")))
  }
  
  # Calculate beta diversity if not already present on the full dataset
  # This ensures the PCoA is performed once on the entire dataset for consistent axes if desired,
  # though faceting implies the primary comparison is within facets.
  # For truly independent PCoAs per facet, one would filter then ordinate.
  # However, for faceting a global ordination, this is correct.
  if (is.null(dataset$beta_diversity[[distance_metric]])) {
    message(paste0("Calculating ", distance_metric, " beta diversity for the full dataset..."))
    dataset$cal_betadiv(method = distance_metric)
  }
  
  beta_obj <- microeco::trans_beta$new(dataset = dataset,
                                       measure = distance_metric,
                                       group = group_var) # group for potential ellipse logic if ever added
  
  beta_obj$cal_ordination(method = ordination_method, ncomp = max(axes_choice))
  
  # Extract ordination data
  scores_df_full <- as.data.frame(beta_obj$res_ordination$scores)
  if (ncol(scores_df_full) < max(axes_choice)) {
    stop(paste0("Requested axes (", paste(axes_choice, collapse=", "), ") exceed available ordination components (", ncol(scores_df_full), ")."))
  }
  scores_df <- scores_df_full[, axes_choice, drop = FALSE]
  colnames(scores_df) <- c("AxisX", "AxisY") # Generic names
  
  scores_df$SampleID <- rownames(scores_df_full) # Use original rownames
  
  # Add metadata
  metadata <- dataset$sample_table
  # Ensure SampleID is the first column or rownames for join
  if (!"SampleID" %in% colnames(metadata) && rownames(metadata)[1] != metadata[[1]][1]) {
    if (length(unique(metadata[[1]])) == nrow(metadata)) {
      rownames(metadata) <- as.character(metadata[[1]])
      message("Setting rownames of metadata from its first column for joining.")
    } else {
      stop("First column of metadata does not contain unique SampleIDs and rownames are not set as SampleIDs.")
    }
  }
  if (!"SampleID" %in% colnames(metadata)) {
    metadata$SampleID <- rownames(metadata)
  }
  
  
  # Convert SampleID columns to character for robust join
  scores_df$SampleID <- as.character(scores_df$SampleID)
  metadata$SampleID <- as.character(metadata$SampleID)
  
  plot_data <- dplyr::left_join(scores_df, metadata, by = "SampleID")
  
  if(nrow(plot_data) != nrow(scores_df)){
    warning("Not all samples in ordination scores were matched in metadata. Check SampleID consistency.")
  }
  if(any(is.na(plot_data[[facet_var]]))){
    warning(paste0("NA values found in the facet_var column '", facet_var, "'. These samples will be omitted or grouped into an NA facet."))
    plot_data <- plot_data[!is.na(plot_data[[facet_var]]),] # Optionally remove them
  }
  
  
  # Convert time_var and group_var to factor for correct plotting
  plot_data[[time_var]] <- as.factor(plot_data[[time_var]])
  plot_data[[group_var]] <- as.factor(plot_data[[group_var]])
  plot_data[[facet_var]] <- as.factor(plot_data[[facet_var]])
  
  
  # Extract variance explained for axis labels
  eig_values <- NULL
  if (!is.null(beta_obj$res_ordination$values) && !is.null(beta_obj$res_ordination$values$Eigenvalues)) {
    eig_values <- beta_obj$res_ordination$values$Eigenvalues
  } else if (!is.null(beta_obj$res_ordination$eig)) {
    eig_values <- beta_obj$res_ordination$eig
  } else if (!is.null(beta_obj$res_ordination$sdev)) {
    eig_values <- beta_obj$res_ordination$sdev^2
  }
  
  axis_label_x <- colnames(scores_df_full)[axes_choice[1]]
  axis_label_y <- colnames(scores_df_full)[axes_choice[2]]
  
  if (!is.null(eig_values) && length(eig_values) >= max(axes_choice)) {
    explained_variance <- eig_values / sum(eig_values) * 100
    axis_label_x <- paste0(colnames(scores_df_full)[axes_choice[1]], " [", round(explained_variance[axes_choice[1]], 1), "%]")
    axis_label_y <- paste0(colnames(scores_df_full)[axes_choice[2]], " [", round(explained_variance[axes_choice[2]], 1), "%]")
  }
  
  # Auto-generate title if not provided
  if (is.null(title)) {
    title <- paste("Longitudinal", toupper(ordination_method), "by", group_var, "Faceted by", facet_var)
  }
  
  # Create the plot
  p <- ggplot2::ggplot(plot_data, ggplot2::aes(x = AxisX, y = AxisY, group = !!sym(subject_var), color = !!sym(group_var))) +
    ggplot2::geom_line(alpha = line_alpha) +
    ggplot2::geom_point(aes(shape = !!sym(time_var)), size = point_size) +
    ggplot2::facet_wrap(stats::as.formula(paste("~", facet_var)), scales = facet_scales, nrow = facet_nrow, ncol = facet_ncol) +
    ggplot2::labs(
      x = axis_label_x,
      y = axis_label_y,
      title = title,
      color = group_var,
      shape = time_var
    ) +
    ggplot2::theme_bw() +
    ggplot2::theme(
      legend.position = "right",
      plot.title = ggplot2::element_text(hjust = 0.5),
      strip.background = ggplot2::element_rect(fill="grey90"),
      strip.text = ggplot2::element_text(size = 10)
    )
  
  return(p)
}
