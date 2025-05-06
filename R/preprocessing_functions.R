#' Preprocessing Functions for Microbiome Analysis
#' This file contains functions for data preprocessing, batch correction,
#' and data preparation for NM microbiome analysis
#' 
#' @author Orr Tobaly & Nikol Elyashov
#' @date April 2025
#' 
#' dmckdfnkdfnmdf,dchanges

# Required packages --------------------------------------------------------
# library(vegan)
# library(microbiome)
# library(phyloseq)
# library(qiime2R)
# library(microeco)
# library(dplyr)
# library(tidyverse)
# library(magrittr)
# library(sva)
# library(GUniFrac)

#' Import Raw Data from QIIME2 and Create Dataset
#'
#' Imports metadata, taxonomy, feature table and phylogenetic tree from QIIME2
#' outputs and creates a microeco dataset object
#'
#' @param metadata_path Path to the QIIME2 metadata file
#' @param feature_table_path Path to the QIIME2 feature table
#' @param taxonomy_path Path to the QIIME2 taxonomy file
#' @param tree_path Path to the QIIME2 phylogenetic tree
#' @return A microtable object with the imported data
#' @export
import_qiime2_data <- function(metadata_path, feature_table_path, taxonomy_path, tree_path) {
  # Import data
  metadata <- read_q2metadata(metadata_path)
#  metadata <- metadata %>% filter(time == 2) # time filter
  taxa <- read_qza(taxonomy_path)$data %>% parse_taxonomy()
  feature_table <- as.data.frame(read_qza(feature_table_path)$data)
  tree <- read_qza(tree_path)$data
  demodata <- as.data.frame(read_excel(file.path(raw_data_dir, "demographic data _Jan2022_final.YM.20241230.xlsx"))) %>% 
    select(c('pid','cohort'))
  demodata$pid <- as.factor(demodata$pid)
  demodata$cohort <- as.factor(demodata$cohort)
  metadata <- left_join(metadata,demodata,by = 'pid')
  
  # Set sample IDs as rownames for metadata
  rownames(metadata) <- metadata[,1]
  
  # Create microeco dataset
  dataset <- microtable$new(sample_table = metadata, 
                            otu_table = feature_table, 
                            tax_table = taxa,
                            phylo_tree = tree)
  
  # Basic data tidying
  dataset$tidy_dataset()
  dataset$tax_table %<>% base::subset(Kingdom == "Archaea" | Kingdom == "Bacteria")
  dataset$filter_pollution(taxa = c("mitochondria", "chloroplast"))
  
  return(dataset)
}

#' Create Enhanced Time Variables
#'
#' Creates meaningful time variables based on cohort and time offset information
#' using vectorized operations and working directly with year-month combinations
#'
#' @param dataset A microtable object
#' @param cohort_var The name of the cohort variable (default: "cohort")
#' @param time_var The name of the time offset variable (default: "time")
#' @return The dataset with enhanced time variables added
#' @export
create_enhanced_time_variables <- function(dataset, cohort_var = "cohort", time_var = "time") {
  # Clone the dataset to avoid modifying the original
  dataset_enhanced <- dataset$clone()
  
  # Define cohort base months and years
  cohort_base_months <- c("1" = 3, "2" = 8, "3" = 11)  # Mar, Aug, Nov
  cohort_base_year <- 2019
  
  # Store current locale
  current_locale <- Sys.getlocale("LC_TIME")
  
  # Temporarily set locale to English
  Sys.setlocale("LC_TIME", "C")
  
  # Add enhanced time variables using vectorized operations
  dataset_enhanced$sample_table <- dataset_enhanced$sample_table %>%
    mutate(
      # Create base month and year based on cohort
      base_month = unname(cohort_base_months[as.character(!!sym(cohort_var))]),
      
      # Calculate actual collection month (base month + time offset)
      collection_month = base_month + as.numeric(as.character(!!sym(time_var))),
      collection_year = cohort_base_year + floor((collection_month - 1) / 12),  # Adjust year if month > 12
      collection_month = ((collection_month - 1) %% 12) + 1,  # Adjust month if > 12
      
      # Create full date object (first day of the month)
      collection_date_full = as.Date(paste0(
        collection_year, "-", 
        formatC(collection_month, width = 2, flag = "0"), 
        "-01"
      )),
      
      # Format as MMM-YY string for display - using English month abbreviations
      collection_date_str = paste0(
        month.abb[collection_month], "-", 
        substr(as.character(collection_year), 3, 4)
      ),
      
      # Calculate months since study start (March 2019)
      months_since_start = (collection_year - cohort_base_year) * 12 + 
        (collection_month - cohort_base_months["1"]),
      
      # Determine season
      season = factor(case_when(
        collection_month %in% c(3, 4, 5) ~ "Spring",
        collection_month %in% c(6, 7, 8) ~ "Summer",
        collection_month %in% c(9, 10, 11) ~ "Fall",
        collection_month %in% c(12, 1, 2) ~ "Winter"
      ), levels = c("Winter", "Spring", "Summer", "Fall"))
    )
  
  # Get unique dates and sort them chronologically
  unique_dates <- dataset_enhanced$sample_table %>%
    select(collection_date_str, collection_date_full) %>%
    distinct() %>%
    arrange(collection_date_full) %>%
    pull(collection_date_str)
  
  # Create the ordered factor for collection_date
  dataset_enhanced$sample_table <- dataset_enhanced$sample_table %>%
    mutate(
      collection_date = factor(
        collection_date_str, 
        levels = unique_dates, 
        ordered = TRUE
      )
    ) %>%
    # Remove temporary columns
    select(-base_month, -collection_month, -collection_year, -collection_date_str)
  
  # Restore original locale
  Sys.setlocale("LC_TIME", current_locale)
  
  return(dataset_enhanced)
}

#' Rarefy Samples to Equal Depth
#'
#' Rarefies all samples to the same sequencing depth
#'
#' @param dataset A microtable object
#' @param sample_size Target sequencing depth (default: 5000)
#' @param method Rarefaction method (default: "rarefy")
#' @return The dataset with rarefied counts
#' @export
rarefy_dataset <- function(dataset, sample_size = 5000, method = "rarefy") {
  dataset$rarefy_samples(sample.size = sample_size, method = method)
  return(dataset)
}

#' Create Beta Diversity Ordination Plot
#'
#' Creates a beta diversity ordination plot using a specified distance metric and method
#'
#' @param dataset A microtable object
#' @param factor_name Factor to color/group by in the plot
#' @param measure Beta diversity measure (default: "bray")
#' @param method Ordination method (default: "PCoA")
#' @param plot_type Plot type (default: c("point", "ellipse"))
#' @param title Plot title (optional)
#' @param color_values Color values for plotting (optional)
#' @return A list containing the beta diversity object and the plot
#' @export
create_beta_ordination <- function(dataset, 
                                   factor_name, 
                                   measure = "bray", 
                                   method = "PCoA",
                                   plot_type = c("point", "ellipse"),
                                   title = NULL,
                                   color_values = NULL) {
  
  # Create the trans_beta object
  beta_obj <- trans_beta$new(dataset = dataset, 
                             group = factor_name, 
                             measure = measure)
  
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
    plot_shape = factor_name,
  ) + 
    labs(title = title)
  
  # Return both the object and the plot
  return(list(
    beta_object = beta_obj,
    plot = beta_plot
  ))
}

#' Batch Effect Correction Using ComBat
#'
#' Applies ComBat batch correction to microbiome data to remove batch effects.
#' Offers two approaches: a simple log transformation (recommended) or the CLR approach.
#'
#' @param dataset A microtable object containing microbiome data
#' @param batch_var The variable in sample_table indicating batch (e.g., "year", "batch")
#' @param method Which method to use for transformation: "simple_log" (recommended) or "clr"
#' @param round_counts Whether to round the corrected values to integers (default: TRUE)
#' @param pseudocount Value to add to zeros before transformation (default: 0.5)
#' @return A new dataset with batch-corrected abundances
#' @export
correct_batch_effect <- function(dataset, batch_var = "year", 
                                 method = c("simple_log", "clr"),
                                 round_counts = TRUE,
                                 pseudocount = 0.5) {
  # Match method argument
  method <- match.arg(method)
  
  # Extract the abundance matrix
  abundance_matrix <- dataset$otu_table
  
  # Get the batch information
  batch <- dataset$sample_table[[batch_var]]
  
  # Make sure the data is in the right format for ComBat
  # ComBat expects samples as columns, features as rows
  if(nrow(abundance_matrix) > ncol(abundance_matrix)) {
    # If we have more rows than columns, the data is likely already in the right format
    # (features as rows, samples as columns)
    abundance_matrix_for_combat <- abundance_matrix
  } else {
    # If we have more columns than rows, transpose
    abundance_matrix_for_combat <- t(abundance_matrix)
  }
  
  # Add pseudocount to handle zeros
  abundance_matrix_for_combat[abundance_matrix_for_combat == 0] <- pseudocount
  
  # Apply transformation based on selected method
  if (method == "simple_log") {
    # Simple log2 transformation (recommended for ComBat)
    transformed_abundance <- log2(abundance_matrix_for_combat)
  } else if (method == "clr") {
    # CLR transformation - compute geometric mean for each sample
    geo_means <- apply(abundance_matrix_for_combat, 1, function(x) exp(mean(log(x))))
    transformed_abundance <- log2(sweep(abundance_matrix_for_combat, 1, geo_means, "/"))
  }
  
  # Run ComBat - ensure batch is a factor
  batch_factor <- as.factor(batch)
  corrected_abundance <- ComBat(dat = transformed_abundance, 
                                batch = batch_factor, 
                                mod = NULL,
                                par.prior = TRUE)
  
  # Convert back from log2
  corrected_counts <- 2^corrected_abundance
  
  # Handle any values below pseudocount
  corrected_counts[corrected_counts < pseudocount] <- 0
  
  # Round to integers if requested (makes sense for count data)
  if (round_counts) {
    corrected_counts <- round(corrected_counts)
  }
  
  # Ensure the corrected data is in the right format for the dataset object
  if(nrow(abundance_matrix) > ncol(abundance_matrix)) {
    # If the original data had features as rows, no need to transpose back
    final_corrected_counts <- corrected_counts
  } else {
    # If we transposed earlier, transpose back
    final_corrected_counts <- t(corrected_counts)
  }
  
  # Create a new dataset with corrected abundances
  dataset_corrected <- dataset$clone()
  dataset_corrected$otu_table <- final_corrected_counts
  
  return(dataset_corrected)
}

#' Classify Longitudinal Carriage Status
#'
#' Classifies subjects based on NM carriage dynamics across time points
#'
#' @param dataset A microtable object with longitudinal microbiome data
#' @return A dataframe with carriage status classifications for each subject
#' @export
classify_carriage_status <- function(dataset) {
  # Create a wide-format table showing carriage status at each time point
  carriage_patterns <- dataset$sample_table %>%
    select(pid, time, case_control) %>%
    # Convert to wide format to see patterns across time
    tidyr::pivot_wider(
      id_cols = pid,
      names_from = time,
      values_from = case_control,
      names_prefix = "time_"
    )
  
  # Classify each subject based on their carriage pattern with missing data awareness
  carriage_patterns <- carriage_patterns %>%
    mutate(
      available_timepoints = (!is.na(time_0)) + (!is.na(time_2)) + (!is.na(time_4)),
      carriage_status = case_when(
        # Subjects with all three time points
        !is.na(time_0) & !is.na(time_2) & !is.na(time_4) & 
          time_0 == "case" & time_2 == "case" & time_4 == "case" ~ "Persistent_carrier",
        
        !is.na(time_0) & !is.na(time_2) & !is.na(time_4) & 
          time_0 == "control" & (time_2 == "case" | time_4 == "case") ~ "Acquisition_carrier",
        
        !is.na(time_0) & !is.na(time_2) & !is.na(time_4) & 
          time_0 == "case" & (time_2 == "control" | time_4 == "control") ~ "Clearance_carrier",
        
        !is.na(time_0) & !is.na(time_2) & !is.na(time_4) & 
          time_0 == "control" & time_2 == "control" & time_4 == "control" ~ "Non_carrier",
        
        # Most common pattern: time_0 and time_2 available, time_4 missing
        !is.na(time_0) & !is.na(time_2) & is.na(time_4) & 
          time_0 == "case" & time_2 == "case" ~ "Persistent_carrier_incomplete",
        
        !is.na(time_0) & !is.na(time_2) & is.na(time_4) & 
          time_0 == "control" & time_2 == "case" ~ "Acquisition_carrier_incomplete",
        
        !is.na(time_0) & !is.na(time_2) & is.na(time_4) & 
          time_0 == "case" & time_2 == "control" ~ "Clearance_carrier_incomplete",
        
        !is.na(time_0) & !is.na(time_2) & is.na(time_4) & 
          time_0 == "control" & time_2 == "control" ~ "Non_carrier_incomplete",
        
        # Pattern with missing time_2
        !is.na(time_0) & is.na(time_2) & !is.na(time_4) & 
          time_0 == "case" & time_4 == "case" ~ "Likely_persistent_carrier",
        
        !is.na(time_0) & is.na(time_2) & !is.na(time_4) & 
          time_0 == "control" & time_4 == "case" ~ "Acquisition_carrier",
        
        !is.na(time_0) & is.na(time_2) & !is.na(time_4) & 
          time_0 == "case" & time_4 == "control" ~ "Clearance_carrier",
        
        !is.na(time_0) & is.na(time_2) & !is.na(time_4) & 
          time_0 == "control" & time_4 == "control" ~ "Non_carrier",
        
        # Missing baseline but have other timepoints
        is.na(time_0) & !is.na(time_2) & (!is.na(time_4)) ~ "Insufficient_baseline_data",
        
        # Insufficient data (only one timepoint)
        available_timepoints == 1 ~ "Insufficient_data",
        
        # Catch-all for any other patterns
        TRUE ~ "Undetermined"
      )
    )
  
  # Create a simplified categorization for analysis
  carriage_patterns <- carriage_patterns %>%
    mutate(
      carriage_group = case_when(
        grepl("Persistent", carriage_status) ~ "Persistent_carrier",
        grepl("Acquisition", carriage_status) ~ "Acquisition_carrier",
        grepl("Clearance", carriage_status) ~ "Clearance_carrier",
        grepl("Non_carrier", carriage_status) ~ "Non_carrier",
        TRUE ~ "Insufficient_data"
      )
    )
  
  return(carriage_patterns)
}

#' Create Expanded Susceptibility Group
#'
#' Creates an expanded group definition for susceptibility analysis
#'
#' @param carriage_patterns The carriage patterns dataframe from classify_carriage_status
#' @return A dataframe with expanded susceptibility grouping
#' @export
create_susceptibility_groups <- function(carriage_patterns) {
  # Define expanded acquisition status
  expanded_carriage_patterns <- carriage_patterns %>%
    mutate(
      shows_acquisition = case_when(
        # Traditional acquisition (control → case)
        !is.na(time_0) & time_0 == "control" & 
          ((!is.na(time_2) & time_2 == "case") | (!is.na(time_4) & time_4 == "case")) ~ TRUE,
        
        # Reacquisition (case → control → case)
        !is.na(time_0) & !is.na(time_2) & !is.na(time_4) & 
          time_0 == "case" & time_2 == "control" & time_4 == "case" ~ TRUE,
        
        # Alternatively, capture pattern with missing middle time point
        !is.na(time_0) & is.na(time_2) & !is.na(time_4) & 
          time_0 == "case" & time_4 == "case" & grepl("Clearance", carriage_status) ~ TRUE,
        
        # Default
        TRUE ~ FALSE
      )
    )
  
  return(expanded_carriage_patterns)
}

#' Update Dataset with Carriage Classifications
#'
#' Updates the sample_table in the dataset with carriage classifications
#'
#' @param dataset A microtable object
#' @param carriage_patterns The carriage patterns dataframe from classify_carriage_status
#' @param expanded_patterns The expanded carriage patterns from create_susceptibility_groups (optional)
#' @return Updated dataset with carriage classifications
#' @export
update_dataset_with_carriage_status <- function(dataset, carriage_patterns, expanded_patterns = NULL) {
  # Create a lookup table with pid and carriage_status/group
  carriage_lookup <- carriage_patterns %>%
    select(pid, carriage_status, carriage_group)
  
  # Join this information back to the sample data
  dataset$sample_table <- dataset$sample_table %>%
    left_join(carriage_lookup, by = "pid")
  
  # Add expanded susceptibility grouping if provided
  if(!is.null(expanded_patterns)) {
    # Get list of PIDs showing acquisition
    acquisition_pids <- expanded_patterns %>%
      filter(shows_acquisition == TRUE) %>%
      pull(pid)
    
    # Create a new column in the sample table for the expanded grouping
    dataset$sample_table <- dataset$sample_table %>%
      mutate(
        susceptibility_group = case_when(
          pid %in% acquisition_pids ~ "Ever_acquired_NM",
          carriage_group == "Non_carrier" ~ "Never_acquired_NM",
          TRUE ~ "Other"  # All other patterns
        )
      )
  }
  
  return(dataset)
}

#' Convert Microeco Dataset to Phyloseq
#'
#' Converts a microeco dataset to a phyloseq object for compatibility with other packages
#'
#' @param dataset A microtable object
#' @return A phyloseq object
#' @export
convert_to_phyloseq <- function(dataset) {
  otu_table <- dataset$otu_table
  tax_table <- as.matrix(dataset$tax_table)
  rownames(dataset$sample_table) <- dataset$sample_table[,1]
  sample_data <- dataset$sample_table
  
  ps <- phyloseq(
    otu_table(otu_table, taxa_are_rows = TRUE),
    tax_table(tax_table),
    sample_data(sample_data)
  )
  
  return(ps)
}