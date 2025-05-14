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

#' Classify Longitudinal NM Carriage Status and Susceptibility
#'
#' Classifies subjects based on Neisseria meningitidis (NM) carriage dynamics
#' across three time points (time_0, time_2, time_4) from the 'case_control'
#' column in the dataset. It also determines a binary susceptibility group
#' (1 if ever a 'case', 0 otherwise).
#'
#' @param dataset A microtable object containing a sample_table with 'pid',
#'   'time', and 'case_control' columns.
#' @return A data frame with columns: 'pid', 'carriage_group',
#'   'susceptibility_group', 'time_0', 'time_2', 'time_4' (carriage status
#'   at each timepoint), and 'available_timepoints'.
#' @export
#' @examples
#' # Assuming 'my_microbiome_dataset' is a microtable object:
#' # carriage_df <- classify_longitudinal_carriage(my_microbiome_dataset)
classify_longitudinal_carriage <- function(dataset) {
  if (!"sample_table" %in% names(dataset)) {
    stop("Dataset must have a sample_table component.")
  }
  if (!all(c("pid", "time", "case_control") %in% names(dataset$sample_table))) {
    stop("dataset$sample_table must contain 'pid', 'time', and 'case_control' columns.")
  }
  
  # Create a wide-format table showing carriage status at each time point
  carriage_wide_df <- dataset$sample_table %>%
    dplyr::select(pid, time, case_control) %>%
    dplyr::mutate(time = as.character(time)) %>% # Ensure time is character for correct pivoting
    tidyr::pivot_wider(
      id_cols = pid,
      names_from = time,
      values_from = case_control,
      names_prefix = "time_"
    )
  
  # Ensure all expected time columns exist, fill with NA if not present for some pids
  # This handles cases where a pid might not have an entry for e.g. time_4 at all
  expected_time_cols <- c("time_0", "time_2", "time_4")
  for (col_name in expected_time_cols) {
    if (!col_name %in% names(carriage_wide_df)) {
      carriage_wide_df[[col_name]] <- NA_character_
    }
  }
  
  # Define carriage patterns
  classified_df <- carriage_wide_df %>%
    dplyr::mutate(
      available_timepoints = (!is.na(time_0)) + (!is.na(time_2)) + (!is.na(time_4)),
      
      # Susceptibility Group Logic (1 if ever 'case', 0 otherwise)
      susceptibility_group = dplyr::case_when(
        (!is.na(time_0) & time_0 == "case") |
          (!is.na(time_2) & time_2 == "case") |
          (!is.na(time_4) & time_4 == "case") ~ 1,
        TRUE ~ 0
      ),
      
      # Carriage Group Logic
      carriage_group = dplyr::case_when(
        available_timepoints <= 1 ~ "Undetermined",
        
        # Persistent Carrier: >=2 available timepoints, and all available are "case"
        available_timepoints >= 2 &
          (is.na(time_0) | time_0 == "case") &
          (is.na(time_2) | time_2 == "case") &
          (is.na(time_4) | time_4 == "case") &
          ( (time_0 == "case" & !is.na(time_0)) | # Ensure at least one actual 'case' observation
              (time_2 == "case" & !is.na(time_2)) |
              (time_4 == "case" & !is.na(time_4)) ) ~ "Persistent_carrier",
        
        # Non-Carrier: >=2 available timepoints, and all available are "control"
        available_timepoints >= 2 &
          (is.na(time_0) | time_0 == "control") &
          (is.na(time_2) | time_2 == "control") &
          (is.na(time_4) | time_4 == "control") &
          ( (time_0 == "control" & !is.na(time_0)) | # Ensure at least one actual 'control' observation
              (time_2 == "control" & !is.na(time_2)) |
              (time_4 == "control" & !is.na(time_4)) ) ~ "Non_carrier",
        
        # Acquisition Carrier
        ( (!is.na(time_0) & time_0 == "control") &
            ((!is.na(time_2) & time_2 == "case") | (!is.na(time_4) & time_4 == "case"))
        ) |
          ( is.na(time_0) & (!is.na(time_2) & time_2 == "control") & (!is.na(time_4) & time_4 == "case")
          ) ~ "Acquisition_carrier",
        
        # Clearance Carrier
        ( (!is.na(time_0) & time_0 == "case") &
            ((!is.na(time_2) & time_2 == "control") | (!is.na(time_4) & time_4 == "control"))
        ) |
          ( is.na(time_0) & (!is.na(time_2) & time_2 == "case") & (!is.na(time_4) & time_4 == "control")
          ) ~ "Clearance_carrier",
        
        TRUE ~ "Undetermined" # Fallback for any other complex patterns not explicitly covered
      )
    ) %>%
    dplyr::select(pid, carriage_group, susceptibility_group, time_0, time_2, time_4, available_timepoints)
  
  return(classified_df)
}


#' Add Carriage Columns to Dataset Sample Table
#'
#' Merges participant-level carriage group and susceptibility group classifications
#' into the sample_table of a microtable object. Ensures rownames (SampleIDs)
#' of the sample_table are preserved.
#'
#' @param dataset A microtable object.
#' @param carriage_classification_df A data frame typically generated by
#'   `classify_longitudinal_carriage`, containing 'pid', 'carriage_group',
#'   and 'susceptibility_group'.
#' @return The input microtable object with 'carriage_group' and
#'   'susceptibility_group' columns added to its sample_table.
#' @export
#' @examples
#' # Assuming 'my_dataset' is a microtable object and 'carriage_info' is from classify_longitudinal_carriage
#' # my_dataset_updated <- add_carriage_columns_to_dataset(my_dataset, carriage_info)
add_carriage_columns_to_dataset <- function(dataset, carriage_classification_df) {
  if (!"sample_table" %in% names(dataset)) {
    stop("Dataset must have a sample_table component.")
  }
  if (!"pid" %in% names(dataset$sample_table)) {
    stop("dataset$sample_table must contain a 'pid' column for joining.")
  }
  if (!all(c("pid", "carriage_group", "susceptibility_group") %in% names(carriage_classification_df))) {
    stop("carriage_classification_df must contain 'pid', 'carriage_group', and 'susceptibility_group' columns.")
  }
  
  # Ensure 'pid' types are consistent for the join
  current_sample_table <- dataset$sample_table
  current_sample_table$pid <- as.character(current_sample_table$pid)
  
  temp_carriage_lookup <- carriage_classification_df %>%
    dplyr::select(pid, carriage_group, susceptibility_group)
  temp_carriage_lookup$pid <- as.character(temp_carriage_lookup$pid)
  
  # Store original rownames (SampleIDs) of the sample_table.
  # Assumes the first column of sample_table contains the unique SampleIDs which also serve as rownames.
  # This should have been set during the initial data import (e.g., by import_qiime2_data).
  if (nrow(current_sample_table) > 0 && !is.null(rownames(current_sample_table)) &&
      all(rownames(current_sample_table) == as.character(current_sample_table[[1]]))) {
    original_sample_ids_in_order <- rownames(current_sample_table)
  } else if (nrow(current_sample_table) > 0 && !is.null(current_sample_table[[1]])) {
    original_sample_ids_in_order <- as.character(current_sample_table[[1]])
    if (length(unique(original_sample_ids_in_order)) != nrow(current_sample_table)) {
      stop("The first column of sample_table does not appear to contain unique SampleIDs. Cannot reliably preserve row order for rownames.")
    }
    if (is.null(rownames(current_sample_table)) || !all(rownames(current_sample_table) == original_sample_ids_in_order)) {
      rownames(current_sample_table) <- original_sample_ids_in_order
      # message("Rownames of sample_table were reset to values from its first column before the join.")
    }
  } else {
    stop("Cannot determine unique SampleIDs for preserving rowname order. Ensure sample_table has valid rownames or a unique SampleID in its first column.")
  }
  
  # Remove old/conflicting classification columns from sample_table before join
  cols_to_remove <- c("carriage_group", "susceptibility_group", "carriage_status")
  for (col_name in cols_to_remove) {
    if (col_name %in% names(current_sample_table)) {
      current_sample_table[[col_name]] <- NULL
    }
  }
  
  # Perform the left join
  updated_sample_table <- current_sample_table %>%
    dplyr::left_join(temp_carriage_lookup, by = "pid")
  
  # CRITICAL STEP: Restore/Set rownames using the original SampleID order
  if (nrow(updated_sample_table) == length(original_sample_ids_in_order)) {
    rownames(updated_sample_table) <- original_sample_ids_in_order
  } else {
    # This case should ideally not happen if 'pid' correctly maps classifications to samples.
    warning("Number of rows in sample_table changed unexpectedly after join. Rowname restoration might be incorrect. Further investigation needed.")
    # Attempt a fallback if the first column still represents unique SampleIDs and matches the new row count
    if (length(unique(as.character(updated_sample_table[[1]]))) == nrow(updated_sample_table)) {
      rownames(updated_sample_table) <- as.character(updated_sample_table[[1]])
      # message("Attempted to restore rownames from the first column after row number change during join.")
    }
  }
  
  dataset$sample_table <- updated_sample_table
  return(dataset)
}

#' Remove Batch 4 Samples from a Microtable Dataset
#'
#' Creates a new microtable object that specifically excludes samples
#' belonging to batch "4". It subsets the sample_table and otu_table
#' and then tidies the new dataset to ensure consistency.
#'
#' @param dataset_original A microtable object. Assumes this object has a
#'   sample_table with a column named "batch" and that the first column
#'   of sample_table contains unique SampleIDs that are also its rownames
#'   and correspond to otu_table colnames.
#' @return A new microtable object with batch "4" samples removed.
#'   Returns the original dataset with a warning if the 'batch' column is missing
#'   or if no batch "4" samples are found.
#' @export
#' @examples
#' # Assuming 'my_dataset' is a microtable object:
#' # dataset_no_b4 <- remove_batch4(my_dataset)
#' # plot_pcoa(dataset = remove_batch4(my_dataset), ...)
remove_batch4 <- function(dataset_original) {
  if (!inherits(dataset_original, "microtable")) {
    stop("Input 'dataset_original' must be a microtable object.")
  }
  if (!"batch" %in% colnames(dataset_original$sample_table)) {
    warning("Column 'batch' not found in sample_table. Returning original dataset.")
    return(dataset_original)
  }
  
  # Create a true clone to work on, so the original dataset is not modified
  dataset_filtered <- dataset_original$clone(deep = TRUE)
  
  # Ensure the 'batch' column is character for consistent filtering
  dataset_filtered$sample_table$batch <- as.character(dataset_filtered$sample_table$batch)
  
  # Identify the SampleIDs to keep (those NOT in Batch "4")
  # Assumes the first column of sample_table holds the SampleIDs and these are also the rownames.
  # This should be ensured by your import_qiime2_data function.
  
  # Defensive check for rownames based on first column
  sample_ids_from_col1 <- as.character(dataset_filtered$sample_table[[1]])
  if (is.null(rownames(dataset_filtered$sample_table)) || 
      !all(rownames(dataset_filtered$sample_table) == sample_ids_from_col1)) {
    if (length(unique(sample_ids_from_col1)) == nrow(dataset_filtered$sample_table)) {
      rownames(dataset_filtered$sample_table) <- sample_ids_from_col1
      # message("Rownames of sample_table in the clone were aligned with its first column for filtering.")
    } else {
      stop("First column of sample_table does not contain unique SampleIDs. Cannot reliably filter by batch.")
    }
  }
  
  samples_to_keep_ids <- dataset_filtered$sample_table %>%
    dplyr::filter(batch != "4") %>%
    rownames() # Pull the rownames, which are the SampleIDs
  
  if (length(samples_to_keep_ids) == nrow(dataset_filtered$sample_table)) {
    message("No samples from batch '4' found to remove. Returning a clone of the original dataset.")
    return(dataset_filtered) # Or return(dataset_original) if no clone is needed in this case
  }
  
  if (length(samples_to_keep_ids) == 0 && nrow(dataset_filtered$sample_table) > 0) {
    warning("Filtering for batch '4' resulted in zero samples to keep. This might indicate an issue or that all samples were batch 4.")
    # Return an empty-like structure or stop, based on desired behavior
    dataset_filtered$sample_table <- dataset_filtered$sample_table[0, , drop = FALSE]
    dataset_filtered$otu_table <- dataset_filtered$otu_table[, 0, drop = FALSE]
    dataset_filtered$tidy_dataset() # Tidy will handle taxa consistency
    return(dataset_filtered)
  }
  
  
  # Subset the sample_table using the identified SampleIDs (which are rownames)
  dataset_filtered$sample_table <- dataset_filtered$sample_table[samples_to_keep_ids, , drop = FALSE]
  
  # Subset the otu_table (samples are typically columns in microeco's otu_table)
  otu_samples_as_colnames <- colnames(dataset_filtered$otu_table)
  samples_in_otu_to_keep <- otu_samples_as_colnames[otu_samples_as_colnames %in% samples_to_keep_ids]
  
  if (length(samples_in_otu_to_keep) < length(samples_to_keep_ids) && length(samples_in_otu_to_keep) > 0) {
    warning("Not all selected samples (from sample_table after filtering batch 4) were found in otu_table. Subsetting otu_table with available matches.")
  } else if (length(samples_in_otu_to_keep) == 0 && length(samples_to_keep_ids) > 0) {
    warning("None of the selected samples (from sample_table after filtering batch 4) were found in otu_table. Otu_table will be empty for these samples.")
  }
  dataset_filtered$otu_table <- dataset_filtered$otu_table[, samples_in_otu_to_keep, drop = FALSE]
  
  # Tidy the dataset_filtered to ensure consistency across all components
  dataset_filtered$tidy_dataset()
  
  message(paste0("Successfully created a new dataset excluding batch '4'. Original samples: ", 
                 nrow(dataset_original$sample_table), 
                 ", Filtered samples: ", nrow(dataset_filtered$sample_table), "."))
  
  return(dataset_filtered)
}
