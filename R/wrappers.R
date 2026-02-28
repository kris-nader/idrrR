#' Analyze User Data
#'
#' Takes a dataframe of user data, fits curves, calculates metrics, and optionally plots the result.
#'
#' @param data Data frame containing columns: Drug_Name, Cell_Line, Dose, Inhibition.
#' @param plot Logical, whether to generate plots.
#' @return A list containing metrics and plots.
#'
#' @export
analyze_own_data <- function(data, plot = TRUE) {
    # Strict validation: Require "Inhibition" column
    if (!"Inhibition" %in% colnames(data)) {
        stop("Data must contain 'Inhibition' column (Percentage 0-100 or Fraction 0-1)")
    }

    # Validation for other required columns
    req_cols <- c("Drug_Name", "Cell_Line", "Dose")
    if (!all(req_cols %in% colnames(data))) {
        stop(
            "Data must contain columns: ", paste(req_cols, collapse = ", "),
            ", and 'Inhibition'"
        )
    }

    results <- list()
    plots <- list()

    if (!"Dataset" %in% colnames(data)) {
        data$Dataset <- "User_Upload"
    } else {
        # Ensure it's character
        data$Dataset <- as.character(data$Dataset)
    }

    unique_pairs <- unique(data[, c("Drug_Name", "Cell_Line")])

    for (i in 1:nrow(unique_pairs)) {
        d <- unique_pairs$Drug_Name[i]
        c <- unique_pairs$Cell_Line[i]

        sub_data <- data[data$Drug_Name == d & data$Cell_Line == c, ]

        # Use Inhibition column
        inhibition <- sub_data$Inhibition

        # Auto-detect scale (0-1 vs 0-100)
        # If max value is <= 1, assume it's fractional (0-1) and convert to percentage
        if (max(inhibition, na.rm = TRUE) <= 1) {
            inhibition <- inhibition * 100
        }

        # Calculate Metrics
        # Use the "Overlay" function which handles both cases:
        # 1. Drug found in DB -> Overlays user data + builds combined consensus
        # 2. Drug NOT found -> Uses user data (replicates) to build consensus

        analysis_res <- tryCatch(
            {
                plot_custom_on_consensus(
                    drug_name = d,
                    cell_line = c,
                    dose = sub_data$Dose,
                    inhibition = inhibition,
                    dataset_name = sub_data$Dataset
                )
            },
            error = function(e) {
                # Fallback if something goes wrong (e.g. nls error)
                warning("Fitting failed for ", d, " x ", c, ": ", e$message)
                return(NULL)
            }
        )

        if (!is.null(analysis_res)) {
            metrics <- analysis_res$metrics_consensus
            # Add identifiers back to metrics
            metrics$Drug_Name <- d
            metrics$Cell_Line <- c

            pair_id <- paste(d, c, sep = "_")
            results[[pair_id]] <- metrics

            if (plot) {
                plots[[pair_id]] <- analysis_res$plot
            }
        }
    }

    return(list(metrics = do.call(rbind, results), plots = plots))
}

#' Get Metrics for a Drug-Cell Line Pair
#'
#' Wraps the analysis function for a single pair input.
#'
#' @export
get_metrics <- function(drug, cell_line, dose, inhibition) {
    df <- data.frame(
        Drug_Name = drug,
        Cell_Line = cell_line,
        Dose = dose,
        Inhibition = inhibition
    )
    res <- analyze_own_data(df, plot = FALSE)
    return(res$metrics)
}

#' Plot Consensus Curve from Database
#'
#' Search for a drug and cell line in the provided database and plot the consensus curve.
#' Uses optimized caching and indexing for sub-millisecond lookup.
#'
#' @param drug_name Character string, drug name to search (case-insensitive).
#' @param cell_line Character string, cell line name to search (case-insensitive).
#' @return A list containing the ggplot object, metrics, and cell line name.
#'
#' @import data.table
#' @export
plot_consensus_from_db <- function(drug_name, cell_line) {
    # 1. Access Cache (Loaded on first query by get_consensus_data in zzz.R)
    cached_db <- get_consensus_data()

    # 2. Fast Lookup
    # Standardize input
    norm <- function(x) gsub("[^[:alnum:]]", "", toupper(x))
    target_drug <- norm(drug_name)
    target_cell <- norm(cell_line)

    # Use explicit vector subsetting since keys may not be formally set to clean_drug/clean_cell
    # This maintains data.table speed while avoiding key mapping errors
    db_subset <- cached_db[grepl(target_drug, clean_drug) & clean_cell == target_cell]

    if (nrow(db_subset) == 0) {
        stop("No data found for Drug: ", drug_name, " and Cell Line: ", cell_line)
    }

    unique_dd <- unique(db_subset$DD)
    if (length(unique_dd) > 1) {
        stop("Search matched multiple unique drugs (", paste(unique_dd, collapse = ", "), "). Please refine your search query for greater specificity.")
    }

    # 3. Proceed with Calculation
    # First, summarize technical replicates by taking the median viability per dose & dataset
    dr_data <- db_subset
    if (!data.table::is.data.table(dr_data)) data.table::setDT(dr_data)

    # User requested data.table grouping specifically:
    med_data_og <- dr_data[, .(med_viability = median(viability, na.rm = TRUE)), by = .(DD, dataset, dose)]
    med_data_og_sorted <- med_data_og[order(dataset, dose)]

    b <- calculate_all(
        dose = med_data_og_sorted$dose,
        inhibition_percent = 100 - med_data_og_sorted$med_viability,
        dataset = med_data_og_sorted$dataset,
        drug_name = med_data_og_sorted$DD[1]
    )

    metrics <- b
    doses <- med_data_og_sorted$dose
    inhibition <- 100 - med_data_og_sorted$med_viability
    datasets <- med_data_og_sorted$dataset

    p <- plot_consensus_curve(
        dose = doses,
        response = inhibition,
        ic50 = metrics$IC50,
        slope = metrics$SLOPE,
        max_response = metrics$Max_Response,
        dataset = datasets,
        title = paste("Consensus:", drug_name, "x", cell_line)
    )

    return(list(
        plot = p,
        metrics = metrics,
        cell_line = cell_line,
        drug_name = drug_name,
        n_points = nrow(med_data_og_sorted) # Now returns the number of summarized points
    ))
}

#' Plot Custom Data on top of Consensus
#'
#' Overlay custom dose-response points on the existing consensus data from the database.
#' Re-calculates the consensus curve including the new data.
#'
#' @param drug_name Character string, drug name to search.
#' @param cell_line Character string, cell line name to search.
#' @param dose Numeric vector of custom doses.
#' @param inhibition Numeric vector of custom inhibition. Accepts both 0-1 and 0-100 scales (auto-detected).
#' @param dataset_name Character string for the custom dataset label.
#' @param background_alpha Numeric transparency for existing DB points (0-1). Default 0.3.
#' @return A list containing the combined plot and metrics.
#'
#' @export
plot_custom_on_consensus <- function(drug_name, cell_line, dose, inhibition, dataset_name = "User_Upload", background_alpha = 0.3) {
    # 1. Access Cache
    cached_db <- get_consensus_data()

    # 2. Get Database Data
    norm <- function(x) gsub("[^[:alnum:]]", "", toupper(x))
    target_drug <- norm(drug_name)
    target_cell <- norm(cell_line)

    db_subset <- cached_db[grepl(target_drug, clean_drug) & clean_cell == target_cell]

    unique_dd <- unique(db_subset$DD)
    if (length(unique_dd) > 1) {
        stop("Search matched multiple unique drugs (", paste(unique_dd, collapse = ", "), "). Please refine your search query for greater specificity.")
    }

    # Handle case where DB has no data -> just plot custom
    if (nrow(db_subset) == 0) {
        warning(paste("No data found in database for", drug_name, "x", cell_line, "- plotting custom data only."))
        # Create minimal DB structure for combination logic
        db_dose <- numeric(0)
        db_inhib <- numeric(0)
        db_dataset <- character(0)
        db_alpha <- numeric(0)
    } else {
        db_dose <- db_subset$dose
        db_inhib <- 100 - db_subset$viability
        db_dataset <- db_subset$dataset
        db_alpha <- rep(background_alpha, length(db_dose))
    }

    # 3. Combine Data
    # Auto-detect scale (0-1 vs 0-100)
    custom_inhib <- inhibition
    if (max(custom_inhib, na.rm = TRUE) <= 1) {
        custom_inhib <- custom_inhib * 100
    }
    custom_alpha <- rep(1, length(dose)) # Custom points fully opaque

    final_dose <- c(db_dose, dose)
    final_inhib <- c(db_inhib, custom_inhib)

    # Handle dataset_name (vector vs single)
    if (length(dataset_name) == length(dose)) {
        custom_datasets <- dataset_name
    } else {
        custom_datasets <- rep(dataset_name, length(dose))
    }
    final_dataset <- c(db_dataset, custom_datasets)
    final_alpha <- c(db_alpha, custom_alpha)

    # 4a. Calculate Combined Metrics (Consensus Fit)
    combined_df <- data.frame(dose = final_dose, inhibition = final_inhib, dataset = final_dataset, DD = drug_name)
    if (!data.table::is.data.table(combined_df)) data.table::setDT(combined_df)

    agg_combined_df <- combined_df[, .(inhibition = median(inhibition, na.rm = TRUE)), by = .(DD, dataset, dose)]

    metrics_consensus <- calculate_all(
        dose = agg_combined_df$dose,
        inhibition_percent = agg_combined_df$inhibition,
        dataset = "Combined_Consensus",
        drug_name = drug_name
    )

    # 4b. Calculate Custom-Only Metrics (Gray Curve Fit)
    custom_df <- data.frame(dose = dose, inhibition = custom_inhib, dataset = custom_datasets, DD = drug_name)
    if (!data.table::is.data.table(custom_df)) data.table::setDT(custom_df)

    agg_custom_df <- custom_df[, .(inhibition = median(inhibition, na.rm = TRUE)), by = .(DD, dataset, dose)]

    metrics_custom <- calculate_all(
        dose = agg_custom_df$dose,
        inhibition_percent = agg_custom_df$inhibition,
        dataset = custom_datasets[1], # Use a single label
        drug_name = drug_name
    )

    # 4c. Calculate Standard Error of Estimate (SEE) of Custom Points vs Consensus Curve
    # Consensus Parameters
    ic50_cons <- metrics_consensus$IC50
    slope_cons <- metrics_consensus$SLOPE
    max_cons <- metrics_consensus$Max_Response

    # Predict values for custom points using Consensus Model
    # y = min + (max - min) / (1 + 10^(slope * (log_ic50 - x)))
    predicted_y <- 0 + (max_cons - 0) / (1 + 10^(slope_cons * (log10(ic50_cons) - log10(dose))))

    residuals <- custom_inhib - predicted_y
    # SEE Formula: sqrt(sum(residuals^2) / (n - p))
    # p = 3 (Slope, Max, IC50) since Min is fixed to 0
    df_resid <- max(1, length(residuals) - 3)
    see_val <- sqrt(sum(residuals^2, na.rm = TRUE) / df_resid)

    # 5. Plot
    p <- plot_consensus_curve(
        dose = final_dose,
        response = final_inhib,
        ic50 = metrics_consensus$IC50,
        slope = metrics_consensus$SLOPE,
        max_response = metrics_consensus$Max_Response,
        dataset = final_dataset,
        point_alpha = final_alpha,
        custom_fit = list(
            ic50 = metrics_custom$IC50,
            slope = metrics_custom$SLOPE,
            max = metrics_custom$Max_Response,
            dose_min = min(dose),
            dose_max = max(dose)
        ),
        title = paste("Consensus + Custom:", drug_name, "x", cell_line)
    )

    return(list(
        plot = p,
        metrics_consensus = metrics_consensus,
        metrics_custom = metrics_custom,
        SEE_custom_vs_consensus = see_val,
        combined_data = data.frame(dose = final_dose, inhibition = final_inhib, dataset = final_dataset)
    ))
}
