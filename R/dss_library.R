#' Calculate DRS (formerly MSF4)
#'
#' @export
dss_library_all <- function(ic50, slope, max_response, min_conc, max_conc, residuals = NULL, observed_values = NULL,
                            ic50_abs = NULL, nls_model = NULL, min_activity = 10,
                            min_dose_range = NULL, max_dose_range = NULL) {
  # Use provided dose range if specified, otherwise use default range
  min_conc_calc <- if (!is.null(min_dose_range)) min_dose_range else min_conc
  max_conc_calc <- if (!is.null(max_dose_range)) max_dose_range else max_conc

  # Basic validation
  if (is.na(ic50) || is.na(slope) || is.na(max_response)) {
    return(data.frame(DRS = 0))
  }

  # Handle infinite/very large absIC50 values
  if (!is.null(ic50_abs) && (is.infinite(ic50_abs) || ic50_abs >= max_conc * 10)) {
    ic50_abs <- max_conc * 10 # Use double the max_conc as a proxy
  }

  # Set parameters
  max_response <- min(as.numeric(max_response), 100)[[1]]
  slope <- abs(as.numeric(slope))[[1]]
  min_response <- 0

  # Convert to log scale - use the calculation range here
  log_min_conc <- log10(min_conc_calc)
  log_max_conc <- log10(max_conc_calc)
  log_ic50 <- log10(ic50)

  # Define consistent dose-response function
  dose_response <- function(x) {
    min_response + (max_response - min_response) / (1 + 10^(slope * (log_ic50 - x)))
  }

  # Basic inverse linear weighting (original method)
  weight_1 <- function(x, center) {
    1 / (1 + 2 * abs(x - center))
  }

  # MSF4: Linear weight, IC50 center, normalized by log of max response
  # This corresponds to the original MSF1 calculation logic regarding weights and steps

  steps <- 100
  log_range <- seq(log_min_conc, log_max_conc, length.out = steps)
  max_slope_point_1 <- log_ic50 # Standard IC50 (default)

  weights <- sapply(log_range, weight_1, center = max_slope_point_1)
  weights <- weights / sum(weights) # Normalize weights

  responses <- sapply(log_range, dose_response)
  weighted_responses <- (responses - min_activity) * weights * (responses > min_activity)
  weighted_sum <- sum(weighted_responses)

  # MSF4 Calculation
  # Normalized by log10 of max response (avoid division by zero if max_response <= 1)
  norm_factor <- if (max_response > 1) log10(max_response) else 0.001
  drs_value <- weighted_sum / norm_factor

  # Validate
  if (is.na(drs_value) || !is.finite(drs_value)) drs_value <- 0
  if (drs_value < 0) drs_value <- 0
  drs_value <- round(drs_value, 4)

  return(data.frame(DRS = drs_value))
}
