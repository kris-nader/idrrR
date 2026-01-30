#' Plot Consensus Dose-Response Curve
#'
#' Creates a nature-style dose-response plot with confidence bands and data points.
#'
#' @param dose Numeric vector of doses.
#' @param response Numeric vector of responses (Inhibition %).
#' @param ic50 Numeric IC50 value (absolute).
#' @param slope Numeric slope value.
#' @param max_response Numeric max response.
#' @param dataset Vector of dataset names corresponding to each dose/response point.
#' @param point_alpha Numeric vector (or single value) for point transparency (0-1).
#' @param custom_fit List containing optional 'ic50', 'slope', 'max' for a secondary custom curve (default NULL).
#' @param title Optional title for the plot.
#'
#' @importFrom ggplot2 ggplot geom_ribbon geom_line geom_vline geom_hline geom_point scale_x_continuous scale_y_continuous ggtitle theme_minimal aes element_blank element_rect element_line element_text scale_color_gradientn scale_fill_manual scale_alpha_identity
#' @export
plot_consensus_curve <- function(dose, response, ic50, slope, max_response, dataset = NULL, point_alpha = 1, custom_fit = NULL, title = NULL) {
    if (is.null(dataset)) {
        dataset <- rep("Unknown", length(dose))
    }

    # Ensure point_alpha is same length if it's a single value
    if (length(point_alpha) == 1) {
        point_alpha <- rep(point_alpha, length(dose))
    }

    # Data frame for points
    data_points <- data.frame(
        logconc = log10(dose),
        inhibition_percent = response,
        dataset = dataset,
        alpha_val = point_alpha
    )

    # Generate fitted curve points (Consensus)
    x <- seq(log10(min(dose)), log10(max(dose)), length.out = 100)
    ic50_log <- log10(ic50)

    # 4PL Function: y = min + (max - min) / (1 + 10^(slope * (log_ic50 - x)))
    # Min fixed at 0
    get_y <- function(x_vec, ic50_val, slope_val, max_val) {
        ic50_l <- log10(ic50_val)
        0 + (max_val - 0) / (1 + 10^(slope_val * (ic50_l - x_vec)))
    }

    y_fitted <- get_y(x, ic50, slope, max_response)

    df <- data.frame(x = x, y = y_fitted)
    df <- df[order(df$x), ]

    width_factor <- 0.07
    y_rng <- max(df$y, na.rm = TRUE) - min(df$y, na.rm = TRUE)
    if (is.infinite(y_rng) || is.na(y_rng)) y_rng <- 100

    df$upper <- df$y + y_rng * width_factor
    df$lower <- df$y - y_rng * width_factor
    df$pos <- seq(0, 1, length.out = nrow(df))

    # Generate custom curve if parameters provided
    df_custom <- NULL
    if (!is.null(custom_fit)) {
        if (all(c("ic50", "slope", "max", "dose_min", "dose_max") %in% names(custom_fit))) {
            # Use custom data dose range only
            x_custom <- seq(log10(custom_fit$dose_min), log10(custom_fit$dose_max), length.out = 100)
            y_custom <- get_y(x_custom, custom_fit$ic50, custom_fit$slope, custom_fit$max)
            df_custom <- data.frame(x = x_custom, y = y_custom)
        }
    }

    # Dataset palette for consistency
    dataset_pal <- c(
        ASTRAZENICA = "#fb9a99",
        CCLE = "#e31a1c",
        CTRPv2 = "#a6cee3",
        gCSI = "#1f78b4",
        GDSC1 = "#b2df8a",
        GDSC2 = "#33a02c",
        NAIR = "#fdbf6f",
        oneil = "#ff7f00",
        PRISM = "#cab2d6",
        NCI60 = "#6a3d9a",
        FIMM = "#1f78b4",
        BeatAML = "#b15928",
        GBM2021 = "#ffff99",
        UHNBreast = "#a6cee3",
        GRAY = "#fb9a99",
        jaaks = "#cab2d6",
        Consensus = "#333333",
        User_Upload = "#333333",
        Unknown = "#999999"
    )

    p <- ggplot2::ggplot() +
        ggplot2::theme_minimal() +
        ggplot2::geom_ribbon(data = df, ggplot2::aes(x = x, ymin = lower, ymax = upper), fill = "#3A8BC9", alpha = 0.15) +
        ggplot2::geom_line(data = df, ggplot2::aes(x = x, y = upper), color = "#3A8BC9", linewidth = 0.5, alpha = 0.6) +
        ggplot2::geom_line(data = df, ggplot2::aes(x = x, y = lower), color = "#3A8BC9", linewidth = 0.5, alpha = 0.6) +
        ggplot2::geom_line(data = df, ggplot2::aes(x = x, y = y, color = pos), linewidth = 1.4) +
        ggplot2::scale_color_gradientn(colors = c("#4A90E2", "#E74C3C"), guide = "none")

    # Add Custom Fit Curve (Gray) if exists
    if (!is.null(df_custom)) {
        p <- p + ggplot2::geom_line(data = df_custom, ggplot2::aes(x = x, y = y), color = "gray50", size = 1.0, linetype = "solid")
    }

    p <- p +
        ggplot2::geom_vline(xintercept = ic50_log, color = "grey40", size = 0.6, linetype = "dashed") +
        ggplot2::geom_hline(yintercept = 50, linetype = "dotted", color = "grey60", size = 0.4) +
        ggplot2::geom_point(data = data_points, ggplot2::aes(x = logconc, y = inhibition_percent, fill = dataset, alpha = alpha_val), size = 2, shape = 21, color = "black", stroke = 0.4)

    # Handle separate legend/colors for new datasets
    # Find datasets not in palette
    used_datasets <- unique(dataset)
    unknown_datasets <- setdiff(used_datasets, names(dataset_pal))

    if (length(unknown_datasets) > 0) {
        # Assign colors for unknowns (simple hues)
        # We can use hcl checks or just a few fixed fallback colors or simple rainbow
        # For simplicity, create a function to generate distinctive colors not in palette would be overkill
        # Just use rainbow-like colors
        new_colors <- stats::setNames(rainbow(length(unknown_datasets)), unknown_datasets)
        dataset_pal <- c(dataset_pal, new_colors)
    }

    p <- p +
        ggplot2::scale_fill_manual(values = dataset_pal, name = "Dataset", drop = FALSE) +
        ggplot2::scale_alpha_identity() +
        ggplot2::scale_x_continuous(breaks = log10(c(1, 10, 100, 1000, 10000)), labels = c("1", "10", "10^2", "10^3", "10^4"), name = "Dose (nM)") +
        ggplot2::scale_y_continuous(breaks = seq(0, 100, by = 25), limits = c(-25, 125), name = "Response (%)") +
        ggplot2::ggtitle(if (is.null(title)) "Dose-Response" else title) +
        ggplot2::theme(
            panel.grid.major.y = ggplot2::element_line(color = "grey95", size = 0.3),
            panel.grid.major.x = ggplot2::element_blank(),
            plot.background = ggplot2::element_rect(fill = "white", colour = NA),
            panel.background = ggplot2::element_rect(fill = "white", colour = NA),
            panel.border = ggplot2::element_rect(color = "black", fill = NA, size = 0.5),
            plot.title = ggplot2::element_text(hjust = 0.5, size = 12, face = "bold"),
            legend.position = "right"
        )

    return(p)
}
