#' Calculate All Metrics
#'
#' @param dose Numeric vector of doses
#' @param inhibition_percent Numeric vector of inhibition
#' @param dataset Dataset name
#' @param drug_name Drug name
#' @param DSS_type DSS type (default 2)
#' @param min_dose_range Optional min dose
#' @param max_dose_range Optional max dose
#'
#' @importFrom drc drm LL.4 L.4 drmc

#' @importFrom minpack.lm nlsLM

#' @import dplyr
#' @export
calculate_all <- function(dose, inhibition_percent, dataset, drug_name, DSS_type = 2,
                          min_dose_range = NULL, max_dose_range = NULL) {
  if (all(inhibition_percent <= 0)) inhibition_percent <- rep(0, length(inhibition_percent))
  if (any(duplicated(inhibition_percent))) inhibition_percent <- seq(from = 0, length.out = length(inhibition_percent), by = 0.005) + inhibition_percent


  # DSS FUNCTION
  dss <- function(ic50, slope, max, min.conc.tested, max.conc.tested, y = 10, DSS.type = 2, concn_scale = 1e-9,
                  min_dose_range = NULL, max_dose_range = NULL) {
    # Use provided dose range if specified
    min.conc.tested <- if (!is.null(min_dose_range)) min_dose_range else min.conc.tested
    max.conc.tested <- if (!is.null(max_dose_range)) max_dose_range else max.conc.tested

    a <- as.numeric(unname(max))
    b <- as.numeric(unname(slope))
    d <- 0
    ic50 <- as.numeric(unname(ic50))
    min.conc.tested <- as.numeric(unname(min.conc.tested))
    max.conc.tested <- as.numeric(unname(max.conc.tested))
    Min.Conc <- log10(min.conc.tested * concn_scale)
    Max.Conc <- max.conc.tested
    x2 <- log10(Max.Conc * concn_scale)

    if (is.na(ic50) || is.na(b) || is.na(a) || is.na(Min.Conc) || is.na(Max.Conc)) {
      dss <- NA
    } else if (isTRUE(ic50 >= Max.Conc)) {
      dss <- 0
    } else if (isTRUE(b == 0)) {
      dss <- 0
    } else {
      if (a > 100) {
        a <- 100
      }
      if (isTRUE(b < 0)) {
        b <- -b
      }
      c <- log10(ic50 * concn_scale)
      if (a > y) {
        if (y != 0) {
          x1 <- (c - ((log(a - y) - log(y - d)) / (b * log(10))))
          if (isTRUE(x1 < Min.Conc)) {
            x1 <- Min.Conc
          } else if (isTRUE(x1 > x2)) {
            x1 <- x2
          }
        } else {
          x1 <- Min.Conc
        }
        int_y <- (((((a - d) * log(1 + 10^(b * (c - x2)))) / (b * log(10))) + a * x2) - ((((a - d) * log(1 + 10^(b * (c - x1)))) / (b * log(10))) + a * x1)) - (y * (x2 - x1))
        total_area <- (x2 - Min.Conc) * (100 - y)

        if (DSS.type == 1) {
          norm_area <- ((int_y / total_area) * 100)
        } # DSS1
        if (DSS.type == 2) {
          norm_area <- ((int_y / total_area) * 100) / log10(a) # DSS2 #AUC1
          if (isTRUE(norm_area > 50)) {
            norm_area <- 0
          }
        }
        if (DSS.type == 3) {
          norm_area <- ((int_y / total_area) * 100) * (log10(100) / log10(a)) * ((x2 - x1) / (x2 - Min.Conc))
        } # DSS3 #AUC5
        if (isTRUE(norm_area < 0 | norm_area > 100)) {
          dss <- 0
        } else {
          dss <- round(norm_area, digits = 4)
        }
      } else {
        dss <- 0
      }
    }
    return(dss)
  }


  # PLOTTING AND CALCULATIONS
  # combine the data and sort by dose.
  mat_tbl <- data.frame(inhibition_percent, dose, logconc = log10(dose), 100 - inhibition_percent, dataset = dataset)
  mat_tbl <- mat_tbl[order(mat_tbl[, "dose"]), ]

  #############    IC50
  estimate_param <- tryCatch(
    {
      drm(inhibition_percent ~ logconc, data = mat_tbl, fct = LL.4(fixed = c(NA, NA, NA, NA), names = c("SLOPE", "MIN", "MAX", "IC50")), logDose = 10, control = drmc(errorm = F))
    },
    warning = function(w) {
      drm(inhibition_percent ~ logconc, data = mat_tbl, fct = L.4(fixed = c(NA, NA, NA, NA), names = c("SLOPE", "MIN", "MAX", "IC50")), logDose = 10)
    },
    error = function(e) {
      drm(inhibition_percent ~ logconc, data = mat_tbl, fct = L.4(fixed = c(NA, NA, NA, NA), names = c("SLOPE", "MIN", "MAX", "IC50")), logDose = 10)
    }
  )
  # (extract and name coefficients)
  coef_estim <- coef(estimate_param)
  names(coef_estim) <- c("SLOPE", "MIN", "MAX", "IC50")
  # see http://www.ncbi.nlm.nih.gov/pmc/articles/PMC4696819/
  coef_estim["SLOPE"] <- coef_estim["SLOPE"] * -1

  # if curve decreases or IC50 is higher than max (i.e. IC50 is "outlier"), set IC50 to max conc.
  coef_estim["IC50"] <- ifelse(coef_estim["MAX"] <= coef_estim["MIN"] | coef_estim["IC50"] > max(mat_tbl$dose, na.rm = T), max(mat_tbl$dose, na.rm = T), coef_estim["IC50"])
  # if IC50 is less than 0 set it to min. conc. and if even min. conc. < 0, then set IC50 to mean of all conc.
  coef_estim["IC50"] <- ifelse(coef_estim["IC50"] < 0, min(mat_tbl$dose, na.rm = T), coef_estim["IC50"])
  coef_estim["IC50"] <- ifelse(coef_estim["IC50"] < 0, mean(mat_tbl$dose, na.rm = T), coef_estim["IC50"])
  # similar to previous step but now compare log10(IC50) with log(min. conc.).
  coef_estim["IC50"] <- log10(coef_estim["IC50"])
  coef_estim["IC50"] <- ifelse(coef_estim["IC50"] < min(mat_tbl$logconc), max(mat_tbl$logconc), coef_estim["IC50"])
  # if all inhib. < 0 set IC50 to max. log. conc !!!!! not obvious why!
  coef_estim["IC50"] <- ifelse(all(mat_tbl$inhibition_percent < 0), max(mat_tbl$logconc, na.rm = T), coef_estim["IC50"])
  # (Trying to fix curves that need outlier kickout)
  coef_estim["MIN"] <- 0
  coef_estim["MAX"] <- max(mat_tbl$inhibition_percent, na.rm = T)
  # (Fix off minimums) Find lowest inhibition_percent value. If it is not in (0:100), fix it whether to 0 or 99.
  min_lower <- ifelse(min(mat_tbl$inhibition_percent, na.rm = T) > 0, min(mat_tbl$inhibition_percent, na.rm = T), 0)
  min_lower <- ifelse(min_lower >= 100, 99, min_lower)
  # similar to previous step but for MAX
  coef_estim["MAX"] <- ifelse(coef_estim["MAX"] > 100, 100, coef_estim["MAX"])
  coef_estim["MAX"] <- ifelse(coef_estim["MAX"] < 0, 100, coef_estim["MAX"])
  # max_lower and max_upper - lower and upper bounds for 'nl2sol' algorithm in nonlinear least-squares
  max_lower <- ifelse(max(mat_tbl$inhibition_percent, na.rm = T) > 100, coef_estim["MAX"], max(mat_tbl$inhibition_percent, na.rm = T))
  max_lower <- ifelse(max_lower < 0, coef_estim["MAX"], max(mat_tbl$inhibition_percent, na.rm = T))
  max_lower <- ifelse(max_lower < 0, 0, max_lower)
  max_lower <- ifelse(max_lower > 100, 100, max_lower)
  # (Fix upper maximum for negative slopes)
  # Base R replacement for caTools::runmean
  runmean_base <- function(x, k) {
    # Implement a simple sliding window mean
    n <- length(x)
    if (n < k) {
      return(mean(x, na.rm = TRUE))
    }
    res <- numeric(n)

    # Fill ends with standard means or padding similar to caTools
    half_k <- floor(k / 2)
    for (i in seq_len(n)) {
      start_i <- max(1, i - half_k)
      end_i <- min(n, i + half_k)
      res[i] <- mean(x[start_i:end_i], na.rm = TRUE)
    }
    return(res)
  }
  run_avg <- runmean_base(mat_tbl$inhibition_percent, 10)
  max_upper <- ifelse(any(run_avg[-nrow(mat_tbl)] > run_avg[nrow(mat_tbl)]), max(mat_tbl$inhibition_percent[run_avg > run_avg[nrow(mat_tbl)]]), coef_estim["MAX"])
  max_upper <- ifelse(any(mat_tbl$inhibition_percent > max_upper), mean(mat_tbl$inhibition_percent[mat_tbl$inhibition_percent > max_upper]) + 5, max_upper)
  max_upper <- ifelse(max_upper < 0, coef_estim["MAX"], max_upper)
  max_upper <- ifelse(max_upper > 100, 100, max_upper) # coef_estim["MAX"]
  max_upper <- ifelse(max_lower > max_upper, coef_estim["MAX"], max_upper)
  # left it as it was, just rewritten a bit (ALEKS). not clear how values 25, 60 and 5 are chosen.
  mean_inh_last <- mean(tail(mat_tbl$inhibition_percent, 2), na.rm = T)
  if (mean_inh_last < 60) {
    if (mean_inh_last > 25) {
      coef_estim["IC50"] <- mean(mat_tbl$logconc, na.rm = T)
    } else if (mean_inh_last < 25) coef_estim["IC50"] <- max(mat_tbl$logconc, na.rm = T)
  }
  if (mean(mat_tbl$inhibition_percent[1:3], na.rm = T) < 5) coef_estim["IC50"] <- max(mat_tbl$logconc, na.rm = T)
  # add a bit of positive noise to MAX if it is the same as MIN.
  if (unname(coef_estim["MIN"]) == unname(coef_estim["MAX"])) coef_estim["MAX"] <- coef_estim["MAX"] + 0.001

  # adaptive nonlinear Least-Squares algorithm NL2SOL to estimate parameters.
  nls_result_ic50_old <- function() {
    tryCatch(
      {
        nls(inhibition_percent ~ MIN + (MAX - MIN) / (1 + (10^(SLOPE * (IC50 - logconc)))), data = mat_tbl, algorithm = "port", start = list(SLOPE = 1, MIN = coef_estim["MIN"][[1]], MAX = coef_estim["MAX"][[1]], IC50 = coef_estim["IC50"][[1]]), lower = list(SLOPE = 0.1, MIN = 0, MAX = max_lower, IC50 = min(mat_tbl$logconc)), upper = list(SLOPE = 3.5, MIN = 0, MAX = max_upper, IC50 = max(mat_tbl$logconc)), control = list(warnOnly = T, minFactor = 1 / 2048))
      },
      error = function(e) {
        # allows higher residual sum-of-squares
        minpack.lm::nlsLM(inhibition_percent ~ MIN + (MAX - MIN) / (1 + (10^(SLOPE * (IC50 - logconc)))),
          data = mat_tbl,
          start = list(SLOPE = 1, MIN = coef_estim["MIN"][[1]], MAX = coef_estim["MAX"][[1]], IC50 = coef_estim["IC50"][[1]]),
          lower = c(SLOPE = 0.1, MIN = 0, MAX = max_lower, IC50 = min(mat_tbl$logconc)),
          upper = c(SLOPE = 3.5, MIN = 0, MAX = max_upper, IC50 = max(mat_tbl$logconc))
        )
      }
    )
  }

  # IC50 first
  nls_result_ic50 <- nls_result_ic50_old()

  # IC50 second
  nls_result_ic50_2 <- tryCatch(
    {
      # allows higher residual sum-of-squares
      nls(inhibition_percent ~ MIN + (MAX - MIN) / (1 + (10^(SLOPE * (IC50 - logconc)))), data = mat_tbl, algorithm = "port", start = list(SLOPE = 1, MIN = coef_estim["MIN"][[1]], MAX = coef_estim["MAX"][[1]], IC50 = median(mat_tbl$logconc)), lower = list(SLOPE = 0.1, MIN = 0, MAX = max_lower, IC50 = min(mat_tbl$logconc)), upper = list(SLOPE = 3.5, MIN = 0, MAX = max_upper, IC50 = max(mat_tbl$logconc)), control = list(warnOnly = T, minFactor = 1 / 2048))
    },
    warning = function(w) {
      nls_result_ic50_old()
    },
    error = function(e) {
      nls_result_ic50_old()
    }
  )

  # element (4, 4) is zero, so the inverse cannot be computed
  nls_result_ic50 <- tryCatch(
    {
      summary(nls_result_ic50)
      nls_result_ic50
    },
    error = function(e) {
      nls_result_ic50_2
    }
  )

  # Calculate the standard error scores
  sumIC50 <- list(summary(nls_result_ic50), summary(nls_result_ic50_2))

  ic50std_resid <- round(sqrt(sum((sumIC50[[1]]$residuals)^2) / (length(sumIC50[[1]]$residuals) - 1)), 1)
  ic50std_resid2 <- round(sqrt(sum((sumIC50[[2]]$residuals)^2) / (length(sumIC50[[2]]$residuals) - 1)), 1)

  # continue with the best
  switch_ <- which.min(c(ic50std_resid, ic50std_resid2))
  nls_result_ic50 <- list(nls_result_ic50, nls_result_ic50_2)[[switch_]]


  # Calculate the standard error scores
  sumIC50 <- summary(nls_result_ic50)
  ic50std_Error <- sumIC50$coefficients["IC50", "Std. Error"] # tec50std_Error <- sumTEC50$coefficients["TEC50","Std. Error"]
  ic50std_resid <- round(sqrt(sum((sumIC50$residuals)^2) / (length(sumIC50$residuals) - 1)), 1)
  max_signal <- max(mat_tbl$dose, na.rm = T)
  min_signal <- min(mat_tbl$dose, na.rm = T)
  mat_tbl$residuals <- sumIC50$residuals


  se_data <- mat_tbl %>%
    group_by(dataset) %>%
    summarise(standard_error = round(sqrt(sum(residuals^2) / (n() - 1)), 1))
  # Compute some summary metric, the mean standard error
  ic50std_resid <- mean(se_data$standard_error)


  #############################
  #############   Final modification & STD error

  # prepare final data and convert IC50 back from log scale (inverse)
  coef_ic50 <- coef(nls_result_ic50)[c("IC50", "SLOPE", "MAX", "MIN")]
  coef_ic50["IC50"] <- 10^coef_ic50["IC50"]
  # (Fix ic50 for curves in wrong direction)
  coef_ic50["IC50"] <- ifelse(coef_ic50["SLOPE"] < 0, max_signal, coef_ic50["IC50"])
  # (Fix based on MAX)
  coef_ic50["IC50"] <- ifelse(coef_ic50["MAX"] < 0, max_signal, coef_ic50["IC50"])
  coef_ic50["IC50"] <- ifelse(coef_ic50["MAX"] < 10, max_signal, coef_ic50["IC50"])
  coef_ic50["MAX"] <- ifelse(coef_ic50["MAX"] < 0, 0, coef_ic50["MAX"])
  # (Fix over sensitive drugs)
  coef_ic50["IC50"] <- ifelse(all(c(max(mat_tbl$inhibition_percent, na.rm = T), min(mat_tbl$inhibition_percent, na.rm = T)) > 50), min_signal, coef_ic50["IC50"])

  # for ploting
  x_range_min <- if (!is.null(min_dose_range)) log10(min_dose_range) else min(mat_tbl$logconc)
  x_range_max <- if (!is.null(max_dose_range)) log10(max_dose_range) else max(mat_tbl$logconc)
  x <- seq(x_range_min, x_range_max, length = 100)
  yic <- predict(nls_result_ic50, data.frame(logconc = x))


  ####
  # Absolute IC50
  xIC50ABS <- seq(min(mat_tbl$logconc), max(mat_tbl$logconc) * 15, length = 5000)
  yicIC50ABS <- predict(nls_result_ic50, data.frame(logconc = xIC50ABS))
  if (all(yicIC50ABS < 50)) coef_ic50ABS <- Inf else coef_ic50ABS <- 10**xIC50ABS[which.min(abs(yicIC50ABS - 50))]
  ####


  #############################
  #############    DSS
  dss_more <- dss_library_all(as.numeric(coef_ic50["IC50"]), as.numeric(coef_ic50["SLOPE"]), as.numeric(coef_ic50["MAX"]),
    min_conc = min_signal, max_conc = max_signal, residuals = mat_tbl$residuals,
    nls_model = nls_result_ic50, observed_values = mat_tbl$inhibition_percent, ic50_abs = coef_ic50ABS, min_dose_range = min_dose_range, max_dose_range = max_dose_range
  )

  # dataframe for IC50
  IC50_dataframe <- data.frame(
    DRUG_NAME = drug_name,
    IC50ABS = coef_ic50ABS,
    IC50 = coef_ic50["IC50"],
    SLOPE = coef_ic50["SLOPE"],
    Max_Response = coef_ic50["MAX"],
    SE_of_estimate = as.numeric(ic50std_resid),
    DRS = dss_more$DRS
  )


  #  #round by 2 dex. all the numeric colums
  numeric_cols <- sapply(IC50_dataframe, is.numeric)
  IC50_dataframe[, numeric_cols] <- round(IC50_dataframe[, numeric_cols], 2)
  #
  # plot IC50
  # mat_tbl$inhibition_percent = xpr_tbl$inhibition_percent_percent[idx_filt]; # if we have all values < 0, they will be replaced
  # mat_tbl$viability = 100 - mat_tbl$inhibition_percent;  # we are replacing them back here.

  # Modified function that uses the existing fitted logistic curve values
  create_nature_style_plot <- function(x, y_fitted, ic50_value, data_points, d_name, d_score) {
    df <- data.frame(x = x, y = y_fitted)
    df <- df[order(df$x), ]
    width_factor <- 0.07
    upper <- y_fitted + (max(y_fitted) - min(y_fitted)) * width_factor
    lower <- y_fitted - (max(y_fitted) - min(y_fitted)) * width_factor
    curve_data <- data.frame(
      x = x, y = y_fitted, upper = upper, lower = lower,
      pos = seq(0, 1, length.out = length(x))
    )
    ic50_x <- log10(ic50_value)

    p <- ggplot2::ggplot() +
      ggplot2::theme_minimal() +
      ggplot2::geom_ribbon(data = curve_data, ggplot2::aes(x = x, ymin = lower, ymax = upper), fill = "#3A8BC9", alpha = 0.15) +
      ggplot2::geom_line(data = curve_data, ggplot2::aes(x = x, y = upper), color = "#3A8BC9", linewidth = 0.5, alpha = 0.6) +
      ggplot2::geom_line(data = curve_data, ggplot2::aes(x = x, y = lower), color = "#3A8BC9", linewidth = 0.5, alpha = 0.6) +
      ggplot2::geom_line(data = curve_data, ggplot2::aes(x = x, y = y, color = pos), linewidth = 1.4) +
      ggplot2::scale_color_gradientn(colors = c("#4A90E2", "#E74C3C"), guide = "none") +
      ggplot2::geom_vline(xintercept = ic50_x, color = "grey40", linewidth = 0.6, linetype = "dashed") +
      ggplot2::geom_hline(yintercept = 50, linetype = "dotted", color = "grey60", linewidth = 0.4) +
      ggplot2::geom_point(data = data_points, ggplot2::aes(x = logconc, y = inhibition_percent, fill = dataset), size = 2, shape = 21, color = "black", stroke = 0.4) +
      ggplot2::scale_fill_manual(
        values = c("#a6cee3", "#1f78b4", "#b2df8a", "#33a02c", "#fb9a99", "#e31a1c", "#fdbf6f", "#ff7f00", "#cab2d6", "#6a3d9a", "#ffff99", "#b15928"),
        name = "Dataset"
      ) +
      ggplot2::scale_x_continuous(breaks = log10(c(1, 10, 100, 1000, 10000)), labels = c("1", "10", "10²", "10³", "10⁴"), minor_breaks = NULL, name = "Dose (nM)", expand = c(0.02, 0)) +
      ggplot2::scale_y_continuous(breaks = seq(0, 100, by = 25), limits = c(-25, 125), name = "Response (%)", expand = c(0.02, 0)) +
      ggplot2::ggtitle(paste0(d_name, " (DRS: ", d_score, ")")) +
      ggplot2::theme(
        panel.grid.major.y = ggplot2::element_line(color = "grey95", linewidth = 0.3),
        panel.grid.major.x = ggplot2::element_blank(),
        panel.grid.minor = ggplot2::element_blank(),
        plot.background = ggplot2::element_rect(fill = "white", colour = NA),
        panel.background = ggplot2::element_rect(fill = "white", colour = NA),
        panel.border = ggplot2::element_rect(color = "black", fill = NA, linewidth = 0.5),
        plot.title = ggplot2::element_text(hjust = 0.5, size = 12, face = "bold"),
        axis.title = ggplot2::element_text(size = 10, face = "bold"),
        axis.title.x = ggplot2::element_text(margin = ggplot2::margin(t = 10)),
        axis.title.y = ggplot2::element_text(margin = ggplot2::margin(r = 10)),
        axis.text = ggplot2::element_text(size = 9, color = "black"),
        legend.position = c(0.9, 0.25),
        legend.background = ggplot2::element_rect(fill = NA, color = NA),
        legend.key = ggplot2::element_rect(fill = NA),
        legend.title = ggplot2::element_text(face = "bold", size = 9),
        legend.text = ggplot2::element_text(size = 8),
        legend.key.size = ggplot2::unit(0.8, "lines"),
        plot.margin = ggplot2::margin(10, 10, 10, 10)
      )
    return(p)
  }

  nature_plot <- create_nature_style_plot(
    x = x,
    y_fitted = yic,
    ic50_value = coef_ic50["IC50"],
    data_points = mat_tbl,
    d_name = drug_name,
    d_score = IC50_dataframe$DRS[1]
  )

  grDevices::pdf(paste0("./", drug_name, "_nature_style_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".pdf"), width = 3.5, height = 2.8)
  print(nature_plot)
  grDevices::dev.off()

  grDevices::png(paste0("./", drug_name, "_nature_style_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".png"), width = 5, height = 4.5, units = "in", res = 600, bg = "white")
  print(nature_plot)
  grDevices::dev.off()

  print(nature_plot)


  return(IC50_dataframe)
}
