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
#' @importFrom caTools runmean
#' @importFrom minpack.lm nlsLM
#' @importFrom MESS auc
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
  run_avg <- caTools::runmean(mat_tbl$inhibition_percent, 10)
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
  auc <- MESS::auc(x, yic)


  # Modify the calculate_AA and calculate_AAC functions to use the custom range
  calculate_AA <- function(dose, inhibition_percent, min_dose = NULL, max_dose = NULL) {
    # Sort data by dose
    df <- data.frame(dose = dose, inhibition = inhibition_percent)
    df <- df[order(df$dose), ]

    # Filter by dose range if specified
    if (!is.null(min_dose)) {
      df <- df[df$dose >= min_dose, ]
    }
    if (!is.null(max_dose)) {
      df <- df[df$dose <= max_dose, ]
    }

    # Calculate area using rectangle method
    n <- length(df$dose)
    if (n <= 1) {
      return(0)
    }

    # Convert doses to log scale (as shown in the paper)
    log_doses <- log10(df$dose)

    # Calculate width of each rectangle in log space
    widths <- diff(log_doses)

    # Calculate sum of areas (rectangle method as described in the paper)
    aa <- sum(df$inhibition[-1] * widths) # Using inhibition values as heights

    return(aa)
  }
  AA <- calculate_AA(mat_tbl$dose, mat_tbl$inhibition_percent, min_dose = min_dose_range, max_dose = max_dose_range)


  # Similarly update the AAC calculation
  calculate_AAC <- function(concentrations, responses, min_dose = NULL, max_dose = NULL) {
    # Filter by dose range if specified
    if (!is.null(min_dose) || !is.null(max_dose)) {
      indices <- rep(TRUE, length(concentrations))
      if (!is.null(min_dose)) {
        indices <- indices & (concentrations >= min_dose)
      }
      if (!is.null(max_dose)) {
        indices <- indices & (concentrations <= max_dose)
      }
      concentrations <- concentrations[indices]
      responses <- responses[indices]
    }

    # Integrate response curve over log concentration
    log_conc <- log10(concentrations)
    # Area calculation with trapezoidal rule
    areas <- 0
    for (i in 2:length(log_conc)) {
      areas <- areas + (log_conc[i] - log_conc[i - 1]) *
        (responses[i] + responses[i - 1]) / 2
    }
    return(areas / (log_conc[length(log_conc)] - log_conc[1]))
  }
  AAC <- calculate_AAC(mat_tbl$dose, mat_tbl$inhibition_percent, min_dose = min_dose_range, max_dose = max_dose_range)
  if (length(AAC) == 0) {
    AAC <- 0
  }


  # OAUC - Outlier-Adjusted AUC
  OAUC <- function(auc, residuals, threshold = 2) {
    # Penalize AUC based on presence of outliers
    outlier_penalty <- sum(abs(residuals) > (threshold * sd(residuals))) / length(residuals)
    return(auc * (1 - outlier_penalty))
  }
  auc2 <- OAUC(auc, mat_tbl$residuals)

  ## average replicates
  mat_tblCp <- mat_tbl[, c("inhibition_percent", "dose")]
  cols_ <- colnames(mat_tblCp)[!grepl("inhibition_percent", colnames(mat_tblCp))] # columns which should be equal to average PI
  X <- as.data.table(mat_tblCp)
  mat_tblCp <- as.data.frame(X[, list(inhibition_percent = mean(inhibition_percent)), cols_], stringAsFactors = !1)


  perInh <- t(matrix(mat_tblCp[, "inhibition_percent"],
    dimnames =
      list(paste0(rep("D", length(mat_tblCp[, "inhibition_percent"])), 1:length(mat_tblCp[, "inhibition_percent"])))
  ))

  coef_tec50 <- coef_ic50
  coef_tec50["IC50"] <- ifelse(coef_tec50["MAX"] > 25, coef_tec50["IC50"], max(mat_tbl$dose, na.rm = T))
  names(coef_tec50) <- c("EC50", "SLOPE", "MAX", "MIN")
  coef_tec50["SLOPE"] <- -1 * coef_tec50["SLOPE"] # min - 0, max - 77 in ec50 it is max - 100, min - 23
  tmp <- coef_tec50["MAX"]
  coef_tec50["MAX"] <- 100 - coef_tec50["MIN"]
  coef_tec50["MIN"] <- 100 - tmp
  ytec <- 100 - yic
  perViaTox <- 100 - perInh


  ####
  # Absolute IC50
  xIC50ABS <- seq(min(mat_tbl$logconc), max(mat_tbl$logconc) * 15, length = 5000)
  yicIC50ABS <- predict(nls_result_ic50, data.frame(logconc = xIC50ABS))
  if (all(yicIC50ABS < 50)) coef_ic50ABS <- Inf else coef_ic50ABS <- 10**xIC50ABS[which.min(abs(yicIC50ABS - 50))]
  ####


  # Absolute IC20
  if (all(yicIC50ABS < 10)) coef_ic10ABS <- Inf else coef_ic10ABS <- 10**xIC50ABS[which.min(abs(yicIC50ABS - 10))]
  ####

  # Absolute IC20
  if (all(yicIC50ABS < 20)) coef_ic20ABS <- Inf else coef_ic20ABS <- 10**xIC50ABS[which.min(abs(yicIC50ABS - 20))]
  ####

  ##  Decay parameter - related to how quickly response decays with concentration
  Decay <- (coef_ic50["SLOPE"] * (coef_ic50["MAX"] / 100))[[1]]


  # Potency Interval (PI)
  # Range between ic20 and IC10 - shows steepness of dose-response
  PI <- function(ic20, ic10) {
    return(log10(ic20 / ic10))
  }

  # Measures inequality in drug response across cell lines
  Gini <- function(responses) {
    responses <- sort(responses)
    n <- length(responses)
    index <- 1:n
    return((2 * sum(index * responses) / (n * sum(responses))) - (n + 1) / n)
  }

  activity_at_conc <- function(nls_model, conc) {
    predict(nls_model, data.frame(logconc = log10(conc)))
  }

  ACT1 <- activity_at_conc(nls_result_ic50, 1)
  ACT10 <- activity_at_conc(nls_result_ic50, 10)
  ACT100 <- activity_at_conc(nls_result_ic50, 100)

  hill_coefficient <- abs(coef_ic50["IC50"]) * log(10)

  #############################
  #############    DSS

  # Calculate specific DSS variants requested by user
  dss_score1 <- round(as.numeric(dss(coef_ic50["IC50"], coef_ic50["SLOPE"], coef_ic50["MAX"], min_signal, max_signal, DSS.type = as.integer(1), y = 10, min_dose_range = min_dose_range, max_dose_range = max_dose_range)), 2)
  dss_score2 <- round(as.numeric(dss(coef_ic50["IC50"], coef_ic50["SLOPE"], coef_ic50["MAX"], min_signal, max_signal, DSS.type = as.integer(2), y = 10, min_dose_range = min_dose_range, max_dose_range = max_dose_range)), 2)
  dss_score3 <- round(as.numeric(dss(coef_ic50["IC50"], coef_ic50["SLOPE"], coef_ic50["MAX"], min_signal, max_signal, DSS.type = as.integer(3), y = 10, min_dose_range = min_dose_range, max_dose_range = max_dose_range)), 2)
  dss_score4 <- round(as.numeric(dss(coef_ic50["IC50"], coef_ic50["SLOPE"], coef_ic50["MAX"], min_signal, max_signal, DSS.type = as.integer(1), y = 5, min_dose_range = min_dose_range, max_dose_range = max_dose_range)), 2)
  dss_score5 <- round(as.numeric(dss(coef_ic50["IC50"], coef_ic50["SLOPE"], coef_ic50["MAX"], min_signal, max_signal, DSS.type = as.integer(2), y = 5, min_dose_range = min_dose_range, max_dose_range = max_dose_range)), 2)
  dss_score6 <- round(as.numeric(dss(coef_ic50["IC50"], coef_ic50["SLOPE"], coef_ic50["MAX"], min_signal, max_signal, DSS.type = as.integer(3), y = 5, min_dose_range = min_dose_range, max_dose_range = max_dose_range)), 2)
  dss_score10 <- round(as.numeric(dss(coef_ic50["IC50"], coef_ic50["SLOPE"], coef_ic50["MAX"], min_signal, max_signal, DSS.type = as.integer(1), y = 15, min_dose_range = min_dose_range, max_dose_range = max_dose_range)), 2)
  dss_score11 <- round(as.numeric(dss(coef_ic50["IC50"], coef_ic50["SLOPE"], coef_ic50["MAX"], min_signal, max_signal, DSS.type = as.integer(2), y = 15, min_dose_range = min_dose_range, max_dose_range = max_dose_range)), 2)
  dss_score12 <- round(as.numeric(dss(coef_ic50["IC50"], coef_ic50["SLOPE"], coef_ic50["MAX"], min_signal, max_signal, DSS.type = as.integer(3), y = 15, min_dose_range = min_dose_range, max_dose_range = max_dose_range)), 2)

  dss_scores_df <- data.frame(
    DSS_Type1_y10 = dss_score1, DSS_Type2_y10 = dss_score2, DSS_Type3_y10 = dss_score3,
    DSS_Type1_y5 = dss_score4, DSS_Type2_y5 = dss_score5, DSS_Type3_y5 = dss_score6,
    DSS_Type1_y15 = dss_score10, DSS_Type2_y15 = dss_score11, DSS_Type3_y15 = dss_score12
  )

  dss_more <- dss_library_all(as.numeric(coef_ic50["IC50"]), as.numeric(coef_ic50["SLOPE"]), as.numeric(coef_ic50["MAX"]),
    min_conc = min_signal, max_conc = max_signal, residuals = mat_tbl$residuals,
    nls_model = nls_result_ic50, observed_values = mat_tbl$inhibition_percent, ic50_abs = coef_ic50ABS, min_dose_range = min_dose_range, max_dose_range = max_dose_range
  )


  # dataframe for IC50
  IC50_dataframe <- data.frame(
    DRUG_NAME = drug_name, ANALYSIS_NAME = "IC50", IC50ABS = coef_ic50ABS, t(as.matrix(coef_ic50)),
    SE_of_estimate = as.numeric(ic50std_resid), AUC = auc, IC10ABS = coef_ic10ABS, IC20ABS = coef_ic20ABS
  )

  # Include MSF4 and DSS variants
  IC50_dataframe$MSF4 <- dss_more$MSF4
  IC50_dataframe <- cbind(IC50_dataframe, dss_scores_df)

  IC50_dataframe$ACT1 <- ACT1
  IC50_dataframe$ACT10 <- ACT10
  IC50_dataframe$ACT100 <- ACT100
  IC50_dataframe$AA <- AA
  IC50_dataframe$AAC <- AAC
  IC50_dataframe$OAUC <- auc2
  IC50_dataframe$Decay <- Decay
  IC50_dataframe$pIC50 <- -log10(IC50_dataframe$IC50)
  IC50_dataframe$pIC50ABS <- -log10(IC50_dataframe$IC50ABS)
  IC50_dataframe$pIC20ABS <- -log10(IC50_dataframe$IC20ABS)
  IC50_dataframe$pIC10ABS <- -log10(IC50_dataframe$IC10ABS)
  IC50_dataframe$PI <- PI(IC50_dataframe$pIC20ABS, IC50_dataframe$pIC10ABS)
  IC50_dataframe$Gini <- Gini(mat_tbl$inhibition_percent)
  IC50_dataframe$Hill <- hill_coefficient

  # Compatibility columns for wrappers.R
  IC50_dataframe$IC50_absolute <- IC50_dataframe$IC50
  IC50_dataframe$Slope <- IC50_dataframe$SLOPE
  IC50_dataframe$Max_Response <- IC50_dataframe$MAX


  #  #round by 2 dex. all the numeric colums
  numeric_cols <- sapply(IC50_dataframe, is.numeric)
  IC50_dataframe[, numeric_cols] <- round(IC50_dataframe[, numeric_cols], 2)
  #
  #  # plot IC50
  #  #mat_tbl$inhibition_percent = xpr_tbl$inhibition_percent_percent[idx_filt]; # if we have all values < 0, they will be replaced
  #  #mat_tbl$viability = 100 - mat_tbl$inhibition_percent;  # we are replacing them back here.


  return(IC50_dataframe)
}
