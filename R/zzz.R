# Package internal environment for caching data
.idrr_env <- new.env(parent = emptyenv())

.onLoad <- function(libname, pkgname) {
    # Use delayedAssign to load data only when accessed
    delayedAssign("consensus_data",
        {
            # Check for RDA file first (which we know exists)
            data_path_rda <- system.file("extdata", "consensus_data.rda", package = "idrrR", lib.loc = libname)
            data_path_rds <- system.file("extdata", "consensus_data.rds", package = "idrrR", lib.loc = libname)

            cached_db <- NULL

            if (data_path_rda != "") {
                tryCatch(
                    {
                        e <- new.env()
                        load(data_path_rda, envir = e)
                        if (exists("consensus_data", envir = e)) {
                            cached_db <- e$consensus_data
                        } else {
                            # If variable name is different, take the first one
                            vars <- ls(envir = e)
                            if (length(vars) > 0) {
                                cached_db <- get(vars[1], envir = e)
                            } else {
                                stop("RDA file is empty")
                            }
                        }
                    },
                    error = function(e) {
                        stop(paste("Failed to load consensus data from RDA:", e$message))
                    }
                )
            } else if (data_path_rds != "") {
                tryCatch(
                    {
                        cached_db <- readRDS(data_path_rds)
                    },
                    error = function(e) {
                        stop(paste("Failed to load consensus data from RDS:", e$message))
                    }
                )
            } else {
                stop("Consensus data file (RDA or RDS) not found in package.")
            }

            if (!is.null(cached_db)) {
                # Convert to data.table if needed (should be already from optimize step)
                if (!requireNamespace("data.table", quietly = TRUE)) {
                    warning("data.table package is required for optimal performance.")
                } else {
                    if (!data.table::is.data.table(cached_db)) {
                        data.table::setDT(cached_db)

                        # Create keys if missing
                        norm <- function(x) gsub("-", "", toupper(x))
                        # Check checks for existing columns before assigning
                        if ("treatmentid" %in% names(cached_db)) {
                            cached_db[, clean_drug := norm(treatmentid)]
                        }
                        if ("sampleid" %in% names(cached_db)) {
                            cached_db[, clean_cell := norm(sampleid)]
                        }

                        if ("clean_drug" %in% names(cached_db) && "clean_cell" %in% names(cached_db)) {
                            data.table::setkey(cached_db, clean_drug, clean_cell)
                        }
                    }
                }
                cached_db
            } else {
                NULL
            }
        },
        assign.env = .idrr_env
    )
}
