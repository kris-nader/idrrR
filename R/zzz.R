# Package internal environment for caching data
.idrr_env <- new.env(parent = emptyenv())

# Zenodo URL
DEFAULT_ZENODO_URL <- "https://zenodo.org/records/18775863/files/iDRR_raw_data.rds?download=1"

.onLoad <- function(libname, pkgname) {
    # Do nothing on load to avoid crashing R's internal grepRaw during installation
}

get_consensus_data <- function() {
    if (is.null(.idrr_env$consensus_data)) {
        # Setup Persistent Cache Directory
        cache_dir <- tools::R_user_dir("idrrR", which = "cache")
        if (!dir.exists(cache_dir)) {
            dir.create(cache_dir, recursive = TRUE)
        }

        db_path <- file.path(cache_dir, "iDRR_raw_data.rds")

        # Download if data is missing from cache
        if (!file.exists(db_path)) {
            message("iDRR Database not found locally. Downloading from Zenodo (this may take a while)...")
            options(timeout = max(3000, getOption("timeout"))) # Prevent timeout on large file
            tryCatch(
                {
                    utils::download.file(DEFAULT_ZENODO_URL, db_path, mode = "wb", quiet = FALSE)
                    message("Download complete! Database safely cached to: ", db_path)
                },
                error = function(e) {
                    stop(paste("Failed to download database from Zenodo:", e$message))
                }
            )
        }

        tryCatch(
            {
                m <- readRDS(db_path)
                cached_db <- m
                message("Loaded consensus_data from ", db_path)

                if (!requireNamespace("data.table", quietly = TRUE)) {
                    warning("data.table package is required for optimal performance.")
                } else {
                    if (!data.table::is.data.table(cached_db)) {
                        data.table::setDT(cached_db)
                    }

                    norm <- function(x) gsub("[^[:alnum:]]", "", toupper(x))

                    # Only add clean_drug and clean_cell if they don't already exist
                    if ("treatmentid" %in% names(cached_db) && !"clean_drug" %in% names(cached_db)) {
                        cached_db[, clean_drug := norm(treatmentid)]
                    }
                    if ("sampleid" %in% names(cached_db) && !"clean_cell" %in% names(cached_db)) {
                        cached_db[, clean_cell := norm(sampleid)]
                    }

                    if ("clean_drug" %in% names(cached_db) && "clean_cell" %in% names(cached_db)) {
                        data.table::setkey(cached_db, clean_drug, clean_cell)
                    }
                }

                # Save into package local environment
                assign("consensus_data", cached_db, envir = .idrr_env)
            },
            error = function(e) {
                stop(paste("Failed to load or parse database natively:", e$message))
            }
        )
    }
    return(.idrr_env$consensus_data)
}

.onAttach <- function(libname, pkgname) {
    if (interactive()) {
        packageStartupMessage("Initializing idrrR and loading database...")
        # Eagerly load the data into the specific environment exactly once on library startup
        # Guarded by interactive() to avoid breaking R CMD INSTALL grepRaw memory limit
        get_consensus_data()
        packageStartupMessage("Database loaded successfully.")
    } else {
        packageStartupMessage("idrrR loaded in non-interactive mode. Database will be loaded on first use.")
    }
}
