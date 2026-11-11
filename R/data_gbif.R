# prepare_gbif_data -------

#' @title Download and Process GBIF Occurrence Data
#'
#' @description Downloads, combines, cleans, and saves occurrence records from the
#'   Global Biodiversity Information Facility ([GBIF](https://gbif.org/)) for
#'   the supplied taxon keys.
#'
#' @param gbif_ids (character) Vector of one or more GBIF (COL XR) taxon keys.
#'   Data for multiple keys are combined. If `NULL` (default),
#'   keys are retrieved from the `onesdm_gbif_ids` option. \strong{Required}.
#' @param model_dir (character) Directory for saving model outputs.
#'   A `data` subdirectory is created automatically to store raw and processed data.
#'   If `NULL` (default), path is retrieved from the `onesdm_model_dir` option.
#'   \strong{Required}.
#' @param r_environ (character) Path to an `.Renviron` file containing GBIF
#'   credentials. If `GBIF_USER`, `GBIF_PWD`, and `GBIF_EMAIL` are not set in the
#'   current R session, values are read from this file. If `NULL`, path is retrieved
#'   from the `onesdm_r_environ` option. Defaults to \emph{.Renviron}.
#' @param start_year (integer) Earliest year for which records are included.
#'   Defaults to `1981L`. Can also be set via the `onesdm_start_year` option.
#' @param boundaries (numeric) Vector of four spatial boundaries in decimal
#'   degrees, specified as (left, right, bottom, top). Defaults to
#'   \emph{c(-180L, 180L, -90L, 90L)} (global extent). Can also be set via the
#'   `onesdm_gbif_boundaries` option.
#' @param max_uncertainty (numeric) Maximum allowed spatial uncertainty in
#'   kilometers. Defaults to `10L`. Can also be set via the
#'   `onesdm_gbif_max_uncertainty` option.
#' @param overwrite (logical) If `TRUE`, overwrites existing GBIF data file.
#'    Defaults to `FALSE`. Can also be set via `onesdm_gbif_overwrite` option.
#' @param return_data (logical) If `TRUE`, returns processed GBIF data as an
#'   `sf` object in addition to saving it. Defaults to `FALSE`.
#' @param verbose (logical) If `TRUE` (default), prints progress messages during
#'   execution. Can also be set via the `onesdm_verbose` option.
#'
#' @details The function performs the following steps:
#' - Validates input parameters and environment variables.
#' - Checks for existing cleaned GBIF data and reuses it if available,
#'   unless `overwrite = TRUE`.
#' - Requests occurrence data from GBIF using the `rgbif` package, applying
#'   filters for coordinates, geospatial issues, year, occurrence status, basis
#'   of record, and spatial boundaries.
#' - Cleans downloaded data by:
#'   - Removing records with missing or low-precision coordinates.
#'   - Excluding records exceeding maximum spatial uncertainty.
#'   - Applying additional checks using the CoordinateCleaner package.
#' - Converts cleaned records to an `sf` object and saves it as an `.RData` file.
#' - Supports global configuration of default arguments via [base::options()]
#'   at the beginning of an R session, avoiding repeated specification in
#'   function calls. Available options and corresponding arguments:
#'   - `onesdm_gbif_ids`: `gbif_ids`
#'   - `onesdm_model_dir`: `model_dir`
#'   - `onesdm_gbif_verbose`: `verbose`
#'   - `onesdm_start_year`: `start_year`
#'   - `onesdm_r_environ`: `r_environ`
#'   - `onesdm_gbif_boundaries`: `boundaries`
#'   - `onesdm_gbif_max_uncertainty`: `max_uncertainty`
#'   - `onesdm_gbif_overwrite`: `overwrite`
#'
#' Example of setting [base::options()]:
#'   ```r
#'   options(
#'     onesdm_gbif_ids = c("793ZB", "6WRZD"),
#'     onesdm_model_dir = "path/to/model_dir",
#'     onesdm_gbif_verbose = TRUE,
#'     onesdm_start_year = 1981L,
#'     onesdm_r_environ = ".Renviron",
#'     onesdm_gbif_boundaries = c(-180L, 180L, -90L, 90L),
#'     onesdm_gbif_max_uncertainty = 10L,
#'     onesdm_gbif_overwrite = FALSE
#'   )
#'   ```
#'
#' @return If `return_data = TRUE`, returns cleaned GBIF data as an sf
#'   object. Otherwise, invisibly returns a named list of paths to saved files:
#' - `gbif_request`: GBIF request object (`data/gbif_request.RData`)
#' - `gbif_status`: GBIF download status (`data/gbif_status.RData`)
#' - `gbif_data_raw`: Raw GBIF download (`data/gbif_data_raw.zip`)
#' - `gbif_data`: Processed and cleaned GBIF data (`data/gbif_data.RData`).
#'
#' @examples
#' \dontrun{
#' prepare_gbif_data(
#'   gbif_ids = c("793ZB", "6WRZD"),
#'   model_dir = "path/to/model_dir"
#' )
#' }
#'
#' @export
#' @author Ahmed El-Gabbas, Maryna Golivets

prepare_gbif_data <- function(
  gbif_ids = NULL,
  model_dir = NULL,
  r_environ = ".Renviron",
  start_year = 1981L,
  boundaries = c(-180, 180, -90, 90),
  max_uncertainty = 10,
  overwrite = FALSE,
  return_data = FALSE,
  verbose = TRUE
) {
  .start_time <- lubridate::now(tzone = "CET")

  # Variables used in dplyr/data.table-like NSE expressions
  latitude <- longitude <- n_dec_long <- n_dec_lat <- NULL
  decimalLongitude <- decimalLatitude <- n_obs <- NULL
  year <- uncertain_km <- NULL
  id <- speciesKey <- name <- matched_name <- NULL
  occurrenceStatus <- coordinatePrecision <- coordinateUncertaintyInMeters <- NULL

  # ---------------------------------------------------------------------------
  # Packages
  # ---------------------------------------------------------------------------

  ecokit::check_packages(c(
    "CoordinateCleaner",
    "crayon",
    "curl",
    "fs",
    "lubridate",
    "purrr",
    "rcol",
    "readr",
    "rgbif",
    "sf",
    "stringr",
    "tidyselect",
    "dplyr",
    "tibble"
  ))

  # ---------------------------------------------------------------------------
  # Arguments from options
  # ---------------------------------------------------------------------------

  gbif_ids <- ecokit::assign_from_options(
    gbif_ids,
    "onesdm_gbif_ids",
    "character"
  )
  model_dir <- ecokit::assign_from_options(
    model_dir,
    "onesdm_model_dir",
    "character"
  )
  r_environ <- ecokit::assign_from_options(
    r_environ,
    "onesdm_r_environ",
    "character"
  )
  start_year <- ecokit::assign_from_options(
    start_year,
    "onesdm_start_year",
    c("numeric", "integer")
  )
  boundaries <- ecokit::assign_from_options(
    boundaries,
    "onesdm_gbif_boundaries",
    c("numeric", "integer")
  )
  max_uncertainty <- ecokit::assign_from_options(
    max_uncertainty,
    "onesdm_gbif_max_uncertainty",
    c("numeric", "integer")
  )
  overwrite <- ecokit::assign_from_options(
    overwrite,
    "onesdm_gbif_overwrite",
    "logical"
  )
  verbose <- ecokit::assign_from_options(
    verbose,
    "onesdm_gbif_verbose",
    "logical"
  )

  # ---------------------------------------------------------------------------
  # Validate arguments
  # ---------------------------------------------------------------------------

  if (
    is.null(model_dir) ||
      length(model_dir) != 1L ||
      is.na(model_dir) ||
      !nzchar(trimws(model_dir))
  ) {
    ecokit::stop_ctx(
      paste0(
        "The model_dir argument must be provided either directly or via the ",
        "`onesdm_model_dir` option."
      ),
      cat_timestamp = FALSE
    )
  }

  if (is.null(gbif_ids)) {
    ecokit::stop_ctx(
      "gbif_ids must contain at least one non-missing GBIF taxon key.",
      gbif_ids = gbif_ids,
      cat_timestamp = FALSE
    )
  }

  gbif_ids <- gbif_ids[!is.na(gbif_ids)]

  if (length(gbif_ids) == 0L) {
    ecokit::stop_ctx(
      "gbif_ids must contain at least one non-missing GBIF taxon key.",
      gbif_ids = gbif_ids,
      cat_timestamp = FALSE
    )
  }

  if (
    !is.numeric(start_year) ||
      length(start_year) != 1L ||
      is.na(start_year) ||
      !is.finite(start_year) ||
      start_year != floor(start_year) ||
      start_year < 1900 ||
      start_year > lubridate::year(lubridate::now())
  ) {
    ecokit::stop_ctx(
      "start_year must be one integer between 1900 and the current year.",
      start_year = start_year,
      cat_timestamp = FALSE
    )
  }
  start_year <- as.integer(start_year)

  if (
    !is.numeric(boundaries) ||
      length(boundaries) != 4L ||
      anyNA(boundaries) ||
      any(!is.finite(boundaries)) ||
      boundaries[1L] >= boundaries[2L] ||
      boundaries[3L] >= boundaries[4L] ||
      boundaries[1L] < -180 ||
      boundaries[2L] > 180 ||
      boundaries[3L] < -90 ||
      boundaries[4L] > 90
  ) {
    ecokit::stop_ctx(
      paste0(
        "boundaries must be finite numeric values (left, right, bottom, top), ",
        "within longitude [-180, 180] and latitude [-90, 90], with ",
        "left < right and bottom < top."
      ),
      boundaries = boundaries,
      cat_timestamp = FALSE
    )
  }

  if (
    !is.numeric(max_uncertainty) ||
      length(max_uncertainty) != 1L ||
      is.na(max_uncertainty) ||
      !is.finite(max_uncertainty) ||
      max_uncertainty < 0
  ) {
    ecokit::stop_ctx(
      "max_uncertainty must be one finite, non-negative numeric value.",
      max_uncertainty = max_uncertainty,
      cat_timestamp = FALSE
    )
  }

  for (arg in c("overwrite", "return_data", "verbose")) {
    value <- get(arg)
    if (!is.logical(value) || length(value) != 1L || is.na(value)) {
      ecokit::stop_ctx(
        paste0(arg, " must be one non-missing logical value."),
        cat_timestamp = FALSE
      )
    }
  }

  if (
    !is.null(r_environ) &&
      (!is.character(r_environ) ||
        length(r_environ) != 1L ||
        is.na(r_environ) ||
        !nzchar(r_environ))
  ) {
    ecokit::stop_ctx(
      "r_environ must be NULL or one non-empty character string.",
      cat_timestamp = FALSE
    )
  }

  # ---------------------------------------------------------------------------
  # GBIF IDs and taxon names
  # ---------------------------------------------------------------------------

  is_valid_gbif_id <- function(id) {
    result <- tryCatch(
      rcol::col_usage(id = id),
      error = function(e) NULL
    )
    !is.null(result) &&
      !is.null(result$id) &&
      identical(result$id, id)
  }

  gbif_ids_valid <- purrr::map_lgl(gbif_ids, is_valid_gbif_id)

  if (!all(gbif_ids_valid)) {
    ecokit::stop_ctx(
      paste0(
        "The following GBIF taxon keys do not exist: ",
        toString(gbif_ids[!gbif_ids_valid])
      ),
      gbif_ids = gbif_ids[!gbif_ids_valid],
      cat_timestamp = FALSE
    )
  }

  gbif_taxa <- purrr::map_dfr(
    sort(gbif_ids),
    function(id) {
      x <- rcol::col_usage(id)
      if (is.null(x) || nrow(x) == 0L) {
        return(tibble::tibble(id = id, name = NA_character_))
      }
      if (!"id" %in% names(x)) {
        x$id <- id
      }
      if (!"label" %in% names(x)) {
        x$name <- NA_character_
      }
      tibble::tibble(
        id = as.character(x$id[1L]),
        name = as.character(x$label[1L])
      )
    }
  )

  gbif_ids_names <- gbif_taxa |>
    dplyr::mutate(
      matched_name = paste0(id, " (", name, ")")
    ) |>
    dplyr::pull(matched_name) |>
    toString()

  ecokit::cat_time(
    "The following GBIF taxon keys will be processed:",
    cat_timestamp = FALSE,
    level = 1L,
    verbose = verbose
  )
  ecokit::cat_time(
    gbif_ids_names,
    cat_timestamp = FALSE,
    level = 2L,
    verbose = verbose
  )

  # ---------------------------------------------------------------------------
  # Paths and processed-data cache
  # ---------------------------------------------------------------------------

  fs::dir_create(model_dir)
  path_data <- fs::path(model_dir, "data")
  fs::dir_create(path_data)

  path_gbif_request <- fs::path(path_data, "gbif_request.RData")
  path_gbif_status <- fs::path(path_data, "gbif_status.RData")
  path_gbif_data_raw <- fs::path(path_data, "gbif_data_raw.zip")
  path_gbif_data <- fs::path(path_data, "gbif_data.RData")

  output_list <- list(
    gbif_request = path_gbif_request,
    gbif_status = path_gbif_status,
    gbif_data_raw = path_gbif_data_raw,
    gbif_data = path_gbif_data
  )

  if (ecokit::check_data(path_gbif_data, warning = FALSE) && !overwrite) {
    ecokit::cat_time(
      paste0(
        "Processed GBIF data already exist at: ",
        crayon::blue(path_gbif_data)
      ),
      cat_timestamp = FALSE,
      verbose = verbose
    )
    if (return_data) {
      return(invisible(ecokit::load_as(path_gbif_data)))
    }
    return(invisible(output_list))
  }

  if (ecokit::check_data(path_gbif_data, warning = FALSE) && overwrite) {
    ecokit::cat_time(
      "Existing GBIF data will be overwritten.",
      cat_timestamp = FALSE,
      verbose = verbose
    )
  }

  # ---------------------------------------------------------------------------
  # Print parameters
  # ---------------------------------------------------------------------------

  if (verbose) {
    ecokit::cat_time(
      "\nGBIF data extraction parameters",
      cat_timestamp = FALSE,
      cat_bold = TRUE,
      cat_red = TRUE
    )

    gbif_user <- Sys.getenv("GBIF_USER", unset = "")
    gbif_email <- Sys.getenv("GBIF_EMAIL", unset = "")
    if (!nzchar(gbif_user)) {
      gbif_user <- getOption("GBIF_USER", "")
    }
    if (!nzchar(gbif_email)) {
      gbif_email <- getOption("GBIF_EMAIL", "")
    }

    parameters <- c(
      "GBIF account" = if (nzchar(gbif_user)) gbif_user else "(not available)",
      "GBIF email" = if (nzchar(gbif_email)) gbif_email else "(not available)",
      "GBIF taxon key(s)" = gbif_ids_names,
      "Modelling directory" = model_dir,
      "Modelling directory (absolute)" = fs::path_abs(model_dir),
      "Start year" = as.character(start_year),
      ".Renviron file" = if (is.null(r_environ)) {
        "(not specified)"
      } else {
        r_environ
      },
      "Geographic boundaries" = toString(boundaries),
      "Maximum spatial uncertainty (km)" = as.character(max_uncertainty),
      "Overwrite" = as.character(overwrite),
      "Return processed data" = as.character(return_data)
    )

    purrr::iwalk(parameters, function(value, name) {
      ecokit::cat_time(
        paste0(crayon::italic(paste0(name, ": ")), crayon::blue(value)),
        level = 1L,
        cat_timestamp = FALSE
      )
    })

    ecokit::cat_time(
      "\nExtracting GBIF data",
      cat_timestamp = FALSE,
      cat_bold = TRUE,
      cat_red = TRUE
    )
  }

  ecokit::check_gbif(r_environ = r_environ)

  # ---------------------------------------------------------------------------
  # Reuse or submit GBIF request
  # ---------------------------------------------------------------------------

  request_status_exist <-
    ecokit::check_data(path_gbif_request, warning = FALSE) &&
    ecokit::check_data(path_gbif_status, warning = FALSE)

  if (request_status_exist && !overwrite) {
    gbif_request <- ecokit::load_as(path_gbif_request)
    gbif_status <- ecokit::load_as(path_gbif_status)
    ecokit::cat_time(
      "Reusing the saved GBIF request and status.",
      cat_timestamp = FALSE,
      verbose = verbose
    )
  } else {
    if (file.exists(path_gbif_data_raw)) {
      fs::file_delete(path_gbif_data_raw)
    }

    .request_start <- lubridate::now(tzone = "UTC")
    ecokit::cat_time(
      "Submitting GBIF occurrence download request.",
      cat_timestamp = FALSE,
      verbose = verbose
    )

    excluded_basis_of_record <- c("FOSSIL_SPECIMEN", "LIVING_SPECIMEN")
    data_boundary <- ecokit::boundary_to_wkt(
      left = boundaries[1L],
      right = boundaries[2L],
      bottom = boundaries[3L],
      top = boundaries[4L]
    )

    gbif_request <- rgbif::occ_download(
      rgbif::pred_in("taxonKey", gbif_ids),
      rgbif::pred("hasCoordinate", TRUE),
      rgbif::pred("hasGeospatialIssue", FALSE),
      rgbif::pred_gte("year", start_year),
      rgbif::pred("occurrenceStatus", "DETECTED"),
      rgbif::pred_not(
        rgbif::pred_in("BASIS_OF_RECORD", excluded_basis_of_record)
      ),
      rgbif::pred_within(data_boundary),
      format = "SIMPLE_CSV"
    )

    save(gbif_request, file = path_gbif_request)

    ecokit::cat_time(
      "Waiting for the GBIF download to be prepared.",
      cat_timestamp = FALSE,
      verbose = verbose
    )
    gbif_status <- rgbif::occ_download_wait(gbif_request, quiet = TRUE)
    save(gbif_status, file = path_gbif_status)

    ecokit::cat_diff(
      init_time = .request_start,
      prefix = "Submitting and preparing the GBIF request took ",
      level = 1L,
      verbose = verbose,
      cat_timestamp = FALSE
    )
  }

  # Check that GBIF reports a successful download and provides a URL
  if (
    !is.null(gbif_status$status) &&
      !identical(toupper(gbif_status$status), "SUCCEEDED")
  ) {
    ecokit::stop_ctx(
      paste0(
        "GBIF download is not ready; reported status: ",
        gbif_status$status
      ),
      cat_timestamp = FALSE
    )
  }
  if (
    is.null(gbif_status$downloadLink) ||
      length(gbif_status$downloadLink) != 1L ||
      is.na(gbif_status$downloadLink) ||
      !nzchar(gbif_status$downloadLink)
  ) {
    gbif_meta <- tryCatch(
      rgbif::occ_download_meta(key = gbif_status$key),
      error = function(e) NULL
    )
    if (!is.null(gbif_meta$downloadLink)) {
      gbif_status$downloadLink <- gbif_meta$downloadLink
    }
  }
  if (
    is.null(gbif_status$downloadLink) ||
      length(gbif_status$downloadLink) != 1L ||
      is.na(gbif_status$downloadLink) ||
      !nzchar(gbif_status$downloadLink)
  ) {
    ecokit::stop_ctx(
      "GBIF did not provide a valid download URL.",
      cat_timestamp = FALSE
    )
  }

  # ---------------------------------------------------------------------------
  # Download raw GBIF data
  # ---------------------------------------------------------------------------

  .download_start <- NULL
  raw_zip_valid <- ecokit::check_zip(path_gbif_data_raw, warning = FALSE)

  if (!raw_zip_valid || overwrite) {
    .download_start <- lubridate::now(tzone = "UTC")
    ecokit::cat_time(
      paste0(
        "Downloading raw GBIF data to: ",
        crayon::blue(path_gbif_data_raw)
      ),
      cat_timestamp = FALSE,
      verbose = verbose
    )

    curl::curl_download(
      url = gbif_status$downloadLink,
      destfile = path_gbif_data_raw,
      mode = "wb",
      quiet = TRUE
    )

    if (!ecokit::check_zip(path_gbif_data_raw, warning = FALSE)) {
      if (file.exists(path_gbif_data_raw)) {
        fs::file_delete(path_gbif_data_raw)
      }
      ecokit::stop_ctx(
        "The downloaded GBIF archive is not a valid ZIP file.",
        cat_timestamp = FALSE
      )
    }
  } else {
    ecokit::cat_time(
      paste0("Reusing raw GBIF archive: ", crayon::blue(path_gbif_data_raw)),
      cat_timestamp = FALSE,
      verbose = verbose
    )
  }

  # ---------------------------------------------------------------------------
  # Read raw GBIF data
  # ---------------------------------------------------------------------------

  ecokit::cat_time(
    "\nReading and processing GBIF occurrence data.",
    cat_timestamp = FALSE,
    verbose = verbose
  )

  gbif_data <- readr::read_tsv(
    file = path_gbif_data_raw,
    progress = FALSE,
    show_col_types = FALSE,
    col_types = readr::cols(
      gbifID = readr::col_double(),
      datasetKey = readr::col_character(),
      occurrenceID = readr::col_character(),
      kingdom = readr::col_character(),
      phylum = readr::col_character(),
      class = readr::col_character(),
      order = readr::col_character(),
      family = readr::col_character(),
      genus = readr::col_character(),
      species = readr::col_character(),
      infraspecificEpithet = readr::col_character(),
      taxonRank = readr::col_character(),
      scientificName = readr::col_character(),
      verbatimScientificName = readr::col_character(),
      verbatimScientificNameAuthorship = readr::col_character(),
      countryCode = readr::col_character(),
      locality = readr::col_character(),
      stateProvince = readr::col_character(),
      occurrenceStatus = readr::col_character(),
      individualCount = readr::col_double(),
      publishingOrgKey = readr::col_character(),
      decimalLatitude = readr::col_character(),
      decimalLongitude = readr::col_character(),
      coordinateUncertaintyInMeters = readr::col_double(),
      coordinatePrecision = readr::col_double(),
      elevation = readr::col_double(),
      elevationAccuracy = readr::col_double(),
      depth = readr::col_double(),
      depthAccuracy = readr::col_double(),
      eventDate = readr::col_character(),
      day = readr::col_integer(),
      month = readr::col_integer(),
      year = readr::col_integer(),
      taxonKey = readr::col_character(),
      speciesKey = readr::col_character(),
      basisOfRecord = readr::col_character(),
      institutionCode = readr::col_character(),
      collectionCode = readr::col_character(),
      catalogNumber = readr::col_character(),
      recordNumber = readr::col_character(),
      identifiedBy = readr::col_character(),
      dateIdentified = readr::col_datetime(format = ""),
      license = readr::col_character(),
      rightsHolder = readr::col_character(),
      recordedBy = readr::col_character(),
      typeStatus = readr::col_character(),
      establishmentMeans = readr::col_character(),
      lastInterpreted = readr::col_datetime(format = ""),
      mediaType = readr::col_character(),
      issue = readr::col_character()
    )
  )

  n_rows_raw <- nrow(gbif_data)
  ecokit::cat_time(
    paste0("Read ", ecokit::format_number(n_rows_raw), " raw GBIF records."),
    cat_timestamp = FALSE,
    verbose = verbose
  )

  # ---------------------------------------------------------------------------
  # Clean GBIF records
  # ---------------------------------------------------------------------------

  if (n_rows_raw > 0L) {
    gbif_data <- gbif_data |>
      dplyr::mutate(
        n_dec_long = ecokit::n_decimals(decimalLongitude),
        n_dec_lat = ecokit::n_decimals(decimalLatitude),
        uncertain_km = coordinateUncertaintyInMeters / 1000,
        longitude = suppressWarnings(as.numeric(decimalLongitude)),
        latitude = suppressWarnings(as.numeric(decimalLatitude)),
        .keep = "unused"
      )

    ecokit::cat_time(
      "\nCleaning GBIF records",
      cat_timestamp = FALSE,
      cat_bold = TRUE,
      cat_red = TRUE,
      verbose = verbose
    )

    # Both coordinates are required. Missing uncertainty is retained
    if (nrow(gbif_data) > 0L) {
      gbif_data <- CoordinateCleaner::cc_val(
        gbif_data,
        lon = "longitude",
        lat = "latitude",
        verbose = FALSE
      ) |>
        CoordinateCleaner::cc_equ(
          lon = "longitude",
          lat = "latitude",
          test = "identical",
          verbose = FALSE
        )
      ecokit::cat_time(
        paste0(
          "Coordinate validity: removed ",
          ecokit::format_number(n_rows_raw - nrow(gbif_data)),
          " records; ",
          ecokit::format_number(nrow(gbif_data)),
          " remain."
        ),
        level = 1L,
        cat_timestamp = FALSE,
        verbose = verbose
      )
    }

    report_filter <- function(data, label, filter_expr) {
      n_before <- nrow(data)
      data <- dplyr::filter(data, {{ filter_expr }})
      ecokit::cat_time(
        paste0(
          label,
          ": removed ",
          ecokit::format_number(n_before - nrow(data)),
          " records; ",
          ecokit::format_number(nrow(data)),
          " remain."
        ),
        level = 1L,
        cat_timestamp = FALSE,
        verbose = verbose
      )
      data
    }

    if (nrow(gbif_data) > 0L) {
      gbif_data <- report_filter(
        gbif_data,
        "Coordinate precision",
        n_dec_long > 1L & n_dec_lat > 1L
      )
    }
    if (nrow(gbif_data) > 0L) {
      gbif_data <- report_filter(
        gbif_data,
        "Spatial uncertainty",
        is.na(uncertain_km) | uncertain_km <= max_uncertainty
      )
    }

    # -------------------------------------------------------------------------
    # CoordinateCleaner additional filters
    # -------------------------------------------------------------------------

    if (nrow(gbif_data) > 0L) {
      n_before <- nrow(gbif_data)
      gbif_data <- gbif_data |>
        CoordinateCleaner::cc_cen(
          buffer = 100,
          lon = "longitude",
          lat = "latitude",
          verbose = FALSE
        ) |>
        CoordinateCleaner::cc_cap(
          buffer = 100,
          lon = "longitude",
          lat = "latitude",
          verbose = FALSE
        ) |>
        CoordinateCleaner::cc_inst(
          buffer = 100,
          geod = TRUE,
          lon = "longitude",
          lat = "latitude",
          verbose = FALSE
        ) |>
        CoordinateCleaner::cc_gbif(
          buffer = 100,
          lon = "longitude",
          lat = "latitude",
          verbose = FALSE
        )

      ecokit::cat_time(
        paste0(
          "CoordinateCleaner: removed ",
          ecokit::format_number(n_before - nrow(gbif_data)),
          " records due to their proximity to capitals, ",
          "centroids of countries/provinces, and biodiversity institutions; ",
          ecokit::format_number(nrow(gbif_data)),
          " remain."
        ),
        level = 1L,
        cat_timestamp = FALSE,
        verbose = verbose
      )
    }
  }

  # ---------------------------------------------------------------------------
  # Convert to sf and save, including empty results
  # ---------------------------------------------------------------------------

  if (nrow(gbif_data) > 0L) {
    gbif_data <- sf::st_as_sf(
      gbif_data,
      coords = c("longitude", "latitude"),
      crs = 4326,
      remove = FALSE
    )
  } else {
    ecokit::cat_time(
      "No GBIF records remain after filtering.",
      cat_timestamp = FALSE,
      verbose = verbose
    )
  }

  ecokit::cat_time(
    paste0("Saving processed GBIF data to: ", crayon::blue(path_gbif_data)),
    level = 1L,
    cat_timestamp = FALSE,
    verbose = verbose
  )
  save(gbif_data, file = path_gbif_data)

  if (!is.null(.download_start)) {
    ecokit::cat_diff(
      init_time = .download_start,
      prefix = "Downloading and processing GBIF data took ",
      level = 1L,
      verbose = verbose,
      cat_timestamp = FALSE
    )
  }

  # ---------------------------------------------------------------------------
  # Summary and return
  # ---------------------------------------------------------------------------

  if (nrow(gbif_data) > 0L && "speciesKey" %in% names(gbif_data)) {
    n_obs_per_species <- sf::st_drop_geometry(gbif_data) |>
      dplyr::count(speciesKey, name = "n_obs") |>
      dplyr::mutate(
        summary = paste0(
          speciesKey,
          " (",
          format(n_obs, big.mark = ",", scientific = FALSE),
          ")"
        )
      ) |>
      dplyr::pull(summary) |>
      toString()
  } else {
    n_obs_per_species <- "none"
  }

  ecokit::cat_time(
    paste0(
      "DONE! A total of ",
      ecokit::format_number(nrow(gbif_data), underline = TRUE),
      " filtered observations were extracted for GBIF taxon key(s): ",
      crayon::blue(n_obs_per_species)
    ),
    cat_timestamp = FALSE,
    level = 1L,
    verbose = verbose
  )

  ecokit::cat_diff(
    init_time = .start_time,
    prefix = "Extracting GBIF data took ",
    level = 1L,
    verbose = verbose,
    cat_timestamp = FALSE
  )

  if (return_data) {
    return(invisible(gbif_data))
  }
  invisible(output_list)
}
