# prepare_gbif_data -------

#' @title Download and Process GBIF Occurrence Data
#'
#' @description Downloads, filters, and processes occurrence records from the
#'   Global Biodiversity Information Facility ([GBIF](https://gbif.org/)) for
#'   specified taxon keys (GBIF IDs). The function handles data requests,
#'   downloading, cleaning, and conversion to an `sf` object, with optional
#'   spatial and temporal filters, coordinate uncertainty thresholds, and
#'   geographic boundary constraints.
#'
#' @param gbif_ids \emph{(character or numeric)} A vector of one or more GBIF
#'   taxon keys. When multiple IDs are supplied, data are combined across all keys.
#'   If `NULL` (default), the function attempts to retrieve IDs from the
#'   `onesdm_gbif_ids` option. \strong{Required}.
#' @param model_dir \emph{(character)} Path to the directory where model outputs
#'   will be saved. A subdirectory named `data` is automatically created within this
#'   directory to store processed species data. When modelling multiple species,
#'   it is recommended to use a separate directory for each run to avoid overwriting
#'   or mixing data files. This path can also be set via the `onesdm_model_dir`
#'   option.  Default is `NULL`. \strong{Required}.
#' @param r_environ \emph{(character)} Path to an `.Renviron` file containing GBIF
#'   credentials (default: ".Renviron"). The function uses
#'   [ecokit::check_gbif()] to verify the presence of GBIF credentials. If the
#'   `GBIF_USER`, `GBIF_PWD`, and `GBIF_EMAIL` environment variables are not set
#'   in the current R session, `check_gbif` attempts to read them from the
#'   specified `.Renviron` file. See [ecokit::check_gbif()] and
#'   this
#'   [`rgbif` documentation](https://docs.ropensci.org/rgbif/articles/gbif_credentials.html)
#'   for details on setting GBIF credentials. If `r_environ = NULL`, the path to
#'   the `.Renviron` file is retrieved from the "`onesdm_r_environ`" option.
#' @param start_year \emph{(integer)} Include only records from this year onward. The
#'   default is `1981L`, to match the temporal coverage of CHELSA climate data.
#'   Can also be set via the `onesdm_start_year` option.
#' @param boundaries \emph{(numeric)} A vector of size 4 containing spatial boundaries
#'   as (left, right, bottom, top) in decimal degrees (default:
#'   `c(-180L, 180L, -90L, 90L)` for global extent). Can also be set via
#'   `onesdm_gbif_boundaries` option.
#' @param max_uncertainty \emph{(numeric)} Maximum allowed spatial uncertainty in
#'   kilometres. Default is `10L`. Can also be set via
#'   `onesdm_gbif_max_uncertainty` option.
#' @param overwrite \emph{(logical)} If `TRUE`, overwrites existing cleaned GBIF data
#'   file in the model directory. Default is `FALSE`. Can also be set via
#'   `onesdm_gbif_overwrite` option.
#' @param return_data \emph{(logical)} If `TRUE`, returns the processed GBIF data as an
#'   `sf` object in addition to saving it. Default is `FALSE`.
#' @param verbose \emph{(logical)} If `TRUE` (default), prints progress messages during
#'   execution. Can also be set via the `onesdm_verbose` option. Default is `"TRUE"`.
#'
#' @details The function performs the following steps:
#' - Validates input parameters and environment.
#' - Checks for existing processed data to avoid redundant downloads/requests.
#' - Requests occurrence data from GBIF using the `rgbif` package, applying
#'   filters for coordinates, geospatial issues, year, occurrence status, basis
#'   of record, and spatial boundaries.
#' - Downloads and reads the raw data, applying further cleaning:
#'   - Removes records with missing or low-precision coordinates.
#'   - Filters by spatial uncertainty and accepted taxonomic ranks.
#'   - Applies additional cleaning using the `CoordinateCleaner` package.
#' - Converts the cleaned data to an `sf` object and saves it.
#'
#' - Function default arguments can be set globally using the `options()`
#'   function. Users can set these options at the start of their R session to
#'   avoid repeatedly specifying them in function calls. The following options
#'   correspond to the function arguments:
#'   - `onesdm_gbif_ids`: `gbif_ids`
#'   - `onesdm_model_dir`: `model_dir`
#'   - `onesdm_gbif_verbose`: `verbose`
#'   - `onesdm_start_year`: `start_year`
#'   - `onesdm_r_environ`: `r_environ`
#'   - `onesdm_gbif_boundaries`: `boundaries`
#'   - `onesdm_gbif_max_uncertainty`: `max_uncertainty`
#'   - `onesdm_gbif_overwrite`: `overwrite`
#' - Example of setting options:
#'   ```r
#'   options(
#'     onesdm_gbif_ids = c("1234567", "7654321"),
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
#' @return If `return_data` is `TRUE`, return the cleaned GBIF data as an `sf`
#'   object. Otherwise, returns (invisibly) a named list with paths to the saved
#'   GBIF data files:
#' - `gbif_request`: Path to the saved GBIF request object
#'   (`data/gbif_request.RData`).
#' - `gbif_status`: Path to the saved GBIF status object
#'   (`data/gbif_status.RData`).
#' - `gbif_data_raw`: Path to the downloaded raw GBIF data
#'   (`data/gbif_data_raw.zip`).
#' - `gbif_data`: Path to the processed and cleaned GBIF data
#'   (`data/gbif_data.RData`).
#'
#' @examples
#' \dontrun{
#' prepare_gbif_data(
#'   gbif_ids = c("1234567", "7654321"), model_dir = "path/to/model_dir")
#' }
#'
#' @export
#' @author Ahmed El-Gabbas

prepare_gbif_data <- function(
  gbif_ids = NULL,
  model_dir = NULL,
  r_environ = ".Renviron",
  start_year = 1981L,
  boundaries = c(-180L, 180L, -90L, 90L),
  max_uncertainty = 10L,
  overwrite = FALSE,
  return_data = FALSE,
  verbose = TRUE
) {
  .start_gbif_time <- lubridate::now(tzone = "CET")

  latitude <- longitude <- year <- uncertain_km <- n_dec_long <- n_dec_lat <-
    key <- matched_name <- decimalLongitude <- decimalLatitude <- #nolint
      occurrenceStatus <- taxonRank <- coordinatePrecision <- #nolint
        scientificName <- coordinateUncertaintyInMeters <- NULL #nolint

  ecokit::check_packages(
    c(
      "CoordinateCleaner",
      "crayon",
      "curl",
      "fs",
      "lubridate",
      "purrr",
      "readr",
      "rgbif",
      "sf",
      "stringr",
      "tidyselect"
    )
  )

  # # ********************************************************************** #
  # Check inputs and environment ------
  # # ********************************************************************** #

  gbif_ids <- ecokit::assign_from_options(
    gbif_ids,
    "onesdm_gbif_ids",
    c("character", "numeric", "integer")
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

  # # ||||||||||||||||||||||||||||||||||||||||| #
  ## Check boundaries ------
  # # ||||||||||||||||||||||||||||||||||||||||| #

  if (
    length(boundaries) != 4L || anyNA(boundaries) || !is.numeric(boundaries)
  ) {
    ecokit::stop_ctx(
      "The boundaries argument must be a numeric vector of length 4",
      boundaries = boundaries
    )
  }

  if (
    boundaries[1L] >= boundaries[2L] ||
      boundaries[3L] >= boundaries[4L] ||
      boundaries[1L] < -180L ||
      boundaries[2L] > 180L ||
      boundaries[3L] < -90L ||
      boundaries[4L] > 90L
  ) {
    ecokit::stop_ctx(
      paste0(
        "The boundaries argument must be a numeric vector of length 4 with ",
        "values (Left, Right, Bottom, Top) within the ranges: ",
        "Left [-180, 180], Right [-180, 180], Bottom [-90, 90],",
        "Top [-90, 90], and Left < Right and Bottom < Top."
      ),
      boundaries = boundaries,
      cat_timestamp = FALSE
    )
  }

  # # ||||||||||||||||||||||||||||||||||||||||| #
  ## Check r_environ / GBIF credentials ------
  # # ||||||||||||||||||||||||||||||||||||||||| #

  ecokit::check_gbif(r_environ = r_environ)

  # # ||||||||||||||||||||||||||||||||||||||||| #
  ## Check start_year ------
  # # ||||||||||||||||||||||||||||||||||||||||| #

  current_year <- lubridate::year(lubridate::now())
  if (
    !is.numeric(start_year) ||
      length(start_year) != 1L ||
      is.na(start_year) ||
      start_year < 1900L ||
      start_year > current_year
  ) {
    ecokit::stop_ctx(
      paste0(
        "The start_year argument must be a single numeric value between",
        "1900 and ",
        current_year,
        "."
      ),
      start_year = start_year,
      cat_timestamp = FALSE
    )
  }

  # # ||||||||||||||||||||||||||||||||||||||||| #
  ## Check model_dir ------
  # # ||||||||||||||||||||||||||||||||||||||||| #

  if (is.null(model_dir)) {
    ecokit::stop_ctx(
      paste0(
        "The model_dir argument must be provided either directly or via the ",
        "onesdm_model_dir option."
      ),
      cat_timestamp = FALSE
    )
  }

  # # ||||||||||||||||||||||||||||||||||||||||| #
  ## Check gbif_ids ------
  # # ||||||||||||||||||||||||||||||||||||||||| #

  if (is.character(gbif_ids)) {
    gbif_ids <- as.integer(gbif_ids)
  }

  if (!all(stringr::str_detect(gbif_ids, "^\\d+$"))) {
    #nolint
    ecokit::stop_ctx(
      "Each gbif_id can contain only numbers.",
      gbif_ids = gbif_ids,
      cat_timestamp = FALSE
    )
  }

  gbif_ids_valid <- purrr::map_lgl(gbif_ids, is_valid_gbif_id)
  if (!all(gbif_ids_valid)) {
    ecokit::stop_ctx(
      paste0(
        "The following gbif_ids do not exist in GBIF: ",
        toString(gbif_ids[!gbif_ids_valid])
      ),
      gbif_ids = gbif_ids,
      cat_timestamp = FALSE
    )
  }

  gbif_ids_names <- purrr::map(
    .x = sort(gbif_ids),
    .f = ~ rgbif::name_usage(.x)$data
  ) %>%
    dplyr::bind_rows() %>%
    dplyr::select(tidyselect::all_of(c("key", "scientificName"))) %>%
    dplyr::mutate(matched_name = paste0(key, " (", scientificName, ")")) %>%
    dplyr::pull(matched_name) %>%
    toString()

  ecokit::cat_time(
    "The following GBIF ID(s) will be processed:",
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

  # # ********************************************************************** #
  # Print function arguments ------
  # # ********************************************************************** #

  if (verbose) {
    ecokit::cat_time(
      crayon::italic("\nGBIF data extraction parameters:"),
      cat_timestamp = FALSE,
      cat_bold = TRUE,
      cat_red = TRUE
    )

    gbif_user <- Sys.getenv("GBIF_USER")
    if (is.null(gbif_user)) {
      gbif_user <- getOption("GBIF_USER")
    }
    gbif_email <- Sys.getenv("GBIF_EMAIL")
    if (is.null(gbif_email)) {
      gbif_email <- getOption("GBIF_EMAIL")
    }
    gbif_details <- paste0(gbif_user, " - ", gbif_email)
    ecokit::cat_time(
      paste0(crayon::italic("GBIF user details: "), crayon::blue(gbif_details)),
      level = 1L,
      cat_timestamp = FALSE
    )

    ecokit::cat_time(
      paste0(crayon::italic("GBIF ID(s): "), crayon::blue(gbif_ids_names)),
      level = 1L,
      cat_timestamp = FALSE
    )
    ecokit::cat_time(
      paste0(crayon::italic("Modelling directory: "), crayon::blue(model_dir)),
      level = 1L,
      cat_timestamp = FALSE
    )
    ecokit::cat_time(
      paste0(
        crayon::italic("Modelling directory (absolute): "),
        crayon::blue(fs::path_abs(model_dir))
      ),
      level = 1L,
      cat_timestamp = FALSE
    )
    ecokit::cat_time(
      paste0(crayon::italic("Start year: "), crayon::blue(start_year)),
      level = 1L,
      cat_timestamp = FALSE
    )
    ecokit::cat_time(
      paste0(crayon::italic(".Renviron file: "), crayon::blue(r_environ)),
      level = 1L,
      cat_timestamp = FALSE
    )
    ecokit::cat_time(
      paste0(
        crayon::italic("Geographic boundaries: "),
        crayon::blue(toString(boundaries))
      ),
      level = 1L,
      cat_timestamp = FALSE
    )
    ecokit::cat_time(
      paste0(
        crayon::italic("Overwrite exisiting data: "),
        crayon::blue(overwrite)
      ),
      level = 1L,
      cat_timestamp = FALSE
    )
    ecokit::cat_time(
      paste0(
        crayon::italic("Return processed data: "),
        crayon::blue(return_data)
      ),
      level = 1L,
      cat_timestamp = FALSE
    )
    ecokit::cat_time(
      paste0(
        crayon::italic("Maximum allowed spatial uncertainty: "),
        crayon::blue(max_uncertainty),
        " km"
      ),
      level = 1L,
      cat_timestamp = FALSE
    )
    ecokit::cat_time(
      "\nExtracting GBIF data",
      cat_timestamp = FALSE,
      cat_bold = TRUE,
      cat_red = TRUE
    )
  }

  # # ********************************************************************** #
  # Check if data exist ------
  # # ********************************************************************** #

  path_data <- fs::path(model_dir, "data")
  path_gbif_request <- fs::path(path_data, "gbif_request.RData")
  path_gbif_status <- fs::path(path_data, "gbif_status.RData")
  path_gbif_data <- fs::path(path_data, "gbif_data.RData")
  path_gbif_data_raw <- fs::path(path_data, "gbif_data_raw.zip")

  output_list <- list(
    gbif_request = path_gbif_request,
    gbif_status = path_gbif_status,
    gbif_data_raw = path_gbif_data_raw,
    gbif_data = path_gbif_data
  )
  all_okay <- all(
    ecokit::check_data(path_gbif_request, warning = FALSE),
    ecokit::check_data(path_gbif_status, warning = FALSE),
    ecokit::check_data(path_gbif_data, warning = FALSE)
  )

  if (all_okay) {
    if (!overwrite) {
      ecokit::cat_time(
        paste0(
          "GBIF data already exist at: ",
          crayon::blue(path_gbif_data)
        ),
        cat_timestamp = FALSE,
        verbose = verbose
      )

      if (return_data) {
        ecokit::cat_time(
          "Loading GBIF data",
          cat_timestamp = FALSE,
          verbose = verbose
        )
        return(invisible(ecokit::load_as(path_gbif_data)))
      } else {
        return(invisible(output_list))
      }
    }
    ecokit::cat_time(
      paste0(
        "GBIF data already exist. ",
        "Overwriting as per the `overwrite = TRUE` argument."
      ),
      cat_timestamp = FALSE,
      verbose = verbose
    )
  }

  fs::dir_create(path_data)

  # # ********************************************************************** #
  # Request GBIF data ------
  # # ********************************************************************** #

  request_status_okay <- ecokit::check_data(
    path_gbif_request,
    warning = FALSE
  ) &&
    ecokit::check_data(path_gbif_status, warning = FALSE)

  if (request_status_okay && !overwrite) {
    ecokit::cat_time(
      "Loading GBIF data",
      verbose = verbose,
      cat_timestamp = FALSE
    )
    gbif_request <- ecokit::load_as(path_gbif_request)
    gbif_status <- ecokit::load_as(path_gbif_status)
  } else {
    .start_time_request <- lubridate::now(tzone = "CET")

    ecokit::cat_time(
      "Requesting GBIF data (this may take some time, depending on the data volume)",
      verbose = verbose,
      cat_timestamp = FALSE
    )

    excluded_basis_of_record <- c("FOSSIL_SPECIMEN", "LIVING_SPECIMEN")
    data_boundary <- ecokit::boundary_to_wkt(
      left = boundaries[1L],
      right = boundaries[2L],
      bottom = boundaries[3L],
      top = boundaries[4L]
    )

    gbif_request <- rgbif::occ_download(
      # Species keys
      rgbif::pred_in("taxonKey", gbif_ids),
      # Only with coordinates & no spatial issues
      rgbif::pred("hasCoordinate", TRUE),
      rgbif::pred("hasGeospatialIssue", FALSE),
      # Only after (>=) a certain year
      rgbif::pred_gte("year", start_year),
      # Only "PRESENCE" data (i.e. exclude ABSENCE)
      rgbif::pred("occurrenceStatus", "PRESENT"),
      # Only specific basis of record
      rgbif::pred_not(
        rgbif::pred_in("BASIS_OF_RECORD", excluded_basis_of_record)
      ),
      # Only within specific boundaries
      rgbif::pred_within(data_boundary),
      # File size of "SIMPLE_CSV" format can be much smaller than "DWCA"
      format = "SIMPLE_CSV"
    )

    # # ||||||||||||||||||||||||||||||||||||||||| #
    ## Save data request -------
    # # ||||||||||||||||||||||||||||||||||||||||| #

    ecokit::cat_time(
      "Saving data request",
      level = 1L,
      verbose = verbose,
      cat_timestamp = FALSE
    )
    save(gbif_request, file = path_gbif_request)

    # # ||||||||||||||||||||||||||||||||||||||||| #
    # Wait for data to be ready ------
    # # ||||||||||||||||||||||||||||||||||||||||| #

    ecokit::cat_time(
      "Waiting for data to be ready .... ",
      level = 1L,
      verbose = verbose,
      cat_timestamp = FALSE
    )
    gbif_status <- rgbif::occ_download_wait(gbif_request, quiet = TRUE)

    ecokit::cat_diff(
      init_time = .start_time_request,
      prefix = "Requesting GBIF data took ",
      level = 1L,
      verbose = verbose
    )

    # # ||||||||||||||||||||||||||||||||||||||||| #
    # Save status details ------
    # # ||||||||||||||||||||||||||||||||||||||||| #

    ecokit::cat_time(
      "Saving status details",
      level = 1L,
      cat_timestamp = FALSE,
      verbose = verbose
    )
    save(gbif_status, file = fs::path(path_data, "gbif_status.RData"))

    ecokit::cat_time(
      "Data are ready - status summary:",
      ... = "\n",
      level = 1L,
      verbose = verbose,
      cat_timestamp = FALSE
    )
    print(rgbif::occ_download_meta(key = gbif_status$key))
  }

  # # ********************************************************************** #
  # Download GBIF data -------
  # # ********************************************************************** #

  if (ecokit::check_zip(path_gbif_data_raw, warning = FALSE) && !overwrite) {
    ecokit::cat_time(
      paste0(
        "Raw GBIF raw data already exist at: ",
        crayon::blue(path_gbif_data_raw)
      ),
      cat_timestamp = FALSE,
      verbose = verbose
    )
  } else {
    ecokit::cat_time(
      "\nDownload GBIF data",
      cat_timestamp = FALSE,
      verbose = verbose
    )
    .start_time_download <- lubridate::now(tzone = "CET")

    # download file to path_gbif_data_raw
    ecokit::cat_time(
      paste0(
        "Downloading raw GBIF data to: ",
        crayon::blue(path_gbif_data_raw)
      ),
      verbose = verbose,
      level = 1L,
      cat_timestamp = FALSE
    )
    curl::curl_download(
      url = gbif_status$downloadLink,
      destfile = path_gbif_data_raw,
      mode = "wb",
      quiet = TRUE
    )
  }

  # # ********************************************************************** #
  # Process GBIF data -------
  # # ********************************************************************** #

  if (!ecokit::check_data(path_gbif_data, warning = FALSE) || overwrite) {
    ecokit::cat_time(
      "Preparing GBIF data as an `sf` object",
      verbose = verbose,
      cat_timestamp = FALSE
    )

    # # ||||||||||||||||||||||||||||||||||||||||| #
    ## Read GBIF data -------
    # # ||||||||||||||||||||||||||||||||||||||||| #

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
        # Read coordinates as character to keep the original precision
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
        taxonKey = readr::col_double(),
        speciesKey = readr::col_double(),
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
    if (n_rows_raw == 0L) {
      ecokit::cat_time(
        "No GBIF data were extracted after reading the raw data",
        cat_timestamp = FALSE,
        level = 1L,
        verbose = verbose
      )
      save(gbif_data, file = path_gbif_data)
      return(invisible(output_list))
    }

    # # ||||||||||||||||||||||||||||||||||||||||| #
    # Clean GBIF data -------
    # # ||||||||||||||||||||||||||||||||||||||||| #

    gbif_data <- gbif_data %>%
      dplyr::rename(
        longitude = decimalLongitude,
        latitude = decimalLatitude,
        uncertain_km = coordinateUncertaintyInMeters
      ) %>%
      dplyr::mutate(
        # Number of decimal places for longitude / latitude
        n_dec_long = ecokit::n_decimals(longitude),
        n_dec_lat = ecokit::n_decimals(latitude),
        # Convert uncertainty to kilometers
        uncertain_km = uncertain_km / 1000L,
        # Convert coordinates to numeric
        longitude = as.numeric(longitude),
        latitude = as.numeric(latitude)
      ) %>%
      dplyr::filter(
        # Exclude occurrences with empty coordinates
        !is.na(longitude) | !is.na(latitude),
        # Exclude high spatial uncertainty (keep empty uncertainty values)
        uncertain_km <= max_uncertainty | is.na(uncertain_km),
        # Exclude occurrences with low precision (keep empty precision values)
        # coordinatePrecision <= 0.05 | is.na(coordinatePrecision),
        # Exclude occurrences if either latitude/longitude has 1 or 0 decimals
        (n_dec_long > 1L & n_dec_lat > 1L),
        # Exclude equal coordinates
        longitude != latitude,
        # Only occurrences recorded after specific year
        year >= start_year,
        # Only "PRESENT" data (i.e. exclude ABSENT)
        occurrenceStatus == "PRESENT" | is.na(occurrenceStatus),
        # Only accepted taxonomy ranks (species level and below)
        taxonRank %in% c("FORM", "SPECIES", "SUBSPECIES", "VARIETY")
      )

    # # ||||||||||||||||||||||||||||||||||||||||| #
    # Clean data further using CoordinateCleaner package ------
    # # ||||||||||||||||||||||||||||||||||||||||| #

    if (nrow(gbif_data) > 0L) {
      gbif_data <- gbif_data %>%
        # Exclude coordinates in the vicinity of country and province centroids
        CoordinateCleaner::cc_cen(
          buffer = 100L,
          lon = "longitude",
          lat = "latitude",
          verbose = FALSE
        ) %>%
        # Exclude coordinates in the vicinity of country capitals
        CoordinateCleaner::cc_cap(
          buffer = 100L,
          lon = "longitude",
          lat = "latitude",
          verbose = FALSE
        ) %>%
        # Exclude records in the vicinity of biodiversity institutions
        CoordinateCleaner::cc_inst(
          lon = "longitude",
          lat = "latitude",
          verbose = FALSE
        ) %>%
        # Exclude records assigned to GBIF Headquarters
        CoordinateCleaner::cc_gbif(
          buffer = 100L,
          lon = "longitude",
          lat = "latitude",
          verbose = FALSE
        ) %>%
        # Exclude records with identical lat/lon
        CoordinateCleaner::cc_equ(
          lon = "longitude",
          lat = "latitude",
          verbose = FALSE
        )
    }

    n_rows_cleaned <- nrow(gbif_data)

    if (n_rows_cleaned == 0L) {
      ecokit::cat_time(
        "No GBIF data were extracted after cleaning the raw data",
        cat_timestamp = FALSE,
        level = 1L,
        verbose = verbose
      )
      save(gbif_data, file = path_gbif_data)
      return(invisible(output_list))
    }

    ecokit::cat_time(
      paste0(
        "A total of ",
        ecokit::format_number(n_rows_cleaned, underline = TRUE),
        " observations were kept after cleaning from the initial ",
        ecokit::format_number(n_rows_raw, underline = TRUE),
        " observations."
      ),
      cat_timestamp = FALSE,
      level = 1L,
      verbose = verbose
    )

    # # ||||||||||||||||||||||||||||||||||||||||| #
    # Convert to sf object ------
    # # ||||||||||||||||||||||||||||||||||||||||| #

    gbif_data <- sf::st_as_sf(
      gbif_data,
      coords = c("longitude", "latitude"),
      crs = 4326L,
      remove = FALSE
    )

    # # ||||||||||||||||||||||||||||||||||||||||| #
    # Save data ------
    # # ||||||||||||||||||||||||||||||||||||||||| #

    ecokit::cat_time(
      paste0("Saving GBIF data to: ", crayon::blue(path_gbif_data)),
      verbose = verbose,
      level = 1L,
      cat_timestamp = FALSE
    )
    save(gbif_data, file = path_gbif_data)

    ecokit::cat_diff(
      init_time = .start_time_download,
      prefix = "\nDownloading GBIF data took ",
      verbose = verbose,
      level = 1L,
      cat_timestamp = FALSE
    )
  }

  # # ||||||||||||||||||||||||||||||||||||||||| #
  ## Prepare data summary -----
  # # ||||||||||||||||||||||||||||||||||||||||| #

  ecokit::cat_diff(
    init_time = .start_gbif_time,
    prefix = "\nExtracting GBIF data took ",
    verbose = verbose,
    level = 1L,
    cat_timestamp = FALSE
  )

  ecokit::cat_time(
    paste0(
      "A total of ",
      ecokit::format_number(nrow(gbif_data), underline = TRUE),
      " filtered observations were extracted for GBIF ID(s): ",
      crayon::blue(toString(sort(unique(gbif_data$speciesKey))))
    ),
    cat_timestamp = FALSE,
    level = 1L,
    verbose = verbose
  )

  # # ********************************************************************** #

  if (return_data) {
    return(invisible(gbif_data))
  } else {
    return(invisible(output_list))
  }
}
