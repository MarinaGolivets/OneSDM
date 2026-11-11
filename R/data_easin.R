# prepare_easin_data -------
#
#' @title Download and Process EASIN Occurrence Data
#'
#' @description Downloads, combines, cleans, and saves occurrence records
#'   from the European Alien Species Information Network
#'   ([EASIN](https://easin.jrc.ec.europa.eu/easin)) Geodatabase for the
#'   supplied EASIN IDs.
#'
#' @param easin_ids (character) Vector of one or more EASIN IDs,
#'   each starting with "R" followed by five digits (e.g., "R00544"). IDs can
#'   be obtained by searching for a taxon name on the
#'   [EASIN website](https://easin.jrc.ec.europa.eu/spexplorer/search/)
#'   by searching for a taxon name. Data for multiple IDs are combined.
#'   If `NULL` (default), IDs are retrieved from the `onesdm_easin_ids` option.
#'   \strong{Required}.
#' @param model_dir (character) Directory for saving model outputs.
#'   A `data` subdirectory is created automatically to store raw and processed data.
#'   If `NULL` (default), path is retrieved from the `onesdm_model_dir` option.
#'   \strong{Required}.
#' @param timeout (integer) Timeout (seconds) per download attempt and chunk.
#'   Default is `600L`. Does not limit total download duration. Can also be set via
#'   the `onesdm_easin_timeout` option.
#' @param n_search (integer) Number of records requested per API call (chunk
#'   size). Defaults to `1000L` (EASIN API maximum). Can also be set via the
#'   `onesdm_easin_n_search` option.
#' @param n_attempts (integer) Maximum number of download attempts per chunk.
#'   Defaults to `10L`. Can also be set via the `onesdm_easin_n_attempts`
#'   option.
#' @param sleep_time (integer) Delay (seconds) between successive
#'   chunk downloads to prevent overloading the EASIN server.
#'   Defaults to `5L`. Can also be set via the `onesdm_easin_sleep_time` option.
#' @param exclude_gbif (logical) If `TRUE` (default), excludes records
#'   originating from [GBIF](https://www.gbif.org/). Can also be set via the
#'   `onesdm_easin_exclude_gbif` option.
#' @param start_year (integer) Earliest year for which records are included.
#'   Defaults to `1981L`. Can also be set via the `onesdm_start_year` option.
#' @param overwrite (logical) If `TRUE`, overwrites existing cleaned EASIN data.
#'   Defaults to `FALSE`. Can also be set via `onesdm_easin_overwrite` option.
#' @param return_data (logical) If `TRUE`, returns the processed EASIN data as an
#'   `sf` object in addition to saving it. Defaults to `FALSE`.
#' @param verbose (logical) If `TRUE` (default), prints progress and
#'   informational messages, including URL of currently processed chunk.
#'   Can also be set via the `onesdm_easin_verbose` option.
#'
#' @details
#' - Validates input parameters.
#' - Supports chunked downloads and retries.
#' - Checks for existing cleaned EASIN data and reuses it if available,
#'   unless `overwrite = TRUE`.
#' - Cleans downloaded data by:
#'   - Removing records with missing or low-precision coordinates.
#'   - Excluding records exceeding maximum spatial uncertainty.
#'   - Applying additional checks using the CoordinateCleaner package.
#' - Converts cleaned records to an `sf` object and saves it as an `.RData` file.
#' - Supports global configuration of default arguments via [base::options()]
#'   at the beginning of an R session, avoiding repeated specification in
#'   function calls. Available options and corresponding arguments:
#'   - `onesdm_easin_ids`: `easin_ids`
#'   - `onesdm_model_dir`: `model_dir`
#'   - `onesdm_easin_timeout`: `timeout`
#'   - `onesdm_easin_n_search`: `n_search`
#'   - `onesdm_easin_n_attempts`: `n_attempts`
#'   - `onesdm_easin_sleep_time`: `sleep_time`
#'   - `onesdm_easin_exclude_gbif`: `exclude_gbif`
#'   - `onesdm_start_year`: `start_year`
#'   - `onesdm_easin_verbose`: `verbose`
#'   - `onesdm_easin_overwrite`: `overwrite`
#'
#' Example of setting [base::options()]:
#'   ```r
#'   options(
#'     onesdm_easin_ids = c("R00042", "R00544"),
#'     onesdm_model_dir = "path/to/model_dir",
#'     onesdm_easin_exclude_gbif = TRUE,
#'     onesdm_easin_overwrite = TRUE
#'   )
#'   ```
#'
#' @return Invisibly returns path to saved cleaned EASIN data (`data/easin_data.RData`
#'   within `model_dir`). If `return_data = TRUE`, returns cleaned EASIN data as an
#'   `sf` object instead.
#'
#' @examples
#' \dontrun{
#'   temp_model_dir_1 <- fs::path_temp("onesdm_model_dir")
#'   easin_data <- prepare_easin_data(
#'     easin_ids = "R00042",
#'     model_dir = temp_model_dir,
#'     return_data = TRUE
#'   )
#'   print(easin_data)
#'
#'   temp_model_dir_2 <- fs::path_temp("onesdm_model_dir")
#'   options(
#'     onesdm_easin_ids = c("R00042", "R00544"),
#'     onesdm_easin_exclude_gbif = TRUE,
#'     onesdm_easin_overwrite = TRUE,
#'     onesdm_model_dir = temp_model_dir_2
#'   )
#'   prepare_easin_data()
#' }
#'
#' @export
#' @author Ahmed El-Gabbas, Maryna Golivets
#' @references  EASIN Geospatial Web Service:
#'   <https://easin.jrc.ec.europa.eu/apixg/home/geoqueries/>

prepare_easin_data <- function(
  easin_ids = NULL,
  model_dir = NULL,
  timeout = 600L,
  n_search = 1000L,
  n_attempts = 10L,
  sleep_time = 5L,
  exclude_gbif = TRUE,
  verbose = TRUE,
  start_year = 1981L,
  overwrite = FALSE,
  return_data = FALSE
) {
  start_time <- lubridate::now(tzone = "CET")

  # Variables used in dplyr/data.table-like NSE expressions
  WKT <- Year <- EASINID <- Name <- Authorship <- SpeciesId <- NULL
  longitude <- latitude <- point_coords <- NULL
  n_dec_long <- n_dec_lat <- n_obs <- NULL
  matched <- taxon_fact_sheet <- matched_taxon <- NULL

  # ---------------------------------------------------------------------------
  # Packages
  # ---------------------------------------------------------------------------

  ecokit::check_packages(
    c(
      "cli",
      "CoordinateCleaner",
      "crayon",
      "dplyr",
      "fs",
      "httr",
      "jsonlite",
      "lubridate",
      "purrr",
      "RCurl",
      "sf",
      "stringr",
      "tibble",
      "tidyr",
      "tidyselect",
      "withr"
    )
  )

  # ---------------------------------------------------------------------------
  # Arguments from options
  # ---------------------------------------------------------------------------

  easin_ids <- ecokit::assign_from_options(
    easin_ids,
    "onesdm_easin_ids",
    "character"
  )

  model_dir <- ecokit::assign_from_options(
    model_dir,
    "onesdm_model_dir",
    "character"
  )

  timeout <- ecokit::assign_from_options(
    timeout,
    "onesdm_easin_timeout",
    c("numeric", "integer")
  )

  n_search <- ecokit::assign_from_options(
    n_search,
    "onesdm_easin_n_search",
    c("numeric", "integer")
  )

  n_attempts <- ecokit::assign_from_options(
    n_attempts,
    "onesdm_easin_n_attempts",
    c("numeric", "integer")
  )

  sleep_time <- ecokit::assign_from_options(
    sleep_time,
    "onesdm_easin_sleep_time",
    c("numeric", "integer")
  )

  exclude_gbif <- ecokit::assign_from_options(
    exclude_gbif,
    "onesdm_easin_exclude_gbif",
    "logical"
  )

  start_year <- ecokit::assign_from_options(
    start_year,
    "onesdm_start_year",
    c("numeric", "integer")
  )

  overwrite <- ecokit::assign_from_options(
    overwrite,
    "onesdm_easin_overwrite",
    "logical"
  )

  verbose <- ecokit::assign_from_options(
    verbose,
    "onesdm_easin_verbose",
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

  if (is.null(easin_ids)) {
    ecokit::stop_ctx(
      "easin_ids must contain at least one EASIN ID.",
      cat_timestamp = FALSE
    )
  }

  easin_ids <- easin_ids[!is.na(easin_ids)]

  if (length(easin_ids) == 0L) {
    ecokit::stop_ctx(
      "easin_ids must contain at least one EASIN ID.",
      cat_timestamp = FALSE
    )
  }

  if (!all(stringr::str_detect(easin_ids, "^R\\d{5}$"))) {
    ecokit::stop_ctx(
      paste0(
        "easin_ids must be in the format 'RXXXXX', where X is an integer.\n",
        "Please provide correct EASIN ID(s)."
      ),
      easin_ids = easin_ids,
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

  if (timeout <= 0) {
    ecokit::stop_ctx(
      "timeout must be greater than zero.",
      cat_timestamp = FALSE
    )
  }

  if (n_search <= 0) {
    ecokit::stop_ctx(
      "n_search must be greater than zero.",
      cat_timestamp = FALSE
    )
  }

  if (n_attempts <= 0) {
    ecokit::stop_ctx(
      "n_attempts must be greater than zero.",
      cat_timestamp = FALSE
    )
  }

  if (sleep_time < 0) {
    ecokit::stop_ctx(
      "sleep_time cannot be negative.",
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

  # ---------------------------------------------------------------------------
  # Match EASIN IDs to taxa
  # ---------------------------------------------------------------------------

  easin_ids <- unique(easin_ids)

  ecokit::cat_time(
    "Match EASIN IDs to taxon names and fact sheets",
    cat_timestamp = FALSE,
    cat_bold = TRUE,
    cat_red = TRUE,
    verbose = verbose
  )

  easin_taxon_url <- "https://easin.jrc.ec.europa.eu/apixg/catxg"

  matched_taxa <- purrr::map_dfr(
    easin_ids,
    function(easin_id) {
      url <- stringr::str_glue(
        "{easin_taxon_url}/easinid/{easin_id}"
      )

      taxon_data <- tryCatch(
        RCurl::getURL(
          url,
          .mapUnicode = FALSE,
          timeout = timeout
        ),
        error = function(e) NULL
      )

      # Request failed
      if (is.null(taxon_data)) {
        return(
          tibble::tibble(
            EASINID = easin_id,
            Name = NA_character_,
            matched_taxon = NA_character_,
            taxon_fact_sheet = NA_character_
          )
        )
      }

      # EASIN returned no match
      if (
        stringr::str_detect(
          taxon_data,
          "There are no results"
        )
      ) {
        return(
          tibble::tibble(
            EASINID = easin_id,
            Name = NA_character_,
            matched_taxon = NA_character_,
            taxon_fact_sheet = NA_character_
          )
        )
      }

      taxon_data <- tryCatch(
        jsonlite::fromJSON(
          taxon_data,
          flatten = TRUE
        ) |>
          tibble::as_tibble(),
        error = function(e) NULL
      )

      if (is.null(taxon_data) || nrow(taxon_data) == 0L) {
        return(
          tibble::tibble(
            EASINID = easin_id,
            Name = NA_character_,
            matched_taxon = NA_character_,
            taxon_fact_sheet = NA_character_
          )
        )
      }

      # Make sure expected fields exist
      if (!"Name" %in% names(taxon_data)) {
        taxon_data$Name <- NA_character_
      }

      if (!"Authorship" %in% names(taxon_data)) {
        taxon_data$Authorship <- NA_character_
      }

      taxon_data <- taxon_data |>
        dplyr::mutate(
          EASINID = easin_id,
          matched_taxon = dplyr::if_else(
            is.na(Name),
            NA_character_,
            stringr::str_squish(
              paste(Name, Authorship)
            )
          )
        ) |>
        dplyr::select(
          tidyselect::any_of(
            c("Name", "EASINID", "matched_taxon")
          )
        )

      fact_sheet_url <- stringr::str_glue(
        "https://easin.jrc.ec.europa.eu/",
        "spexplorer/species/factsheet/{easin_id}"
      )

      fact_sheet_ok <- tryCatch(
        !httr::http_error(fact_sheet_url) &&
          ecokit::check_url(fact_sheet_url),
        error = function(e) FALSE
      )

      taxon_data |>
        dplyr::mutate(
          taxon_fact_sheet = if (fact_sheet_ok) {
            fact_sheet_url
          } else {
            NA_character_
          }
        )
    }
  ) |>
    dplyr::mutate(
      matched = paste0(
        EASINID,
        dplyr::if_else(
          is.na(matched_taxon),
          ": not matched to EASIN database",
          paste0(": ", matched_taxon)
        ),
        dplyr::if_else(
          is.na(taxon_fact_sheet),
          "",
          paste0(
            " (",
            cli::style_hyperlink(
              text = crayon::blue("fact sheet"),
              url = taxon_fact_sheet
            ),
            ")"
          )
        )
      )
    )

  # ---------------------------------------------------------------------------
  # Stop if no IDs could be matched
  # ---------------------------------------------------------------------------

  if (all(is.na(matched_taxa$matched_taxon))) {
    ecokit::cat_time(
      paste0(
        "None of the provided EASIN IDs could be matched to EASIN database.\n",
        "No EASIN data will be downloaded."
      ),
      cat_timestamp = FALSE,
      level = 1L,
      verbose = verbose
    )

    return(invisible(NULL))
  }

  # Report skipped IDs
  if (anyNA(matched_taxa$matched_taxon)) {
    skipped_ids <- matched_taxa |>
      dplyr::filter(is.na(matched_taxon)) |>
      dplyr::pull(EASINID) |>
      toString()

    ecokit::cat_time(
      paste0(
        "Some EASIN IDs could not be matched to EASIN database.\n",
        "  >>>  These EASIN IDs will be skipped: ",
        crayon::red(skipped_ids)
      ),
      cat_timestamp = FALSE,
      level = 1L,
      verbose = verbose
    )
  }

  matched_taxa |>
    dplyr::filter(!is.na(matched_taxon)) |>
    dplyr::pull(matched) |>
    paste(collapse = "\n  >>  ") |>
    ecokit::cat_time(
      cat_timestamp = FALSE,
      level = 1L,
      verbose = verbose
    )

  # Only process successfully matched IDs
  easin_ids_matched <- matched_taxa |>
    dplyr::filter(!is.na(matched_taxon)) |>
    dplyr::pull(EASINID)

  # ---------------------------------------------------------------------------
  # Paths
  # ---------------------------------------------------------------------------

  path_data <- fs::path(model_dir, "data")
  path_easin_data <- fs::path(path_data, "easin_data.RData")
  path_easin_data_raw <- fs::path(path_data, "easin_data_raw.RData")

  # ---------------------------------------------------------------------------
  # Existing processed data
  # ---------------------------------------------------------------------------

  if (ecokit::check_data(path_easin_data, warning = FALSE)) {
    if (!overwrite) {
      ecokit::cat_time(
        paste0(
          "\nEASIN data already exist at: ",
          crayon::blue(path_easin_data)
        ),
        cat_timestamp = FALSE,
        verbose = verbose
      )

      ecokit::cat_time(
        "  >>>  Use overwrite = TRUE to re-download and re-process the data.",
        cat_timestamp = FALSE,
        verbose = verbose
      )

      if (return_data) {
        return(
          invisible(
            ecokit::load_as(path_easin_data)
          )
        )
      }

      return(invisible(path_easin_data))
    }

    ecokit::cat_time(
      crayon::blue(
        "\nEASIN data already exist and will be overwritten."
      ),
      cat_timestamp = FALSE,
      cat_bold = TRUE,
      verbose = verbose
    )
  }

  fs::dir_create(path_data)

  # ---------------------------------------------------------------------------
  # Print parameters
  # ---------------------------------------------------------------------------

  if (verbose) {
    ecokit::cat_time(
      "\nEASIN data extraction parameters",
      cat_timestamp = FALSE,
      cat_bold = TRUE,
      cat_red = TRUE
    )

    parameters <- c(
      `EASIN ID(s)` = toString(easin_ids_matched),
      `Modelling directory` = model_dir,
      `Modelling directory (absolute)` = fs::path_abs(model_dir),
      `Timeout` = timeout,
      `Number of search results per chunk` = ecokit::format_number(n_search),
      `Number of download attempts` = ecokit::format_number(n_attempts),
      `Sleep time` = sleep_time,
      `Exclude GBIF` = exclude_gbif,
      `Start year` = start_year,
      `Overwrite` = overwrite,
      `Return processed data` = return_data
    )

    purrr::iwalk(
      parameters,
      ~ ecokit::cat_time(
        paste0(
          crayon::italic(paste0(.y, ": ")),
          crayon::blue(.x)
        ),
        level = 1L,
        cat_timestamp = FALSE
      )
    )

    ecokit::cat_time(
      "\nExtracting EASIN data",
      cat_timestamp = FALSE,
      cat_bold = TRUE,
      cat_red = TRUE
    )
  }

  # ---------------------------------------------------------------------------
  # Download settings
  # ---------------------------------------------------------------------------

  withr::local_options(
    list(
      scipen = 999L,
      timeout = timeout
    )
  )

  # ---------------------------------------------------------------------------
  # Download raw EASIN data
  # ---------------------------------------------------------------------------

  if (!ecokit::check_data(path_easin_data_raw, warning = FALSE) || overwrite) {
    easin_data <- purrr::map_dfr(
      easin_ids_matched,
      get_easin_internal,
      path_data = path_data,
      timeout = timeout,
      n_search = n_search,
      n_attempts = n_attempts,
      sleep_time = sleep_time,
      exclude_gbif = exclude_gbif,
      overwrite = overwrite,
      verbose = verbose
    )

    if (nrow(easin_data) == 0L) {
      ecokit::cat_time(
        "\nNo EASIN data were downloaded for the provided EASIN IDs.",
        cat_timestamp = FALSE,
        verbose = verbose
      )

      if (return_data) {
        return(tibble::tibble())
      }

      return(invisible(NA_character_))
    }

    ecokit::cat_time(
      paste0(
        "Saving raw EASIN data to: `",
        crayon::blue(path_easin_data_raw),
        "`"
      ),
      cat_timestamp = FALSE,
      verbose = verbose,
      level = 1L
    )

    ecokit::save_as(
      object = easin_data,
      object_name = "easin_data_raw",
      out_path = path_easin_data_raw
    )
  } else {
    easin_data <- ecokit::load_as(path_easin_data_raw)
  }

  # ---------------------------------------------------------------------------
  # Clean EASIN data
  # ---------------------------------------------------------------------------

  ecokit::cat_time(
    "\nCleaning EASIN data",
    cat_timestamp = FALSE,
    cat_bold = TRUE,
    cat_red = TRUE,
    verbose = verbose
  )

  n_rows_raw <- nrow(easin_data)

  # ---------------------------------------------------------------------------
  # Remove missing coordinates
  # ---------------------------------------------------------------------------

  ecokit::cat_time(
    "Discarding entries without coordinates",
    cat_timestamp = FALSE,
    verbose = verbose
  )

  easin_data <- easin_data |>
    dplyr::filter(!is.na(WKT))

  n_rows_1 <- nrow(easin_data)

  ecokit::cat_time(
    paste0(
      "Filtered out ",
      ecokit::format_number(n_rows_raw - n_rows_1),
      " records without coordinates"
    ),
    cat_timestamp = FALSE,
    level = 1L,
    verbose = verbose
  )

  if (n_rows_1 == 0L) {
    ecokit::cat_time(
      "No EASIN records remain after filtering",
      cat_timestamp = FALSE,
      level = 1L,
      verbose = verbose
    )

    if (return_data) {
      return(tibble::tibble())
    }

    return(invisible(NA_character_))
  }

  # ---------------------------------------------------------------------------
  # Filter by year
  # ---------------------------------------------------------------------------

  ecokit::cat_time(
    paste0(
      "Discarding records older than ",
      start_year
    ),
    cat_timestamp = FALSE,
    verbose = verbose
  )

  easin_data <- easin_data |>
    dplyr::mutate(
      Year = suppressWarnings(as.integer(Year))
    ) |>
    dplyr::filter(
      !is.na(Year),
      Year >= start_year
    )

  n_rows_2 <- nrow(easin_data)

  if (n_rows_2 == 0L) {
    ecokit::cat_time(
      "No EASIN records remain after filtering",
      cat_timestamp = FALSE,
      level = 1L,
      verbose = verbose
    )

    if (return_data) {
      return(tibble::tibble())
    }

    return(invisible(NA_character_))
  }

  ecokit::cat_time(
    paste0(
      "Filtered out ",
      ecokit::format_number(n_rows_1 - n_rows_2),
      " records older than ",
      start_year,
      "."
    ),
    cat_timestamp = FALSE,
    level = 1L,
    verbose = verbose
  )

  # ---------------------------------------------------------------------------
  # Extract coordinates from WKT
  # ---------------------------------------------------------------------------

  ecokit::cat_time(
    "Extracting coordinates from WKT strings",
    cat_timestamp = FALSE,
    verbose = verbose
  )

  easin_data <- easin_data |>
    dplyr::mutate(
      point_coords = purrr::map(
        WKT,
        function(wkt) {
          point <- tryCatch(
            sf::st_as_sfc(wkt, crs = 4326),
            error = function(e) NULL
          )

          if (
            is.null(point) ||
              length(point) == 0L ||
              !inherits(point[[1L]], "POINT")
          ) {
            return(
              tibble::tibble(
                longitude = NA_real_,
                latitude = NA_real_
              )
            )
          }

          coords <- sf::st_coordinates(point)

          tibble::tibble(
            longitude = coords[1L, "X"],
            latitude = coords[1L, "Y"]
          )
        }
      )
    ) |>
    tidyr::unnest(point_coords)

  n_rows_3 <- nrow(easin_data)

  if (n_rows_3 == 0L) {
    ecokit::cat_time(
      "No EASIN records remain after coordinate extraction",
      cat_timestamp = FALSE,
      level = 1L,
      verbose = verbose
    )

    if (return_data) {
      return(tibble::tibble())
    }

    return(invisible(NA_character_))
  }

  # ---------------------------------------------------------------------------
  # Coordinate validity
  # ---------------------------------------------------------------------------

  ecokit::cat_time(
    "Checking coordinate validity",
    cat_timestamp = FALSE,
    verbose = verbose
  )

  easin_data <- easin_data |>
    dplyr::filter(
      !is.na(longitude),
      !is.na(latitude)
    ) |>
    CoordinateCleaner::cc_val(
      lon = "longitude",
      lat = "latitude",
      verbose = FALSE
    )

  n_rows_4 <- nrow(easin_data)

  if (n_rows_4 == 0L) {
    ecokit::cat_time(
      "No EASIN records remain after coordinate validation",
      cat_timestamp = FALSE,
      level = 1L,
      verbose = verbose
    )

    if (return_data) {
      return(tibble::tibble())
    }

    return(invisible(NA_character_))
  }

  # ---------------------------------------------------------------------------
  # Coordinate precision
  # ---------------------------------------------------------------------------

  ecokit::cat_time(
    "Excluding observations with low coordinate precision",
    cat_timestamp = FALSE,
    verbose = verbose
  )

  easin_data <- easin_data |>
    dplyr::mutate(
      n_dec_long = ecokit::n_decimals(longitude),
      n_dec_lat = ecokit::n_decimals(latitude)
    ) |>
    dplyr::filter(
      n_dec_long > 1L,
      n_dec_lat > 1L
    )

  n_rows_5 <- nrow(easin_data)

  if (n_rows_5 == 0L) {
    ecokit::cat_time(
      "No EASIN records remain after precision filtering",
      cat_timestamp = FALSE,
      level = 1L,
      verbose = verbose
    )

    if (return_data) {
      return(tibble::tibble())
    }

    return(invisible(NA_character_))
  }

  ecokit::cat_time(
    paste0(
      "Filtered out ",
      ecokit::format_number(n_rows_4 - n_rows_5),
      " records with low spatial precision"
    ),
    cat_timestamp = FALSE,
    level = 1L,
    verbose = verbose
  )

  # ---------------------------------------------------------------------------
  # CoordinateCleaner
  # ---------------------------------------------------------------------------

  ecokit::cat_time(
    "Cleaning using `CoordinateCleaner`",
    cat_timestamp = FALSE,
    verbose = verbose
  )

  # Explicitly use a 100-m buffer for the institution test.
  # `geod = TRUE` makes the buffer unit metres.
  easin_data <- easin_data |>
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
      .,
      buffer = 100,
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

  n_rows_6 <- nrow(easin_data)

  if (n_rows_6 == 0L) {
    ecokit::cat_time(
      "No EASIN records remain after CoordinateCleaner",
      cat_timestamp = FALSE,
      level = 1L,
      verbose = verbose
    )

    if (return_data) {
      return(tibble::tibble())
    }

    return(invisible(NA_character_))
  }

  ecokit::cat_time(
    paste0(
      "Filtered out ",
      ecokit::format_number(n_rows_5 - n_rows_6),
      " records using `CoordinateCleaner`"
    ),
    cat_timestamp = FALSE,
    level = 1L,
    verbose = verbose
  )

  # ---------------------------------------------------------------------------
  # Convert to sf
  # ---------------------------------------------------------------------------

  easin_data <- easin_data |>
    sf::st_as_sf(
      coords = c("longitude", "latitude"),
      crs = 4326,
      remove = FALSE
    )

  # ---------------------------------------------------------------------------
  # Summary
  # ---------------------------------------------------------------------------

  n_obs_per_id <- sf::st_drop_geometry(easin_data) |>
    dplyr::count(SpeciesId, name = "n_obs") |>
    dplyr::mutate(
      n_obs = paste0(
        SpeciesId,
        " (",
        format(n_obs, big.mark = ",", scientific = FALSE),
        ")"
      )
    ) |>
    dplyr::pull(n_obs) |>
    toString()

  ecokit::cat_time(
    paste0(
      "A total of ",
      ecokit::format_number(nrow(easin_data)),
      " filtered observations were extracted for EASIN ID(s): ",
      crayon::blue(n_obs_per_id)
    ),
    cat_timestamp = FALSE,
    verbose = verbose
  )

  # ---------------------------------------------------------------------------
  # Save
  # ---------------------------------------------------------------------------

  ecokit::cat_time(
    paste0(
      "Saving EASIN data to: `",
      crayon::blue(path_easin_data),
      "`"
    ),
    cat_timestamp = FALSE,
    verbose = verbose
  )

  save(
    easin_data,
    file = path_easin_data
  )

  # ---------------------------------------------------------------------------
  # Finish
  # ---------------------------------------------------------------------------

  ecokit::cat_diff(
    init_time = start_time,
    prefix = "\nExtracting EASIN data was finished in ",
    verbose = verbose
  )

  if (return_data) {
    return(easin_data)
  }

  invisible(path_easin_data)
}


# ---------------------------------------------------------------------------
# Internal function to download EASIN chunk data
# ---------------------------------------------------------------------------

#' @noRd
#' @keywords internal

get_easin_internal <- function(
  easin_id = NULL,
  path_data = NULL,
  timeout = 300L,
  n_search = 1000L,
  n_attempts = 10L,
  sleep_time = 5L,
  exclude_gbif = TRUE,
  overwrite = FALSE,
  verbose = TRUE
) {
  SpeciesId <- NULL

  ecokit::check_args(
    args_to_check = "easin_id",
    args_type = "character"
  )

  easin_url <- "https://easin.jrc.ec.europa.eu/apixg/geoxg"

  easin_data_sub <- list()
  chunk_n <- 0L

  ecokit::cat_time(
    paste0(
      "Processing EASIN ID: ",
      crayon::blue(easin_id)
    ),
    cat_timestamp = FALSE,
    level = 1L,
    verbose = verbose
  )

  repeat {
    chunk_n <- chunk_n + 1L
    skip <- (chunk_n - 1L) * n_search

    url <- if (exclude_gbif) {
      stringr::str_glue(
        "{easin_url}/{easin_id}/exclude/dps/1/{skip}/{n_search}"
      )
    } else {
      stringr::str_glue(
        "{easin_url}/{easin_id}/{skip}/{n_search}"
      )
    }

    chunk_name <- paste0(
      "easin_",
      easin_id,
      "_chunk_",
      chunk_n
    )

    chunk_file <- fs::path(
      path_data,
      paste0(chunk_name, ".RData")
    )

    chunk_data <- NULL
    success <- FALSE
    no_results <- FALSE

    # -----------------------------------------------------------------------
    # Retry download
    # -----------------------------------------------------------------------

    for (attempt in seq_len(n_attempts)) {
      ecokit::cat_time(
        paste0(
          cli::style_hyperlink(
            text = crayon::blue(
              paste0("chunk ", chunk_n)
            ),
            url = url
          ),
          " (attempt ",
          attempt,
          "/",
          n_attempts,
          ")"
        ),
        level = 2L,
        cat_timestamp = FALSE,
        verbose = verbose
      )

      # Use cached chunk only when overwrite = FALSE
      if (
        !overwrite &&
          ecokit::check_data(chunk_file, warning = FALSE)
      ) {
        ecokit::cat_time(
          "Loading cached chunk data from disk",
          level = 3L,
          cat_timestamp = FALSE,
          verbose = verbose
        )

        chunk_data <- ecokit::load_as(chunk_file)
        success <- TRUE
        break
      }

      # Download chunk
      response <- tryCatch(
        RCurl::getURL(
          url,
          .mapUnicode = FALSE,
          timeout = timeout
        ),
        error = function(e) NULL
      )

      if (is.null(response)) {
        if (attempt < n_attempts) {
          Sys.sleep(sleep_time)
        }

        next
      }

      # No more observations
      if (
        stringr::str_detect(
          response,
          "There are no results based on your"
        )
      ) {
        no_results <- TRUE
        chunk_data <- tibble::tibble()
        success <- TRUE

        break
      }

      # Parse JSON
      chunk_data <- tryCatch(
        jsonlite::fromJSON(
          response,
          flatten = TRUE
        ) |>
          tibble::as_tibble() |>
          dplyr::mutate(
            json_url = url
          ),
        error = function(e) NULL
      )

      if (!is.null(chunk_data)) {
        success <- TRUE
        break
      }

      if (attempt < n_attempts) {
        Sys.sleep(sleep_time)
      }
    }

    # -----------------------------------------------------------------------
    # Handle failed chunk
    # -----------------------------------------------------------------------

    if (!success) {
      ecokit::cat_time(
        paste0(
          "Failed to download ",
          "chunk ",
          chunk_n,
          " after ",
          n_attempts,
          " attempts."
        ),
        level = 1L,
        cat_timestamp = FALSE,
        verbose = verbose
      )

      break
    }

    # No results means we have reached the end
    if (no_results) {
      break
    }

    # Save successfully downloaded chunk
    if (inherits(chunk_data, "data.frame")) {
      easin_data_sub[[length(easin_data_sub) + 1L]] <- chunk_data

      # Only save non-empty downloaded chunks
      if (nrow(chunk_data) > 0L) {
        ecokit::save_as(
          object = chunk_data,
          object_name = chunk_name,
          out_path = chunk_file
        )
      }
    }

    # A short chunk means this was the final page
    if (
      !inherits(chunk_data, "data.frame") ||
        nrow(chunk_data) < n_search
    ) {
      break
    }

    Sys.sleep(sleep_time)
  }

  easin_data_sub <- dplyr::bind_rows(
    easin_data_sub
  ) |>
    ecokit::add_missing_columns(
      fill_value = easin_id,
      SpeciesId
    )

  ecokit::cat_time(
    paste0(
      "A total of ",
      ecokit::format_number(
        nrow(easin_data_sub)
      ),
      " observations were extracted for EASIN ID: ",
      crayon::blue(easin_id)
    ),
    cat_timestamp = FALSE,
    level = 2L,
    verbose = verbose
  )

  easin_data_sub
}
