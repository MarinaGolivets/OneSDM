#' @title Download and Clean EASIN Occurrence Data for Given Species IDs
#'
#' @description Downloads, combines, cleans, and saves occurrence data for given
#'   EASIN species IDs from the [EASIN](https://alien.jrc.ec.europa.eu/easin)
#'   (European Alien Species Information Network) API.
#'
#' @param easin_ids Character vector. One or more EASIN species IDs.
#'   Each ID should start with 'R' followed by five digits (e.g. "R00544").
#'   This cannot be `NULL`. Species IDs can be found on the [EASIN
#'   website](https://easin.jrc.ec.europa.eu/apixg/home/geoqueries/) by
#'   searching for a species and locating the "EASIN ID" in the species
#'   details section. Can also be set via the `onestop_easin_ids` option.
#' @param model_dir Character. Path to the modelling directory where data and
#'   fitted models will be saved. This cannot be `NULL` and must the same
#'   directory used for the same species data. This can also be set via the
#'   `onestop_model_dir` option.
#' @param timeout Integer. Timeout (in seconds) for each download attempt.
#'   Default is 300L. Can also be set via the `onestop_easin_timeout` option.
#'   Note: This timeout only applies to the download of each chunk, not the
#'   entire download process.
#' @param n_search Integer. Number of records to request per API call (chunk
#'   size). Default is 1000L, which is the maximum allowed by the EASIN API. Can
#'   also be set via the `onestop_easin_n_search` option.
#' @param n_attempts Integer. Maximum number of download attempts per chunk.
#'   Default is 10L. Can also be set via the `onestop_easin_n_attempts` option.
#' @param sleep_time Integer. Seconds to wait between chunk downloads. This
#'   helps avoid overwhelming the EASIN server. Default is 5L. Can also be set
#'   via the `onestop_easin_sleep_time` option.
#' @param exclude_gbif Logical. If `TRUE` (default), exclude
#'   [GBIF](https://www.gbif.org/) — the Global Biodiversity Information
#'   Facility — records from the download. Can also be set via the
#'   `onestop_easin_exclude_gbif` option.
#' @param verbose Logical. If `TRUE` (default), print progress and information
#'   messages, including the URL of the currently processed chunk. Can also be
#'   set via the `onestop_easin_verbose` option.
#' @param start_year Integer. Include only records from this year onward.
#'   The default is 1981L. Can also be set via the `onestop_start_year` option.
#'
#' @details
#' - The function supports chunked downloads, multiple attempts, and
#' extensive data cleaning, including coordinate precision checks and spatial
#' filtering using the `CoordinateCleaner` package. The cleaned data is saved as
#' an `.RData` file in the specified model directory.
#' - The function checks for an existing cleaned EASIN data file and skips
#' download if found.
#' - The function applies several cleaning steps; e.g. filtering out records
#' with low coordinate precision, equal longitude / latitude, or near centroids,
#' capitals, biodiversity institutions, and GBIF headquarters.
#'
#' @return (Invisibly) the path to the saved `.RData` file containing the
#'   cleaned EASIN data (`data/easin_data.RData`) in the specified `model_dir`.
#'
#' @examples
#' \dontrun{
#'   process_easin_data(easin_ids = "R00042", model_dir = "path/to/model_dir")
#' }
#'
#' @export
#' @author Ahmed El-Gabbas
#' @references  EASIN geospatial Web service:
#'   <https://easin.jrc.ec.europa.eu/apixg/home/geoqueries/>

process_easin_data <- function(
  easin_ids = NULL,
  model_dir = NULL,
  timeout = 300L,
  n_search = 1000L,
  n_attempts = 10L,
  sleep_time = 5L,
  exclude_gbif = TRUE,
  verbose = TRUE,
  start_year = 1981L
) {
  .start_easin_time <- lubridate::now(tzone = "CET")

  longitude <- latitude <- WKT <- Year <- n_obs <- SpeciesId <- #nolint
    point_coords <- n_dec_long <- n_dec_lat <- NULL

  # # ********************************************************************** #
  # Assigning function arguments from options if not provided directly ------
  # # ********************************************************************** #

  easin_ids <- ecokit::assign_from_options(
    easin_ids,
    "onestop_easin_ids",
    "character"
  )
  model_dir <- ecokit::assign_from_options(
    model_dir,
    "onestop_model_dir",
    "character"
  )
  timeout <- ecokit::assign_from_options(
    timeout,
    "onestop_easin_timeout",
    c("numeric", "integer")
  )
  n_search <- ecokit::assign_from_options(
    n_search,
    "onestop_easin_n_search",
    c("numeric", "integer")
  )
  n_attempts <- ecokit::assign_from_options(
    n_attempts,
    "onestop_easin_n_attempts",
    c("numeric", "integer")
  )
  sleep_time <- ecokit::assign_from_options(
    sleep_time,
    "onestop_easin_sleep_time",
    c("numeric", "integer")
  )
  exclude_gbif <- ecokit::assign_from_options(
    exclude_gbif,
    "onestop_easin_exclude_gbif",
    "logical"
  )
  verbose <- ecokit::assign_from_options(
    verbose,
    "onestop_easin_verbose",
    "logical"
  )
  start_year <- ecokit::assign_from_options(
    start_year,
    "onestop_start_year",
    c("numeric", "integer")
  )

  # # ********************************************************************** #
  # Checking function arguments -------
  # # ********************************************************************** #

  # Check easin_ids
  if (is.null(easin_ids)) {
    ecokit::stop_ctx(
      "easin_ids can not be NULL. Please provide at least one EASIN ID.",
      cat_timestamp = FALSE
    )
  }
  if (!all(stringr::str_detect(easin_ids, "^R\\d{5}$"))) {
    #nolint
    ecokit::stop_ctx(
      "All easin_ids must be in the format 'RXXXXX' where X is an integer.",
      easin_ids = easin_ids,
      cat_timestamp = FALSE
    )
  }

  if (is.null(model_dir)) {
    ecokit::stop_ctx(
      paste0(
        "The model_dir argument must be provided either directly or via the ",
        "`onestop_model_dir` option."
      ),
      cat_timestamp = FALSE
    )
  }

  path_data <- fs::path(model_dir, "data")
  path_easin_data <- fs::path(path_data, "easin_data.RData")

  if (ecokit::check_data(path_easin_data, warning = FALSE)) {
    ecokit::cat_time(
      paste0(
        "EASIN data file already exists at: ",
        crayon::blue(path_easin_data)
      ),
      cat_timestamp = FALSE,
      verbose = verbose
    )
    return(invisible(path_easin_data))
  }

  fs::dir_create(path_data)

  # # ********************************************************************** #
  # Print function arguments ------
  # # ********************************************************************** #

  if (verbose) {
    ecokit::cat_time(
      "EASIN data extraction parameters:",
      cat_timestamp = FALSE,
      cat_bold = TRUE,
      cat_red = TRUE
    )
    ecokit::cat_time(
      paste0(crayon::bold("EASIN IDs: "), toString(easin_ids)),
      level = 1L,
      cat_timestamp = FALSE
    )
    ecokit::cat_time(
      paste0(crayon::bold("Modelling directory: "), model_dir),
      level = 1L,
      cat_timestamp = FALSE
    )
    ecokit::cat_time(
      paste0(crayon::bold("Timeout: "), timeout),
      level = 1L,
      cat_timestamp = FALSE
    )
    ecokit::cat_time(
      paste0(crayon::bold("Number of search results: "), n_search),
      level = 1L,
      cat_timestamp = FALSE
    )
    ecokit::cat_time(
      paste0(crayon::bold("Number of attempts: "), n_attempts),
      level = 1L,
      cat_timestamp = FALSE
    )
    ecokit::cat_time(
      paste0(crayon::bold("Sleep time: "), sleep_time),
      level = 1L,
      cat_timestamp = FALSE
    )
    ecokit::cat_time(
      paste0(crayon::bold("Exclude GBIF: "), exclude_gbif),
      level = 1L,
      cat_timestamp = FALSE
    )
    ecokit::cat_time(
      paste0(crayon::bold("Start year: "), start_year),
      level = 1L,
      cat_timestamp = FALSE
    )
    ecokit::cat_time(
      "\nExtracting EASIN data",
      cat_timestamp = TRUE,
      cat_bold = TRUE,
      cat_red = TRUE
    )
  }

  # Temporarily set download time out only within the function
  withr::local_options(list(scipen = 999L, timeout = timeout))

  # # ********************************************************************** #
  # Downloading EASIN data -----
  # # ********************************************************************** #

  easin_data <- purrr::map(
    .x = easin_ids,
    .f = get_easin_internal,
    timeout = timeout,
    n_search = n_search,
    n_attempts = n_attempts,
    sleep_time = sleep_time,
    exclude_gbif = exclude_gbif,
    verbose = verbose
  ) %>%
    dplyr::bind_rows()

  # # ********************************************************************** #
  # Cleaning EASIN data -----
  # # ********************************************************************** #

  # # ||||||||||||||||||||||||||||||||||||||||| #
  # Initial filtering of EASIN data
  # # ||||||||||||||||||||||||||||||||||||||||| #

  if (nrow(easin_data) > 0L) {
    easin_data <- easin_data %>%
      dplyr::mutate(Year = as.integer(Year)) %>%
      dplyr::filter(!is.na(WKT), Year >= start_year)
  }

  # # ||||||||||||||||||||||||||||||||||||||||| #
  # Extracting coordinates from WKT strings and filtering observations without
  # coordinates
  # # ||||||||||||||||||||||||||||||||||||||||| #

  if (nrow(easin_data) > 0L) {
    easin_data <- easin_data %>%
      dplyr::mutate(
        point_coords = purrr::map(
          .x = WKT,
          .f = ~ {
            # Extract POINT coordinates from WKT string
            point_coords_0 <- stringr::str_extract_all(
              .x,
              "POINT\\s*\\(\\s*-?\\d+\\.\\d+\\s+-?\\d+\\.\\d+\\s*\\)"
            )[[1L]]

            if (length(point_coords_0) > 0L) {
              purrr::map(
                .x = point_coords_0,
                .f = ecokit::text_to_coordinates,
                name_x = "longitude",
                name_y = "latitude"
              ) %>%
                dplyr::bind_rows() %>%
                dplyr::mutate(
                  dplyr::across(
                    .cols = c("longitude", "latitude"),
                    .fns = ~ round(.x, 5L)
                  )
                )
            } else {
              tibble::tibble(longitude = NA_real_, latitude = NA_real_)
            }
          }
        )
      ) %>%
      tidyr::unnest(point_coords) %>%
      dplyr::filter(!is.na(longitude) & !is.na(latitude))
  }

  # # ||||||||||||||||||||||||||||||||||||||||| #
  # Exclude observations with low coordinate precision
  # # ||||||||||||||||||||||||||||||||||||||||| #

  if (nrow(easin_data) > 0L) {
    easin_data <- easin_data %>%
      dplyr::mutate(
        # number of decimal places for longitude / latitude
        n_dec_long = ecokit::n_decimals(longitude),
        n_dec_lat = ecokit::n_decimals(latitude)
      ) %>%
      dplyr::filter(
        # exclude occurrences if either latitude/longitude has 1 or 0 decimals
        (n_dec_long > 1L & n_dec_lat > 1L),
        # exclude equal coordinates
        longitude != latitude
      )
  }

  # # ||||||||||||||||||||||||||||||||||||||||| #
  # Further cleaning of EASIN data using CoordinateCleaner package
  # # ||||||||||||||||||||||||||||||||||||||||| #

  if (nrow(easin_data) > 0L) {
    easin_data <- easin_data %>%
      # Clean coordinates further using CoordinateCleaner package
      # exclude coordinates in vicinity of country and province centroids
      CoordinateCleaner::cc_cen(
        buffer = 100L,
        lon = "longitude",
        lat = "latitude",
        verbose = FALSE
      ) %>%
      # exclude coordinates in vicinity of country capitals
      CoordinateCleaner::cc_cap(
        buffer = 100L,
        lon = "longitude",
        lat = "latitude",
        verbose = FALSE
      ) %>%
      # exclude records in the vicinity of biodiversity institutions
      CoordinateCleaner::cc_inst(
        lon = "longitude",
        lat = "latitude",
        verbose = FALSE
      ) %>%
      # Identify Records Assigned to GBIF Headquarters
      CoordinateCleaner::cc_gbif(
        buffer = 100L,
        lon = "longitude",
        lat = "latitude",
        verbose = FALSE
      ) %>%
      CoordinateCleaner::cc_equ(
        lon = "longitude",
        lat = "latitude",
        verbose = FALSE
      )
  }

  # # ||||||||||||||||||||||||||||||||||||||||| #
  # Converting to sf object
  # # ||||||||||||||||||||||||||||||||||||||||| #

  if (nrow(easin_data) > 0L) {
    easin_data <- easin_data %>%
      # convert to sf object, while keeping original coordinates as columns
      sf::st_as_sf(
        coords = c("longitude", "latitude"),
        crs = 4326L,
        remove = FALSE
      )
  }

  # # ||||||||||||||||||||||||||||||||||||||||| #
  # Saving EASIN data
  # # ||||||||||||||||||||||||||||||||||||||||| #

  ecokit::cat_time(
    paste0("Saving EASIN data to ", crayon::blue(path_easin_data)),
    cat_timestamp = FALSE,
    level = 1L,
    verbose = verbose
  )
  save(easin_data, file = path_easin_data)

  # # ||||||||||||||||||||||||||||||||||||||||| #
  # Printing summary of extracted EASIN data
  # # ||||||||||||||||||||||||||||||||||||||||| #

  n_rows <- nrow(easin_data)
  if (n_rows > 0L) {
    n_obs_per_id <- sf::st_drop_geometry(easin_data) %>%
      dplyr::group_by(SpeciesId) %>%
      dplyr::tally(name = "n_obs") %>%
      dplyr::mutate(n_obs = paste0(SpeciesId, " (", n_obs, ")")) %>%
      dplyr::pull(n_obs) %>%
      toString()

    ecokit::cat_time(
      paste0(
        "A total of ",
        n_rows,
        " filtered observations were extracted for EASIN IDs: ",
        n_obs_per_id
      ),
      cat_timestamp = FALSE,
      level = 1L,
      verbose = verbose
    )
  } else {
    ecokit::cat_time(
      "No EASIN data were extracted",
      cat_timestamp = FALSE,
      level = 1L,
      verbose = verbose
    )
  }

  ecokit::cat_diff(
    init_time = .start_easin_time,
    prefix = "\nExtracting EASIN data was finished in ",
    verbose = verbose
  )

  # # ********************************************************************** #
  # # ********************************************************************** #

  return(invisible(path_easin_data))
}


#' @noRd
#' @keywords internal

get_easin_internal <- function(
  easin_id = NULL,
  timeout = 300L,
  n_search = 1000L,
  n_attempts = 10L,
  sleep_time = 5L,
  exclude_gbif = TRUE,
  verbose = TRUE
) {
  SpeciesId <- NULL #nolint
  ecokit::check_args(args_to_check = "easin_id", args_type = "character")

  easin_url <- "https://easin.jrc.ec.europa.eu/apixg/geoxg"
  easin_data_sub <- list()
  chunk_n <- 0L

  ecokit::cat_time(
    paste0("Processing EASIN ID: ", crayon::blue(easin_id)),
    cat_timestamp = FALSE,
    level = 1L,
    verbose = verbose
  )

  repeat {
    download_try <- 0L
    chunk_n <- chunk_n + 1L
    skip <- (chunk_n - 1L) * n_search

    if (exclude_gbif) {
      # `exclude/dps/1` to excludes GBIF observations
      url <- stringr::str_glue(
        "{easin_url}/{easin_id}/exclude/dps/1/{skip}/{n_search}"
      ) #nolint
    } else {
      url <- stringr::str_glue("{easin_url}/{easin_id}/{skip}/{n_search}") #nolint
    }

    while (download_try < n_attempts) {
      download_try <- download_try + 1L

      ecokit::cat_time(
        paste0(
          cli::style_hyperlink(
            text = crayon::blue(paste0("chunk ", chunk_n)),
            url = url
          ),
          " (attempt ",
          download_try,
          ")"
        ),
        level = 2L,
        cat_timestamp = FALSE,
        verbose = verbose
      )

      chunk_data <- try(
        RCurl::getURL(url, .mapUnicode = FALSE, timeout = timeout),
        silent = TRUE
      )
      if (inherits(chunk_data, "try-error")) {
        easin_data_sub[[chunk_n]] <- tibble::tibble()
        break
      }

      no_obs <- stringr::str_detect(
        chunk_data,
        pattern = "There are no results based on your"
      )

      if (no_obs) {
        easin_data_sub[[chunk_n]] <- tibble::tibble()
        break
      }

      chunk_data <- jsonlite::fromJSON(chunk_data, flatten = TRUE) %>%
        tibble::tibble() %>%
        dplyr::mutate(json_url = url)

      if (inherits(chunk_data, "data.frame")) {
        break
      }
    }

    if (inherits(chunk_data, "data.frame")) {
      easin_data_sub[[chunk_n]] <- chunk_data
      if (nrow(chunk_data) < n_search) {
        break
      }
    } else {
      break
    }

    # sleep at each chunk download
    Sys.sleep(sleep_time)
  }

  rm(easin_url, skip, envir = environment())

  easin_data_sub <- dplyr::bind_rows(easin_data_sub) %>%
    ecokit::add_missing_columns(fill_value = easin_id, SpeciesId)

  ecokit::cat_time(
    paste0(
      "A total of ",
      nrow(easin_data_sub),
      " observations were extracted for EASIN ID: ",
      easin_id
    ),
    cat_timestamp = FALSE,
    level = 2L,
    verbose = verbose
  )

  easin_data_sub
}
