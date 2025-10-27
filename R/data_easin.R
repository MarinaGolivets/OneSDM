#' @title Download and Clean EASIN Occurrence Data for Given Species IDs
#'
#' @description Downloads, combines, cleans, and saves occurrence data for given
#'   EASIN species IDs from the [EASIN](https://alien.jrc.ec.europa.eu/easin)
#'   (European Alien Species Information Network) API.
#'
#' @param easin_ids Character vector. One or more EASIN species IDs. Each ID
#'   should start with 'R' followed by five digits (e.g. "R00544"). This cannot
#'   be `NULL`. Species IDs can be found on the [EASIN
#'   website](https://easin.jrc.ec.europa.eu/spexplorer/search/) by searching
#'   for a species and locating the "EASIN ID" in the species details section.
#'   Can also be set via the `onesdm_easin_ids` option.
#' @param model_dir Character. Path to the modelling directory where data and
#'   fitted models will be saved. This cannot be `NULL` and must the same
#'   directory used for the same species data. This can also be set via the
#'   `onesdm_model_dir` option.
#' @param timeout Integer. Timeout (in seconds) for each download attempt.
#'   Default is `600L`. Can also be set via the `onesdm_easin_timeout` option.
#'   Note: This timeout only applies to the download of each chunk, not the
#'   entire download process.
#' @param n_search Integer. Number of records to request per API call (chunk
#'   size). Default is `1000L`, which is the maximum allowed by the EASIN API.
#'   Can also be set via the `onesdm_easin_n_search` option.
#' @param n_attempts Integer. Maximum number of download attempts per chunk.
#'   Default is `10L`. Can also be set via the `onesdm_easin_n_attempts` option.
#' @param sleep_time Integer. Seconds to wait between chunk downloads. This
#'   helps avoid overwhelming the EASIN server. Default is 5L. Can also be set
#'   via the `onesdm_easin_sleep_time` option.
#' @param exclude_gbif Logical. If `TRUE` (default), exclude
#'   [GBIF](https://www.gbif.org/) — the Global Biodiversity Information
#'   Facility — records from the download. Can also be set via the
#'   `onesdm_easin_exclude_gbif` option.
#' @param verbose Logical. If `TRUE` (default), print progress and information
#'   messages, including the URL of the currently processed chunk. Can also be
#'   set via the `onesdm_easin_verbose` option.
#' @param start_year Integer. Include only records from this year onward. The
#'   default is `1981L`, to match the temporal coverage of CHELSA climate data.
#'   Can also be set via the `onesdm_start_year` option.
#' @param overwrite Logical. If `TRUE`, overwrite existing cleaned EASIN data
#'   file in the model directory. Default is `FALSE`. Can also be set via
#'   `onesdm_easin_overwrite` option.
#' @param return_data Logical. If `TRUE`, returns the processed EASIN data as an
#'   `sf` object in addition to saving it. Default is `FALSE`.
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
#' - Function default arguments can be set globally using the `options()`
#' function. Users can set these options at the start of their R session to
#' avoid repeatedly specifying them in function calls. The following options
#' correspond to the function arguments:
#'   - `onesdm_easin_ids`: Character vector of EASIN species IDs.
#'   - `onesdm_model_dir`: Character. Path to the modelling directory.
#'   - `onesdm_easin_timeout`: Integer. Timeout (in seconds) for each download
#' attempt.
#'   - `onesdm_easin_n_search`: Integer. Number of records to request per API
#' call.
#'   - `onesdm_easin_n_attempts`: Integer. Maximum number of download attempts
#' per chunk.
#'   - `onesdm_easin_sleep_time`: Integer. Seconds to wait between chunk
#' downloads.
#'   - `onesdm_easin_exclude_gbif`: Logical. Whether to exclude GBIF records.
#'   - `onesdm_easin_verbose`: Logical. Whether to print progress messages.
#'   - `onesdm_start_year`: Integer. Minimum year for records to include.
#'   - `onesdm_easin_overwrite`: Logical. Whether to overwrite existing cleaned
#' data.
#' - Example of setting options:
#'   ```r
#'   options(
#'     onesdm_easin_ids = c("R00042", "R00544"),
#'     onesdm_easin_exclude_gbif = TRUE,
#'     onesdm_easin_overwrite = TRUE,
#'     onesdm_model_dir = "path/to/model_dir"
#'   )
#'   ```
#'
#' @return Invisibly returns the file path to the saved cleaned EASIN data
#'   `easin_data.RData` file in the `data` subdirectory of `model_dir`, unless
#'   `return_data` is `TRUE`, in which case it returns the cleaned EASIN data as
#'   an `sf` object.
#'
#' @examples
#' \dontrun{
#'  require(fs)
#'
#'  # Prepare EASIN data for species with EASIN ID "R00042": Acacia karroo Hayne
#'  temp_model_dir <- fs::path_temp("onesdm_model_dir")
#'  easin_data <- prepare_easin_data(
#'    easin_ids = "R00042", model_dir = temp_model_dir, return_data = TRUE)
#'
#'  print(easin_data)
#'
#'  # # ||||||||||||||||||||||||||||||||||||||||||||||||||| #
#'
#'  # Prepare EASIN data, using argument values from `options`
#'  temp_model_dir_2 <- fs::path_temp("onesdm_model_dir")
#'  options(
#'     onesdm_easin_ids = c("R00042", "R00544"),
#'     onesdm_easin_exclude_gbif = TRUE,
#'     onesdm_easin_overwrite = TRUE,
#'     onesdm_model_dir = temp_model_dir_2)
#'
#'  prepare_easin_data()
#'
#' }
#'
#' @export
#' @author Ahmed El-Gabbas
#' @references  EASIN geospatial Web service:
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
  .start_easin_time <- lubridate::now(tzone = "CET")

  WKT <- Year <- EASINID <- Name <- SpeciesId <- longitude <- latitude <- #nolint
    point_coords <- n_dec_long <- n_dec_lat <- n_obs <- matched <-
    species_fact_sheet <- matched_species <- NULL

  ecokit::check_packages(
    c(
      "cli", "CoordinateCleaner", "crayon", "fs", "httr", "jsonlite",
      "lubridate", "purrr", "RCurl", "sf", "stringr", "tibble", "tidyr",
      "tidyselect", "withr"))

  # # ********************************************************************** #
  # Assigning function arguments from options if not provided directly ------
  # # ********************************************************************** #

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
  verbose <- ecokit::assign_from_options(
    verbose,
    "onesdm_easin_verbose",
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

  # # ********************************************************************** #
  # Checking function arguments -------
  # # ********************************************************************** #

  # Check easin_ids
  if (is.null(easin_ids)) {
    ecokit::stop_ctx(
      "easin_ids cannot be NULL. Please provide at least one EASIN ID.",
      cat_timestamp = FALSE
    )
  }
  if (!all(stringr::str_detect(easin_ids, "^R\\d{5}$"))) { #nolint
    ecokit::stop_ctx(
      "All easin_ids must be in the format 'RXXXXX', where X is an integer.",
      easin_ids = easin_ids,
      cat_timestamp = FALSE
    )
  }
  # Match EASIN IDs to species names and fact sheets
  ecokit::cat_time(
    "Match EASIN IDs to species names and fact sheets",
    cat_timestamp = FALSE,
    cat_bold = TRUE,
    cat_red = TRUE
  )

  matched_taxa <- purrr::map_dfr(
    .x = easin_ids,
    .f = ~ {

      easin_taxononly_url <- "https://easin.jrc.ec.europa.eu/apixg/catxg"
      # Extract species data as tibble
      url <- stringr::str_glue("{easin_taxononly_url}/easinid/{.x}") #nolint
      taxa_data <- try(RCurl::getURL(url, .mapUnicode = FALSE), silent = TRUE)
      if (inherits(taxa_data, "try-error")) {
        break
      }
      if (stringr::str_detect(taxa_data, "There are no results")) {
        return(
          dplyr::tibble(
            EASINID = .x,
            matched_species = NA_character_,
            species_fact_sheet = NA_character_
          ))
      }
      taxa_data <- jsonlite::fromJSON(taxa_data, flatten = TRUE) %>%
        dplyr::tibble() %>%
        dplyr::mutate(matched_species = paste0(Name, " ", Authorship)) %>%
        dplyr::select(
          tidyselect::all_of(c("Name", "EASINID", "matched_species")))

      species_fact_sheet <- paste0(
        "https://easin.jrc.ec.europa.eu/spexplorer/species/factsheet/", .x)
      # Check that URL is valid
      if (httr::http_error(species_fact_sheet) ||
          !ecokit::check_url(species_fact_sheet)) {
        taxa_data <- dplyr::mutate(
          taxa_data, species_fact_sheet = NA_character_)
      } else {
        taxa_data <- dplyr::mutate(
          taxa_data, species_fact_sheet = species_fact_sheet)
      }
      taxa_data
    }) %>%
    dplyr::mutate(
      matched = paste0(
        EASINID,
        dplyr::if_else(
          is.na(Name),
          ": not matched with EASIN database",
          paste0(": ", matched_species)),

        dplyr::if_else(
          is.na(species_fact_sheet),
          "",
          paste0(
            " (",
            cli::style_hyperlink(
              text = crayon::blue("fact sheet"),
              url = species_fact_sheet
            ),
            ")"
          ))
      )
    )

  if (all(is.na(matched_taxa$matched_species))) {
    ecokit::cat_time(
      paste0(
        "All provided EASIN IDs were not matched with the EASIN database.\n",
        "No EASIN data will be downloaded."),
      cat_timestamp = FALSE, level = 1L
    )
    return(invisible(NULL))
  }

  if (anyNA(matched_taxa$matched_species)) {
    skipped_ids <- matched_taxa %>%
      dplyr::filter(is.na(matched_species)) %>%
      dplyr::pull(EASINID) %>%
      toString()
    ecokit::cat_time(
      paste0(
        "Some provided EASIN ID(s) were not matched with the EASIN database.\n",
        "  >>>  These EASIN ID(s) will be skipped: ",
        crayon::red(skipped_ids)),
      cat_timestamp = FALSE, level = 1L
    )
  }

  matched_taxa %>%
    dplyr::filter(!is.na(matched_species)) %>%
    dplyr::pull(matched) %>%
    paste(collapse = "\n  >>>  ") %>%
    ecokit::cat_time(cat_timestamp = FALSE, level = 1L)


  if (is.null(model_dir)) {
    ecokit::stop_ctx(
      paste0(
        "The model_dir argument must be provided either directly or via the ",
        "`onesdm_model_dir` option."
      ),
      cat_timestamp = FALSE
    )
  }

  path_data <- fs::path(model_dir, "data")
  path_easin_data <- fs::path(path_data, "easin_data.RData")

  if (ecokit::check_data(path_easin_data, warning = FALSE)) {
    if (!overwrite) {
      ecokit::cat_time(
        paste0(
          "\nEASIN data file already exists at: ",
          crayon::blue(path_easin_data)),
        cat_timestamp = FALSE, verbose = verbose)
      ecokit::cat_time(
        paste0(
          "  >>>  Use overwrite = TRUE to re-download ",
          "and re-process the data."),
        cat_timestamp = FALSE, verbose = verbose)

      if (return_data) {
        ecokit::cat_time(
          "Loading EASIN data file", cat_timestamp = FALSE, verbose = verbose)
        return(invisible(ecokit::load_as(path_easin_data)))
      } else {
        return(invisible(path_easin_data))
      }
    }

    ecokit::cat_time(
      crayon::blue(
        "\nEASIN data file already exists and will be overwritten"),
      cat_timestamp = FALSE, cat_bold = TRUE, verbose = verbose)

  }

  fs::dir_create(path_data)

  # # ********************************************************************** #
  # Print function arguments ------
  # # ********************************************************************** #

  if (verbose) {
    ecokit::cat_time(
      "\nEASIN data extraction parameters:",
      cat_timestamp = FALSE,
      cat_bold = TRUE,
      cat_red = TRUE
    )
    ecokit::cat_time(
      paste0(
        crayon::italic("EASIN ID(s): "), crayon::blue(toString(easin_ids))),
      level = 1L,
      cat_timestamp = FALSE
    )
    ecokit::cat_time(
      paste0(
        crayon::italic("Modelling directory: "), crayon::blue(model_dir)),
      level = 1L, cat_timestamp = FALSE)
    ecokit::cat_time(
      paste0(
        crayon::italic("Modelling directory (absolute): "),
        crayon::blue(fs::path_abs(model_dir))),
      level = 1L, cat_timestamp = FALSE)
    ecokit::cat_time(
      paste0(crayon::italic("Timeout: "), crayon::blue(timeout)),
      level = 1L,
      cat_timestamp = FALSE
    )
    ecokit::cat_time(
      paste0(
        crayon::italic("Number of search results per chunk: "),
        crayon::blue(format(n_search, big.mark = ",", scientific = FALSE))),
      level = 1L,
      cat_timestamp = FALSE
    )
    ecokit::cat_time(
      paste0(
        crayon::italic("Number of download attempts: "),
        crayon::blue(n_attempts)),
      level = 1L,
      cat_timestamp = FALSE
    )
    ecokit::cat_time(
      paste0(crayon::italic("Sleep time: "), crayon::blue(sleep_time)),
      level = 1L,
      cat_timestamp = FALSE
    )
    ecokit::cat_time(
      paste0(crayon::italic("Exclude GBIF: "), crayon::blue(exclude_gbif)),
      level = 1L,
      cat_timestamp = FALSE
    )
    ecokit::cat_time(
      paste0(crayon::italic("Start year: "), crayon::blue(start_year)),
      level = 1L,
      cat_timestamp = FALSE
    )
    ecokit::cat_time(
      paste0(crayon::italic("Overwrite: "), crayon::blue(overwrite)),
      level = 1L,
      cat_timestamp = FALSE
    )
    ecokit::cat_time(
      paste0(
        crayon::italic("Return processed data: "), crayon::blue(return_data)),
      level = 1L,
      cat_timestamp = FALSE
    )
    ecokit::cat_time(
      "\nExtracting EASIN data",
      cat_timestamp = FALSE, cat_bold = TRUE, cat_red = TRUE)
  }

  # Temporarily set download time out only within the function
  withr::local_options(list(scipen = 999L, timeout = timeout))

  # # ********************************************************************** #
  # Downloading EASIN data -----
  # # ********************************************************************** #

  file_data_raw <- fs::path(path_data, "easin_data_raw.RData")

  if (!ecokit::check_data(file_data_raw, warning = FALSE) || overwrite) {

    # Download EASIN data for each provided EASIN ID
    easin_data <- purrr::map_dfr(
      .x = easin_ids,
      .f = get_easin_internal,
      path_data = path_data,
      timeout = timeout,
      n_search = n_search,
      n_attempts = n_attempts,
      sleep_time = sleep_time,
      exclude_gbif = exclude_gbif,
      verbose = verbose
    )

    # Check if any data were downloaded
    if (nrow(easin_data) == 0L) {
      ecokit::cat_time(
        "\nNo EASIN data were downloaded for the provided EASIN IDs.",
        cat_timestamp = FALSE,
        verbose = verbose)

      if (return_data) {
        return(tibble::tibble())
      } else {
        return(invisible(NA_character_))
      }
    }

    # Save raw EASIN data before filtering
    ecokit::cat_time(
      paste0("Saving raw EASIN data to: `", crayon::blue(file_data_raw), "`"),
      cat_timestamp = FALSE, verbose = verbose, level = 1L)
    ecokit::save_as(
      object = easin_data, object_name = "easin_data_raw",
      out_path = file_data_raw)

  } else {
    easin_data <- ecokit::load_as(file_data_raw)
  }

  # # ********************************************************************** #
  # Cleaning EASIN data -----
  # # ********************************************************************** #

  ecokit::cat_time(
    "\nCleaning EASIN data", cat_timestamp = FALSE,
    cat_bold = TRUE, cat_red = TRUE, verbose = verbose)

  n_rows_raw <- nrow(easin_data)

  # # ||||||||||||||||||||||||||||||||||||||||| #
  # Initial filtering of EASIN data
  # # ||||||||||||||||||||||||||||||||||||||||| #

  ecokit::cat_time(
    paste0("Discarding empty coordinates and data older than ", start_year),
    cat_timestamp = FALSE, verbose = verbose)

  easin_data <- easin_data %>%
    dplyr::mutate(Year = as.integer(Year)) %>%
    dplyr::filter(!is.na(WKT), Year >= start_year)

  n_rows_1 <- nrow(easin_data)
  if (n_rows_1 == 0L) {
    ecokit::cat_time(
      "No EASIN data were extracted after filtering.",
      cat_timestamp = FALSE, level = 1L, verbose = verbose)

    if (return_data) {
      return(tibble::tibble())
    } else {
      return(invisible(NA_character_))
    }
  }

  ecokit::cat_time(
    paste0(
      "Filtered out ",
      format(n_rows_raw - n_rows_1, big.mark = ",", scientific = FALSE),
      " records with empty coordinates or older than ", start_year, "."),
    cat_timestamp = FALSE, level = 1L, verbose = verbose)

  # # ||||||||||||||||||||||||||||||||||||||||| #
  # Extracting coordinates from WKT strings and filtering observations without
  # coordinates
  # # ||||||||||||||||||||||||||||||||||||||||| #

  ecokit::cat_time(
    paste0(
      "Extracting coordinates from WKT strings and filtering ",
      "observations without coordinates"),
    cat_timestamp = FALSE, verbose = verbose)

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
              dplyr::bind_rows()
          } else {
            tibble::tibble(longitude = NA_real_, latitude = NA_real_)
          }
        }
      )
    ) %>%
    tidyr::unnest(point_coords) %>%
    dplyr::filter(!is.na(longitude) & !is.na(latitude))

  n_rows_2 <- nrow(easin_data)
  if (n_rows_2 == 0L) {
    ecokit::cat_time(
      "No EASIN data were extracted after filtering.",
      cat_timestamp = FALSE, level = 1L, verbose = verbose)

    if (return_data) {
      return(tibble::tibble())
    } else {
      return(invisible(NA_character_))
    }
  }

  ecokit::cat_time(
    paste0(
      "filltered out ",
      crayon::blue(
        format(n_rows_1 - n_rows_2, big.mark = ",", scientific = FALSE)),
      " records with invalid coordinates."),
    cat_timestamp = FALSE, level = 1L, verbose = verbose)

  # # ||||||||||||||||||||||||||||||||||||||||| #
  # Exclude observations with low coordinate precision
  # # ||||||||||||||||||||||||||||||||||||||||| #

  ecokit::cat_time(
    "Exclude observations with low coordinate precision",
    cat_timestamp = FALSE, verbose = verbose)

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

  n_rows_3 <- nrow(easin_data)
  if (n_rows_3 == 0L) {
    ecokit::cat_time(
      "No EASIN data were extracted after filtering.",
      cat_timestamp = FALSE, level = 1L, verbose = verbose)

    if (return_data) {
      return(tibble::tibble())
    } else {
      return(invisible(NA_character_))
    }
  }

  ecokit::cat_time(
    paste0(
      "filltered out ",
      crayon::blue(
        format(n_rows_2 - n_rows_3, big.mark = ",", scientific = FALSE)),
      " records with low spatial precision or equal longitude and latitude."),
    cat_timestamp = FALSE, level = 1L, verbose = verbose)

  # # ||||||||||||||||||||||||||||||||||||||||| #
  # Further cleaning of EASIN data using CoordinateCleaner package
  # # ||||||||||||||||||||||||||||||||||||||||| #

  ecokit::cat_time(
    "cleaning using `CoordinateCleaner` package",
    cat_timestamp = FALSE, verbose = verbose)

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

  n_rows_4 <- nrow(easin_data)
  if (n_rows_4 == 0L) {
    ecokit::cat_time(
      "No EASIN data were extracted after filtering.",
      cat_timestamp = FALSE, level = 1L, verbose = verbose)

    if (return_data) {
      return(tibble::tibble())
    } else {
      return(invisible(NA_character_))
    }
  }

  ecokit::cat_time(
    paste0(
      "filltered out ",
      crayon::blue(
        format(n_rows_3 - n_rows_4, big.mark = ",", scientific = FALSE)),
      " records using `CoordinateCleaner`."),
    cat_timestamp = FALSE, level = 1L, verbose = verbose)

  # # ||||||||||||||||||||||||||||||||||||||||| #
  # Converting to sf object
  # # ||||||||||||||||||||||||||||||||||||||||| #

  easin_data <- easin_data %>%
    # convert to sf object, while keeping original coordinates as columns
    sf::st_as_sf(
      coords = c("longitude", "latitude"), crs = 4326L, remove = FALSE)

  # # ||||||||||||||||||||||||||||||||||||||||| #
  # Printing summary of extracted EASIN data
  # # ||||||||||||||||||||||||||||||||||||||||| #

  n_obs_per_id <- sf::st_drop_geometry(easin_data) %>%
    dplyr::group_by(SpeciesId) %>%
    dplyr::tally(name = "n_obs") %>%
    dplyr::mutate(n_obs = paste0(SpeciesId, " (", n_obs, ")")) %>%
    dplyr::pull(n_obs) %>%
    format(big.mark = ",", scientific = FALSE) %>%
    toString()

  ecokit::cat_time(
    paste0(
      "A total of ",
      crayon::blue(format(n_rows_4, big.mark = ",", scientific = FALSE)),
      " filtered observations were extracted for EASIN ID(s): ",
      crayon::blue(n_obs_per_id)),
    cat_timestamp = FALSE, verbose = verbose)

  # # ||||||||||||||||||||||||||||||||||||||||| #
  # Saving EASIN data
  # # ||||||||||||||||||||||||||||||||||||||||| #

  ecokit::cat_time(
    paste0(
      "Saving EASIN data to: `", crayon::blue(path_easin_data), "`"),
    cat_timestamp = FALSE, verbose = verbose)
  save(easin_data, file = path_easin_data)

  # # ********************************************************************** #
  # # ********************************************************************** #

  ecokit::cat_diff(
    init_time = .start_easin_time,
    prefix = "\nExtracting EASIN data was finished in ", verbose = verbose)

  # # ********************************************************************** #
  # # ********************************************************************** #

  if (return_data) {
    return(easin_data)
  } else {
    return(invisible(path_easin_data))
  }

}



# Internal function to download EASIN chunk data

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
        "{easin_url}/{easin_id}/exclude/dps/1/{skip}/{n_search}"  #nolint
      )
    } else {
      url <- stringr::str_glue("{easin_url}/{easin_id}/{skip}/{n_search}") #nolint
    }

    chunk_name <- paste0("easin_", easin_id, "_chunk_", chunk_n)
    chunk_file <- fs::path(path_data, paste0(chunk_name, ".RData"))

    while (download_try < n_attempts) {
      download_try <- download_try + 1L

      ecokit::cat_time(
        paste0(
          cli::style_hyperlink(
            text = crayon::blue(paste0("chunk ", chunk_n)), url = url),
          " (attempt ", download_try, ")"),
        level = 2L, cat_timestamp = FALSE, verbose = verbose)

      if (ecokit::check_data(chunk_file, warning = FALSE)) {
        ecokit::cat_time(
          "loading chunk data from disk",
          level = 3L, cat_timestamp = FALSE, verbose = verbose)
        easin_data_sub[[chunk_n]] <- chunk_data <- ecokit::load_as(chunk_file)
        break
      }

      chunk_data <- try(
        RCurl::getURL(url, .mapUnicode = FALSE, timeout = timeout),
        silent = TRUE)

      if (inherits(chunk_data, "try-error")) {
        chunk_data <- tibble::tibble()
        break
      }

      no_obs <- stringr::str_detect(
        chunk_data,
        pattern = "There are no results based on your"
      )

      if (no_obs) {
        chunk_data <- tibble::tibble()
        break
      }

      chunk_data <- jsonlite::fromJSON(chunk_data, flatten = TRUE) %>%
        tibble::tibble() %>%
        dplyr::mutate(json_url = url)

      if (inherits(chunk_data, "data.frame")) {
        easin_data_sub[[chunk_n]] <- chunk_data
        ecokit::save_as(
          object = chunk_data, object_name = chunk_name, out_path = chunk_file)
        break
      }
    }

    if (inherits(chunk_data, "data.frame")) {
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
      crayon::blue(
        format(nrow(easin_data_sub), big.mark = ",", scientific = FALSE)),
      " observations were extracted for EASIN ID: ", crayon::blue(easin_id)),
    cat_timestamp = FALSE, level = 2L, verbose = verbose)

  return(easin_data_sub)
}
