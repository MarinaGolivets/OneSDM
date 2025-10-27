# prepare_species_data -------

#' Prepare Species Occurrence Data for SDM Modelling using OneSDM
#'
#' @description This function prepares species occurrence data from multiple
#'   sources (GBIF, EASIN, and user-provided coordinates) for species
#'   distribution modelling. It processes, validates, rasterizes the data, and
#'   optionally excludes spatial outliers and specified geographic extents.
#'   Internally, the function uses the [prepare_gbif_data()],
#'   [prepare_easin_data()], and [prepare_user_data()] functions to retrieve and
#'   process data from each source.
#'
#' @param easin_ids Character vector. One or more EASIN species IDs. Each ID
#'   should start with 'R' followed by five digits (e.g. "R00544"). This cannot
#'   be `NULL`. Species IDs can be found on the [EASIN
#'   website](https://easin.jrc.ec.europa.eu/spexplorer/search/) by searching
#'   for a species and locating the "EASIN ID" in the species details section.
#'   Can also be set via the `onesdm_easin_ids` option. For more details, see
#'   [prepare_easin_data()].
#' @param gbif_ids Character or numeric vector of GBIF taxon keys (as numeric
#'   strings) to query. If `NULL`, attempts to retrieve from the
#'   `onesdm_gbif_ids` option. For more details, see [prepare_gbif_data()].
#' @param coordinates A data frame or matrix containing user-provided longitude
#'   and latitude values. Must have exactly two columns. If not provided
#'   directly, the function will attempt to retrieve it from the
#'   `onesdm_coordinates` option. See [prepare_user_data()] for more details.
#' @param model_dir Character. Directory path where model outputs will be saved.
#'   Can be set via the `onesdm_model_dir` option. A subdirectory `data` will be
#'   created within this directory to store processed species data. A separate
#'   directory is expected in case of modelling multiple species, otherwise data
#'   files may be overwritten or mixed up. Default is `NULL`.
#' @param climate_dir Character. Directory path for climate data files. Can be
#'   set via the `onesdm_climate_dir` option. This directory is used to load
#'   mask layers matching the specified resolution. The same directory should be
#'   used in case of modelling multiple species to ensure consistency. Default
#'   is `NULL`.
#' @param resolution Numeric. Spatial resolution. valid values are 5, 10, or 20
#'   for resolutions of approximately 10km, 20km, and 40km (2.5, 5, and 10
#'   arc-minutes) respectively. Can be set via the `onesdm_resolution` option.
#'   Default is `NULL`.
#' @param verbose Logical. If `TRUE` (default), prints progress messages. Can be
#'   set via the `onesdm_verbose` option.
#' @param exclude_extents List. A list of `SpatExtent` objects (created using
#'   [terra::ext] function), defining geographic areas to exclude from the
#'   species data. Default is an empty list, meaning no areas are excluded. Can
#'   be set via the `onesdm_exclude_extents` option.
#' @param species_name Character. Name of the species for labelling outputs.
#'   This will not be used to retrieve data from `GBIF` or `EASIN`. Can be set
#'   via the `onesdm_species_name` option. Default is `"species"`.
#' @param outlier_dist_km Numeric. Distance threshold in kilometers for
#'   identifying spatial outliers. If `0L` (default), no outlier detection is
#'   performed. Can be set via the `onesdm_outlier_dist_km` option. For more
#'   details, see [ecokit::nearest_dist_sf()].
#' @param outlier_resolution Numeric. Spatial resolution for outlier detection
#'   calculations. Can be set via the `onesdm_outlier_resolution` option. Only
#'   used if `outlier_dist_km` is greater than 0. A coarser resolution can speed
#'   up calculations. See [ecokit::nearest_dist_sf()] for more details. Default
#'   is 0.125.
#' @param outlier_n_cores Integer. Number of CPU cores to use for outlier
#'   detection. Can be set via the `onesdm_outlier_n_cores` option. Only used if
#'   `outlier_dist_km` is greater than 0. Default is `6L`.
#' @param plot_distribution Logical. If `TRUE` (default), generates a JPEG map
#'   of species distribution. Can be set via the `onesdm_plot_distribution`
#'   option.
#'
#' @details
#' - The function creates a `data` subdirectory in `model_dir` and saves
#' multiple outputs including raw data, processed data, rasterized data, and
#' optional plots.
#' - The function performs the following steps:
#'   - Validates input parameters
#'   - Retrieves species data from GBIF, EASIN, and/or user coordinates using
#' [prepare_gbif_data()], [prepare_easin_data()], and [prepare_user_data()].
#' **Note** that
#'     - at least one of `easin_ids`, `gbif_ids`, or `coordinates` must be
#' provided.
#'     - not all arguments of these functions can be set via this wrapper
#' function. Please use respective options to set those arguments and refer to
#' the respective function documentation for details.
#'
#'   - Merges data from all sources
#'   - Optionally removes spatial outliers based on nearest neighbour distances
#'   - Rasterizes occurrence data to match the specified resolution
#'   - Excludes specified geographic extents if provided
#'   - Saves processed data and generates visualizations if requested. The map
#' shows the rasterized species distribution overlaid on a world map. If
#' `exclude_extents` is provided, these areas are highlighted on the map in red.
#' If a specific extent is used to retrieve GBIF data (via the
#' `onesdm_gbif_boundaries` option in [prepare_gbif_data()]), it is also
#' highlighted on the map in green
#' - Function default arguments can be set globally using the `options()`
#' function. Users can set these options at the start of their R session to
#' avoid repeatedly specifying them in function calls. The following options
#' correspond to the function arguments:
#'   - `onesdm_easin_ids`: for the `easin_ids` argument.
#'   - `onesdm_gbif_ids`: for the `gbif_ids` argument.
#'   - `onesdm_coordinates`: for the `coordinates` argument.
#'   - `onesdm_model_dir`: for the `model_dir` argument.
#'   - `onesdm_climate_dir`: for the `climate_dir` argument.
#'   - `onesdm_resolution`: for the `resolution` argument.
#'   - `onesdm_verbose`: for the `verbose` argument.
#'   - `onesdm_exclude_extents`: for the `exclude_extents` argument.
#'   - `onesdm_species_name`: for the `species_name` argument.
#'   - `onesdm_outlier_dist_km`: for the `outlier_dist_km` argument.
#'   - `onesdm_outlier_resolution`: for the `outlier_resolution` argument.
#'   - `onesdm_outlier_n_cores`: for the `outlier_n_cores` argument.
#'   - `onesdm_plot_distribution`: for the `plot_distribution` argument.
#'   - This in addition to options used in the internal functions called.
#'
#' @return Invisibly returns `NULL`. The function is called for its side effects
#'   (saving processed data files and generating plots).
#'
#' @examples
#' \dontrun{
#' options(
#'   onesdm_gbif_verbose = TRUE,
#'   onesdm_start_year = 1990L,
#'   onesdm_gbif_overwrite = FALSE,
#'   onesdm_gbif_boundaries = c(-170L, 170L, -80L, 80L),
#'   onesdm_gbif_max_uncertainty = 5L)
#'
#' prepare_species_data(
#'   easin_ids = "R00042", gbif_ids = 2979128, coordinates = NULL,
#'   model_dir = "test/Acacia_karroo", climate_dir = "test/climate_data",
#'   resolution = 10L, verbose = TRUE,
#'   exclude_extents = list(terra::ext(110, 180, -50, -10)),
#'   species_name = "Acacia_karroo",
#'   outlier_dist_km = 200L, outlier_resolution = 0.25, outlier_n_cores = 4L)
#' }
#'
#' @export
#' @author Ahmed El-Gabbas

prepare_species_data <- function(
    easin_ids = NULL, gbif_ids = NULL, coordinates = NULL,
    model_dir = NULL, climate_dir = NULL, resolution = NULL, verbose = TRUE,
    exclude_extents = list(), species_name = "species",
    outlier_dist_km = 0L, outlier_resolution = 0.125, outlier_n_cores = 6L,
    plot_distribution = TRUE) {

  valid <- xmin <- xmax <- ymin <- ymax <- nearest_dist <- value <- NULL

  ecokit::info_chunk(
    "Preparing species occurrence data", line_char = "=", cat_bold =  TRUE,
    cat_red = TRUE, line_char_rep  = 65L, verbose = verbose, cat_date = FALSE)

  # # ********************************************************************** #
  # Check inputs ------
  # # ********************************************************************** #

  easin_ids <- ecokit::assign_from_options(
    easin_ids, "onesdm_easin_ids", "character", allow_null = TRUE)
  gbif_ids <- ecokit::assign_from_options(
    gbif_ids, "onesdm_gbif_ids",
    c("numeric", "character", "integer"), allow_null = TRUE)
  coordinates <- ecokit::assign_from_options(
    coordinates, "onesdm_coordinates",
    c("data.frame", "matrix"), allow_null = TRUE)
  model_dir <- ecokit::assign_from_options(
    model_dir, "onesdm_model_dir", "character")
  climate_dir <- ecokit::assign_from_options(
    climate_dir, "onesdm_climate_dir", "character")
  resolution <- ecokit::assign_from_options(
    resolution, "onesdm_resolution", c("numeric", "integer"))
  verbose <- ecokit::assign_from_options(verbose, "onesdm_verbose", "logical")
  species_name <- ecokit::assign_from_options(
    species_name, "onesdm_species_name", "character")
  outlier_dist_km <- ecokit::assign_from_options(
    outlier_dist_km, "onesdm_outlier_dist_km", c("numeric", "integer"))
  outlier_resolution <- ecokit::assign_from_options(
    outlier_resolution, "onesdm_outlier_resolution", c("numeric", "integer"))
  outlier_n_cores <- ecokit::assign_from_options(
    outlier_n_cores, "onesdm_outlier_n_cores", c("numeric", "integer"))
  plot_distribution <- ecokit::assign_from_options(
    plot_distribution, "onesdm_plot_distribution", "logical")

  # # ||||||||||||||||||||||||||||||||||||||||| #

  ## easin_ids - gbif_ids - coordinates ----

  # at least one of easin_ids, gbif_ids, coordinates must be provided
  if (is.null(easin_ids) && is.null(gbif_ids) && is.null(coordinates)) {
    ecokit::stop_ctx(
      paste0(
        "At least one of the arguments `easin_ids`, `gbif_ids`, ",
        "or `coordinates` must be provided either directly or via the ",
        "corresponding options."),
      easin_ids = easin_ids, gbif_ids = gbif_ids, coordinates = coordinates,
      cat_timestamp = FALSE)
  }

  # # ||||||||||||||||||||||||||||||||||||||||| #

  ## model_dir -----

  if (is.null(model_dir)) {
    ecokit::stop_ctx(
      paste0(
        "The model_dir argument must be provided either directly or via the ",
        "onesdm_model_dir option."),
      cat_timestamp = FALSE)
  }

  # # ||||||||||||||||||||||||||||||||||||||||| #

  ## climate_dir -----

  if (is.null(climate_dir)) {
    ecokit::stop_ctx(
      paste0(
        "`climate_dir` must be provided as a character string ",
        "indicating the directory to save the climate data files."),
      cat_timestamp = FALSE)
  }

  if (!fs::dir_exists(climate_dir)) {
    ecokit::cat_time(
      paste0(
        "Creating climate data directory at: ",
        crayon::blue(climate_dir)),
      cat_timestamp = FALSE, verbose = verbose)
    fs::dir_create(climate_dir)
  }

  # # ||||||||||||||||||||||||||||||||||||||||| #

  ## resolution -----

  valid_resolutions <- c(5L, 10L, 20L)
  if (length(resolution) != 1L || !is.numeric(resolution)) {
    ecokit::stop_ctx(
      paste0(
        "`resolution` must be a single numeric value. ",
        "Provided: ", toString(resolution), "."),
      cat_timestamp = FALSE)
  }
  if (!resolution %in% valid_resolutions) {
    ecokit::stop_ctx(
      paste0(
        "Invalid `resolution`: ", resolution, ".\n",
        "Valid options are: ", toString(valid_resolutions), "."),
      cat_timestamp = FALSE)
  }

  # # ||||||||||||||||||||||||||||||||||||||||| #

  ## exclude_extents ------

  # exclude_extents must be a list of SpatExtent objects
  if (!is.list(exclude_extents)) {
    ecokit::stop_ctx(
      "The exclude_extents argument must be a list.",
      exclude_extents = exclude_extents, cat_timestamp = FALSE)
  }

  if (length(exclude_extents) > 0L) {
    valid_extents <- purrr::map_lgl(exclude_extents, inherits, "SpatExtent")
    if (!all(valid_extents)) {
      ecokit::stop_ctx(
        paste0(
          "All elements in the `exclude_extents` argument must ",
          "be SpatExtent objects."),
        exclude_extents = exclude_extents,
        class_exclude_extents = purrr::map(exclude_extents, class),
        cat_timestamp = FALSE)
    }

    # check that all extents have valid coordinate ranges
    extents_vector <- purrr::map_dfr(exclude_extents, as.vector)
    invalid_extents_n <- extents_vector %>%
      dplyr::mutate(
        ext_id = seq_along(exclude_extents),
        valid = dplyr::case_when(
          xmin >= -180L & xmax <= 180L & ymin >= -90L & ymax <= 90L ~ TRUE,
          .default = FALSE)) %>%
      dplyr::filter(!valid)

    if (nrow(invalid_extents_n) > 0L) {
      ecokit::stop_ctx(
        paste0(
          "All extents in the `exclude_extents` argument must have valid ",
          "coordinate ranges: xmin >= -180, xmax <= 180, ",
          "ymin >= -90, ymax <= 90."),
        invalid_extents_n = invalid_extents_n,
        invalid_extents_ids = invalid_extents_n$ext_id,
        cat_timestamp = FALSE)
    }
  }

  # # ********************************************************************** #
  # File paths ------
  # # ********************************************************************** #

  path_data <- fs::path(model_dir, "data")
  file_data_raw <- fs::path(path_data, "species_data_raw.RData")
  file_data_outliers <- fs::path(path_data, "species_data_outliers.RData")
  file_data <- fs::path(path_data, "species_data.RData")
  file_data_raster <- fs::path(
    path_data, paste0("species_data_r_", resolution, "_km.tif"))
  file_data_final <- fs::path(
    path_data, paste0("species_data_r_", resolution, "_km_PA.tif"))
  file_plot <- fs::path(
    path_data, paste0("species_data_", resolution, "_km_plot.jpeg"))
  fs::dir_create(path_data)

  # # ********************************************************************** #
  # Print function arguments ------
  # # ********************************************************************** #

  if (verbose) {

    ecokit::cat_time(
      crayon::italic("Species data preparation parameters:"),
      cat_timestamp = FALSE, cat_bold = TRUE, cat_red = TRUE)

    ecokit::cat_time(
      paste0(crayon::italic("Species name: "), crayon::blue(species_name)),
      level = 1L, cat_timestamp = FALSE)
    ecokit::cat_time(
      paste0(
        crayon::italic("GBIF ID(s): "),
        crayon::blue(
          dplyr::if_else(length(gbif_ids) > 0L, toString(gbif_ids), "None"))),
      level = 1L, cat_timestamp = FALSE)
    ecokit::cat_time(
      paste0(
        crayon::italic("EASIN ID(s): "),
        crayon::blue(
          dplyr::if_else(length(easin_ids) > 0L, toString(easin_ids), "None"))),
      level = 1L, cat_timestamp = FALSE)
    ecokit::cat_time(
      paste0(
        crayon::italic("User coordinates: "),
        crayon::blue(
          dplyr::if_else(
            (!is.null(coordinates) && nrow(coordinates) > 0L),
            paste0(nrow(coordinates), " unique coordinates"), "None"))),
      level = 1L, cat_timestamp = FALSE)
    ecokit::cat_time(
      paste0(crayon::italic("Modelling directory: "), crayon::blue(model_dir)),
      level = 1L, cat_timestamp = FALSE)
    ecokit::cat_time(
      paste0(
        crayon::italic("Modelling directory (absolute): "),
        crayon::blue(fs::path_abs(model_dir))),
      level = 1L, cat_timestamp = FALSE)
    ecokit::cat_time(
      paste0(
        crayon::italic("Climate data directory: "), crayon::blue(climate_dir)),
      level = 1L, cat_timestamp = FALSE)
    ecokit::cat_time(
      paste0(
        crayon::italic("Climate data directory (absolute): "),
        crayon::blue(fs::path_abs(climate_dir))),
      level = 1L, cat_timestamp = FALSE)
    ecokit::cat_time(
      paste0(
        crayon::italic("Modelling resolution: "), crayon::blue(resolution)),
      level = 1L, cat_timestamp = FALSE)

    # exclude_extents <- list(ext(), ext())
    if (length(exclude_extents) > 0L) {
      ext_txt <- dplyr::mutate(
          extents_vector,
          ext_txt = paste0(
            "xmin: ", round(xmin, 2L), ", xmax: ", round(xmax, 2L),
            ", ymin: ", round(ymin, 2L), ", ymax: ", round(ymax, 2L))
        ) %>%
        dplyr::pull(ext_txt)

      ecokit::cat_time(
        crayon::italic("Excluded extents: "), level = 1L, cat_timestamp = FALSE)
      purrr::walk(ext_txt, ecokit::cat_time, level = 2L, cat_timestamp = FALSE)

    } else {
      ecokit::cat_time(
        paste0(
          crayon::italic("Excluded extents: "), crayon::blue("None")),
        level = 1L, cat_timestamp = FALSE)
    }

    if (outlier_dist_km > 0L) {
      ecokit::cat_time(
        paste0(
          crayon::italic("Outlier distance threshold: "),
          crayon::blue(outlier_dist_km), " km"),
        level = 1L, cat_timestamp = FALSE)
    } else {
      ecokit::cat_time(
        paste0(
          crayon::italic("Outlier distance threshold: "),
          crayon::blue("None")),
        level = 1L, cat_timestamp = FALSE)
    }
  }

  # # ********************************************************************** #
  # GBIF -----
  # # ********************************************************************** #

  if (is.null(gbif_ids)) {

    gbif_data <- tibble::tibble()

  } else {

    ecokit::info_chunk(
      paste0(
        "Preparing GBIF data for taxon ID(s): ",
        crayon::blue(toString(gbif_ids))),
      line_char_rep  = 65L, verbose = verbose, cat_date = FALSE,
      cat_bold = TRUE)

    start_year <- ecokit::get_option_with_default(
      "onesdm_start_year", "OneSDM::prepare_gbif_data", "start_year")
    r_environ <- ecokit::get_option_with_default(
      "onesdm_r_environ", "OneSDM::prepare_gbif_data", "r_environ")
    gbif_boundaries <- ecokit::get_option_with_default(
      "onesdm_gbif_boundaries", "OneSDM::prepare_gbif_data", "boundaries")
    max_uncertainty <- ecokit::get_option_with_default(
      "onesdm_gbif_max_uncertainty", "OneSDM::prepare_gbif_data",
      "max_uncertainty")
    overwrite_gbif <- ecokit::get_option_with_default(
      "onesdm_gbif_overwrite", "OneSDM::prepare_gbif_data", "overwrite")
    verbose_gbif <- ecokit::get_option_with_default(
      "onesdm_gbif_verbose", "OneSDM::prepare_gbif_data", "verbose")

    gbif_data <- OneSDM::prepare_gbif_data(
      gbif_ids = gbif_ids, model_dir = model_dir, verbose = verbose_gbif,
      start_year = start_year, r_environ = r_environ,
      boundaries = gbif_boundaries, max_uncertainty = max_uncertainty,
      overwrite = overwrite_gbif, return_data = TRUE)

    if (nrow(gbif_data) > 0L) {
      gbif_data <- dplyr::select(
        gbif_data, tidyselect::all_of(c("longitude", "latitude"))) %>%
        dplyr::mutate(source = "gbif")
    }
  }

  # # ********************************************************************** #
  # EASIN -----
  # # ********************************************************************** #

  if (is.null(easin_ids)) {

    easin_data <-  tibble::tibble()

  } else {

    ecokit::info_chunk(
      paste0(
        "Preparing EASIN data for taxon ID(s): ",
        crayon::blue(toString(easin_ids))),
      line_char_rep  = 65L, verbose = verbose, cat_date = FALSE,
      cat_bold = TRUE)

    timeout <- ecokit::get_option_with_default(
      "onesdm_easin_timeout", "OneSDM::prepare_easin_data", "timeout")
    n_search <- ecokit::get_option_with_default(
      "onesdm_easin_n_search", "OneSDM::prepare_easin_data", "n_search")
    n_attempts <- ecokit::get_option_with_default(
      "onesdm_easin_n_attempts", "OneSDM::prepare_easin_data", "n_attempts")
    sleep_time <- ecokit::get_option_with_default(
      "onesdm_easin_sleep_time", "OneSDM::prepare_easin_data", "sleep_time")
    exclude_gbif <- ecokit::get_option_with_default(
      "onesdm_easin_exclude_gbif", "OneSDM::prepare_easin_data", "exclude_gbif")
    start_year <- ecokit::get_option_with_default(
      "onesdm_start_year", "OneSDM::prepare_easin_data", "start_year")
    overwrite_easin <- ecokit::get_option_with_default(
      "onesdm_easin_overwrite", "OneSDM::prepare_easin_data", "overwrite")
    verbose_easin <- ecokit::get_option_with_default(
      "onesdm_easin_verbose", "OneSDM::prepare_easin_data", "verbose")

    easin_data <- OneSDM::prepare_easin_data(
      easin_ids = easin_ids, model_dir = model_dir, timeout = timeout,
      n_search = n_search, n_attempts = n_attempts, sleep_time = sleep_time,
      exclude_gbif = exclude_gbif, verbose = verbose_easin,
      start_year = start_year, overwrite = overwrite_easin, return_data = TRUE)

    if (nrow(easin_data) > 0L) {
      easin_data <- dplyr::select(
        easin_data, tidyselect::all_of(c("longitude", "latitude"))) %>%
        dplyr::mutate(source = "easin")
    }
  }

  # # ********************************************************************** #
  # User coordinates -------
  # # ********************************************************************** #

  if (is.null(coordinates)) {
    coordinates_data <-  tibble::tibble()
  } else {
    ecokit::info_chunk(
      "Processing user-defined coordinates",
      line_char_rep  = 65L, verbose = verbose, cat_date = FALSE,
      cat_bold = TRUE)

    coordinates_data <- OneSDM::prepare_user_data(
      coordinates = coordinates, model_dir = model_dir, return_data = TRUE) %>%
      dplyr::mutate(source = "user")
  }

  # # ********************************************************************** #
  # Merge all data sources ------
  # # ********************************************************************** #

  ecokit::info_chunk(
    "Processing combined species data",
    line_char_rep  = 65L, verbose = verbose, cat_date = FALSE,
    cat_bold = TRUE)

  # Initialize an empty sf data frame
  species_data <- sf::st_sf(
    longitude = numeric(0L), latitude = numeric(0L),
    geometry = sf::st_sfc(crs = 4326L))

  if (nrow(gbif_data) > 0L) {
    species_data <- dplyr::bind_rows(species_data, gbif_data)
  }
  if (nrow(easin_data) > 0L) {
    species_data <- dplyr::bind_rows(species_data, easin_data)
  }
  if (nrow(coordinates_data) > 0L) {
    species_data <- dplyr::bind_rows(species_data, coordinates_data)
  }

  if (nrow(species_data) == 0L) {
    ecokit::stop_ctx(
      "No species occurrence data available from the provided sources.",
      easin_ids = easin_ids, gbif_ids = gbif_ids,
      coordinates = coordinates, cat_timestamp = FALSE)
  }

  # Save the total number of records before excluding outliers

  ecokit::cat_time(
    paste0(
      "Saving raw merged species data to: ", crayon::blue(file_data_raw)),
    verbose = verbose, cat_timestamp = FALSE)
  ecokit::save_as(
    object = species_data, object_name = "species_data_raw",
    out_path = file_data_raw)

  n_records <- nrow(species_data)
  ecokit::cat_time(
    paste0(
      "Total number of species records: ",
      crayon::blue(format(n_records, big.mark = ",", scientific = FALSE))),
    verbose = verbose, cat_timestamp = FALSE)

  # Ensure species_data2 is always defined
  species_data2 <- species_data

  # # ********************************************************************** #
  # Exclude spatial outliers ------
  # # ********************************************************************** #

  if (outlier_dist_km > 0L) {

    if (!requireNamespace("nngeo", quietly = TRUE)) {
      ecokit::stop_ctx(
        paste0(
          "The nngeo package is required to exclude spatial ",
          "outliers. Please install it using: `install.packages('nngeo')`"),
        cat_timestamp = FALSE)
    }

    # CoordinateCleaner::cc_outl does not work with large datasets
    # https://github.com/ropensci/CoordinateCleaner/issues/110
    # species_data2 <- CoordinateCleaner::clean_coordinates(
    #   x = sf::st_drop_geometry(species_data),
    #   lon = "longitude", lat = "latitude", species = "species",
    #   tests = "outliers", outliers_method = "distance",
    #   outliers_td = outlier_dist_km)

    ecokit::cat_time(
      paste0(
        "Excluding spatial outliers beyond ",
        crayon::blue(outlier_dist_km), " km distance threshold using ",
        crayon::blue(outlier_n_cores), " cores"),
      verbose = verbose, cat_timestamp = FALSE)

    species_data2 <- ecokit::nearest_dist_sf(
      sf_object = species_data, resolution = outlier_resolution,
      n_cores = outlier_n_cores)

    excluded_data <- dplyr::filter(
      species_data2, nearest_dist >= outlier_dist_km)
    species_data2 <- dplyr::filter(
      species_data2, nearest_dist < outlier_dist_km)

    n_excluded <- n_records - nrow(species_data2)
    ecokit::cat_time(
      paste0(
        "Number of records excluded as spatial outliers: ",
        crayon::blue(format(n_excluded, big.mark = ",", scientific = FALSE)),
        " (", round((n_excluded / n_records) * 100L, 2L), "%)"),
      verbose = verbose, cat_timestamp = FALSE, level = 1L)

    if (nrow(excluded_data) > 0L) {
      ecokit::cat_time(
        paste0(
          "Saving excluded outlier species data to: ",
          crayon::blue(file_data_outliers)),
        verbose = verbose, cat_timestamp = FALSE)
      ecokit::save_as(
        object = excluded_data, object_name = "species_data_outliers",
        out_path = file_data_outliers)
    }

  }

  # # ********************************************************************** #
  # Saving species data ------
  # # ********************************************************************** #

  ecokit::cat_time(
    paste0(
      "Saving processed species data to: ", crayon::blue(file_data)),
    verbose = verbose, cat_timestamp = FALSE)

  ecokit::save_as(
    object = species_data2, object_name = "species_data",
    out_path = file_data)


  # # ********************************************************************** #
  # Load mask layer ------
  # # ********************************************************************** #

  ecokit::cat_time(
    paste0(
      "Loading mask layer at resolution of ", crayon::blue(resolution), " km"),
    verbose = verbose, cat_timestamp = FALSE)

  mask_layer <- OneSDM::load_mask_layer(
    resolution = resolution, climate_dir = climate_dir, verbose = verbose,
    overwrite = FALSE, return_spatraster = TRUE, wrap = FALSE)

  # # ********************************************************************** #
  # Rasterize species observations ------
  # # ********************************************************************** #

  ecokit::cat_time(
    "Rasterizing species occurrence data to match mask layer resolution",
    verbose = verbose, cat_timestamp = FALSE)

  species_data_r <- terra::rasterize(
    x = species_data2, y = mask_layer, field = 1L, background = 0L)
  species_data_r <- species_data_r * mask_layer

  # Number of grid cells with species presence
  n_cells_presence <- terra::global(
    species_data_r, fun = "sum", na.rm = TRUE)[1L, 1L]
  ecokit::cat_time(
    paste0(
      "Number of grid cells at resolution of ", crayon::blue(resolution),
      " km with species presence: ",
      crayon::blue(format(
        n_cells_presence, big.mark = ",", scientific = FALSE))),
    verbose = verbose, cat_timestamp = FALSE, level = 1L)

  # Save the number of presence cells before exclusion

  ecokit::cat_time(
    paste0(
      "Saving rasterized species data to: ",
      crayon::blue(file_data_raster)),
    verbose = verbose, cat_timestamp = FALSE)
  terra::writeRaster(
    species_data_r, filename = file_data_raster, overwrite = TRUE,
    gdal = c("COMPRESS=ZSTD", "ZSTD_LEVEL=22", "TILED=YES"))

  # # ********************************************************************** #
  # Exclude spatial extents ------
  # # ********************************************************************** #

  if (length(exclude_extents) > 0L) {

    ecokit::cat_time(
      "Excluding provided spatial extent(s) from species data",
      verbose = verbose, cat_timestamp = FALSE)

    for (ext in exclude_extents) {
      species_data_r <- terra::mask(
        species_data_r, ext, inverse = TRUE, updatevalue = 0L)
    }
    species_data_r <- species_data_r * mask_layer

    n_cells_presence_2 <- terra::global(
      species_data_r, fun = "sum", na.rm = TRUE)[1L, 1L]
    n_cells_excluded <- n_cells_presence - n_cells_presence_2
    ecokit::cat_time(
      paste0(
        "Number of grid cells excluded: ",
        crayon::blue(format(
          n_cells_excluded, big.mark = ",", scientific = FALSE)),
        " (", round((n_cells_excluded / n_cells_presence) * 100L, 2L), "%)"),
      verbose = verbose, cat_timestamp = FALSE, level = 1L)
    ecokit::cat_time(
      paste0(
        "Number of presence grid cells after exclusion: ",
        crayon::blue(
          format(n_cells_presence_2, big.mark = ",", scientific = FALSE))),
      verbose = verbose, cat_timestamp = FALSE, level = 1L)

    # Save the final rasterized species data after exclusion
    ecokit::cat_time(
      paste0(
        "Saving final rasterized species data to: ",
        crayon::blue(file_data_final)),
      verbose = verbose, cat_timestamp = FALSE)

    terra::writeRaster(
      species_data_r, filename = file_data_final, overwrite = TRUE,
      gdal = c("COMPRESS=ZSTD", "ZSTD_LEVEL=22", "TILED=YES"))

  } else {
    ecokit::cat_time(
      "No spatial extents provided for exclusion",
      verbose = verbose, cat_timestamp = FALSE)

    # copy file_data_raster to file_data_final
    fs::file_copy(
      path = file_data_raster, new_path = file_data_final, overwrite = TRUE)

  }

  # # ********************************************************************** #
  # Plotting ------
  # # ********************************************************************** #

  n_cells_presence <- terra::global(
    species_data_r, fun = "sum", na.rm = TRUE)[1L, 1L]

  if (plot_distribution && n_cells_presence > 0L) {

    ecokit::cat_time(
      paste0(
        "Generating species distribution plot at: ", crayon::blue(file_plot)),
      verbose = verbose, cat_timestamp = FALSE)

    ecokit::check_packages(c("rworldmap", "tidyterra", "ragg"))

    # aggregate species distribution to 0.5 degree resolution for plotting
    species_data_r_4_plot <- species_data_r %>%
      terra::aggregate(
        fact = round(0.5 / (resolution / 111.32)),
        fun = "max", na.rm = TRUE) %>%
      terra::classify(cbind(0L, NA)) %>%
      terra::as.factor()

    ecokit::quietly({
      species_data_plot <- ggplot2::ggplot() +
        ggplot2::geom_sf(
          data = sf::st_as_sf(rworldmap::getMap(resolution = "low")),
          fill = "lightgray", color = "black", linewidth = 0.1) +
        tidyterra::geom_spatraster(
          data = species_data_r_4_plot, maxcell = 25e4L,
          mapping = ggplot2::aes(fill = ggplot2::after_stat(value)),
          na.rm = TRUE, alpha = 0.8, show.legend = FALSE) +
        ggplot2::scale_fill_manual(values = "blue", na.value = "transparent") +
        ggplot2::coord_sf(
          xlim = c(-180L, 180L), ylim = c(-60L, 85L), expand = FALSE) +
        ggplot2::annotate(
          "text", x = -175L, y = -57L,
          label = paste0("Distribution of ", species_name, " data"),
          hjust = 0L, vjust = 0L, size = 2L, fontface = "bold")
    },
    "resampled to")
    rm(species_data_r_4_plot, species_data_r, envir = environment())

    if (length(exclude_extents) > 0L) {
      # Plot excluded extents if provided
      species_data_plot <- species_data_plot +
        ggplot2::geom_rect(
          data = extents_vector,
          mapping = ggplot2::aes(
            xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
          fill = "red", color = "red", alpha = 0.1, linewidth = 0.1,
          linetype = "dashed")
    }

    if (!is.null(gbif_ids) &&
        !identical(gbif_boundaries, c(-180L, 180L, -90L, 90L))) {
      # Plot GBIF gbif_boundaries if not global
      boundaries_data <- as.data.frame(matrix(gbif_boundaries, nrow = 1L)) %>%
        stats::setNames(c("xmin", "xmax", "ymin", "ymax"))
      species_data_plot <- species_data_plot +
        ggplot2::geom_rect(
          data = tibble::tibble(boundaries_data),
          mapping = ggplot2::aes(
            xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
          fill = "green", color = "green", alpha = 0.1, linewidth = 0.1,
          linetype = "dashed")
    }

    species_data_plot <- species_data_plot +
      ggplot2::theme(
        text = ggplot2::element_text(size = 6L),
        axis.title = ggplot2::element_blank(),
        axis.ticks = ggplot2::element_blank(),
        axis.text.x = ggplot2::element_text(
          size = 4L, margin = ggplot2::margin(t = -1L, b = 0L)),
        axis.text.y = ggplot2::element_text(
          size = 4L, margin = ggplot2::margin(r = -1L, l = 0L)),
        panel.grid = ggplot2::element_line(
          linewidth = 0.1, color = "gray80"),
        panel.background = ggplot2::element_rect(fill = "white"),
        plot.margin = grid::unit(c(0L, 0L, 0L, 0L), "cm"),
        plot.background = ggplot2::element_blank(),
        legend.position = "none")

    ragg::agg_jpeg(
      filename = file_plot, width = 18L, height = 7.5,
      res = 600L, quality = 100L, units = "cm")
    print(species_data_plot)
    grDevices::dev.off()
  }

  return(invisible(NULL))
}
