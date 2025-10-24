# load_mask_layer -----

#' Load or Create a Mask Layer at Specified Resolution
#'
#' Loads a mask layer (raster) at a specified model resolution from the OneSDM
#' Open Science Framework (OSF) storage. It can either return a `SpatRaster`
#' object directly or saved to a specified path.
#'
#' @param resolution Numeric. Spatial resolution. Valid values are 5, 10, or 20
#'   for resolutions of approximately 10km, 20km, and 40km (2.5, 5, and 10
#'   arc-minutes) respectively. Can be set via the `onesdm_resolution` option.
#'   Default is `NULL`.
#' @param climate_dir Character. Directory where climate data and mask layers
#'   are stored. If `NULL` (default), a temporary file is created.
#' @param verbose Logical. Should a download progress bar be shown? Default is
#'   `FALSE`.
#' @param overwrite Logical. Should existing files be overwritten? Default is
#'   `FALSE`.
#' @param return_spatraster Logical. If `TRUE` (default), returns a `SpatRaster`
#'   object. If `FALSE`, saves the raster to a file and returns the file path
#'   (invisibly).
#' @param wrap Logical. Should the resulting `SpatRaster` be wrapped using
#'   [terra::wrap]? Default is `FALSE`.
#'
#' @return A `SpatRaster` object representing the mask layer at the specified
#'   resolution, or the file path (invisibly) if `return_spatraster` is `FALSE`.
#'
#' @examples
#' \dontrun{
#' # Load mask layer as `SpatRaster` object
#' mask <- load_mask_layer(resolution = 10L)
#'
#' # Save mask layer to file
#' load_mask_layer(resolution = 20L, climate_dir = NULL)
#' }
#'
#' @export
#' @author Ahmed El-Gabbas

load_mask_layer <- function(
    resolution = NULL, climate_dir = NULL, verbose = FALSE,
    overwrite = FALSE, return_spatraster = TRUE, wrap = FALSE) {

  name <- NULL

  # # ********************************************************************** #
  # Validate inputs ------
  # # ********************************************************************** #

  resolution <- ecokit::assign_from_options(
    resolution, "onesdm_resolution", c("numeric", "integer"))
  climate_dir <- ecokit::assign_from_options(
    climate_dir, "onesdm_climate_dir", "character", allow_null = TRUE)
  verbose <- ecokit::assign_from_options(
    verbose, "onesdm_verbose", "logical")

  ecokit::check_args(
    c("overwrite", "wrap", "return_spatraster"), "logical")

  ## resolution -----

  if (length(resolution) != 1L || !is.numeric(resolution)) {
    ecokit::stop_ctx(
      paste0(
        "`resolution` must be a single numeric value. ",
        "Provided: ", toString(resolution), "."),
      cat_timestamp = FALSE)
  }

  valid_resolutions <- c(5L, 10L, 20L)
  if (!resolution %in% valid_resolutions) {
    ecokit::stop_ctx(
      paste0(
        "Invalid `resolution`: ", resolution, ".\n",
        "Valid options are: ", toString(valid_resolutions), "."),
      cat_timestamp = FALSE)
  }

  ## climate_dir -----

  if (is.null(climate_dir)) {
    save_path <- fs::file_temp(
      pattern = paste0("mask_agg", resolution, "_"), ext = ".tif")
  } else {
    if (!fs::dir_exists(climate_dir)) {
      ecokit::cat_time(
        paste0(
          "Creating climate data directory at: ", crayon::blue(climate_dir)),
        cat_timestamp = FALSE)
      fs::dir_create(climate_dir)
    }
    save_path <- fs::path(climate_dir, paste0("mask_agg_", resolution, ".tif"))
  }

  # # ********************************************************************** #
  # Get mask links and download mask file ------
  # # ********************************************************************** #

  if (!ecokit::check_tiff(save_path, warning = FALSE) || overwrite) {
    mask_file <- osfr::osf_retrieve_node("zarx2") %>%
      osfr::osf_ls_files(n_max = 20L) %>%
      dplyr::filter(name == paste0("mask_agg_", resolution, ".tif"))

    if (nrow(mask_file) == 0L) {
      ecokit::stop_ctx(
        paste0(
          "Mask file for `resolution` ", resolution, " not found on OSF."),
        cat_timestamp = FALSE)
    }
    if (nrow(mask_file) > 1L) {
      ecokit::stop_ctx(
        paste0(
          "Multiple mask files for `resolution` ", resolution,
          " found on OSF."),
        mask_file = mask_file, cat_timestamp = FALSE)
    }

    file_down <- httr::GET(
      mask_file$meta[[1L]]$links$download,
      httr::write_disk(save_path, overwrite = overwrite),
      httr::timeout(300L), if (verbose) httr::progress())
    rm(file_down, envir = environment())
  }

  # # ********************************************************************** #

  if (return_spatraster) {
    mask_r <- terra::toMemory(terra::rast(save_path))
    # Set variable name to empty string for downstream compatibility
    terra::varnames(mask_r) <- ""
    if (wrap) {
      mask_r <- terra::wrap(mask_r)
    }
    mask_r
  } else {
    invisible(save_path)
  }
}


# get_climate_data ------

#' Download climate data files from the OneSDM OSF project
#'
#' Downloads selected climate raster files (GeoTIFF) from the OneSDM Open
#' Science Framework (OSF) storage based on spatial resolution, climate
#' scenario, climate model, time period, and variable names. The function
#' validates inputs, locates the corresponding files in the OSF project, creates
#' local directories as needed, downloads the files, and verifies their
#' integrity.
#'
#' @param climate_dir Character. Destination directory where climate data files
#'   will be saved. Can be set via the `onesdm_climate_dir` option. The same
#'   directory should be used in case of modelling multiple species to ensure
#'   consistency. Default is `NULL`.
#' @param resolution Numeric. Spatial resolution. Valid values are 5, 10, or 20
#'   for resolutions of approximately 10km, 20km, and 40km (2.5, 5, and 10
#'   arc-minutes) respectively. Can be set via the `onesdm_resolution` option.
#'   Default is `10L`.
#' @param climate_scenario Character scalar. Climate scenario; one of `current`,
#'   `ssp126`, `ssp370`, `ssp585`.
#' @param climate_model Character scalar. Abbreviation of Global Circulation
#'   Models; one of `current`, `gfdl`, `ipsl`, `mpi`, `mri`, `ukesm1`.
#' @param year Character scalar. Time period; one of `1981-2010`, `2011-2040`,
#'   `2041-2070`, `2071-2100`.
#' @param var_names Character vector of climate variable codes to download. See
#'   [OneSDM::climate_data] for the list of valid variable names.
#' @param verbose Logical scalar. If `TRUE`, prints progress and informative
#'   messages.
#'
#' @details
#' - The function filters the [OneSDM::climate_data] data to identify files
#' matching the requested parameter combination, ensures a single OSF directory,
#' lists files via the `osfr` package, and joins metadata to obtain download
#' links.
#' - GeoTIFF files are downloaded, if not already available and valid, using
#' [httr::GET()] (with a maximum timeout of 5 minutes per file) into a mirrored
#' subdirectory structure under the `climate_dir`. Output files follow the
#' pattern: `<climate_dir>/res_<resolution>/1981_2010/var_name.tif` for current
#' climate, and `<climate_dir>/res_<resolution>/<year>_<climate_scenario>_`
#' `<climate_model>/var_name.tif` for future climate.
#' - Each file is validated as a GeoTIFF using [ecokit::check_tiff()]. Failed
#' downloads are reported.
#' - Input validation is strict; any invalid values, missing combinations, or
#' inconsistent OSF directory layouts will trigger informative errors.
#'
#' @return Invisibly returns `NULL`. On success, the requested `.tif` files are
#'   written to disk.
#'
#' @seealso
#' - [OneSDM::climate_data] for the metadata table used to resolve files.
#' - [ecokit::chelsa_var_info] and [CHELSA](https://chelsa-climate.org/)
#' webpage for details on climatic variables description.
#'
#' @examples
#' \dontrun{
#' ecokit::load_packages(fs, terra)
#'
#' # Create a temporary directory for climate data
#' tmp_dir <- fs::path_temp("onesdm_climate123")
#' fs::dir_create(tmp_dir)
#'
#' # |||||||||||||||||||||||||
#'
#' # Example usage of get_climate_data
#' get_climate_data(
#'   climate_dir = tmp_dir,
#'   resolution = 20L,
#'   climate_scenario = "current",
#'   climate_model = "current",
#'   year = "1981-2010",
#'   var_names = c("bio1", "bio12"),
#'   verbose = TRUE)
#'
#' # List downloaded files in the subdirectory of combination of parameters
#' print(list.dirs(tmp_dir))
#'
#' list_files <- list.files(
#'   fs::path(tmp_dir, "res_20", "1981_2010"), full.names = TRUE)
#' print(list_files)
#'
#' # Load the downloaded files as a SpatRaster stack
#' terra::rast(list_files)
#'
#' # |||||||||||||||||||||||||
#'
#' # Example using future climate data set using options to set defaults
#'
#' options(
#'   onesdm_climate_dir = tmp_dir, onesdm_resolution = 20L,
#'   onesdm_verbose = TRUE)
#' get_climate_data(
#'   climate_scenario = "ssp585", climate_model = "ukesm1",
#'   year = "2071-2100", var_names = c("bio5", "bio6", "npp"))
#'
#' # List downloaded files in the subdirectory of combination of parameters
#' print(list.dirs(tmp_dir))
#'
#' # Load the downloaded files as a SpatRaster stack
#' list_files_future <- list.files(
#'   fs::path(tmp_dir, "res_20", "2071_2100_ssp585_ukesm1"),
#'  full.names = TRUE)
#'
#' print(list_files_future)
#'
#' # Load the downloaded files as a SpatRaster stack
#' terra::rast(list_files_future)
#'
#' # |||||||||||||||||||||||||
#'
#' # Clean up temporary directory after use
#' fs::dir_delete(tmp_dir)
#' }
#'
#' @author Ahmed El-Gabbas
#' @export

get_climate_data <- function(
    climate_dir = NULL, resolution = 10L,
    climate_scenario = "current", climate_model = "current",
    year = "1981-2010", var_names = NULL, verbose = TRUE) {

  osf_path <- download_link <- out_dir <- meta <- out_file <-
    name <- climate_model_abb <- var_name <- NULL

  ecokit::check_args("verbose", "logical")

  ecokit::info_chunk(
    "Download climate data files from OneSDM OSF project",
    line_char_rep = 65L, verbose = verbose)

  # # ********************************************************************** #
  # Validate inputs ------
  # # ********************************************************************** #

  climate_dir <- ecokit::assign_from_options(
    climate_dir, "onesdm_climate_dir", "character")
  resolution <- ecokit::assign_from_options(
    resolution, "onesdm_resolution", c("numeric", "integer"))
  var_names <- ecokit::assign_from_options(
    var_names, "onesdm_var_names", "character")

  ecokit::cat_time(
    "Validate inputs for climate data download",
    cat_timestamp = FALSE, verbose = verbose)

  # # |||||||||||||||||||||||||||||||||||||||||||||||||||||| #

  ## climate_dir -----

  if (is.null(climate_dir)) {
    ecokit::stop_ctx(
      paste0(
        "`climate_dir` must be provided as a character string ",
        "indicating the directory to save the climate data files."),
      cat_timestamp = FALSE)
  }
  # ecokit::check_args("climate_dir", "character")
  if (!fs::dir_exists(climate_dir)) {
    ecokit::cat_time(
      paste0(
        "Creating climate data directory at: ",
        crayon::blue(climate_dir)),
      cat_timestamp = FALSE, verbose = verbose)
    fs::dir_create(climate_dir)
  }

  # # |||||||||||||||||||||||||||||||||||||||||||||||||||||| #

  ## resolution ----
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

  # # |||||||||||||||||||||||||||||||||||||||||||||||||||||| #

  ## climate models ----
  valid_models <- c("current", "gfdl", "ipsl", "mpi", "mri", "ukesm1")
  if (length(climate_model) != 1L || !is.character(climate_model)) {
    ecokit::stop_ctx(
      paste0(
        "`climate_model` must be a single character string. ",
        "Provided: ", toString(climate_model), "."),
      cat_timestamp = FALSE)
  }
  if (!climate_model %in% valid_models) {
    ecokit::stop_ctx(
      paste0(
        "Invalid `climate_model`: ", climate_model, ".\n",
        "Valid options are: ", toString(valid_models), "."),
      cat_timestamp = FALSE)
  }

  # # |||||||||||||||||||||||||||||||||||||||||||||||||||||| #

  ## climate scenarios ----
  valid_scenarios <- c("current", "ssp126", "ssp370", "ssp585")
  if (length(climate_scenario) != 1L || !is.character(climate_scenario)) {
    ecokit::stop_ctx(
      paste0(
        "`climate_scenario` must be a single character string. ",
        "Provided: ", toString(climate_scenario), "."),
      cat_timestamp = FALSE)
  }
  if (!climate_scenario %in% valid_scenarios) {
    ecokit::stop_ctx(
      paste0(
        "Invalid `climate_scenario`: ", climate_scenario, ".\n",
        "Valid options are: ", toString(valid_scenarios), "."),
      cat_timestamp = FALSE)
  }

  # # |||||||||||||||||||||||||||||||||||||||||||||||||||||| #

  ## years ----
  valid_years <- c("1981-2010", "2011-2040", "2041-2070", "2071-2100")
  if (length(year) != 1L || !is.character(year)) {
    ecokit::stop_ctx(
      paste0(
        "`year` must be a single character string. ",
        "Provided: ", toString(year), "."),
      cat_timestamp = FALSE)
  }
  if (!year %in% valid_years) {
    ecokit::stop_ctx(
      paste0(
        "Invalid `year`: ", year, ".\nValid options are: ",
        toString(valid_years), "."),
      cat_timestamp = FALSE)
  }

  # # |||||||||||||||||||||||||||||||||||||||||||||||||||||| #

  ## var_names ----
  valid_var_names <- c(
    "bio1", "bio2", "bio3", "bio4", "bio5", "bio6", "bio7", "bio8", "bio9",
    "bio10", "bio11", "bio12", "bio13", "bio14", "bio15", "bio16", "bio17",
    "bio18", "bio19", "fcf", "fgd", "gdd0", "gdd10", "gdd5", "gddlgd0",
    "gddlgd10", "gddlgd5", "gdgfgd0", "gdgfgd10", "gdgfgd5", "gsl", "gsp",
    "gst", "lgd", "ngd0", "ngd10", "ngd5", "npp", "scd", "swe")
  if (is.null(var_names)) {
    ecokit::stop_ctx(
      paste0(
        "`var_names` must be provided as a character vector.\n",
        "Valid options are: ", toString(valid_var_names), "."),
      cat_timestamp = FALSE)
  }
  if (!all(var_names %in% valid_var_names)) {
    ecokit::stop_ctx(
      paste0(
        "Invalid `var_names` values: ",
        toString(var_names[!var_names %in% valid_var_names]), ".\n",
        "Valid options are: ", toString(valid_var_names), "."),
      cat_timestamp = FALSE)
  }

  # # ********************************************************************** #
  # Print function arguments ------
  # # ********************************************************************** #

  if (verbose) {
    ecokit::cat_time(
      "Function arguments for climate data download", cat_timestamp = FALSE)
    ecokit::cat_time(
      paste0(
        crayon::italic("Climate directory: "), crayon::blue(climate_dir)),
      level = 1L, cat_timestamp = FALSE)
    ecokit::cat_time(
      paste0(
        crayon::italic("Climate directory (absolute): "),
        crayon::blue(fs::path_abs(climate_dir))),
      level = 1L, cat_timestamp = FALSE)
    ecokit::cat_time(
      paste0(crayon::italic("resolution: "), crayon::blue(resolution)),
      level = 1L, cat_timestamp = FALSE)
    ecokit::cat_time(
      paste0(
        crayon::italic("climate scenario: "), crayon::blue(climate_scenario)),
      level = 1L, cat_timestamp = FALSE)
    ecokit::cat_time(
      paste0(crayon::italic("climate model: "), crayon::blue(climate_model)),
      level = 1L, cat_timestamp = FALSE)
    ecokit::cat_time(
      paste0(crayon::italic("climate year: "), crayon::blue(year)),
      level = 1L, cat_timestamp = FALSE)
    ecokit::cat_time(
      crayon::italic("Climate variables: "),
      level = 1L, cat_timestamp = FALSE)
    stringr::str_wrap(toString(var_names), 50L) %>%
      stringr::str_split("\n", simplify = TRUE) %>%
      purrr::walk(
        ~ ecokit::cat_time(crayon::blue(.x), level = 2L, cat_timestamp = FALSE))
  }

  # # ********************************************************************** #
  # Filter climate_data to get files to download ------
  # # ********************************************************************** #

  ecokit::cat_time(
    "Preparing climate data download links",
    cat_timestamp = FALSE, verbose = verbose)

  download_links <- dplyr::filter(
    OneSDM::climate_data,
    resolution == !!resolution,
    climate_scenario == !!climate_scenario,
    climate_model_abb == !!climate_model,
    year == !!year,
    var_name %in% !!var_names)

  if (nrow(download_links) == 0L) {
    ecokit::stop_ctx(
      "No climate data files found for the specified parameter combination.",
      cat_timestamp = FALSE,
      resolution = resolution, climate_scenario = climate_scenario,
      climate_model_abb = climate_model, year = year,
      var_names = toString(var_names))
  }

  download_links <- download_links %>%
    dplyr::mutate(
      # directory name at the OSF project
      osf_dir = dirname(osf_path),
      # local directory to save the downloaded file
      out_dir = fs::path(climate_dir, out_dir),
      # full local path to save the downloaded file
      out_file = fs::path(climate_dir, out_file),
      # file name to match with OSF files
      name = paste0(var_name, ".tif"),
      out_okay = purrr::map_lgl(out_file, ecokit::check_tiff, warning = FALSE))

  if (!all(var_names %in% download_links$var_name)) {
    ecokit::stop_ctx(
      paste0(
        "Some requested climate variable files were not found for the ",
        "specified parameter combination: ",
        toString(var_names[!var_names %in% download_links$var_name]), "."),
      cat_timestamp = FALSE)
  }

  if (all(download_links$out_okay)) {
    ecokit::cat_time(
      "All requested climate data files are already downloaded and valid.",
      cat_timestamp = FALSE, verbose = verbose)
    return(invisible(NULL))
  }

  # Create download directories
  fs::dir_create(unique(download_links$out_dir))

  # # ********************************************************************** #
  # Extract download links ------
  # # ********************************************************************** #

  # Ensure a single unique OSF directory
  osf_dir_unique <- unique(download_links$osf_dir)
  if (length(osf_dir_unique) != 1L) {
    ecokit::stop_ctx(
      paste0(
        "Expected a single unique OSF directory but found: ",
        toString(osf_dir_unique), "."),
      cat_timestamp = FALSE)
  }

  # Get list of files in the specified OSF directory
  osf_l_files <- osfr::osf_retrieve_node("tuyxh") %>%
    osfr::osf_ls_files(n_max = 1000L, type = "folder") %>%
    dplyr::filter(name == osf_dir_unique) %>%
    osfr::osf_ls_files(n_max = 1000L, type = "file")

  ecokit::cat_time(
    "Downloading climate data files", cat_timestamp = FALSE, verbose = verbose)

  # Ensure unique names before joining
  if (any(duplicated(download_links$name))) {
    ecokit::stop_ctx(
      paste0(
        "Duplicate names found in requested download links: ",
        toString(download_links$name[duplicated(download_links$name)]), "."),
      cat_timestamp = FALSE)
  }

  if (any(duplicated(osf_l_files$name))) {
    ecokit::stop_ctx(
      paste0(
        "Duplicate names found in OSF files: ",
        toString(osf_l_files$name[duplicated(osf_l_files$name)]), "."),
      cat_timestamp = FALSE)
  }

  download_links <- download_links %>%
    # left join to get the metadata including download links
    dplyr::left_join(osf_l_files, by = "name") %>%
    dplyr::mutate(
      download_link = purrr::map_chr(meta, ~ .x$links$download),
      download_check = purrr::map2_lgl(
        .x = download_link, .y = out_file,
        .f = ~ {
          ecokit::cat_time(
            crayon::blue(stringr::str_remove(basename(.y), ".tif")),
            cat_timestamp = FALSE, level = 1L, verbose = verbose)

          if (!ecokit::check_tiff(.y, warning = FALSE)) {
            httr::GET(
              .x, httr::write_disk(path = .y, overwrite = TRUE),
              httr::timeout(300L))
          }
          if (!ecokit::check_tiff(.y, warning = FALSE)) {
            return(FALSE)
          }
          return(TRUE)
        },
        .progress = verbose))

  # # ********************************************************************** #
  # Check download results ------
  # # ********************************************************************** #

  if (!all(download_links$download_check)) {
    ecokit::stop_ctx(
      paste0(
        "Some climate data files failed to download: ",
        toString(download_links$name[!download_links$download_check]), "."),
      cat_timestamp = FALSE)
  }

  ecokit::cat_time(
    "All climate data files downloaded successfully.",
    cat_timestamp = FALSE, verbose = verbose)

  return(invisible(NULL))

}


# download_landuse_data ------

# download_landuse_data <- function(
    #
#   climate_dir = NULL, resolution = 10L,
#   climate_scenario = "current",
#   # pft
#   # year = "1981-2010", verbose = TRUE
#   # climate_model = "current",
#   # , var_names = NULL,
# ) {
#
# }
