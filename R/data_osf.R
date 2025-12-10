# get_mask_layer -----

#' Load or Create a Mask Layer at Specified Resolution
#'
#' Loads a mask layer (raster) at a specified model resolution from the OneSDM
#' Open Science Framework (OSF) storage. It can either return a `SpatRaster`
#' object directly or saved to a specified path.
#'
#' @param resolution Numeric. Spatial resolution. Valid values are 5, 10, or 20
#'   for resolutions of approximately 5, 10, and 20 km (2.5, 5, and 10
#'   arc-minutes) respectively. Can be set via the "`onesdm_resolution`" option.
#'   Default is `NULL`.
#' @param climate_dir Character. Directory where climate data and mask layers
#'   are stored. If `NULL` (default), a temporary file is created.
#' @param verbose Logical. Should a download progress bar be shown? Default is
#'   `FALSE`.
#' @param europe_only Logical. If `TRUE`, downloads the Europe-only mask layer.
#'   Default is `FALSE`.
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
#'   require(terra)
#'   # Load mask layer as `SpatRaster` object
#'   mask <- get_mask_layer(resolution = 10L)
#'   print(mask)
#' }
#'
#' @export
#' @author Ahmed El-Gabbas

get_mask_layer <- function(
    resolution = NULL, climate_dir = NULL, verbose = FALSE, europe_only = FALSE,
    overwrite = FALSE, return_spatraster = TRUE, wrap = FALSE) {

  name <- NULL

  ecokit::check_packages(
    c("crayon", "dplyr", "fs", "httr", "osfr", "purrr", "stringr", "terra"))

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
    c("overwrite", "wrap", "return_spatraster", "europe_only"), "logical")

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
    if (europe_only) {
      save_path <- fs::file_temp(
        pattern = paste0("mask_europe_agg", resolution, "_"), ext = ".tif")
    } else {
      save_path <- fs::file_temp(
        pattern = paste0("mask_agg", resolution, "_"), ext = ".tif")
    }
  } else {
    fs::dir_create(climate_dir)
    if (europe_only) {
      save_path <- fs::path(
        climate_dir, paste0("mask_europe_agg_", resolution, ".tif"))
    } else {
      save_path <- fs::path(
        climate_dir, paste0("mask_agg_", resolution, ".tif"))
    }
  }

  # # ********************************************************************** #
  # Get mask links and download mask file ------
  # # ********************************************************************** #

  if (!ecokit::check_tiff(save_path, warning = FALSE) || overwrite) {

    if (europe_only) {
      mask_file <- osfr::osf_retrieve_node("zarx2") %>%
        osfr::osf_ls_files(n_max = 20L) %>%
        dplyr::filter(name == "europe") %>%
        osfr::osf_ls_files(n_max = 20L) %>%
        dplyr::filter(name == paste0("europe_mask_", resolution, ".tif"))
    } else {
      mask_file <- osfr::osf_retrieve_node("zarx2") %>%
        osfr::osf_ls_files(n_max = 20L) %>%
        dplyr::filter(name == paste0("mask_agg_", resolution, ".tif"))
    }

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

# # ********************************************************************** #
# # ********************************************************************** #


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
#' @param climate_dir Character. Destination directory where climate and land
#'   use data files will be saved. Can be set via the "`onesdm_climate_dir`"
#'   option. The same directory should be used in case of modelling multiple
#'   species to ensure consistency. Default is `NULL`.
#' @param resolution Numeric. Spatial resolution. Valid values are 5, 10, or 20
#'   for resolutions of approximately 5, 10, and 20 km (2.5, 5, and 10
#'   arc-minutes) respectively. Can be set via the "`onesdm_resolution`" option.
#'   Default is `10L`.
#' @param climate_scenario Character scalar. Climate scenario; one of `current`,
#'   `ssp126`, `ssp370`, `ssp585`. Default is `"current"`.
#' @param climate_model Character scalar. Abbreviation of Global Circulation
#'   Models; one of `current`, `gfdl`, `ipsl`, `mpi`, `mri`, `ukesm1`. Default
#'   is `"current"`.
#' @param year Character scalar. Time period; one of `1981_2010`, `2011_2040`,
#'   `2041_2070`, `2071_2100`. Default is `"1981_2010"`.
#' @param var_names Character vector of climate variable codes to download. See
#'   [OneSDM::climate_data] for the list of valid climate variable names. Can be
#'   set via the "`onesdm_var_names`" option. This parameter is *required* and
#'   cannot be `NULL`.
#' @param verbose Logical scalar. If `TRUE` (default), prints progress and
#'   informative messages.
#' @param sleep_time Numeric scalar. Number of seconds to pause after
#'   downloading files to avoid overwhelming the server. Default is `1L`.
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
#' @return Returns invisibly a tibble with metadata of the downloaded climate
#'   data files, including their local paths.
#'
#' @seealso
#' - [OneSDM::climate_data] for the metadata table used to resolve files.
#' - [ecokit::chelsa_var_info] and [CHELSA](https://chelsa-climate.org/)
#' webpage for details on climatic variables description.
#'
#' @examples
#' \dontrun{
#'   ecokit::load_packages(fs, terra)
#'
#'   # Create a temporary directory for climate data
#'   tmp_dir <- fs::path_temp("onesdm_climate")
#'   fs::dir_create(tmp_dir)
#'
#'   # |||||||||||||||||||||||||
#'
#'   # Example usage of get_climate_data
#'   get_climate_data(
#'     climate_dir = tmp_dir,
#'     resolution = 20L,
#'     climate_scenario = "current",
#'     climate_model = "current",
#'     year = "1981_2010",
#'     var_names = c("bio1", "bio12"),
#'     verbose = TRUE)
#'
#'   # List downloaded files in the subdirectory of combination of parameters
#'   print(list.dirs(tmp_dir))
#'
#'   list_files <- list.files(
#'     fs::path(tmp_dir, "res_20", "1981_2010"), full.names = TRUE)
#'   print(list_files)
#'
#'   # Load the downloaded files as a SpatRaster stack
#'   terra::rast(list_files)
#'
#'   # |||||||||||||||||||||||||
#'
#'   # Example using future climate data set using options to set defaults
#'
#'   options(
#'     onesdm_climate_dir = tmp_dir, onesdm_resolution = 20L,
#'     onesdm_verbose = TRUE)
#'   get_climate_data(
#'     climate_scenario = "ssp585", climate_model = "ukesm1",
#'     year = "2071_2100", var_names = c("bio5", "bio6", "npp"))
#'
#'   # List downloaded files in the subdirectory of combination of parameters
#'   print(list.dirs(tmp_dir))
#'
#'   # Load the downloaded files as a SpatRaster stack
#'   list_files_future <- list.files(
#'     fs::path(tmp_dir, "res_20", "2071_2100_ssp585", "ukesm1"),
#'     recursive = TRUE, pattern = "\\.tif$", full.names = TRUE)
#'
#'   print(list_files_future)
#'
#'   # Load the downloaded files as a SpatRaster stack
#'   terra::rast(list_files_future)
#'
#'   # |||||||||||||||||||||||||
#'
#'   # Clean up temporary directory after use
#'   fs::dir_delete(tmp_dir)
#' }
#'
#' @author Ahmed El-Gabbas
#' @export

get_climate_data <- function(
    climate_dir = NULL, resolution = 10L,
    climate_scenario = "current", climate_model = "current",
    year = "1981_2010", var_names = NULL, verbose = TRUE,
    sleep_time = 1L) {

  osf_path <- download_link <- out_dir <- meta <- out_file <-
    climate_model_abb <- var_name <- NULL

  ecokit::check_args("verbose", "logical")

  ecokit::info_chunk(
    "Download climate data files from OneSDM OSF project",
    line_char_rep = 65L, verbose = verbose, cat_date = FALSE)

  ecokit::check_packages(
    c("crayon", "dplyr", "fs", "httr", "osfr", "purrr", "stringr", "terra"))

  # # |||||||||||||||||||||||||||||||||||||||||||||||||||||| #

  # Load climate data tibble
  climate_data_0 <- OneSDM::climate_data

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
        "indicating the directory to save the climate and land use ",
        "data files."),
      cat_timestamp = FALSE)
  }

  fs::dir_create(climate_dir)

  # # |||||||||||||||||||||||||||||||||||||||||||||||||||||| #

  ## resolution ----

  valid_resolutions <- unique(climate_data_0$resolution)

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

  valid_models <- unique(climate_data_0$climate_model_abb)
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

  valid_scenarios <- unique(climate_data_0$climate_scenario)

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

  valid_years <- unique(climate_data_0$year)

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

  valid_var_names <- unique(climate_data_0$var_name)

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

  # # |||||||||||||||||||||||||||||||||||||||||||||||||||||| #

  # sleep_time -----

  if (length(sleep_time) != 1L || !is.numeric(sleep_time) || sleep_time < 0L) {
    ecokit::stop_ctx(
      paste0(
        "`sleep_time` must be a single non-negative numeric value. ",
        "Provided: ", toString(sleep_time), "."),
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
    climate_data_0,
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
    return(invisible(download_links))
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
    # Without using pattern osfr may return incomplete list of directories, with
    # some duplicates
    osfr::osf_ls_files(
      n_max = 1000L, type = "folder", pattern = osf_dir_unique) %>%
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
        .progress = FALSE))

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

  if (sleep_time > 0L) {
    Sys.sleep(sleep_time)
  }

  return(invisible(download_links))

}

# # ********************************************************************** #
# # ********************************************************************** #

# get_landuse_data ------

#' Download Land Use Data from OneSDM OSF Project
#'
#' @description Downloads selected land use data files from the OneSDM Open
#'   Science Framework (OSF) storage based on spatial resolution, climate
#'   scenario, climate model, time period, and Plant Functional Type (PFT)
#'   names. The function validates input parameters, checks for existing files,
#'   and downloads only missing or corrupted files.
#'
#' @param climate_dir Character. Destination directory where climate and land
#'   use data files will be saved. Can be set via the "`onesdm_climate_dir`"
#'   option. The same directory should be used in case of modelling multiple
#'   species to ensure consistency. Default is `NULL`.
#' @param resolution Numeric. Spatial resolution. Valid values are 5, 10, or 20
#'   for resolutions of approximately 5, 10, and 20 km (2.5, 5, and 10
#'   arc-minutes) respectively. Can be set via the "`onesdm_resolution`" option.
#'   Default is `10L`.
#' @param climate_scenario Character scalar. Climate scenario; one of `current`,
#'   `ssp126`, `ssp370`, `ssp585`. Default is `"current"`.
#' @param year Character scalar. Time period; one of `1981_2010`, `2011_2040`,
#'   `2041_2070`, `2071_2100`. Default is `"1981_2010"`.
#' @param pft_type Character string. Plant functional type category. Has to be
#'   one of "`cross-walk`" (default) "`original`"; see [OneSDM::landuse_data]
#'   for details. This can be set via the "`onesdm_pft_type`" option.
#' @param pft_id Numeric vector. One or more plant functional type identifiers
#'   to download. Must be valid for the specified `pft_type`: 1-20 for `pft_type
#'   == "original"`, and 1-12 for `pft_type == "cross-walk"`; See
#'   [OneSDM::landuse_data] for details. Can be set via the "`onesdm_pft_id`"
#'   option. This parameter has no default and must be provided.
#' @param verbose Logical scalar. If `TRUE` (default), prints progress and
#'   informative messages.
#' @param sleep_time Numeric scalar. Number of seconds to pause after
#'   downloading files to avoid overwhelming the server. Default is `1L`.
#'
#' @return Invisibly returns a tibble containing metadata about the downloaded
#'   files, including local file paths, OSF paths, download links, and
#'   validation status.
#'
#' @details The function performs the following steps:
#' - Validates all input parameters against available options
#' - Checks if requested files already exist locally and are valid
#' - Creates necessary output directories
#' - Connects to the OneSDM OSF project and downloads only missing or corrupted
#'   files
#' - Validates downloaded files using TIFF format checks
#'
#' @examples
#' \dontrun{
#'   ecokit::load_packages(fs, terra)
#'
#'   # Create a temporary directory for land use data
#'   tmp_dir <- fs::path_temp("onesdm_landuse")
#'   fs::dir_create(tmp_dir)
#'
#'   # |||||||||||||||||||||||||
#'
#'   # Example usage of get_landuse_data
#'   get_landuse_data(
#'     climate_dir = tmp_dir,
#'     resolution = 20L,
#'     climate_scenario = "current",
#'     year = "1981_2010",
#'     pft_type = "cross-walk",
#'     pft_id = c(2, 4, 12),
#'     verbose = TRUE)
#'
#'   # List downloaded files in the subdirectory of combination of parameters
#'   print(list.dirs(tmp_dir))
#'
#'   list_files <- list.files(
#'     fs::path(tmp_dir, "res_20", "1981_2010"), full.names = TRUE)
#'   print(list_files)
#'
#'   # Load the downloaded files as a SpatRaster stack
#'   terra::rast(list_files)
#'
#'   # |||||||||||||||||||||||||
#'
#'   # Example using future land use data set using options to set defaults
#'
#'   options(
#'     onesdm_climate_dir = tmp_dir, onesdm_resolution = 20L,
#'     onesdm_verbose = TRUE,
#'     onesdm_pft_type = "original",
#'     onesdm_pft_id = c(1, 14, 18))
#'
#'   get_landuse_data(climate_scenario = "ssp585", year = "2071_2100")
#'
#'   # List downloaded files in the subdirectory of combination of parameters
#'   print(list.dirs(tmp_dir))
#'
#'   # Load the downloaded files as a SpatRaster stack
#'   list_files_future <- list.files(
#'     fs::path(tmp_dir, "res_20", "2071_2100_ssp585"),
#'     recursive = TRUE, pattern = "\\.tif$", full.names = TRUE)
#'
#'   print(list_files_future)
#'
#'   # Load the downloaded files as a SpatRaster stack
#'   terra::rast(list_files_future)
#'
#'   # |||||||||||||||||||||||||
#'
#'   # Clean up temporary directory after use
#'   fs::dir_delete(tmp_dir)
#' }
#'
#' @author Ahmed El-Gabbas
#' @export

get_landuse_data <- function(
    climate_dir = NULL, resolution = 10L, climate_scenario = "current",
    year = "1981_2010", pft_type = "cross-walk", pft_id = NULL,
    verbose = TRUE, sleep_time = 1L) {

  osf_path <- download_link <- out_dir <- meta <- out_file <- name <- NULL

  ecokit::check_args("verbose", "logical")

  ecokit::info_chunk(
    "Download land use files from OneSDM OSF project",
    line_char_rep = 65L, verbose = verbose, cat_date = FALSE)

  ecokit::check_packages(
    c("crayon", "dplyr", "fs", "httr", "osfr", "purrr", "stringr", "terra"))

  # # |||||||||||||||||||||||||||||||||||||||||||||||||||||| #

  # Load land use data tibble
  landuse_data_0 <- OneSDM::landuse_data

  # # ********************************************************************** #
  # Validate inputs ------
  # # ********************************************************************** #

  climate_dir <- ecokit::assign_from_options(
    climate_dir, "onesdm_climate_dir", "character")
  resolution <- ecokit::assign_from_options(
    resolution, "onesdm_resolution", c("numeric", "integer"))
  pft_type <- ecokit::assign_from_options(
    pft_type, "onesdm_pft_type", "character")
  pft_id <- ecokit::assign_from_options(
    pft_id, "onesdm_pft_id", c("numeric", "integer"))

  ecokit::cat_time(
    "Validate inputs for land use data download",
    cat_timestamp = FALSE, verbose = verbose)

  # # |||||||||||||||||||||||||||||||||||||||||||||||||||||| #

  ## climate_dir -----

  if (is.null(climate_dir)) {
    ecokit::stop_ctx(
      paste0(
        "`climate_dir` must be provided as a character string ",
        "indicating the directory to save the climate and land ",
        "use data files."),
      cat_timestamp = FALSE)
  }

  fs::dir_create(climate_dir)

  # # |||||||||||||||||||||||||||||||||||||||||||||||||||||| #

  ## resolution ----

  valid_resolutions <- unique(landuse_data_0$resolution)

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

  ## climate scenarios ----

  valid_scenarios <- unique(landuse_data_0$climate_scenario)

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

  valid_years <- unique(landuse_data_0$year)

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

  ## pft_type ----

  valid_pft_types <- unique(landuse_data_0$pft_type)

  if (length(pft_type) != 1L || !is.character(pft_type)) {
    ecokit::stop_ctx(
      paste0(
        "`pft_type` must be a single character string. ",
        "Provided: ", toString(pft_type), "."),
      cat_timestamp = FALSE)
  }

  if (!pft_type %in% valid_pft_types) {
    ecokit::stop_ctx(
      paste0(
        "Invalid `pft_type`: ", pft_type, ".\nValid options are: ",
        toString(valid_pft_types), "."),
      cat_timestamp = FALSE)
  }

  # # |||||||||||||||||||||||||||||||||||||||||||||||||||||| #

  ## pft_id ----

  valid_pft_ids <- dplyr::filter(landuse_data_0, pft_type == !!pft_type) %>%
    dplyr::pull(pft_id) %>%
    unique()

  if (is.null(pft_id)) {
    ecokit::stop_ctx(
      paste0(
        "`pft_id` must be provided as a numeric vector.\n",
        "Valid options are: ", toString(valid_pft_ids), "."),
      cat_timestamp = FALSE)
  }

  if (!all(pft_id %in% valid_pft_ids)) {
    ecokit::stop_ctx(
      paste0(
        "Invalid `pft_id` values: ",
        toString(pft_id[!pft_id %in% valid_pft_ids]), ".\n",
        "Valid options are: ", toString(valid_pft_ids), "."),
      cat_timestamp = FALSE)
  }

  # # ********************************************************************** #

  # Get filtered land use data for printing
  download_links <- dplyr::filter(
    landuse_data_0, pft_type == !!pft_type, pft_id %in% !!pft_id,
    climate_scenario == !!climate_scenario, year == !!year,
    resolution == !!resolution)

  if (nrow(download_links) == 0L) {
    ecokit::stop_ctx(
      "No land use data files found for the specified parameter combination.",
      cat_timestamp = FALSE,
      resolution = resolution, climate_scenario = climate_scenario,
      year = year, pft = toString(download_links$pft_name),
      pft_id = toString(pft_id), pft_type = pft_type)
  }

  # # |||||||||||||||||||||||||||||||||||||||||||||||||||||| #

  # sleep_time -----

  if (length(sleep_time) != 1L || !is.numeric(sleep_time) || sleep_time < 0L) {
    ecokit::stop_ctx(
      paste0(
        "`sleep_time` must be a single non-negative numeric value. ",
        "Provided: ", toString(sleep_time), "."),
      cat_timestamp = FALSE)
  }

  # # ********************************************************************** #
  # Print function arguments ------
  # # ********************************************************************** #

  if (verbose) {
    ecokit::cat_time(
      "Function arguments for land use data download", cat_timestamp = FALSE)
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
      paste0(crayon::italic("climate year: "), crayon::blue(year)),
      level = 1L, cat_timestamp = FALSE)
    ecokit::cat_time(
      paste0(crayon::italic("Land use type: "), crayon::blue(pft_type)),
      level = 1L, cat_timestamp = FALSE)
    ecokit::cat_time(
      paste0(
        crayon::italic("land use id(s): "), crayon::blue(toString(pft_id))),
      level = 1L, cat_timestamp = FALSE)
    ecokit::cat_time(
      crayon::italic("Land use pft(s): "),
      level = 1L, cat_timestamp = FALSE)
    stringr::str_wrap(toString(download_links$pft_name), 50L) %>%
      stringr::str_split("\n", simplify = TRUE) %>%
      purrr::walk(
        ~ ecokit::cat_time(crayon::blue(.x), level = 2L, cat_timestamp = FALSE))
  }

  # # ********************************************************************** #
  # Filter land use data to get files to download ------
  # # ********************************************************************** #

  ecokit::cat_time(
    "Preparing land use data download links",
    cat_timestamp = FALSE, verbose = verbose)

  download_links <- download_links %>%
    dplyr::mutate(
      # directory name at the OSF project
      osf_dir = dirname(osf_path),
      # local directory to save the downloaded file
      out_dir = fs::path(climate_dir, out_dir),
      # full local path to save the downloaded file
      out_file = fs::path(climate_dir, out_file),
      # file name to match with OSF files
      name = basename(osf_path),
      out_okay = purrr::map_lgl(out_file, ecokit::check_tiff, warning = FALSE))

  if (all(download_links$out_okay)) {
    ecokit::cat_time(
      "All requested land use data files are already downloaded and valid.",
      cat_timestamp = FALSE, verbose = verbose)
    return(invisible(download_links))
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
  osf_l_files <- osfr::osf_retrieve_node("jvns4") %>%
    osfr::osf_ls_files(
      n_max = 1000L, type = "folder", pattern = osf_dir_unique) %>%
    osfr::osf_ls_files(n_max = 1000L, type = "file") %>%
    dplyr::filter(name %in% download_links$name)

  ecokit::cat_time(
    "Downloading land use data files",
    cat_timestamp = FALSE, verbose = verbose)

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
    dplyr::left_join(
      dplyr::select(osf_l_files, -tidyselect::all_of("id")), by = "name") %>%
    dplyr::mutate(
      download_link = purrr::map_chr(meta, ~ .x$links$download),
      download_check = purrr::map2_lgl(
        .x = download_link, .y = out_file,
        .f = ~ {

          map_name <- stringr::str_remove(basename(.y), ".tif")

          ecokit::cat_time(
            crayon::blue(map_name), cat_timestamp = FALSE,
            level = 1L, verbose = verbose)

          if (!ecokit::check_tiff(.y, warning = FALSE)) {
            # Write to temp file to change variable name before saving
            temp_tif <- fs::file_temp(pattern = "landuse_", ext = ".tif")
            httr::GET(
              .x, httr::write_disk(path = temp_tif, overwrite = TRUE),
              httr::timeout(300L))

            # Rename the raster layer to match pft_name
            lu_rast <- stats::setNames(terra::rast(temp_tif), map_name)
            # save the renamed raster to the final location
            terra::writeRaster(lu_rast, filename = .y, overwrite = TRUE)

            # Clean up temp file
            try(fs::file_delete(temp_tif), silent = TRUE)
          }

          if (!ecokit::check_tiff(.y, warning = FALSE)) {
            return(FALSE)
          }
          return(TRUE)
        },
        .progress = FALSE))

  # # ********************************************************************** #
  # Check download results ------
  # # ********************************************************************** #

  if (!all(download_links$download_check)) {
    ecokit::stop_ctx(
      paste0(
        "Some land use data files failed to download: ",
        toString(download_links$name[!download_links$download_check]), "."),
      cat_timestamp = FALSE)
  }

  ecokit::cat_time(
    "All land use data files downloaded successfully.",
    cat_timestamp = FALSE, verbose = verbose)

  if (sleep_time > 0L) {
    Sys.sleep(sleep_time)
  }

  return(invisible(download_links))

}


# # ********************************************************************** #
# # ********************************************************************** #

# get_sampling_efforts ------

#' Get sampling effort (bias) raster for a taxonomic group
#'
#' Download and prepare a sampling-effort (bias) raster for a specified
#' taxonomic group. If a precomputed bias TIFF already exists in the provided
#' climate directory it will be returned; otherwise the function downloads a 1°
#' bias grid from Davis *et al.* 2023, available at
#' [Zenodo](https://zenodo.org/records/7556851), projects/disaggregates it to
#' the requested resolution using the study mask, masks out invalid cells, names
#' the layer, and writes a compressed GeoTIFF to disk.
#'
#' @param bias_group character scalar. Taxonomic group for which to obtain the
#'   bias raster. Valid values are `"amphibians"`, `"birds"`, `"mammals"`,
#'   `"molluscs"`, `"plants"` and `"reptiles"`.
#' @param resolution integer scalar. Target spatial resolution (see
#'   [OneSDM::get_mask_layer()]) used to disaggregate the 1° bias grid. Defaults
#'   to `10L`.
#' @param climate_dir character. Path to the directory where the output bias
#'   TIFF should be stored (and where existing bias TIFFs are checked).
#' @param return_spatraster logical. If `TRUE` (default) the function returns a
#'   `SpatRaster` object (invisibly). If `FALSE` the function returns
#'   (invisibly) the path to the written TIFF file.
#' @param verbose logical. If `TRUE`, prints progress and informative messages.
#'   Defaults to `FALSE`.
#'
#' @details The function:
#' - constructs a target output filename of the form
#'   `"bias_<bias_group>_res_<resolution>.tif"` in `climate_dir`;
#' - if the file exists, returns it (as a `SpatRaster` or path depending on
#'   `return_spatraster`);
#' - otherwise downloads the appropriate 1° bias raster from Zenodo, projects
#'   it to the study mask resolution via [OneSDM::get_mask_layer()], converts
#'   `NA` values to zeros, applies the mask, sets a descriptive layer name, and
#'   writes a compressed GeoTIFF (ZSTD compression, tiled) to disk.
#'
#' @return Invisibly returns either:
#' - a `SpatRaster` (if `return_spatraster = TRUE`), or
#' - a character scalar with the path to the written TIFF (if
#'   `return_spatraster = FALSE`).
#'
#' @examples
#' \dontrun{
#' # Return a `SpatRaster` for birds at resolution 10
#' r <- get_sampling_efforts(
#'   bias_group = "birds", resolution = 10L, climate_dir = "climate_data")
#' r
#'
#' # Only return path to file after creation
#' path <- get_sampling_efforts(
#'   bias_group = "mammals", resolution = 5L, climate_dir = "climate_data",
#'   return_spatraster = FALSE)
#' path
#'
#' terra::rast(path)
#' }
#'
#' @author Ahmed El-Gabbas
#' @export

get_sampling_efforts <- function(
    bias_group = NULL, resolution = 10L, climate_dir = NULL,
    return_spatraster = TRUE, verbose = FALSE) {

  valid_groups <- c(
    "amphibians", "birds", "mammals", "molluscs", "reptiles", "plants")

  if (is.null(bias_group) || length(bias_group) != 1L || !nzchar(bias_group)) {
    ecokit::stop_ctx(
      "The `bias_group` argument must be a single non-empty string.",
      bias_group = bias_group, cat_timestamp = FALSE)
  }

  if (!(bias_group %in% valid_groups)) {
    ecokit::stop_ctx(
      paste0(
        "Invalid `bias_group` value. Valid options are: ",
        toString(valid_groups), "."),
      bias_group = bias_group, cat_timestamp = FALSE)
  }

  bias_file <- dplyr::case_when(
    bias_group == "amphibians" ~ "amphib_1deg_grid.tif",
    bias_group == "birds" ~ "birds_1deg_grid.tif",
    bias_group == "mammals" ~ "mammals_1deg_grid.tif",
    bias_group == "molluscs" ~ "molluscs_1deg_grid.tif",
    bias_group == "reptiles" ~ "reptiles_1deg_grid.tif",
    bias_group == "plants" ~ "plants_1deg_min5.tif",
    .default = NA_character_)

  out_bias_file <- fs::path(
    climate_dir, paste0("bias_", bias_group, "_res_", resolution, ".tif"))

  if (ecokit::check_tiff(out_bias_file, warning = FALSE)) {
    ecokit::cat_time(
      paste0(
        "Bias file for group ", crayon::blue(bias_group),
        " already exists at: ", crayon::blue(out_bias_file)),
      cat_timestamp = FALSE, verbose = verbose)
    if (return_spatraster) {
      return(invisible(terra::rast(out_bias_file)))
    } else {
      return(invisible(out_bias_file))
    }
  }

  out_bias_temp <- fs::file_temp(
    pattern = paste0("bias_", bias_group, "_"), ext = ".tif")

  mask_r <- OneSDM::get_mask_layer(
    resolution = resolution, climate_dir = climate_dir,
    verbose = verbose, overwrite = FALSE, return_spatraster = TRUE)

  bias_r <- ecokit::zenodo_download_file(
    record_id = "7556851", file = bias_file, dest_file = out_bias_temp,
    read_func = terra::rast, verbose = verbose) %>%
    # disaggregate from 1 degree to current resolution
    terra::project(mask_r) %>%
    terra::classify(cbind(NA, 0L)) %>%
    terra::mask(mask_r) %>%
    stats::setNames(paste0("bias_", bias_group, "_res_", resolution))

  terra::writeRaster(
    x = bias_r, filename = out_bias_file, overwrite = TRUE,
    gdal = c("COMPRESS=ZSTD", "ZSTD_LEVEL=22", "TILED=YES"))

  if (return_spatraster) {
    return(invisible(bias_r))
  } else {
    return(invisible(out_bias_file))
  }

}
