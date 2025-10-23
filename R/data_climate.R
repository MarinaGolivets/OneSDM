#' Prepare CHELSA Climate Data
#'
#' Prepare global terrestrial [CHELSA](https://chelsa-climate.org/) climate data
#' as TIFF files at different spatial resolutions and future projections for use
#' in species distribution modelling.
#'
#' @param agg_factors Integer vector. Aggregation factors for spatial
#'   resolution. Original CHELSA data is at 30 arc-seconds (~1km at the
#'   equator). Aggregation factors are applied to both latitude and longitude
#'   dimensions. Defaults to `c(5L, 10L, 20L)` for 2.5, 5, and 10 arc-minutes
#' @param n_cores Integer. Number of CPU cores to use for parallel processing;
#'   defaults to `20L`.
#' @param min_land_percent Numeric. Minimum percentage of land required in an
#'   aggregated grid cell to be considered as land. Value must be between 0 and
#'   100; defaults to 20L for 20% land.
#' @param path_chelsa_tif Character or `NULL`. Path to directory containing
#'   pre-downloaded CHELSA `.tif` files. If `NULL`, a temporary directory is
#'   used and files are downloaded/removed as needed. If a directory path is
#'   provided, it must exist and contain the original CHELSA `.tif` files.
#' @param temp_dir Character or `NULL`. Path to temporary directory for
#'   intermediate files, including downloaded CHELSA `.tif` files if
#'   `path_chelsa_tif` is `NULL` and temporary files created during processing.
#'   If `NULL`, a new temporary directory is created. The directory and its
#'   contents are deleted after processing.
#' @param predictor_pattern Character. Regular expression pattern that is
#'   compatible with [stringr::str_detect] to filter CHELSA predictor variables.
#'   See CHELSA [Technical
#'   Specifications](https://chelsa-climate.org/downloads/) for a list of
#'   available variables and the `var_name` column of [ecokit::get_chelsa_links]
#'   for the list of supported variables. Defaults to `"bio|npp"` for 19
#'   bioclimatic and net primary productivity variables. Use `".+"` to
#'   process all available variables.
#'
#' @return The function invisibly returns `NULL`. Side effects include writing
#'   processed data to the `climate` sub-directory:
#' - Summary data frame
#'   (`climate/climate_data.RData`) contains metadata for each processed raster.
#' - Tiff files are saved in respective sub-directories:
#'   `1981_2010_res_<agg_factor>` for historic climates (1981-2010) or future
#'   projections with the structure:
#'   `<time_period>_<climate_scenario>_<climate_model>_res_<agg_factor>/`.
#'   - `<time_period>`: `2011_2040`, `2041_2070`, or `2071_2100`.
#'   - `<climate_scenario>`: ssp126, ssp370, or ssp585.
#'   - `<climate_model>`: e.g., `gfdl-esm4`, `ipsl_cm6a_lr`, `mpi_esm1_2_hr`,
#'   `mpi_esm2_0`, and `ukesm1_0_LL`.
#'   - `<agg_factor>`: Aggregation factor (`5`, `10`, or `20`).
#' @details
#' - *This function is not intended to be called by the end user of the `OneSDM`
#' package, but to prepare necessary climate data for use in `OneSDM` species
#' distribution modelling workflow.*
#' - The function loads [land masks](https://doi.org/10.48364/ISIMIP.836809.3/)
#' and prepares reference masks at various resolutions, retrieves CHELSA data
#' links from [ecokit::get_chelsa_links], and processes climate rasters in
#' parallel. It handles downloading, cropping, masking, aggregating, and writing
#' output rasters, as well as attaching metadata.

#'
#' @examples
#' \dontrun{
#'   prepare_climate()
#' }
#'
#' @author Ahmed El-Gabbas
#' @export

prepare_climate <- function(
    agg_factors = c(5L, 10L, 20L), n_cores = 20L, min_land_percent = 20L,
    path_chelsa_tif = NULL, temp_dir = NULL, predictor_pattern = "bio|npp") {

  ecokit::info_chunk("Prepare CHELSA climate data")

  year <- climate_scenario <- climate_model <- climate_name <- projected_map <-
    var_name <- climate_model_abb <- . <- NULL

  # # ********************************************************************** #
  # Setup -----
  # # ********************************************************************** #

  # Temporary directory for intermediate files
  if (is.null(temp_dir)) {
    temp_dir <- fs::path(
      fs::path_temp(),
      paste0("onesdm_chelsa_", as.integer(Sys.time())))
    on.exit(try(fs::dir_delete(temp_dir), silent = TRUE), add = TRUE)
  }
  fs::dir_create(temp_dir)

  terra_options(temp_dir = temp_dir)

  # GDAL settings for writing GeoTIFF files
  # using ZSTD compression with level 22 (maximum compression)
  gdal_settings <- c("COMPRESS=ZSTD", "ZSTD_LEVEL=22", "TILED=YES")

  # # ********************************************************************** #
  # Argument checking -----
  # # ********************************************************************** #

  ecokit::cat_time("Argument checking")

  # # ||||||||||||||||||||||||||||||||||||||||| #
  # agg_factors
  # # ||||||||||||||||||||||||||||||||||||||||| #

  if (is.null(agg_factors) || length(agg_factors) == 0L) {
    ecokit::stop_ctx(
      "Argument `agg_factors` must be a non-empty vector of positive integers.",
      agg_factors = agg_factors, class_agg_factors = class(agg_factors),
      length_agg_factors = length(agg_factors))
  }

  if (!is.numeric(agg_factors) || any(agg_factors < 1L)) {
    ecokit::stop_ctx(
      "Argument `agg_factors` must be a vector of positive integers.",
      agg_factors = agg_factors, class_agg_factors = class(agg_factors))
  }

  # # ||||||||||||||||||||||||||||||||||||||||| #
  # n_cores
  # # ||||||||||||||||||||||||||||||||||||||||| #

  ecokit::check_args(args_to_check = "n_cores", args_type = "numeric")
  if (!is.integer(n_cores)) {
    n_cores <- as.integer(n_cores)
  }
  if (n_cores > parallelly::availableCores()) {
    warning(
      "Argument `n_cores` exceeds the number of available CPU cores.\n",
      "n_cores: ", n_cores, "; available cores: ", parallel::detectCores(),
      call. = FALSE)
    n_cores <- parallelly::availableCores()
  }

  # # ||||||||||||||||||||||||||||||||||||||||| #
  # min_land_percent
  # # ||||||||||||||||||||||||||||||||||||||||| #

  ecokit::check_args(args_to_check = "min_land_percent", args_type = "numeric")
  if (min_land_percent < 0.0 || min_land_percent > 100.0) {
    ecokit::stop_ctx(
      "Argument `min_land_percent` must be a single number between 0 and 100.",
      min_land_percent = min_land_percent)
  }

  # # ||||||||||||||||||||||||||||||||||||||||| #
  # path_chelsa_tif
  # # ||||||||||||||||||||||||||||||||||||||||| #

  if (is.null(path_chelsa_tif)) {
    path_chelsa_tif <- fs::path(temp_dir, "chelsa_raw")
    fs::dir_create(path_chelsa_tif)
  } else if (!is.character(path_chelsa_tif) || length(path_chelsa_tif) != 1L) {
    ecokit::stop_ctx(
      "Argument `path_chelsa_tif` must be `NULL` or a single character string.",
      path_chelsa_tif = path_chelsa_tif,
      class_path_chelsa_tif = class(path_chelsa_tif),
      length_path_chelsa_tif = length(path_chelsa_tif))
  } else if (!(fs::dir_exists(path_chelsa_tif) &&
               fs::is_dir(path_chelsa_tif))) {
    ecokit::stop_ctx(
      "Argument `path_chelsa_tif` must be `NULL` or a valid directory path.",
      path_chelsa_tif = path_chelsa_tif,
      dir_exists = fs::dir_exists(path_chelsa_tif),
      is_dir = fs::is_dir(path_chelsa_tif))
  }

  # # ********************************************************************** #
  # Prepare land masks -----
  # # ********************************************************************** #

  ecokit::cat_time("Prepare land masks")
  path_mask <- "land_mask"
  path_climate <- "climate"
  if (fs::dir_exists(path_mask)) {
    warning(
      "The land_mask directory `",
      crayon::blue(ecokit::normalize_path(path_mask)),
      "` already exists. Existing files may be overwritten.", call. = FALSE)
  }
  if (fs::dir_exists(path_climate)) {
    warning(
      "The climate directory `",
      crayon::blue(ecokit::normalize_path(path_climate)),
      "` already exists. Existing files may be overwritten.", call. = FALSE)
  }
  fs::dir_create(c(path_mask, path_climate))

  # # ||||||||||||||||||||||||||||||||||||||||| #
  # download original land mask
  # # ||||||||||||||||||||||||||||||||||||||||| #

  # file to store original mask as delivered by CHELSA
  file_mask_orig_nc <- fs::path(path_mask, "land_mask_original.nc")

  if (!ecokit::check_tiff(file_mask_orig_nc, warning = FALSE)) {

    ecokit::cat_time("download original land mask", level = 1L)

    mask_orig_url <- stringr::str_glue(
      "https://files.isimip.org/ISIMIP3a/InputData/climate/\\
atmosphere/obsclim/global/daily/historical/CHELSA-W5E5/\\
chelsa-w5e5_obsclim_mask_30arcsec_global.nc")

    file_down <- httr::GET(
      mask_orig_url, httr::write_disk(file_mask_orig_nc, overwrite = TRUE),
      httr::timeout(600L))
    rm(file_down)
  }

  if (!ecokit::check_tiff(file_mask_orig_nc, warning = FALSE)) {
    ecokit::stop_ctx(
      "The original land mask file could not be found or downloaded.",
      file_mask_orig_nc = file_mask_orig_nc,
      file_exists = fs::file_exists(file_mask_orig_nc),
      file_size = if (fs::file_exists(file_mask_orig_nc)) {
        fs::file_size(file_mask_orig_nc)
      } else {
        NA_integer_
      })
  }

  # # ||||||||||||||||||||||||||||||||||||||||| #
  # Crop and save cropped land mask as tiff file
  # # ||||||||||||||||||||||||||||||||||||||||| #

  # file to store original land mask (cropped and classified)
  file_mask_orig <- fs::path(path_mask, "land_mask.tif")

  if (!ecokit::check_tiff(file_mask_orig, warning = FALSE)) {

    ecokit::cat_time("processing original land mask", level = 1L)
    mask_orig <- terra::rast(file_mask_orig_nc) %>%
      ecokit::set_raster_crs("epsg:4326") %>%
      terra::classify(cbind(0L, NA)) %>%
      # exclude Antarctica
      terra::crop(terra::ext(-180L, 180L, -60L, 84L))

    terra::writeRaster(
      x = mask_orig, filename = file_mask_orig, overwrite = TRUE,
      gdal = gdal_settings)

    rm(mask_orig, envir = environment())
  }


  # # ||||||||||||||||||||||||||||||||||||||||| #
  # prepare aggregated land masks
  # # ||||||||||||||||||||||||||||||||||||||||| #

  # files for aggregated land mask(s)
  file_mask_agg <- fs::path(path_mask, paste0("mask_agg_", agg_factors, ".tif"))

  ecokit::cat_time(
    "Prepare aggregated land masks", level = 1L, cat_timestamp = FALSE)

  purrr::walk(
    .x = seq_along(agg_factors),
    .f = ~ {

      file_mask <- file_mask_agg[[.x]]
      agg_factor <- agg_factors[[.x]]

      if (!ecokit::check_tiff(file_mask, warning = FALSE)) {

        ecokit::cat_time(
          paste0("aggregation factor: ", agg_factor),
          level = 2L, cat_timestamp = FALSE)

        # minimum number of land pixels required in the aggregated pixel
        min_land <- (agg_factor * agg_factor) * (min_land_percent / 100L)

        mask_agg <- terra::aggregate(
          x = terra::rast(file_mask_orig), fact = agg_factor, fun = "sum",
          cores = n_cores, na.rm = TRUE) %>%
          terra::classify(cbind(0L, min_land, NA)) %>%
          terra::classify(cbind(min_land, Inf, 1L)) %>%
          stats::setNames(paste0("mask_agg_", agg_factor))

        terra::writeRaster(
          x = mask_agg, filename = file_mask, overwrite = TRUE,
          gdal = gdal_settings)
        rm(mask_agg, envir = environment())
      }
      invisible(gc())
    })

  # # ********************************************************************** #
  # Processing CHELSA data -----
  # # ********************************************************************** #

  # # ||||||||||||||||||||||||||||||||||||||||| #
  # Retrieve CHELSA data links
  # # ||||||||||||||||||||||||||||||||||||||||| #

  ecokit::cat_time("Retrieve CHELSA data links")

  chelsa_links <- ecokit::get_chelsa_links()
  if (!inherits(chelsa_links, "data.frame") || nrow(chelsa_links) == 0L) {
    ecokit::stop_ctx(
      "Error in retrieving CHELSA data links.",
      class_chelsa_links = class(chelsa_links),
      nrow_chelsa_links = nrow(chelsa_links))
  }

  chelsa_links <- dplyr::mutate(
    chelsa_links,
    # extract climate model abbreviations
    climate_model_abb = purrr::map_chr(
      climate_model, stringr::str_extract, "^[A-Za-z0-9]+"),
    climate_name = paste0(year, "_", climate_scenario, "_", climate_model_abb),
    climate_name = stringr::str_replace_all(climate_name, "-", "_"),
    climate_name = stringr::str_replace_all(
      climate_name, "1981_2010_current_current", "1981_2010"),
    climate_name = stringr::str_to_lower(climate_name),
    output_files = purrr::map2(
      .x = climate_name, .y = var_name,
      .f = ~ {
        fs::path(
          path_climate, paste0(.x, "_res_", agg_factors),
          paste0(.y, ".tif")) %>%
          matrix(nrow = 1L) %>%
          as.data.frame() %>%
          setNames(paste0("file_res_", agg_factors)) %>%
          tibble::as_tibble()
      })) %>%
    ecokit::arrange_alphanum(climate_name, var_name) %>%
    # only process bioclimatic (bio1-19) and npp variables (default value)
    dplyr::filter(stringr::str_detect(var_name, predictor_pattern)) %>%
    # Exclude kg0-5 as they are for categorical variables
    dplyr::filter(!stringr::str_detect(var_name, "^kg[0-5]$")) %>%
    tidyr::unnest("output_files")

  # # ||||||||||||||||||||||||||||||||||||||||| #
  # create directories for storing processed climate maps
  # # ||||||||||||||||||||||||||||||||||||||||| #

  ecokit::cat_time("create directories for storing processed climate maps")

  dir_names <- dplyr::select(chelsa_links, dplyr::starts_with("file_res_")) %>%
    unlist() %>%
    dirname() %>%
    unique()
  fs::dir_create(dir_names)

  # # ||||||||||||||||||||||||||||||||||||||||| #
  # Parallel processing of CHELSA data
  # # ||||||||||||||||||||||||||||||||||||||||| #

  ecokit::cat_time("Parallel processing of CHELSA data")

  .start_chelsa_time <- lubridate::now(tzone = "CET")

  if (n_cores == 1L) {
    future::plan("sequential", gc = TRUE)
  } else {
    ecokit::set_parallel(n_cores, show_log = FALSE)
    withr::defer(future::plan("sequential", gc = TRUE))
  }

  par_packages <- c("fs", "terra", "ecokit", "dplyr", "magrittr", "stringr")
  par_globals <- c(
    "path_chelsa_tif", "agg_factors", "file_mask_agg", "terra_options",
    "temp_dir", "chelsa_links", "gdal_settings")

  climate_data <- future.apply::future_lapply(
    X = seq_len(nrow(chelsa_links)),
    FUN = function(map_id) {

      terra_options(temp_dir = temp_dir)

      subset_data <- dplyr::slice(chelsa_links, map_id)

      output_tibble <- dplyr::select(
        subset_data, dplyr::starts_with("file_res_"))

      maps_metatags <- tibble::tibble(
        names = c(
          "variable", "description", "unit", "time",
          "climate_model", "climate_scenario"),
        values = c(
          subset_data$var_name, subset_data$long_name,
          subset_data$unit, subset_data$year,
          stringr::str_replace(subset_data$climate_model, "current", "-"),
          stringr::str_replace(subset_data$climate_scenario, "current", "-")))

      output_exist <- purrr::map_lgl(
        .x = unlist(output_tibble), .f = ecokit::check_tiff, warning = FALSE)

      if (all(output_exist)) {
        return(output_tibble)
      }

      temp_file <- fs::path(
        temp_dir,
        paste0(subset_data$var_name, "_", subset_data$climate_name, ".tif"))

      map_file <-  dplyr::if_else(
        is.null(path_chelsa_tif),
        NA_character_, fs::path(path_chelsa_tif, basename(subset_data$url)))
      map_okay <- !is.na(map_file) &&
        ecokit::check_tiff(map_file, warning = FALSE)

      if (!map_okay) {
        # Check file accessibility with HEAD request
        head_response <- httr::HEAD(subset_data$url)
        if (httr::status_code(head_response) != 200L) {
          return(NULL)
        }
        on.exit(try(fs::file_delete(temp_file), silent = TRUE))
        file_down <- httr::GET(
          subset_data$url, httr::write_disk(temp_file, overwrite = TRUE),
          httr::timeout(600L))
        rm(file_down)
        map_file <- temp_file
      }

      # terra package by default considers the scale and offset information
      # stored in the tiff files. Here I disable this to read the raw values
      # as it is and later consider the scale and offset information
      # manually. This is more safe as I found that some of the future
      # projections do not include such information in the tiff files.

      # For `npp` layers, all tiff maps except for current climate does have
      # a scaling factor all scale and offset information were set manually

      map <- terra::rast(map_file, raw = TRUE) %>%
        # `gsp` maps contains extremely high values instead of NA; the
        # following replace extreme values with NA
        terra::classify(cbind(420000000L, Inf, NA))

      # Manually considering offset and scale
      if (!is.na(subset_data$scale) && subset_data$scale != 1L) {
        map <- map * subset_data$scale
      }
      if (!is.na(subset_data$offset) && subset_data$offset != 0L) {
        map <- map + subset_data$offset
      }

      purrr::walk(
        .x = seq_along(agg_factors),
        .f = function(i) {

          if (ecokit::check_tiff(
            dplyr::pull(output_tibble, i), warning = FALSE)) {
            return(NULL)
          }

          land_mask <- terra::rast(file_mask_agg[[i]])

          projected_map <- terra::aggregate(
            x = map, fact = agg_factors[[i]], na.rm = TRUE) %>%
            terra::crop(land_mask) %>%
            terra::mask(land_mask) %>%
            stats::setNames(subset_data$var_name)

          # set metadata
          maps_metatags0 <- dplyr::bind_rows(
            as.data.frame(maps_metatags),
            data.frame(
              names = "resolution", values = as.character(agg_factors[[i]]),
              stringsAsFactors = FALSE))
          terra::metags(projected_map) <- maps_metatags0

          terra::writeRaster(
            x = projected_map, filename = dplyr::pull(output_tibble, i),
            overwrite = TRUE, gdal = gdal_settings)
          rm(land_mask, projected_map, envir = environment())
          invisible(gc())
        })

      if (fs::file_exists(temp_file)) {
        try(fs::file_delete(temp_file), silent = TRUE)
      }

      rm(map, envir = environment())
      invisible(gc())

      output_tibble
    },
    future.scheduling = Inf, future.seed = TRUE,
    future.packages = par_packages,
    future.globals = par_globals) %>%
    dplyr::bind_rows() %>%
    dplyr::bind_cols(chelsa_links, .) %>%
    tidyr::pivot_longer(
      cols = tidyselect::starts_with("file_res_"),
      names_to = "resolution", values_to = "projected_map",
      names_transform = list(
        resolution = ~ as.integer(stringr::str_remove(.x, "file_res_")))) %>%
    dplyr::mutate(
      file_size_mb = purrr::map_dbl(
        .x = projected_map, .f = ~ round(file.size(.x) / (1024L * 1024L), 2L)),
      projected_exists = purrr::map_lgl(projected_map, fs::file_exists))

  future::plan("sequential", gc = TRUE)

  ecokit::cat_diff(
    init_time = .start_chelsa_time,
    prefix = "\nProcessing CHELSA data was finished in ")

  # # ||||||||||||||||||||||||||||||||||||||||| #
  # Save climate data summary
  # # ||||||||||||||||||||||||||||||||||||||||| #

  save(climate_data, file = fs::path(path_climate, "climate_data.RData"))

  # # ********************************************************************** #

  return(invisible(NULL))

}
