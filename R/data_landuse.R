# prepare_landuse -------

#'Prepare Land-Use Data for `OneSDM` Framework
#'
#'This function downloads, processes, and aggregates global land-use projection;
#'See Chen et al. ([2022](https://doi.org/10.1038/s41597-022-01208-6)) for
#'details. The function handles PFT (Plant Functional Type) data under different
#'Shared Socioeconomic Pathways (SSPs) and Representative Concentration Pathways
#'(RCPs) scenarios. It creates aggregated rasters at multiple resolutions, and
#'generates both majority class maps and percentage coverage maps for plant
#'functional types (PFTs) and their cross-walked groups, and saves the results
#'for future use.
#'
#'@param agg_factors Integer vector. Aggregation factors (resolutions) to apply
#'  when resampling the land-use data. Each value determines a spatial
#'  resolution of output rasters. Must be positive integers. Defaults to c(`5L`,
#'  `10L`, `20L`) for `2.5`, `5`, and `10` arc-minutes resolutions.
#'@param n_cores Integer. Number of CPU cores to use for parallel processing.
#'  Default: 20L.
#'@param temp_dir Character or `NULL`. Path to temporary directory for
#'  intermediate files. If `NULL`, a new temporary directory is created. The
#'  directory and its contents are deleted after processing.
#'
#'@return Invisibly returns `NULL`. The function creates a directory structure
#'  under `"landuse/"` containing:
#'   - Majority class rasters for each climate scenario and resolution: most
#'  common land-use class per grid cell;
#'   - Percentage coverage maps for each of `20` original PFT class;
#'   - Percentage coverage maps for `12` grouped PFT classes (cross-walk, see
#'  below);
#'   - An RData file (`landuse_results.RData`) with metadata tibbles.
#'
#'@details This function is not intended to be called by the end user of the
#'  `OneSDM` package, but to prepare landuse data for use in the workflow.
#'
#'  The function performs the following steps:
#' - Downloads land-use data from [Zenodo](https://10.0.20.161/zenodo.4584775)
#' - Extracts files for climate scenarios (SSP1-RCP26, SSP3-RCP70, SSP5-RCP85)
#'  and time periods (2025-2100) matching CHELSA climate data v2.1. Original
#'  landuse data is available at 30-arcsecond resolution and temporal resolution
#'  of 5 years from 2015 to 2100. To match CHELSA data, the function processes
#'  data for:
#'   - `Current`: 2015
#'   - `2021-2040`: mode of data for 2025, 2030, 2035, 2040
#'   - `2041-2070`: mode of data for 2045, 2050, 2055, 2060, 2065, 2070
#'   - `2071-2100`: mode of data for 2075, 2080, 2085, 2090, 2095, 2100
#' - Creates cross-walk tables to group similar PFT classes
#' - Processes data at multiple resolutions using parallel computation
#' - Generates majority class maps using modal aggregation
#' - Creates percentage coverage maps for each PFT class (original and grouped)
#' - Saves results as TIFF files and a summary `RData` file
#'
#'@note
#' - Requires corresponding mask layers to exist in `"land_mask/"` directory;
#' - If results already exist, the function returns immediately without
#'reprocessing
#' - GDAL compression settings use ZSTD level 22 for efficient storage
#'
#' @section Plant Functional Types and Land Cover Classes Table:
#' This table provides a mapping between unique ID (`ID`) of original PFT classes
#' (`pft`) and the IDs and names of their corresponding cross-walk (`cw_id` and
#' `pft_cw`).
#' | **ID** | **pft**                   | **cw_id** | **pft_cw**               |
#' |----|-----------------------------------|-------|-----------------------------|
#' | 1  | Water                             | 1     | Water                       |
#' | 2  | Broadleaf_evergreen_tree_tropical | 2     | Broadleaf_evergreen_forest  |
#' | 3  | Broadleaf_evergreen_tree_temperate| 2     | Broadleaf_evergreen_forest  |
#' | 4  | Broadleaf_deciduous_tree_tropical | 3     | Broadleaf_deciduous_forest  |
#' | 5  | Broadleaf_deciduous_tree_temperate| 3     | Broadleaf_deciduous_forest  |
#' | 6  | Broadleaf_deciduous_tree_boreal   | 3     | Broadleaf_deciduous_forest  |
#' | 7  | Needleleaf_evergreen_tree_temperate|4     | Needleleaf_evergreen_forest |
#' | 8  | Needleleaf_evergreen_tree_boreal  | 4     | Needleleaf_evergreen_forest |
#' | 9  | Needleleaf_deciduous_tree         | 5     | Needleleaf_deciduous_forest |
#' | 10 | Broadleaf_evergreen_shrub, temperate|6    | Broadleaf_evergreen_shrubland|
#' | 11 | Broadleaf_deciduous_shrub, temperate|7    | Broadleaf_deciduous_shrubland|
#' | 12 | Broadleaf_deciduous_shrub, boreal | 7     | Broadleaf_deciduous_shrubland|
#' | 13 | C3_grass, arctic                  | 8     | Grasslands                  |
#' | 14 | C3_grass                          | 8     | Grasslands                  |
#' | 15 | C4_grass                          | 8     | Grasslands                  |
#' | 16 | Mixed_c3_c4_grass                 | 8     | Grasslands                  |
#' | 17 | Barren                            | 9     | Barren                      |
#' | 18 | Cropland                          | 10    | Cropland                    |
#' | 19 | Urban                             | 11    | Urban                       |
#' | 20 | Permanent_snow_and_ice            | 12    | Permanent_snow_and_ice      |
#'
#'@references
#' - Chen G, Li X, Liu X (2022). Global land projection based on plant
#'functional types with a 1-km resolution under socio-climatic scenarios.
#'Scientific Data, 9, 125. <https://doi.org/10.1038/s41597-022-01208-6>
#' @examples
#' \dontrun{
#'   prepare_landuse()
#' }
#'
#'@author Ahmed El-Gabbas
#'@export

prepare_landuse <- function(
    agg_factors = c(5L, 10L, 20L), n_cores = 20L, temp_dir = NULL) {

  cw_id <- ID <- resolution <- year <- climate_scenario <- cw_id <- #nolint
    climate_name <- path <- ids <- pft <- pft_cw <- NULL

  gdal_settings <- c("COMPRESS=ZSTD", "ZSTD_LEVEL=22", "TILED=YES")

  ecokit::check_packages(
    c(
      "archive", "fs", "future", "future.apply", "parallelly", "purrr",
      "stringr", "terra", "tibble", "tidyr", "tidyselect", "withr"))

  # # ********************************************************************** #
  # Argument checking -----
  # # ********************************************************************** #

  ecokit::cat_time("Argument checking")

  # # ||||||||||||||||||||||||||||||||||||||||| #
  # temp_dir
  # # ||||||||||||||||||||||||||||||||||||||||| #

  # Temporary directory for intermediate files
  if (is.null(temp_dir)) {
    temp_dir <- fs::path(
      fs::path_temp(), paste0("onesdm_landuse_", as.integer(Sys.time())))
    on.exit(try(fs::dir_delete(temp_dir), silent = TRUE), add = TRUE)
  }
  fs::dir_create(temp_dir)

  terra_options(temp_dir = temp_dir)

  # # ||||||||||||||||||||||||||||||||||||||||| #
  # agg_factors
  # # ||||||||||||||||||||||||||||||||||||||||| #

  if (!is.numeric(agg_factors) || any(agg_factors < 1L) ||
      any(agg_factors != as.integer(agg_factors))) {
    ecokit::stop_ctx(
      "Argument `agg_factors` must be a vector of positive integers.",
      agg_factors = agg_factors, class_agg_factors = class(agg_factors))
  }
  agg_factors <- as.integer(agg_factors)

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
      "n_cores: ", n_cores, "; available cores: ", parallelly::availableCores(),
      call. = FALSE)
    n_cores <- parallelly::availableCores()
  }

  # # ********************************************************************** #
  # Check if results already exist -----
  # # ********************************************************************** #
  path_landuse <- fs::path("landuse")
  fs::dir_create(path_landuse)

  file_landuse_results <- fs::path(path_landuse, "landuse_results.RData")
  if (ecokit::check_data(file_landuse_results, warning = FALSE)) {
    return(ecokit::load_as(file_landuse_results))
  }

  # # ********************************************************************** #
  # Download land-use data from Zenodo -----
  # # ********************************************************************** #

  ecokit::cat_time("Download land-use data from Zenodo")

  file_landuse_raw <- fs::path(path_landuse, "landuse.zip")

  if (!ecokit::check_zip(file_landuse_raw, warning = FALSE)) {

    # Zenodo record: https://zenodo.org/record/4584775
    landuse_file <- paste0(
      "Global PFT-based land projection dataset under SSPs-RCPs.zip")
    landuse_file <- ecokit::zenodo_download_file(
      "4584775", file_name = landuse_file, dest_file = file_landuse_raw,
      timeout = 1800L)

    ecokit::cat_time(
      "Checking that downloaded land-use zip file is valid", level = 1L)
    if (!ecokit::check_zip(file_landuse_raw, warning = FALSE)) {
      ecokit::stop_ctx(
        "Land-use zip file does not exist or is not a valid zip file.",
        file_landuse_raw = file_landuse_raw)
    }
  }

  # # ********************************************************************** #
  # Extract specific files from the zip file -----
  # # ********************************************************************** #

  ecokit::cat_time("Extract specific files from the zip file")

  # File name pattern to match
  year_regex <- paste0(
    "(_", paste(seq(2025L, 2100L, by = 5L), collapse = "|_"), ")")
  climate_scenarios <- c("SSP1_RCP26", "SSP3_RCP70", "SSP5_RCP85")
  regex_to_match <- paste0(
    "^global_PFT_2015.+|(",
    paste0("^", climate_scenarios, collapse = ".+|"),
    ").+", year_regex, "\\.+")

  dir_landuse_extract <- fs::path(temp_dir, "landuse_extract")
  fs::dir_create(dir_landuse_extract)

  # List all files in the zip that match the pattern
  landuse_data <- archive::archive(file_landuse_raw) %>%
    dplyr::select(path) %>%
    dplyr::filter(stringr::str_detect(path, regex_to_match)) %>%
    dplyr::mutate(
      file_path = fs::path(dir_landuse_extract, path),
      # Match with year ranges used in CHELSA climate data
      year = purrr::map_chr(
        .x = path,
        .f = ~ {
          base_name <- basename(.x)
          year <- stringr::str_extract(base_name, "\\d{4}")
          future_time_1 <- c("2025", "2030", "2035", "2040")
          future_time_2 <- c("2045", "2050", "2055", "2060", "2065", "2070")
          future_time_3 <- c("2075", "2080", "2085", "2090", "2095", "2100")

          dplyr::case_when(
            year == "2015" ~ "current",
            year %in% future_time_1 ~ "2021_2040",
            year %in% future_time_2 ~ "2041_2070",
            year %in% future_time_3 ~ "2071_2100",
            .default = NA_character_)
        }),
      # Match with climate scenario names
      climate_scenario = purrr::map_chr(
        .x = path,
        .f = ~ {
          base_name <- basename(.x)
          dplyr::case_when(
            stringr::str_detect(base_name, "2015") ~ "current",
            stringr::str_detect(base_name, "SSP1_RCP26") ~ "ssp126",
            stringr::str_detect(base_name, "SSP3_RCP70") ~ "ssp370",
            stringr::str_detect(base_name, "SSP5_RCP85") ~ "ssp585",
            .default = NA_character_)
        })
    )

  if (nrow(landuse_data) == 0L) {
    ecokit::stop_ctx("No matched land-use files were found in the zip file. ")
  }

  # Extract only matching files; suppress print progress
  suppressMessages(
    suppressWarnings(
      archive::archive_extract(
        archive = file_landuse_raw, dir = dir_landuse_extract,
        files = landuse_data$path)
    ))

  # # ********************************************************************** #
  # PFT cross-walks table -------
  # # ********************************************************************** #

  classes_tbl <- tibble::tribble(
    ~ID, ~pft, ~cw_id, ~pft_cw,
    1L, "Water", 1L, "Water",
    2L, "Broadleaf_evergreen_tree_tropical", 2L, "Broadleaf_evergreen_forest",
    3L, "Broadleaf_evergreen_tree_temperate", 2L, "Broadleaf_evergreen_forest",
    4L, "Broadleaf_deciduous_tree_tropical", 3L, "Broadleaf_deciduous_forest",
    5L, "Broadleaf_deciduous_tree_temperate", 3L, "Broadleaf_deciduous_forest",
    6L, "Broadleaf_deciduous_tree_boreal", 3L, "Broadleaf_deciduous_forest",
    7L, "Needleleaf_evergreen_tree_temperate", 4L,
    "Needleleaf_evergreen_forest",
    8L, "Needleleaf_evergreen_tree_boreal", 4L, "Needleleaf_evergreen_forest",
    9L, "Needleleaf_deciduous_tree", 5L, "Needleleaf_deciduous_forest",
    10L, "Broadleaf_evergreen_shrub, temperate", 6L,
    "Broadleaf_evergreen_shrubland",
    11L, "Broadleaf_deciduous_shrub, temperate", 7L,
    "Broadleaf_deciduous_shrubland",
    12L, "Broadleaf_deciduous_shrub, boreal", 7L,
    "Broadleaf_deciduous_shrubland",
    13L, "C3_grass, arctic", 8L, "Grasslands",
    14L, "C3_grass", 8L, "Grasslands",
    15L, "C4_grass", 8L, "Grasslands",
    16L, "Mixed_c3_c4_grass", 8L, "Grasslands",
    17L, "Barren", 9L, "Barren",
    18L, "Cropland", 10L, "Cropland",
    19L, "Urban", 11L, "Urban",
    20L, "Permanent_snow_and_ice", 12L, "Permanent_snow_and_ice")

  # # ********************************************************************** #
  # Processing land-use data at different mask grids -----
  # # ********************************************************************** #

  landuse_data <- landuse_data %>%
    tidyr::nest(file_info = tidyselect::all_of(c("path", "file_path"))) %>%
    dplyr::mutate(
      climate_name = paste0(climate_scenario, "_", year),
      climate_name = stringr::str_replace_all(
        climate_name, "current_current", "current")) %>%
    tidyr::expand_grid(agg_factor = agg_factors)

  if (n_cores == 1L) {
    future::plan("sequential", gc = TRUE)
  } else {
    ecokit::set_parallel(min(n_cores, nrow(landuse_data)), show_log = FALSE)
    withr::defer(future::plan("sequential", gc = TRUE))
  }

  par_packages <- c(
    "fs", "terra", "ecokit", "dplyr", "magrittr", "purrr", "tibble",
    "stringr", "tidyselect")

  landuse_data0 <- future.apply::future_lapply(
    X = seq_len(nrow(landuse_data)),
    FUN = function(map_id) {

      terra_options(temp_dir = temp_dir)

      # # ||||||||||||||||||||||||||||||||||||||||| #
      # Load the original raster and add metadata
      # # ||||||||||||||||||||||||||||||||||||||||| #

      landuse_sub <- dplyr::slice(landuse_data, map_id)
      agg_factor <- landuse_sub$agg_factor

      orig <- landuse_sub$file_info[[1L]]$file_path %>%
        stringr::str_subset(".tif$") %>%
        terra::rast()
      if (terra::nlyr(orig) > 1L) {
        orig <- terra::app(orig, "modal")
      }
      orig <- terra::as.factor(orig)

      rat <- terra::levels(orig)[[1L]] %>%
        dplyr::left_join(classes_tbl, by = "ID") %>%
        dplyr::select(tidyselect::all_of(c("ID", "pft", "pft_cw")))
      levels(orig)[[1L]] <- rat
      # Set the active attribute table to the one containing the cross-walked
      # PFT classes
      terra::activeCat(orig) <- 2L

      # # ||||||||||||||||||||||||||||||||||||||||| #
      # Create output directory for this climate and aggregation factor
      # # ||||||||||||||||||||||||||||||||||||||||| #

      dir_out <- fs::path(
        path_landuse, paste0(landuse_sub$climate_name, "_res_", agg_factor))
      fs::dir_create(dir_out)

      # # |||||||||||||||||||||||||||||||||||||||||||||||||| #
      # Load the corresponding mask layer
      # # |||||||||||||||||||||||||||||||||||||||||||||||||| #

      map_mask <- fs::path("land_mask", paste0("mask_agg_", agg_factor, ".tif"))
      if (!ecokit::check_tiff(map_mask, warning = FALSE)) {
        ecokit::stop_ctx(
          "Mask layer file does not exist", map_mask = map_mask)
      }
      map_mask <- terra::rast(map_mask)

      # # |||||||||||||||||||||||||||||||||||||||||||||||||| #
      # Project to get the most common class (mode), preserving metadata
      # # |||||||||||||||||||||||||||||||||||||||||||||||||| #

      path_majority <- fs::path(dir_out, "landuse_majority.tif")
      majority_tibble <- tibble::tibble(
        resolution = as.integer(agg_factor),
        climate_scenario = landuse_sub$climate_scenario,
        year = landuse_sub$year, climate_name = landuse_sub$climate_name,
        file_majority = path_majority)

      if (!ecokit::check_tiff(path_majority, warning = FALSE)) {
        map_majority <- terra::project(
          x = orig, y = map_mask, method = "mode", threads = TRUE) %>%
          terra::mask(map_mask) %>%
          stats::setNames(
            paste0("majority_", landuse_sub$climate_name, "_res_", agg_factor))
        invisible(gc())
        terra::writeRaster(
          x = map_majority, filename = path_majority, overwrite = TRUE,
          gdal = gdal_settings)
        rm(map_majority, envir = environment())
        invisible(gc())
      }

      # # |||||||||||||||||||||||||||||||||||||||||||||||||| #
      # Percentage coverage maps for each original ptf classes
      # # |||||||||||||||||||||||||||||||||||||||||||||||||| #

      percent_pft <- purrr::map_dfr(
        .x = seq_len(nrow(rat)),
        .f = ~ {

          file_percent <- fs::path(
            dir_out,
            paste0("landuse_pft_", .x, "_", tolower(rat$pft[[.x]]), ".tif"))
          out_tibble <- tibble::tibble(
            resolution = as.integer(agg_factor),
            climate_scenario = landuse_sub$climate_scenario,
            year = landuse_sub$year, climate_name = landuse_sub$climate_name,
            ID = as.integer(rat$ID[[.x]]),
            pft = rat$pft[[.x]], file_pft = file_percent)

          if (ecokit::check_tiff(file_percent, warning = FALSE)) {
            return(out_tibble)
          }

          prop_map <- terra::project(
            x = (orig == .x), y = map_mask,
            method = "average", threads = TRUE) %>%
            terra::mask(map_mask) %>%
            magrittr::multiply_by(100L) %>%
            stats::setNames(paste0("pft_", .x, "_", tolower(rat$pft[[.x]])))

          terra::writeRaster(
            x = prop_map, filename = file_percent, overwrite = TRUE,
            gdal = gdal_settings)

          rm(prop_map, envir = environment())
          invisible(gc())
          out_tibble
        })

      # # |||||||||||||||||||||||||||||||||||||||||||||||||| #
      # Percentage coverage maps for ptf cross-walk
      # # |||||||||||||||||||||||||||||||||||||||||||||||||| #

      percent_groups <- dplyr::group_by(classes_tbl, pft_cw) %>%
        dplyr::select(-pft) %>%
        dplyr::summarise(ids = list(ID)) %>%
        dplyr::pull(ids, name = pft_cw)

      percent_pft_cw <- purrr::map_dfr(
        .x = seq_len(length(percent_groups)),
        .f = ~ {

          cw_name <- names(percent_groups)[[.x]]
          cw_id <- unique(dplyr::filter(classes_tbl, pft_cw == cw_name)$cw_id)

          file_percent <- fs::path(
            dir_out,
            paste0("landuse_pft_cw_", cw_id, "_", tolower(cw_name), ".tif"))

          out_tibble <- tibble::tibble(
            resolution = as.integer(agg_factor),
            climate_scenario = landuse_sub$climate_scenario,
            year = landuse_sub$year, climate_name = landuse_sub$climate_name,
            cw_id = as.integer(cw_id),
            pft_cw = cw_name, file_pft_cw = file_percent)

          if (ecokit::check_tiff(file_percent, warning = FALSE)) {
            return(out_tibble)
          }

          prop_map <- percent_pft %>%
            dplyr::filter(ID %in% percent_groups[[.x]]) %>%
            dplyr::pull(file_pft) %>%
            purrr::map(terra::rast) %>%
            terra::rast() %>%
            terra::app(fun = sum, na.rm = TRUE) %>%
            stats::setNames(paste0("pft_cw_", cw_id, "_", tolower(cw_name)))

          terra::writeRaster(
            x = prop_map, filename = file_percent, overwrite = TRUE,
            gdal = gdal_settings)

          rm(prop_map, envir = environment())
          invisible(gc())

          out_tibble
        })

      # # |||||||||||||||||||||||||||||||||||||||||||||||||| #

      # Return a tibble with all output file paths
      tibble::tibble(
        majority = list(majority_tibble),
        percent_pft = list(percent_pft),
        percent_pft_cw = list(percent_pft_cw))

    },
    future.scheduling = Inf, future.seed = TRUE, future.packages = par_packages,
    future.globals = c(
      "classes_tbl", "landuse_data", "agg_factors", "temp_dir",
      "path_landuse", "gdal_settings", "terra_options")) %>%
    dplyr::bind_rows()

  ecokit::set_parallel(stop_cluster = TRUE, show_log = FALSE)

  landuse_majority <- dplyr::bind_rows(landuse_data0$majority) %>%
    dplyr::arrange(year, climate_scenario, resolution)
  landuse_percent_pft <- dplyr::bind_rows(landuse_data0$percent_pft) %>%
    dplyr::arrange(year, climate_scenario, resolution, ID)
  landuse_percent_pft_cw <- dplyr::bind_rows(landuse_data0$percent_pft_cw) %>%
    dplyr::arrange(year, climate_scenario, resolution, cw_id)

  landuse_results <- list(
    majority = landuse_majority,
    percent_pft = landuse_percent_pft,
    percent_pft_cw = landuse_percent_pft_cw)

  save(landuse_results, file = file_landuse_results)

  # # ********************************************************************** #

  return(invisible(NULL))

}
