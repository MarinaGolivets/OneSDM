# prepare_model_data -------

#' Prepare Modelling Data for SDM Modelling using OneSDM
#'
#' High-level workflow to prepare presence/pseudo-absence, predictors and
#' cross-validation datasets required for fitting species distribution models
#' using the `OneSDM` workflow. This function orchestrates species data
#' preparation, predictor assembly (climate and land use data under current and
#' optional future scenarios), variable selection (using the variance inflation
#' factor, VIF), spatial block creation for cross-validation, pseudo-absence
#' sampling, and persistent saving of intermediate raster/data objects to disk.
#'
#' @param model_dir Character (required). Path to the modelling directory where
#'   data and fitted models will be saved. This should be a single directory for
#'   all models of a given species (a separate directory is expected in case of
#'   modelling multiple species, otherwise data files may be overwritten or
#'   mixed up). This can also be set via the "`onesdm_model_dir`" option.
#' - The function creates a subdirectory "`data`" to store processed species
#'   data.
#' - The function creates a subdirectory "`models_res_<resolution>`" to store
#'   model data and outputs for the specified resolution.
#' @param climate_dir Character (required). Directory path for climate data
#'   files. Can be set via the "`onesdm_climate_dir`" option. This directory is
#'   used to load/download mask layers matching the specified resolution and
#'   current/future climate and land use data (see below). The same directory
#'   should be used in case of modelling multiple species to ensure consistency
#'   and redundant re-download of large geospatial data.
#' @param resolution Integer. Spatial resolution used to prepare data for
#'   analysis and model fitting. valid values are 5, 10 (default), or 20,
#'   corresponding to approximate spatial resolutions of 5, 10, and 20 km (2.5,
#'   5, and 10 arc-minutes) respectively. This value can also be set via the
#'   "`onesdm_resolution`" option.
#' @param climate_scenario,climate_model,climate_year Character. Future climate
#'   scenario(s), model(s), and year(s). These three parameters control which
#'   future climate data will be prepared for future predictions (data for
#'   selected variables at current climate conditions are always prepared):
#'   - **`climate_scenario`**:
#'   Shared Socioeconomic Pathways (SSPs). Valid values are: `"ssp126"`,
#'   `"ssp370"`, `"ssp585"`, `"all"` (default), or `"none"`. This can also be
#'   set via the "`onesdm_climate_scenario`" option.
#'   - **`climate_model`**: abbreviation for future climate
#'   model(s). Valid values are: `"gfdl"`, `"ipsl"`, `"mpi"`, `"mri"`,
#'   `"ukesm1"`, `"all"` (default), or `"none"`. This can also be set via the
#'   "`onesdm_climate_model`" option.
#'   - **`climate_year`**: future year range(s). Valid values are:
#'   `"2011_2040"`, `"2041_2070"`, `"2071_2100"`, `"all"` (default), `"none"`.
#'   This can also be set via the "`onesdm_climate_year`" option.
#'
#'   Multiple values are possible for each of the three parameters, each
#'   provided as a character vector. Combinations of scenario × model × year
#'   will be prepared. This is unless any of the three parameters is set to
#'   `"none"` or `NULL`, in which case no future climate data will be prepared.
#'   For more details, see [OneSDM::climate_data] and
#'   [OneSDM::get_climate_data].
#'
#' @param pft_type,pft_id Information on the plant functional type (PFT)
#'   parameters for land-use predictors. See [OneSDM::landuse_data] and
#'   [OneSDM::get_landuse_data] for details.
#'   - **`pft_type`**: Character vector of length 1. Plant functional type
#'   category to download land-use predictors for. Must be one of "`cross-walk`"
#'   (default) or "`original`". This can be set via the "`onesdm_pft_type`"
#'   option.
#'   - **`pft_id`**: Numeric vector. One or more plant functional type
#'   identifiers to download. Must be valid for the specified `pft_type`: 1-20
#'   for `pft_type == "original"`, and 1-12 for `pft_type == "cross-walk"`. If
#'   `NULL` (default), no land-use predictors are used. This can be set via the
#'   "`onesdm_pft_id`" option.
#'
#' @param bias_group Character scalar. Taxonomic group used to compute sampling
#'   effort surface. Valid options are `"amphibians"`, `"birds"`, `"mammals"`,
#'   `"molluscs"`, `"plants"`, `"reptiles"`, or `"none"`. Default is `"none"`;
#'   i.e., no sampling-effort predictor is used. If not "none", the function
#'   downloads sampling-effort raster of the specified group from
#'   [Zenodo](https://zenodo.org/records/7556851). See
#'   [OneSDM::get_sampling_efforts] for details. This can also be set via the
#'   "`onesdm_bias_group`" option.
#'   - When not "none", the sampling-effort raster
#'   (log<sub>10</sub> scale) is used as a predictor in the models. The
#'   90<sup>th</sup>-percentile value of the  sampling-effort raster in the
#'   modelling study area is used across the whole study area to correct for
#'   sampling bias in projections (but not during models evaluation). For more
#'   details on model-based sampling bias correction, see Warton et al.
#'   ([2013](https://doi.org/10.1371/journal.pone.0079168)).
#'   - The `bias` predictor, if used, is enforced to be kept during VIF
#'   screening (see below).
#' @param vif_th Numeric. Variance inflation factor (VIF) threshold for
#'   multicollinearity screening. Predictors with VIF larger than this value are
#'   excluded. Default is `10`.
#'
#'   - This argument refers to the `th` parameter of the [usdm::vifstep]
#'   function. This can also be set via the "`onesdm_vif_th`" option.
#'   - The function checks if the option "`onesdm_vif_sample`" is set in the
#'   current R session; if so, its value (numeric integer) is used as the
#'   `sample` argument of [usdm::vifstep]. This option controls the number of
#'   random sample points drawn from the study area raster stack. If not set,
#'   the default value is used (5000).
#'   - Similarly, if the option "`onesdm_vif_keep`" is set in the current R
#'   session, its value (a character vector of predictor names) is used as the
#'   `exclude` argument of [usdm::vifstep]. This option allows specifying
#'   predictors that should always be retained regardless of their VIF values.
#'   If not set, no predictors are forced to be kept. If a valid `bias_group` is
#'   provided, the corresponding `bias` predictor is always included in this
#'   list to ensure it is retained during VIF screening.
#' @param pres_buffer Numeric or `NULL`. Distance in kilometres. Non-presence
#'   grid cells beyond this distance from any presence grid cells are masked out
#'   of the modelling area (i.e., no pseudo-absences will be sampled there).
#'   This helps to exclude areas that are unlikely to be accessible to the
#'   species or are too environmentally distant from the known species
#'   presences. Larger values lead to more inclusive modelling areas. If `NULL`
#'   or \eqn{\leq 0}, no presences-based masking is applied (i.e.,
#'   pseudo-absences can be sampled anywhere in the study area). Default is
#'   `1000` km. This can be set via the "`onesdm_pres_buffer`" option.
#' @param abs_buffer Numeric or `NULL`. Minimum distance in kilometres from
#'   presence points inside which pseudo-absences will not be sampled. This
#'   helps avoid sampling pseudo-absences in nearby areas of presences that are
#'   likely to be suitable (e.g., due to spatial autocorrelation) or are already
#'   occupied by the species. Larger values lead to more conservative
#'   pseudo-absence sampling. If `NULL` or \eqn{\leq 0}, no presence-based
#'   pseudo-absence exclusion buffer is applied (i.e., pseudo-absences can be
#'   sampled anywhere outside presence grid cells, prior to sampling
#'   pseudo-absences using [sdm::background]). Default is `10` km. This can be
#'   set via the "`onesdm_abs_buffer`" option.
#' @param abs_exclude_ext List (Optional). A list of terra `SpatExtent` objects,
#'   typically generated using [terra::ext]. Each extent defines a geographic
#'   area where pseudo-absences will not be sampled. This is useful for
#'   excluding regions such as areas outside the species' native range or
#'   regions with specific environmental conditions that preclude the species'
#'   presence. If an empty list (default), no additional exclusion extents are
#'   applied. This can also be set via the "`onesdm_abs_exclude_ext`" option.
#' @param min_pres_grids Integer. Minimum number of presence grid cells required
#'   to proceed with modelling. If the number of presence grid cells after
#'   filtering is less than this value, the function will stop with an error.
#'   Default is `50`. This can also be set via "`onesdm_min_pres_grids`" option.
#' @param cv_folds,cv_block_size,cv_random Parameters controlling spatial-block
#'   cross-validation.
#'   - **`cv_folds`**: Integer. Number of cross-validation folds. Default is
#'   `5`. This can also be set via the "`onesdm_cv_folds`" option.
#'   - **`cv_block_size`**: Numeric. Approximate cross-validation block
#'   size in kilometres. Default is `500` km. Larger block sizes lead to fewer,
#'   more spatially distinct blocks. This can also be set via the
#'   "`onesdm_cv_block_size`" option.
#'   - **`cv_random`**: Logical. If `TRUE` (default), create random spatial
#'   blocks by aggregating and randomly assigning blocks into `cv_folds` groups.
#'   Although this is done randomly, the function tries to ensure that each fold
#'   has a balanced number of presence grid cells. If `FALSE`,
#'   [blockCV::cv_spatial] is used for block determination. Note that
#'   [blockCV::cv_spatial] can be computationally intensive for large study
#'   areas and small block sizes. This can also be set via the
#'   "`onesdm_cv_random`" option.
#' @param abs_ratio,model_n_reps,n_cores Parameters for sampling
#'   pseudo-absences.
#'   - **`abs_ratio`**: Integer scalar. Desired ratio of pseudo-absences to
#'   training presences used when sampling pseudo-absences grid cells for each
#'   cross-validation fold. Default is `10`, which means that for each training
#'   presence grid cell, 10 pseudo-absence grid cells will be sampled. This can
#'   also be set via the "`onesdm_abs_ratio`" option.
#'   - **`model_n_reps`**: Integer scalar. Number of repetition datasets
#'   (different sets of pseudo-absences) to generate for each cross-validation
#'   fold. Pseudo-absences are sampled using [sdm::background] with method
#'   `"eDist"`, which draws pseudo-absences weighted by environmental distance
#'   from presence grid cells. Default is `5`, for five different pseudo-absence
#'   sample datasets per fold. If the available non-presence grid cells become
#'   limiting (when the number of requested pseudo-absences exceeds 80% of
#'   available non-presence grid cells), the function will reduce the number of
#'   repetitions to `1` as further repetitions would be redundant. This can also
#'   be set via the "`onesdm_model_n_reps`" option.
#'   - **`n_cores`**: Integer. Number of CPU cores to use for parallel
#'   processing of data for cross-validation folds. Default is equal to
#'   `cv_folds`. Parallel processing helps to speed up the preparation of
#'   cross-validation modelling data; however, the `sdm::background` can be
#'   memory-intensive when sampling many pseudo-absences over large study areas,
#'   so care should be taken when setting this value. This can also be set via
#'   the "`onesdm_model_n_cores`" option.
#'
#' @inheritParams prepare_species_data
#' @inheritParams get_climate_data
#'
#' @details The function has several side-effects: it writes GeoTIFFs and RData
#'   files under a subdirectory of "`model_dir`" (folder
#'   "`models_res_<resolution>`"), and download/process climate and land-use
#'   data via OneSDM utilities, if needed.
#'
#'   The function performs the following major steps:
#'
#'   - Validates input arguments and resolves defaults from package options.
#'   - Calls [OneSDM::prepare_species_data] to check/prepare species
#'   distribution data. This includes loading occurrence data from EASIN/GBIF or
#'   user-provided coordinates.
#'   - Masks and filters the study area using optional exclusion extents and
#'   distances (`abs_exclude_ext`, `pres_buffer`, and `abs_buffer`).
#'   -  Loads/download current climate (and optionally land-use) data for
#'   selected predictors via [OneSDM::get_climate_data] and
#'   [OneSDM::get_landuse_data], combines them and masks them to the study area.
#'   - Loads/download future climate (and optionally land-use) predictors for
#'   each requested scenario/model/year combination (if any) via
#'   [OneSDM::get_climate_data] and [OneSDM::get_landuse_data], and masks them
#'   to the study area.
#'   - Optionally adds a sampling-effort (bias) predictor when `bias_group`
#'   is valid and != "none".
#'   - Performs VIF-based predictor selection (using [usdm::vifstep]) and
#'   excludes highly collinear variables. The `bias` predictor, if used, is
#'   always retained.
#'   - Creates spatial cross-validation blocks either with a simple random-block
#'   aggregation or via [blockCV::cv_spatial] (controlled by `cv_random`).
#'   - Builds a modelling table (raster + predictors), saves wrapped rasters and
#'   per-fold data to disk and samples pseudo-absences per cross-validation
#'   fold/repetition.
#'   - Saves a tibble describing all cross-validation folds and model
#'   repetitions with file paths to persisted datasets (training/testing
#'   rasters and pseudo-absence objects).
#'
#'   The function also saves a summary list for the parameters used and file
#'   paths to important intermediate data objects under
#'   "`<model_dir>/models_res_<resolution>/model_data_summary.RData`". This list
#'   contains the following elements:
#'   - **`model_dir`**: Character. The modelling directory path.
#'   - **`species_data_raw`**: list of GBIF and EASIN IDs and/or user-provided
#'   coordinates used. This also includes the file path to the saved GeoTIFF
#'   file for prepared species  distribution `species_data_tiff`.
#'   - **`climate`**: list. Information on the current and future climate data
#'   used, including climate directory path, climate variable names, scenarios,
#'   models, years, and file paths to the downloaded GeoTIFF files.
#'   - **`landuse`**: list. Information on the land-use data used, including PFT
#'   type, PFT IDs, and file paths to the downloaded GeoTIFF files.
#'   - **`bias`**: list. Information on the sampling-effort (bias) predictor
#'   used, if any, including the taxonomic group and the value of the
#'   90<sup>th</sup>-percentile used for bias correction `bias_fix_value`.
#'   - **`var_selection`**: list. Information on the variable selection process,
#'   including the VIF threshold used, names of always kept or excluded
#'   predictors and file path for the VIF selection summary.
#'   - **`model_options`**: list. Information on modelling options used,
#'   including resolution, presence filtering buffer, absence exclusion buffer,
#'   absence-to-presence ratio, minimum number of presence grid cells,
#'   cross-validation parameters, and number of model repetitions.
#'   - **`model_data`**: list. Information on the prepared modelling data,
#'   including file paths to the saved presence/pseudo-absence GeoTIFF files,
#'   predictor names, number of presence grid cells, and file paths to the saved
#'   model predictors and modelling data RData files. This also includes the
#'   tibble `model_data_cv` with one row per cross-validation fold × repetition,
#'   saved as an RData file (see the `value` section below).
#'
#' @return Invisibly returns a tibble (`model_data_cv`) with one row per
#'   cross-validation fold × repetition. The tibble is saved as `.RData` file
#'   under "`<model_dir>/models_res_<resolution>/model_data_cv.RData`". The
#'   tibble consists of the following columns.
#'   - **`cv`**: Integer. Cross-validation fold index.
#'   - **`training_data`/`training_data_r`**: Character. File path to `.RData`
#'   files for training data table and raster.
#'   - **`training_pres`/`training_abs`**: Character. File path to `.RData` file
#'   containing training presence/non-presence data tables.
#'   - **`n_training_pres` / `n_training_abs`**: Integer. Number of
#'   presence/non-presence grid cells in the training data.
#'   - **`n_model_reps`**: Integer. Number of model repetitions (pseudo-absence
#'   datasets) used for this fold.
#'   - **`n_pseudo_abs`**: Integer. Number of pseudo-absences sampled per
#'   fold/rep.
#'   - **`model_rep_id`** / **`pseudo_abs_files`**: Integer/Character. Model
#'   repetition ID (from 1 to `n_model_reps`) and file paths to the sampled
#'   pseudo-absence objects for this repetition.
#'   - **`testing_data`/`testing_data_r`**: Character. File path to `.RData`
#'   files for testing data table and raster.
#'   - **`n_test_pres` / `n_test_abs`**: Integer. Number of presence/
#'   non-presence grid cells in the testing data.
#'
#' @examples
#' \dontrun{
#'
#'   model_cv_tbl <- prepare_model_data(
#'     model_dir = "test/Acacia_karroo",
#'     climate_dir = "test/climate_data",
#'     easin_ids = "R00042",
#'     gbif_ids = 2979128,
#'     coordinates = NULL,
#'     resolution = 10L,
#'     var_names = c("bio1","bio12","bio15"),
#'     pft_id = NULL,
#'     bias_group = "plants",
#'     cv_folds = 5L,
#'     model_n_reps = 4L,
#'     n_cores = 2L)
#' }
#'
#' @author Ahmed El-Gabbas
#' @export

prepare_model_data <- function(
    model_dir = NULL, climate_dir = NULL, resolution = 10L, verbose = TRUE,
    easin_ids = NULL, gbif_ids = NULL, coordinates = NULL,
    climate_scenario = "all", climate_model = "all", climate_year = "all",
    var_names = NULL, pft_type = "cross-walk", pft_id = NULL,
    bias_group = "none", vif_th = 10L, pres_buffer = 1000L, abs_buffer = 10L,
    abs_exclude_ext = list(), min_pres_grids = 50L,
    cv_folds = 5L, cv_block_size = 500L,
    cv_random = TRUE, abs_ratio = 10L, model_n_reps = 5L, n_cores = cv_folds) {

  down_clim <- down_lu <- species <- block_group <- cv <- n_pseudo_abs <-
    training_abs <- training_pres <- test_abs <- test_pres <- valid <- NULL

  ecokit::check_packages(
    c("usdm", "sdm", "parallelly", "future.apply", "future", "gtools", "fs",
      "tidyr", "withr", "sf", "terra", "crayon", "purrr", "tibble", "dplyr"))

  # # ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||| #

  # # ********************************************************************** #
  # Check inputs ------
  # # ********************************************************************** #

  # Directories
  model_dir <- ecokit::assign_from_options(
    model_dir, "onesdm_model_dir", "character")
  climate_dir <- ecokit::assign_from_options(
    climate_dir, "onesdm_climate_dir", "character")

  resolution <- ecokit::assign_from_options(
    resolution, "onesdm_resolution", c("numeric", "integer"))
  verbose <- ecokit::assign_from_options(verbose, "onesdm_verbose", "logical")

  # Species data
  easin_ids <- ecokit::assign_from_options(
    easin_ids, "onesdm_easin_ids", "character", allow_null = TRUE)
  gbif_ids <- ecokit::assign_from_options(
    gbif_ids, "onesdm_gbif_ids",
    c("numeric", "character", "integer"), allow_null = TRUE)
  coordinates <- ecokit::assign_from_options(
    coordinates, "onesdm_coordinates",
    c("data.frame", "matrix"), allow_null = TRUE)

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

  # Climate and land use data
  climate_scenario <- ecokit::assign_from_options(
    climate_scenario, "onesdm_climate_scenario", "character")
  climate_model <- ecokit::assign_from_options(
    climate_model, "onesdm_climate_model", "character")
  climate_year <- ecokit::assign_from_options(
    climate_year, "onesdm_climate_year", "character")

  var_names <- ecokit::assign_from_options(
    var_names, "onesdm_var_names", "character")
  pft_id <- ecokit::assign_from_options(
    pft_id, "onesdm_pft_id",
    c("numeric", "character", "integer"), allow_null = TRUE)
  pft_type <- ecokit::assign_from_options(
    pft_type, "onesdm_pft_type", "character")
  vif_th <- ecokit::assign_from_options(
    vif_th, "onesdm_vif_th", c("numeric", "integer"))
  bias_group <- ecokit::assign_from_options(
    bias_group, "onesdm_bias_group", "character")

  # Modelling parameters
  pres_buffer <- ecokit::assign_from_options(
    pres_buffer, "onesdm_pres_buffer", c("numeric", "integer"),
    allow_null = TRUE)
  abs_ratio <- ecokit::assign_from_options(
    abs_ratio, "onesdm_abs_ratio", c("numeric", "integer"))
  abs_buffer <- ecokit::assign_from_options(
    abs_buffer, "onesdm_abs_buffer", c("numeric", "integer"), allow_null = TRUE)
  abs_exclude_ext <- ecokit::assign_from_options(
    abs_exclude_ext, "onesdm_abs_exclude_ext", "list")
  min_pres_grids <- ecokit::assign_from_options(
    min_pres_grids, "onesdm_min_pres_grids", c("numeric", "integer"))

  model_n_reps <- ecokit::assign_from_options(
    model_n_reps, "onesdm_model_n_reps", c("numeric", "integer"))
  cv_folds <- ecokit::assign_from_options(
    cv_folds, "onesdm_cv_folds", c("numeric", "integer"))
  cv_block_size <- ecokit::assign_from_options(
    cv_block_size, "onesdm_cv_block_size", c("numeric", "integer"))
  cv_random <- ecokit::assign_from_options(
    cv_random, "onesdm_cv_random", "logical")
  n_cores <- ecokit::assign_from_options(
    n_cores, "onesdm_model_n_cores", c("numeric", "integer"))

  # # ********************************************************************** #

  dir_sub <- fs::path(model_dir, paste0("models_res_", resolution))
  dir_model_data <- fs::path(dir_sub, "modelling_data")
  fs::dir_create(c(dir_sub, dir_model_data))

  # # ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||| #

  # # ********************************************************************** #
  # Validate values ------
  # # ********************************************************************** #

  ## `climate_scenario` -----

  valid_clim_scenarios <- c("ssp126", "ssp370", "ssp585", "all", "none")
  if (!all(climate_scenario %in% valid_clim_scenarios)) {
    ecokit::stop_ctx(
      paste0(
        "Invalid climate_scenario value(s). Valid options are: ",
        toString(valid_clim_scenarios), "."),
      climate_scenario = climate_scenario, cat_timestamp = FALSE)
  }

  # # ********************************************************************** #

  ## `climate_model` -----

  valid_clim_models <- c("gfdl", "ipsl", "mpi", "mri", "ukesm1", "all", "none")
  if (!all(climate_model %in% valid_clim_models)) {
    ecokit::stop_ctx(
      paste0(
        "Invalid climate_model value(s). Valid options are: ",
        toString(valid_clim_models), "."),
      climate_model = climate_model, cat_timestamp = FALSE)
  }

  # # ********************************************************************** #

  ## `climate_year` -----

  valid_clim_years <- c("2011_2040", "2041_2070", "2071_2100", "all", "none")
  if (!all(climate_year %in% valid_clim_years)) {
    ecokit::stop_ctx(
      paste0(
        "Invalid climate_year value(s). Valid options are: ",
        toString(valid_clim_years), "."),
      climate_year = climate_year, cat_timestamp = FALSE)
  }

  # # ********************************************************************** #

  make_future_predictions <- purrr::none(
    list(climate_scenario, climate_model, climate_year),
    ~ any(.x == "none") || is.null(.x))

  if (make_future_predictions) {

    if (climate_scenario == "all") {
      climate_scenario <- c("ssp126", "ssp370", "ssp585")
    }
    if (climate_model == "all") {
      climate_model <- c("gfdl", "ipsl", "mpi", "mri", "ukesm1")
    }
    if (climate_year == "all") {
      climate_year <- c("2011_2040", "2041_2070", "2071_2100")
    }

  } else {

    ecokit::cat_time(
      paste0(
        "No climate data will be used for future predictions as one of the ",
        "climate_scenario, climate_model, or climate_year arguments is ",
        "set to `none` or `NULL`."),
      cat_timestamp = FALSE, verbose = verbose)

  }

  # # ********************************************************************** #

  ## `bias_group` ----

  valid_bias_groups <- c(
    "amphibians", "birds", "mammals", "molluscs", "plants", "reptiles", "none")

  if (is.null(bias_group) || length(bias_group) != 1L || !nzchar(bias_group)) {
    ecokit::stop_ctx(
      "The `bias_group` argument must be a single non-empty string.",
      bias_group = bias_group, cat_timestamp = FALSE)
  }

  if (!(bias_group %in% valid_bias_groups)) {
    ecokit::stop_ctx(
      paste0(
        "Invalid `bias_group` value. Valid options are: ",
        toString(valid_bias_groups), "."),
      bias_group = bias_group, cat_timestamp = FALSE)
  }

  # # ********************************************************************** #

  ## `abs_exclude_ext` ------

  # abs_exclude_ext must be a list of SpatExtent objects
  if (!is.list(abs_exclude_ext)) {
    ecokit::stop_ctx(
      "The `abs_exclude_ext` argument must be a list.",
      abs_exclude_ext = abs_exclude_ext,
      class_abs_exclude_ext = class(abs_exclude_ext), cat_timestamp = FALSE)
  }

  if (length(abs_exclude_ext) > 0L) {
    valid_extents <- purrr::map_lgl(abs_exclude_ext, inherits, "SpatExtent")
    if (!all(valid_extents)) {
      ecokit::stop_ctx(
        paste0(
          "All elements in the `abs_exclude_ext` argument must ",
          "be SpatExtent objects."),
        abs_exclude_ext = abs_exclude_ext,
        class_abs_exclude_ext = purrr::map(abs_exclude_ext, class),
        cat_timestamp = FALSE)
    }

    # check that all extents have valid coordinate ranges
    extents_vector <- purrr::map_dfr(abs_exclude_ext, as.vector)
    invalid_extents_n <- extents_vector %>%
      dplyr::mutate(
        ext_id = seq_along(abs_exclude_ext),
        valid = dplyr::case_when(
          xmin >= -180L & xmax <= 180L & ymin >= -90L & ymax <= 90L ~ TRUE,
          .default = FALSE)) %>%
      dplyr::filter(!valid)

    if (nrow(invalid_extents_n) > 0L) {
      ecokit::stop_ctx(
        paste0(
          "All extents in the `abs_exclude_ext` argument must have valid ",
          "coordinate ranges: xmin >= -180, xmax <= 180, ",
          "ymin >= -90, ymax <= 90."),
        invalid_extents_n = invalid_extents_n,
        invalid_extents_ids = invalid_extents_n$ext_id,
        cat_timestamp = FALSE)
    }
  }

  # # ********************************************************************** #

  ## `vif_th` -----

  if (!is.numeric(vif_th) || length(vif_th) != 1L || vif_th <= 0L) {
    ecokit::stop_ctx(
      "The `vif_th` argument must be a single positive numeric or `NULL`.",
      vif_th = vif_th, cat_timestamp = FALSE)
  }

  # # ********************************************************************** #

  ## `var_names` ----

  all_vars_valid <- all(var_names %in% unique(OneSDM::climate_data$var_name))
  if (!all_vars_valid) {
    invalid_var_names <- setdiff(
      var_names, unique(OneSDM::climate_data$var_name))
    ecokit::stop_ctx(
      paste0(
        "One or more variable names in `var_names` are invalid. Please ",
        "refer to `OneSDM::climate_data` for valid variable names."),
      var_names = var_names, invalid_var_names = invalid_var_names,
      cat_timestamp = FALSE)
  }

  # # ********************************************************************** #

  ## `pft_id` / `pft_type` ----

  if (!is.null(pft_id)) {

    valid_pft_types <- c("original", "cross-walk")
    if (length(pft_type) != 1L || !(pft_type %in% valid_pft_types)) {
      ecokit::stop_ctx(
        paste0(
          "The `pft_type` argument must be a single string with one of the ",
          "following values: 'original' or 'cross-walk'."),
        pft_type = pft_type, cat_timestamp = FALSE)
    }

    if (pft_type == "original") {
      valid_pft_ids <- seq_len(20L)
    } else if (pft_type == "cross-walk") {
      valid_pft_ids <- seq_len(12L)
    }

    if (!all(pft_id %in% valid_pft_ids)) {
      invalid_pft_ids <- setdiff(pft_id, valid_pft_ids)
      ecokit::stop_ctx(
        paste0(
          "One or more `pft_id` values are invalid for the specified ",
          "`pft_type`."),
        pft_type = pft_type, pft_id = pft_id,
        invalid_pft_ids = invalid_pft_ids, cat_timestamp = FALSE)
    }
  }

  # # ********************************************************************** #

  ## `cv_folds` ----

  if (!is.numeric(cv_folds) || length(cv_folds) != 1L || cv_folds <= 0L) {
    ecokit::stop_ctx(
      "The `cv_folds` argument must be a single positive numeric.",
      cv_folds = cv_folds, cat_timestamp = FALSE)
  }

  if (!is.integer(cv_folds)) {
    if (cv_folds != as.integer(cv_folds)) {
      ecokit::stop_ctx(
        "The `cv_folds` argument must be an integer value.",
        cv_folds = cv_folds, cat_timestamp = FALSE)
    }
    cv_folds <- as.integer(cv_folds)
  }

  if (cv_folds > 10L || cv_folds < 1L) {
    ecokit::stop_ctx(
      "The `cv_folds` argument must be between 1 and 10.",
      cv_folds = cv_folds, cat_timestamp = FALSE)
  }

  # # ********************************************************************** #

  ## `n_cores` ----

  if (!is.numeric(n_cores) || length(n_cores) != 1L || n_cores <= 0L) {
    ecokit::stop_ctx(
      "The `n_cores` argument must be a single positive numeric.",
      n_cores = n_cores, cat_timestamp = FALSE)
  }

  if (!is.integer(n_cores)) {
    if (n_cores != as.integer(n_cores)) {
      ecokit::stop_ctx(
        "The `n_cores` argument must be an integer value.",
        n_cores = n_cores, cat_timestamp = FALSE)
    }
    n_cores <- as.integer(n_cores)
  }

  if (n_cores > parallelly::availableCores()) {
    warning(
      "Argument `n_cores` exceeds the number of available CPU cores.\n",
      "n_cores: ", n_cores, "; available cores: ", parallelly::availableCores(),
      call. = FALSE)
    n_cores <- parallelly::availableCores()
  }

  if (n_cores > cv_folds) {
    warning(
      "Argument `n_cores` exceeds the number of CV folds (`cv_folds`).\n",
      "`n_cores`: ", n_cores, "; cv_folds: ", cv_folds,
      ". Setting `n_cores` to ", cv_folds,
      call. = FALSE)
    n_cores <- cv_folds
  }

  # # ********************************************************************** #

  ## `cv_block_size` ------

  if (!is.numeric(cv_block_size) || length(cv_block_size) != 1L ||
      cv_block_size <= 0L) {
    ecokit::stop_ctx(
      "The `cv_block_size` argument must be a single positive numeric.",
      cv_block_size = cv_block_size, cat_timestamp = FALSE)
  }

  # # ********************************************************************** #

  ## `abs_ratio` -----

  if (!is.numeric(abs_ratio) || length(abs_ratio) != 1L || abs_ratio <= 0L) {
    ecokit::stop_ctx(
      "The `abs_ratio` argument must be a single positive numeric.",
      abs_ratio = abs_ratio, cat_timestamp = FALSE)
  }

  if (!is.integer(abs_ratio)) {
    if (abs_ratio != as.integer(abs_ratio)) {
      ecokit::stop_ctx(
        "The `abs_ratio` argument must be an integer value.",
        abs_ratio = abs_ratio, cat_timestamp = FALSE)
    }
    abs_ratio <- as.integer(abs_ratio)
  }

  # # ********************************************************************** #

  ## `model_n_reps` -----

  if (!is.numeric(model_n_reps) || length(model_n_reps) != 1L ||
      model_n_reps <= 0L) {
    ecokit::stop_ctx(
      "The `model_n_reps` argument must be a single positive numeric.",
      model_n_reps = model_n_reps, cat_timestamp = FALSE)
  }

  if (!is.integer(model_n_reps)) {
    if (model_n_reps != as.integer(model_n_reps)) {
      ecokit::stop_ctx(
        "The `model_n_reps` argument must be an integer value.",
        model_n_reps = model_n_reps, cat_timestamp = FALSE)
    }
    model_n_reps <- as.integer(model_n_reps)
  }

  # # ********************************************************************** #

  ## `min_pres_grids` -----

  if (!is.numeric(min_pres_grids) || length(min_pres_grids) != 1L ||
      min_pres_grids <= 10L) {
    ecokit::stop_ctx(
      "The `min_pres_grids` argument must be a single numeric value > 10.",
      min_pres_grids = min_pres_grids, cat_timestamp = FALSE)
  }

  if (!is.integer(min_pres_grids)) {
    if (min_pres_grids != as.integer(min_pres_grids)) {
      ecokit::stop_ctx(
        "The `min_pres_grids` argument must be an integer value.",
        min_pres_grids = min_pres_grids, cat_timestamp = FALSE)
    }
    min_pres_grids <- as.integer(min_pres_grids)
  }


  # # ********************************************************************** #

  ## parameters set via options -----

  vif_sample <- getOption("onesdm_vif_sample")

  if (is.null(vif_sample)) {
    vif_sample <- 5000L
  } else {
    if (!is.numeric(vif_sample) || length(vif_sample) != 1L ||
        vif_sample <= 0L) {
      ecokit::stop_ctx(
        paste0(
          "The `onesdm_vif_sample` option must be a single positive numeric ",
          "value."),
        vif_sample = vif_sample, cat_timestamp = FALSE)
    }
    vif_sample <- as.integer(vif_sample)
  }

  # # ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||| #

  # # ********************************************************************** #
  # Prepare species data ------
  # # ********************************************************************** #

  ecokit::info_chunk(
    "Preparing species data",
    line_char_rep = 65L, verbose = verbose, cat_date = FALSE)

  exclude_extents <- ecokit::get_option_with_default(
    "onesdm_exclude_extents", "OneSDM::prepare_species_data", "exclude_extents")
  species_name <- ecokit::get_option_with_default(
    "onesdm_species_name", "OneSDM::prepare_species_data", "species_name")
  onesdm_outlier_dist_km <- ecokit::get_option_with_default(
    "onesdm_outlier_dist_km", "OneSDM::prepare_species_data", "outlier_dist_km")
  outlier_resolution <- ecokit::get_option_with_default(
    "onesdm_outlier_resolution", "OneSDM::prepare_species_data",
    "outlier_resolution")
  outlier_n_cores <- ecokit::get_option_with_default(
    "onesdm_outlier_n_cores", "OneSDM::prepare_species_data", "outlier_n_cores")
  plot_distribution <- ecokit::get_option_with_default(
    "onesdm_plot_distribution", "OneSDM::prepare_species_data",
    "plot_distribution")

  species_pa_file <- OneSDM::prepare_species_data(
    easin_ids = easin_ids, gbif_ids = gbif_ids, coordinates = coordinates,
    model_dir = model_dir, climate_dir = climate_dir, resolution = resolution,
    verbose = verbose, exclude_extents = exclude_extents,
    species_name = species_name, outlier_dist_km = onesdm_outlier_dist_km,
    outlier_resolution = outlier_resolution, outlier_n_cores = outlier_n_cores,
    plot_distribution = plot_distribution)

  if (length(species_pa_file) != 1L || !nzchar(species_pa_file)) {
    ecokit::stop_ctx(
      "The output of `prepare_species_data` is invalid.",
      species_pa_file = species_pa_file, cat_timestamp = FALSE)
  }

  if (!ecokit::check_tiff(species_pa_file, warning = FALSE)) {
    ecokit::stop_ctx(
      paste0(
        "The species presence-absence file does not exist or is not a valid",
        "GeoTIFF file."),
      species_pa_file = species_pa_file, cat_timestamp = FALSE)
  }

  species_pa_r <- terra::rast(species_pa_file)

  # # ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||| #

  # # ********************************************************************** #
  # Check number of presence grid cells -----
  # # ********************************************************************** #

  n_pres_grids <- (terra::classify(species_pa_r, cbind(0L, NA)) == 1L) %>%
    terra::global(fun = "sum", na.rm = TRUE) %>%
    as.numeric() %>%
    as.integer()

  if (n_pres_grids < min_pres_grids) {
    ecokit::stop_ctx(
      paste0(
        "The number of presence grid cells (",
        ecokit::format_number(n_pres_grids), ") ",
        "is less than the minimum required (",
        ecokit::format_number(min_pres_grids), "). ",
        "Consider adjusting the exclusion extents (if used), or ",
        "lowering the `min_pres_grids` parameter."),
      n_pres_grids = n_pres_grids, min_pres_grids = min_pres_grids,
      cat_timestamp = FALSE)
  }

  ecokit::cat_time(
    paste0(
      "Number of presence grid cells: ",
      ecokit::format_number(n_pres_grids)),
    cat_timestamp = FALSE, verbose = verbose)

  # # ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||| #

  # # ********************************************************************** #
  # Defining study area and exclude pseudo-absences ------
  # # ********************************************************************** #

  # Exclude pseudo-absences based on provided extents

  if (length(abs_exclude_ext) > 0L) {

    ecokit::cat_time(
      "\nExcluding pseudo-absences based on provided extents",
      cat_timestamp = FALSE, verbose = verbose)

    for (ext in abs_exclude_ext) {
      ecokit::cat_time(
        paste0("Extent: ", crayon::blue(as.character(ext))),
        cat_timestamp = FALSE, level = 1L, verbose = verbose)
      species_pa_r <- terra::mask(species_pa_r, ext, updatevalue = 100L) %>%
        terra::classify(cbind(0L, NA)) %>%
        terra::mask(x = species_pa_r, mask = ., updatevalue = NA) %>%
        stats::setNames(names(species_pa_r))
    }
    rm(ext, envir = environment())

  } else {
    ecokit::cat_time(
      paste0(
        "No pseudo-absences were excluded based on provided extents as ",
        "`abs_exclude_ext` is NULL or invalid."),
      cat_timestamp = FALSE, verbose = verbose)
  }

  invisible(gc())

  # # |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||| #

  # Exclude grid cells that are too far from presence locations ( >
  # pres_buffer km)

  if (!is.null(pres_buffer) && is.numeric(pres_buffer) &&
      length(pres_buffer) == 1L && pres_buffer > 0L) {

    ecokit::cat_time(
      paste0(
        "Excluding grid cells beyond ", ecokit::format_number(pres_buffer),
        " km from presence locations"),
      cat_timestamp = FALSE, verbose = verbose)

    species_pa_r <- terra::classify(species_pa_r, cbind(0L, NA)) %>%
      # Aggregate to coarser resolution to speed up distance calculation
      terra::aggregate(fact = 5L, fun = "max", na.rm = TRUE) %>%
      # Calculate distance to nearest presence cell
      terra::buffer(width = pres_buffer * 1000L) %>%
      # Resample back to original resolution
      terra::disagg(fact = 5L, method = "near") %>%
      # Ensure that buffered areas align with original raster
      terra::resample(species_pa_r, method = "near") %>%
      # Mask original raster to keep only areas within buffer
      terra::classify(cbind(0L, NA)) %>%
      terra::mask(x = species_pa_r, mask = ., updatevalue = NA) %>%
      stats::setNames(names(species_pa_r))

    invisible(gc())
  } else {
    ecokit::cat_time(
      paste0(
        "No grid cells were excluded based on distance from presence ",
        "locations as `pres_buffer` is NULL or non-positive."),
      cat_timestamp = FALSE, verbose = verbose)
  }

  # # |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||| #

  # Exclude grid cells that are too close to presence locations from
  # pseudo-absences (< `abs_buffer` km). This is needed for sampling
  # pseudo-absences only beyond a certain distance from presences, which is
  # relevant for model types other than MaxEnt.

  if (!is.null(abs_buffer) && is.numeric(abs_buffer) &&
      length(abs_buffer) == 1L && abs_buffer > 0L) {

    ecokit::cat_time(
      paste0(
        "Excluding pseudo-absences within ",
        ecokit::format_number(abs_buffer),
        " km of presence locations"),
      cat_timestamp = FALSE, verbose = verbose)

    buffered_pres <- terra::classify(species_pa_r, cbind(0L, NA)) %>%
      terra::as.polygons(aggregate = TRUE) %>%
      terra::buffer(width = abs_buffer * 1000L) %>%
      terra::rasterize(species_pa_r, field = 1L, background = NA) %>%
      stats::setNames(names(species_pa_r))

    species_pa_r <- terra::ifel(
      species_pa_r == 1L, 1L,
      terra::ifel(species_pa_r == 0L & is.na(buffered_pres), 0L, NA))

    rm(buffered_pres, envir = environment())
    invisible(gc())

  } else {
    ecokit::cat_time(
      paste0(
        "No pseudo-absences were excluded based on distance from presence ",
        "locations as `abs_buffer` is NULL or non-positive."),
      cat_timestamp = FALSE, verbose = verbose)
  }

  # # ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||| #

  # # ********************************************************************** #
  # Save species presence-absence raster ------
  # # ********************************************************************** #

  species_pa_r_file <- fs::path(
    dir_sub, paste0("species_pa_res_", resolution, ".tif"))
  terra::writeRaster(
    species_pa_r, filename = species_pa_r_file,
    overwrite = TRUE, gdal = c("COMPRESS=ZSTD", "ZSTD_LEVEL=22", "TILED=YES"))

  # # ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||| #

  # # ********************************************************************** #
  # Predictors - current ------
  # # ********************************************************************** #

  climate_preds <- OneSDM::get_climate_data(
    climate_dir = climate_dir, resolution = resolution,
    climate_scenario = "current", climate_model = "current", year = "1981_2010",
    var_names = var_names, verbose = verbose)

  if (!is.null(pft_id)) {

    lu_preds <- OneSDM::get_landuse_data(
      climate_dir = climate_dir, resolution = resolution,
      climate_scenario = "current", year = "1981_2010", pft_type = pft_type,
      pft_id = pft_id, verbose = verbose)

  } else {

    ecokit::cat_time(
      "No land use predictors will be used as `pft_id` is NULL.",
      cat_timestamp = FALSE, verbose = verbose)
    lu_preds <- tibble::tibble()

  }

  ecokit::cat_sep(
    sep_lines_before = 1L, sep_lines_after = 2L, verbose = verbose,
    line_char_rep = 65L)

  ecokit::cat_time(
    "Combining and masking predictors to study area",
    cat_timestamp = FALSE, verbose = verbose)

  model_predictors <- terra::rast(
    c(climate_preds$out_file, lu_preds$out_file)) %>%
    terra::mask(species_pa_r)

  # # ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||| #

  # # ********************************************************************** #
  # Sampling efforts ------
  # # ********************************************************************** #

  if (bias_group == "none") {

    bias_as_predictor <- FALSE
    bias_fix_value <- NA_real_
    ecokit::cat_time(
      "No sampling effort predictor will be used.",
      cat_timestamp = FALSE, verbose = verbose)

  } else {

    bias_as_predictor <- TRUE
    ecokit::cat_time(
      paste0(
        "Using sampling efforts for '", crayon::blue(bias_group),
        "' group as predictor to correct for sampling bias"),
      cat_timestamp = FALSE, verbose = verbose)

    bias_r <- OneSDM::get_sampling_efforts(
      bias_group = bias_group, resolution = resolution,
      climate_dir = climate_dir, return_spatraster = TRUE) %>%
      terra::mask(species_pa_r) %>%
      magrittr::add(1L) %>%
      log10() %>%
      stats::setNames("bias")

    # calculate 90 percentile value for later use in bias correction
    bias_fix_value <- terra::global(
      x = bias_r,
      fun = function(x, ...) {
        stats::quantile(x, probs = 0.9, na.rm = TRUE)})[[1L]]

    bias_fixed_r <- stats::setNames(bias_r, "bias_fixed")
    bias_fixed_r[bias_fixed_r <= bias_fixed_r] <- bias_fix_value

    # Add bias raster to predictors
    model_predictors <- c(model_predictors, bias_r, bias_fixed_r)

    rm(bias_fixed_r, bias_r, envir = environment())
    invisible(gc())

  }

  ecokit::cat_time(
    "Saving masked predictors", cat_timestamp = FALSE,
    verbose = verbose)
  model_predictors_f <- fs::path(
    dir_sub, paste0("model_predictors_res_", resolution, ".RData"))
  ecokit::save_as(
    object = terra::wrap(model_predictors), object_name = "model_predictors",
    out_path = model_predictors_f)

  # # ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||| #

  # # ********************************************************************** #
  # Predictor selection ------
  # # ********************************************************************** #

  vif_keep <- ecokit::get_option_with_default(
    "onesdm_vif_keep", "usdm::vifstep", "keep")

  # check that all values in `vif_keep` overlaps with predictors
  vif_keep_null <- is.null(vif_keep) || length(vif_keep) == 0L ||
    all(nzchar(vif_keep))

  if (!vif_keep_null && (!all(vif_keep %in% names(model_predictors)))) {
    ecokit::stop_ctx(
      paste0(
        "The `vif_keep` option contains variable names that are not in the ",
        "predictor set"),
      vif_keep = vif_keep, predictor_names = names(model_predictors),
      cat_timestamp = FALSE)
  }

  if (bias_group != "none") {
    vif_keep <- unique(c(vif_keep, "bias"))
  }

  ecokit::cat_time(
    paste0(
      "Performing VIF-based predictor selection (VIF threshold = ",
      ecokit::format_number(vif_th), ")"),
    cat_timestamp = FALSE, verbose = verbose)

  vif_results <- usdm::vifstep(
    model_predictors, th = vif_th, keep = vif_keep, size = vif_sample,
    method = "pearson")
  vif_file <- fs::path(dir_sub, "vif_results.RData")
  save(vif_results, file = vif_file)
  excluded_vars <- vif_results@excluded

  if (length(excluded_vars) > 0L) {

    ecokit::cat_time(
      paste0(
        "The following predictor variables were excluded based on VIF > ",
        ecokit::format_number(vif_th),
        ": ", crayon::blue(toString(excluded_vars))),
      cat_timestamp = FALSE, verbose = verbose, level = 1L)

    model_predictors <- terra::subset(
      x = model_predictors, subset = excluded_vars, negate = TRUE)

  } else {

    ecokit::cat_time(
      "No predictor variables were excluded based on VIF.",
      cat_timestamp = FALSE, verbose = verbose, level = 1L)

  }

  predictor_names <- names(model_predictors) %>%
    stringr::str_subset("bias_fixed", negate = TRUE) %>%
    gtools::mixedsort()

  ecokit::cat_time(
    paste0(
      "Final set of predictor variables: ",
      crayon::blue(toString(predictor_names))),
    cat_timestamp = FALSE, verbose = verbose, level = 1L)

  climate_predictor_names <- setdiff(climate_preds$var_name, excluded_vars)

  if (nrow(lu_preds) == 0L) {
    lu_predictor_names <- NA_character_
  } else {
    lu_predictor_names <- setdiff(
      fs::path_ext_remove(basename(lu_preds$out_file)), excluded_vars)
  }

  # # ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||| #

  # # ********************************************************************** #
  # Predictors - future ------
  # # ********************************************************************** #

  if (make_future_predictions) {

    ecokit::info_chunk(
      "Preparing future climate predictors",
      line_char_rep = 65L, verbose = verbose, cat_date = FALSE)

    climate_future <- tidyr::expand_grid(
      climate_scenario = climate_scenario,
      climate_model = climate_model,
      climate_year = climate_year)

    ecokit::cat_time(
      paste0(
        "Total climate scenario combinations to process: ",
        nrow(climate_future)),
      cat_timestamp = FALSE, verbose = verbose, level = 1L)

    ecokit::cat_time(
      "Downloading and processing future climate predictors",
      cat_timestamp = FALSE, verbose = verbose, level = 1L)

    climate_future <- climate_future %>%
      dplyr::mutate(
        down_clim = purrr::pmap(
          .l = list(climate_scenario, climate_model, climate_year),
          .f = function(climate_scenario, climate_model, climate_year) {

            ecokit::cat_time(
              paste0(
                "scenario: ", crayon::blue(climate_scenario),
                "; model: ", crayon::blue(climate_model),
                "; year: ", crayon::blue(climate_year)),
              cat_timestamp = FALSE, verbose = verbose, level = 2L)

            OneSDM::get_climate_data(
              climate_dir = climate_dir, resolution = resolution,
              climate_scenario = climate_scenario,
              climate_model = climate_model, year = climate_year,
              var_names = climate_predictor_names,
              verbose = FALSE, sleep_time = 5L)
          }
        )) %>%
      dplyr::pull(down_clim) %>%
      dplyr::bind_rows()


    if (!is.null(pft_id) && !all(is.na(lu_predictor_names))) {

      ecokit::info_chunk(
        "Preparing future landuse predictors",
        line_char_rep = 65L, verbose = verbose, cat_date = FALSE)

      lu_future <- tidyr::expand_grid(
        climate_scenario = climate_scenario,
        climate_year = climate_year)

      ecokit::cat_time(
        paste0(
          "Total land use scenario combinations to process: ", nrow(lu_future)),
        cat_timestamp = FALSE, verbose = verbose, level = 1L)

      ecokit::cat_time(
        "Downloading and processing future land use predictors",
        cat_timestamp = FALSE, verbose = verbose, level = 1L)

      lu_future <- lu_future %>%
        dplyr::mutate(
          down_lu = purrr::map2(
            .x = climate_scenario, .y = climate_year,
            .f = ~ {

              ecokit::cat_time(
                paste0(
                  "scenario: ", crayon::blue(.x), "; year: ", crayon::blue(.y)),
                cat_timestamp = FALSE, verbose = verbose, level = 2L)

              OneSDM::get_landuse_data(
                climate_dir = climate_dir, resolution = resolution,
                climate_scenario = .x, year = .y,
                pft_type = pft_type,
                pft_id = as.integer(
                  stringr::str_extract(lu_predictor_names, "[0-9]+")),
                verbose = FALSE, sleep_time = 5L)
            }
          )) %>%
        dplyr::pull(down_lu) %>%
        dplyr::bind_rows()
    }

  } else {
    climate_future <- lu_future <- tibble::tibble()

  }

  invisible(gc())

  # # ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||| #

  # # ********************************************************************** #
  # Spatial blocks ------
  # # ********************************************************************** #

  ecokit::cat_sep(
    sep_lines_before = 1L, sep_lines_after = 2L, verbose = verbose,
    line_char_rep = 65L)

  if (cv_random) {

    ecokit::cat_time(
      "Creating spatial cross-validation blocks using random blocks",
      cat_timestamp = FALSE, verbose = verbose)

    # Find aggregation factor based on cv_block_size and map resolution
    resol <- terra::res(species_pa_r)[[1L]]
    agg_factor <- ceiling(cv_block_size / (resol * 111.320))

    # aggregate species_pa_r to larger blocks, calculating sum of presences
    ecokit::cat_time(
      paste0(
        "Aggregating presence-absence raster to blocks of approx. ",
        ecokit::format_number(cv_block_size),
        " km (aggregation factor: ", ecokit::format_number(agg_factor), ")"),
      cat_timestamp = FALSE, verbose = verbose, level = 1L)

    species_blocks <- terra::aggregate(
      x = species_pa_r, fact = agg_factor, fun = "sum", na.rm = TRUE)

    # convert aggregated blocks to polygons and assign cv folds. Folds are
    # assigned randomly but maintaining balance of presences across folds
    ecokit::cat_time(
      "Assigning cross-validation folds to spatial blocks",
      cat_timestamp = FALSE, verbose = verbose, level = 1L)

    species_blocks <- species_blocks %>%
      terra::as.polygons(aggregate = FALSE) %>%
      sf::st_as_sf() %>%
      dplyr::arrange(dplyr::desc(species)) %>%
      dplyr::mutate(block_group = ceiling(dplyr::row_number() / cv_folds)) %>%
      dplyr::group_by(block_group) %>%
      dplyr::mutate(cv = sample.int(cv_folds, size = dplyr::n())) %>%
      dplyr::ungroup() %>%
      dplyr::select(-block_group) %>%
      # Convert back to raster
      terra::vect() %>%
      terra::rasterize(y = species_blocks, field = "cv", background = NA) %>%
      # convert back to original resolution
      terra::disagg(fact = agg_factor, method = "near", progress = FALSE) %>%
      # mask by original buffered raster
      terra::mask(species_pa_r, progress = FALSE) %>%
      stats::setNames("cv")

    species_block_cv_file <- NA_character_

    rm(agg_factor, resol, envir = environment())
    invisible(gc())

  } else {

    ecokit::cat_time(
      "Creating spatial cross-validation blocks using `blockCV` package",
      cat_timestamp = FALSE, verbose = verbose)

    ecokit::check_packages("blockCV")

    ecokit::cat_time(
      "Converting presence-absence raster to points",
      cat_timestamp = FALSE, verbose = verbose, level = 1L)
    species_pa_sf <- terra::as.points(species_pa_r) %>%
      sf::st_as_sf() %>%
      stats::setNames(c("species", "geometry"))

    ecokit::cat_time(
      paste0(
        "Creating spatial cross-validation blocks using `blockCV` package.",
        "\n  >>>  This could take a while, depending on the size ",
        "of the study area"),
      cat_timestamp = FALSE, verbose = verbose, level = 1L)

    species_blocks <- blockCV::cv_spatial(
      x = species_pa_sf, column = "species", r = species_pa_r, hexagon = FALSE,
      # selection = "random", iteration = 100L,
      k = cv_folds, size = cv_block_size * 1000L, plot = FALSE,
      progress = verbose, report = verbose)

    ecokit::cat_time(
      "Save output of `blockCV` spatial blocks",
      cat_timestamp = FALSE, verbose = verbose, level = 1L)
    species_block_cv_file <- fs::path(
      dir_sub, paste0("species_cv_blockCV_", resolution, ".RData"))
    ecokit::save_as(
      object = species_blocks, object_name = "species_blocks",
      out_path = species_block_cv_file)

    # Extract raster of spatial blocks
    ecokit::cat_time(
      "Extracting spatial blocks raster from `blockCV` output",
      cat_timestamp = FALSE, verbose = verbose, level = 1L)

    species_blocks <- dplyr::select(
      species_blocks$blocks, tidyselect::all_of("folds")) %>%
      terra::vect() %>%
      terra::rasterize(y = species_pa_r, field = "folds", background = NA) %>%
      terra::mask(species_pa_r, progress = FALSE) %>%
      stats::setNames("cv")
    terra::varnames(species_blocks) <- ""

    rm(species_pa_sf, envir = environment())
    invisible(gc())

  }

  ecokit::cat_time(
    "Save spatial blocks to GeoTIFF file",
    cat_timestamp = FALSE, verbose = verbose, level = 1L)
  species_blocks_file <- fs::path(
    dir_sub, paste0("species_cv_blocks_", resolution, ".tif"))
  terra::writeRaster(
    species_blocks, filename = species_blocks_file,
    overwrite = TRUE, gdal = c("COMPRESS=ZSTD", "ZSTD_LEVEL=22", "TILED=YES"))

  # # ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||| #

  # # ********************************************************************** #
  # Prepare modelling data ------
  # # ********************************************************************** #

  ecokit::info_chunk(
    "Preparing modelling data",
    line_char_rep = 65L, verbose = verbose, cat_date = FALSE)

  model_data_r <- c(species_pa_r, model_predictors, species_blocks)
  model_data <- as.data.frame(model_data_r, xy = TRUE, na.rm = TRUE) %>%
    tibble::tibble()
  model_data_f <- fs::path(dir_model_data, "model_data_res_.RData")
  ecokit::save_as(
    object = model_data, object_name = "model_data", out_path = model_data_f)

  # save model_data_r raster as wrapped terra SpatRaster, and load as needed to
  # free up memory
  model_data_r_file <- fs::path(dir_model_data, "model_data_r.RData")
  ecokit::save_as(
    object = terra::wrap(model_data_r), object_name = "model_data_r",
    out_path = model_data_r_file)
  rm(species_pa_r, model_data_r, species_blocks, envir = environment())
  invisible(gc())

  if (model_n_reps >= 1L) {
    pseudo_abs_check <- model_data  %>%
      dplyr::select(tidyselect::all_of(c("species", "cv"))) %>%
      dplyr::group_by(cv) %>%
      dplyr::summarise(
        test_pres = sum(species == 1L), test_abs = sum(species == 0L)) %>%
      dplyr::mutate(
        training_pres = sum(test_pres) - test_pres,
        training_abs = sum(test_abs) - test_abs,
        n_pseudo_abs = pmin(training_pres * abs_ratio, training_abs),
        n_reps = ifelse(
          n_pseudo_abs >= (0.80 * training_abs), 1L, model_n_reps))

    if (all(pseudo_abs_check$n_reps == model_n_reps)) {
      ecokit::cat_time(
        paste0(
          "Sampling pseudo-absences with ", ecokit::format_number(model_n_reps),
          " model repetitions"),
        cat_timestamp = FALSE, verbose = verbose, level = 1L)
    } else {
      ecokit::cat_time(
        paste0(
          "Adjusting number of model fitting repetitions to ",
          ecokit::format_number(1L), ", as the number of\n  >>>  potential ",
          "pseudo-absences for some cross-validation folds is limited,",
          "\n  >>>  exceeding 80% of available absences."),
        cat_timestamp = FALSE, verbose = verbose, level = 1L)
      model_n_reps <- 1L
    }
  }

  if (n_cores == 1L) {
    future::plan("sequential", gc = TRUE)
  } else {
    ecokit::set_parallel(
      min(n_cores, cv_folds), show_log = FALSE, future_max_size = 2000L)
    withr::defer(future::plan("sequential", gc = TRUE))
  }

  par_packages <- c(
    "fs", "terra", "ecokit", "dplyr", "magrittr", "sdm", "tibble", "purrr")
  par_globals <- c(
    "model_data", "dir_model_data", "verbose", "abs_ratio",
    "bias_fix_value", "model_n_reps",  "model_data_r_file", "bias_as_predictor")

  model_data_cv <- future.apply::future_lapply(
    X = seq_len(cv_folds),
    FUN = function(cv_fold) {

      model_rep_id <- NULL

      # Training data - tibble -----
      training_data_file <- fs::path(
        dir_model_data, paste0("training_data_cv_", cv_fold, ".RData"))
      if (ecokit::check_data(training_data_file, warning = FALSE)) {
        training_data <- ecokit::load_as(training_data_file)
      } else {
        training_data <- dplyr::filter(model_data, cv != cv_fold)
        ecokit::save_as(
          object = training_data, object_name = "training_data",
          out_path = training_data_file)
      }

      # Training presences - tibble  ------
      training_pres_file <- fs::path(
        dir_model_data, paste0("training_pres_cv_", cv_fold, ".RData"))
      if (ecokit::check_data(training_pres_file, warning = FALSE)) {
        training_pres <- ecokit::load_as(training_pres_file)
      } else {
        training_pres <- dplyr::filter(training_data, species == 1L)
        ecokit::save_as(
          object = training_pres, object_name = "training_pres",
          out_path = training_pres_file)
      }
      n_training_pres <- nrow(training_pres)
      rm(training_pres, envir = environment())
      invisible(gc())

      # Training absences - tibble ------
      training_abs_file <- fs::path(
        dir_model_data, paste0("training_abs_cv_", cv_fold, ".RData"))
      if (ecokit::check_data(training_abs_file, warning = FALSE)) {
        training_abs <- ecokit::load_as(training_abs_file)
      } else {
        training_abs <- dplyr::filter(training_data, species == 0L)
        ecokit::save_as(
          object = training_abs, object_name = "training_abs",
          out_path = training_abs_file)
      }
      n_training_abs <- nrow(training_abs)

      rm(training_data, training_abs, envir = environment())
      invisible(gc())

      # Training data - raster -----
      training_r_file <- fs::path(
        dir_model_data, paste0("training_r_cv_", cv_fold, ".RData"))
      if (ecokit::check_data(training_r_file, warning = FALSE)) {
        training_r <- ecokit::load_as(training_r_file, unwrap_r = TRUE)
      } else {
        training_r <- ecokit::load_as(model_data_r_file, unwrap_r = TRUE)
        training_r[training_r$cv == cv_fold] <- NA
        ecokit::save_as(
          object = terra::wrap(training_r),
          object_name = paste0("training_r_cv_", cv_fold),
          out_path = training_r_file)
      }

      # Pseudo-absences -------
      n_pseudo_abs <- pmin(n_training_pres * abs_ratio, n_training_abs)

      training_r_sub <- terra::subset(
        x = training_r, subset = "species", negate = TRUE)
      if (bias_as_predictor) {
        training_r_sub$bias_fixed <- NULL
      }
      training_r_pres <- training_r
      training_r_pres[training_r_pres$species == 0L] <- NA
      training_r_pres <- terra::subset(
        x = training_r_pres, subset = "species") %>%
        terra::as.points()
      rm(training_r, envir = environment())
      invisible(gc())

      pseudo_abs <- tibble::tibble(model_rep_id = seq_len(model_n_reps)) %>%
        dplyr::mutate(
          pseudo_abs_file = purrr::map_chr(
            .x = model_rep_id,
            .f = function(rep_id) {

              pseudo_abs_file <- fs::path(
                dir_model_data,
                paste0("pseudo_abs_cv_", cv_fold, "_rep_", rep_id, ".RData"))

              if (!ecokit::check_data(pseudo_abs_file, warning = FALSE)) {
                pseudo_abs <- sdm::background(
                  x = training_r_sub, n = n_pseudo_abs,
                  method = "eDist", sp = training_r_pres)

                ecokit::save_as(
                  object = pseudo_abs,
                  object_name = paste0(
                    "pseudo_abs_cv_", cv_fold, "_rep_", rep_id),
                  out_path = pseudo_abs_file)
              }
              pseudo_abs_file
            }
          ))

      rm(training_r_sub, training_r_pres, envir = environment())
      invisible(gc())


      # Testing data -----
      testing_data_file <- fs::path(
        dir_model_data, paste0("testing_data_cv_", cv_fold, ".RData"))
      if (ecokit::check_data(testing_data_file, warning = FALSE)) {
        testing_data <- ecokit::load_as(testing_data_file)
      } else {
        testing_data <- dplyr::filter(model_data, cv != cv_fold)
        # if bias correction is used, add bias_fixed column
        if (bias_as_predictor) {
          testing_data <- dplyr::mutate(
            testing_data, bias_fixed = bias_fix_value)
        }
        ecokit::save_as(
          object = testing_data, object_name = "testing_data",
          out_path = testing_data_file)
      }
      n_test_pres <- sum(testing_data$species == 1L)
      n_test_abs <- sum(testing_data$species == 0L)

      rm(testing_data, envir = environment())
      invisible(gc())


      # Testing raster -----
      testing_r_file <- fs::path(
        dir_model_data, paste0("testing_r_cv_", cv_fold, ".RData"))
      if (!ecokit::check_data(testing_r_file, warning = FALSE)) {
        testing_r <- ecokit::load_as(model_data_r_file, unwrap_r = TRUE)
        testing_r[testing_r$cv == cv_fold] <- NA

        # if bias correction is used, set bias values to fixed value
        if (bias_as_predictor) {
          testing_r$bias_fixed <- terra::classify(
            testing_r$bias, rbind(cbind(-Inf, Inf, bias_fix_value)))
        }

        ecokit::save_as(
          object = terra::wrap(testing_r),
          object_name = paste0("testing_r_cv_", cv_fold),
          out_path = testing_r_file)
        rm(testing_r, envir = environment())
        invisible(gc())
      }

      tibble::tibble(
        cv = cv_fold,
        training_data = training_data_file,
        training_data_r = training_r_file,
        training_pres = training_pres_file,
        n_training_pres = n_training_pres,
        training_abs = training_abs_file,
        n_training_abs = n_training_abs,
        n_model_reps = model_n_reps,
        n_pseudo_abs = n_pseudo_abs,
        pseudo_abs = pseudo_abs,
        testing_data = testing_data_file,
        testing_data_r = testing_r_file,
        n_test_pres = n_test_pres,
        n_test_abs = n_test_abs) %>%
        tidyr::unnest(col = "pseudo_abs")

    },
    future.scheduling = Inf, future.seed = TRUE,
    future.packages = par_packages,
    future.globals = par_globals) %>%
    dplyr::bind_rows()

  future::plan("sequential", gc = TRUE)

  ## Saving model data
  ecokit::cat_time(
    "Saving model data", cat_timestamp = FALSE, verbose = verbose, level = 1L)

  ecokit::save_as(
    object = model_data_cv, object_name = "model_data_cv",
    out_path = fs::path(dir_model_data, "model_data_cv.RData"))

  # # ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||| #

  # # ********************************************************************** #
  # Summary of modelling data ------
  # # ********************************************************************** #

  if (length(abs_exclude_ext) > 0L) {
    abs_exclude_ext <- purrr::map(abs_exclude_ext, terra::wrap)
  }

  model_data_summary <- list(
    # model output directory
    model_dir = model_dir,
    # species PA data
    species_data_raw = list(
      gbif_ids = gbif_ids,
      easin_ids = easin_ids,
      coordinates = coordinates,
      # path to species presence-absence GeoTIFF file
      species_data_tiff = species_pa_file),
    # climate options and data
    climate = list(
      # directory with climate data
      climate_dir = climate_dir,
      # input climate variable names
      var_names = var_names,
      # whether to make future predictions
      future_predictions = make_future_predictions,
      # future climate options
      climate_scenario = climate_scenario,
      climate_model = climate_model,
      climate_year = climate_year,
      # tibble with current and future climate data
      climate_current = climate_preds,
      climate_future = climate_future),
    # land use options and data
    land_use = list(
      # land use PFT type and IDs
      pft_type = pft_type, pft_id = pft_id,
      # tibble with current and future land use data
      lu_current = lu_preds, lu_future = lu_future),
    # sampling bias options
    bias = list(
      bias_group = bias_group,
      bias_as_predictor = bias_as_predictor,
      bias_fix_value = bias_fix_value),
    # variable selection options and results
    var_selection = list(
      vif_th = vif_th,
      vif_keep = vif_keep,
      excluded_vars = excluded_vars,
      vif_results_file = vif_file),
    model_options = list(
      # modelling resolution
      resolution = resolution,
      # minimum required presence grid cells
      min_pres_grids = min_pres_grids,
      # pseudo-absence sampling options
      pres_buffer = pres_buffer, abs_buffer = abs_buffer,
      abs_exclude_ext = abs_exclude_ext, abs_ratio = abs_ratio,
      # cross-validation options
      cv_folds = cv_folds, cv_block_size = cv_block_size,
      cv_random = cv_random,
      # spatial blocks files (if used blockCV package)
      species_block_cv_file = species_block_cv_file,
      # GeoTIFF file for spatial blocks
      species_blocks_file = species_blocks_file,
      # number of model fitting repetitions
      model_n_reps = model_n_reps),
    model_data = list(
      # final presence-absence GeoTIFF file
      species_pa_file  = species_pa_r_file,
      # final predictor names
      predictor_names = predictor_names,
      climate_predictor_names = climate_predictor_names,
      lu_predictor_names = lu_predictor_names,
      # number of presence grid cells
      n_pres_grids = n_pres_grids,
      # model predictors raster file
      model_predictors = model_predictors_f,
      # model predictors raster file (include species and blocks)
      model_data_r_file = model_data_r_file,
      # model data tibble file
      model_data_f = model_data_f,
      # model data for cross-validation folds and repetitions
      model_data_cv = model_data_cv)
  )

  ecokit::save_as(
    object = model_data_summary, object_name = "model_data_summary",
    out_path = fs::path(dir_model_data, "model_data_summary.RData"))

  return(invisible(model_data_cv))

}
