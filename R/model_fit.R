# process_models -------

#' Prepare Modelling Data, Fit Models, and Summarise Model Outputs
#'
#' Orchestrates the end-to-end species distribution modelling (SDM) workflow.
#' This includes preparing occurrence data, downloading and preparing predictors
#' (climate and land use), handling spatial sampling bias, variable selection
#' (collinearity filtering via VIF), building spatial-blocks for
#' cross-validation (CV), assembling modelling datasets, fitting multiple SDMs
#' using different methods with repetitions across CV folds, generating current and 
#' future projections, and producing comprehensive summaries at repetition, CV-fold,
#' and overall levels. All intermediate artefacts and summaries are saved to a
#' structured output directory.
#'
#' @param model_dir,models_prefix Character.
#'   - `model_dir` (required) represents the path to the root modelling
#'   directory where data and fitted models will be saved. This should be a
#'   single directory for all models of a given species (a separate directory is
#'   expected in case of modelling multiple species, otherwise data files may be
#'   overwritten or mixed up). It can also be set via the "`onesdm_model_dir`"
#'   option.
#'   - `models_prefix` (defaults to "`models`") is the character prefix for
#'   model subdirectory under `model_dir`. This is helpful when fitting multiple
#'   sets of models for the same species (e.g., using different SDM methods or
#'   settings) to avoid overwriting files.
#'   - The function creates two subdirectories: "`data`" to store processed
#'   species data, and "`<models_prefix>_res_<resolution>`" to store model data
#'   and outputs for the specified resolution.
#' @param climate_dir Character (required). Directory path for climate data
#'   files. This directory is used to load/download mask layers matching the
#'   specified resolution and current/future climate and land use data (see
#'   below). The same directory should be used across `OneSDM` workflows for
#'   different models (including different species, predictors, resolutions,
#'   etc.) to ensure consistency and redundant re-download of large geospatial
#'   data. This can be set via the "`onesdm_climate_dir`" option.
#' @param resolution Integer. Spatial resolution used to prepare data for
#'   analysis and model fitting. Valid values are 5, 10 (default), or 20,
#'   corresponding to approximate spatial resolutions of 5, 10, and 20 km (2.5,
#'   5, and 10 arc-minutes) respectively. This value can also be set via the
#'   "`onesdm_resolution`" option.
#' @param climate_scenario,climate_model,climate_year Character vectors. Future
#'   climate scenario(s), model(s), and year(s). These three parameters control
#'   which future climate data will be prepared for future predictions (data for
#'   selected variables at current climate conditions will be always prepared):
#'   - **`climate_scenario`**: Shared Socioeconomic Pathways (SSPs). Valid
#'   values are: `"ssp126"`, `"ssp370"`, `"ssp585"`, `"all"` (default), or
#'   `"none"`. This can also be set via the "`onesdm_climate_scenario`" option.
#'   - **`climate_model`**: abbreviation for future climate model(s). Valid
#'   values are: `"gfdl"`, `"ipsl"`, `"mpi"`, `"mri"`, `"ukesm1"`, `"all"`
#'   (default), or `"none"`. This can also be set via the
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
#'   parameters for land-use predictors (optional). See [OneSDM::landuse_data]
#'   and [OneSDM::get_landuse_data] for details.
#'   - **`pft_type`**: PFT to download land-use predictors for. Must be one of
#'   "`cross-walk`" (default) or "`original`". This can be set via the
#'   "`onesdm_pft_type`" option.
#'   - **`pft_id`**: Numeric vector. One or more PFT identifiers to download.
#'   Must be valid for the specified `pft_type`: 1-20 for `"original"`, and 1-12
#'   for `"cross-walk"`. If `NULL` (default), no land-use predictors are used.
#'   This can be set via the "`onesdm_pft_id`" option.
#'
#' @param bias_group Character scalar. Taxonomic group used to compute sampling
#'   effort surface. Valid options are `"amphibians"`, `"birds"`, `"mammals"`,
#'   `"molluscs"`, `"plants"`, `"reptiles"`, or `"none"`; all in small case.
#'   Default is `"none"`; i.e., no sampling-effort predictor is used. If not
#'   "none", the function downloads sampling-effort raster of the specified
#'   group from [Zenodo](https://zenodo.org/records/7556851). See
#'   [OneSDM::get_sampling_efforts] for details. This can also be set via the
#'   "`onesdm_bias_group`" option.
#'   - When not "none", the sampling-effort raster
#'   (log<sub>10</sub> scale) is used as a predictor in the models. The
#'   90<sup>th</sup>-percentile value of the  sampling-effort raster in the
#'   modelling study area (after filtering; i.e., presence grid cells +
#'   candidate grid cells for pseudo-absence sampling, see blow) is used across
#'   the whole study area to correct for sampling bias in projections (but not
#'   during models evaluation). For more details on model-based sampling bias
#'   correction, see Warton et al.
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
#'   sampled anywhere in the study area except at presence grid cells). Default
#'   is `10` km. This can be set via the "`onesdm_abs_buffer`" option.
#' @param abs_exclude_ext List (Optional). A list of terra `SpatExtent` objects,
#'   typically generated using [terra::ext]. Each extent defines a geographic
#'   area (rectangles) where pseudo-absences will not be sampled. This is useful
#'   for excluding regions such as areas outside the species' native range or
#'   regions with specific environmental conditions that preclude the species'
#'   presence. If an empty list (default), no additional exclusion extents are
#'   applied. This can also be set via the "`onesdm_abs_exclude_ext`" option.
#' @param min_pres_grids Integer. Minimum number of presence grid cells required
#'   to proceed with modelling. If the number of presence grid cells after
#'   filtering is less than this value, the function will stop with an error.
#'   Default is **`50`** grids. This can also be set via
#'   "`onesdm_min_pres_grids`" option.
#' @param cv_folds,cv_block_size,cv_random Parameters controlling spatial-block
#'   cross-validation.
#'   - **`cv_folds`**: Integer. Number of cross-validation folds. Default is
#'   `5`. This can also be set via the "`onesdm_cv_folds`" option.
#'   - **`cv_block_size`**: Numeric. Approximate cross-validation block
#'   size in kilometres. Larger block sizes lead to fewer, more spatially
#'   distinct blocks. The function checks if the provided value is feasible
#'   given the study area extent and number of CV folds; if not, it adjusts the
#'   block size to ensure at least `cv_folds` blocks can fit within the study
#'   area. This helps to avoid situations where the block size is too large
#'   relative to the study area, leading to insufficient blocks for
#'   cross-validation. Default is `500` km. This can also be set via the
#'   "`onesdm_cv_block_size`" option.
#'   - **`cv_random`**: Logical. If `TRUE` (default), create random spatial
#'   blocks by aggregating and randomly assigning blocks into `cv_folds` groups.
#'   Although this is done randomly, the function tries to ensure that each fold
#'   has a balanced number of presence grid cells. If `FALSE`,
#'   [blockCV::cv_spatial] is used for block determination. Note that
#'   [blockCV::cv_spatial] can be computationally intensive for large study
#'   areas and small block sizes. This can also be set via the
#'   "`onesdm_cv_random`" option.
#' @param abs_ratio,model_n_reps Parameters for sampling pseudo-absences.
#'   - **`abs_ratio`**: Integer determining the ratio of sampled pseudo-absences
#'   to training presences for each cross-validation fold. Default is `20`,
#'   which means that for each training presence grid cell, 20 pseudo-absence
#'   grid cells will be sampled. This can also be set via the
#'   "`onesdm_abs_ratio`" option.
#'   - **`model_n_reps`**: Integer representing the number of repetition
#'   datasets (different sets of pseudo-absences) to generate for each
#'   cross-validation fold. Default is `5`, for five different pseudo-absence
#'   sample datasets per CV fold (i.e., 5 fitted models per CV). If the
#'   available non-presence grid cells become limiting (when the number of
#'   requested pseudo-absences exceeds 80% of available non-presence grid
#'   cells), the function reduces the number of repetitions to `1` as further
#'   repetitions would be redundant. This can also be set via the
#'   "`onesdm_model_n_reps`" option.
#'   - Pseudo-absences are sampled using [sdm::background] with method `"eDist"`
#'   , which draws pseudo-absences weighted by environmental distance from
#'   presence grid cells.
#'   - No pseudo-absences are sampled for maxent models; instead, all
#'   available non-presence grid cells are used as background points (capped to
#'   1,000,000 points for memory efficiency; can also be set using the
#'   "`onesdm_max_n_bg_maxent`") option. Therefore, for maxent models,
#'   `model_n_reps` is set to `1` regardless of the provided value.
#' @param n_cores Integer. Number of CPU cores to use for parallel processing of
#'   data for cross-validation folds and for model fitting (capped to
#'   [parallelly::availableCores]). This can also be set via the
#'   "`onesdm_model_n_cores`" option. Default is `4`.
#' @param sdm_methods Character vector. List of SDM methods to be used for model
#'   fitting.
#'   - Default is c("glm", "glmnet", "gam", "maxent", "rf").
#'   - Valid methods are those supported by the `sdm` R package: `glm`,
#'   `glmpoly`, `gam`, `glmnet`, `mars`, `gbm`, `rf`, `ranger`, `cart`, `rpart`,
#'   `maxent`, `mlp`, `svm`, `mda`, and `fda`. See documentation of the `sdm` R
#'   package for more details.
#'   - The `svm` method refers to the [kernlab::ksvm] implementation of SVM.
#'   However, the function also support the SVM implementation from the
#'   [e1071::svm] R package, referred to as `svm2` here. This can be useful if
#'   users prefer this implementation or have specific reasons to use it.
#'   - This can also be set via the "`onesdm_sdm_methods`" option.
#' @param sdm_settings List or `NULL`. A named list of lists, each containing
#'   method-specific settings for SDM methods. Each list element should be named
#'   after a valid SDM method (see "`sdm_methods`" above) and contain a list of
#'   settings specific to that method.
#'   - For example: `list(glm = list(control = list(maxit = 200L)))`.
#'   - If `NULL` (default), pre-defined default settings are used; see
#'   [sdm_model_settings].
#'   - If an empty list, no additional method-specific settings are applied, and
#'   default settings of the `sdm` R package are used.
#'   - This can also be set via the "`onesdm_sdm_settings`" option.
#' @param proj_extent Character scalar or list of `terra::SpatExtent` objects.
#'   This parameter controls the geographic extent for saving model projections
#'   under current and future climate scenarios.
#'   - Projections are made for all combinations of SDM model types, model
#'   repetitions, cross-validation folds, and climate options.
#'   - Projections are run in parallel if `n_cores` is set to > 1.
#'   - Processing time and memory consumption during projections can be high
#'   depending on the extent of the study area, modelling resolution, and
#'   modelling options used.
#'   - Caution should be taken when projecting to areas very far from the model
#'   fitting area due to high expected environmental extrapolation.
#'   - Valid options are:
#'     - `"modelling_area"` (*default*): predict only at grid cells within the
#'   modelling study area (i.e., after filtering, where presence/non-presences
#'   data are prepared).
#'     - `"europe"`: predict across Europe (using a predefined extent).
#'     - `"world"`: predict across the whole world land. This can highly slow
#'   down the projection step (high memory consumption and processing time) and
#'   should be chosen carefully only if needed.
#'     - A list of `terra::SpatExtent` objects, each defining a custom
#'   geographic extent for projections. In this case, global land grid cells
#'   overlap with these extents will be used for projections (irrespective of
#'   the modelling area).
#'
#' @inheritParams prepare_species_data
#' @inheritParams get_climate_data
#'
#' @details The function performs the following major steps:
#'
#'   - Validates input arguments and resolves defaults from package options.
#'   - Calls [OneSDM::prepare_species_data] to check/prepare species
#'   distribution data. This includes processing/loading occurrence data from
#'   EASIN/GBIF or user-provided coordinates.
#'   - Masks and filters the study area using optional exclusion extents and
#'   distances (`abs_exclude_ext`, `pres_buffer`, and `abs_buffer`).
#'   -  Loads/download current climate (and optionally land-use) data for
#'   selected predictors via [OneSDM::get_climate_data] and
#'   [OneSDM::get_landuse_data], combines them and masks them to the study area.
#'   - Loads/download future climate (and optionally land-use) predictors for
#'   each requested scenario \eqn{\times} model \eqn{\times} year combination
#'   (if any) via [OneSDM::get_climate_data] and [OneSDM::get_landuse_data], and
#'   masks them to the study area.
#'   - Optionally adds a sampling-effort (bias) predictor when `bias_group`
#'   is valid.
#'   - Performs VIF-based predictor selection (using [usdm::vifstep]) and
#'   excludes highly collinear variables. The `bias` predictor, if used, is
#'   always retained.
#'   - Creates spatial cross-validation blocks either with a simple random-block
#'   aggregation or via [blockCV::cv_spatial] (controlled by `cv_random`).
#'   - Samples pseudo-absences per cross-validation fold/repetition.
#'   - Assembles modelling data (species, predictors, CV fold) and persists both
#'   data.frame and wrapped `SpatRaster` for memory efficiency.
#'   - Determines feasible repetitions per CV fold based on available absences
#'   and requested `abs_ratio`.
#'   - Enforces minimum presence count after filtering; otherwise, stops with
#'   an error; controlled by `min_pres_count`.
#'   - Saves a tibble describing all model input data, including
#'   cross-validation folds and model repetitions with file paths to persisted
#'   datasets (training/testing rasters and pseudo-absence objects).
#'   - Fits models in parallel across tasks, capturing: fitted model path,
#'   training and testing evaluation metrics, variable importance, response
#'   curve data, and projection outputs for current and future scenarios.
#'   - Summarises results:
#'     - Per-CV-fold: aggregates repetitions (mean and sd where applicable) for
#'   evaluations, variable importance, response curves; compiles projections.
#'     - Overall (per method): aggregates across folds for the same metrics and
#'   summarizes projections.
#'   - Plots maps for species distribution data to be used in the models,
#'   sampling efforts data used, and cross-validation spatial blocks used.
#'   - Saves all relevant intermediate files and final summaries to structured
#'   directories: `modelling_data` (model inputs), `fitted_models` (fitted model
#'   objects and results), "`projections_reps`" (raw projection maps per model
#'   repetition), "`projections_cv`" (summary projection maps per CV fold), and
#'   "`projections`" (overall summary maps), and metadata `RData` files.
#'
#' @return A a summary list for the parameters used and file paths to important
#'   intermediate data objects under
#'   "`<model_dir>/<models_prefix>_res_<resolution>/model_summary.RData`". This
#'   list contains the following elements:
#'   - **`dir_models_root`**: Character. The root modelling directory path
#'   defined via the "`<model_dir>`" argument.
#'   - **`dir_models`**: Character. The modelling directory path:
#'   "`<model_dir>/<models_prefix>_res_<resolution>`".
#'   - **`species_data_raw`**: List of GBIF and EASIN IDs and/or user-provided
#'   coordinates used. This also includes the file path to the saved GeoTIFF
#'   file for prepared species distribution (`species_data_tiff`).
#'   - **`climate`**: List. Information on the current and future climate data
#'   used, including climate directory path (`<climate_dir>`), climate variable
#'   names (`<var_names>`), whether to make future projections
#'   (`future_predictions`), future climate options (`<climate_scenario>`,
#'   `<climate_model>`, and `<climate_year>`), and tibbles for processed current
#'   and future climate data (`climate_current` and `climate_future`).
#'   - **`landuse`**: List. Information on the land-use data used, including PFT
#'   type (`<pft_type>`), PFT IDs (`<pft_id>`), and tibbles for processed
#'   current and future land use data (`lu_current` and `lu_future`).
#'   - **`bias`**: List. Information on the sampling-effort (bias) predictor
#'   used, if any, including the taxonomic group `<bias_group>`, whether
#'   sampling-effort predictor is used (`bias_as_predictor`), and the value of
#'   the 90<sup>th</sup>-percentile used for bias correction (`bias_fix_value`).
#'   - **`var_selection`**: List. Information on the variable selection process,
#'   including the VIF threshold used (`<vif_th>`), names of always kept
#'   (`vif_keep`) or excluded (`excluded_vars`) predictors and file path for the
#'   VIF selection summary (`vif_results_file`).
#'   - **`model_options`**: List. Information on modelling options used,
#'   including resolution (`<resolution>`), minimum required presence grid cells
#'   (`<min_pres_grids>`), pseudo-absence sampling options (`<pres_buffer>` ,
#'   `<abs_buffer>`, `<abs_exclude_ext>`, and `<abs_ratio>`), cross-validation
#'   options (`<cv_folds>`, `<cv_block_size>`, and `<cv_random>`), spatial
#'   blocks files if used blockCV package (`species_block_cv_file`), GeoTIFF
#'   file for spatial blocks (`species_blocks_file`), number of model fitting
#'   repetitions (`<model_n_reps>`), and SDM methods and settings used
#'   (`<sdm_methods>` and `<sdm_settings>`).
#'   - **`model_data`**: List. Information on the prepared modelling data and
#'   raw outputs, including final presence-absence GeoTIFF file
#'   (`species_pa_file`), final filtered predictor names (`predictor_names`,
#'   `climate_predictor_names`, and `lu_predictor_names`), number of presence
#'   grid cells (`n_pres_grids`), predictors as `PackedSpatRaster` file
#'   (`model_predictors`, and `model_data_r_file`; the latter also includes
#'   species distribution and CV blocks), modelling data as tibble
#'   (`model_data_file`, including paths to cv \eqn{\times} repetition training
#'   and testing datasets and sampled pseudo-absences), paths to climate and
#'   land use data needed for making spatial projections (`projection_inputs`),
#'   projection extent (`<proj_extent>` and `proj_mask_file` for projection
#'   extent as `PackedSpatRaster`), and tibble for fitted models, raw model
#'   repetition outputs, and spatial projection paths (`models_cv`).
#'   - **`models_cv_summary`/`models_cv_summary_file`**: A tibble summarising
#'   model outputs (evaluation metrics, variable importance, response curves,
#'   and projections) per cross-validation fold and method (i.e., summarised
#'   over model repetitions). Saved as an `RData` file under
#'   "`<model_dir>/<models_prefix>_res_<resolution>/models_cv_summary.RData`".
#'  - **`models_summary`/`models_summary_file`**: A tibble summarising the
#'   overall model outputs (evaluation metrics, variable importance, response
#'   curves, and projections) per SDM method. Saved as an `RData` file under
#'   "`<model_dir>/<models_prefix>_res_<resolution>/models_summary.RData`".
#'
#' @examples
#' \dontrun{
#'
#'   OneSDM::process_models(
#'     model_dir = "test/Acacia_karroo",
#'     models_prefix = "model_test2",
#'     climate_dir = "test/climate_data",
#'     # resolution = 10L,
#'     # verbose = TRUE,
#'     easin_ids = "R00042",
#'     gbif_ids = 2979128,
#'     # coordinates = NULL,
#'     climate_scenario = c("ssp370", "ssp585"),
#'     climate_model = c("mri", "ukesm1"),
#'     climate_year = c("2041_2070", "2071_2100"),
#'     var_names = c("bio1", "bio2", "bio12", "bio15", "bio18"),
#'     # pft_type = "cross-walk",
#'     pft_id = c(5L, 8L, 10L),
#'     bias_group = "plants",
#'     # vif_th = 10L
#'     # pres_buffer = 1000L,
#'     # abs_buffer = 10L,
#'     # abs_exclude_ext = list(),
#'     # min_pres_grids = 50L,
#'     # cv_folds = 5L,
#'     # cv_block_size = 500L,
#'     # cv_random = TRUE,
#'     # abs_ratio = 20L,
#'     # model_n_reps = 5L,
#'     # n_cores = 4L,
#'     sdm_methods = c("glm", "glmnet", "gam", "maxent"),
#'     # sdm_settings = NULL,
#'     proj_extent = "europe"
#'   )
#' }
#'
#' @author Ahmed El-Gabbas
#' @export

process_models <- function(
    model_dir = NULL, models_prefix = "models",
    climate_dir = NULL, resolution = 10L, verbose = TRUE,
    easin_ids = NULL, gbif_ids = NULL, coordinates = NULL,
    climate_scenario = "all", climate_model = "all", climate_year = "all",
    var_names = NULL, pft_type = "cross-walk", pft_id = NULL,
    bias_group = "none", vif_th = 10L, pres_buffer = 1000L, abs_buffer = 10L,
    abs_exclude_ext = list(), min_pres_grids = 50L,
    cv_folds = 5L, cv_block_size = 500L,
    cv_random = TRUE, abs_ratio = 20L, model_n_reps = 5L, n_cores = 4L,
    sdm_methods = c("glm", "glmnet", "gam", "maxent", "rf"),
    sdm_settings = NULL, proj_extent = "modelling_area") {

  down_clim <- down_lu <- species <- block_group <- cv <- n_pseudo_abs <- x <-
    training_abs <- training_pres <- test_abs <- test_pres <- valid <- label <-
    var_name <- out_file <- out_file_name <- year <- lu_file <- climate_file <-
    mod_method <- climate_model_abb <- climate_option <- packages <- y <-
    response_curves <- n_model_reps <- cv_summary <- model_rep_id <- value <-
    model_rep_results <- sdm_method <- variable_importance <- NULL

  ecokit::check_packages(
    c("usdm", "sdm", "parallelly", "future.apply", "future", "gtools", "fs",
      "tidyr", "withr", "sf", "terra", "crayon", "purrr", "tibble", "dplyr",
      "rworldmap", "rworldxtra", "tidyterra", "ragg", "ggplot2", "stringr",
      "grid", "ggtext", "scales"))

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
  bias_group <- tolower(bias_group)

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
  # Creating directories ------
  # # ********************************************************************** #

  if (is.null(models_prefix) || length(models_prefix) != 1L ||
      !nzchar(models_prefix)) {
    models_prefix <- "models"
  }
  models_prefix <- stringr::str_remove(models_prefix, "^_|_$") %>%
    stringr::str_trim() %>%
    stringr::str_to_lower()

  dir_models <- fs::path(
    model_dir, paste0(models_prefix, "_res_", resolution))
  dir_model_data <- fs::path(dir_models, "modelling_data")
  dir_fit <- fs::path(dir_models, "fitted_models")
  # directory path for projections of model repetitions
  dir_proj_reps <- fs::path(dir_models, "projections_reps")
  # directory path for summary of cross-validated projections
  dir_proj_cv <- fs::path(dir_models, "projections_cv")
  # directory path for final summary of projections
  dir_proj <- fs::path(dir_models, "projections")
  fs::dir_create(
    c(dir_models, dir_model_data, dir_fit, dir_proj_reps,
      dir_proj_cv, dir_proj))

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

    if (all(climate_scenario == "all")) {
      climate_scenario <- c("ssp126", "ssp370", "ssp585")
    }
    if (all(climate_model == "all")) {
      climate_model <- c("gfdl", "ipsl", "mpi", "mri", "ukesm1")
    }
    if (all(climate_year == "all")) {
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

  ## `sdm_methods` -----

  sdm_method_valid <- any(
    is.null(sdm_methods), length(sdm_methods) == 0L, !is.character(sdm_methods))
  if (sdm_method_valid) {
    ecokit::stop_ctx(
      "Invalid model method",
      sdm_methods = sdm_methods, class_sdm_method = class(sdm_methods),
      valid_sdm_methods = valid_sdm_methods)
  }

  sdm_methods <- sort(unique(tolower(sdm_methods)))
  valid_sdm_methods <- c(
    "glm", "glmpoly", "gam", "glmnet", "mars", "gbm", "rf", "ranger",
    "cart", "rpart", "maxent", "mlp", "svm", "svm2", "mda", "fda")
  if (!all(sdm_methods %in% valid_sdm_methods)) {
    invalid_sdm_methods <- setdiff(sdm_methods, valid_sdm_methods)
    ecokit::stop_ctx(
      paste0(
        "One or more model methods in `sdm_methods` are invalid.\n",
        "Invalid methods: ", toString(invalid_sdm_methods), "\n",
        "Valid methods are: ", toString(valid_sdm_methods)),
      sdm_methods = sdm_methods, invalid_sdm_methods = invalid_sdm_methods,
      valid_sdm_methods = valid_sdm_methods, cat_timestamp = FALSE)
  }

  sdm_packages <- tibble::tribble(
    ~mod_method, ~packages,
    "glm", NA_character_,
    "glmpoly", NA_character_,
    "gam", "mgcv",
    "glmnet", "glmnet",
    "mars", "earth",
    "gbm", "gbm",
    "rf", "randomForest",
    "ranger", "ranger",
    "cart", "tree",
    "rpart", "rpart",
    "maxent", c("dismo", "rJava"),
    "mlp", "RSNNS",
    # "rbf", "RSNNS",
    "svm", "kernlab",
    "svm", "e1071",
    "mda", "mda",
    "fda", "mda") %>%
    tidyr::unnest_longer("packages") %>%
    dplyr::filter(mod_method %in% sdm_methods, !is.na(packages))
  ecokit::check_packages(sdm_packages$packages)

  # Ensure that svm2 method is registered
  if ("svm2" %in% sdm_methods) {
    copy_svm2()
  }

  # # ********************************************************************** #

  ## `sdm_settings` -----

  if (is.null(sdm_settings)) {
    ecokit::cat_time("Loading default model settings")
    sdm_settings <- OneSDM::sdm_model_settings
  }

  if (all(sdm_methods %in% names(sdm_settings))) {
    sdm_settings <- sdm_settings[sdm_methods]
  } else {
    sdm_settings <- list()
  }


  # # ********************************************************************** #

  ## `proj_extent` -----

  if (is.null(proj_extent) ||
      !(is.character(proj_extent) || is.list(proj_extent))) {
    ecokit::stop_ctx(
      paste0(
        "The `proj_extent` argument must be either a character string of ",
        "length 1 or a list object of SpatExtent object(s)."),
      proj_extent = proj_extent, class_proj_extent = class(proj_extent),
      cat_timestamp = FALSE)
  }

  if (is.character(proj_extent)) {
    proj_extent <- tolower(proj_extent)
    valid_proj_extents <- c(
      "modelling_area", "modeling_area", "global", "europe")
    if (length(proj_extent) != 1L || !(proj_extent %in% valid_proj_extents)) {
      ecokit::stop_ctx(
        paste0(
          "Invalid `proj_extent` value. Valid options are: ",
          toString(valid_proj_extents), "."),
        proj_extent = proj_extent, cat_timestamp = FALSE)
    }
    if (proj_extent == "modeling_area") {
      proj_extent <- "modelling_area"
    }
  } else if (is.list(proj_extent)) {
    valid_extents <- purrr::map_lgl(proj_extent, inherits, "SpatExtent")
    if (!all(valid_extents)) {
      ecokit::stop_ctx(
        paste0(
          "All elements in the `proj_extent` argument must ",
          "be SpatExtent objects."),
        proj_extent = proj_extent,
        class_proj_extent = purrr::map(proj_extent, class),
        cat_timestamp = FALSE)
    }
  }

  # # ********************************************************************** #

  ## Parameters set via options -----

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

  species_pa_r <- terra::rast(species_pa_file) %>%
    stats::setNames("species")

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
    cat_timestamp = FALSE, verbose = verbose, level = 1L)

  # # ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||| #

  # # ********************************************************************** #
  # Defining study area and exclude pseudo-absences ------
  # # ********************************************************************** #

  ecokit::cat_time(
    "\nDefining study area and exclude pseudo-absences",
    cat_timestamp = FALSE, verbose = verbose)


  # Exclude pseudo-absences based on provided extents
  if (length(abs_exclude_ext) > 0L) {

    ecokit::cat_time(
      "Excluding pseudo-absences based on provided extents",
      cat_timestamp = FALSE, verbose = verbose, level = 1L)

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
      cat_timestamp = FALSE, verbose = verbose, level = 1L)
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
      cat_timestamp = FALSE, verbose = verbose, level = 1L)

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
      cat_timestamp = FALSE, verbose = verbose, level = 1L)
  }

  # # |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||| #

  # Exclude grid cells that are too close to presence locations from
  # pseudo-absences (< `abs_buffer` km). This is needed for sampling
  # pseudo-absences only beyond a certain distance from presences, which is
  # relevant for model types other than maxent.

  if (!is.null(abs_buffer) && is.numeric(abs_buffer) &&
      length(abs_buffer) == 1L && abs_buffer > 0L) {

    ecokit::cat_time(
      paste0(
        "Excluding pseudo-absences within ", ecokit::format_number(abs_buffer),
        " km of presence locations"),
      cat_timestamp = FALSE, verbose = verbose, level = 1L)

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
      cat_timestamp = FALSE, verbose = verbose, level = 1L)
  }

  # # ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||| #

  # # ********************************************************************** #
  # Save species presence-non-presence raster ------
  # # ********************************************************************** #

  ecokit::cat_time(
    "\nSaving species presence-non-presence raster",
    cat_timestamp = FALSE, verbose = verbose)
  species_pa_r_file <- fs::path(
    dir_models, paste0("species_pa_res_", resolution, ".tif"))
  terra::writeRaster(
    species_pa_r, filename = species_pa_r_file,
    overwrite = TRUE, gdal = c("COMPRESS=ZSTD", "ZSTD_LEVEL=22", "TILED=YES"))

  # # ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||| #

  # # ********************************************************************** #
  # Plotting species presence-non-presence raster ------
  # # ********************************************************************** #

  ecokit::cat_time(
    "Plotting species distribution", cat_timestamp = FALSE, verbose = verbose)

  species_pa_4_plot <- species_pa_r %>%
    terra::aggregate(
      fact = round(0.3333 / (resolution / 111.32)),
      fun = "max", na.rm = TRUE) %>%
    terra::as.factor() %>%
    terra::trim()
  p_ext <- terra::ext(species_pa_4_plot)
  file_plot <- fs::path(dir_model_data, "species_pa_distribution.jpeg")
  plot_width <- 18L
  aspect_ratio <- (p_ext[2L] - p_ext[1L]) / (p_ext[4L] - p_ext[3L])
  plot_height <- plot_width / (aspect_ratio * 0.975)
  world_outline <- sf::st_as_sf(rworldmap::getMap(resolution = "high"))

  ecokit::quietly({
    species_pa_plot <- ggplot2::ggplot() +
      ggplot2::geom_sf(
        data = world_outline, fill = "gray95", color = "black",
        linewidth = 0.05, show.legend = FALSE) +
      tidyterra::geom_spatraster(
        data = species_pa_4_plot, maxcell = 25e4L, na.rm = TRUE,
        mapping = ggplot2::aes(fill = ggplot2::after_stat(value))) +
      ggplot2::scale_fill_manual(
        values = c("0" = "grey70", "1" = "blue"), na.translate = FALSE,
        labels = c("0" = "potential pseudo-absences", "1" = "presences"),
        name = NULL, na.value = "transparent", drop = TRUE) +
      ggtext::geom_richtext(
        inherit.aes = FALSE, hjust = 0L, vjust = 1L, size = 3L,
        data = tibble::tibble(
          x = -Inf, y = Inf,
          label = paste0(
            "**Species presence-potential pseudo-absences**<br/>",
            "(aggregated to 1/3 degree resolution for visualization)")), #nolint
        mapping = ggplot2::aes(x = x, y = y, label = label),
        fill = scales::alpha("lightblue", 0.4), label.color = NA) +
      ggplot2::geom_sf(
        data = world_outline, fill = "transparent", color = "black",
        linewidth = 0.05, show.legend = FALSE) +
      ggplot2::coord_sf(
        xlim = p_ext[1L:2L], ylim = p_ext[3L:4L], expand = FALSE) +
      ggplot2::theme(
        text = ggplot2::element_text(size = 6L),
        axis.title = ggplot2::element_blank(),
        axis.ticks = ggplot2::element_blank(),
        axis.text.x = ggplot2::element_text(
          size = 5L, margin = ggplot2::margin(t = 0L, b = -0.5), vjust = 0L),
        axis.text.y = ggplot2::element_text(
          size = 5L, margin = ggplot2::margin(r = -0.25, l = -0.5), hjust = 0L),
        panel.grid = ggplot2::element_line(linewidth = 0.1, color = "gray85"),
        panel.background = ggplot2::element_rect(fill = "white"),
        plot.margin = grid::unit(c(-2L, 1L, 1L, 1L), "mm"),
        plot.background = ggplot2::element_blank(),
        legend.position = "bottom",
        legend.key.size = ggplot2::unit(0.2, "cm"),
        legend.margin = ggplot2::margin(t = -7L, b = -7L),
        legend.text = ggplot2::element_text(size = 6L))
  },
  "resampled to")

  ragg::agg_jpeg(
    filename = file_plot, width = plot_width, height = plot_height,
    res = 600L, quality = 100L, units = "cm")
  print(species_pa_plot)
  grDevices::dev.off()

  rm(
    species_pa_4_plot, p_ext, species_pa_plot, file_plot, plot_width,
    aspect_ratio, plot_height, envir = environment())
  invisible(gc())

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

  ecokit::cat_time(
    "Processing sampling effort predictor",
    cat_timestamp = FALSE, verbose = verbose)

  if (bias_group == "none") {

    bias_as_predictor <- FALSE
    bias_fix_value <- NA_real_
    ecokit::cat_time(
      "No sampling effort predictor will be used.",
      cat_timestamp = FALSE, verbose = verbose, level = 1L)

  } else {

    bias_as_predictor <- TRUE
    ecokit::cat_time(
      paste0(
        "Using sampling efforts for '", crayon::blue(bias_group),
        "' group as predictor"),
      cat_timestamp = FALSE, verbose = verbose, level = 1L)

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
    ecokit::cat_time(
      paste0(
        "90th percentile value of sampling effort predictor: ",
        ecokit::format_number(bias_fix_value)),
      cat_timestamp = FALSE, verbose = verbose, level = 1L)

    bias_fixed_r <- stats::setNames(bias_r, "bias_fixed")
    bias_fixed_r[bias_fixed_r <= bias_fixed_r] <- bias_fix_value

    # Add bias raster to predictors
    model_predictors <- c(model_predictors, bias_r, bias_fixed_r)

    # |||||||||||||||||||||||||||||||||||||||||||| #
    ## Plotting sampling effort predictor ----
    # |||||||||||||||||||||||||||||||||||||||||||| #

    ecokit::cat_time(
      "Plotting sampling effort predictor", cat_timestamp = FALSE,
      verbose = verbose, level = 1L)

    bias_r_4_plot <- terra::trim(bias_r)
    p_ext <- terra::ext(bias_r_4_plot)
    file_plot <- fs::path(dir_model_data, "sampling_efforts_log.jpeg")
    plot_width <- 18L
    aspect_ratio <- (p_ext[2L] - p_ext[1L]) / (p_ext[4L] - p_ext[3L])
    plot_height <- plot_width / (aspect_ratio * 0.97)

    ecokit::quietly({
      bias_plot <- ggplot2::ggplot() +
        ggplot2::geom_sf(
          data = world_outline, fill = "gray95", color = "black",
          linewidth = 0.05, show.legend = FALSE) +
        tidyterra::geom_spatraster(
          data = bias_r_4_plot, maxcell = 25e4L, na.rm = TRUE) +
        ggplot2::scale_fill_viridis_c(
          option = "magma", na.value = "transparent", name = NULL) +
        ggplot2::geom_sf(
          data = world_outline, fill = "transparent", color = "black",
          linewidth = 0.05, show.legend = FALSE) +
        ggplot2::coord_sf(
          xlim = p_ext[1L:2L], ylim = p_ext[3L:4L], expand = FALSE) +
        ggtext::geom_richtext(
          inherit.aes = FALSE, hjust = 0L, vjust = 1L, size = 3L,
          data = tibble::tibble(
            x = -Inf, y = Inf,
            label = "Sampling efforts (log<sub>10</sub> scale)"),
          mapping = ggplot2::aes(x = x, y = y, label = label),
          fill = scales::alpha("lightblue", 0.4), label.color = NA) +
        ggplot2::theme(
          text = ggplot2::element_text(size = 6L),
          axis.title = ggplot2::element_blank(),
          axis.ticks = ggplot2::element_blank(),
          axis.text.x = ggplot2::element_text(
            size = 5L, margin = ggplot2::margin(t = 0L, b = -0.5), vjust = 0L),
          axis.text.y = ggplot2::element_text(
            size = 5L, margin = ggplot2::margin(r = -0.25, l = -0.5),
            hjust = 0L),
          panel.grid = ggplot2::element_line(linewidth = 0.1, color = "gray85"),
          panel.background = ggplot2::element_rect(fill = "white"),
          plot.margin = grid::unit(c(-2L, 1L, 1L, 1L), "mm"),
          plot.background = ggplot2::element_blank(),
          legend.position = "bottom",
          legend.key.width = ggplot2::unit(0.4, "cm"),
          legend.key.height = ggplot2::unit(0.2, "cm"),
          legend.key.spacing = ggplot2::unit(0.05, "cm"),
          legend.margin = ggplot2::margin(t = -7L, b = -7L),
          legend.text = ggplot2::element_text(size = 4L),
          legend.title = ggtext::element_markdown(
            size = 6L, margin = ggplot2::margin(r = 4L)))
    },
    "resampled to")

    ragg::agg_jpeg(
      filename = file_plot, width = plot_width, height = plot_height,
      res = 600L, quality = 100L, units = "cm")
    print(bias_plot)
    grDevices::dev.off()

    rm(
      bias_r_4_plot, p_ext, bias_plot, file_plot, plot_width,
      bias_fixed_r, bias_r, aspect_ratio, plot_height, envir = environment())
    invisible(gc())

  }

  # |||||||||||||||||||||||||||||||||||||||||||| #
  # Saving masked predictors
  # |||||||||||||||||||||||||||||||||||||||||||| #

  ecokit::cat_time(
    "Saving masked predictors", cat_timestamp = FALSE, verbose = verbose)
  model_predictors_f <- fs::path(
    dir_model_data, paste0("model_predictors_res_", resolution, ".RData"))
  model_predictors <- terra::toMemory(model_predictors)
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
    "VIF-based predictor selection", cat_timestamp = FALSE, verbose = verbose)
  ecokit::cat_time(
    paste0("threshold = ", ecokit::format_number(vif_th)),
    cat_timestamp = FALSE, verbose = verbose, level = 1L)

  if (bias_as_predictor) {
    model_predictors_4vif <- terra::subset(
      x = model_predictors, subset = "bias_fixed", negate = TRUE)
  } else {
    model_predictors_4vif <- model_predictors
  }

  vif_results <- usdm::vifstep(
    x = model_predictors_4vif, th = vif_th, keep = vif_keep,
    size = vif_sample, method = "pearson")

  vif_file <- fs::path(dir_models, "vif_results.RData")
  save(vif_results, file = vif_file)
  excluded_vars <- vif_results@excluded

  rm(model_predictors_4vif, vif_results, envir = environment())
  invisible(gc())

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
      "No predictors were excluded based on VIF.",
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

  ecokit::info_chunk(
    "Future predictors", line_char_rep = 65L, verbose = verbose,
    cat_date = FALSE)

  if (make_future_predictions) {

    ecokit::cat_time(
      "future climate predictors", cat_timestamp = FALSE, verbose = verbose)

    climate_future <- tidyr::expand_grid(
      climate_scenario = climate_scenario,
      climate_model = climate_model,
      climate_year = climate_year)

    ecokit::cat_time(
      paste0(
        "Total climate scenario combinations to process: ",
        nrow(climate_future)),
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
              cat_timestamp = FALSE, verbose = verbose, level = 1L)

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

      ecokit::cat_time(
        "Future land use predictors", cat_timestamp = FALSE, verbose = verbose)

      lu_future <- tidyr::expand_grid(
        climate_scenario = climate_scenario,
        climate_year = climate_year)

      ecokit::cat_time(
        paste0(
          "Total land use scenario combinations to process: ", nrow(lu_future)),
        cat_timestamp = FALSE, verbose = verbose, level = 1L)

      lu_future <- lu_future %>%
        dplyr::mutate(
          down_lu = purrr::map2(
            .x = climate_scenario, .y = climate_year,
            .f = ~ {

              ecokit::cat_time(
                paste0(
                  "scenario: ", crayon::blue(.x), "; year: ", crayon::blue(.y)),
                cat_timestamp = FALSE, verbose = verbose, level = 1L)

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

  ecokit::info_chunk(
    "spatial block cross-validation", line_char_rep = 65L, verbose = verbose,
    cat_date = FALSE)

  # |||||||||||||||||||||||||||||||||||||||||||| #
  # Check if `cv_block_size` is appropriate for the filtered species data
  # |||||||||||||||||||||||||||||||||||||||||||| #

  ecokit::cat_time(
    "Checking spatial cross-validation block size",
    cat_timestamp = FALSE, verbose = verbose)

  # resolution in degrees
  resol <- terra::res(species_pa_r)[[1L]]
  # resolution in km
  resol_km <- resol * 111.320

  # extent width and height in km
  ext_dim <- ecokit::raster_dims_km(species_pa_r, exclude_na = TRUE)
  ext_w <- ext_dim$extent_width_occupied
  ext_h <- ext_dim$extent_height_occupied

  # Estimate number of blocks, ensure at least 1 block per direction
  n_blocks_x <- max(1L, floor(ext_w / cv_block_size))
  n_blocks_y <- max(1L, floor(ext_h / cv_block_size))
  potential_blocks <- n_blocks_x * n_blocks_y

  # If current `cv_block_size` cannot provide enough folds, reduce block size
  # appropriately
  if (potential_blocks < cv_folds) {

    ecokit::cat_time(
      paste0(
        "`cv_block_size` (", ecokit::format_number(cv_block_size),
        " km) is large to allow ", ecokit::format_number(cv_folds),
        " spatial folds"),
      cat_timestamp = FALSE, verbose = verbose, level = 1L)
    ecokit::cat_time(
      paste0(
        "Occupied extent (excluding NA rows and columns): ",
        ecokit::format_number(round(ext_w, 1L)),
        " [width] \u00D7 ", ecokit::format_number(round(ext_h, 1L)),
        " [height] km."),
      cat_timestamp = FALSE, verbose = verbose, level = 1L)

    # Compute block size so that area is evenly split among folds
    cv_block_size <- sqrt((ext_w * ext_h) / cv_folds)

    # Ensure block size is at least one raster cell (avoid too small block)
    if (cv_block_size < resol_km) {
      cv_block_size <- resol_km
    }

    ecokit::cat_time(
      paste0(
        "`cv_block_size` is set to ", ecokit::format_number(cv_block_size),
        " km to accommodate ", ecokit::format_number(cv_folds),
        " spatial folds."),
      cat_timestamp = FALSE, verbose = verbose, level = 1L)
  }

  # Further check if enough spatial blocks can be created with chosen block
  # size, by aggregating raster cells into blocks (to account for blocks with
  # only NA values)

  # Find aggregation factor (number of raster cells merged per block) based on
  # chosen block size and raster resolution
  agg_factor <- ceiling(cv_block_size / resol_km)

  # Aggregate raster to spatial blocks according to aggregation factor
  species_blocks <- terra::aggregate(
    x = species_pa_r, fact = agg_factor, fun = "sum", na.rm = TRUE)

  # Count how many aggregated blocks have data
  n_agg_cells <- terra::global(x = species_blocks, fun = "notNA")[[1L]] %>%
    as.integer()

  # If not enough spatial blocks remain after aggregating, recompute aggregation
  # factor and block size to allow required number of folds
  if (n_agg_cells < cv_folds) {

    ecokit::cat_time(
      paste0(
        "`cv_block_size` (", ecokit::format_number(cv_block_size),
        " km) is large to allow ", ecokit::format_number(cv_folds),
        " spatial folds"),
      cat_timestamp = FALSE, verbose = verbose, level = 1L)
    ecokit::cat_time(
      paste0(
        "Occupied extent (excluding NA rows and columns): ",
        ecokit::format_number(round(ext_w, 1L)),
        " [width] \u00D7 ", ecokit::format_number(round(ext_h, 1L)),
        " [height] km."),
      cat_timestamp = FALSE, verbose = verbose, level = 1L)

    # Use geometric means to set aggregation so that blocks = cv_folds
    agg_factor <- ceiling(sqrt((ext_w * ext_h) / (cv_folds * resol_km^2L)))

    # Update block size to match aggregation factor
    cv_block_size <- agg_factor * resol_km

    ecokit::cat_time(
      paste0(
        "`cv_block_size` is set to ", ecokit::format_number(cv_block_size),
        " km to accommodate ", ecokit::format_number(cv_folds),
        " spatial folds."),
      cat_timestamp = FALSE, verbose = verbose, level = 1L)
  }

  # |||||||||||||||||||||||||||||||||||||||||||| #
  ## Create spatial blocks -----
  # |||||||||||||||||||||||||||||||||||||||||||| #

  if (cv_random) {

    ecokit::cat_time(
      "Creating blocks using random blocks",
      cat_timestamp = FALSE, verbose = verbose)

    # aggregate species_pa_r to larger blocks, calculating sum of presences
    agg_factor <- ceiling(cv_block_size / (resol_km))
    ecokit::cat_time(
      paste0(
        "Aggregating presence-non-presence raster to blocks of ~",
        ecokit::format_number(cv_block_size), " km"),
      cat_timestamp = FALSE, verbose = verbose, level = 1L)
    ecokit::cat_time(
      paste0("Aggregation factor: ", ecokit::format_number(agg_factor)),
      cat_timestamp = FALSE, verbose = verbose, level = 1L)

    # convert aggregated blocks to polygons and assign cv folds. Folds are
    # assigned randomly but maintaining balance of presences across folds
    ecokit::cat_time(
      "Assigning CV folds to spatial blocks",
      cat_timestamp = FALSE, verbose = verbose, level = 1L)

    species_blocks <- terra::aggregate(
      x = species_pa_r, fact = agg_factor, fun = "sum", na.rm = TRUE) %>%
      terra::as.polygons(aggregate = FALSE) %>%
      sf::st_as_sf() %>%
      dplyr::arrange(dplyr::desc(species)) %>%
      dplyr::mutate(block_group = ceiling(dplyr::row_number() / cv_folds)) %>%
      dplyr::group_by(block_group) %>%
      dplyr::mutate(cv = sample.int(cv_folds, size = dplyr::n())) %>%
      dplyr::ungroup() %>%
      dplyr::select(-tidyselect::all_of("block_group")) %>%
      # Convert back to raster
      terra::vect() %>%
      terra::rasterize(y = species_blocks, field = "cv", background = NA) %>%
      # convert back to original resolution
      terra::disagg(fact = agg_factor, method = "near", progress = FALSE) %>%
      # mask by original buffered raster
      terra::mask(species_pa_r, progress = FALSE) %>%
      stats::setNames("cv")

    species_block_cv_file <- NA_character_

  } else {

    ecokit::cat_time(
      "Creating blocks using `blockCV` package",
      cat_timestamp = FALSE, verbose = verbose)

    ecokit::check_packages("blockCV")

    ecokit::cat_time(
      "Converting presence-non-presence raster to points",
      cat_timestamp = FALSE, verbose = verbose, level = 1L)
    species_pa_sf <- terra::as.points(species_pa_r) %>%
      sf::st_as_sf() %>%
      stats::setNames(c("species", "geometry"))

    ecokit::cat_time(
      paste0(
        "Creating blocks using `blockCV` package.",
        "\n  >>>  This could take a while, depending on the size ",
        "of the study area"),
      cat_timestamp = FALSE, verbose = verbose, level = 1L)

    species_blocks <- blockCV::cv_spatial(
      x = species_pa_sf, column = "species", r = species_pa_r, hexagon = FALSE,
      k = cv_folds, size = cv_block_size * 1000L, plot = FALSE,
      progress = verbose, report = verbose)

    species_block_cv_file <- fs::path(
      dir_models, paste0("species_cv_blockCV_", resolution, ".RData"))
    ecokit::cat_time(
      paste0(
        "Save output of `blockCV` blocks to: ",
        crayon::blue(species_block_cv_file)),
      cat_timestamp = FALSE, verbose = verbose, level = 1L)
    ecokit::save_as(
      object = species_blocks, object_name = "species_blocks",
      out_path = species_block_cv_file)

    # Extract raster of spatial blocks
    ecokit::cat_time(
      "Extracting blocks raster from `blockCV` output",
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

  # |||||||||||||||||||||||||||||||||||||||||||| #
  ## Plotting spatial blocks -----
  # |||||||||||||||||||||||||||||||||||||||||||| #

  ecokit::cat_time(
    "Plotting cross-validation blocks",
    cat_timestamp = FALSE, verbose = verbose)

  cv_r_4_plot <- terra::as.factor(terra::trim(species_blocks))
  p_ext <- terra::ext(cv_r_4_plot)
  file_plot <- fs::path(dir_model_data, "cv_blocks.jpeg")
  plot_width <- 18L
  aspect_ratio <- (p_ext[2L] - p_ext[1L]) / (p_ext[4L] - p_ext[3L])
  plot_height <- plot_width / (aspect_ratio * 0.97)

  species_pa_4_plot <- terra::classify(species_pa_r, cbind(0L, NA)) %>%
    terra::trim() %>%
    terra::as.points() %>%
    sf::st_as_sf()

  ecokit::quietly({
    cv_plot <- ggplot2::ggplot() +
      ggplot2::geom_sf(
        data = world_outline, fill = "gray95", color = "black",
        linewidth = 0.05, show.legend = FALSE) +
      tidyterra::geom_spatraster(
        data = cv_r_4_plot, maxcell = 25e4L, na.rm = TRUE,
        mapping = ggplot2::aes(fill = ggplot2::after_stat(value))) +
      ggplot2::scale_fill_discrete(
        na.value = "transparent", name = NULL, na.translate = FALSE) +
      ggplot2::geom_sf(
        data = species_pa_4_plot, fill = "gray97", color = "black",
        size = 0.125, shape = 19L, linewidth = 0.05, show.legend = FALSE) +
      ggplot2::geom_sf(
        data = world_outline, fill = "transparent", color = "black",
        linewidth = 0.05, show.legend = FALSE) +
      ggplot2::coord_sf(
        xlim = p_ext[1L:2L], ylim = p_ext[3L:4L], expand = FALSE) +
      ggtext::geom_richtext(
        inherit.aes = FALSE, hjust = 0L, vjust = 1L, size = 3L,
        data = tibble::tibble(
          x = -Inf, y = Inf, label = "Cross-validation folds"),
        mapping = ggplot2::aes(x = x, y = y, label = label),
        fill = scales::alpha("lightblue", 0.4), label.color = NA) +
      ggplot2::theme(
        text = ggplot2::element_text(size = 6L),
        axis.title = ggplot2::element_blank(),
        axis.ticks = ggplot2::element_blank(),
        axis.text.x = ggplot2::element_text(
          size = 5L, margin = ggplot2::margin(t = 0L, b = -0.5), vjust = 0L),
        axis.text.y = ggplot2::element_text(
          size = 5L, margin = ggplot2::margin(r = -0.25, l = -0.5), hjust = 0L),
        panel.grid = ggplot2::element_line(linewidth = 0.1, color = "gray85"),
        panel.background = ggplot2::element_rect(fill = "white"),
        plot.margin = grid::unit(c(-2L, 1L, 1L, 1L), "mm"),
        plot.background = ggplot2::element_blank(),
        legend.position = "bottom",
        legend.key.size = ggplot2::unit(0.35, "cm"),
        legend.key.spacing = ggplot2::unit(0.05, "cm"),
        legend.margin = ggplot2::margin(t = -7L, b = -7L),
        legend.text = ggplot2::element_text(size = 6L),
        legend.title = ggtext::element_markdown(
          size = 6L, margin = ggplot2::margin(r = 4L)))
  },
  "resampled to")

  ecokit::cat_time(
    paste0("Saving cross-validation blocks plot to: ", crayon::blue(file_plot)),
    cat_timestamp = FALSE, verbose = verbose, level = 1L)

  ragg::agg_jpeg(
    filename = file_plot, width = plot_width, height = plot_height,
    res = 600L, quality = 100L, units = "cm")
  print(cv_plot)
  grDevices::dev.off()

  rm(
    cv_r_4_plot, p_ext, cv_plot, file_plot, plot_width, species_pa_4_plot,
    ext_dim, ext_w, ext_h, n_blocks_x, n_blocks_y, n_agg_cells, resol, resol_km,
    aspect_ratio, plot_height, envir = environment())
  invisible(gc())

  # |||||||||||||||||||||||||||||||||||||||||||| #
  # Save spatial blocks to GeoTIFF file
  # |||||||||||||||||||||||||||||||||||||||||||| #

  species_blocks_file <- fs::path(
    dir_models, paste0("species_cv_blocks_", resolution, ".tif"))
  ecokit::cat_time(
    paste0(
      "Save spatial blocks to GeoTIFF file: ",
      crayon::blue(species_blocks_file)),
    cat_timestamp = FALSE, verbose = verbose)
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

  ecokit::cat_time(
    "Extracting modelling data from rasters",
    cat_timestamp = FALSE, verbose = verbose)

  model_data_r <- c(species_pa_r, model_predictors, species_blocks) %>%
    terra::toMemory()
  model_data <- as.data.frame(model_data_r, xy = TRUE, na.rm = TRUE) %>%
    tibble::tibble()
  model_data_file <- fs::path(
    dir_model_data, paste0("model_data_res_", resolution, ".RData"))
  ecokit::save_as(
    object = model_data, object_name = "model_data", out_path = model_data_file)

  # save model_data_r raster as wrapped terra SpatRaster, and load as needed to
  # free up memory
  model_data_r_file <- fs::path(dir_model_data, "model_data_r.RData")
  ecokit::save_as(
    object = terra::wrap(model_data_r), object_name = "model_data_r",
    out_path = model_data_r_file)
  rm(
    species_pa_r, model_data_r, species_blocks, model_predictors,
    envir = environment())
  invisible(gc())

  if (model_n_reps >= 1L) {
    pseudo_abs_check <- model_data %>%
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

  n_cores_modelling <- min(n_cores, cv_folds)
  ecokit::cat_time(
    paste0(
      "Preparing modelling data in parallel using ",
      ecokit::format_number(n_cores_modelling), " core(s)"),
    cat_timestamp = FALSE, verbose = verbose)

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
    "model_data", "dir_model_data", "abs_ratio", "bias_fix_value",
    "model_n_reps",  "model_data_r_file", "bias_as_predictor",
    "prep_model_data")

  models_cv <- future.apply::future_lapply(
    X = seq_len(cv_folds), FUN = prep_model_data,
    dir_model_data = dir_model_data, model_data = model_data,
    abs_ratio = abs_ratio, bias_fix_value = bias_fix_value,
    model_n_reps = model_n_reps, model_data_r_file = model_data_r_file,
    bias_as_predictor = bias_as_predictor,
    future.scheduling = Inf, future.seed = TRUE,
    future.packages = par_packages, future.globals = par_globals) %>%
    dplyr::bind_rows()

  # expand models_cv for all sdm_methods, except for maxent which will be
  # run only for a single repetition per CV fold
  models_cv_0 <- models_cv %>%
    tidyr::crossing(
      sdm_method = stringr::str_subset(sdm_methods, "maxent", negate = TRUE))

  if ("maxent" %in% sdm_methods) {
    models_cv_maxent <- models_cv %>%
      dplyr::filter(model_rep_id == 1L) %>%
      dplyr::mutate(
        n_pseudo_abs = NA_integer_, pseudo_abs_file = NA_character_,
        n_model_reps = 1L, sdm_method = "maxent")
    models_cv <- dplyr::bind_rows(models_cv_0, models_cv_maxent)
  } else {
    models_cv <- models_cv_0
  }

  models_cv <- dplyr::arrange(models_cv, cv, sdm_method, model_rep_id)

  future::plan("sequential", gc = TRUE)
  rm(
    models_cv_0, models_cv_maxent, par_globals, par_packages,
    envir = environment())
  invisible(gc())

  ## Saving model data
  models_cv_file <- fs::path(dir_models, "models_cv.RData")
  ecokit::cat_time(
    paste0("Saving model data to: ", crayon::blue(models_cv_file)),
    cat_timestamp = FALSE, verbose = verbose, level = 1L)
  ecokit::save_as(
    object = models_cv, object_name = "models_cv", out_path = models_cv_file)

  # # ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||| #

  # # ||||||||||||||||||||||||||||||||||||||||||||||||||||||| #
  # Projection data -------
  # # ||||||||||||||||||||||||||||||||||||||||||||||||||||||| #

  ecokit::info_chunk(
    "Preparing projection data",
    line_char_rep = 65L, verbose = verbose, cat_date = FALSE)

  # projection mask
  if (is.character(proj_extent)) {

    proj_mask <- switch(
      proj_extent,
      modelling_area = {
        ecokit::cat_time(
          "Models will be projected to the `modelling area`",
          cat_timestamp = FALSE, verbose = verbose)
        ecokit::load_as(model_data_r_file) %>%
          terra::unwrap() %>%
          terra::subset("species") %>%
          terra::classify(cbind(0L, 1L)) %>%
          stats::setNames("proj_extent")
      },
      global = {
        ecokit::cat_time(
          "Models will be projected `globally`",
          cat_timestamp = FALSE, verbose = verbose)
        OneSDM::get_mask_layer(
          resolution = resolution, climate_dir = climate_dir,
          europe_only = FALSE, return_spatraster = TRUE, wrap = FALSE) %>%
          stats::setNames("proj_extent")
      },
      europe = {
        ecokit::cat_time(
          "Models will be projected to `Europe`",
          cat_timestamp = FALSE, verbose = verbose)
        proj_mask <- OneSDM::get_mask_layer(
          resolution = resolution, climate_dir = climate_dir,
          europe_only = TRUE, return_spatraster = TRUE, wrap = FALSE) %>%
          stats::setNames("proj_extent")
      }
    )
  }

  if (is.list(proj_extent)) {
    ecokit::cat_time(
      "Models will be projected to custom extent(s)",
      cat_timestamp = FALSE, verbose = verbose)
    proj_mask <- OneSDM::get_mask_layer(
      resolution = resolution, climate_dir = climate_dir, europe_only = FALSE,
      return_spatraster = TRUE, wrap = FALSE)
    proj_mask <- purrr::map(
      .x = proj_extent,
      .f = ~ {
        proj_mask_0 <- proj_mask
        proj_mask_0[.x] <- 1000L
        proj_mask_0
      }) %>%
      terra::rast() %>%
      terra::app(fun = max, na.rm = TRUE) %>%
      terra::classify(cbind(1L, NA)) %>%
      terra::classify(cbind(1000L, 1L)) %>%
      stats::setNames("proj_extent") %>%
      terra::mask(proj_mask) %>%
      terra::trim()
  }

  proj_mask_file <- fs::path(dir_model_data, "projection_mask.tif")
  ecokit::cat_time(
    paste0("Saving projection mask to: ", crayon::blue(proj_mask_file)),
    cat_timestamp = FALSE, verbose = verbose, level = 1L)
  terra::writeRaster(
    proj_mask, filename = proj_mask_file,
    overwrite = TRUE, gdal = c("COMPRESS=ZSTD", "ZSTD_LEVEL=22", "TILED=YES"))
  rm(proj_mask, envir = environment())
  invisible(gc())

  ## Current ----
  ecokit::cat_time(
    "Preparing projection inputs for current climate",
    cat_timestamp = FALSE, verbose = verbose)
  ecokit::cat_time(
    "climate data", cat_timestamp = FALSE, verbose = verbose, level = 1L)
  projection_inputs <- dplyr::filter(
    climate_preds, var_name %in% predictor_names) %>%
    dplyr::pull(out_file)

  if (length(pft_id) > 0L) {
    ecokit::cat_time(
      "Land use data", cat_timestamp = FALSE, verbose = verbose, level = 1L)
    projection_inputs <- dplyr::filter(
      lu_preds, out_file_name %in% predictor_names) %>%
      dplyr::pull(out_file) %>%
      c(projection_inputs, .)
  }
  projection_inputs <- tibble::tibble(
    climate_model = "current", climate_scenario = "current",
    year = "current", map_paths = list(projection_inputs))

  ## Future -----
  if (make_future_predictions) {

    ecokit::cat_time(
      "Preparing projection inputs for future scenarios",
      cat_timestamp = FALSE, verbose = verbose)
    ecokit::cat_time(
      "climate data", cat_timestamp = FALSE, verbose = verbose, level = 1L)

    predictors_future <- dplyr::filter(
      climate_future, var_name %in% predictor_names) %>%
      dplyr::summarise(
        climate_file = list(out_file),
        .by = c(climate_model_abb, climate_scenario, year)) %>%
      dplyr::rename(climate_model = climate_model_abb)

    if (length(pft_id) > 0L) {
      ecokit::cat_time(
        "Land use data", cat_timestamp = FALSE, verbose = verbose, level = 1L)
      lu_future_files <- dplyr::filter(
        lu_future, out_file_name %in% predictor_names) %>%
        dplyr::summarise(
          lu_file = list(out_file), .by = c(climate_scenario, year))
      predictors_future <- dplyr::left_join(
        predictors_future, lu_future_files,
        by = c("climate_scenario", "year"))
    } else {
      predictors_future <- dplyr::mutate(
        predictors_future, lu_file = NA_character_)
    }

    predictors_future <- predictors_future %>%
      dplyr::mutate(
        map_paths = purrr::map2(.x = climate_file, .y = lu_file, .f = c),
        climate_file = NULL, lu_file = NULL) %>%
      dplyr::arrange(climate_model, climate_scenario, year)

    projection_inputs <- dplyr::bind_rows(
      projection_inputs, predictors_future) %>%
      dplyr::mutate(
        climate_option = paste0(
          climate_model, "_", climate_scenario, "_", year),
        climate_option = stringr::str_replace_all(
          climate_option, "current_current_current", "current"))

    rm(predictors_future, lu_future_files, envir = environment())
    invisible(gc())
  }

  # # ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||| #

  # # ********************************************************************** #
  # Model fitting and projections ------
  # # ********************************************************************** #

  ecokit::info_chunk(
    "Model fitting and projections",
    line_char_rep = 65L, verbose = verbose, cat_date = FALSE)

  ecokit::cat_time(
    paste0(
      "Total model fitting tasks: ",
      ecokit::format_number(nrow(models_cv), bold = TRUE)),
    cat_timestamp = FALSE, verbose = verbose)
  ecokit::cat_time(
    paste0(
      ecokit::format_number(length(sdm_methods)), " modelling methods \u00D7 ",
      ecokit::format_number(cv_folds), " CV folds \u00D7 ",
      ecokit::format_number(model_n_reps), " repetitions"),
    cat_timestamp = FALSE, verbose = verbose, level = 1L)

  if ("maxent" %in% sdm_methods) {
    ecokit::cat_time(
      "`Maxent` models will be fitted only for a single repetition per CV fold",
      cat_timestamp = FALSE, verbose = verbose, level = 1L)
  }
  ecokit::cat_time(
    paste0(
      "Modelling methods: ", crayon::blue(toString(sdm_methods))),
    level = 1L, cat_timestamp = FALSE, verbose = verbose)
  ecokit::cat_time(
    paste0(
      "Predictor variables: ",
      crayon::blue(toString(predictor_names))),
    level = 1L, cat_timestamp = FALSE, verbose = verbose)

  ecokit::cat_time(
    "Model projections", cat_timestamp = FALSE, verbose = verbose)
  n_projections <- ecokit::format_number(nrow(projection_inputs))
  ecokit::cat_time(
    paste0(n_projections, " projection options for each model fitting task"),
    cat_timestamp = FALSE, verbose = verbose, level = 1L)
  n_tiff_files <- ecokit::format_number(
    nrow(models_cv) * nrow(projection_inputs))
  ecokit::cat_time(
    paste0("Total projection GeoTIFF files: ", n_tiff_files),
    cat_timestamp = FALSE, verbose = verbose, level = 1L)

  n_cores_fit <- min(n_cores, nrow(models_cv))
  ecokit::cat_time(
    paste0(
      "Model fitting and post-processing using ",
      ecokit::format_number(n_cores_fit), " parallel core(s)"),
    cat_timestamp = FALSE, verbose = verbose)

  if (n_cores == 1L) {
    future::plan("sequential", gc = TRUE)
  } else {
    ecokit::set_parallel(n_cores_fit, show_log = FALSE, future_max_size = 2000L)
    withr::defer(future::plan("sequential", gc = TRUE))
  }

  par_packages <- c(
    "fs", "terra", "ecokit", "dplyr", "magrittr", "tidyselect",
    "sdm", "tibble", "purrr", "tidyr", "withr")
  par_globals <- c(
    "bias_fix_value", "bias_as_predictor", "projection_inputs",
    "predictor_names", "models_cv", "sdm_settings", "reduce_sdm_formulas",
    "dir_fit", "extract_sdm_info", "dir_proj_reps", "sdm_packages",
    "fit_predict", "proj_mask_file")

  model_fitting_pp <- future.apply::future_lapply(
    X = seq_len(nrow(models_cv)),
    FUN = function(cv_rep_id) {
      withr::with_envvar(
        new = c(NOAWT = "TRUE"),
        code = {
          suppressPackageStartupMessages(
            ecokit::quietly(
              fit_predict(
                cv_rep_id = cv_rep_id, bias_fix_value = bias_fix_value,
                bias_as_predictor = bias_as_predictor,
                projection_inputs = projection_inputs,
                predictor_names = predictor_names,
                models_cv = models_cv, sdm_settings = sdm_settings,
                dir_fit = dir_fit, dir_proj_reps = dir_proj_reps,
                sdm_packages = sdm_packages, proj_mask_file = proj_mask_file),
              "The following object is masked from",
              "Attaching package: ", "Loading required package",
              "Loaded ", "This version of "))
        })
    },
    future.scheduling = Inf, future.seed = TRUE,
    future.packages = par_packages, future.globals = par_globals)

  future::plan("sequential", gc = TRUE)
  invisible(gc())

  # # ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||| #

  # # ********************************************************************** #
  # Summarising model repetition outputs per modelling method and CV fold------
  # # ********************************************************************** #

  ecokit::info_chunk(
    "Summarise model outputs",
    line_char_rep = 65L, verbose = verbose, cat_date = FALSE)

  ecokit::cat_time(
    "Summarising model repetition outputs per modelling method and CV fold",
    cat_timestamp = FALSE, verbose = verbose)

  summ_fun <- function(.x, n_reps) { #nolint
    if (n_reps > 1L) {
      list(mean = ~ mean(.x, na.rm = TRUE), sd = ~ stats::sd(.x, na.rm = TRUE))
    } else {
      list(mean = ~ mean(.x, na.rm = TRUE), sd = ~ NA_real_)
    }
  }

  models_cv_summary <- models_cv %>%
    dplyr::mutate(
      model_rep_results = as.character(model_fitting_pp),
      model_fitting = purrr::map(model_rep_results, ecokit::load_as)) %>%
    tidyr::unnest(col = "model_fitting") %>%
    tidyr::nest(cv_summary = -c("cv", "sdm_method", "n_model_reps")) %>%
    dplyr::mutate(

      # Model repetition data
      reps_data = purrr::map(
        .x = cv_summary,
        .f = function(cv_data) {

          summary_columns_1 <- c(
            "variable_importance", "response_curves", "projections",
            "evaluation_training", "evaluation_testing", "pseudo_abs_file",
            "model_rep_id", "model_rep_results", "model_fit_file")
          df_1 <- cv_data %>%
            dplyr::select(-tidyselect::all_of(summary_columns_1)) %>%
            dplyr::distinct()

          summary_columns_2 <- c(
            "model_rep_id", "pseudo_abs_file", "model_rep_results",
            "model_fit_file", "evaluation_training", "evaluation_testing",
            "variable_importance", "response_curves", "projections")
          df_2 <- dplyr::rename(cv_data, rep_id = model_rep_id) %>%
            dplyr::summarise(
              dplyr::across(
                .cols = tidyselect::all_of(summary_columns_2[-1L]),
                .fns = ~ {
                  tibble::tibble(rep_id = rep_id, .x) %>%
                    tidyr::unnest(cols = .x) %>%
                    list()
                }))

          dplyr::bind_cols(df_1, df_2) %>%
            dplyr::rename_with(~paste0("rep_", .x))
        }),

      # Clean up some columns
      cv_summary = purrr::map(
        .x = cv_summary,
        .f = ~ {
          summary_columns_3 <- c(
            "evaluation_training", "evaluation_testing",
            "variable_importance", "response_curves", "projections")
          dplyr::select(.x, tidyselect::all_of(summary_columns_3))
        }),

      rep_summary = purrr::map2(
        .x = cv_summary, .y = n_model_reps,
        .f = function(cv_summary, n_reps) {

          # |||||||||||||||||||||||||||||||||||||||||||||||||
          # Train and test evaluation metrics
          # |||||||||||||||||||||||||||||||||||||||||||||||||

          cols_eval <- c("evaluation_training", "evaluation_testing")
          eval_train_test <- cv_summary %>%
            dplyr::summarize(
              dplyr::across(
                .cols = tidyselect::all_of(cols_eval),
                .fns = ~ {
                  dplyr::bind_rows(.x) %>%
                    dplyr::summarize(
                      dplyr::across(
                        .cols = tidyselect::everything(),
                        .fns = summ_fun(.x, n_reps),
                        .names = "{.col}_{.fn}")) %>%
                    list()
                }))

          # |||||||||||||||||||||||||||||||||||||||||||||||||
          # Variable importance
          # |||||||||||||||||||||||||||||||||||||||||||||||||

          var_imp <- cv_summary %>%
            dplyr::summarize(
              dplyr::across(
                .cols = tidyselect::all_of("variable_importance"),
                .fns = ~ {
                  dplyr::bind_rows(.x) %>%
                    dplyr::group_by(variable) %>%
                    dplyr::summarize(
                      dplyr::across(
                        .cols = tidyselect::all_of(c("cor_test", "auc_test")),
                        .fns = summ_fun(.x, n_reps),
                        .names = "{.col}_{.fn}"),
                      .groups = "drop") %>%
                    list()
                }))

          # |||||||||||||||||||||||||||||||||||||||||||||||||
          # Response curves
          # |||||||||||||||||||||||||||||||||||||||||||||||||

          # For response curves, there might be different x_value points across
          # repetitions. To properly summarize, we need to interpolate all reps
          # to a common grid per variable before calculating mean and sd.

          if (n_reps > 1L) {
            resp_curv <- cv_summary %>%
              dplyr::summarize(
                dplyr::across(
                  .cols = tidyselect::all_of("response_curves"),
                  .fns = ~ {

                    # Bind all response curve data
                    resp <- dplyr::bind_rows(.x)

                    # Build a per-variable grid from the overlap of all
                    # repetitions’ ranges
                    resp_grid <- dplyr::group_by(resp, variable) %>%
                      dplyr::summarise(
                        x_min = min(x_value), x_max = max(x_value),
                        .groups = "drop") %>%
                      dplyr::group_by(variable) %>%
                      dplyr::summarise(
                        x_min = x_min, x_max = x_max, n = 100L,
                        grid = list(seq(x_min, x_max, length.out = n)),
                        .groups = "drop") %>%
                      dplyr::select(tidyselect::all_of(c("variable", "grid")))

                    tidyr::nest(resp, rc_data = -variable) %>%
                      dplyr::left_join(resp_grid, by = "variable") %>%
                      tidyr::unnest(cols = "rc_data") %>%
                      dplyr::mutate(
                        x_value = purrr::map2_dbl(
                          .x = x_value, .y = grid,
                          .f = function(x_vals, grid) {
                            grid[which.min(abs(grid - x_vals))[1L]]
                          })) %>%
                      dplyr::select(-tidyselect::all_of("grid")) %>%
                      dplyr::group_by(variable, x_value) %>%
                      dplyr::summarize(
                        dplyr::across(
                          .cols = tidyselect::all_of("prediction"),
                          .fns = list(
                            mean = ~ mean(.x, na.rm = TRUE),
                            sd = ~ stats::sd(.x, na.rm = TRUE)),
                          .names = "{.col}_{.fn}"),
                        .groups = "drop") %>%
                      dplyr::mutate(
                        prediction_plus = pmin(
                          1.0, prediction_mean + prediction_sd),
                        prediction_minus = pmax(
                          0.0, prediction_mean - prediction_sd)) %>%
                      list()
                  }))
          } else {
            resp_curv <- dplyr::select(
              cv_summary, tidyselect::all_of("response_curves")) %>%
              dplyr::mutate(
                response_curves = purrr::map(
                  .x = response_curves,
                  .f = ~{

                    dplyr::mutate(
                      .x,
                      prediction_mean = prediction, prediction_sd = NA_real_,
                      prediction_plus = NA_real_,
                      prediction_minus = NA_real_) %>%
                      dplyr::select(-tidyselect::all_of("prediction"))
                  }))
          }
          dplyr::bind_cols(eval_train_test, var_imp, resp_curv)
        })) %>%
    tidyr::unnest("rep_summary")

  # |||||||||||||||||||||||||||||||||||||||||||||||||
  # Projections
  # |||||||||||||||||||||||||||||||||||||||||||||||||

  if (n_cores == 1L) {
    future::plan("sequential", gc = TRUE)
  } else {
    ecokit::set_parallel(
      min(n_cores, nrow(models_cv_summary)), show_log = FALSE,
      future_max_size = 2000L)
    withr::defer(future::plan("sequential", gc = TRUE))
  }

  par_packages <- c(
    "fs", "terra", "ecokit", "dplyr", "magrittr", "tidyselect",
    "sdm", "tibble", "purrr", "tidyr", "withr")
  par_globals <- c("summarise_preds_cv", "models_cv_summary", "dir_proj_cv")

  models_cv_summary_preds <- future.apply::future_lapply(
    X = seq_len(nrow(models_cv_summary)),
    FUN = summarise_preds_cv,
    models_cv_summary = models_cv_summary, dir_proj_cv = dir_proj_cv,
    future.scheduling = Inf, future.seed = TRUE,
    future.packages = par_packages, future.globals = par_globals)

  models_cv_summary <- models_cv_summary %>%
    dplyr::mutate(projections = models_cv_summary_preds) %>%
    dplyr::select(-tidyselect::all_of("cv_summary"))

  future::plan("sequential", gc = TRUE)
  rm(
    models_cv_summary_preds, par_globals, par_packages, model_fitting_pp,
    envir = environment())
  invisible(gc())

  # Saving cv summary data
  models_cv_summary_file <- fs::path(dir_models, "models_cv_summary.RData")
  ecokit::cat_time(
    paste0("Saving summary data to: ", crayon::blue(models_cv_summary_file)),
    cat_timestamp = FALSE, verbose = verbose, level = 1L)
  ecokit::save_as(
    object = models_cv_summary, object_name = "models_cv_summary",
    out_path = models_cv_summary_file)

  # # ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||| #

  # # ********************************************************************** #
  # Overall summary of cross-validated models ------
  # # ********************************************************************** #

  ecokit::cat_time(
    "Summarising model cross-validated models per modelling method",
    cat_timestamp = FALSE, verbose = verbose)

  models_summary <- models_cv_summary %>%
    dplyr::select(-tidyselect::all_of(c("cv", "n_model_reps", "reps_data"))) %>%
    dplyr::group_by(sdm_method) %>%
    dplyr::summarise(
      dplyr::across(
        .cols = tidyselect::everything(), .fns = ~ list(bind_rows(.x)))) %>%
    dplyr::mutate(

      # |||||||||||||||||||||||||||||||||||||||||||||||||
      # Train and test evaluation metrics
      # |||||||||||||||||||||||||||||||||||||||||||||||||

      dplyr::across(
        .cols = tidyselect::all_of(
          c("evaluation_training", "evaluation_testing")),
        .fns = ~ {
          dplyr::bind_rows(.x) %>%
            dplyr::summarize(
              dplyr::across(
                .cols = tidyselect::everything(),
                .fns = list(
                  mean = ~ mean(.x, na.rm = TRUE),
                  sd = ~ stats::sd(.x, na.rm = TRUE)),
                .names = "{.col}_{.fn}")) %>%
            dplyr::rename_with(~stringr::str_replace(.x, "_mean_", "_")) %>%
            list()
        }
      ),

      # |||||||||||||||||||||||||||||||||||||||||||||||||
      # Variable importance
      # |||||||||||||||||||||||||||||||||||||||||||||||||

      variable_importance = purrr::map(
        .x = variable_importance,
        .f = ~{

          dplyr::bind_rows(.x) %>%
            ecokit::arrange_alphanum(variable) %>%
            dplyr::group_by(variable) %>%
            dplyr::summarize(
              dplyr::across(
                .cols = tidyselect::ends_with("_mean"),
                .fns = list(
                  mean = ~ mean(.x, na.rm = TRUE),
                  sd = ~ stats::sd(.x, na.rm = TRUE)),
                .names = "{.col}_{.fn}"),
              .groups = "drop") %>%
            dplyr::rename_with(~stringr::str_replace(.x, "_mean_", "_"))
        }),

      # |||||||||||||||||||||||||||||||||||||||||||||||||
      # Response curves
      # |||||||||||||||||||||||||||||||||||||||||||||||||

      # For response curves, there might be different x_value points across
      # repetitions. To properly summarize, we need to interpolate all reps
      # to a common grid per variable before calculating mean and sd.

      response_curves = purrr::map(
        .x = response_curves,
        .f = ~ {

          # Bind all response curve data
          resp <- dplyr::bind_rows(.x)

          # Build a per-variable grid from the overlap of all CV fold ranges
          resp_grid <- dplyr::group_by(resp, variable) %>%
            dplyr::summarise(
              x_min = min(x_value), x_max = max(x_value),
              .groups = "drop") %>%
            dplyr::group_by(variable) %>%
            dplyr::summarise(
              x_min = x_min, x_max = x_max, n = 100L,
              grid = list(seq(x_min, x_max, length.out = n)),
              .groups = "drop") %>%
            dplyr::select(tidyselect::all_of(c("variable", "grid")))

          tidyr::nest(resp, rc_data = -variable) %>%
            dplyr::left_join(resp_grid, by = "variable") %>%
            tidyr::unnest(cols = "rc_data") %>%
            dplyr::mutate(
              x_value = purrr::map2_dbl(
                .x = x_value, .y = grid,
                .f = function(x_vals, grid) {
                  grid[which.min(abs(grid - x_vals))[1L]]
                })) %>%
            dplyr::group_by(variable, x_value) %>%
            dplyr::summarize(
              dplyr::across(
                .cols = tidyselect::all_of("prediction_mean"),
                .fns = list(
                  mean = ~ mean(.x, na.rm = TRUE),
                  sd = ~ stats::sd(.x, na.rm = TRUE)),
                .names = "{.col}_{.fn}"),
              .groups = "drop") %>%
            dplyr::rename_with(~stringr::str_replace(.x, "_mean_", "_")) %>%
            dplyr::mutate(
              prediction_plus = pmin(
                1.0, prediction_mean + prediction_sd),
              prediction_minus = pmax(
                0.0, prediction_mean - prediction_sd)) %>%
            ecokit::arrange_alphanum(variable, x_value)
        })
    )

  # |||||||||||||||||||||||||||||||||||||||||||||||||
  # Projections
  # |||||||||||||||||||||||||||||||||||||||||||||||||

  if (n_cores == 1L) {
    future::plan("sequential", gc = TRUE)
  } else {
    ecokit::set_parallel(
      min(n_cores, nrow(models_summary)), show_log = FALSE,
      future_max_size = 2000L)
    withr::defer(future::plan("sequential", gc = TRUE))
  }

  par_packages <- c(
    "fs", "terra", "ecokit", "dplyr", "magrittr", "tidyselect",
    "tibble", "purrr", "tidyr")
  par_globals <- c("summarise_preds", "models_summary", "dir_proj")

  models_summary_preds <- future.apply::future_lapply(
    X = seq_len(nrow(models_summary)),
    FUN = summarise_preds, models_summary = models_summary, dir_proj = dir_proj,
    future.scheduling = Inf, future.seed = TRUE,
    future.packages = par_packages, future.globals = par_globals)

  models_summary <- dplyr::mutate(
    models_summary, projections = models_summary_preds)

  future::plan("sequential", gc = TRUE)
  rm(models_summary_preds, par_globals, par_packages, envir = environment())
  invisible(gc())

  # Saving summary data
  models_summary_file <- fs::path(dir_models, "models_summary.RData")
  ecokit::cat_time(
    paste0("Saving summary data to: ", crayon::blue(models_summary_file)),
    cat_timestamp = FALSE, verbose = verbose, level = 1L)
  ecokit::save_as(
    object = models_summary, object_name = "models_summary",
    out_path = models_summary_file)

  # # ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||| #

  # # ********************************************************************** #
  # Overall modelling summary ------
  # # ********************************************************************** #

  ecokit::cat_time(
    "Prepare list of overall modelling summary (inputs and outputs)",
    cat_timestamp = FALSE, verbose = verbose)

  if (length(abs_exclude_ext) > 0L) {
    abs_exclude_ext <- purrr::map(abs_exclude_ext, terra::wrap)
  }
  if (is.list(proj_extent) && length(proj_extent) > 0L) {
    proj_extent <- purrr::map(proj_extent, terra::wrap)
  }

  model_summary <- list(
    # root directory for species models
    dir_models_root = model_dir,
    # model output directory
    dir_models = dir_models,

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

    # modelling options
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
      model_n_reps = model_n_reps,
      # sdm methods
      sdm_methods = sdm_methods,
      # sdm settings
      sdm_settings = sdm_settings),

    # modelling data
    model_data = list(
      # final presence-absence GeoTIFF file
      species_pa_file = species_pa_r_file,
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
      model_data_file = model_data_file,
      # paths for projection data
      projection_inputs = projection_inputs,
      # input projection extent
      proj_extent = proj_extent,
      # projection extent as PackedSpatRaster
      proj_mask_file = proj_mask_file,
      # model data for cross-validation folds and repetitions
      models_cv = models_cv),

    # Summary of model repetitions
    models_cv_summary = models_cv_summary,
    models_cv_summary_file = models_cv_summary_file,

    # Overall summary of models
    models_summary = models_summary,
    models_summary_file = models_summary_file
  )

  overall_summary_file <- fs::path(dir_models, "model_summary.RData")
  ecokit::cat_time(
    paste0(
      "Saving overall summary data to: ", crayon::blue(overall_summary_file)),
    cat_timestamp = FALSE, verbose = verbose, level = 1L)
  ecokit::save_as(
    object = model_summary, object_name = "model_summary",
    out_path = overall_summary_file)

  return(invisible(model_summary))

}
