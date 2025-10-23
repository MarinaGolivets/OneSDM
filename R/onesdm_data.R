#' Climate Data
#'
#' A dataset containing information about climate data files available for
#' download from OSF.
#'
#' @format ## `climate_data` A tibble 9 columns:
#' - `climate_scenario`: Character. The climate scenario. Valid values:
#'   `current`, `ssp126`, `ssp370`, `ssp585`.
#' - `climate_model`: Character. The climate model. Valid values: `current`,
#'   `GFDL-ESM4`, `IPSL-CM6A-LR`, `MPI-ESM1-2-HR`, `MRI-ESM2-0`, and
#'   `UKESM1-0-LL`.
#' - `climate_model_abb`: Character. The abbreviated climate model name. Valid
#'   values are: `current`, `gfdl`, `ipsl`, `mpi`, `mri`, and `ukesm1`.
#' - `year`: Character. The year ranges. Valid values are: `1981-2010`,
#'   `2011-2040`, `2041-2070`, and `2071-2100`.
#' - `var_name`: Character. The variable name. Valid values: `bio1`, `bio2`,
#'   `bio3`, `bio4`, `bio5`, `bio6`, `bio7`, `bio8`, `bio9`, `bio10`, `bio11`,
#'   `bio12`, `bio13`, `bio14`, `bio15`, `bio16`, `bio17`, `bio18`, `bio19`,
#'   `fcf`, `fgd`, `gdd0`, `gdd5`, `gdd10`, `gddlgd0`, `gddlgd5`, `gddlgd10`,
#'   `gdgfgd0`, `gdgfgd5`, `gdgfgd10`, `gsl`, `gsp`, `gst`, `lgd`, `ngd0`,
#'   `ngd5`, `ngd10`, `npp`, `scd`, and `swe`.
#' - `climate_name`: Character. The climate name (e.g., `current` and
#'   <`year`>_<`climate_scenario`>_<`climate_model_abb`> for future options (all
#'   small letters)).
#' - `osf_path`: Character. The OSF path to the climate data file (e.g.,
#'   "1981_2010_res_10/bio1.tif").
#' - `file_size_mb`: Double. The size of the climate data file in megabytes.
#' - `resolution`: Integer. The resolution of the climate data file. Valid
#'   values are 5, 10, and 20.
#' @examples
#' #' # View the climate data
#' data("climate_data")
#' head(climate_data)

"climate_data"
