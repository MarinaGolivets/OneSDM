#' Climate Data
#'
#' A dataset containing information about climate data files available for
#' download from OSF. See [prepare_climate()] for more details.
#'
#' @format ## `climate_data` A tibble 11 columns:
#' - `climate_scenario`: Character. The climate scenario: "`current`",
#'   "`ssp126`", "`ssp370`", and "`ssp585`".
#' - `climate_model`: Character. The climate model: "`current`", "`GFDL-ESM4`",
#'   "`IPSL-CM6A-LR`", "`MPI-ESM1-2-HR`", "`MRI-ESM2-0`", and "`UKESM1-0-LL`".
#' - `year`: Character. The year ranges: "`1981_2010`", "`2011_2040`",
#'   "`2041_2070`", and "`2071_2100`".
#' - `var_name`: Character. The variable name: "`bio1`", "`bio2`", "`bio3`",
#'   "`bio4`", "`bio5`", "`bio6`", "`bio7`", "`bio8`", "`bio9`", "`bio10`",
#'   "`bio11`", "`bio12`", "`bio13`", "`bio14`", "`bio15`", "`bio16`",
#'   "`bio17`", "`bio18`", "`bio19`", "`fcf`", "`fgd`", "`gdd0`", "`gdd5`",
#'   "`gdd10`", "`gddlgd0`", "`gddlgd5`", "`gddlgd10`", "`gdgfgd0`",
#'   "`gdgfgd5`", "`gdgfgd10`", "`gsl`", "`gsp`", "`gst`", "`lgd`", "`ngd0`",
#'   "`ngd5`", "`ngd10`", "`npp`", "`scd`", and "`swe`".
#' - `climate_name`: Character. The climate name; e.g., "`current`" and
#'   "<`year`>_<`climate_scenario`>_<`climate_model_abb`>" for future options
#'   (all small letters).
#' - `climate_model_abb`: Character. The abbreviated climate model name:
#'   "`current`", "`gfdl`", "`ipsl`", "`mpi`", "`mri`", and "`ukesm1`".
#' - `osf_path`: Character. The OSF path to the climate data file; e.g.,
#'   "`1981_2010_res_10/bio1.tif`".
#' - `file_size_mb`: Double. The size of the climate data file in megabytes.
#'  - `resolution`: Integer. The resolution of the climate data. Valid values
#'   are `5`, `10`, and `20`, for resolutions of approximately 10, 20, and 40 km
#'   (2.5, 5, and 10 arc-minutes) respectively.
#' - `out_dir`: Character. The output directory for the climate data file; e.g.
#'   "`res_5/1981_2010`". This would be relative to the `climate_dir` argument.
#' - `out_file`: Character. The output file path for the climate data file;
#'   e.g., "`res_5/1981_2010/bio1.tif`". This would be relative to the
#'   `climate_dir` argument.
#' @docType data
#' @examples
#' # View the climate data
#' data("climate_data")
#' head(climate_data)

"climate_data"


#' landuse_data Data
#'
#' A dataset containing information about land use data files available for
#' download from OSF. For more details, see [prepare_landuse()].
#'
#' @format ## `landuse_data` A tibble 10 columns:
#' - `pft_type`: Character. The type of plant functional type (PFT); one of
#'   "`original`" for the original PFTs as in the original land use data and
#'   "`cross-walk`" for the cross-walked PFTs as per the cross-walk table.
#'  - `resolution`: Integer. The resolution of the land use data: `5`, `10`,
#'   and `20`, for resolutions of approximately 10, 20, and 40 km (2.5, 5, and
#'   10 arc-minutes) respectively.
#'  - `climate_scenario`: Character. The climate scenario: "`current`",
#'   "`ssp126`", "`ssp370`", and "`ssp585`".
#'  - `year`: Character. The year ranges: "`1981_2010`", "`2011_2040`",
#'   "`2041_2070`", and "`2071_2100`".
#'  - `climate_name`: Character. The climate name: "`current`" and
#'   `"<year`>_<climate_scenario>"` for future options.
#'  - `id`: Integer. The PFT ID.
#'  - `pft`: Character. The plant functional type name.
#'  - `out_file`: Character. The output file path for the land use data
#'   files; e.g., "`res_20/ssp585_2071_2100/landuse_pft_cw_8_grasslands.tif`".
#'   This would be relative to the `landuse_dir` argument.
#'  - `osf_path`: Character. The OSF path to the land use data file;
#'   e.g., "`current_res_5/landuse_pft_1_water.tif`".
#'  - `out_dir`: Character. The output directory for the land use data file;
#'   e.g. "`res_5/current`". This would be relative to the `landuse_dir`
#'   argument.
#' @docType data
#' @examples
#' # View the land use data
#' data("landuse_data")
#' head(landuse_data)

"landuse_data"
