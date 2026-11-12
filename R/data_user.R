#' @title Process and Validate User-Provided Occurrence Data
#'
#' @description Validates user-provided geographic coordinates, converts them to an `sf`
#' object using WGS84 (EPSG:4326), and saves the result in the specified model directory.
#'
#' @param coordinates #' @param coordinates (data.frame or matrix) Two-column object
#'   containing longitude and latitude values. If NULL, coordinates are retrieved from
#'   `onesdm_coordinates`option. \strong{Required}.
#' @param model_dir (character) Directory for saving model outputs.
#'   A `data` subdirectory is created automatically to store raw and processed data.
#'   If `NULL` (default), path is retrieved from the `onesdm_model_dir` option.
#'   \strong{Required}.
#' @param return_data (logical) If `TRUE`, returns processed GBIF data as an
#'   `sf` object in addition to saving it. Defaults to `FALSE`.
#'
#' @details The function performs the following steps:
#' - Validates that `coordinates` contains at least one row and exactly two columns.
#' - Checks that longitude and latitude columns are numeric and within valid
#'   ranges (longitude: `[-180, 180]`, latitude: `[-90, 90]`).
#' - Converts input to a tibble with standardized column names.
#' - Removes rows with missing or out-of-range values.
#' - Converts cleaned coordinates to an `sf` object using WGS84 (EPSG:4326).
#' - Saves the spatial object as `data/user_coordinates.RData` in `model_dir`.
#' - Supports global configuration of default arguments via [base::options].
#'   Available options:
#'   - `onesdm_coordinates`: `coordinates`
#'   - `onesdm_model_dir`: `model_dir`.
#'
#' @return If `return_data = TRUE`, returns processed coordinates as an `sf`
#'   object. Otherwise, invisibly returns path to the saved data file.
#'
#' @examples
#' \dontrun{
#'   coords <- data.frame(lon = c(10, 20), lat = c(50, 60))
#'   prepare_user_data(coords, model_dir = "path/to/model_dir")
#' }
#'
#' @export
#' @author Ahmed El-Gabbas, Maryna Golivets

prepare_user_data <- function(
  coordinates = NULL,
  model_dir = NULL,
  return_data = FALSE
) {
  longitude <- latitude <- NULL

  ecokit::check_packages(c("fs", "sf", "tibble"))

  # ---------------------------------------------------------------------------
  # Arguments from options
  # ---------------------------------------------------------------------------

  coordinates <- ecokit::assign_from_options(
    coordinates,
    "onesdm_coordinates",
    c("data.frame", "matrix")
  )
  model_dir <- ecokit::assign_from_options(
    model_dir,
    "onesdm_model_dir",
    "character"
  )

  # ---------------------------------------------------------------------------
  # Validate arguments
  # ---------------------------------------------------------------------------

  if (is.null(coordinates)) {
    return(invisible(NULL))
  }

  if (!is.data.frame(coordinates) && !is.matrix(coordinates)) {
    ecokit::stop_ctx(
      "The `coordinates` argument must be a data frame or matrix.",
      coordinates = coordinates,
      cat_timestamp = FALSE
    )
  }

  if (ncol(coordinates) != 2L) {
    ecokit::stop_ctx(
      "The `coordinates` argument must have exactly two columns: longitude and latitude.",
      coordinates = coordinates,
      cat_timestamp = FALSE
    )
  }

  if (nrow(coordinates) == 0L) {
    ecokit::stop_ctx(
      "The `coordinates` argument must contain at least one row.",
      coordinates = coordinates,
      cat_timestamp = FALSE
    )
  }

  if (
    !is.logical(return_data) ||
      length(return_data) != 1L ||
      is.na(return_data)
  ) {
    ecokit::stop_ctx(
      "The `return_data` argument must be a single, non-missing logical value.",
      return_data = return_data,
      cat_timestamp = FALSE
    )
  }

  if (
    is.null(model_dir) ||
      !is.character(model_dir) ||
      length(model_dir) != 1L ||
      is.na(model_dir) ||
      !nzchar(trimws(model_dir))
  ) {
    ecokit::stop_ctx(
      paste0(
        "The `model_dir` argument must be provided as a non-empty character ",
        "path, either directly or via the `onesdm_model_dir` option."
      ),
      model_dir = model_dir,
      cat_timestamp = FALSE
    )
  }

  # Standardize the coordinate column names.
  coordinates <- as.data.frame(coordinates)
  names(coordinates) <- c("longitude", "latitude")

  if (
    !is.numeric(coordinates$longitude) ||
      !is.numeric(coordinates$latitude) ||
      is.complex(coordinates$longitude) ||
      is.complex(coordinates$latitude)
  ) {
    ecokit::stop_ctx(
      paste0(
        "The longitude and latitude columns in `coordinates` must be ",
        "numeric (integer or double) values."
      ),
      coordinates = coordinates,
      cat_timestamp = FALSE
    )
  }

  valid_rows <- (is.finite(coordinates$longitude) &
    is.finite(coordinates$latitude) &
    coordinates$longitude >= -180 &
    coordinates$longitude <= 180 &
    coordinates$latitude >= -90 &
    coordinates$latitude <= 90)

  coordinates <- tibble::as_tibble(
    coordinates[valid_rows, , drop = FALSE]
  )

  if (nrow(coordinates) == 0L) {
    ecokit::stop_ctx(
      paste0(
        "The `coordinates` argument does not contain any valid rows. ",
        "Longitude must be within [-180, 180] and latitude within [-90, 90]."
      ),
      coordinates = coordinates,
      cat_timestamp = FALSE
    )
  }

  # ---------------------------------------------------------------------------
  # Paths
  # ---------------------------------------------------------------------------

  path_data <- fs::path(model_dir, "data")
  fs::dir_create(path_data)
  path_user_coordinates <- fs::path(path_data, "user_coordinates.RData")

  # ---------------------------------------------------------------------------
  # # Convert to sf and save
  # ---------------------------------------------------------------------------

  user_coordinates <- sf::st_as_sf(
    coordinates,
    coords = c("longitude", "latitude"),
    crs = 4326,
    remove = FALSE
  )

  save(user_coordinates, file = path_user_coordinates)

  if (return_data) {
    return(user_coordinates)
  }
  invisible(path_user_coordinates)
}
