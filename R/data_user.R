#' Process and Validate User Observation Coordinates
#'
#' This function processes a set of user-provided geographic coordinates,
#' validates their structure and values, converts them into a spatial object,
#' and saves the result to a specified model directory.
#'
#' @param coordinates A data frame or matrix containing longitude and latitude
#'   values. Must have exactly two columns. If not provided directly, the
#'   function will attempt to retrieve it from the `onestop_coordinates` option.
#' @param model_dir Character. Path to the modelling directory where data and
#'   fitted models will be saved. This can not be `NULL` and should be the same
#'   directory used for the same species data. This can also be set via the
#'   `onestop_model_dir` option.
#'
#' @details The function performs the following steps:
#' - Validates that `coordinates` has at least one row and exactly two columns.
#' - Ensures that the columns represent numeric longitude and latitude values
#'   within valid ranges (longitude: `[-180, 180]`, latitude: `[-90, 90]`).
#' - Converts the coordinates to a tibble with standardized column names.
#' - Removes rows with missing or out-of-range values.
#' - Converts the cleaned coordinates to an `sf` spatial object (WGS84,
#'   EPSG:4326).
#' - Saves the resulting spatial object as `data/user_coordinates.RData` in
#'   `model_dir`.
#'
#' @return (Invisibly) The file path to the saved `user_coordinates.RData` file.
#'
#' @examples
#' \dontrun{
#'   coords <- data.frame(lon = c(10, 20), lat = c(50, 60))
#'   process_user_data(coords, model_dir = "path/to/model_dir")
#' }
#'
#' @export
#' @author Ahmed El-Gabbas

process_user_data <- function(coordinates = NULL, model_dir = NULL) {

  longitude <- latitude <- NULL

  # # ********************************************************************** #
  # Assigning function arguments from options if not provided directly ------
  # # ********************************************************************** #

  coordinates <- ecokit::assign_from_options(
    coordinates, "onestop_coordinates", c("data.frame", "matrix"))
  model_dir <- ecokit::assign_from_options(
    model_dir, "onestop_model_dir", "character")

  # # ********************************************************************** #
  # Checking function arguments -------
  # # ********************************************************************** #

  # # ||||||||||||||||||||||||||||||||||||||||| #
  ## coordinates ------
  # # ||||||||||||||||||||||||||||||||||||||||| #

  if (is.null(coordinates)) {
    return(invisible(NULL))
  }

  if (!inherits(coordinates, "data.frame")) {
    coordinates <- as.data.frame(coordinates)
  }

  if (nrow(coordinates) == 0L) {
    ecokit::stop_ctx(
      "The coordinates argument must contain at least one row.",
      coordinates = coordinates, cat_timestamp = FALSE)
  }
  if (ncol(coordinates) < 2L) {
    ecokit::stop_ctx(
      "The coordinates argument must contain at least two columns.",
      coordinates = coordinates, cat_timestamp = FALSE)
  }
  if (ncol(coordinates) > 2L) {
    ecokit::stop_ctx(
      "The coordinates argument contains more than two columns. ",
      coordinates = coordinates, cat_timestamp = FALSE)
  }

  coordinates <- tibble::tibble(coordinates) %>%
    stats::setNames(c("longitude", "latitude"))

  if (!is.numeric(coordinates$longitude) || !is.numeric(coordinates$latitude)) {
    ecokit::stop_ctx(
      paste0(
        "The longitude and latitude columns in the coordinates argument",
        "must be numeric."),
      coordinates = coordinates, cat_timestamp = FALSE)
  }

  coordinates <- dplyr::filter(
    coordinates, !is.na(longitude) & !is.na(latitude),
    longitude >= -180L, longitude <= 180L,
    latitude >= -90L, latitude <= 90L)

  if (nrow(coordinates) == 0L) {
    ecokit::stop_ctx(
      paste0(
        "The `coordinates` argument does not contain any valid rows. ",
        "Please ensure that the longitude and latitude columns contain ",
        "numeric values within the ranges: longitude [-180, 180], latitude ",
        "[-90, 90]."),
      coordinates = coordinates, cat_timestamp = FALSE)
  }

  # # ||||||||||||||||||||||||||||||||||||||||| #
  ## model_dir ------
  # # ||||||||||||||||||||||||||||||||||||||||| #

  if (is.null(model_dir)) {
    ecokit::stop_ctx(
      paste0(
        "The model_dir argument must be provided either directly or via the ",
        "`onestop_model_dir` option."),
      cat_timestamp = FALSE)
  }

  path_data <- fs::path(model_dir, "data")
  fs::dir_create(path_data)
  path_user_coordinates <- fs::path(path_data, "user_coordinates.RData")

  # # ********************************************************************** #
  # Processing user coordinate data
  # # ********************************************************************** #

  user_coordinates <- sf::st_as_sf(
    coordinates, coords = c("longitude", "latitude"),
    crs = 4326L, remove = FALSE)

  # Save processed coordinates
  save(user_coordinates, file = path_user_coordinates)

  return(invisible(path_user_coordinates))

}
