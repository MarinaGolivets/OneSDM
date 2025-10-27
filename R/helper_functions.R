# terra package options

#' @noRd
#' @keywords internal

terra_options <- function(temp_dir) {
  terra::terraOptions(
    # fraction of RAM terra may use (0-0.9)
    memfrac = 0.5,
    # (GB) below which mem is assumed available
    memmin = 1L,
    # (GB) cap for terra on this node (set to node RAM - margin) or NA
    # prefer on-disk intermediate files
    tempdir = temp_dir,
    # silence per-worker progress bars
    progress = 0L,
    # write temporary files to disk directly
    todisk = TRUE)
}

# # ********************************************************************** #
# # ********************************************************************** #


# is_valid_gbif_id ------

#' Check the validity of a GBIF taxon ID
#'
#' This internal function verifies whether a given GBIF taxon ID is valid by
#' attempting to retrieve its usage information via the [rgbif::name_usage]
#' function. If the ID is valid and recognized by GBIF, the function returns
#' `TRUE`; otherwise, it returns `FALSE`.
#'
#' @param gbif_id An integer or character string representing a GBIF taxon ID.
#'
#' @return Logical value: `TRUE` if the GBIF ID is valid, `FALSE` otherwise.
#'
#' @keywords internal
#' @noRd

is_valid_gbif_id <- function(gbif_id) {
  tryCatch({
    rgbif::name_usage(key = gbif_id)
    # If no error, it's valid
    TRUE
  }, error = function(e) {
    # If error, it's invalid
    FALSE
  })
}

# # ********************************************************************** #
# # ********************************************************************** #
