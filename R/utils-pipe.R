## quiets concerns of R CMD check re: the .'s that appear in pipelines
if (getRversion() >= "2.15.1") {
  utils::globalVariables(".")
}

#' Pipe operator
#'
#' See \code{magrittr::\link[magrittr:pipe]{\%>\%}} for details.
#'
#' @name %>%
#' @rdname pipe
#' @keywords internal
#' @export
#' @importFrom magrittr %>%
#' @usage lhs \%>\% rhs
#' @param lhs A value or the magrittr placeholder.
#' @param rhs A function call using the magrittr semantics.
#' @return The result of calling `rhs(lhs)`.
NULL


#' @noRd
.onAttach <- function(...) {    #nolint

  # Retrieve the package version and date dynamically
  package_version_info <- utils::packageVersion("OneSDM")
  package_date_info <- utils::packageDescription("OneSDM")$Date

  # Display the startup message
  packageStartupMessage(
    "OneSDM v", package_version_info,
    " - Last updated on ", package_date_info, "\n")
}
