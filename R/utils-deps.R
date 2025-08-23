#' Silently Load a Package
#'
#' It attempts to load a specified package silently, without printing any startup messages.
#' If the package is not installed, it will stop the execution and provide an error message
#' instructing the user to install the required package.
#'
#' @param pkg A character string specifying the name of the package to load.
#'
#' @details
#' - This function suppresses the package startup messages using `suppressPackageStartupMessages()`.
#' - If the package is not installed, it triggers an error with a message to guide the user to install the package.
#'
#' @return NULL
#'
#' @export
.silently_require <- function(pkg) {
  suppressPackageStartupMessages({
    if (!requireNamespace(pkg, quietly = TRUE)) {
      stop(sprintf("This feature requires '%s'. Please install it.", pkg), call. = FALSE)
    }
  })
}
