#' @keywords internal
"_PACKAGE"

## usethis namespace: start
#' @import S7
#' @importFrom gdalraster GDALRaster
## usethis namespace: end
NULL

.onLoad <- function(libname, pkgname) {
  S7::methods_register()
}