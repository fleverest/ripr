.onLoad <- function(libname, pkgname) {
  # Register the S7 classes and methods defined across the package.
  S7::methods_register()
}
