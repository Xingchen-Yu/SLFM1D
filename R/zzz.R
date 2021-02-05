.onUnload <- function (libpath) {
  library.dynam.unload("SLFM1D", libpath)
}
.onAttach <- function(libname, pkgname) {
  packageStartupMessage("Hello")
}
