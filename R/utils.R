.onUnload <- function (libpath) {
  library.dynam.unload("fastWavelets", libpath)
}
