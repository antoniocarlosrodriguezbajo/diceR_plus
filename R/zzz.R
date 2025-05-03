# Clean up when package is unloaded
.onUnload <- function(libpath) {
  library.dynam.unload("diceRplus", libpath)
}
