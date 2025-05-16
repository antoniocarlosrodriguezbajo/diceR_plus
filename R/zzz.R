# Clean up when package is unloaded
.onUnload <- function(libpath) {
  library.dynam.unload("diceRplus", libpath)
}
.onLoad <- function(libname, pkgname) {
  options(
    ufs.ahk_executable_path = "C:/Program Files/AutoHotkey/v2/AutoHotkey.exe",
    ufs.ahk_script_path = system.file("scripts", "send_input.ahk", package = "diceRplus")
  )
}

