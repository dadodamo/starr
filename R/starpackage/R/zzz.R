.onAttach <- function(...) {
  readRenviron(system.file("starpackage.Renviron", package = "starr"))
}

# Unset BAYESMIX_EXE variable on detaching
.onDetach <- function(...) {
  Sys.unsetenv("STARPACKAGE_HOME")
  Sys.unsetenv("STARPACKAGE_EXE")
}
