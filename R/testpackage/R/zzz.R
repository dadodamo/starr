.onAttach <- function(...) {
  readRenviron(system.file("testpackage.Renviron", package = "testpackage"))
}

# Unset BAYESMIX_EXE variable on detaching
.onDetach <- function(...) {
  Sys.unsetenv("TESTPACKAGE_HOME")
  Sys.unsetenv("TESTPACKAGE_EXE")
}
