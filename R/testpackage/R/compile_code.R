#' @export
compile_code <- function() {
  build_subdir <- "build"
  renviron = system.file("testpackage.Renviron", package = "testpackage") # get Renviron file from package
  readRenviron(renviron) # read keys/values
  home_dir = Sys.getenv("TESTPACKAGE_HOME") #directory in which package is saved (get key PWD)
  if(home_dir == ""){
    stop("TESTPACKAGE_HOME environment variable is not set")
  }
  testpackage_home = dirname(dirname(home_dir))
  build_dir = sprintf("%s/src/%s", testpackage_home, build_subdir)
  dir.create(build_dir, showWarnings = FALSE)
  #
  print(build_dir)
  # # Configure testpackage
  cat("*** Configuring testpackage ***\n")
  flags = '-DCMAKE_BUILD_TYPE=Release -DCMAKE_CXX_FLAGS="-O3"'
  CONFIGURE = sprintf('cd %s && cmake .. %s && make ', build_dir, flags)
  system(CONFIGURE)
  cat("\n")
  cat("\n")

  # Set TESTPACKAGE_EXE environment variable
  #cat("*** Setting TESTPACKAGE_EXE environment variable ***\n")
  # write(x = sprintf('TESTPACKAGE_EXE=%s/ar_gibbs', build_dir), file = renviron, append = TRUE)
  #  cat("\n")
  # #
   cat("Successfully installed Testpackage and compiled C++ code \n")
}


