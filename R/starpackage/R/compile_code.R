#' Builds the starr executable
#'
#' \code{compile_code()} will execute \code{cmake} and \code{make} commands to be able to execute
#' the C++ code in the \code{src} directory.
#'
#' @return No output if build is successful, it raises errors otherwise
#'
#' @export
compile_code <- function() {
  build_subdir <- "build"
  renviron = system.file("starpackage.Renviron", package = "starr") # get Renviron file from package
  readRenviron(renviron) # read keys/values
  home_dir = Sys.getenv("STARPACKAGE_HOME") #directory in which package is saved (get key PWD)
  if(home_dir == ""){
    stop("STARPACKAGE_HOME environment variable is not set")
  }
  starpackage_home = dirname(dirname(home_dir))
  build_dir = sprintf("%s/src/%s", starpackage_home, build_subdir)
  dir.create(build_dir, showWarnings = FALSE)
  #
  print(build_dir)
  # # Configure testpackage
  cat("*** Configuring starr package ***\n")
  flags = '-DCMAKE_BUILD_TYPE=Release -DCMAKE_CXX_FLAGS="-O3"'
  CONFIGURE = sprintf('cd %s && cmake .. %s && make ', build_dir, flags)
  system(CONFIGURE)
  cat("\n")
  cat("\n")
  cat("Successfully installed starr package and compiled C++ code \n")
}


