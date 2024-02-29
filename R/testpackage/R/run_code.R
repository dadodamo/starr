#' @export
run_code <- function(data, response, covariates, IDStations, x,y ,Time , nu = 0.5, phi_prior = c(2, 4) ) {
  ## error messages for input


  ##

  # get Renviron file from package
  renviron = system.file("testpackage.Renviron", package = "testpackage")
  # read keys/values
  readRenviron(renviron)
  #directory in which package is saved (get key PWD)
  home_dir = Sys.getenv("TESTPACKAGE_HOME")
  if(home_dir == ""){
    stop("TESTPACKAGE_HOME environment variable is not set")
  }

  ## make binary file via protobuf
  #source of proto files

  testpackage_home = dirname(dirname(home_dir))
  proto_path <- sprintf("%s/src/proto/parsedata.proto", testpackage_home)
  y_proto <- RProtoBuf::readProtoFiles(files = proto_path)

  ## N and T from data
  data$IDStations <- as.factor(data$IDStations)
  N <- length(levels(data$IDStations))
  time_interval <- length(unique(data$Time))

  ## sort data according to Time
  sorted_data <- data[order(data$Time),]


  ## creating message for data parsing
  msg_parsedata <- new(parsedata.input_data)

  ##hyperparam
  msg_parsedata$N <- N
  msg_parsedata$T <- time_interval

  ## response vector
  msg_parsedata$y <- new(parsedata.vector, vec_value =  as.vector(sorted_data[,response]))

  ## covariate vector
  for(i in 1:length(covariates)){
    temp <- new(parsedata.vector, vec_value = as.vector(sorted_data[,covariates[i]]))
    msg_parsedata$x$m_vec[[i]]<- temp
  }

  ## location vector
  for(i in 1:N){
    location <- new(parsedata.location)
    location$lat <- sorted_data$x[i];
    location$long <- sorted_data$y[i];

    msg_parsedata$loc[[i]] <- location;
  }

   #create temporary binary file and save it to working src directory
  DATAPARSED_PATH <- sprintf("%s/src/dataparsed.bin", testpackage_home)
  file.create(DATAPARSED_PATH)
  file <- file(DATAPARSED_PATH, open = "wb")
  serialize(msg_parsedata, file)
  close(file)

  # test if it works
  #binary_data <- readBin(DATAPARSED_PATH, "raw", file.info(DATAPARSED_PATH)$size)
  #msg_data <- read(parsedata.input_data, binary_data)
  #print(msg_data$y$vec_value)

   EXECUTE <- sprintf("cd %s/src/build && ./main %s/src/dataparsed.bin %s %s %s", testpackage_home, testpackage_home, nu, phi_prior[1], phi_prior[2])
   system(EXECUTE)
}

