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

  # EXECUTE <- sprintf("cd %s/src/build && ./ar_gibbs %s/src/dataparsed.bin %s %s %s", testpackage_home, testpackage_home, nu, phi_prior[1], phi_prior[2])
  # system(EXECUTE)


  PROTO_PATH <- sprintf("%s/src/proto/paramdata.proto", testpackage_home)
  samples_proto <- RProtoBuf::readProtoFiles(files = PROTO_PATH)

  READ_PATH <- sprintf("%s/src/build/samples_serialized.bin", testpackage_home)
  binary_data_samples <- readBin(READ_PATH, "raw", file.info(READ_PATH)$size)

  msg_samples <- read(sampler_data.samples, binary_data_samples)
  list_rho <- as.list(msg_samples$rho)
  list_phi <- as.list(msg_samples$phi)
  list_sigma_w <- as.list(msg_samples$sigma_w)
  list_sigma_eps <- as.list(msg_samples$sigma_eps)
  list_sigma_0 <- as.list(msg_samples$sigma_0)
  list_beta <- as.list(msg_samples$beta)
  list_mu0 <- as.list(msg_samples$mu_0)
  #list_o <- as.list(msg_samples$o)

  beta <- sapply(list_beta, function(x){x$vec_value})
  mu0 <- sapply(list_mu0, function(x){x$vec_value})
  rho <- sapply(list_rho, function(x){x});
  phi <- sapply(list_phi, function(x){x});
  sig_w <- sapply(list_sigma_w, function(x){x});
  sig_eps <- sapply(list_sigma_eps, function(x){x});
  sig_0 <- sapply(list_sigma_0, function(x){x});
  #o <- sapply(list_o, function(x){x$vec});
  #o <- sapply(o, function(x){x$vec_value});
}

