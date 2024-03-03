#' Run starr MCMC sampler
#'
#' This function calls the starr executable that is build through \code{compile_code()} function, and print
#' the results on the terminal for the parameters of the spatio-temporal AR(1) model. Additionally, it returns the
#' MCMC chains of the parameters.
#'
#' @param data a data frame that contains spatial and temporal information, with columns corresponding to a response variable, covariate variables, 2D geographical information (such as latitude, longitude), and location specification (can be numerical).
#' @param response A string containing the column name in \code{data} corresponding to the response variable.
#' @param covariates A vector of strings containing the column names in \code{data} that shall be included as covariates.
#' @param IDStations A string containing the column name in \code{data} that specifies the location corresponding to each row.
#' @param Time A string containing the column name in \code{data} that specifies the time variable. Can be in any sortable format.
#' @param x A string containing the column name in \code{data} that specifies the x-coordinate.
#' @param y A string containing the column name in \code{data} that specifies the y-coordinate.
#' @param nu A numeric value specifying nu for the Matérn correlation function. If not set, the default value is set to 0.5.
#' @param phi_prior a numeric vector of size 2 specifying the prior hyperparameters to be set for phi. If not set, the default value is
#' set to \code{c(2, 4)}.
#' @param n_iter An integer value specifying the number of iterations to run. If not set, the default is 5000.
#' @param burn_in A numeric value specifying number of burn-in iterations. If not set, the default  is 1000.
#'
#' @return A list with the following components:
#' \itemize{
#'   \item{\strong{beta_samples} returns a data matrix of size \code{n_covariates} x \code{n_iter}. The order of the rows correspond to the order specified by the user in \code{covariates}.}
#'   \item{\strong{rho_samples} returns a data vector of size \code{n_iter}.}
#'   \item{\strong{sigw_samples} returns a data vector of size \code{n_iter}.}
#'   \item{\strong{sigeps_samples} returns a data vector of size \code{n_iter}.}
#'   \item{\strong{sig0_samples} returns a data vector of size \code{n_iter}.}
#'   \item{\strong{phi_samples} returns a data vector of size \code{n_iter}.}
#'   \item{\strong{mu0_samples} returns a data matrix of size \code{number_locations} x \code{n_iter}. The rows are sorted by time, i.e., after sorting \code{data} by time, \code{mu0_samples[1]} will correspond to the first location
#'    appearing at time t = 1, \code{mu0_samples[2]} to the second location at time t=1, and so on.}
#'    \item{\strong{y_fitted} returns a data matrix of size \code{number_locations} x \code{T} of the fitted values. The rows are sorted by time, i.e., after sorting \code{data} by time, \code{y_fitted[1]} will correspond to the first location
#'    appearing at time t = 1, \code{y_fitted[2]} to the second location at time t=1, and so on.}
#' }
#'
#'@export
run_code <- function(data, response, covariates, IDStations, x,y ,Time , nu = 0.5, phi_prior = c(2, 4), n_iter = 5000, burn_in = 1000) {
  ## error messages for input
  for( i in 1:length(covariates)){
    if(!is.character(covariates[i])){stop("covariates must be a vector of strings!")}
    if(!is.numeric(data[,covariates[i]])){stop("Either covariate selection is wrong or data is not numeric!")}
  }
  if(!is.character(IDStations)){stop("IDStations must be a string!")}
  if(!is.character(x)){stop("x must be a string!")}
  if(!is.character(y)){stop("x must be a string!")}
  if(!is.character(Time)){stop("Time must be a string!")}
  if(!is.numeric(nu)){stop("nu must be numeric!")}
  if(length(phi_prior) != 2){stop("phi_prior must be vector of length 2!")}
  if(!is.numeric(phi_prior)){stop("phi_prior must be numeric!")}


  ##

  # get Renviron file from package
  renviron = system.file("starpackage.Renviron", package = "starr")
  # read keys/values
  readRenviron(renviron)
  #directory in which package is saved (get key PWD)
  home_dir = Sys.getenv("STARPACKAGE_HOME")
  if(home_dir == ""){
    stop("STARPACKAGE_HOME environment variable is not set")
  }

  ## make binary file via protobuf
  #source of proto files

  starpackage_home = dirname(dirname(home_dir))
  proto_path <- sprintf("%s/src/proto/parsedata.proto", starpackage_home)
  y_proto <- RProtoBuf::readProtoFiles(files = proto_path)

  ## N and T from data
  data[,IDStations] <- as.factor(data[,IDStations])
  N <- length(levels(data[,IDStations]))
  time_interval <- length(unique(data[,Time]))

  ## sort data according to Time
  sorted_data <- data[order(data[,Time]),]


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
    location$lat <- sorted_data[i, x];
    location$long <- sorted_data[i, y];

    msg_parsedata$loc[[i]] <- location;
  }

   #create binary file and write serialize msg_parsedata into it
  DATAPARSED_PATH <- sprintf("%s/src/dataparsed.bin", starpackage_home)
  writeBin(RProtoBuf::serialize(msg_parsedata, NULL), DATAPARSED_PATH)


  EXECUTE <- sprintf("cd %s/src/build && ./main %s/src/dataparsed.bin %s %s %s %s %s", starpackage_home,
                     starpackage_home, nu, phi_prior[1], phi_prior[2], n_iter, burn_in)
  system(EXECUTE)

  proto_path <- sprintf("%s/src/proto/paramdata.proto", starpackage_home)
  samples_proto <- RProtoBuf::readProtoFiles(files = proto_path)

  file_path_samples <- sprintf("%s/src/build/samples_serialized.bin", starpackage_home)
  binary_data_samples <- readBin(file_path_samples, "raw", file.info(file_path_samples)$size)

  msg_samples <- RProtoBuf::read(sampler_data.samples, binary_data_samples)
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

  ### fitted value parsing

  y_proto_path <- sprintf("%s/src/proto/ydata.proto", starpackage_home)
  y_proto <- RProtoBuf::readProtoFiles(files = y_proto_path)
  file_path_y <- sprintf("%s/src/build/y_serialized.bin", starpackage_home)
  binary_data_y <- readBin(file_path_y, "raw", file.info(file_path_y)$size)
  msg_y <- RProtoBuf::read(y_data.full_y, binary_data_y)
  ## COMPARING FITTED VALUES WITH DATA
  list_fitted_y <- as.list(msg_y$fitted_values)
  y_fitted <- list_fitted_y$vec_value
  y_fitted_mat <- matrix(y_fitted, nrow = N)


  alpha <- 0.05;

  ## beta CI and mean calculation
  beta_low_bound <- c();
  beta_high_bound <- c();
  beta_means <- c();

  for(i in 1:length(covariates)){
    beta_low_bound <- c(beta_low_bound, sort(beta[i,])[n_iter*alpha/2])
    beta_means <- c(beta_means, mean(beta[i,]))
    beta_high_bound <- c(beta_high_bound, sort(beta[i,])[n_iter*(1-alpha/2)])
  }


  CI <- data.frame(CI = paste("(", round(beta_low_bound,3), ", ", round(beta_high_bound,3), ")", sep = ""))
  beta_df <- data.frame(round(beta_means,3), CI)
  names(beta_df) <- c("Mean", "95% CI")
  row.names(beta_df) <- covariates

  ## variance components DF
  var_medians <- c(median(sig_w), median(sig_eps), median(sig_0))
  var_low_bound <- c(sort(sig_w)[n_iter*alpha/2],sort(sig_eps)[n_iter*alpha/2], sort(sig_0)[n_iter*alpha/2])
  var_high_bound <- c(sort(sig_w)[n_iter*(1-alpha/2)]  ,sort(sig_eps)[n_iter*(1-alpha/2)], sort(sig_0)[n_iter*(1-alpha/2)])

  CI_var <- data.frame(CI = paste("(", round(var_low_bound,3), ", ", round(var_high_bound,3), ")", sep = ""))
  var_df <- data.frame(round(var_medians,3), CI_var)
  names(var_df) <- c("Median", "95% CI")
  row.names(var_df) <- c("sigw^2", "sigeps^2", "sig0^2")

   # beta_conf_int <- data.frame(beta_means, CI)
  cat("Covariate effects: \n")
  cat("================================================\n")
  print(beta_df)
  cat("================================================\n")
  cat("Variance components:\n")
  cat("================================================\n")
  print(var_df)
  cat("================================================\n")
  cat(paste("mean rho: ", round(mean(rho),3),"  mean phi: ", round(mean(phi),3), "\n", sep = ""))

  ## return MCMC chains list objects for access

  return(list(beta_samples = beta, rho_samples = rho, sigw_samples = sig_w,
              sigeps_samples = sig_eps, sig0_samples = sig_0,
              phi_samples = phi, mu0_samples = mu0, y_fitted = y_fitted_mat))

}

