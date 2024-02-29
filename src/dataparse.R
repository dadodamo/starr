library( "RProtoBuf")
setwd("/users/daniel/desktop")
data <- read.csv("pm10_dataset.csv")

setwd("/users/daniel/desktop/stAR-sampler")

colnames(data)
## covariate names for later + selection
covariates <- c("WE_temp_2m", "WE_tot_precipitation", "WE_surface_pressure", "WE_wind_speed_100m_mean", "WE_blh_layer_max")
response <- "AQ_pm10"
y_proto <- RProtoBuf::readProtoFiles(files = "/users/daniel/desktop/stAR-sampler/proto/parsedata.proto")

oneloc <- data[data$IDStations == '1264',]
plot(scale(oneloc[,'WE_blh_layer_max']), type = 'l')
  #factor ID
data$IDStations <- as.factor(data$IDStations)
data$NomeStazione <- as.factor(data$NomeStazione)

## sorting data wrt time
data$Time <- as.Date(data$Time)
sorted_data <- data[order(data$Time),]
sorted_data <- sorted_data[1:(64*365),]
scaled_matrix <- scale(sorted_data[, c(response, "Latitude", "Longitude", covariates)])

sorted_data[, c(response, "Latitude", "Longitude", covariates)] <- scaled_matrix

## N and T
N <- length(levels(data$IDStations))
time_interval <- length(unique(data$Time))



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
  location$lat <- sorted_data[i,"Latitude"];
  location$long <- sorted_data[i,"Longitude"];

  msg_parsedata$loc[[i]] <- location;
}



#tf2 <- tempfile()
#con <- file( tf2, open = "wb" )
# serialize( msg_parsedata, con )
#file.copy(tf2, target_file_path2, overwrite = TRUE)
#close( con )

target <- '/Users/daniel/Desktop/rpackage/src/dataparsed.bin'
file.create(target)
con <- file(target, "wb")
ser_msg <- serialize(msg_parsedata, con)
# Write the serialized message to the file
close(con)

target_file_path <- "/users/daniel/Desktop/stAR-sampler/dataparsed.bin"
target_file_path2 <- "/Users/daniel/Desktop/rpackage/src/dataparsed.bin"


binary_data <- readBin(target, "raw", file.info(target)$size)
msg_data <- read(parsedata.input_data, binary_data)
msg_data$y$vec_value
msg_data$loc[[3]]
list <- as.list(msg_data$loc)
list <- sapply(list, function(x){x$lat})
