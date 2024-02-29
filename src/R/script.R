library( "RProtoBuf")
setwd("/users/daniel/desktop/stAR-sampler")


##### gamma plotting #####
x = seq(0,100,0.01);
dev.off()

plot(x, dgamma(x,2,rate = 1), type = 'l', c(0,10))

##### new proto files ##### 

samples_proto <- RProtoBuf::readProtoFiles(files = "proto/paramdata.proto")

file_path_samples <- "cmake-build-debug/samples_serialized.bin"  
binary_data_samples <- readBin(file_path_samples, "raw", file.info(file_path_samples)$size)

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



beta_true <- c(
  0.52218,
  -1.09668,
  -0.415183,
  -1.27539,
  -0.032302
  
)
n_iter= 5000;
alpha <- 0.05;
beta_low_bound <- c();
beta_high_bound <- c();
beta_means <- c();
for(i in 1:5){
beta_low_bound <- c(beta_low_bound, sort(beta[i,])[n_iter*alpha/2]) 
beta_means <- c(beta_means, mean(beta[i,]))
beta_high_bound <- c(beta_high_bound, sort(beta[i,])[n_iter*(1-alpha/2)])
}
beta_conf_int <- data.frame(beta_low_bound, beta_means, beta_high_bound)
beta_conf_int
## beta only sample  
dev.off()
par(mfrow = c(5, 1))
par(mar = c(2, 4, 2, 1))
for (i in 1:5) {
  # Create the plot
  plot(beta[i, 1000:5000], type = 'l', main = "", col = 'grey', ylab = bquote(beta[.(i)]))
  mean(beta[i,])
  # Add the true beta line
  #lines(1:4000, rep(beta_true[i], 4000), col = "green", type = 'l', lwd = 2)  
}
hist(beta[1,], breaks = 500)
abline( v = beta_low_bound[1], col = "blue", lty = 2)
abline( v= beta_high_bound[1], col = "blue", lty = 2)
abline( v = beta_true[1], col = 'red', lty = 2)
dev.off()


###### rho



## read rho file
dev.off()
plot(rho[1000:5000], type = 'l', col = 'grey', ylim = c(0.45, 0.55), ylab = "", main = expression(rho), xlab = "Iteration")
lines(1:4000, rep(0.5, 4000), type = 'l', col = 'green', lwd = 2)
plot(rho, type = 'l')
sort(rho)[n_iter*alpha/2]

#mu_0


## read file

mu_0_true <- c(
  1.44256,
  0.540784,
  -0.00554277,
  -1.18828,
  -0.340807,
  0.815028,
  -1.37403,
  -0.486838,
  -0.272475,
  1.639
)
dev.off()
par(mfrow = c(5, 2))
par(mar = c(2, 4, 2, 1))
for (i in 1:10) {
  plot(mu0[i,1000:5000], ylab = bquote(paste(mu[0], "[", .(i), "]")) ,  type = 'l', col = 'grey');
  lines(1:4000, rep(mu_0_true[i], 4000), type = 'l',col = 'green', lwd = 2);
}
mu_low_bound <- c();
mu_high_bound <- c();
means_mu <- c();
for(i in 1:10){
  mu_low_bound <- c(mu_low_bound, sort(mu0[i,])[n_iter*alpha/2])
  means_mu <- c(means_mu, mean(mu0[i,]));
  mu_high_bound <- c(mu_high_bound, sort(mu0[i,])[n_iter*(1-alpha/2)])
}
mu_conf_int <- data.frame(mu_low_bound, means_mu, mu_high_bound)
mu_conf_int

##### plot of variance components

## sig_eps
dev.off()
plot(sig_eps[1000:5000], type = 'l', col = 'grey', main =  bquote(sigma[epsilon]^2), ylab = "", xlab = "Iteration")
lines(1:4000, rep(0.9, 4000), type = 'l', col = "green", lwd = 2)
mean(sig_eps)
hist(sig_eps, steps = 500)

## sig_w
dev.off()
plot(sig_w[1000:5000], type = 'l',main =  bquote(sigma[w]^2), col ='grey', ylim = c(0, 1) ,xlab = "Iteration" , ylab = "")
lines(1:4000, rep(0.1, 4000), type = 'l', col = "green")
plot(sig_w, type = 'l')
mean(sig_w)

## sig_0
dev.off()
plot(sig_0,col = 'grey' ,type = 'l', main = bquote(sigma[0]^2), xlab = "Iteration", ylab = "")
lines(1:4000, rep(1, 4000), type = 'l', col = "green")
mean(sig_0)
median(sig_0)
plot(sig_0, type = 'l')
plot(density(sig_0), main = bquote(sigma[0]^2), xlab = "", xlim = c(0,10), ylab = '', bty = "n", lwd = 1)
polygon(density(sig_0), col = rgb(0, 0, 0, alpha = 0.2))
abline( v = 1, col = "black", lty = 1, lwd = 2, )
sig0_low_bound <- sort(sig_0)[n_iter*alpha/2];
sig0_high_bound <- sort(sig_0)[n_iter*(1-alpha/2)]
abline(v = sig0_low_bound, col = "blue", lty  = 2, lwd = 2)
abline( v = sig0_high_bound, col = "blue", lty = 2, lwd = 2)

### phi
dev.off()
plot(phi, type = 'l', col = 'grey')
mean(phi)
phi_low_bound <- sort(phi)[n_iter*alpha/2];
phi_high_bound <- sort(phi)[n_iter*(1-alpha/2)]
hist(phi, breaks = 500)
abline(v = phi_low_bound, col = 'blue', lty = 2)
abline(v = phi_high_bound, col = 'blue', lty = 2)
#### o's
o0 <- matrix(0, ncol = 5000, nrow = 10);
col = 1;
for(i in seq(11,999811 , 200)){
  o0[,col] = o[,i];
  col = col + 1;
}
dev.off()
plot(o0[1,], type = "l")  
mean(o0[1,])
plot(o0[2,], type = "l")
mean(o0[2,])
plot(o0[3,], type = "l")
mean(o0[3,])
plot(o0[4,], type = "l")
mean(o0[4,])
plot(o0[5,], type = "l")
mean(o0[5,])
plot(o0[6,], type = "l")
mean(o0[6,])
plot(o0[7,], type = "l")
mean(o0[7,])
plot(o0[8,], type = "l")
mean(o0[8,])
plot(o0[9,], type = "l")
mean(o0[9,])
plot(o0[10,], type = "l")
mean(o0[10,])

mean <- c(
  0.790282,
  2.44915,
  3.55486,
  5.77145,
  5.61658,
  5.07802,
  7.88179,
  8.32403,
  9.91118,
  10.2209
)
mean2 <- c(
  0.799816,
  2.05917,
  0.767974,
  0.722412,
  0.146886,
  -1.35127,
  2.82488,
  -0.877237,
  -2.75092,
  1.8037
)
par(mfrow = c(5, 2))
par(mar = c(2, 4, 2, 1))
for (i in 1:10) {
  plot(o0[i,1000:5000], main = paste("o10 ", i) ,  type = 'l', col = 'grey');
  lines(1:4000, rep(mean2[i], 4000), type = 'l', col = "green");
}

