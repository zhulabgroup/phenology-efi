path<-paste0("archive/",today,"/")
dir.create(path)

cl <- makeCluster(5, outfile = "")
registerDoSNOW(cl)

epoch=1

# neighbors and lags
var_list <- c("gcc" 
              , "tmean"
              , "prcp"
              )
vars <- 1:length(var_list)

var_oi<-1

neighbors <- vector(mode = "list", length(vars))
for (i in 1:length(vars)) {
  if (var_list[i] == "doy") {
    neighbors[[i]] <- 1:1
  } # for doy
  else {
    neighbors[[i]] <- 1:1
  }
}

lags <- vector(mode = "list", length(vars))
for (i in 1:length(vars)) {
  if (var_list[vars[i]]  %in% c("doy", "gccmean")) {
    lags[[i]] <- 1:1
  } # for doy
  if (var_list[vars[i]]  =="gcc") {
    lags[[i]] <- list()#list(1,2,3,4,5)
    for (period in 1:23) {
      lags[[i]] <-rlist::list.append(lags[[i]] ,(period-1)*16+1:16)
    }
  }
  else {
    lags[[i]] <- list()#list(1,2,3,4,5)
    for (period in 1:8) {
      lags[[i]] <-rlist::list.append(lags[[i]] ,(period-1)*8+1:8)
    }
  }
}

# weeks <- vector(mode = "list", length(vars))
# for (i in 1:length(vars)) {
#   if (var_list[vars[i]]  %in% c("doy", "gccmean")) {
#     weeks[[i]] <- 1:1
#   } # for doy
#   if (var_list[vars[i]]  =="gcc") {
#     weeks[[i]] <- 1:23
#   }
#   else {
#     weeks[[i]] <- 1:4
#   }
# }

ndim <- 0
for (i in 1:length(vars)) {
  ndim <- ndim + length(neighbors[[i]]) * length(lags[[i]])
}

range_list<-vector(mode="list", length=length(var_list))
for (i in 1:length(var_list)) {
  var<-var_list[i]
  if (var=="gcc") {range_list[[i]]<-c(0.3, 0.5)}
  if (var=="gccmean") {range_list[[i]]<-c(0.3, 0.5)}
  if (var=="tmax") {range_list[[i]]<-c(-20, 40)}
  if (var=="tmin") {range_list[[i]]<-c(-30, 20)}
  if (var=="humi") {range_list[[i]]<-c(0, 1)}
  # if (var=="srad") {range_list[[i]]<-c(0, 700)}
  if (var=="prcp") {range_list[[i]]<-c(0, 100)}
}

# Initial  hyperparameters
phimin <- 1e-50
phimax <- 0.99
vemin <- 0.001
vemax <- 0.999
taumin <- 0.001
taumax <- 0.999#4.99
gammamin <- 1/30^2 # exp(-d^2*gamma) gamma smaller -> u more similar over space/time
gammamax <- 1/1^2
rhomin <- .0001 
rhomax <- 0.999

V_list<-vector(mode="list")
for (v in 1:length(vars)){
  for (i in 1:length(lags[[v]])) {
    V_list <- rlist::list.append(V_list,0.05*exp(-(max(lags[[v]][[i]])/16/23)^2)) 
  }
}
V_phi<-unlist(V_list)
priors <- list(
  E_phi = matrix(c(0 * rep(1, length = ndim))),
  # V_phi = matrix(0.05, ncol = 1, nrow = ndim),
  V_phi = V_phi, #informed prior
  # Other priors
  E_ve = 0,
  V_ve = 1,
  E_tau = 0,
  V_tau = 0.01,
  E_gamma = 1/10^2,
  V_gamma = 100,
  E_rho = 0.5,
  V_rho = 5
)

num_part <- 20

num_epoch <- 1 #20

basisnumber <-500 # I tried 500 but exceeded MEM capacity
