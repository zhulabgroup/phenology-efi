update<-F
source("code/11 get NEON phenology data.R")
source("code/12 get NEON weather data api.R")
source("code/13 get NOAA weather data.R")

source("code/21 preprocess data.R")

###
dir.create(paste0(path,"train_out"))
res<-PrepareEmbedding(x,start=max(unlist(lags))+1,end=length(date_list)-365, focalsites =1:(nrow(coord_df)/2), lags=lags, neighbors=neighbors,vars=vars, distMat = distMat)
X_train<-res$X
Y_train<-res$Y
D_train<-res$D
P_train<-res$P

write_csv(as.data.frame(X_train),paste0(path,"train_out/X.csv"))
write_csv(as.data.frame(Y_train),paste0(path,"train_out/Y.csv"))
write_csv(as.data.frame(D_train),paste0(path,"train_out/D.csv"))
write_csv(as.data.frame(P_train),paste0(path,"train_out/P.csv"))

######
path_data<-paste0(path,"train_out")
path_results<-paste0(path, "results_out")
dir.create(path_results, recursive = T)

########################
# Get new data
########################
X_all<-read_csv(list.files(path_data, pattern="X.csv", full.names = T)) %>% as.matrix()
Y_all<-read_csv(list.files(path_data, pattern="Y.csv", full.names = T))%>% as.matrix()
D_all<-read_csv(list.files(path_data, pattern="D.csv", full.names = T)) %>% as.matrix()
P_all<-read_csv(list.files(path_data, pattern="P.csv", full.names = T)) %>% as.matrix()

missing_id<-which(rowSums(is.na(cbind(Y_all,X_all)))!=0)
if (length(missing_id)>0) {
  X_all<-X_all[-missing_id,,drop=F]
  Y_all<-Y_all[-missing_id,,drop=F]
  P_all<-P_all[-missing_id,,drop=F]
  D_all<-D_all[-missing_id,,drop=F]
}

train_id <- sample(1:nrow(X_all), round(nrow(X_all) * 10 / 10))
# valid_id <- setdiff(1:nrow(X_all), train_id)

X_train <- X_all[train_id, , drop = F]
Y_train <- Y_all[train_id, , drop = F]
P_train <- P_all[train_id, , drop = F]
D_train <- D_all[train_id, , drop = F]

# X_valid <- X_all[valid_id, , drop = F]
# Y_valid <- Y_all[valid_id, , drop = F]
# P_valid <- P_all[valid_id, , drop = F]
# D_valid <- D_all[valid_id, , drop = F]

basisnumber<-min(basisnumber, nrow(X_all))
if (basisnumber==nrow(X_all)) {
  X_basis <- X_all
  Y_basis <- Y_all
  P_basis <- P_all
  D_basis <- D_all
} else {
  set.seed(42)
  cluster <- kmeans(cbind(P_all,X_all), basisnumber )
  cluster_id_sort <- seq(basisnumber )[order(cluster$size, decreasing = T)]
  basis_id <- rep(NA, basisnumber)
  for (i in 1:basisnumber) {
    # basis_id<-which(cl$cluster==cluster_id_sort[i])[length(which(cl$cluster==cluster_id_sort[i]))]
    basis_id[i] <- sample(which(cluster$cluster == cluster_id_sort[i]), 1)
  }
  X_basis <- X_all[basis_id, , drop = F]
  Y_basis <- Y_all[basis_id, , drop = F]
  P_basis <- P_all[basis_id, , drop = F]
  D_basis <- D_all[basis_id, , drop = F]
}

write_csv(as.data.frame(X_basis), paste0(path_results, "/", "X_basis", ".csv"))
write_csv(as.data.frame(Y_basis), paste0(path_results, "/", "Y_basis", ".csv"))
write_csv(as.data.frame(P_basis), paste0(path_results, "/", "P_basis", ".csv"))
write_csv(as.data.frame(D_basis), paste0(path_results, "/", "D_basis", ".csv"))

if (epoch == 1) {
  pars_id_prev <- numeric()
  pars_prev_updated <- matrix(NA, nrow = ndim + 4, ncol = 0)
  pars_var_prev_updated <- vector(mode = "list", length = 0)
  log_p_prev_updated <- numeric(0)
}
if (epoch > 1) {
  ########################
  # Get previous hyperparameters
  ########################
  # pars_results_prev<-paste0(path, "/results/", "LC", lc_code, "/", "DOY", doy_oi, "/", "EPOCH",epoch-1)
  # change this
  
  pars_prev <- as.matrix(read_csv(list.files(pars_results_prev, pattern="pars", recursive = T)))
  pars_id_prev <- colnames(pars_prev)
  pars_var_prev <- vector(mode = "list", length = num_part)
  for (i in 1:num_part) {
    pars_var_prev[[i]] <- as.matrix(nearPD(as.matrix(read_csv(list.files(pars_results_prev, pattern=paste0("pars_var_", i), recursive = T))))$mat)
  }
  
  ########################
  # Learn hyperparameters
  ########################
  Rprop_res <-
    foreach(
      i = 1:num_part,
      # .export = ls(globalenv()),
      .export = c(
        "phimax", "phimin", "vemax", "vemin", "taumax", "taumin", "rhomax", "rhomin",
        "pars_prev", "pars_var_prev", "distMat", "priors",
        "X_train", "P_train", "Y_train",
        "X_basis", "P_basis", "Y_basis",
        "GPSDM", "fmingrad_Rprop"
      ),
      .packages = c("Matrix", "LaplacesDemon", "RhpcBLASctl")
    ) %dopar% {
      blas_set_num_threads(1)
      omp_set_num_threads(1)
      
      # subset<-sample(1:nrow(X_new), nrow(X_new))
      # X_new<-X_new[subset,,drop=F]
      # P_new<-P_new[subset,,drop=F]
      # Y_new<-Y_new[subset,,drop=F]
      
      lpost <- function(p,sample_n) {
        minibatch <- sample(1:nrow(X_train), sample_n)
        return(GPSDM(pars = p, distMat = distMat, basisX = X_basis, basisP = P_basis, basisY = Y_basis, newX = X_train[minibatch, , drop = F], newP = P_train[minibatch, , drop = F], newY = Y_train[minibatch, , drop = F], priors_tr = list(E = pars_prev[, i, drop = F], V = pars_var_prev[[i]]), mode = c("optimize")))
      }
      print("start")
      
      fmingrad_Rprop(lpost, pars_prev[, i, drop = F] + matrix(rnorm(ndim + 3, 0, 0.1)), sample_n = min(nrow(X_train),100), maxcount = 200)
    }
  
  ########################
  # Estimate hyperparameter uncertainty
  ########################
  pars_prev_updated <- matrix(NA, nrow = ndim + 3, ncol = num_part)
  log_p_prev_updated <- rep(NA, num_part)
  for (i in 1:num_part) {
    pars_prev_updated[, i] <- Rprop_res[[i]]$xopt
    log_p_prev_updated[i] <- -Rprop_res[[i]]$fopt
  }
  pars_var_prev_updated <-
    foreach(
      i = 1:num_part,
      .export = c(
        "phimax", "phimin", "vemax", "vemin", "taumax", "taumin", "gammamax", "gammamin",
        "pars_prev", "pars_var_prev", "pars_prev_updated", "distMat", "priors",
        "X_train", "P_train", "Y_train",
        "X_basis", "P_basis", "Y_basis",
        "GPSDM", "fmingrad_Rprop"
      ),
      .packages = c("Matrix", "LaplacesDemon", "RhpcBLASctl")
    ) %dopar% {
      blas_set_num_threads(1)
      omp_set_num_threads(1)
      
      num_sample <- 50
      pars_sample <- pars_prev_updated[, i, drop = F] %*% matrix(1, nrow = 1, ncol = num_sample) + t(rmvn(num_sample, mu = rep(0, ndim + 3), Sigma = diag(1, nrow = ndim + 3, ncol = ndim + 3)))
      loglik_res <- rep(NA, num_sample)
      
      for (k in 1:num_sample) {
        lpost <- function(p, frac) {
          minibatch <- sample(1:nrow(X_train), round(nrow(X_train) * frac))
          return(GPSDM(pars = p, distMat = distMat, basisX = X_basis, basisP = P_basis, basisY = Y_basis, newX = X_train[minibatch, , drop = F], newP = P_train[minibatch, , drop = F], newY = Y_train[minibatch, , drop = F], priors_tr = list(E = pars_prev[, i, drop = F], V = pars_var_prev[[i]]), mode = c("optimize")))
        }
        
        res <- lpost(pars_sample[, k, drop = F], frac = 1)
        loglik_res[k] <- -res$neglpost
        print(k)
      }
      loglik_res <- loglik_res - mean(loglik_res)
      loglik_res[loglik_res < -50] <- -50
      loglik_res[loglik_res > 50] <- 50
      lik_res <- exp(loglik_res)
      pars_var_prev_updated_i <- cov.wt(t(pars_sample),
                                        wt = as.numeric(exp(loglik_res)),
                                        cor = FALSE, center = t(pars_prev_updated[, i, drop = F])
      )$cov
      pars_var_prev_updated_i
    }
}

########################
# Add new sets of hyperparameters
########################
particles <- matrix(NA, nrow = ndim + 4, ncol = num_part)
for (i in 1:ndim) {
  particles[i, ] <- rtrunc(n = num_part, spec = "norm", a = phimin, b = phimax, mean = priors$E_phi[i], sd = sqrt(priors$V_phi[i]))
  particles[i, ] <- -log((phimax - phimin) / (particles[i, ] - phimin) - 1) # transform
}

particles[ndim + 1, ] <- rtrunc(n = num_part, spec = "norm", a = vemin, b = vemax, mean = priors$E_ve, sd = sqrt(priors$V_ve))
particles[ndim + 1, ] <- -log((vemax - vemin) / (particles[ndim + 1, ] - vemin) - 1) # transform

particles[ndim + 2, ] <- rtrunc(n = num_part, spec = "norm", a = taumin, b = taumax, mean = priors$E_tau, sd = sqrt(priors$V_tau))
particles[ndim + 2, ] <- -log((taumax - taumin) / (particles[ndim + 2, ] - taumin) - 1) # transform

particles[ndim + 3, ] <- rtrunc(n = num_part, spec = "norm", a = rhomin, b = rhomax, mean = priors$E_rho, sd = sqrt(priors$V_rho))
particles[ndim + 3, ] <- -log((rhomax - rhomin) / (particles[ndim + 3, ] - rhomin) - 1) # transform

particles[ndim + 4, ] <- rtrunc(n = num_part, spec = "norm", a = gammamin, b = gammamax, mean = priors$E_gamma, sd = sqrt(priors$V_gamma))
particles[ndim + 4, ] <- -log((gammamax - gammamin) / (particles[ndim + 4, ] - gammamin) - 1) # transform

pars_add <- particles
pars_id_add <- ((epoch - 1) * num_part + 1):(epoch * num_part)

Rprop_res_add <-
  foreach(
    i = 1:20,
    # .export = ls(globalenv()),
    .export = c(
      "phimax", "phimin", "vemax", "vemin", "taumax", "taumin", "gammamax", "gammamin",
      "pars_add", "distMat", "priors",
      "X_train", "P_train", "D_train","Y_train",
      "X_basis", "P_basis", "D_basis","Y_basis",
      "GPSDM", "fmingrad_Rprop"
    ),
    .packages = c("Matrix", "LaplacesDemon", "RhpcBLASctl")
  ) %dopar% {
    blas_set_num_threads(1)
    omp_set_num_threads(1)
    
    # subset<-sample(1:nrow(X_new), nrow(X_new))
    # X_new<-X_new[subset,,drop=F]
    # P_new<-P_new[subset,,drop=F]
    # Y_new<-Y_new[subset,,drop=F]
    
    lpost <- function(p, sample_n, seed) {
      set.seed(seed)
      minibatch <- sample(1:nrow(X_train), sample_n)
      return(GPSDM(pars = p, distMat = distMat, basisX = X_basis, basisP = P_basis, basisD=D_basis,basisY = Y_basis, newX = X_train[minibatch, , drop = F], newP = P_train[minibatch, , drop = F],newD = D_train[minibatch, , drop = F], newY = Y_train[minibatch, , drop = F], priors = priors, mode = c("optimize")))
    }
    print("start")
    
    fmingrad_Rprop(lpost, pars_add[, i, drop = F] + matrix(rnorm(ndim + 4, 0, 0.1)), sample_n = min(nrow(X_train),100), maxcount = 200)
  }

pars_add_updated <- matrix(NA, nrow = ndim + 4, ncol = num_part)
log_p_add_updated <- rep(NA, num_part)
for (i in 1:num_part) {
  pars_add_updated[, i] <- Rprop_res_add[[i]]$xopt
  log_p_add_updated[i] <- -Rprop_res_add[[i]]$fopt
}

pars_var_add_updated <-
  foreach(
    i = 1:num_part,
    .export = c( 
      "phimax", "phimin", "vemax", "vemin", "taumax", "taumin", "gammamax", "gammamin",
      "pars_add_updated", "distMat", "priors",
      "X_train", "P_train", "Y_train",
      "X_basis", "P_basis", "Y_basis",
      "GPSDM", "fmingrad_Rprop"
    ),
    .packages = c("Matrix", "LaplacesDemon", "RhpcBLASctl")
  ) %dopar% {
    blas_set_num_threads(1)
    omp_set_num_threads(1)
    
    num_sample <- 50
    pars_sample <- pars_add_updated[, i, drop = F] %*% matrix(1, nrow = 1, ncol = num_sample) + t(rmvn(num_sample, mu = rep(0, ndim + 4), Sigma = diag(1, nrow = ndim + 4, ncol = ndim + 4)))
    loglik_res <- rep(NA, num_sample)
    
    for (k in 1:num_sample) {
      lpost <- function(p, sample_n) {
        minibatch <- sample(1:nrow(X_train), sample_n)
        return(GPSDM(pars = p, distMat = distMat, basisX = X_basis, basisP = P_basis,basisD=D_basis, basisY = Y_basis, newX = X_train[minibatch, , drop = F], newP = P_train[minibatch, , drop = F],newD = D_train[minibatch, , drop = F], newY = Y_train[minibatch, , drop = F], priors = priors, mode = c("optimize")))
      }
      
      res <- lpost(pars_sample[, k, drop = F], sample_n = min(nrow(X_train),100))
      loglik_res[k] <- -res$neglpost
      print(k)
    }
    loglik_res <- loglik_res - mean(loglik_res)
    loglik_res[loglik_res < -50] <- -50
    loglik_res[loglik_res > 50] <- 50
    lik_res <- exp(loglik_res)
    pars_var_add_updated_i <- cov.wt(t(pars_sample),
                                     wt = as.numeric(exp(loglik_res)),
                                     cor = FALSE, center = t(pars_add_updated[, i, drop = F])
    )$cov
    pars_var_add_updated_i
  }

pars_new <- cbind(pars_prev_updated, pars_add_updated)
pars_var_new <- c(pars_var_prev_updated, pars_var_add_updated)
log_p_new <- c(log_p_prev_updated, log_p_add_updated)
pars_id_new <- c(pars_id_prev, pars_id_add)

########################
# Choose sets of hyperparameters
########################
# valid_id_subset<-vector(mode="list")
# for (n in 1:10) {
#   valid_id_subset[[n]]<-sample(1:nrow(X_valid),100)
# }

# log_p_valid_new <-
#   foreach(
#     i = 1:ncol(pars_new),
#     # .export = ls(globalenv()),
#     .export = c(
#       "phimax", "phimin", "vemax", "vemin", "taumax", "taumin", "gammamax", "gammamin",
#       "pars_new", "distMat", "priors",
#       "X_valid", "P_valid", "Y_valid",
#       "X_basis", "P_basis", "Y_basis",
#       "GPSDM", "fmingrad_Rprop"
#     ),
#     .combine = cbind
#   ) %dopar% {
#     blas_set_num_threads(1)
#     omp_set_num_threads(1)
#     
#     
#     num_sample <- 50
#     loglik_res <- rep(NA, num_sample)
#     
#     for (k in 1:num_sample) {
#       minibatch <- sample(1:nrow(X_valid), min(nrow(X_valid),100))
#       res <- GPSDM(pars = pars_new[, i, drop = F], distMat = distMat, basisX = X_basis, basisP = P_basis, basisY = Y_basis, newX = X_valid[minibatch, , drop = F], newP = P_valid[minibatch, , drop = F], newY = Y_valid[minibatch, , drop = F],priors=priors, mode = c("optimize"))
#       loglik_res[k] <- -res$neglpost
#       print(k)
#     }
#     
#     median(loglik_res)
#   }

sort_id <- rev(order(log_p_new))[1:num_part]
pars <- pars_new[, sort_id, drop = F]
pars_var <- pars_var_new[sort_id]
log_p <- log_p_new[sort_id]
pars_id <- pars_id_new[sort_id]
# log_p_valid <- log_p_valid_new[sort_id]
colnames(pars) <- pars_id

write_csv(as.data.frame(pars), paste0(path_results, "/", "pars", ".csv"))
for (i in 1:num_part) {
  write_csv(as.data.frame(pars_var[[i]]), paste0(path_results, "/", "pars_var_", i, ".csv"))
}
write_csv(as.data.frame(log_p), paste0(path_results, "/", "log_p_train", ".csv"))
# write_csv(as.data.frame(log_p_valid), paste0(path_results, "/", "log_p_valid", ".csv"))



###
path_results<-paste0(path, "/results_out")
# X_basis<-as.matrix(read_csv(paste0(path_results, "/", "X_basis", ".csv")))
# Y_basis<-as.matrix(read_csv(paste0(path_results, "/", "Y_basis", ".csv")))
# P_basis<-as.matrix(read_csv(paste0(path_results, "/", "P_basis", ".csv")))
# D_basis<-as.matrix(read_csv(paste0(path_results, "/", "D_basis", ".csv")))
pars <- as.matrix(read_csv(paste0(path_results, "/", "pars", ".csv")))
log_p <- as.matrix(read_csv(paste0(path_results, "/", "log_p_train", ".csv")))

# fill gaps in gcc
res<-PrepareEmbedding(x,start=max(unlist(lags))+1,end=length(date_list), focalsites =1:nrow(coord_df), lags=lags, neighbors=neighbors,vars=vars, distMat = distMat)

X_all<-res$X
Y_all<-res$Y
D_all<-res$D
P_all<-res$P

missing_id<-which(rowSums(is.na(cbind(Y_all,X_all)))!=0)
if (length(missing_id)>0) {
  X_all<-X_all[-missing_id,,drop=F]
  Y_all<-Y_all[-missing_id,,drop=F]
  P_all<-P_all[-missing_id,,drop=F]
  D_all<-D_all[-missing_id,,drop=F]
}

today_id<-which (date_list==today)
for(site in 1:nrow(coord_df)) {
  all_filled=F
  start_search<-max(unlist(lags))+1
  while(!all_filled) {
    na_id<-which(is.na(x[site,,1])) 
    # date_list[forecast_start]
    if (length(na_id[which(na_id>=start_search&na_id<today_id)])>0) {
      forecast_start<-min(na_id[which(na_id>=start_search&na_id<today_id)])-1
      steps=length(date_list)-forecast_start
      D_forecast <- date_list[(forecast_start + 1):(forecast_start + steps)] %>% 
        as.character() %>% 
        as.matrix()
      
      filled=F
      t=1
      while(!filled) {
        date<-as.Date(D_forecast[t,1])
        doy<-as.integer(format(date, "%j"))
        if(doy==366) {doy<-365}
        # moy<-as.integer(format(date, "%m"))
        # 
        # X_basis<-X_basis_list[[moy]]
        # Y_basis<-Y_basis_list[[moy]]
        # P_basis<-P_basis_list[[moy]]
        # # D_basis<-D_basis_list[[moy]]
        # pars<-pars_list[[moy]]
        # log_p<-log_p_list[[moy]]
        
        select_date<-data.frame(date_all=as.Date(D_all),
                                doy_all=date2doy(D_all)) %>% 
          mutate(a=abs(doy_all-doy),b=365-abs(doy_all-doy)) %>% 
          rowwise() %>% 
          mutate(diff=min(a,b))
        window_id<-which(select_date$diff<=30
                         # &select_date$date_all<date
                         )
        X_window <- X_all[window_id, , drop = F]
        Y_window <- Y_all[window_id, , drop = F]
        P_window <- P_all[window_id, , drop = F]
        D_window <- D_all[window_id, , drop = F]
        
        basisnumber<-min(basisnumber, nrow(X_window))
        if (basisnumber==nrow(X_all)) {
          X_basis <- X_window
          Y_basis <- Y_window
          P_basis <- P_window
          D_basis <- D_window
        } else {
          set.seed(42)
          cluster <- kmeans(cbind(P_window,X_window), basisnumber )
          cluster_id_sort <- seq(basisnumber )[order(cluster$size, decreasing = T)]
          basis_id <- rep(NA, basisnumber)
          for (i in 1:basisnumber) {
            # basis_id<-which(cl$cluster==cluster_id_sort[i])[length(which(cl$cluster==cluster_id_sort[i]))]
            basis_id[i] <- sample(which(cluster$cluster == cluster_id_sort[i]), 1)
          }
          X_basis <- X_window[basis_id, , drop = F]
          Y_basis <- Y_window[basis_id, , drop = F]
          P_basis <- P_window[basis_id, , drop = F]
          D_basis <- D_window[basis_id, , drop = F]
        }
        
        res<-PrepareEmbedding(x,start=forecast_start+t,end=forecast_start+t, focalsites = site, lags=lags, neighbors=neighbors,vars=vars, distMat = distMat)
        newX<-res$X
        newP<-res$P
        newD<-res$D
        res<-PrepareEmbedding(Sigma,start=forecast_start+t,end=forecast_start+t, focalsites = site, lags=lags, neighbors=neighbors,vars=vars, distMat = distMat)
        newSigma<-res$X
        
        if (sum(is.na(newX))>0) {
          x[site,(forecast_start+t),1]<-NA
          Sigma[site,(forecast_start+t),1]<-NA
        } else {
          Y_forecast_all<-
            foreach (i=1:num_part,
                     .export = c( # "num_part","h","lags","neighbors","vars","ndim",
                       "phimax", "phimin", "vemax", "vemin", "taumax", "taumin", "gammamax", "gammamin",
                       "distMat",
                       "GPSDM"
                     ),
                     .packages = c("tidyverse","RhpcBLASctl","MASS"),
                     .combine = rbind
            ) %dopar% {
              blas_set_num_threads(1)
              omp_set_num_threads(1)
              
              particle_pick <- pars[, i, drop = F]
              res <- GPSDM(pars = particle_pick, distMat = distMat, basisX = X_basis, basisP = P_basis, basisD = D_basis, basisY = Y_basis, newX = newX, newSigma=newSigma,newP = newP,newD=newD, mode = c("predict"))
              #store prediction
              
              c(res$mt, res$Ct)
            }
          x[site,(forecast_start+t),1]<-weighted.mean(Y_forecast_all[,1], w=exp(log_p-mean(log_p)))
          Sigma[site,(forecast_start+t),1]<-weighted.mean(Y_forecast_all[,2], w=exp(log_p-mean(log_p)))+Hmisc::wtd.var(Y_forecast_all[,1],weights=exp(log_p-mean(log_p)), normwt = T)
        }
        
        print(paste0(site,",",t))
        t<-t+1
        if(forecast_start+t==today_id) {
          all_filled<-T
          filled<-T
        } else {
          filled<-!is.na(x[site,(forecast_start+t),1])
          all_filled<-sum(is.na(x[site,(forecast_start+t):(today_id-1),1]))==0 
          start_search<-forecast_start+t
        }
      }
    } else {
      all_filled=T
    }
  }
}

#####
P_fit<-P_forecast<-as.matrix(1:nrow(coord_df))

forecast_start<-length(date_list)-365*2#max(unlist(lags))+1
date_list[forecast_start]

steps=length(date_list)-forecast_start
D_forecast <- date_list[(forecast_start + 1):(forecast_start + steps)] %>% 
  as.character() %>% 
  as.matrix()

xnew <- x[P_forecast,,,drop=F]
# xnew[, (forecast_start + 1):length(date_list), 1] <- NA
# xnew[,,3]<-xnew[,,3]+0.1 #warming
Sigmanew <- Sigma[P_forecast,,,drop=F]


Y_forecast<-Var_forecast<-matrix(NA, nrow=nrow(P_forecast), ncol = steps )      
for (t in 1:steps) {
  date<-as.Date(D_forecast[t,1])
  doy<-as.integer(format(date, "%j"))
  if(doy==366) {doy<-365}
  # moy<-as.integer(format(date, "%m"))
  # 
  # X_basis<-X_basis_list[[moy]]
  # Y_basis<-Y_basis_list[[moy]]
  # P_basis<-P_basis_list[[moy]]
  # # D_basis<-D_basis_list[[moy]]
  # pars<-pars_list[[moy]]
  # log_p<-log_p_list[[moy]]
  
  select_date<-data.frame(date_all=as.Date(D_all),
                          doy_all=date2doy(D_all)) %>% 
    mutate(a=abs(doy_all-doy),b=365-abs(doy_all-doy)) %>% 
    rowwise() %>% 
    mutate(diff=min(a,b))
  window_id<-which(select_date$diff<=30&select_date$date_all<date)
  X_window <- X_all[window_id, , drop = F]
  Y_window <- Y_all[window_id, , drop = F]
  P_window <- P_all[window_id, , drop = F]
  D_window <- D_all[window_id, , drop = F]
  
  basisnumber<-min(basisnumber, nrow(X_window))
  if (basisnumber==nrow(X_all)) {
    X_basis <- X_window
    Y_basis <- Y_window
    P_basis <- P_window
    D_basis <- D_window
  } else {
    set.seed(42)
    cluster <- kmeans(cbind(P_window,X_window), basisnumber )
    cluster_id_sort <- seq(basisnumber )[order(cluster$size, decreasing = T)]
    basis_id <- rep(NA, basisnumber)
    for (i in 1:basisnumber) {
      # basis_id<-which(cl$cluster==cluster_id_sort[i])[length(which(cl$cluster==cluster_id_sort[i]))]
      basis_id[i] <- sample(which(cluster$cluster == cluster_id_sort[i]), 1)
    }
    X_basis <- X_window[basis_id, , drop = F]
    Y_basis <- Y_window[basis_id, , drop = F]
    P_basis <- P_window[basis_id, , drop = F]
    D_basis <- D_window[basis_id, , drop = F]
  }
  
  res<-PrepareEmbedding(xnew,start=forecast_start+t,end=forecast_start+t, focalsites = P_forecast, lags=lags, neighbors=neighbors,vars=vars, distMat = distMat)
  newX<-res$X
  newP<-res$P
  newD<-res$D
  
  missing_id<-which(rowSums(is.na(newX))!=0)
  if (length(missing_id)>0) {
    newX<-newX[-missing_id,,drop=F]
    newP<-newP[-missing_id,,drop=F]
    newD<-newD[-missing_id,,drop=F]
  }
  res<-PrepareEmbedding(Sigmanew,start=forecast_start+t,end=forecast_start+t, focalsites = P_forecast, lags=lags, neighbors=neighbors,vars=vars, distMat = distMat)
  newSigma<-res$X
  if (length(missing_id)>0) {
    newSigma<-newSigma[-missing_id,,drop=F]
  }
  
  Y_forecast_all<-
    foreach (i=1:num_part,
             .export = c( # "num_part","h","lags","neighbors","vars","ndim",
               "phimax", "phimin", "vemax", "vemin", "taumax", "taumin", "gammamax", "gammamin",
               "distMat",
               "xnew","Sigmanew",
               "GPSDM","res"
             ),
             .packages = c("tidyverse","RhpcBLASctl")
    ) %dopar% {
      blas_set_num_threads(1)
      omp_set_num_threads(1)
      particle_pick <- pars[, i, drop = F]
      res <- GPSDM(pars = particle_pick, distMat = distMat, basisX = X_basis, basisP = P_basis,basisD = D_basis, basisY = Y_basis, newX = newX, newP = newP,newD = newD, newSigma=newSigma, mode = c("predict"))
      cbind(res$mt, res$Ct)
    }
  for(i in 1:nrow(newP)) {
    means<-variances<-c()
    for (j in 1:num_part) {
      means<-cbind(means,Y_forecast_all[[j]][i,1])
      variances<-cbind(variances,Y_forecast_all[[j]][i,2])
    }
    means<-matrix(means)
    variances<-matrix(variances)
    weighted_mean<-weighted.mean(means,w=exp(log_p-mean(log_p)))
    weighted_variance<-weighted.mean(variances, w=exp(log_p-mean(log_p)))+Hmisc::wtd.var(means,weights=exp(log_p-mean(log_p)), normwt = T)
    #store prediction
    Y_forecast[newP[i,1],t]<-weighted_mean
    Var_forecast[newP[i,1],t]<-weighted_variance
    
    xnew[newP[i,1],(forecast_start+t),1]<-weighted_mean
    Sigmanew[newP[i,1],(forecast_start+t),1]<-weighted_variance
  }
  
  print(t)
}

forecast_df<-as.data.frame(Y_forecast) %>% 
  `names<-`(D_forecast) %>% 
  mutate(site=P_forecast) %>% 
  gather(key="date", value="value",-site) %>% 
  mutate(date=as.Date(date))
var_forecast_df<-as.data.frame(Var_forecast) %>% 
  `names<-`(D_forecast) %>% 
  mutate(site=P_forecast) %>% 
  gather(key="date", value="value",-site) %>% 
  mutate(date=as.Date(date))
dir.create(paste0(path,"analyses_out"))
write_csv(forecast_df,paste0(path,"analyses_out/forecast.csv"))
write_csv(var_forecast_df,paste0(path,"analyses_out/var_forecast.csv"))


#### visualize
obs_df<-as.data.frame(x[P_fit,,1, drop=F]) %>%
  as_tibble() %>% 
  mutate(site=P_fit) %>%
  gather(key="date", value="y", -site) %>% 
  mutate(date=as.Date(date)) #%>%
  # mutate(site=as.character(site)) %>% 
  # drop_na() %>% 
  # mutate(year=as.integer(format(date, "%Y")))
var_obs_df<-as.data.frame(Sigma[P_fit,,1, drop=F]) %>%
  as_tibble() %>% 
  mutate(site=P_fit) %>%
  gather(key="date", value="y", -site) %>% 
  mutate(date=as.Date(date)) #%>%
  # mutate(site=as.character(site)) %>% 
  # drop_na() %>% 
  # mutate(year=as.integer(format(date, "%Y")))
obs_df<-left_join(obs_df,var_obs_df, by=c("site", "date")) %>% 
  dplyr::rename(y=y.x, variance=y.y) %>% 
  mutate(lower=y-1.96*sqrt(variance),
         upper=y+1.96*sqrt(variance))
obs_df_ori<-obs_df %>% 
  left_join(df_upper_lower[[1]], by="site")%>% 
  mutate(y=y*range+lower.y,
         upper=upper.x*range+lower.y,
         lower=lower.x*range+lower.y) #%>% 
# mutate(y=y*(range_list[[1]][2]-range_list[[1]][1])+range_list[[1]][1])

# p<-ggplot()+
#   geom_line(data=obs_df_ori, aes(x=date, y=y), alpha=0.5)+
#   geom_ribbon(data=obs_df_ori, aes(x=date, ymin=lower, ymax=upper),fill="blue", alpha=0.5)+
#   ylab("GCC")+
#   facet_wrap(~site, ncol = 3)+
#   theme_classic()
# p

###
forecast_df<-left_join(forecast_df,var_forecast_df, by=c("site", "date")) %>% 
  dplyr::rename(value=value.x, variance=value.y) %>% 
  mutate(lower=value-1.96*sqrt(variance),
         upper=value+1.96*sqrt(variance))

forecast_df_ori<-forecast_df %>% 
  left_join(df_upper_lower[[1]], by="site",suffix = c("", ".scale"))%>% 
  mutate(value=value*range+lower.scale,
         upper=upper*range+lower.scale,
         lower=lower*range+lower.scale) #%>% 
# mutate(median=median*(range_list[[1]][2]-range_list[[1]][1])+range_list[[1]][1],
#        upper=upper*(range_list[[1]][2]-range_list[[1]][1])+range_list[[1]][1],
#        lower=lower*(range_list[[1]][2]-range_list[[1]][1])+range_list[[1]][1])

#####
site_view<-forecast_df_ori %>%
  # group_by(site) %>%
  # dplyr::summarize(view=sum(!is.na(value))) %>%
  # filter(view!=0) %>%
  dplyr::select(site) %>%
  unlist()
# site_view<-which(rowSums(!is.na(Y_forecast_all[[1]]))!=0)
site_label<-unique(ts_all$siteID) %>% unlist()
names(site_label)<-1:length(site_label)
p<-ggplot()+
  geom_line(data=obs_df_ori  %>% filter(site%in%site_view), aes(x=date, y=y))+
  geom_ribbon(data=obs_df_ori%>% filter(site%in%site_view), aes(x=date, ymax=upper, ymin=lower), alpha=0.25)+
  geom_line(data=forecast_df_ori %>% filter(site%in%site_view), aes(x=date, y=value), col="blue")+
  geom_ribbon(data=forecast_df_ori%>% filter(site%in%site_view), aes(x=date, ymax=upper, ymin=lower), fill="blue", alpha=0.25)+
  geom_vline(xintercept = date_list[forecast_start+1], col="blue", alpha=0.5)+
  geom_vline(xintercept = date_list[forecast_start+1-365], col="blue", alpha=0.5)+
  ylab("GCC_90")+
  facet_wrap(~site, ncol = 2,labeller = labeller(site=site_label))+
  theme_classic()
p


###
combined_df<-obs_df %>% 
  left_join(forecast_df, by=c("site", "date")) %>% 
  drop_na() %>% 
  filter(site %in% 1:4,
         date<date_list[length(date_list)-365])
combined_df_ori<-obs_df_ori %>% 
  dplyr::select(-lower, -upper,-range)%>% 
  left_join(forecast_df_ori, by=c("site", "date")) %>% 
  dplyr::select(-lower.scale, -upper.scale,-range) %>% 
  drop_na() %>% 
  filter(site %in% 1:4,
         date<date_list[length(date_list)-365])
compare_stats( obs_ori=combined_df_ori$y, pred_ori=combined_df_ori$value,obs=combined_df$y,pred=combined_df$value)

### out-of-sample 1
combined_df<-obs_df %>% 
  left_join(forecast_df, by=c("site", "date")) %>% 
  drop_na() %>% 
  filter(site %in% 1:4,
         date>=date_list[length(date_list)-365])
combined_df_ori<-obs_df_ori %>% 
  dplyr::select(-lower, -upper,-range)%>% 
  left_join(forecast_df_ori, by=c("site", "date")) %>% 
  dplyr::select(-lower.scale, -upper.scale,-range) %>% 
  drop_na() %>% 
  filter(site %in% 1:4,
         date>=date_list[length(date_list)-365])
compare_stats( obs_ori=combined_df_ori$y, pred_ori=combined_df_ori$value,obs=combined_df$y,pred=combined_df$value)

### out-of-sample 2
combined_df<-obs_df %>% 
  left_join(forecast_df, by=c("site", "date")) %>% 
  drop_na() %>% 
  filter(site %in% 5:8)
combined_df_ori<-obs_df_ori %>% 
  dplyr::select(-lower, -upper,-range)%>% 
  left_join(forecast_df_ori, by=c("site", "date")) %>% 
  dplyr::select(-lower.scale, -upper.scale,-range) %>% 
  drop_na() %>% 
  filter(site %in% 5:8)
compare_stats( obs_ori=combined_df_ori$y, pred_ori=combined_df_ori$value,obs=combined_df$y,pred=combined_df$value)


