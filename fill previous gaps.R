X_basis_list<-Y_basis_list<-P_basis_list<-D_basis_list<-pars_list<-log_p_list<-vector(mode="list", length=12)
for (moy_oi in 1:12) {
  path_results<-paste0(path, "/results/MOY", moy_oi)
  
  X_basis_list[[moy_oi]]<-as.matrix(read_csv(paste0(path_results, "/", "X_basis", ".csv")))
  Y_basis_list[[moy_oi]]<-as.matrix(read_csv(paste0(path_results, "/", "Y_basis", ".csv")))
  P_basis_list[[moy_oi]]<-as.matrix(read_csv(paste0(path_results, "/", "P_basis", ".csv")))
  D_basis_list[[moy_oi]]<-as.matrix(read_csv(paste0(path_results, "/", "D_basis", ".csv")))
  pars_list[[moy_oi]] <- as.matrix(read_csv(paste0(path_results, "/", "pars", ".csv")))
  log_p_list[[moy_oi]] <- as.matrix(read_csv(paste0(path_results, "/", "log_p_train", ".csv")))
}



# fill gaps in gcc
today_id<-which (date_list==today)
for(site in 1:nrow(coord_df)) {
  all_filled=F
  while(!all_filled) {
    na_id<-which(is.na(x[site,,1])) 
    # date_list[forecast_start]
    if (length(na_id[which(na_id>=1115&na_id<today_id)])>0) {
      forecast_start<-min(na_id[which(na_id>=1115&na_id<today_id)])-1
      steps=length(date_list)-forecast_start
      D_forecast <- date_list[(forecast_start + 1):(forecast_start + steps)] %>% 
        as.character() %>% 
        as.matrix()
      
      filled=F
      t=1
      while(!filled) {
        date<-as.Date(D_forecast[t,1])
        moy<-as.integer(format(date, "%m"))
        
        X_basis<-X_basis_list[[moy]]
        Y_basis<-Y_basis_list[[moy]]
        P_basis<-P_basis_list[[moy]]
        # D_basis<-D_basis_list[[moy]]
        pars<-pars_list[[moy]]
        log_p<-log_p_list[[moy]]
        
        res<-PrepareEmbedding(x,start=forecast_start+t,end=forecast_start+t, focalsites = site, lags=lags, neighbors=neighbors,vars=vars, distMat = distMat)
        newX<-res$X
        newP<-res$P
        res<-PrepareEmbedding(Sigma,start=forecast_start+t,end=forecast_start+t, focalsites = site, lags=lags, neighbors=neighbors,vars=vars, distMat = distMat)
        newSigma<-res$X
        
        if (sum(is.na(newX))>0) {
          x[site,(forecast_start+t),1]<-NA
          Sigma[,(forecast_start+t),1]<-NA
        } else {
          Y_forecast_all<-
            foreach (i=1:num_part,
                     .export = c( # "num_part","h","lags","neighbors","vars","ndim",
                       "phimax", "phimin", "vemax", "vemin", "taumax", "taumin", "gammamax", "gammamin",
                       "distMat",
                       "GPSDM"
                     ),
                     .packages = c("tidyverse","RhpcBLASctl"),
                     .combine = rbind
            ) %dopar% {
              blas_set_num_threads(1)
              omp_set_num_threads(1)
              
              particle_pick <- pars[, i, drop = F]
              res <- GPSDM(pars = particle_pick, distMat = distMat, basisX = X_basis, basisP = P_basis, basisY = Y_basis, newX = newX, newSigma=newSigma,newP = newP, mode = c("predict"))
              #store prediction
              
              c(res$mt, res$Ct)
            }
          x[site,(forecast_start+t),1]<-weighted.mean(Y_forecast_all[,1], w=exp(log_p-mean(log_p)))
          Sigma[,(forecast_start+t),1]<-weighted.mean(Y_forecast_all[,2], w=exp(log_p-mean(log_p)))+Hmisc::wtd.var(Y_forecast_all[,1],weights=exp(log_p-mean(log_p)))
        }
        
        print(paste0(site,",",t))
        t<-t+1
        if(forecast_start+t==today_id) {
          all_filled<-T
          filled<-T
        } else {
          filled<-!is.na(x[site,(forecast_start+t),1])
          all_filled<-sum(is.na(x[site,(forecast_start+t):(today_id-1),1]))==0 
        }
      }
    } else {
      all_filled=T
    }
  }
  
}


