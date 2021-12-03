
# fill gaps in gcc
# res<-PrepareEmbedding(x,start=max(unlist(lags))+1,end=length(date_list), focalsites =coord_df_site$id, lags=lags, neighbors=neighbors,vars=vars, distMat = distMat)
# 
# X_all<-res$X
# Y_all<-res$Y
# D_all<-res$D
# P_all<-res$P
# 
# missing_id<-which(rowSums(is.na(cbind(Y_all,X_all)))!=0)
# if (length(missing_id)>0) {
#   X_all<-X_all[-missing_id,,drop=F]
#   Y_all<-Y_all[-missing_id,,drop=F]
#   P_all<-P_all[-missing_id,,drop=F]
#   D_all<-D_all[-missing_id,,drop=F]
# }

today_id<-which (date_list==today)
for(site in 1:nrow(coord_df)) {
  path_results<-paste0(path, focal_var,"/results/", site)
  X_basis<-as.matrix(read_csv(paste0(path_results,"/X_basis.csv")))
  Y_basis<-as.matrix(read_csv(paste0(path_results,"/Y_basis.csv")))
  P_basis<-as.matrix(read_csv(paste0(path_results,"/P_basis.csv")))
  D_basis<-as.matrix(read_csv(paste0(path_results,"/D_basis.csv")))
  pars <- as.matrix(read_csv(paste0(path_results,"/pars.csv")))
  log_p <- as.matrix(read_csv(paste0(path_results,"/log_p_train.csv")))
  
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