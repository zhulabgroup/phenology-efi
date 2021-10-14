path_results<-paste0(path, "/results")
X_basis<-as.matrix(read_csv(paste0(path_results, "/", "X_basis", ".csv")))
Y_basis<-as.matrix(read_csv(paste0(path_results, "/", "Y_basis", ".csv")))
P_basis<-as.matrix(read_csv(paste0(path_results, "/", "P_basis", ".csv")))
D_basis<-as.matrix(read_csv(paste0(path_results, "/", "D_basis", ".csv")))
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
        # date<-as.Date(D_forecast[t,1])
        # doy<-as.integer(format(date, "%j"))
        # if(doy==366) {doy<-365}
        # # moy<-as.integer(format(date, "%m"))
        # # 
        # # X_basis<-X_basis_list[[moy]]
        # # Y_basis<-Y_basis_list[[moy]]
        # # P_basis<-P_basis_list[[moy]]
        # # # D_basis<-D_basis_list[[moy]]
        # # pars<-pars_list[[moy]]
        # # log_p<-log_p_list[[moy]]
        # 
        # select_date<-data.frame(date_all=as.Date(D_all),
        #                         doy_all=date2doy(D_all)) %>% 
        #   mutate(a=abs(doy_all-doy),b=365-abs(doy_all-doy)) %>% 
        #   rowwise() %>% 
        #   mutate(diff=min(a,b))
        # window_id<-which(select_date$diff<=30
        #                  # &select_date$date_all<date
        # )
        # X_window <- X_all[window_id, , drop = F]
        # Y_window <- Y_all[window_id, , drop = F]
        # P_window <- P_all[window_id, , drop = F]
        # D_window <- D_all[window_id, , drop = F]
        # 
        # basisnumber<-min(basisnumber, nrow(X_window))
        # if (basisnumber==nrow(X_all)) {
        #   X_basis <- X_window
        #   Y_basis <- Y_window
        #   P_basis <- P_window
        #   D_basis <- D_window
        # } else {
        #   set.seed(42)
        #   cluster <- kmeans(cbind(P_window,X_window), basisnumber )
        #   cluster_id_sort <- seq(basisnumber )[order(cluster$size, decreasing = T)]
        #   basis_id <- rep(NA, basisnumber)
        #   for (i in 1:basisnumber) {
        #     # basis_id<-which(cl$cluster==cluster_id_sort[i])[length(which(cl$cluster==cluster_id_sort[i]))]
        #     basis_id[i] <- sample(which(cluster$cluster == cluster_id_sort[i]), 1)
        #   }
        #   X_basis <- X_window[basis_id, , drop = F]
        #   Y_basis <- Y_window[basis_id, , drop = F]
        #   P_basis <- P_window[basis_id, , drop = F]
        #   D_basis <- D_window[basis_id, , drop = F]
        # }
        
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