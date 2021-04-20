P_fit<-P_forecast<-as.matrix(1:nrow(coord_df))

forecast_start<-which(date_list==today)-1#length(date_list)-365*2#max(unlist(lags))+1
date_list[forecast_start]

steps=35#length(date_list)-forecast_start
D_forecast <- date_list[(forecast_start + 1):(forecast_start + steps)] %>% 
  as.character() %>% 
  as.matrix()

xnew <- x[P_forecast,,,drop=F]
# xnew[, (forecast_start + 1):length(date_list), 1] <- NA
# xnew[,,3]<-xnew[,,3]+0.1 #warming
Sigmanew <- Sigma[P_forecast,,,drop=F]


Y_forecast<-Var_forecast<-matrix(NA, nrow=nrow(P_forecast), ncol = steps )      
for (t in 1:steps) {
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
  # window_id<-which(select_date$diff<=180&select_date$date_all<date)
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
dir.create(paste0(path,"analyses"))
write_csv(forecast_df,paste0(path,"analyses/forecast.csv"))
write_csv(var_forecast_df,paste0(path,"analyses/var_forecast.csv"))


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
  mutate(y=(y+0.5)*range+lower.y,
         upper=(upper.x+0.5)*range+lower.y,
         lower=(lower.x+0.5)*range+lower.y) #%>% 
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
  mutate(value=(value+0.5)*range+lower.scale,
         upper=(upper+0.5)*range+lower.scale,
         lower=(lower+0.5)*range+lower.scale) #%>% 
# mutate(median=median*(range_list[[1]][2]-range_list[[1]][1])+range_list[[1]][1],
#        upper=upper*(range_list[[1]][2]-range_list[[1]][1])+range_list[[1]][1],
#        lower=lower*(range_list[[1]][2]-range_list[[1]][1])+range_list[[1]][1])
