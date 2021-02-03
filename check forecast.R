P_fit<-P_forecast<-as.matrix(1:nrow(coord_df))

forecast_start<-which(date_list==max(ts_all$time))
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
  moy<-as.integer(format(date, "%m"))
  
  X_basis<-X_basis_list[[moy]]
  Y_basis<-Y_basis_list[[moy]]
  P_basis<-P_basis_list[[moy]]
  # D_basis<-D_basis_list[[moy]]
  pars<-pars_list[[moy]]
  log_p<-log_p_list[[moy]]
  
  res<-PrepareEmbedding(xnew,start=forecast_start+t,end=forecast_start+t, focalsites = P_forecast, lags=lags, neighbors=neighbors,vars=vars, distMat = distMat)
  newX<-res$X
  newP<-res$P
  res<-PrepareEmbedding(Sigmanew,start=forecast_start+t,end=forecast_start+t, focalsites = P_forecast, lags=lags, neighbors=neighbors,vars=vars, distMat = distMat)
  newSigma<-res$X
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
      res <- GPSDM(pars = particle_pick, distMat = distMat, basisX = X_basis, basisP = P_basis, basisY = Y_basis, newX = newX, newP = newP, newSigma=newSigma, mode = c("predict"))
      cbind(res$mt, res$Ct)
    }
  for(i in 1:nrow(P_forecast)) {
    means<-variances<-c()
    for (j in 1:num_part) {
      means<-cbind(means,Y_forecast_all[[j]][i,1])
      variances<-cbind(variances,Y_forecast_all[[j]][i,2])
    }
    means<-matrix(means)
    variances<-matrix(variances)
    weighted_mean<-weighted.mean(means,w=exp(log_p-mean(log_p)))
    weighted_variance<-weighted.mean(variances, w=exp(log_p-mean(log_p)))+Hmisc::wtd.var(means,weights=exp(log_p-mean(log_p)))
    #store prediction
    Y_forecast[i,t]<-weighted_mean
    Var_forecast[,t]<-weighted_variance
    
    xnew[i,(forecast_start+t),1]<-weighted_mean
    Sigmanew[,(forecast_start+t),1]<-weighted_variance
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
  mutate(date=as.Date(date)) %>%
  # mutate(site=as.character(site)) %>% 
  # drop_na() %>% 
  mutate(year=as.integer(format(date, "%Y")))
var_obs_df<-as.data.frame(Sigma[P_fit,,1, drop=F]) %>%
  as_tibble() %>% 
  mutate(site=P_fit) %>%
  gather(key="date", value="y", -site) %>% 
  mutate(date=as.Date(date)) %>%
  # mutate(site=as.character(site)) %>% 
  # drop_na() %>% 
  mutate(year=as.integer(format(date, "%Y")))
obs_df<-left_join(obs_df,var_obs_df, by=c("site", "date","year")) %>% 
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


df_submit<-forecast_df_ori %>% 
  mutate(sd=sqrt(variance)) %>% 
  dplyr::select(site, time=date, gcc_90=value, gcc_sd=sd) %>% 
  left_join(as.data.frame(site_list) %>% mutate(site=1:8) %>% rename(siteID=site_list), by="site") %>% 
  dplyr::select(time, siteID,gcc_90, gcc_sd) %>% 
  arrange(siteID)
df_submit

year<-format(min(df_submit$time), "%Y")
month<-format(min(df_submit$time), "%m")
day<-format(min(df_submit$time), "%d")
write_csv(df_submit, paste0(path,"phenocam-",year,"-",month,"-",day,"-UCSC_P_EDM.csv"))

cairo_pdf(paste0(path,"phenocam-",year,"-",month,"-",day,"-UCSC_P_EDM.pdf"), width = 16, height = 8)
site_view<-forecast_df_ori %>% 
  group_by(site) %>% 
  dplyr::summarize(view=sum(!is.na(value))) %>% 
  filter(view!=0) %>% 
  dplyr::select(site) %>%
  unlist()
# site_view<-which(rowSums(!is.na(Y_forecast_all[[1]]))!=0)
site_label<-unique(ts_all$siteID) %>% unlist()
names(site_label)<-1:length(site_label)
ggplot()+
  geom_line(data=obs_df_ori  %>% filter(site%in%site_view), aes(x=date, y=y))+
  geom_ribbon(data=obs_df_ori%>% filter(site%in%site_view), aes(x=date, ymax=upper, ymin=lower), alpha=0.25)+
  geom_line(data=forecast_df_ori %>% filter(site%in%site_view), aes(x=date, y=value), col="blue")+
  geom_ribbon(data=forecast_df_ori%>% filter(site%in%site_view), aes(x=date, ymax=upper, ymin=lower), fill="blue", alpha=0.25)+
  geom_vline(xintercept = date_list[forecast_start+1], col="blue", alpha=0.5)+
  ylab("GCC")+
  facet_wrap(~site, ncol = 2,labeller = labeller(site=site_label))+
  theme_classic()
dev.off()
