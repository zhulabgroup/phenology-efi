today<-as.Date("2021-02-02")

source("EFI utils.R")
source("EFI settings.R")

update<-T
source("get NEON phenology data.R")
update<-F
source("get NEON weather data.R")
update<-T
source("get NOAA weather data.R")

source("preprocess data.R")
source("prepare embeddings.R")
source("train GP model.R")
source("fill previous gaps.R")
source("check forecast.R")

cairo_pdf(paste0(path,"phenology-",year,"-",month,"-",day,"-UCSC_P_EDM.pdf"), width = 16, height = 8)
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

