library(tidyverse)
files<-list.files("G:\\My Drive\\GitHub\\EFI\\archive",pattern="-UCSC_P_EDM.csv", recursive=T, full.names = T) %>% sort()
date_list<-seq(as.Date("2021-02-01"),as.Date("2021-06-14"), by="day")
df_list<-vector(mode="list")
for (i in 1:length(files)) {
  df_day<-read_csv(files[i]) %>% 
    mutate(date_start=date_list[i]) %>% 
    mutate(date_ahead=(time-date_start) %>% as.integer)
    df_list[[i]]<- df_day
}
df<-bind_rows(df_list)%>% 
  dplyr::select(-obs_flag, -forecast, -data_assimilation) %>% 
  filter(date_ahead==30) %>% 
  dplyr::select(-rcc_90) %>%
  spread(key=statistic, value="gcc_90")

latest<-read_csv("G:\\My Drive\\GitHub\\EFI\\archive\\2021-07-05\\phenology-targets.csv")

latest<-latest %>% 
  dplyr::select(-rcc_90, -rcc_sd) %>% 
  rename(mean=gcc_90,sd=gcc_sd) 
  

mean <-latest %>% 
  mutate(doy=format(time, "%j") %>% as.integer()) %>% 
  group_by(siteID, doy) %>% 
  summarise (mean=mean(mean, na.rm=T))

climatology<-latest %>% 
  filter (time>=as.Date("2021-02-01")) %>% 
  mutate(doy=format(time, "%j") %>% as.integer()) %>%
  dplyr::select(siteID, time,doy) %>% 
  full_join(mean, by=c("siteID", "doy"))%>% 
  dplyr::select(-doy)

df_compare<-latest %>%
  filter (time>=as.Date("2021-02-01")) %>% 
  dplyr::select(-sd) %>% 
  rename(observed=mean) %>% 
  left_join(climatology %>% 
              rename(climatology="mean"),
            by=c("siteID", "time")) %>% 
  left_join(df %>% 
              dplyr::select(-date_start, -sd, -date_ahead) %>% 
              rename(forecasted=mean),
            by=c("siteID", "time")) %>% 
  mutate(forecast_error=abs(forecasted-observed)) %>%
  mutate(climatology_error=abs(climatology-observed)) %>%
  gather(key="group", value="value",-siteID, -time)
ggplot(df_compare %>% filter(!group %in% c("forecast_error", "climatology_error")))+
  geom_line(aes(x=time, y=value, col=group))+
  theme_classic()+
  facet_wrap(.~siteID)

ggplot(df_compare %>% filter(group %in% c("forecast_error", "climatology_error")))+
  geom_line(aes(x=time, y=value, col=group))+
  theme_classic()+
  facet_wrap(.~siteID)

df_compare %>%
  filter(group %in% c("forecast_error", "climatology_error")) %>%
  spread(key="group", value="value") %>% 
  drop_na() %>% 
  group_by(siteID) %>% 
  summarize (climatology_error=mean(climatology_error, na.rm=T),
             forecast_error=mean(forecast_error, na.rm=T))




####
df<-bind_rows(df_list)%>% 
  dplyr::select(-obs_flag, -forecast, -data_assimilation) %>% 
  filter(date_ahead==30) %>% 
  dplyr::select(-gcc_90) %>%
  spread(key=statistic, value="rcc_90")

latest<-read_csv("G:\\My Drive\\GitHub\\EFI\\archive\\2021-07-05\\phenology-targets.csv")

latest<-latest %>% 
  dplyr::select(-gcc_90, -gcc_sd) %>% 
  rename(mean=rcc_90,sd=rcc_sd) 


mean <-latest %>% 
  mutate(doy=format(time, "%j") %>% as.integer()) %>% 
  group_by(siteID, doy) %>% 
  summarise (mean=mean(mean, na.rm=T))

climatology<-latest %>% 
  filter (time>=as.Date("2021-02-01")) %>% 
  mutate(doy=format(time, "%j") %>% as.integer()) %>%
  dplyr::select(siteID, time,doy) %>% 
  full_join(mean, by=c("siteID", "doy"))%>% 
  dplyr::select(-doy)

df_compare<-latest %>%
  filter (time>=as.Date("2021-02-01")) %>% 
  dplyr::select(-sd) %>% 
  rename(observed=mean) %>% 
  left_join(climatology %>% 
              rename(climatology="mean"),
            by=c("siteID", "time")) %>% 
  left_join(df %>% 
              dplyr::select(-date_start, -sd, -date_ahead) %>% 
              rename(forecasted=mean),
            by=c("siteID", "time")) %>% 
  mutate(forecast_error=abs(forecasted-observed)) %>%
  mutate(climatology_error=abs(climatology-observed)) %>%
  gather(key="group", value="value",-siteID, -time)
ggplot(df_compare %>% filter(!group %in% c("forecast_error", "climatology_error")))+
  geom_line(aes(x=time, y=value, col=group))+
  theme_classic()+
  facet_wrap(.~siteID)

ggplot(df_compare %>% filter(group %in% c("forecast_error", "climatology_error")))+
  geom_line(aes(x=time, y=value, col=group))+
  theme_classic()+
  facet_wrap(.~siteID)

df_compare %>%
  filter(group %in% c("forecast_error", "climatology_error")) %>%
  spread(key="group", value="value") %>% 
  drop_na() %>% 
  group_by(siteID) %>% 
  summarize (climatology_error=mean(climatology_error, na.rm=T),
             forecast_error=mean(forecast_error, na.rm=T))
