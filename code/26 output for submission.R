

######

df_submit<-forecast_df_ori %>% 
  mutate(sd=sqrt(variance)) %>% 
  dplyr::select(site, time=date, gcc_90=value, gcc_sd=sd) %>% 
  left_join(as.data.frame(site_list) %>% mutate(site=1:8) %>% rename(siteID=site_list), by="site") %>% 
  dplyr::select(time, siteID,
                mean=gcc_90,
                sd=gcc_sd) %>% 
  gather(key="statistic",value="gcc_90",-siteID,-time) %>% 
  mutate(obs_flag=2,
         forecast=1,
         data_assimilation=1) %>% 
  dplyr::select(time,
                siteID,
                obs_flag,
                forecast,
                data_assimilation,
                statistic,
                gcc_90) %>% 
  arrange(time,siteID)
df_submit

year<-format(min(df_submit$time), "%Y")
month<-format(min(df_submit$time), "%m")
day<-format(min(df_submit$time), "%d")
write_csv(df_submit, paste0(path,"phenology-",year,"-",month,"-",day,"-UCSC_P_EDM.csv"))



########
site_view<-forecast_df_ori %>%
  group_by(site) %>%
  dplyr::summarize(view=sum(!is.na(value))) %>%
  filter(view!=0) %>%
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
  geom_vline(xintercept = date_list[forecast_start+1]-years(1:4), col="blue", alpha=0.5)+
  ylab("GCC_90")+
  ylim(min(obs_df_ori$lower,forecast_df_ori$lower ,na.rm = T),
       max(obs_df_ori$upper,forecast_df_ori$upper ,na.rm = T))+
  facet_wrap(~site, ncol = 2,labeller = labeller(site=site_label))+
  theme_classic()

cairo_pdf(paste0(path,"phenology-",year,"-",month,"-",day,"-UCSC_P_EDM.pdf"), width = 16, height = 8)
print(p)
dev.off()

p<-ggplot()+
  geom_line(data=obs_df_ori  %>% filter(site%in%site_view) %>% filter(date>=date_list[forecast_start+1]), aes(x=date, y=y))+
  geom_ribbon(data=obs_df_ori%>% filter(site%in%site_view) %>% filter(date>=date_list[forecast_start+1]), aes(x=date, ymax=upper, ymin=lower), alpha=0.25)+
  geom_line(data=forecast_df_ori %>% filter(site%in%site_view) %>% filter(date>=date_list[forecast_start+1]), aes(x=date, y=value), col="blue")+
  geom_ribbon(data=forecast_df_ori%>% filter(site%in%site_view) %>% filter(date>=date_list[forecast_start+1]), aes(x=date, ymax=upper, ymin=lower), fill="blue", alpha=0.25)+
  ylab("GCC_90")+
  # ylim(-0.1,1.1)+
  facet_wrap(~site, ncol = 2,labeller = labeller(site=site_label))+
  theme_classic()

cairo_pdf(paste0(path,"phenology-",year,"-",month,"-",day,"-UCSC_P_EDM_fore_only.pdf"), width = 16, height = 8)
print(p)
dev.off()
########

## define variable names, units, etc
## in practice, this might be kept in a spreadsheet
attributes <- tibble::tribble(
  ~attributeName,     ~attributeDefinition,                          ~unit,                  ~formatString, ~numberType, ~definition,
  "time",              "[dimension]{time}",                          "year",                 "YYYY-MM-DD",  "numberType", NA,
  "siteID",             "[dimension]{ID of NEON site}",            "dimensionless",        NA,           "character",       NA,
  "obs_flag",          "[dimension]{observation error}",             "dimensionless",         NA,           "integer",    NA,
  "forecast",          "[flag]{whether time step assimilated data}", "dimensionless",         NA,           "integer",    NA,
  "data_assimilation", "[flag]{whether time step assimilated data}", "dimensionless",         NA,           "integer",    NA,
  "statistic",     "[dimension]{descriptive statistic}",        "dimensionless",         NA,           "character",    NA,
  "gcc_90",         "[variable]{Green Chromatic Coordinate}",      "dimensionless", NA,           "real",       NA
) 
factors1<-data.frame(attributeName="siteID",
                     code=site_list,
                     definition=site_list)
factors2<-data.frame(attributeName="statistic",
                     code=c("mean","sd"),
                     definition=c("mean","standard deviation"))
factors<-rbind.data.frame(factors1,factors2)
## note: EML uses a different unit standard than UDUNITS. For now use EML. EFI needs to provide a custom unitList.
attributes
factors
attrList <- set_attributes(attributes, 
                           factors,
                           col_classes = c("Date", "factor", "numeric","numeric", 
                                           "numeric","factor","numeric"))

## sets metadata about the file itself (name, file type, size, MD5, etc)
physical <- set_physical(paste0(path,"phenology-",year,"-",month,"-",day,"-UCSC_P_EDM.csv"),
                         recordDelimiter='\n')

## set metadata for the file as a whole
dataTable <- eml$dataTable(
  entityName = "forecast",  ## this is a standard name to allow us to distinguish this entity from 
  entityDescription = "Forecast of phenology using empirical dynamic modeling",
  physical = physical,
  attributeList = attrList)

me <- list(individualName = list(givenName = "Yiluan", 
                                 surName = "Song"),
           electronicMailAddress = "ysong67@ucsc.edu",
           id = "https://orcid.org/0000-0003-3660-3797")


taxa <- tibble::tribble(
  ~Genus,      ~Species)
coverage <- 
  set_coverage(begin = min(df_submit$time), 
               end = max(df_submit$time),
               sci_names = taxa,
               geographicDescription = "USA",
               west = min(coord_df$lon),
               east = max(coord_df$lon), 
               north = max(coord_df$lat),
               south = min(coord_df$lat)
  )

keywordSet <- list(
  list(
    keywordThesaurus = "EFI controlled vocabulary",
    keyword = list("forecast",
                   "phenology",
                   "timeseries")
  ))

dataset = eml$dataset(
  title = "Forecast of phenology using empirical dynamic modeling",
  creator = me,
  contact = list(references="https://orcid.org/0000-0003-3660-3797"),
  pubDate = today,
  intellectualRights = "http://www.lternet.edu/data/netpolicy.html.",
  abstract =  "Forecast of phenology using empirical dynamic modeling",
  dataTable = dataTable,
  keywordSet = keywordSet,
  coverage = coverage
)


additionalMetadata <- eml$additionalMetadata(
  metadata = list(
    forecast = list(
      ## Basic elements
      timestep = "1 day", ## should be udunits parsable; already in coverage -> temporalCoverage?
      forecast_horizon = "35 days",
      forecast_issue_time = today,
      forecast_iteration_id = forecast_iteration_id,
      forecast_project_id = forecast_project_id,
      metadata_standard_version = "0.3",
      model_description = list(
        forecast_model_id = forecast_model_id,
        name = "Gaussian Process empirical dynamic modeling",
        type = "machine learning",
        repository = "https://github.com/yiluansong/EFI"
      ),
      ## MODEL STRUCTURE & UNCERTAINTY CLASSES
      initial_conditions = list(
        # Possible values: absent, present, data_driven, propagates, assimilates
        status = "present",
        # Number of parameters / dimensionality
        complexity = 1  
      ),
      drivers = list(
        status = "present",
        complexity = length(var_list)-1
      ),
      parameters = list(
        status = "present",
        complexity = ndim*12  ## [r, K, alpha] x 2 spp
      ),
      random_effects = list(
        status = "absent"
      ),
      process_error = list(
        status = "propagates",
        propagation = list(
          type = "analytic" # ensemble vs analytic
        ),
        complexity = 1,   
        covariance = FALSE
      ),
      obs_error = list(
        status = "present",
        complexity = 1,   
        covariance = FALSE
      )
    ) # forecast
  ) # metadata
) # eml$additionalMetadata

my_eml <- eml$eml(dataset = dataset,
                  additionalMetadata = additionalMetadata,
                  packageId = forecast_iteration_id , 
                  system = "datetime"  ## system used to generate packageId
)

## check base EML
eml_validate(my_eml)

## check that the EML is also a valid EFI forecast
EFIstandards::forecast_validator(my_eml)

write_eml(my_eml, paste0(path,"phenology-",year,"-",month,"-",day,"-UCSC_P_EDM.eml"))

