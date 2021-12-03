today_list<-seq(as.Date("2021-02-02"),as.Date("2021-12-01"),by=1)
for (d in 1:length(today_list)) {
  todaydone<-F
  while (!todaydone) {
    tryCatch( {
      today<-today_list[d]
      path<-paste0("./archive/",today,"/")
      dir.create(path, recursive = T)
      
      forecast_project_id<-"UCSC_P_EDM" 
      forecast_model_id   <- "v1.8"
      forecast_iteration_id<- uuid::UUIDgenerate() 
      
      
      source("code/01 EFI utils.R")
      cl <- makeCluster(5, outfile = "")
      registerDoSNOW(cl)
      
      update<-T
      source("code/11 get NEON phenology data.R")
      update<-T
      source("code/12 get NEON weather data api.R")
      update<-T
      source("code/13 get NOAA weather data.R")
      
      forecast_df_ori_allvar<-vector(mode = "list", length=2)
      focal_var<-"gcc"
      source("code/02 EFI settings.R")
      source("code/21 preprocess data.R")
      source("code/22 prepare embeddings.R")
      source("code/23 train GP model.R")
      source("code/24 fill previous gaps.R")
      source("code/25 forecast.R")
      obs_df_ori_allvar[[1]]<-obs_df_ori
      forecast_df_ori_allvar[[1]]<-forecast_df_ori
      
      focal_var<-"rcc"
      source("code/02 EFI settings.R")
      source("code/21 preprocess data.R")
      source("code/22 prepare embeddings.R")
      source("code/23 train GP model.R")
      source("code/24 fill previous gaps.R")
      source("code/25 forecast.R")
      obs_df_ori_allvar[[2]]<-obs_df_ori
      forecast_df_ori_allvar[[2]]<-forecast_df_ori
      
      forecast_df_ori<-bind_rows(forecast_df_ori_allvar)
      
      source("code/26 output for submission.R")
      
      closeAllConnections()
      
      todaydone<-T
    },
    error = function(e) {
      todaydone<-F
      message(conditionMessage(e))
    })
  }
  
}


