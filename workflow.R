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


