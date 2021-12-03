for(site in 1:nrow(coord_df)) {
  path_train<-paste0(path, focal_var,"/train/",site)
  dir.create(path_train, recursive=T)
  res<-PrepareEmbedding(x,start=max(unlist(lags))+1,end=length(date_list), focalsites =site, lags=lags, neighbors=neighbors,vars=vars, distMat = distMat)
  X_train<-res$X
  Y_train<-res$Y
  D_train<-res$D
  P_train<-res$P
  
  write_csv(as.data.frame(X_train),paste0(path_train,"/X.csv"))
  write_csv(as.data.frame(Y_train),paste0(path_train,"/Y.csv"))
  write_csv(as.data.frame(D_train),paste0(path_train,"/D.csv"))
  write_csv(as.data.frame(P_train),paste0(path_train,"/P.csv"))
}

