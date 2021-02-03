res<-PrepareEmbedding(x,start=max(unlist(lags))+1,end=length(date_list), focalsites = 1:nrow(coord_df), lags=lags, neighbors=neighbors,vars=vars, distMat = distMat)
dir.create(paste0(path,"train"))

for (moy in 1:12) {
  train_id<-(as.integer(format(as.Date(res$D), "%m")) ==moy)  & (res$D<= max(ts_all$time))
  X_train<-res$X[train_id,, drop=F]
  Y_train<-res$Y[train_id,, drop=F]
  D_train<-res$D[train_id,, drop=F]
  P_train<-res$P[train_id,, drop=F]
  dir.create(paste0(path,"train/MOY", moy))
  write_csv(as.data.frame(X_train),paste0(path,"train/MOY", moy,"/X.csv"))
  write_csv(as.data.frame(Y_train),paste0(path,"train/MOY", moy,"/Y.csv"))
  write_csv(as.data.frame(D_train),paste0(path,"train/MOY", moy,"/D.csv"))
  write_csv(as.data.frame(P_train),paste0(path,"train/MOY", moy,"/P.csv"))
  
  print(moy)
}
