dir.create(paste0(path,"train"))
res<-PrepareEmbedding(x,start=max(unlist(lags))+1,end=length(date_list), focalsites =1:nrow(coord_df), lags=lags, neighbors=neighbors,vars=vars, distMat = distMat)
X_train<-res$X
Y_train<-res$Y
D_train<-res$D
P_train<-res$P

write_csv(as.data.frame(X_train),paste0(path,"train/X.csv"))
write_csv(as.data.frame(Y_train),paste0(path,"train/Y.csv"))
write_csv(as.data.frame(D_train),paste0(path,"train/D.csv"))
write_csv(as.data.frame(P_train),paste0(path,"train/P.csv"))