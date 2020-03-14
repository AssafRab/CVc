### An R code for the California housing example. The dataset is available in scikit-learn package.
### This code was executed 50 times (with differnt seeds). 

#Import relevant library
library(dplyr)

## Path
Data_folder<-"your data folder path"

## Function for estimating the covariance parameters
Loss <- function(par,X,Coor,Y,Block) {
  signal<-par[1]
  eta<-par[2]
  lambda1<-par[3]
  lambda2<-par[4]
  SigSq<-par[5]
  BlockSig<-par[6]
  
  rotCos<-cos(eta)
  rotSin<-sin(eta)
  rotation<-matrix(c(rotCos,rotSin,-rotSin,rotCos),2,2)
  scaling<-matrix(c(lambda1,0,0,lambda2),2,2)
  
  K<-matrix(NA,length(Y),length(Y))
  for(col in 1:length(Y)){
    s<-cbind(Coor[,1]-Coor[col,1],Coor[,2]-Coor[col,2])
    K[,col]<-signal*exp(-sqrt(diag(s%*%chol2inv(chol(rotation%*%scaling%*%t(rotation)))%*%t(s))))
  }
  Var<-K+SigSq*diag(length(Y))
  Var_inv<-chol2inv(chol(Var))
  mean<-X%*%chol2inv(chol(t(X)%*%Var_inv%*%X))%*%t(X)%*%Var_inv%*%Y
  return(t(Y-mean)%*%Var_inv%*%(Y-mean)+sum(diag(log(Var))))
}


#Prepering the data
{

  housing_price<-read.csv(file = paste0(Data_folder,"/housing_price.csv",sep=""))[,-1]
  housing_price<-housing_price[,c(9,7:8,1:6)]
  housing_price<-housing_price %>%mutate(ID = group_indices_(housing_price, .dots=c("Latitude","Longitude"))) 

  # Removing missing values and selecting clusters with 2 observations or more
  housing_price<-housing_price%>%group_by(ID)%>%mutate(n = sum(!is.na(MedInc)))
  housing_price<-housing_price[housing_price$n>1,]
  Blocks<-unique(housing_price[housing_price$n>1,]$ID) 
  
  # Sampling training and test
  set.seed(as.numeric(1))
  sample_cluster<-sample(Blocks, size=700, replace = F)
  tr<-housing_price[housing_price$ID %in% sample_cluster,]
  test<-housing_price[!(housing_price$ID %in% sample_cluster),]
  tr<-tr[order(tr$ID),]
  test<-test[order(test$ID),]
  
  BlockDiag<-model.matrix(~as.factor(ID)-1,tr)%*%t(model.matrix(~as.factor(ID)-1,tr))
  y_test<-unlist(test[,1])
  Coor_test<-as.matrix(test[,c(2,3)])
  y_tr<-unlist(tr[,1])
  Coor_tr<-as.matrix(tr[,c(2,3)])
}


#Defining the models
{
  model_i<-list(c(4),c(4:6),c(4:9))
  model_no<-length(model_i)
  n_tr<-nrow(tr)
  n_test<-nrow(test)
}

#Defining the output matrices
{
  Error_tr<-matrix(NA,1,model_no)
  Correction<-matrix(NA,1,model_no)
  Error_Gen<-matrix(NA,1,model_no)
}


for(j in 1:model_no){
  
  X_test<-as.matrix(cbind(rep(1,n_test),scale(test[,unlist(model_i[j])])))
  X_tr<-as.matrix(cbind(rep(1,n_tr),scale(tr[,unlist(model_i[j])])))
 
  #Estimating the cov matrix
  X_tr_model<-as.matrix(scale(tr[,unlist(model_i[j])]))
  result<-optim(par=c(1,1,1,1,1,1), Loss,X=X_tr, Coor=Coor_tr,
                Y=y_tr,Block=BlockDiag,control = list(maxit=10),method="L-BFGS-B",lower=c(0.001,0,0.001,0.001,0.001,0.001), upper=c(10,pi/2,10,10,10,10))
  
  signal<-result$par[1]
  eta<-result$par[2]
  lambda1<-result$par[3]
  lambda2<-result$par[4]
  SigSq<-result$par[5]
  BlockSig<-result$par[6]
  
  rotCos<-cos(eta)
  rotSin<-sin(eta)
  rotation<-matrix(c(rotCos,rotSin,-rotSin,rotCos),2,2)
  scaling<-matrix(c(lambda1,0,0,lambda2),2,2)
  K<-matrix(NA,length(y_tr),length(y_tr))
  for(col in 1:length(y_tr)){
    s<-cbind(Coor_tr[,1]-Coor_tr[col,1],Coor_tr[,2]-Coor_tr[col,2])
    K[,col]<-signal*exp(-sqrt(diag(s%*%chol2inv(chol(rotation%*%scaling%*%t(rotation)))%*%t(s))))
  }
  V<-K+SigSq*diag(length(y_tr))
 
 # Creating the covaiance matrix between the training and test 
  K_star<-matrix(NA,length(y_test),length(y_tr))
  
  for(col in 1:n_tr){
    ss<-cbind(Coor_test[,1]-Coor_tr[col,1],Coor_test[,2]-Coor_tr[col,2])
    K_star[,col]<-signal*exp(-sqrt(diag(ss%*%chol2inv(chol(rotation%*%scaling%*%t(rotation)))%*%t(ss))))
  }
  
  
  #Creating H_cv
  H_CV<-matrix(0,n_tr,n_tr)
  
  for(i in 1:n_tr){
    X_tr_i<-X_tr[i,]
    X_tr_minus_i<-X_tr[-i,]
    V_minus_i<-as.matrix(V[-i,-i])
    V_minus_i_inv<-chol2inv(chol(V_minus_i))
    H_CV[i,-i]<-X_tr_i%*%chol2inv(chol(t(X_tr_minus_i)%*%V_minus_i_inv%*%X_tr_minus_i))%*%t(X_tr_minus_i)%*%V_minus_i_inv
  }
  
  H_CVTmp<-matrix(0,n_tr,(n_tr-1))
  
  for(i in 1:n_tr){
    H_CVTmp[i,]<-H_CV[i,-i]
  }
  for(i in 1:n_tr){
    H_CV[i,-i]<-H_CV[i,-i]+K[i,-i]%*%chol2inv(chol(V[-i,-i]))%*%(diag(1,n_tr-1)-H_CVTmp[-i,])
  }
  
  
  #CV +correction resutls
  Error_tr[j]<-1/n_tr*(t(y_tr-H_CV%*%y_tr)%*%(y_tr-H_CV%*%y_tr))
  Correction[j]<-2*(K[1,1]-mean(K-K[1,1]*BlockDiag))/n_tr*(sum(diag(H_CV%*%BlockDiag)))
  

  # Generaliztion Error
  Var_model_j_inv<-chol2inv(chol(V))
  H_gen<-X_test%*%chol2inv(chol(t(X_tr)%*%Var_model_j_inv%*%X_tr))%*%t(X_tr)%*%Var_model_j_inv+K_star%*%Var_model_j_inv%*%(diag(1,n_tr)-X_tr%*%chol2inv(chol(t(X_tr)%*%Var_model_j_inv%*%X_tr))%*%t(X_tr)%*%Var_model_j_inv)
  Error_Gen[j]<-1/n_test*(t(y_test-H_gen%*%y_tr)%*%(y_test-H_gen%*%y_tr))
  
}
