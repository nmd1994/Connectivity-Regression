library('MASS')
library('mSSL')
library('matrixStats')
library('expm')
library('gtools')
library('SLOPE')
library('glmnet')

#----------------------------
#Set Observed Response Matrix (or list of Response Matrices) as Nets
#Set Observed Predictor Matrix (or list of Predictor Matrices) as Covs 
#-----------------------------
#mixed p-values 
mSSL_pvalue_multadj=matrix(0,nrow=p,ncol=100)
Lasso_pvalue_multadj=matrix(0,nrow=p,ncol=100)
SLOPE_pvalue_multadj=matrix(0,nrow=p,ncol=100)

#null p-values-global type I error 
mSSL_pvalue_multadj_null=matrix(0,nrow=p,ncol=100)
Lasso_pvalue_multadj_null=matrix(0,nrow=p,ncol=100)
SLOPE_pvalue_multadj_null=matrix(0,nrow=p,ncol=100)
 
for(j in 1:length(Nets)){
  #Get Predicted Correlation Matrices, simulate from mvn normal and
  #re-estimate Matrix From Random Normal 
  Y=Nets[[j]]
  Gd=matrix(0,nrow=n,ncol=q)
  
  for(i in 1:n){
    u=GFT_inv_map(Y[i,],1e-8)[[1]]
    if(!all(eigen(u)$values>=0)){
      u=nearPD(u)$mat
    }
  
    med=mvrnorm(1000,rep(0,10),u)
    mad=cor(med)
    if(!all(eigen(mad)$values>=0)){
      mad=nearPD(mad)
    }
    eerm=logm(mad)
    Gd[i,]=eerm[upper.tri(eerm)]
    print(i)
  }
  
  Xd=Covs[[j]]

  library(doSNOW)
  NumberOfCluster <- 4
  cl <- makeCluster(NumberOfCluster)
  registerDoSNOW(cl)
  nperm=250
  x0 <- foreach(i = 1:nperm) %dopar% {
    library('mSSL')
    library('matrixStats')
    library('expm')
    library('gtools')
    library('SLOPE')
    library('glmnet')
    eddy=permute(seq(1,n))
    Xd=Xd[eddy,]
    red=mSSL_dpe(as.matrix(Xd),Gd)
    
    mSSL_stab=matrix(0,nrow=p,ncol=q)
    mSSL_stab[which(red$B!=0,arr.ind=TRUE)]=1
    
    B_mSSL=red$B
    
    Lasso_Beta=matrix(0,nrow=p,ncol=q)
    for(k in 1:ncol(Gd)){
      emd=cv.glmnet(Xd,Gd[,k],intercept=FALSE)
      bets=coef(emd,s="lambda.min")
      Lasso_Beta[,k]=bets[2:length(bets)]
      #print(k)
    }
    
    Lasso_stab=matrix(0,nrow=p,ncol=q)
    Lasso_stab[which(Lasso_Beta!=0,arr.ind=TRUE)]=1
    
    
    SLOPE_Beta=matrix(0,nrow=p,ncol=q)
    for(r in 1:ncol(Gd)){
      rc=SLOPE(Xd,Gd[,r],family= "gaussian",lambda="bh",
               alpha="estimate",intercept=FALSE,q=0.01)
      SLOPE_Beta[,r]=rc$coefficients
      #print(j)
    }
    
    SLOPE_stab=matrix(0,nrow=p,ncol=q)
    SLOPE_stab[which(SLOPE_Beta!=0,arr.ind=TRUE)]=1
    
    return(list("mSSL_stab"= mSSL_stab,"B_mSSL"= B_mSSL,"B_SLOPE"= SLOPE_Beta, "SLOPE_stab" = SLOPE_stab,
                "Lasso_Beta" = Lasso_Beta, "Lasso_stab" = Lasso_stab))
  }
  stopCluster(cl)
  
  #COMPUTE TEST STATISTICS 
  Xd=Covs[[j]]
  Gd=Gd
  
  rede=mSSL_dpe(as.matrix(Xd),Gd)
  
  mSSL_stabe=matrix(0,nrow=p,ncol=q)
  mSSL_stabe[which(rede$B!=0,arr.ind=TRUE)]=1
  
  SLOPE_Betae=matrix(0,nrow=p,ncol=q)
  for(r in 1:ncol(Gd)){
    rce=SLOPE(Xd,Gd[,r],family= "gaussian",lambda="bh",
              alpha="estimate",intercept=FALSE,q=0.01)
    SLOPE_Betae[,r]=rce$coefficients
  }
  
  SLOPE_stabe=matrix(0,nrow=p,ncol=q)
  SLOPE_stabe[which(SLOPE_Betae!=0,arr.ind=TRUE)]=1
  
  Lasso_Betae=matrix(0,nrow=p,ncol=q)
  for(k in 1:ncol(Gd)){
    emde=cv.glmnet(Xd,Gd[,k],intercept=FALSE)
    bets=coef(emde,s="lambda.min")
    Lasso_Betae[,k]=bets[2:length(bets)]
  }
  
  Lasso_stabe=matrix(0,nrow=p,ncol=q)
  Lasso_stabe[which(Lasso_Betae!=0,arr.ind=TRUE)]=1
  
  tstatcollectmat_mSSL=c()
  tstatcollectmat_Lasso=c()
  tstatcollectmat_SLOPE=c()
  
  dr=which(mSSL_stabe==1,arr.ind=TRUE)
  for(d in 1:p){
    tstatcollectmat_mSSL[d]=length(dr[,1][dr[,1]==d])
  }
  
  dr=which(Lasso_stabe==1,arr.ind=TRUE)
  for(d in 1:p){
    tstatcollectmat_Lasso[d]=length(dr[,1][dr[,1]==d])
  }
  
  dr=which(SLOPE_stabe==1,arr.ind=TRUE)
  for(d in 1:p){
    tstatcollectmat_SLOPE[d]=length(dr[,1][dr[,1]==d])
  }
  
  gtest2_mSSL_permute=c()
  gtest2_Lasso_permute=c()
  gtest2_SLOPE_permute=c()
  
  #Covariate Wise Test 
  collectmat_mSSL_permute=matrix(0,nrow=p,ncol=nperm)
  collectmat_Lasso_permute=matrix(0,nrow=p,ncol=nperm)
  collectmat_SLOPE_permute=matrix(0,nrow=p,ncol=nperm)
  
  for(f in 1:nperm){
    dr=which(x0[[f]]$mSSL_stab==1,arr.ind=TRUE)
    
    for(d in 1:p){
      collectmat_mSSL_permute[d,f]=length(dr[,1][dr[,1]==d])
    }
    
    dr=which(x0[[f]]$Lasso_stab==1,arr.ind=TRUE)
    for(d in 1:p){
      collectmat_Lasso_permute[d,f]=length(dr[,1][dr[,1]==d])
    }
    
    dr=which(x0[[f]]$SLOPE_stab==1,arr.ind=TRUE)
    for(d in 1:p){
      collectmat_SLOPE_permute[d,f]=length(dr[,1][dr[,1]==d])
    }
    
    #Global Test 2
    gtest2_mSSL_permute[f]=max(collectmat_mSSL_permute[,f])
    gtest2_Lasso_permute[f]=max(collectmat_Lasso_permute[,f])
    gtest2_SLOPE_permute[f]=max(collectmat_SLOPE_permute[,f])
  }
  
  
  #compute and save p-values 
  pvaluemSSL_multadj=c()
  pvalueLasso_multadj=c()
  pvalueSLOPE_multadj=c()
  
  for(d in 1:p){
    pvaluemSSL_multadj[d]=length(which(gtest2_mSSL_permute>=tstatcollectmat_mSSL[d]))/length(gtest2_mSSL_permute)
    pvalueLasso_multadj[d]=length(which(gtest2_Lasso_permute>=tstatcollectmat_Lasso[d]))/length(gtest2_Lasso_permute)
    pvalueSLOPE_multadj[d]=length(which(gtest2_SLOPE_permute>=tstatcollectmat_SLOPE[d]))/length(gtest2_SLOPE_permute)
  }
  
  mSSL_pvalue_multadj[,j]= pvaluemSSL_multadj
  Lasso_pvalue_multadj[,j]= pvalueLasso_multadj
  SLOPE_pvalue_multadj[,j]= pvalueSLOPE_multadj
  
  print(j)
}

#Collect p-values for actual signal and non-signal
truepval_mSSL=c()
truepval_Lasso=c()
truepval_SLOPE=c()
t1error_mSSL=c()
t1error_Lasso=c()
t1error_SLOPE=c()


for(i in 1:100){
  #edum=sort(unique(which(Bman[[i]]!=0,arr.ind=TRUE)[,1]))
  edum=sort(unique(which(B!=0,arr.ind=TRUE)[,1]))
  truepval_mSSL=c(truepval_mSSL,mSSL_pvalue_multadj[edum,i])
  truepval_Lasso=c(truepval_Lasso,Lasso_pvalue_multadj[edum,i])
  truepval_SLOPE=c(truepval_SLOPE,SLOPE_pvalue_multadj[edum,i])
  
  
  
  #edun=which(!seq(1,30,1)%in%sort(unique(which(Bman[[i]]!=0,arr.ind=TRUE)[,1])))
  edun=which(!seq(1,10,1)%in%sort(unique(which(B!=0,arr.ind=TRUE)[,1])))
  t1error_mSSL=c(t1error_mSSL,min(mSSL_pvalue_multadj[edun,i]))
  t1error_Lasso=c(t1error_Lasso,min(Lasso_pvalue_multadj[edun,i]))
  t1error_SLOPE=c(t1error_SLOPE,min(SLOPE_pvalue_multadj[edun,i]))
  print(i)
}

#Just get rid of min to get covariate level type 1 error 

#Covariate wise tests 
#power 
#mSSL
length(which(truepval_mSSL<0.05))/length(truepval_mSSL)

#Lasso
length(which(truepval_Lasso<0.05))/length(truepval_Lasso)

#SLOPE 
length(which(truepval_SLOPE<0.05))/length(truepval_SLOPE)


#type 1 error-experimentwsie for null
#mSSL
length(which(t1error_mSSL<0.05))/length(t1error_mSSL)

#Lasso
length(which(t1error_Lasso<0.05))/length(t1error_Lasso)

#SLOPE 
length(which(t1error_SLOPE<0.05))/length(t1error_SLOPE)

