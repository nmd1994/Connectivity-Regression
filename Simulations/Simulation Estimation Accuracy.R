#SIMULATION STUDY: Connectivity Regression 
#--------------------------------------
#STEP 1: DATA GENERATION PROCESS 
#--------------------------------------
#------------------------------------
#STEP 1a): SELECT HYPERPARAMETERS
#------------------------------------
#load relevant package dependencies
library('MASS')
library('matrixStats')
library('Matrix')
library('glmnet')
#library('mSSL')
library('SLOPE')
library('expm')
#------------------------------------
#sample size
n <- 150
#number of covariates
p <- 30
#number of edges 
q <- (10*9)/2
#underlying level of second order dependence
rho=0.9
#------------------------------------
#Generate Empty Lists To Save Simulated Data sets for this setting
sig=diag(p)
Nets=list()
Covs=list()
Bman=list()

#----------------------------------
#STEP 1b): SIMULATE REPLICATED DATASETS FOR EACH SETTING
#----------------------------------
for(k in 1:100){
X=mvrnorm(n,rep(0,p),sig)
bd <- rho^(abs(outer(1:q, 1:q, "-")))

#THIS IS WEIRD STRUCTURE 
bd=as.matrix(bdiag(list(0.9^(abs(outer(1:6, 1:6, "-"))), 0.7^(abs(outer(1:4, 1:4, "-"))), 0.5^(abs(outer(1:4, 1:4, "-"))),diag(45-14))))
#----------------------------------
#select which edges are sparse (SPARSITY LEVEL 1)
signal=round(0.3*q)
Xchoose=seq(1:signal)

#select sparsity within selected covariates (SPARSITY LEVEL 2)
wcov=list(0)
for(i in 1:signal){
wcov[[i]]=round(sample(1:p,p*0.3))
print(i)
}
#---------------------------------
#-----------------------------------
#select mix of coefficients: low, medium, high
prop=c((1/3),(1/3),(1/3)) 
totnum=length(wcov[[1]])*length(Xchoose)

  d=runif(round(prop[1]*totnum),0.1,0.5)
  neg=sample(c(-1,1),round(prop[1]*totnum),replace=TRUE)
  v1=d*neg

  e=runif(round(prop[2]*totnum),0.5,1)
  neg=sample(c(-1,1),round(prop[2]*totnum),replace=TRUE)
  v2=e*neg

  f=runif(totnum-length(e)-length(d),1,1.5)
  neg=sample(c(-1,1),prop[3]*totnum,replace=TRUE)
  v3=f*neg
  
values=c(v1,v2,v3)
values=sample(values)
#------------------------------------------------
#------------------------------------------------
#put together B
B=matrix(0,nrow=p,ncol=q)

set=seq(1,length(values)+1,length(wcov[[1]]))
for(j in 1:length(wcov)){
B[wcov[[j]],Xchoose[j]]=values[set[j]:(set[j+1]-1)]
print(j)
}
ind=sample(1:p,5)
B[ind,]=0

#------------------------------------------------
#----------------------------------
#Generate E 
#----------------------------------
E = mvrnorm(n,rep(0,q),as.matrix(bd))
#----------------------------------
#Generate Real Data And Save Datasets 
#----------------------------------
Y=X%*%B + E

Nets[[k]]=Y
Covs[[k]]=X
Bman[[k]]=B
print(k)
}

#-------------------------------------------------
#STEP 2: RUNNING SIMULATIONS -- COMPUTE MSE & BOOTSTRAP CALCULATIONS
#-------------------------------------------------
#---------------------------
#step 1a) Generate True Correlation Matrix (off-diagonal) given Gamma
#---------------------------
#FOR CORRELATION
mSSLBMSE=c()
#Cov_mSSL=list(0)
LassoBMSE=c()
OLSBMSE=c()
SLOPEBMSE=c()
LassoOLSBMSE=c()
SLOPEOLSEBMSE=c()
SLOPEOLSGLSBMSE=c()
Emp_Cov_Slope=list(0)
LassoOLSGLSBMSE=c()
Emp_Cov_Lasso=list(0)
prec_mSSL=list(0)
lengthjuju=c()
lengthslup=c()


#READ IN SIMULATED DATA HERE as dat

#dat=readRDS("n150rho02.rds")

for(z in 1:100){

#Y=dat[[1]][[z]]
#X=dat[[2]][[z]]
#B=dat[[3]][[z]]

  
Y=Nets[[z]]
X=Covs[[z]]
B=Bman[[z]]
  
Gd=matrix(0,nrow=n,ncol=q)

for(j in 1:n){
  u=GFT_inv_map(Y[j,],1e-8)[[1]]
  if(!all(eigen(u)$values>=0)){
    u=nearPD(u)$mat
  }
  
  #generate vectors to re-estimate correlation matrices (1000 time points)
  med=mvrnorm(1000,rep(0,10),u)
  mad=cor(med)
  if(!all(eigen(mad)$values>=0)){
    mad=nearPD(mad)
  }
  eerm=logm(mad)
  Gd[j,]=eerm[upper.tri(eerm)]
  print(j)
}

Y=Gd

#Compute Solution for mSSL 
ed=mSSL_dpe(X,Y)

B_MSE_mSSL=sum((B-ed$B)^(2))/(nrow(B)*ncol(B))
mSSLBMSE[z]=B_MSE_mSSL
prec_mSSL[[z]]=ed$Omega
#----------------------------------------
#step 1c) Compute solution for Lasso 
#----------------------------------------
Lasso_Beta=matrix(0,nrow=p,ncol=q)
Lasso_BetaOLS=matrix(0,nrow=p,ncol=q)
for(j in 1:ncol(Y)){
  emd=cv.glmnet(X,Y[,j],intercept=FALSE)
  bets=coef(emd,s="lambda.min")
  Lasso_Beta[,j]=bets[2:length(bets)]
  xind=which(coef(emd)!=0)
  if(length(xind)!=0){
  betsols=lm(Y[,j]~X[,(xind-1)]-1)
  Lasso_BetaOLS[(xind-1),j]=betsols$coefficients
  }
}

B_MSE_Lasso=sum((B-Lasso_Beta)^(2))/(nrow(B)*ncol(B))

LassoBMSE[z]=B_MSE_Lasso

B_MSE_LassoOLS=sum((B-Lasso_BetaOLS)^(2))/(nrow(B)*ncol(B))

LassoOLSBMSE[z]=B_MSE_LassoOLS

#----------------------------------
#COMPUTE SOLUTION FOR SLOPE METHOD (USING BH-PARAMETERS AND ALGO 5)
#----------------------------------
SLOPE_Beta=matrix(0,nrow=p,ncol=q)
SLOPE_BetaOLS=matrix(0,nrow=p,ncol=q)
for(j in 1:ncol(Y)){
  rc=SLOPE(X,Y[,j],family= "gaussian",lambda="bh", intercept=FALSE,alpha="estimate",q=0.01)
  xind=which(rc$coefficients!=0)
  if(length(xind)!=0){
  rcOLS=lm(Y[,j]~X[,xind]-1)
  SLOPE_Beta[,j]=rc$coefficients
  SLOPE_BetaOLS[xind,j]=rcOLS$coefficients
  }
}

B_MSE_SLOPE=sum((B-SLOPE_Beta)^(2))/(nrow(B)*ncol(B))
SLOPEBMSE[z]=B_MSE_SLOPE

B_MSE_SLOPEOLS=sum((B-SLOPE_BetaOLS)^(2))/(nrow(B)*ncol(B))
SLOPEOLSEBMSE[z]=B_MSE_SLOPEOLS

ir=unique(which(SLOPE_BetaOLS!=0,arr.ind=TRUE)[,2])
ip=unique(which(Lasso_BetaOLS!=0,arr.ind=TRUE)[,2])

#Compute Empirical Covariance for Both Methods 
Emp_Cov_Slope[[z]]=solve((t(Y[,ir]-X%*%SLOPE_BetaOLS[,ir])%*%(Y[,ir]-X%*%SLOPE_BetaOLS[,ir]))/(n))
Emp_Cov_Lasso[[z]]=solve((t(Y[,ip]-X%*%Lasso_BetaOLS[,ip])%*%(Y[,ip]-X%*%Lasso_BetaOLS[,ip]))/(n))

#Compute GLS for both slope and lasso 

#SLOPE - OUTCOME 
#For GLS, create block diagnoals of each matrix 
SLOPEgls_cov=kronecker(Emp_Cov_Slope[[z]],bdiag(diag(n)))
Lassogls_cov=kronecker(Emp_Cov_Lasso[[z]],bdiag(diag(n)))


#Create list of X based on the sleleciotn -ERROR IS RIGHT HERE -NEED TO FIX THIS 
X_select_SLOPE=list(0)
lengths=c()
for(i in 1:length(ir)){
  X_select_SLOPE[[i]]=X[,which(SLOPE_BetaOLS[,ir[i]]!=0)]
  lengths[i]=length(which(SLOPE_BetaOLS[,ir[i]]!=0))
}


X_sel_SLOPE=bdiag(X_select_SLOPE)
B_SLOPE_GLS=solve(as.matrix(t(X_sel_SLOPE))%*%SLOPEgls_cov%*%as.matrix(X_sel_SLOPE))%*%t(X_sel_SLOPE)%*%SLOPEgls_cov%*%as.vector(Y[,ir])


juju=matrix(0,nrow=p,ncol=q)
lengthjuju[z]=length(juju[which(SLOPE_BetaOLS!=0,arr.ind=TRUE)])
juju[which(SLOPE_BetaOLS!=0,arr.ind=TRUE)]=as.vector(B_SLOPE_GLS)
lengthslup[z]=length(as.vector(B_SLOPE_GLS))

B_MSE_SLOPEGLS=sum((B-juju)^(2))/(nrow(B)*ncol(B))
SLOPEOLSGLSBMSE[z]=B_MSE_SLOPEGLS

#LASSO - OUTCOME 
#For GLS, create block diagnoals of each matrix 
#Lassogls_cov=bdiag(rep(list(Emp_Cov_Lasso[[z]]),n))
#Create list of X based on the sleleciotn 
X_select_Lasso=list(0)
lengthslass=c()
for(i in 1:length(ip)){
  X_select_Lasso[[i]]=X[,which(Lasso_BetaOLS[,ip[i]]!=0)]
  lengthslass[i]=length(which(Lasso_BetaOLS[,ip[i]]!=0))
}

X_sel_Lasso=bdiag(X_select_Lasso)

B_Lasso_GLS=solve(as.matrix(t(X_sel_Lasso))%*%(Lassogls_cov)%*%as.matrix(X_sel_Lasso))%*%t(X_sel_Lasso)%*%Lassogls_cov%*%as.vector(Y[,ip])

juju=matrix(0,nrow=p,ncol=q)
juju[which(Lasso_BetaOLS!=0,arr.ind=TRUE)]=as.vector(B_Lasso_GLS)

B_MSE_LassoGLS=sum((B-juju)^(2))/(nrow(B)*ncol(B))
LassoOLSGLSBMSE[z]=B_MSE_LassoGLS

#---------------------------------
#step 1e) OLS: Compute MSE for B, Gamma, and Predicted Correlation Matrices (off-diagonal)
#---------------------------------
OLS_Beta=matrix(0,nrow=p,ncol=q)
for(j in 1:ncol(Y)){
  sol=lm(Y[,j]~X-1)
  OLS_Beta[,j]=sol$coefficients
}

B_MSE_OLS=sum((B-OLS_Beta)^(2))/(nrow(B)*ncol(B))
OLSBMSE[z]=B_MSE_OLS

print(z)
}








#Compute Precision Off-Diagonal Estimtes 
precMSE_mSSL=c()
precMSE_SLOPE=c()
precMSE_Lasso=c()
for(i in 1:100){
  precMSE_mSSL[i]=mean((prec_mSSL[[i]][upper.tri(prec_mSSL[[i]])]-solve(bd)[upper.tri(solve(bd))])^(2))
  #precMSE_mSSL[i]=mean((solve(prec_mSSL[[i]])-bd)^(2))
  #precMSE_SLOPE[i]=mean((Emp_Cov_Slope[[i]][upper.tri(Emp_Cov_Slope[[i]])]-solve(bd)[upper.tri(solve(bd))])^(2))
  #precMSE_Lasso[i]=mean((Emp_Cov_Lasso[[i]][upper.tri(Emp_Cov_Lasso[[i]])]-solve(bd)[upper.tri(solve(bd))])^(2))
}







