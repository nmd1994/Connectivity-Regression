#STABILITY SELECTION FOR SIMULATIONS 

#NOTE: IF INTERESTED IN RERUNNING THE DATA USED TO GENERATE THE SOLUTIONS IN OUR PAPER,
#UNCOMMENT THE LINE BELOW. OTHERWISE, 
#SET YOUR SIMULATED DATA AS Y=YA, X=XA, B=BA, AND E=EMA.
#THIS IS ASSUMING DATA FOLLOWS THE FORMAT Y = XB + E. 

#------------------------------------------------
#READ IN DATA FROM THE SIMULATIONS SHOWN IN THE MANUSCRIPT HERE

#Dat=readRDS('n250rho02.rds')
#Ya=Dat[[1]]
#Xa=Dat[[2]]
#Ba=Dat[[3]]
#Ema=Dat[[4]]
#-----------------------------------------------
#OTHERWISE, SET YOUR RESPONSE MATRIX TO YA AND PREDICTOR MATRIX TO XA
#Ya=Nets
#Xa=Covs
#-----------------------------------------------
#Compute Stability Selection To Compute ROC Curves 
stabilitymSSL=list(0)
stabilitySLOPE=list(0)
stabilityLasso=list(0)

for(i in 1:length(Ya)){

Y=Ya[[i]]
X=Xa[[i]]

Gd=matrix(0,nrow=n,ncol=q)

for(j in 1:n){
  u=GFT_inv_map(Y[j,],1e-8)[[1]]
  if(!all(eigen(u)$values>=0)){
    u=nearPD(u)$mat
  }
  
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

Boots=100
library(doSNOW)
NumberOfCluster <- 4
cl <- makeCluster(NumberOfCluster)
registerDoSNOW(cl)

x0 <- foreach(i = 1:Boots) %dopar% {
  library('mSSL')
  library("SLOPE")
  library('glmnet')
  library('matrixStats')
  library('expm')
  
  e=sample.int(n,replace=TRUE)
  Xd=X[e,]
  Gd=Y[e,]
  
  red=mSSL_dpe(Xd,Gd)
  
  mSSL_stab=matrix(0,nrow=p,ncol=q)
  mSSL_stab[which(red$B!=0,arr.ind=TRUE)]=1
  
  B_mSSL=red$B
  
  Lasso_Beta=matrix(0,nrow=p,ncol=q)
  for(j in 1:ncol(Gd)){
    emd=cv.glmnet(Xd,Gd[,j],intercept=FALSE)
    bets=coef(emd,s="lambda.min")
    Lasso_Beta[,j]=bets[2:length(bets)]
    print(j)
  }
  
  Lasso_stab=matrix(0,nrow=p,ncol=q)
  Lasso_stab[which(Lasso_Beta!=0,arr.ind=TRUE)]=1
  
  
  SLOPE_Beta=matrix(0,nrow=p,ncol=q)
  for(j in 1:ncol(Gd)){
    rc=SLOPE(Xd,Gd[,j],family= "gaussian",lambda="bh",
             alpha="estimate",intercept=FALSE,q=0.01)
    SLOPE_Beta[,j]=rc$coefficients
    print(j)
  }
  
  SLOPE_stab=matrix(0,nrow=p,ncol=q)
  SLOPE_stab[which(SLOPE_Beta!=0,arr.ind=TRUE)]=1
  
  
  return(list("mSSL_stab"= mSSL_stab,"B_mSSL"= B_mSSL,
              "Lasso_stab"= Lasso_stab,"Lasso_Beta"= Lasso_Beta,"SLOPE_stab"= SLOPE_stab,
              "SLOPE_Beta"=SLOPE_Beta))
}
stopCluster(cl)


x1=list(0)
for(d in 1:Boots){
  x1[[d]]=x0[[d]]$mSSL_stab
  print(d)
}
stabcompute_mSSL=Reduce("+",x1)/length(x1)

stabilitymSSL[[i]]=stabcompute_mSSL

x1=list(0)
for(d in 1:Boots){
  x1[[d]]=x0[[d]]$SLOPE_stab
  print(d)
}
stabcompute_SLOPE=Reduce("+",x1)/length(x1)

stabilitySLOPE[[i]]=stabcompute_SLOPE

x1=list(0)
for(d in 1:Boots){
  x1[[d]]=x0[[d]]$Lasso_stab
  print(d)
}
stabcompute_Lasso=Reduce("+",x1)/length(x1)

stabilityLasso[[i]]=stabcompute_Lasso

print(i)
}

#If re-running our results, uncomment this portion to recreate 
#stability selection scores seen in paper 

#Sta=readRDS('stabilityselectionn250rho02.rds')
#stabilitymSSL=Sta[[1]]
#stabilitySLOPE=Sta[[2]]
#stabilityLasso=Sta[[3]]

#Collect Stats That You Want For 3 methods 
#NOTE: THIS IS SET UP TO PRODUCE RESULTS FOR FP = 0.01
#TO GET FP = 0.05, REPLACE 0.01 IN FOR LOOP WITH 0.05 AND CHANGE THE MULTIPLICATIVE
#FACTOR TO 20 

pAUC01mssl=c()
pAUC05mssl=c()
TP01mssl=c()
TP05mssl=c()
ssthresh01mssl=c()
ssthresh05mssl=c()

pAUC01lasso=c()
pAUC05lasso=c()
TP01lasso=c()
TP05lasso=c()
ssthresh01lasso=c()
ssthresh05lasso=c()

pAUC01SLOPE=c()
pAUC05SLOPE=c()
TP01SLOPE=c()
TP05SLOPE=c()
ssthresh01SLOPE=c()
ssthresh05SLOPE=c()

for(j in 1:100){
#B=Ba[[j]]
Truth=matrix(0,nrow=nrow(B),ncol=ncol(B))
Truth[which(B!=0,arr.ind=TRUE)]=1

ld=1e4
valROC=seq(0.005,1,length.out=ld)

ROCmssl=matrix(0,nrow=ld,ncol=4)
ROClasso=matrix(0,nrow=ld,ncol=4)
ROCSLOPE=matrix(0,nrow=ld,ncol=4)
colar=c(rep(colors()[26],ld/4),rep(colors()[33],ld/4),rep(colors()[21],ld/4),
        rep(colors()[49],ld/4))

for(i in 1:ld){
  
  TMSSL=matrix(0,nrow(B),ncol=ncol(B))
  TMSSL[which(stabilitymSSL[[j]]>=valROC[i],arr.ind=TRUE)]=1
  

  TLasso=matrix(0,nrow(B),ncol=ncol(B))
  TLasso[which(stabilityLasso[[j]]>valROC[i],arr.ind=TRUE)]=1
  
  TSLOPE=matrix(0,nrow(B),ncol=ncol(B))
  TSLOPE[which(stabilitySLOPE[[j]]>valROC[i],arr.ind=TRUE)]=1
  
  #COMPUTE MATRIX OF TRUEPOSITIVE AND TRUE NEGATIVES
  TindexT=which(Truth==1,arr.ind=TRUE)
  numTP=nrow(TindexT)
  TindexFP=which(Truth==0,arr.ind=TRUE)
  numFP=nrow(TindexFP)
  
  
  FPratemssl=sum(TMSSL[TindexFP])/numFP
  FPratelasso=sum(TLasso[TindexFP])/numFP
  FPrateSLOPE=sum(TSLOPE[TindexFP])/numFP
  
  TPratemssl=sum(TMSSL[TindexT])/numTP
  TPratelasso=sum(TLasso[TindexT])/numTP
  TPrateSLOPE=sum(TSLOPE[TindexT])/numTP
  
  ROCmssl[i,]=c(valROC[i],TPratemssl,FPratemssl,colar[i])
  ROClasso[i,]=c(valROC[i],TPratelasso,FPratelasso,colar[i])
  ROCSLOPE[i,]=c(valROC[i],TPrateSLOPE,FPrateSLOPE,colar[i])
  #print(i)
}

values=as.numeric(ROCmssl[,3])
val2=as.numeric(ROCmssl[,2])

minusmssl=values[which(values<=0.01)]
indmin=which.min(abs(minusmssl-0.01))
lowlow=minusmssl[indmin]
downy=mean(val2[which(values==lowlow)])
ROCAcmin=valROC[which(values==lowlow)]


plusmssl=values[which(values>0.01)]
indplus=which.min(abs(plusmssl-0.01))
upup=plusmssl[indplus]
uppy=mean(val2[which(values==upup)])
ROCAcmax=valROC[which(values==upup)]


FPs=values[max(which(values==lowlow)):length(values)]
TPs=val2[max(which(values==lowlow)):length(values)]
FPs=c(0.01,FPs)

dup=abs(0.01-upup)
ddown=abs(0.01-lowlow)
if(length(dup)!=0){

TP01mssl[j]=((uppy*dup)/(dup+ddown)) + ((downy*ddown)/(dup+ddown))
TPs=c(TP01mssl[j],TPs)

library(zoo)
x <- FPs
y <- TPs
id <- order(x)
AUC <- sum(diff(x[id])*rollmean(y[id],2))
pAUC01mssl[j]=AUC*100

ssthresh01mssl[j]=(sum(ROCAcmax)+sum(ROCAcmin))/(length(ROCAcmax)+length(ROCAcmin))
}

values=as.numeric(ROClasso[,3])
val2=as.numeric(ROClasso[,2])

minusmssl=values[which(values<=0.01)]
indmin=which.min(abs(minusmssl-0.01))
lowlow=minusmssl[indmin]
downy=mean(val2[which(values==lowlow)])
ROCAcmin=valROC[which(values==lowlow)]

plusmssl=values[which(values>0.01)]
indplus=which.min(abs(plusmssl-0.01))
upup=plusmssl[indplus]
uppy=mean(val2[which(values==upup)])
ROCAcmax=valROC[which(values==upup)]

FPs=values[max(which(values==lowlow)):length(values)]
TPs=val2[max(which(values==lowlow)):length(values)]
FPs=c(0.01,FPs)

dup=abs(0.01-upup)
ddown=abs(0.01-lowlow)

TP01lasso[j]=((uppy*dup)/(dup+ddown)) + ((downy*ddown)/(dup+ddown))
TPs=c(TP01lasso[j],TPs)

library(zoo)
x <- FPs
y <- TPs
id <- order(x)
AUC <- sum(diff(x[id])*rollmean(y[id],2))
pAUC01lasso[j]=AUC*100

ssthresh01lasso[j]=(sum(ROCAcmax)+sum(ROCAcmin))/(length(ROCAcmax)+length(ROCAcmin))

values=as.numeric(ROCSLOPE[,3])
val2=as.numeric(ROCSLOPE[,2])

minusmssl=values[which(values<=0.01)]
indmin=which.min(abs(minusmssl-0.01))
lowlow=minusmssl[indmin]
downy=mean(val2[which(values==lowlow)])
ROCAcmin=valROC[which(values==lowlow)]

plusmssl=values[which(values>0.01)]
indplus=which.min(abs(plusmssl-0.01))
upup=plusmssl[indplus]
uppy=mean(val2[which(values==upup)])
ROCAcmax=valROC[which(values==upup)]


FPs=values[max(which(values==lowlow)):length(values)]
TPs=val2[max(which(values==lowlow)):length(values)]
FPs=c(0.01,FPs)

dup=abs(0.01-upup)
ddown=abs(0.01-lowlow)

TP01SLOPE[j]=((uppy*dup)/(dup+ddown)) + ((downy*ddown)/(dup+ddown))
TPs=c(TP01SLOPE[j],TPs)

library(zoo)
x <- FPs
y <- TPs
id <- order(x)
AUC <- sum(diff(x[id])*rollmean(y[id],2))
pAUC01SLOPE[j]=AUC*100
ssthresh01SLOPE[j]=(sum(ROCAcmax)+sum(ROCAcmin))/(length(ROCAcmax)+length(ROCAcmin))

print(j)
}
