#Permutation Test: new revamped edition (now try without div by 10)
p=7
q=105
nperm=250
n=1003

egd=as.matrix(Covs)
egd=egd[,-1]
Gamd=as.matrix(NormGam)
Gamd=Gamd[,-1]
avgd=avgd[,-1]
cnote=cnote[,-1]

library(doSNOW)
NumberOfCluster <- 4
cl <- makeCluster(NumberOfCluster)
registerDoSNOW(cl)
x0 <- foreach(i = 1:nperm) %dopar% {
  library('mSSL')
  library('matrixStats')
  library('expm')
  library('gtools')

  eddy=permute(seq(1,nrow(egd)))
  edr=egd[eddy,]
  Xd=edr/10
  Gd=Gamd
  
  red=mSSL_dpe(as.matrix(Xd),Gd)
  
  mSSL_stab=matrix(0,nrow=p,ncol=q)
  mSSL_stab[which(red$B!=0,arr.ind=TRUE)]=1
  
  B_mSSL=red$B
  
  return(list("B_mSSL"= B_mSSL, "mSSL_stab"= mSSL_stab))
}
stopCluster(cl)

Xd=egd/10
Gd=Gamd

rede=mSSL_dpe(as.matrix(Xd),Gd)

mSSL_stabe=matrix(0,nrow=p,ncol=q)
mSSL_stabe[which(rede$B!=0,arr.ind=TRUE)]=1


#test statistics 
tstatcollectmat_mSSL=c()
tstatcollectmat_Lasso=c()
tstatcollectmat_SLOPE=c()

dr=which(mSSL_stabe==1,arr.ind=TRUE)
for(d in 1:p){
  #vecy[d]=length(dr[,1][dr[,1]==d])
  tstatcollectmat_mSSL[d]=length(dr[,1][dr[,1]==d])
}

gtest2_mSSL_permute=c()

#Covariate Wise Test 
collectmat_mSSL_permute=matrix(0,nrow=p,ncol=nperm)

for(f in 1:nperm){
  dr=which(x0[[f]]$mSSL_stab==1,arr.ind=TRUE)
  for(d in 1:p){
    collectmat_mSSL_permute[d,f]=length(dr[,1][dr[,1]==d])
  }
  gtest2_mSSL_permute[f]=max(collectmat_mSSL_permute[,f])
}


pvaluemSSL_multadj=c()
for(d in 1:p){
  pvaluemSSL_multadj[d]=length(which(gtest2_mSSL_permute>=tstatcollectmat_mSSL[d]))/length(gtest2_mSSL_permute)
}
