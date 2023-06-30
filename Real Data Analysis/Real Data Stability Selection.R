#Permutation Test: new revamped edition (now try without div by 10)
p=7
q=105
BOOT=200
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
x0 <- foreach(i = 1:BOOT) %dopar% {
  library('mSSL')
  library('matrixStats')
  library('expm')
  library('gtools')

  e=sample.int(n,replace=TRUE)
  Xd=egd[e,]/10
  Gd=Gamd[e,]
  red=mSSL_dpe(as.matrix(Xd),Gd)
  
  mSSL_stab=matrix(0,nrow=p,ncol=q)
  mSSL_stab[which(red$B!=0,arr.ind=TRUE)]=1
  
  B_mSSL=red$B
  
  Omega_mSSL=red$Omega
  
  Omega_stab=matrix(0,nrow=q,ncol=q)
  Omega_stab[which(red$Omega!=0,arr.ind=TRUE)]=1
  
  return(list("B_mSSL"= B_mSSL, "mSSL_stab"= mSSL_stab,
              "Omega_mSSL" = Omega_mSSL, "Omega_stab" = Omega_stab))
}
stopCluster(cl)


# mSSL-bootstrap stability selection - BETA MATRIX
x1=list(0)
for(d in 1:BOOT){
  x1[[d]]=x0[[d]]$mSSL_stab
  #print(d)
}
stabcompute_mSSL=Reduce("+",x1)/length(x1)



data=as.data.frame(cbind(c(rep(1,105),rep(2,105),rep(3,105),rep(4,105),rep(5,105),
rep(6,105),rep(7,105)),c(stabcompute_mSSL[1,],stabcompute_mSSL[2,],stabcompute_mSSL[3,],
  stabcompute_mSSL[4,],stabcompute_mSSL[5,],stabcompute_mSSL[6,],
  stabcompute_mSSL[7,]),c(rep(1:105,7))))

derp=c(rep("ReadEng",105), rep("PicVocab",105),
                                          rep("Precuneus",105),
                                          rep("Postcentral",105),
                                          rep("Cuneus",105),
                                          rep("Posteriorcingulate",105),
                                          rep("Temporalpole",105))







# "POSTERIOR" BETA MATRIX MEAN
x2=list(0)
for(d in 1:BOOT){
  x2[[d]]=x0[[d]]$B_mSSL
}
Post_Beta=Reduce("+",x2)/length(x2)


# "POSTERIOR" PRECISION MEAN
x2=list(0)
for(d in 1:BOOT){
  x2[[d]]=x0[[d]]$Omega_mSSL
}
Post_Omega_mSSL=Reduce("+",x2)/length(x2)

# "POSTERIOR" PRECISION STABILITY SELECTION
x2=list(0)
for(d in 1:BOOT){
  x2[[d]]=x0[[d]]$Omega_stab
}
Omega_stability=Reduce("+",x2)/length(x2)

diag(Omega_stability)=0

data=as.data.frame(cbind(c(sort(rep(1:105,105))),c(rep(1:105,105)),
           c(Omega_stability[,1:105])))


library("pheatmap")
library(RColorBrewer)
library(gplots)
pheatmap(stabcompute_mSSL[c(1,3,4),],
         fontsize=20,
         colorRampPalette(rev(brewer.pal(n = 7, name =
                                           "RdYlBu")))(1000),cluster_cols=FALSE,cluster_rows=FALSE)









#Supplemental Table Creation 
indices=matrix(0,nrow=15,ncol=15)
indices[upper.tri(indices)]=seq(1,105,1)
indices=t(indices)+indices

oi=which(Omega_stability>=0.5&abs(Post_Omega_mSSL)>=0.25,
        arr.ind=TRUE)

nodeset1=c()
nodeset2=c()
for(i in 1:nrow(oi)){
    re=which(indices==oi[i,1],arr.ind=TRUE)
    nodeset1[i]=paste0("(",re[1],",",re[2],")")

    re=which(indices==oi[i,2],arr.ind=TRUE)
    nodeset2[i]=paste0("(",re[1],",",re[2],")")
  
}

staboo=c()
for(i in 1:nrow(oi)){
  staboo[i]=Omega_stability[oi[i,1],oi[i,2]]
}

partial=c()
for(i in 1:nrow(oi)){
  partial[i]=Post_Omega_mSSL[oi[i,1],oi[i,2]]
}

indv=rev(order(staboo))

indf=rev(order(abs(partial)))

mutec[indf,]


mutec=cbind(nodeset1,nodeset2,staboo,round(partial,3))

muted=mutec[indf,]

muted=muted[orc,]
orc=seq(1,222,2)

Edge=seq(1,105,1)

ror=cbind(Edge,t(stabcompute_mSSL))

write.csv(ror,file="covedge.csv")

write.csv(muted,
          file="nodezc2.csv")

ror=seq(1,105,1)

roed=matrix(0,nrow=15,ncol=15)
