# Load required libraries
library(MixGHD)
library(mixture)
library(cluster)
library(bclust)
library(mnormt)
library(EMMIXskew)
library(mclust)

# Generate simulation data from gaussian distribution with different parameter values

### Scenario 1 : Clusters with and without correlation

set.seed(123)
resARI1=matrix(0,10,4)
ressil1=matrix(0,10,4)
lab=c(rep(1,180),rep(2,180))

sig=matrix(0.8,15,15)
diag(sig)=1
start<-Sys.time()
for(i in 1:10){
  c1=rmnorm(180,mean=rep(0,15),diag(15))
  c2=rmnorm(180,mean=rep(3,15),varcov=sig)
  data=rbind(c1,c2)
  #plot(data,col=lab)
  t1= gpcm(data,G=2,mnames="VVV")
  t2= pam(data,k=2)
  t3= bclust(data,transformed.par = c(0, -10, log(25), 0, 0, 0))
  t4= kmeans(data,centers=2,nstart=10)
  
  resARI1[i,1]=ARI(t1$map,lab)
  resARI1[i,2]=ARI(t2$clustering,lab)
  resARI1[i,3]=ARI(t3$optim.alloc,lab)
  resARI1[i,4]=ARI(t4$cluster,lab)
  
  sil1<-silhouette(x=t1$map,dist=dist(data))
  ressil1[i,1]=mean(sil1[,3])
  sil2<-silhouette(x=t2$clustering,dist=dist(data))
  ressil1[i,2]=mean(sil2[,3])
  sil3<-silhouette(x=t3$optim.alloc,dist=dist(data))
  ressil1[i,3]=mean(sil3[,3])
  sil4<-silhouette(x=t4$cluster,dist=dist(data))
  ressil1[i,4]=mean(sil4[,3])
  
}
# Compute mean ARI for 10 iterations
res1<-apply(resARI1,2,mean)
resw1<-apply(ressil1,2,mean)
end<-Sys.time()
print(end-start)

### Scenario 2: Clusters with and without overlap

set.seed(123)
resARI2=matrix(0,10,4)
ressil2=matrix(0,10,4)
lab=c(rep(1,180),rep(2,180))

sig=matrix(0.3,15,15)
diag(sig)=1
start<-Sys.time()
for(i in 1:10){
  c1=rmnorm(180,mean=rep(0,15),diag(15))
  c2=rmnorm(180,mean=rep(1.5,15),varcov=sig)
  data=rbind(c1,c2)
  #plot(data,col=lab)
  t1= gpcm(data,G=2,mnames="VVV")
  t2= pam(data,k=2)
  t3=bclust(data,transformed.par = c(0, -10, log(25), 0, 0, 0))
  t4= kmeans(data,centers=2,nstart=10)
  
  resARI2[i,1]=ARI(t1$map,lab)
  resARI2[i,2]=ARI(t2$clustering,lab)
  resARI2[i,3]=ARI(t3$optim.alloc,lab)
  resARI2[i,4]=ARI(t4$cluster,lab)
  
  sil1<-silhouette(x=t1$map,dist=dist(data))
  ressil2[i,1]=mean(sil1[,3])
  sil2<-silhouette(x=t2$clustering,dist=dist(data))
  ressil2[i,2]=mean(sil2[,3])
  sil3<-silhouette(x=t3$optim.alloc,dist=dist(data))
  ressil2[i,3]=mean(sil3[,3])
  sil4<-silhouette(x=t4$cluster,dist=dist(data))
  ressil2[i,4]=mean(sil4[,3])
}
# Compute mean ARI for 10 iterations
res2<-apply(resARI2,2,mean)
resw2<-apply(ressil2,2,mean)
end<-Sys.time()
print(end-start)

### Scenario 3: Clusters with and without skewness

set.seed(123)
resARI3=matrix(0,10,4)
ressil3=matrix(0,10,4)
lab=c(rep(1,180),rep(2,180))

start<-Sys.time()
for(i in 1:10){
  c1=rmnorm(180,mean=rep(0,15),diag(15))
  c2=rdmsn(180,15,mean=rep(1.5,15),cov=diag(15),del=c(5,4,1.5,2 ,-3, -3, 2, 5, 2, 2, -4, -3, 1,5,5))
  data=rbind(c1,c2)
  #plot(data,col=lab)
  t1= gpcm(data,G=2,mnames="VVV")
  t2= pam(data,k=2)
  t3=bclust(data,transformed.par = c(0, -10, log(25), 0, 0, 0))
  t4= kmeans(data,centers=2,nstart=10)
  
  resARI3[i,1]=ARI(t1$map,lab)
  resARI3[i,2]=ARI(t2$clustering,lab)
  resARI3[i,3]=ARI(t3$optim.alloc,lab)
  resARI3[i,4]=ARI(t4$cluster,lab)
  
  sil1<-silhouette(x=t1$map,dist=dist(data))
  ressil3[i,1]=mean(sil1[,3])
  sil2<-silhouette(x=t2$clustering,dist=dist(data))
  ressil3[i,2]=mean(sil2[,3])
  sil3<-silhouette(x=t3$optim.alloc,dist=dist(data))
  ressil3[i,3]=mean(sil3[,3])
  sil4<-silhouette(x=t4$cluster,dist=dist(data))
  ressil3[i,4]=mean(sil4[,3])
}
# Compute mean ARI for 10 iterations
res3<-apply(resARI3,2,mean)
resw3<-apply(ressil3,2,mean)
end<-Sys.time()
print(end-start)

#########Scenario 3 result with sil 

######Scenario3 with silhoutte
> print(end-start)
Time difference of 2.383378 hours

resARI3
[,1]      [,2]      [,3]      [,4]
[1,] 1.0000000 0.6038894 0.6549954 0.4656051
[2,] 0.9888888 0.5781822 0.6409890 0.4887040
[3,] 0.9888888 0.3926475 0.6386620 0.3387556
[4,] 1.0000000 0.6125826 0.6552531 0.3652044
[5,] 1.0000000 0.6569798 0.6524622 0.3996637
[6,] 1.0000000 0.7601143 0.6419494 0.3856935
[7,] 1.0000000 0.3788017 0.6348474 0.2009475
[8,] 1.0000000 0.5781822 0.6310017 0.2653762
[9,] 1.0000000 0.5284443 0.6370184 0.2596652
[10,] 0.9778395 0.4356765 0.6494551 0.2110640

> res3
[1] 0.9955617 0.5525501 0.6436634 0.3380679

> ressil3
[,1]      [,2]      [,3]      [,4]
[1,] 0.4528187 0.5184311 0.3752041 0.5302900
[2,] 0.4500245 0.4975577 0.2514452 0.5045585
[3,] 0.4121952 0.5333078 0.3236789 0.5411425
[4,] 0.4349406 0.5059388 0.3406744 0.5346892
[5,] 0.4371548 0.4865674 0.3180841 0.5007277
[6,] 0.4487403 0.4807949 0.2485128 0.5088112
[7,] 0.3918422 0.5098528 0.2941541 0.5596734
[8,] 0.4290566 0.4947846 0.2477683 0.5421403
[9,] 0.4258700 0.5039845 0.2788152 0.5493957
[10,] 0.4070714 0.4874548 0.2568614 0.5275069

> resw3
[1] 0.4289714 0.5018675 0.2935198 0.5298935


#####  Scenario 4: Clusters with skewness, correlation and overlap

set.seed(123)
resARI4=matrix(0,10,4)
ressil4=matrix(0,10,4)
lab=c(rep(1,180),rep(2,180),rep(3,180))

start<-Sys.time()
for(i in 1:10){
  
  sig=matrix(0.3,15,15)
  diag(sig)=1
  c1=rdmsn(180,15,mean=rep(3,15),cov=sig,del=c(-4,-4,1.5,0 ,0, 0, 2, 2, 0, 0, 2, 0, 0,1,1))
  
  c2=rdmsn(180,15,mean=rep(1.5,15),cov=diag(15),del=c(4,4,1.5,0 ,0, 0, 2, 2, 0, 0, 2, 0, 0,3,2))
  
  sig=matrix(0.7,15,15)
  diag(sig)=1
  c3=rdmsn(180,15,mean=rep(4.5,15),cov=sig,del=c(8,4,1.5,2 ,-3, -3, 2, 5, 2, 2, -4, -7, 1,5,5))
  
  data=rbind(c1,c2,c3)
  plot(data,col=lab)
  t1= gpcm(data,G=2,mnames="VVV")
  t2= pam(data,k=2)
  t3= bclust(data,transformed.par = c(0, -10, log(25), 0, 0, 0))
  t4= kmeans(data,centers=2,nstart=10)
  
  resARI4[i,1]=ARI(t1$map,lab)
  resARI4[i,2]=ARI(t2$clustering,lab)
  resARI4[i,3]=ARI(t3$optim.alloc,lab)
  resARI4[i,4]=ARI(t4$cluster,lab)
  
  sil1<-silhouette(x=t1$map,dist=dist(data))
  ressil4[i,1]=mean(sil1[,3])
  sil2<-silhouette(x=t2$clustering,dist=dist(data))
  ressil4[i,2]=mean(sil2[,3])
  sil3<-silhouette(x=t3$optim.alloc,dist=dist(data))
  ressil4[i,3]=mean(sil3[,3])
  sil4<-silhouette(x=t4$cluster,dist=dist(data))
  ressil4[i,4]=mean(sil4[,3])
}
# Compute mean ARI for 10 iterations
res4<-apply(resARI4,2,mean)
resw4<-apply(ressil4,2,mean)
end<-Sys.time()
print(end-start)


