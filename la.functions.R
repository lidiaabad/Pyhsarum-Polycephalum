library(rpart)
library(cluster)

la.naintro <- function(mt, prop=0.05){
  if (is.data.frame(mt)) {
    mt=as.matrix(mt)
  }
  idx=sample(1:length(mt), prop*nrow(mt)*ncol(mt))
  mt[idx] = NA
  return(mt)
}


la.impute = function (mt, mode="median") {
  mt2=mt  #for rpart, we will see
  for (i in 1:ncol(mt)){
    if (!is.numeric(mt[,i])){
      #if it is not a numeric value go to next iteration
      next
    }
    #  create and index idx which has the positions of missing values
    idx=which(is.na(mt[,na]))
    if (mode == "median") {
      mt[idx,i]=median(mt[,i],na.rm=TRUE)
    }
    else if (mode == "mean") {
      mt[idx,i]=mean(mt[,i],na.rm=TRUE)
    }
    else if (mode == "rpart") {
      train=mt[-idx,] # no NA's. training set to create the model
      test=mt[idx,]   # rows with NA's
      form=formula(paste(colnames(mt)[i],"~.")) # the dynamically created formula 
      model=rpart(form,data=as.data.frame(train))       #we use the data we know to create the model using the formula
      mt2[idx,i]=predict(model,newdata=as.data.frame(test[,-i]))      #we use the model to predic the data of test set
    }
    else {
      stop(paste("mode", mode, "not yet implemented"))
    }
  }
  if (mode == "rpart") {
    return(mt2)
  } else {
    return(mt)
  }
}


la.norm = function (mt, mode="z") {
  if (mode %in% c("fs", "uv", "z")) {
    for (i in 1:ncol(mt)) {
      if (mode == "fs"){
        mt[,i]=(mt[,i]-min(mt[,i],na.rm=TRUE))/(max(mt[,i],na.rm=TRUE)-min(mt[,i],na.rm=TRUE))
      }
      else if (mode == "uv") {
        mt[,i]=mt[,i]/sd(mt[,i],na.rm=TRUE)
      }
      else if (mode == "z") {
        mt[,i]=(mt[,i]-mean(mt[,i],na.rm=TRUE))/sd(mt[,i],na.rm=TRUE)
      }
      else {
        stop(paste("mode", mode, "not yet implemented"))  
      }
    }
  }
  else if (mode == "mp") {
    mt <- la.impute(mt, mode="rpart")
    if (is.data.frame(mt)){
      mt=as.matrix(mt)
    }
    ncol=ncol(mt)
    nrow=nrow(mt)
    mtini=mt
    maxiter=10L #L
    for (iter in 1:maxiter){
      newiter.mt = mt
      for (i in 1:nrow){
        mt[i,] = mt[i,]-median(mt[i,])
      } 
      for (i in 1:ncol){
        mt[,i] = mt[,i]-median(mt[,i])
      }
      if (identical(mt, newiter.mt==TRUE)){break}
    }
    mt=mtini-mt
  }
  else {
    stop(paste("mode", mode, "not yet implemented"))
  }
  return(mt)
}


la.mds.plot = function (x,func="cmdscale",dist="euclidean",col='black',pch=19,names=FALSE,...) {
  plotCmd = function (cmd) {
    if (names==TRUE) {
      plot(cmd,xlab="Dim 1", ylab="Dim 2",type="n",...)
      text(cmd,labels=rownames(x),col=col)
    } else {
      plot(cmd,xlab="Dim 1", ylab="Dim 2",col=col,pch=pch,...)
    }	
  }
  if (class(x) == "dist") {
    cmd=cmdscale(x)
    plotCmd(cmd)
  } else {
    for (dm in dist) {
      D=dist(x,method=dm)
      if (func=="cmdscale") {
        cmd=cmdscale(D)
      } else {
        stop("Not ready yet")
        # MASS::sammon
      }
      plotCmd(cmd)
    }
  }
}


la.dist=function(x, method="euclidean",time=NULL) {
  if (is.data.frame(x)){
    x=as.matrix(x)
  }
  if (method %in% c("pearson","spearman","kendall")) {
    #imp to do the transposable because it does it 
    #the other way around with the cor
    D=1-(cor(t(x),method=method))/2
    rownames(D)=colnames(D)=rownames(x)   #important for 1 exercise
    dm=as.dist(D)
  }else if(method=="sts"){
    #important that ncol=nrow. we are creating a dm for the rows!
    D=matrix(data=0, nrow=nrow(x), ncol=nrow(x), dimnames=dimnames(x)[1])
    for (i in 1:nrow(x)){
      for (j in 1:i){
        D[i,j]=sqrt(sum(diff(x[i,])/diff(time)-diff(x[j,])/diff(time)^2))
      }
    }
    dm=as.dist(D)
  }else if (method %in% c("euclidean", "manhattan")) {
    dm=dist(x,method=method)
  } else {
    stop(paste("method", method, "not yet implemented."))
  }
  return(dm)
}


la.mdscluster.plot=function(x, cluster.method="complete", method="hclust"){
  Deuc=la.dist(x)
  Dspe=la.dist(x,method="spearman")
  #to do it col wise
  par(mfcol=c(2,2),mai=rep(0.4,4))
  for (m in c("euclidean", "spearman")){
    D=la.dist(x, method=m)
    plot(cmdscale(D),main=m,ylab='',type='n')
    text(cmdscale(D),rownames(x))
    if (cluster.method=="agnes"){
      #which.plots=1 only banner, =2 only dendrogram
    #  ag <- agnes(data, FALSE, metric="euclidean", FALSE, method ="single")
      
    #  plot(ag, ask = FALSE, which.plots = NULL)
      plot(agnes(D), which.plots = 2, main="agnes", ask = FALSE)
    }else if(cluster.method=="diana"){
      plot(diana(D), which.plots = 2, main="diana", ask = FALSE)
    }else{
      plot(hclust(D), main="")
    }
  }
  mtext("MDS",side = 4, adj = 1, line = 1)
  mtext("Clustering",side = 4, adj = 1, line = 1)
}


la.rand=function(x){
  for (i in 1:ncol(x)){
    x[,i]=sample(x[,i])
  }
  return(x)
}


la.multicluster = function (x, k=2, nstart=5, method="kmeans"){
  for (ki in k){
    for (ns in nstart){
      if (method=="kmeans"){
        km=kmeans(scale(x),centers=ki,nstart=ns)
        col=km$cluster
      }else if (method == "pam"){
        km=pam(x, k=ki)
        col=km$clustering
      }
      la.mds.plot(x,names=TRUE,col=col+1, main=paste("k=",ki,";nstart=",ns))
    }
  }
}


la.bestK = function (x,k=2:10,method="kmeans",dist.method="euclidean") {
  sils=c()
  for (ki in k) {
    if (method == "kmeans") {
      D=dist(x)
      cl=kmeans(x,centers=ki,nstart=5)
      sili=summary(silhouette(cl$cluster,D))$avg.width
      sils=c(sils,sili)
    }	else {
      stop("not ready yet")
    }
  }
  names(sils)=k
  return(list(bestK=k[which(sils==max(sils))],sil.widths=sils))
}	


la.bestKCap= function (x,k=2:10,method="kmeans",dist.method="euclidean") {
  fun=list()
  fun$pam<-function(x,k){
    return(list(cluster=pam(x,k,cluster.only=TRUE)))
  }
  fun$kmeans<-function(x,k){
    return(list(cluster=kmeans(x,centers=k,nstart=5)$cluster))
  }
  #bootstraping hasta 200
  gp=clusGap(x,FUN=fun[[method]],K.max=max(k),B=200)
  #bestk for maximum gap
  bestK=which(gp$Tab[,3]==max(gp$Tab[,3]))
  print(gp$Tab)
  return(bestK)
}


la.pclust = function (x,k=0,method="kmeans", sampsize=nrow(data)/10,metric="euclidean",...) {
  result=list(clustering=NULL,k=k,method=method,metric=metric,sqilinfo=NULL)
  #methods=c("kmeans","pam","clara", "fanny","fastpam", "fastclara")        
  if (method %in% c("fastpam","fastclara")) {
    if (!require("fastkmedoids")) {
      stop("Error: methods 'fastpam' and 'fastclara' needs package fastkmedoids installed!")
    } 
  }     
  if (k == 0) {
    lb=la.bestK(x, k=2:10, clustering=method, distance=metric)
    k=lb$bestK
  }
  D=dist(x, method=metric, ...) 
  if (method=="kmeans"){
    if (metric != "euclidean"){
      stop("kmeand can only use euclidean distances")
    }
    cl=kmeans(D, centers=k, nstart=5, ...)
    clustering=cl$cluster
    result$metric=metric
    silinfo=silhouette(clustering, D)
  }else if (method=="pam"){
    cl=cluster::pam(D, k=k, diss=TRUE, ...)
    clustering=cl$cluster
    silinfo=cl$silinfo
  }else if (method=="clara"){
    cl=cluster::clara(x, k=k, samples=sampsize, ...)
    clustering=cl$cluster
    silinfo=cl$silinfo
  }else if (method=="fanny"){
    cl=fanny(x, k=k, ...)
    clustering=cl$cluster
    silinfo=cl$silinfo
  }else if (method=="fastpam"){
    cl=fastpam(D, nrow(x), k=k)
    clustering=cl@assignment
    silinfo=silhouette(cl@assignment, D)
  }else if (method=="fastclara"){
    cl=fastclara(D, nrow(x), k=k)
    clustering=cl@assignment
    silinfo=silhouette(cl@assignment, D)
  }else {
    stop("not implemented yet")
  }
  result$clustering=clustering
  result$k=k
  result$silinfo=silinfo
  result$diss=D
  result$method=method
  return(result)
}


la.clusterSample=function(k=2, n=50, sd=0.4){
  # The function returns a matrix with the x and y
  #coordinates of the data. The cluster centers
  #should follow the point layout of a dice.
  #Select the initial x and y positions for the 
  #cluster centers ranges carefully to have good 
  #cluster separation with sd=1 the standard distance
  #in x and y. Space between the centers should be around 
  #3 units providing quite good separation you can hardcode the centers
  #for the clusters 2 and 4 in your function
  #where the horizontal and vertical distances for 
  #the 4 solution should be around 3 in case of
  #the 4 cluster solution.
  if (length(n)==1){
    M=matrix(0,nrow=k*n,ncol=2)
    colnames(M)=c("x","y")
    if(k==2){
      rn1x=rnorm(n,mean=1,sd=sd)
      rn1y=rnorm(n,mean=1,sd=sd)
      rn2x=rnorm(n,mean=3,sd=sd)
      rn2y=rnorm(n,mean=3,sd=sd)
      M[,1]=c(rn1x,rn2x)
      M[,2]=c(rn1y,rn2y)
    } else if(k==4){
      #por que esos valores de la media??
      rn1x=rnorm(n,mean=1,sd=sd)
      rn1y=rnorm(n,mean=1,sd=sd)
      rn2x=rnorm(n,mean=3,sd=sd)
      rn2y=rnorm(n,mean=3,sd=sd)
      rn3x=rnorm(n,mean=1,sd=sd)
      rn3y=rnorm(n,mean=3,sd=sd)
      rn4x=rnorm(n,mean=3,sd=sd)
      rn4y=rnorm(n,mean=1,sd=sd)
      M[,1]=c(rn1x,rn2x,rn3x,rn4x)
      M[,2]=c(rn1y,rn2y,rn3y,rn4y)
    }
  } else{
    stop("Error:only a single n can be currently given!")
  }
  return(M)
}


la.clusterCompare = function (x, y) {
  #result=list(rand=rand, rand.rand=rand.rand, jaccard=jaccard, rand.jaccard=rand.jaccard, rao=rao, rand.rao=rand.rao, tanimoto=tanimoto, rand.tanimoto=rand.tanimoto, kulczynski2=kulczynski2, rand.kulczynski2=rand.kulczynski2, ochiai=ochiai, rand.ochiai=rand.ochiai )
  la.indexcount= function (x, y){
    SS=0
    SD=0
    DS=0
    DD=0
    countlist=list(SS, SD, DS, DD)
    for (i in 1:(length(x)-1)){
      for (j in (i+1): length(x)){
        if (x[i]== x[j] & y[i] == y[j]){
          SS=SS+1 
        } else if (x[i]== x[j] & y[i] != y[j]){
          SD=SD+1 
        } else if (x[i]!= x[j] & y[i] == y[j]){
          DS=DS+1 
        } else if (x[i]!= x[j] & y[i] != y[j]){
          DD=DD+1 
        }
      }
    }
    return(countlist)
  }
  rrand=c()
  rjaccard=c()
  rrao=c()
  rtanimoto=c()
  rkulczynski=c()
  rochiai=c()
  for (k in 1:10){
    xr=la.rand(x)
    yr=la.rand(y)
    crlist=la.indexcount(xr, yr)
    SS=crlist[1]
    SD=crlist[2]
    DS=crlist[3]
    DD=crlist[4]
    randi.rand=(SS+DD)/(SS+DS+SD+DD)
    randi.jaccard=SS/(SS+DS+SD)
    randi.rao=SS/(SS+DS+SD+DD)
    randi.tanimoto=(SS+DD)/(SS+2*SD+2*DS+DD)
    randi.kulczynski2=(SS/(SS+SD) + SS/(SS+DS))/2
    randi.ochiai=sqrt(SS/(SS+SD)*SS/(SS+DS))
    rrand=append(rrand, randi.rand)
    rjaccard=c(rjaccard, randi.jaccard)
    rrao=c(rrao, randi.rao)
    rtanimoto=c(rtanimoto, randi.tanimoto)
    rkulczynski=c(rkulczynski, randi.kulczynski2)
    rochiai=c(rochiai, randi.ochiai)
  }
  rand.rand=mean(rrand)
  rand.jaccard=mean(rjaccard)
  rand.rao=mean(rrao)
  rand.tanimoto=mean(rtanimoto)
  rand.kulczynski2=mean(rkulczynski)
  rand.ochiai=mean(rochiai)
  
  clist=la.indexcount(x, y)
  SS=clist[1]
  SD=clist[2]
  DS=clist[3]
  DD=clist[4]
  rand=(SS+DD)/(SS+DS+SD+DD)
  jaccard=SS/(SS+DS+SD)
  rao=SS/(SS+DS+SD+DD)
  tanimoto=(SS+DD)/(SS+2*SD+2*DS+DD)
  kulczynski2=(SS/(SS+SD) + SS/(SS+DS))/2
  ochiai=sqrt(SS/(SS+SD)*SS/(SS+DS))
  result=list(rand=rand, rand.rand=rand.rand, jaccard=jaccard, rand.jaccard=rand.jaccard, rao=rao, rand.rao=rand.rao, tanimoto=tanimoto, rand.tanimoto=rand.tanimoto, kulczynski2=kulczynski2, rand.kulczynski2=rand.kulczynski2, ochiai=ochiai, rand.ochiai=rand.ochiai )
  return(result)
}



la.NetworkCompare= function (x, y) {
  #result=list(rand=rand, rand.rand=rand.rand, jaccard=jaccard, rand.jaccard=rand.jaccard, rao=rao, rand.rao=rand.rao, tanimoto=tanimoto, rand.tanimoto=rand.tanimoto, kulczynski2=kulczynski2, rand.kulczynski2=rand.kulczynski2, ochiai=ochiai, rand.ochiai=rand.ochiai )
  la.indexcount= function (x, y){
    SS=0
    SD=0
    DS=0
    DD=0
    for (i in 1:(nrow(x)-1)){
      for (j in (i+1):ncol(x)){
        if (x[i,j]==1 & y[i,j]==1){
          SS=SS+1 
        } else if (x[i,j]==1 & y[i,j]==0){
          SD=SD+1 
        } else if (x[i,j]==0 & y[i,j]==1){
          DS=DS+1 
        } else if (x[i,j]==0 & y[i,j]==0){
          DD=DD+1 
        }
      }
    }
    countlist=list(SS=SS, SD=SD, DS=DS, DD=DD)
    return(countlist)
  }
  
  clist=la.indexcount(x, y)
  SS=clist[1]
  SD=clist[2]
  DS=clist[3]
  DD=clist[4]
  SS=as.integer(SS)
  SD=as.integer(SD)
  DS=as.integer(DS)
  DD=as.integer(DD)
  
  rand=(SS+DD)/(SS+DS+SD+DD)
  jaccard=SS/(SS+DS+SD)
  rao=SS/(SS+DS+SD+DD)
  tanimoto=(SS+DD)/(SS+2*SD+2*DS+DD)
  kulczynski2=(SS/(SS+SD) + SS/(SS+DS))/2
  ochiai=sqrt(SS/(SS+SD)*SS/(SS+DS))
  
  
  rrand=c()
  rjaccard=c()
  rrao=c()
  rtanimoto=c()
  rkulczynski=c()
  rochiai=c()
  for (k in 1:10){
    xr=la.rand(x)
    yr=la.rand(y)
    
    crlist=la.indexcount(xr, yr)
    SS=crlist[1]
    SD=crlist[2]
    DS=crlist[3]
    DD=crlist[4]
    SS=as.integer(SS)
    SD=as.integer(SD)
    DS=as.integer(DS)
    DD=as.integer(DD)
    
    randi.rand=(SS+DD)/(SS+DS+SD+DD)
    randi.jaccard=SS/(SS+DS+SD)
    randi.rao=SS/(SS+DS+SD+DD)
    randi.tanimoto=(SS+DD)/(SS+2*SD+2*DS+DD)
    randi.kulczynski2=(SS/(SS+SD) + SS/(SS+DS))/2
    randi.ochiai=sqrt(SS/(SS+SD)*SS/(SS+DS))
    rrand=append(rrand, randi.rand)
    rjaccard=c(rjaccard, randi.jaccard)
    rrao=c(rrao, randi.rao)
    rtanimoto=c(rtanimoto, randi.tanimoto)
    rkulczynski=c(rkulczynski, randi.kulczynski2)
    rochiai=c(rochiai, randi.ochiai)
  }
  rand.rand=mean(rrand)
  rand.jaccard=mean(rjaccard)
  rand.rao=mean(rrao)
  rand.tanimoto=mean(rtanimoto)
  rand.kulczynski2=mean(rkulczynski)
  rand.ochiai=mean(rochiai)
  

  
  result=list(rand=rand, rand.rand=rand.rand, jaccard=jaccard, rand.jaccard=rand.jaccard, rao=rao, rand.rao=rand.rao, tanimoto=tanimoto, rand.tanimoto=rand.tanimoto, kulczynski2=kulczynski2, rand.kulczynski2=rand.kulczynski2, ochiai=ochiai, rand.ochiai=rand.ochiai )
  return(result)
}  

