ptm <- proc.time()
setwd("~/GT/Spring 2016/ISYE6740 Computational Data Analysis/R code")
library("mvtnorm")
library("MASS")

#set.seed(231)

#########################################################################################################
###############################################
################################# Function Definitions 
##################################################################################################
#########################################################################################################
###############################################

#### Compute Gaussian Probability 
ComputeProb<-function(Y,K,N,mu,sigmaList){
  prob1=matrix(0,K,N)
  for (k in 1:K){ prob1[k,]=dmvnorm(Y, mu[k,], sigmaList[[k]])}
  prob1
}


#### Gamma Calculation
ComputeGamma<-function(Y,pik,sigmaList,mu,prob1){
  gam=matrix(0,N,K)
  for (i in 1:N){
    for (k in 1:K){
      gam[i,k]=prob1[k,i]*pik[k]/sum(diag(pik)%*%prob1[,i]) #Comoute Gamma
    }
  }
  gam
}

#### SIgma Calculation
ComputeSigma<-function(fit,gam,Y,mu,q,K) {
  sigmaList=vector("list",K)
  
  for (k in 1:K) {
    covmatrix=matrix(0,Col,Col)
    for (i in 1:N){
      covmatrix=covmatrix + gam[i,k]*(Y[i,]-mu[k,])%*%t(Y[i,]-mu[k,])
    }
    covmatrix=covmatrix/sum(gam[,k])
    eigVal=eigen(covmatrix,symmetric=TRUE)
    sigma2Small=1/(Col-q)*sum(eigVal$values[(q+1):Col])
    if (q>0){
      W=eigVal$vectors[,1:q]%*%diag((eigVal$values[(1:q)]-sigma2Small)^(1/2))
      SigmaBig=W%*%t(W)+sigma2Small*diag(Col)
    }
    else {
      SigmaBig=sigma2Small*diag(Col)
    }
    sigmaList[[k]]=SigmaBig
  }
  sigmaList
}

#### Pi Calculation
ComputePik<-function(fit,gam,N){
  colSums(gam)/N
}

#### Mu Calculation
ComputeMu<-function(fit,Y,gam,N,Col){
  mu=matrix(0,K,Col)
  for (k in 1:K){
    muk=0
    for (i in 1:N){
      muk = muk + gam[i,k]*Y[i,]
    }
    muk = muk/sum(gam[,k])
    mu[k,]=muk
    
  }
  mu
}

####Log-likelihood Calculation
ComputeLoglike<-function(pik,prob1){
  loglike=sum(log(colSums(diag(pik)%*%prob1))) #Compute loglikelihood
}

#### AIC Calculation
AICfun<-function(loglike,q,Col){-2*loglike+2*(Col*q+1-q*(q-1)/2)}

#### Accuracy Assessment
ComputeAccuracy<-function(gam,N,myLabel){
  Cat=vector()
  for (a in 1:N){Cat[a]=which.max(gam[a,])} #get cluster with maximum probability for each element
  missClassMat=matrix(0,10,3)
  OverallMissCount=0
  i=1
  for (label in unique(myLabel)){
    LabelCat=table(Cat[myLabel==label]) #Count occurence of clusters for each digit
    MissCat=sum(LabelCat)-max(LabelCat) #COunt the observation not in the most common category
    MissCatRate=MissCat/sum(myLabel==label) #divide by the number of occurences of each digit 
    missClassMat[i,]= c(label,MissCat,MissCatRate)
    OverallMissCount=OverallMissCount+MissCat #count overall mis categorization
    i=i+1
  }
  Overallrate=OverallMissCount/length(myLabel)
  list(missClassMat,Overallrate)
}
#########################################################################################################
#######################
####################### End of the function Definition#########
#########################################################################################################
######################

############## Main

##### Data Importation
myData=read.csv("semeion.csv",header=FALSE) # Read handwritten digits data
Y=data.matrix(myData[,1:256]) # Build data matrix with (thresholded) pixel
myLabel=apply(myData[,257:266],1,function(xx){return(which(xx=="1")-1)}) #Get Label for accuracy check
N=dim(Y)[1]
Col=dim(Y)[2]
K=10

##### Initialization
PCtested=c(0,2,4,6) #list of numbers of PC that we want to try
#PCtested=6
listMuQ=vector("list")
listAICQ=matrix(0,length(PCtested),1)
listSigmaQ=vector("list")
listloglikeQ=vector("list")

#### kmeans with several random starts for preliminary clustering
fit <- kmeans(Y, K,nstart = 10)

dev.new(width=4,height=4)
par(mai=c(0.8,0.8,0.8,0.8),mfrow=c(2,2))

######################################################################
#Loop for each Q
######################################################################
qindx=1
for (q in PCtested){
  #Initialize Gamma, Pi and Sigma
  gam=matrix(0,N,K)
  for (i in 1:N){ gam[i,fit$cluster[i]]=1 }
  pik=colSums(gam)/N
  mu=fit$centers
  sigmaList=ComputeSigma(fit,gam,Y,mu,q,K)
  Prob=ComputeProb(Y,K,N,mu,sigmaList)
  listloglike=c(ComputeLoglike(pik,Prob)) #init
  conv=1
  iter = 0
  ######################################################################
  #Iteration for a given q
  while ( iter <10){#} | conv>0.5) { #iter<1 to be consistent with the first steps
    iter = iter + 1    
    
    gam=ComputeGamma(Y,pik,sigmaList,mu,Prob)#Compute Gamma
    
    pik=ComputePik(fit,gam,N) #Update Pi
    
    mu=ComputeMu(fit,Y,gam,N,Col) #Update Mu
    
    sigmaList=ComputeSigma(fit,gam,Y,mu,q,K) #Update Sigma
    
    Prob=ComputeProb(Y,K,N,mu,sigmaList) #Calculate probability p
    
    loglike_new=ComputeLoglike(pik,Prob) #Calculate loglike for this step
    
    conv= loglike_new-listloglike[length(listloglike)] #Difference with previous loglike
    
    listloglike[iter]=loglike_new #Add loglike to list
  }
  AIC=AICfun(listloglike[length(listloglike)],q,Col)
  listAICQ[qindx]=AIC
  
  listloglikeQ[[qindx]]=listloglike
  listMuQ[[qindx]]=mu
  listSigmaQ[[qindx]]=sigmaList
  
  plot(c(1:length(listloglikeQ[[qindx]])),listloglikeQ[[qindx]],main=paste(toString(q),"PCs",sep=" "),
       ylab="Log-likelihood", xlab="Iteration")
  
  qindx=qindx+1
}
qBestIndx=which.min(listAICQ) ## Get Best q index

###### Vizualization of Clusters
dev.new(width=7,height=3.5)
par(mai=c(0.05,0.05,0.05,0.05),mfrow=c(10,6))
for(k in 1:10){
  image(t(matrix(listMuQ[[qBestIndx]][k,],byrow=TRUE,16,16)[16:1,]),col 
        =gray(seq(0,1,length=10)),axes=FALSE)
  for (r in 1:5){
    random=mvrnorm(1,matrix(listMuQ[[qBestIndx]][k,]),listSigmaQ[[qBestIndx]][[k]])
    image(t(matrix(random,byrow=TRUE,16,16)[16:1,]),col =gray(seq(0,1,length=1000)),axes=FALSE)
  }
}

#### Accuracy Assessment
Acc=ComputeAccuracy(gam,N,myLabel)
Acc[[1]] #mis-categorization rate for each label class
Acc[[2]] #overall mis-categorization rate
proc.time() - ptm