rm(list=ls())    #clear variables
library("voila")
## Load Data 
xData <- read.csv("vwXBin.txt")
xData=t(xData)
xfitData <- read.csv("vwXfitBin.txt")
xfitData=t(xfitData)

#UData <- read.csv("vwU.txt")
#UData=t(UData)

dt=1 
## Variables  
m=2000               #Matlab using 18 
uncertainty=0.004  #Sigma from matlab model fit 
targetIndex=1
inputDim=1
driftLength=0.0338  #sigmaL from matlab model fit 
diffLength=0.0085    #
eps=1e-12
Tol=1e-10

## will eventually by a loop over U's
k=18; 
xU <- xData[,k]
xU <- xU[!is.na(xU)]
xfitU <- xfitData[,k]
xfitU <- xfitU[!is.na(xfitU)]

diffParam=select_diffusion_parameters(xU,1,uncertainty,targetIndex,varX=NULL)
xm=matrix(seq(min(xU,na.rm=TRUE),max(xU,na.rm=TRUE),length.out = m),ncol=1)

driftKer=sde_kernel("exp_kernel",list('amplitude'=uncertainty,'lengthScales'=driftLength),inputDim,eps)
diffKer=sde_kernel("exp_const_kernel",list('maxAmplitude'=diffParam$kernelAmplitude,'expAmplitude'=diffParam$kernelAmplitude*1e-3,'lengthScales'=diffLength),inputDim,eps)

#fit=sde_vi(inputDim,xU,dt,xm,driftKer,diffKer,diffParam$v,10,2,Tol)

#driftPredict=predict(fit$drift,xfitU,lognormal=FALSE,quantiles = c(0.05,0.95))
#diffPredict=predict(fit$diff,xfitU,lognormal=TRUE, ,quantiles = c(0.05,0.95))

