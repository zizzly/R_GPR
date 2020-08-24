rm(list=ls())    #clear variables
library("voila")
## Load Data 

description <- "Aug24"

filename <- paste0("vwXBin",description,".txt")
xData <- read.csv(filename)
xData=t(xData)
filename <- paste0("vwXfit",description,".txt")
xfit <- read.csv(filename)
xfit=t(xfit)

#UData <- read.csv("vwU.txt")
#UData=t(UData)
Nu <- 15  

Nx <- 1000

## Initialize
DriftFit <- matrix(nrow=Nx, ncol=Nu)
DiffFit <- matrix(nrow=Nx, ncol=Nu)
KDriftSave <- matrix(nrow=2,ncol=Nu)
KDiffSave <- matrix(nrow=2,ncol=Nu)
BasisDiffSave <- matrix(nrow=2,ncol=Nu)

dt=1
## Variables  
m=50               #Matlab using 2000
uncertainty=0.001  #Sigma from matlab model fit 
targetIndex=1
inputDim=1
driftLength=0.03 #sigmaL from matlab model fit 
diffLength=0.01    #
eps=1e-11
Tol=1e-10

MaxIt <-20 #number of itterationns 
HPit <-  5     #number of hyperparameter itts

## will eventually by a loop over U's

# Start the clock!
ptm <- proc.time()

for (k in 1:Nu)
{xU <- xData[,k]
xU <- xU[!is.na(xU)]
indS= which.min(abs(xfit-quantile(xU,0.01)))
indE= which.min(abs(xfit-quantile(xU,0.99)))
xfitU = xfit[indS:indE]
xfitU <- data.matrix(xfitU)

diffParam=select_diffusion_parameters(xU,1,uncertainty,targetIndex,varX=NULL)
xm=matrix(seq(min(xU,na.rm=TRUE),max(xU,na.rm=TRUE),length.out = m),ncol=1)

driftKer=sde_kernel("exp_kernel",list('amplitude'=uncertainty,'lengthScales'=driftLength),inputDim,eps)
diffKer=sde_kernel("exp_const_kernel",list('maxAmplitude'=diffParam$kernelAmplitude,'expAmplitude'=diffParam$kernelAmplitude*1e-3,'lengthScales'=diffLength),inputDim,eps)

fit=sde_vi(inputDim,xU,dt,xm,driftKer,diffKer,diffParam$v,MaxIt,HPit,Tol)

driftPredict=predict(fit$drift,xfitU,lognormal=FALSE, ,quantiles = c(0.05,0.95))
diffPredict=predict(fit$diff,xfitU,lognormal=TRUE, ,quantiles = c(0.05,0.95))

D1U=driftPredict$mean
D2U=diffPredict$mean

D1U=data.matrix(D1U)
D2U=data.matrix(D2U)
nu <- length(D1U)

DriftFit[indS:indE,k] = D1U
DiffFit[indS:indE,k] =  D2U

Basis1=diffParam$v 
Basis2=diffParam$kernelAmplitude


KDriftSave[1,k]=fit$drift$gpKernel$amplitude
KDriftSave[2,k]=fit$driftHyperparams
KDiffSave[,k]=fit$diffHyperparams

BasisDiffSave[1,k]=Basis1
BasisDiffSave[2,k]=Basis2
uncertainty=KDriftSave[1,k] #Sigma from matlab model fit 

driftLength=KDriftSave[2,k] #sigmaL from matlab model fit 
diffLength=KDiffSave[2,k] 
write.csv(DriftFit,'DriftFit.csv')
write.csv(DiffFit,'DiffFit.csv')
write.csv(KDriftSave,'KDrift.csv')
write.csv(KDiffSave,'KDiff.csv')
write.csv(BasisDiffSave,'BasisDiff.csv')

}   #

# Stop the clock
proc.time() - ptm

write.csv(DriftFit,'DriftFit.csv')
write.csv(DiffFit,'DiffFit.csv')
write.csv(KDriftSave,'KDrift.csv')
write.csv(KDiffSave,'KDiff.csv')
write.csv(BasisDiffSave,'BasisDiff.csv')
## Plot Drift
plot(driftPredict, ylab = "drift f(x)", xlab = "x", main = "Drift")
lines(xfitU,driftTrue, col=2)


plot(diffPredict, ylab = "diffusion g(x)", xlab = "x", main = "Diffusion")
