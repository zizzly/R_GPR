rm(list=ls())    #clear variables
library("voila")
#simulate SDE 

drift="-0.5*x"
diffusion="sqrt(2)"
dt=0.01
T=2000
x = simulate_sde(drift, diffusion, dt, T)
plot.ts(x, ylab = "x(t)", xlab = "Time t", main = "Ornstein-Uhlenbeck process")

ptm <- proc.time()
## GP parameters
m=50
uncertainty=3
targetIndex=1
inputDim=1
driftLength=1
diffLength=0.1
eps=1e-5

Tol=1e-6
diffParam=select_diffusion_parameters(x,dt,uncertainty,targetIndex,varX=NULL)

xm=matrix(seq(min(x),max(x),len=m),ncol=1)

## Create Kernel 

driftKer=sde_kernel("exp_kernel",list('amplitude'=uncertainty,'lengthScales'=driftLength),inputDim,eps)
diffKer=sde_kernel("exp_const_kernel",list('maxAmplitude'=diffParam$kernelAmplitude,'expAmplitude'=diffParam$kernelAmplitude*1e-5,'lengthScales'=diffLength),inputDim,eps)

## Fit Model
fit=sde_vi(inputDim,x,dt,xm,driftKer,diffKer,diffParam$v,20,5,Tol)

## Predict
xnew=matrix(seq(quantile(x,0.05),quantile(x,0,95),length=100),ncol=1)
driftPredict=predict(fit$drift,xnew,lognormal=FALSE,quantiles = c(0.05,0.95))
diffPredict=predict(fit$diff,xnew,lognormal=TRUE, ,quantiles = c(0.05,0.95))

driftTrue=eval(parse(text=drift), list(x=xnew))

proc.time() - ptm
## Plot Drift
plot(driftPredict, ylab = "drift f(x)", xlab = "x", main = "Ornstein-Uhlenbeck process Drift")
lines(xnew,driftTrue, col=2)
legend("topright", lty = 1, col = 1:2,
       legend = c("Estimate", "Real"), bty = "n")

plot(diffPredict, ylab = "diffusion g(x)", xlab = "x", main = "Ornstein-Uhlenbeck process Diffusion")
abline(h = eval(parse(text = diffusion), list(x = xnew)) ^ 2,
       col = 2)
legend("topright", lty = 1, col = 1:2,
       legend = c("Estimate", "Real"), bty = "n")





