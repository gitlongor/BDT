
rm(list=ls())
time0=proc.time()

nt = 1000
nv = 500
n  = nv + nt
p  = 2000
p1 = 50 # No. non-zero
vare = 10
nIter = 1000
alpha = 2e6 # Dirichlet-t concentration
verbose=10L

# Simulating data
set.seed(1234)
X = replicate(p,rnorm(n))
beta1 = rnorm(p1,0,2)
beta  = array(0,p)
beta[sample(1:p,p1,FALSE)] = beta1
g = X%*%beta
e = rnorm(n,0,sqrt(vare))
y = g+e
gt = g[1:nt]
gv = g[(nt+1):n]

# Prior
vMean = drop(var(g[1:nt])/(0.25*p))
cat("Mean effect variance ",vMean,"\n",sep="")
nu = 4.2
S  = vMean*(nu-2)/nu

# Gibbs sampling
b = array(0,p)
bSample = array(0,p)
v = array(0,p)
vSample = array(vMean,p)
mu = 0
muSample = 0
vare = 0
vareSample=0

muSample = mean(y[1:nt])
ycorr = y[1:nt] - muSample
Xt = X[1:nt,]

for(i in 1:nIter){

  vareSample = (10+crossprod(ycorr))/rchisq(1,nt+10)
  vare = vare + drop(vareSample)

  ycorr = ycorr + muSample
  muSample = rnorm(1,sum(ycorr)/nt,sqrt(vareSample/nt))
  ycorr = ycorr - muSample
  mu = mu + muSample

  for(j in 1:p){
    xj = Xt[,j]
    ycorr = ycorr + xj*bSample[j]
    rhs = crossprod(xj,ycorr)/vareSample
    lhs = crossprod(xj)/vareSample + 1.0/vSample[j]
    invLhs = 1.0/lhs
    mean = invLhs*rhs
    bSample[j]= rnorm(1,mean,sqrt(invLhs))
    ycorr = ycorr - xj*bSample[j]
  }
  b = b + bSample

  # Sample p variances using Dirichlet-gamma prior
  for(j in 1:p){
    if(sum(vSample[-j]==vSample[j])==0) vSample[j]=0
    uniquev = unique(vSample[-j])
    nc = length(uniquev)
    prob = array(0,nc)
    for(k in 1:nc){
      n_other = sum(vSample[-j]==uniquev[k])
      pr = n_other*dnorm(bSample[j],0,sqrt(uniquev[k]))/(p+alpha-1)
      prob[k] = pr
    }
    pr_new = alpha*(dt(bSample[j],df=nu)*sqrt(nu*S))/(p+alpha-1)
    prob = c(prob,pr_new)
    prob = prob/sum(prob)
    coin = array(rmultinom(1,1,prob=prob))
    newv = sum(c(uniquev,0)*coin)
    if(newv==0) newv = (S*nu+bSample[j]^2)/rchisq(1,nu+1)
    vSample[j] = newv
  }

  newUniquev = unique(vSample)
  for(k in 1:length(newUniquev)){
    ids = which(vSample==newUniquev[k])
    newvSample = (S*nu+sum(bSample[ids]^2))/rchisq(1,nu+length(ids))
    vSample[ids] = newvSample
  }
  v = v + vSample

  if(i%%verbose==0){
      time.pass=proc.time()-time0
    ghat = X%*%(b/i)
    corT = cor(gt,ghat[1:nt])
    corV = cor(gv,ghat[(nt+1):n])
    regT = coef(lm(gt~ghat[1:nt]))[2]
    regV = coef(lm(gv~ghat[(nt+1):n]))[2]
    cat("Iteration ",i,": ",round(corT,4),", ",round(regT,4)," || ",
        round(corV,4),", ",round(regV,4)," || ",length(unique(vSample)),"\t",time.pass[3],"\n",sep="")
    time0=proc.time()
  }
}

vare = vare/nIter
b = b/nIter
v = v/nIter
ghat = X%*%b
corT = cor(gt,ghat[1:nt])
corV = cor(gv,ghat[(nt+1):n])
regT = coef(lm(gt~ghat[1:nt]))[2]
regV = coef(lm(gv~ghat[(nt+1):n]))[2]
cat("Posterior: ",round(corT,4),", ",round(regT,4)," || ",
    round(corV,4),", ",round(regV,4),"\n",sep="")
length(unique(vSample))

plot(beta,b)
plot(b)
points(beta,col=2,pch=16)

hist(vSample,breaks=seq(0,max(vSample)+1,0.5))
hist(v,breaks=seq(0,max(vSample)+1,0.2))

x = seq(-5,5,0.01)
plot(dt(x,4.2),pch=16,cex=0.2)

