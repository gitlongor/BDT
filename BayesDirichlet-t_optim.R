
rm(list=ls())
time0=proc.time()

nt = 1000L
nv = 500L
n  = nv + nt
p  = 2000L
p1 = 50L # No. non-zero
vare = 10
nIter = 1000L
alpha = 100 # Dirichlet-t concentration
invRootAlpha=if(alpha==0) 1 else 1/sqrt(alpha)
rootAlpha =if(is.finite(alpha)) sqrt(alpha) else 1
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
sqrt.S.nu=sqrt(S*nu)

# Gibbs sampling
b = array(0,p)
bSample = array(0,p)
if(TRUE){ ## starting from one cluster
	#v = array(0,p)
	#vSample = array(vMean,p);	
	i2class=rep(1L,p)		# i2class[i] is the cluster label for marker i
	n.class=c(p, rep(0L,p-1L))  # markers per cluster
	nclasses=1L		# number of clusters
	vClass = c(vMean, rep(0, p-1))  # var for each cluster
	class2i=c(list(1:p),rep(list(integer(0L)),p-1L)) # cluster membership
	classList=1:p		# first nclasses elements hold cluster labels
	invClassList=c(1L, rep(0L, p-1L))  # classList[inClassList[i]]=i
}else{ ## starting from p clusters
	i2class=1:p		# i2class[i] is the cluster label for marker i
	n.class=rep(1L, p) # markers per cluster
	nclasses=p		# number of clusters
	vClass = runif(p, vMean/10, vMean*10) # var for each cluster
	class2i=as.list(1:p) # cluster membership
	classList=1:p		# first nclasses elements hold cluster labels
	invClassList=c(1L, rep(0L, p-1L))  # classList[inClassList[i]]=i
}
probs = rep(0, p) ## alloc memory
sqrt.vClass=sqrt(vClass)

mu = 0
muSample = 0
vare = 0
vareSample=0

muSample = mean(y[1:nt])
ycorr = y[1:nt] - muSample
Xt = X[1:nt,]
xmean=colMeans(Xt)
Xt=Xt-rep(xmean,each=nt)
X[1:nv+nt,]=X[1:nv+nt,]-rep(xmean, each=nv)


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
    lhs = crossprod(xj)/vareSample + 1.0/vClass[i2class[j]]
    invLhs = 1.0/lhs
    mean = invLhs*rhs
    bSample[j]= rnorm(1,mean,sqrt(invLhs))
    ycorr = ycorr - xj*bSample[j]
  }
  b = b + bSample

  # Sample p variances using Dirichlet-gamma prior
  for(j in 1:p){
	this.class=i2class[j]
	if(n.class[this.class]==1L){ ## switch
	    classList[invClassList[this.class]]=classList[nclasses]
	    invClassList[classList[nclasses]]=invClassList[this.class]
		
	    classList[nclasses]=this.class
		invClassList[this.class]=nclasses
		
		n.class[this.class]=0L
		nclasses=nclasses-1L
	}
	seq.nc=seq_len(nclasses)
	other.classes=classList[seq.nc]
	n.other=n.class[other.classes]
	tmp=invClassList[this.class]
	if(tmp<=nclasses) n.other[tmp]=n.other[tmp]-1L
	
	probs[seq.nc] = invRootAlpha * n.other * dnorm(bSample[j], 0, sqrt.vClass[other.classes])
	#probs[seq.nc] = n.other * dnorms[j, other.classes]
	probs[nclasses+1L]=rootAlpha * dt(bSample[j], nu)*sqrt.S.nu
	new.class=classList[sample.int(nclasses+1L,1L,prob=probs[seq_len(nclasses+1L)])]
	if(n.class[this.class]>0L) n.class[this.class]=n.class[this.class]-1L
	class2i[[this.class]]=class2i[[this.class]][ class2i[[this.class]]!=j ]

	if(new.class==classList[nclasses+1L]){
		nclasses=nclasses+1L
		invClassList[new.class]=nclasses
		n.class[new.class]=1L
	}else  n.class[new.class]=n.class[new.class]+1L
	class2i[[new.class]]=c(class2i[[new.class]], j)
	i2class[j]=new.class
  }
  
  idx=classList[seq_len(nclasses)]
  vClass[idx]=1 / rchisq(nclasses, nu + n.class[idx])
  vClass[-idx]=0
  for(k in idx)
    #vSample[class2i[[k]]] =  (
	vClass[k] = (S*nu+sum(bSample[class2i[[k]]]^2)) * vClass[k]
	#)
  #v = v + vSample
  sqrt.vClass=sqrt(vClass)
  
  if(i<=10 || i%%verbose==0){
      time.pass=proc.time()-time0
      ghat = X%*%(b/i)
      corT = cor(gt,ghat[1:nt])
      corV = cor(gv,ghat[(nt+1):n])
      regT = coef(lm(gt~ghat[1:nt]))[2]
      regV = coef(lm(gv~ghat[(nt+1):n]))[2]
      cat("Iteration ",i,": ",round(corT,4),", ",round(regT,4)," || ",
          round(corV,4),", ",round(regV,4)," || ",nclasses,"\t",time.pass[3],"\n",sep="")
      time0=proc.time()
  }
  
}

vare = vare/nIter
b = b/nIter
#v = v/nIter
ghat = X%*%b
corT = cor(gt,ghat[1:nt])
corV = cor(gv,ghat[(nt+1):n])
regT = coef(lm(gt~ghat[1:nt]))[2]
regV = coef(lm(gv~ghat[(nt+1):n]))[2]
cat("Posterior: ",round(corT,4),", ",round(regT,4)," || ",
    round(corV,4),", ",round(regV,4),"\n",sep="")
nclasses

plot(beta,b)
plot(b)
points(beta,col=2,pch=16)

#hist(vSample,breaks=seq(0,max(vSample)+1,0.5))
#hist(v,breaks=seq(0,max(vSample)+1,0.2))

x = seq(-5,5,0.01)
plot(dt(x,4.2),pch=16,cex=0.2)

