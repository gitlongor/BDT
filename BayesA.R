
bayesA = function(jobStr,seed){

  set.seed(seed)

  # Parameters
  chainLength = 2100
  burnIn = 100
  nua = 4.2
  nue = 4.2

  job = jobStr
  d = read.table(paste(".//",job,"//",job,".mrk",sep=""),header=TRUE)
  dy = read.table(paste(".//",job,"//",job,".phen",sep=""),header=TRUE)
  dt = read.table(paste(".//",job,"//",job,".tbv",sep=""),header=TRUE)

  Z = as.matrix(d[,-1])
  N = nrow(Z)
  nt = N/2
  cat("Number observation: ",N,"\n",sep="")
  cat("Number training: ",nt,"\n",sep="")

  y    = dy$phenotype[1:nt]
  a    = dt$phenotype[1:nt]
  aV   = dt$phenotype[(nt+1):N]

  # inital values
  Z = (Z+10)/10
  freq = apply(Z,2,mean)
  Z = sweep(Z,2,freq,"-")
  freq = freq/2
#  del  = which(freq>0.95 || freq<0.05)
#  Z = Z[,-del]
#  freq = freq[-del]
  nmarkers = ncol(Z)
  cat("Number markers: ",nmarkers,"\n",sep="")

  vara    = var(a)
  mean2pq = mean(2*freq*(1-freq))
  varEffects  = vara/(nmarkers*mean2pq)
  scalea      = varEffects*(nua-2)/nua
  scalee      = (var(y)-vara)*(nue-2)/nue

  mu = 0
  muSample = 0
  vare = 0
  vareSample = 0
  u = array(0.0,nmarkers)
  uSample = u
  muSample = mean(y)
  varj = array(0.0,nmarkers)
  varjSample = varj

  ycorr = y - muSample # - Z[1:nt,]%*%uSample # adjust y

  # mcmc sampling
  for (iter in 1:chainLength){

    # sample vare
    vareSample = ( nue*scalee + crossprod(ycorr) )/rchisq(1,nt+nue)
    if(iter>burnIn) vare = vare + vareSample

    # sample intercept
    ycorr = ycorr + muSample
    rhs    = sum(ycorr)/vareSample
    invLhs = 1.0/(nt/vareSample)
    mean = rhs*invLhs
    muSample = rnorm(1,mean,sqrt(invLhs))
    ycorr = ycorr - muSample
    if(iter>burnIn) mu = mu + muSample

    # sample variance for each locus
    varjSample = (scalea*nua + uSample^2)/rchisq(nmarkers,nua+1)
    if(iter>burnIn) varj = varj + varjSample

    # sample effect for each locus
    for(locus in 1:nmarkers){
      z = Z[1:nt,locus]
      ycorr = ycorr + z*uSample[locus]
      rhs = crossprod(z,ycorr)/vareSample
      lhs = crossprod(z)/vareSample + 1.0/varjSample[locus]
      invLhs = 1.0/lhs
      mean = invLhs*rhs
      uSample[locus]= rnorm(1,mean,sqrt(invLhs))
      ycorr = ycorr - z*uSample[locus]
      if(iter>burnIn) u[locus] = u[locus] + uSample[locus]
    }

    if(iter%%10==0 && iter>burnIn){
      meanu = u/(iter-burnIn)
      aHat  = Z %*% meanu

      gebvT = aHat[1:nt]
      gebvV = aHat[(nt+1):N]
      regT = coef(lm( a~gebvT))[2]
      regV = coef(lm(aV~gebvV))[2]
      cat("Iteration ",iter,": ",sep="")
      cat(" ",cor(a,gebvT),", ",regT,", || ",cor(aV,gebvV),", ",regV,".\n",sep="")
    }
  }

  vare = vare/(chainLength-burnIn)
  u = u/(chainLength-burnIn)
  aHat  = Z %*% u
  post_cor = cor(aV,aHat[(nt+1):N])

  write.table(post_cor,paste(job,"_BayesA.res",sep=""),quote=F,row.names=F,col.names=F,append=F)
}



