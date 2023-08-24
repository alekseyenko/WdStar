is.dist = function(x) any(class(x)=='dist')

dist.sigma2 = function(dm){
  dd = as.matrix(dm)
  dd[upper.tri(dd)]=0 ##
  sum(dd^2)/nrow(dd)/(nrow(dd)-1)
}

dist.ss2 = function(dm2, f){ #dm2 is matrix of square distances; f factor
  K = sapply(levels(f), function(lev) f==lev)
  t(K)%*%dm2%*%K/2
}

dist.group.sigma2 = function(dm, f){
  diag(dist.ss2(as.matrix(dm)^2, f))/table(f)/(table(f)-1)
}

dist.cohen.d = function(dm, f){
  if(nlevels(f) != 2) return(NULL)
  SS2 = dist.ss2(as.matrix(dm^2),f)
  ns = summary(f)
  N = sum(ns)

  SST = sum(SS2)/N
  SSW = sum(diag(SS2)/ns)

  mean.diff = (sqrt((ns[1]+ns[2])/(ns[1]*ns[2])))*sqrt(SST-SSW)

  sigmas = diag(SS2)/ns/(ns-1)
  s1 = sigmas[1]
  s2 = sigmas[2]

  mean.diff/sqrt(((ns[1]-1)*s1 + (ns[2]-1)*s2)/(sum(ns)-2))
}
