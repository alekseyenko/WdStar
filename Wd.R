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

Tw2 = function(dm, f){
  if(nlevels(f) != 2) return(NULL)
  SS2 = dist.ss2(as.matrix(dm^2),f)
  ns = summary(f)
  N = sum(ns)

  SST = sum(SS2)/N
  SSW1 = SS2[1,1]/ns[1] 
  SSW2 = SS2[2,2]/ns[2]
  SSW = SSW1 + SSW2
  
  s1 = SSW1/(ns[1]-1)
  s2 = SSW2/(ns[2]-1)
  
  t.stat = (ns[1]+ns[2])/(ns[1]*ns[2])*(SST-SSW)/(s1/ns[1] + s2/ns[2])
  t.stat
}

WdS = function(dm, f){
  # This method computes Wd* statistic for distance matrix dm and factor f
  ns = table(f)
  SS2 = dist.ss2(as.matrix(dm)^2, f)
  s2 = diag(SS2)/ns/(ns-1)
  W = sum(ns/s2)
  
  idxs = apply(combn(levels(f), 2),2, function(idx) levels(f) %in% idx)
  
  Ws = sum(apply(idxs, 2, 
                 function(idx) sum(ns[idx])/prod(s2[idx]) * 
                   (sum(SS2[idx, idx])/sum(ns[idx]) - sum(diag(SS2[idx, idx])/ns[idx]))))
  k=nlevels(f)
  h = sum( (1-ns/s2/W)^2/(ns-1))
  Ws/W/(k-1)/(1+(2*(k-2)/(k^2-1))*h)
}

generic.distance.permutation.test = 
  function(test.statistic, dm, f, nrep=999, strata = NULL){
  N = length(f)
  generate.permutation=function(){
    f[sample(N)]
  }
  
  if(!is.null(strata)){
    # map elements of each strata back to their positions in the factor variable
    strata.map = order(unlist(tapply(seq_along(f), strata, identity)))
    generate.permutation=function(){
      p = unlist(tapply(f,strata,sample)) # permute within strata
      p[strata.map]
    }  
  }
  
  stats = c(test.statistic(dm, f), 
            replicate(nrep, 
                      test.statistic(dm, generate.permutation())))
  
  p.value = sum(stats>=stats[1])/(nrep+1)
  statistic = stats[1]
  list(p.value = p.value, statistic = statistic, nrep=nrep)
}

Tw2.test = function(dm, f, nrep=999, strata=NULL){
  generic.distance.permutation.test(Tw2, dm = dm, f = f, nrep = nrep, strata=strata)
}

WdS.test = function(dm, f, nrep=999, strata=NULL){
  generic.distance.permutation.test(WdS, dm = dm, f = f, nrep = nrep, strata=strata)
}

Tw2.posthoc.tests = function(dm, f, nrep=999, strata=NULL){
  dd = as.matrix(dm)
  Tw2.subset.test=function(include.levels){
    subs = f %in% include.levels
    c(include.levels, 
      table(f[subs, drop=T]), 
      Tw2.test(dd[subs, subs], f[subs,drop=T], nrep=nrep, strata=strata[subs]))
  }
  res = t(combn(levels(f), 2, Tw2.subset.test))
  colnames(res) = c("Level1", "Level2", "N1", "N2", "p.value", "tw2.stat", "nrep")
  res
}

Tw2.posthoc.1vsAll.tests = function(dm, f, nrep=999, strata=NULL){
  Tw2.subset.test=function(level){
    fs = factor(f == level)
    c(table(fs), Tw2.test(dm, fs, nrep=nrep, strata=strata))
  }
  res = t(sapply(levels(f), Tw2.subset.test))
  colnames(res) = c("N1", "N2", "p.value", "tw2.stat", "nrep")
  res
}


