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
