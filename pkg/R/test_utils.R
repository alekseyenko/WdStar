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
