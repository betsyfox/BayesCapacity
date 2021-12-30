#created by Betsy Fox (fox.119@wright.edu)
#functions for Ct Across Trials (Fox & Houpt, 2021)
#last updated 9 Dec 2021


#fiti = function(x, sh, sc){
#  options(warn=-1)
#  ll = function(p) {-sum(dweibull(x[x>0], exp(p[1]), exp(p[2]), log=TRUE))}
#  fit2 = optim(c(sh,sc), ll) 
#  sh =  exp(fit2$par[1])
#  llsc = function(p) {-sum(dweibull(x[x>0], sh, exp(p), log=TRUE))}
# fit1 = optim(sc, llsc)
#  fit_sc = sc_b(exp(fit1$par[1]),sh)
#  return(c(sh, fit_sc))
#}

sc_b = function(rsc, rsh){rsc^rsh}

#ppars = function(a, b, ns){
##  plist = list()
#  plist[[1]] = fiti(rweibull(ns, shape=a, scale=b), a, b)
#  plist[[2]] = fiti(rweibull(ns, shape=a, scale=b), a, b)
#  plist[[3]] = fiti(pmin(rweibull(ns, shape=a, scale=b), rweibull(ns, shape=a, scale=b)), a, b)
#  return(plist)
#}

#x = list of RTs to fit, sh=starting shape value it param search over, 
#sc=starting scale value to param search over, 
#ns=number of samples to simulate and fit UCIP model (unlimited assumption for prior)

Weibull.paramfit = function(cap='capacity.or', sh, sc, ns = 600){  
  plist = ppars(cap='unlimited.or', sh=sh, sc=sc, ns=ns)
  return(plist)
}



#minimum negative log likelihood function
fsc = function(x, sh, sc){
  llsc = function(p) {-sum(dweibull(x[x>0], sh, exp(p), log=TRUE))}
  res = optim(sc, llsc)
  fit_b = sc_b(exp(res$par[1]), sh)
  return(c(sh, fit_b))
}

ppars = function(cap = 'unlimited.or', sh, sc, ns){
  plist = list()
  
  if(cap=='unlimited.or'){
    plist[[1]] = fsc(rweibull(ns, sh=sh, sc=sc), sh, sc)#Channel1
    plist[[2]] = fsc(rweibull(ns, sh=sh, sc=sc), sh, sc)#Channel2
    plist[[3]] = fsc(pmin(rweibull(ns, sh=sh, sc=sc), rweibull(ns, sh=sh, sc=sc)), sh, sc)
    return(plist)
  }
}

sc_r = function(bsc, rsh){bsc^(1/rsh)}


post = function(n, s, dat, plist, a, ss, sm){    
  bv = plist[[1]][1] 
  sc1 = plist[[1]][2]; sc2 = plist[[2]][2]; sc12 = plist[[3]][2]
  
  x1 = subset(dat, RT > 0 & Channel1 > 0 & Channel2 == 0); x1s <- x1$RT^bv
  x2 = subset(dat, RT > 0 & Channel1 == 0 & Channel2 > 0); x2s <- x2$RT^bv
  x12 = subset(dat, RT > 0 & Channel1 > 0 & Channel2 > 0); x12s <- x12$RT^bv
  
  b1 = sc1*(a+1); b2 = sc2*(a+1); b12 = sc12*(a+1)
  
  pt1 = list(); pt2 = list(); pt12 = list()
  
  for(j in 1:length(ss)){
    ts <-ss[j]
    
    pt1[[j]] = swp(n, s, ts, x1$Step, x1s, a, b1, sm)
    pt2[[j]] = swp(n, s, ts, x2$Step, x2s, a, b2, sm)
    pt12[[j]] = swp(n, s, ts, x12$Step, x12s, a, b12, sm)
    
    m1 = mean(pt1[[j]]); b1 = m1*(a+1)
    m2 = mean(pt2[[j]]); b2 = m2*(a+1)
    m12 = mean(pt12[[j]]); b12 = m12*(a+1)
    
  }
  return(list(pt1, pt2, pt12))
}


swp <- function(n, s, ts, tr, d2s, a, b, sm) { 
  wt <- c(sqrt(2*pi)*sm*dnorm(ts + s - tr[tr <= ts + s], mean=0, sd=sm), rep(0, sum(tr > ts + s)))
  pa <- a + sum(wt)
  pb <- b + wt %*% d2s
  m <- (pb) / ((pa) + 1)
  ssc <- 1/rgamma(n, shape=pa, rate=pb)
  return(ssc)
}



wcs = function(par, ps){
  nst = length(ps[[1]]); ns = length(ps[[1]][[1]])
  cs = list(); t = .99
  for(i in 1:nst){
    cs[[i]] = vector(length=ns)
    for(j in 1:ns){ #assumes stationary CH1 and CH2
      H1 = (t/par[[1]][2]^(1/par[[1]][1]))^par[[1]][1]; H2 = (t/par[[2]][2]^(1/par[[1]][1]))^par[[1]][1]
      H12 = (t/ps[[3]][[i]][j]^(1/par[[3]][1]))^par[[3]][1]
      mod = H1 + H2; c = H12/mod #ratio Ct model
      cs[[i]][j] = c
    }
  }
  return(cs)
}

hdi <- function(sct){
  ns = length(sct)
  u = c(); l = c(); m = c()
  for(s in 1:ns){
    u[s] = quantile(sct[[s]], .975)
    l[s] = quantile(sct[[s]], .025)
    m[s] = mean(sct[[s]])
  }
  return(list(m, u, l))
}


bayes.capacity = function(n, s, x, plist, a, ss, sm, fun, cap, chans){
  
  ps = post(n, s, x, plist, a, ss, sm)
  sct = wcs(plist, ps)
  cts = hdi(sct)
  return(cts)
}
