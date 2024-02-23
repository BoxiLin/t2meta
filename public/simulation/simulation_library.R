library(dplyr)
library(utils)
library(parallel)
library(pbmcapply)
set.seed(123456)

# Data generator
generater <- function(maf, n.f, n.m, beta.f, beta.m, 
                      r = 1,
                      err_dist = "norm", 
                      model = "additive",
                      response = "linear") {
  sigma = 1
  sigma.f = 1*r
  beta.s = 1
  
  g.f = rep(c(2,1,0), c(round(n.f*c(maf^2, 2*maf*(1-maf))),n.f-sum(round(n.f*c(maf^2, 2*maf*(1-maf))))))
  g.m = rep(c(2,1,0), c(round(n.m*c(maf^2, 2*maf*(1-maf))),n.m-sum(round(n.m*c(maf^2, 2*maf*(1-maf))))))
  
  sex.f = rep(1, n.f)
  sex.m = rep(0, n.m)
  
  # Error distribution 
  if (err_dist=="norm") {
    e.f = rnorm(n = n.f, sd = sigma.f)
    e.m = rnorm(n = n.m, sd = sigma)
  } else if (err_dist=="t") {
    e.f = rt(n = n.f, df = 4)
    e.m = rt(n = n.m, df = 4)
  } else if(err_dist == "chisq") {
    e.f = rchisq(n = n.f, df = 4)
    e.m = rchisq(n = n.m, df = 4)
  } else {
    message("Error distribution should be one of norm, t, or chisq")
    break
  }
  
  # Genetic coding
  if (model == "additive") {
    G.f <- g.f
    G.m <- g.m
  } else if (model == "dominant") {
    G.f <- replace(g.f, g.f==2,1)
    G.m <- replace(g.m, g.m==2,1)
  } else {
    message("Genetic model should be one of additive and dominant")
    break
  }
  
  if (response == "linear") {
  y.f = beta.f*G.f + sex.f*beta.s + e.f
  y.m = beta.m*G.m + e.m
  } else if (response == "binary") {
    p.f = 1/(1+exp(-(beta.f*G.f + sex.f*beta.s)))
    p.m = 1/(1+exp(-(beta.m*G.m)))
    y.f = rbinom(n.f,1,p.f)
    y.m = rbinom(n.m,1,p.m)
  } else {
    message("Response should be one of binary and linear")
    break
  }
  
  X = data.frame(G = c(G.f, G.m), sex = c(sex.f, sex.m))
  dat <- data.frame(X, Y = c(y.f, y.m))
  return(dat)
}


# 1. TSG.f: Female only
TSG.f <- function(dat) {
  f <- lm(Y~G+sex, data = dplyr::filter(dat, sex == 1))
  Z.f <- summary(f)$coefficients["G", "t value"]
  p.f <- summary(f)$coefficients["G", "Pr(>|t|)"]
#  return(c(Z.TSG.f = Z.f, p.TSG.f = p.f, rej.TSG.f = (p.f<alpha)))
  return(c(z.TSG.f = Z.f, p.TSG.f = p.f))
}

# TSG.f(dat)

# 2. TSG.min: Minimum p-value
TSG.min <- function(dat) {
  f <- lm(Y~G, data = dplyr::filter(dat, sex == 1))
  m <- lm(Y~G, data = dplyr::filter(dat, sex == 0))
  Z.f <- summary(f)$coefficients["G", "t value"]
  p.f <- summary(f)$coefficients["G", "Pr(>|t|)"]
  Z.m <- summary(m)$coefficients["G", "t value"]
  p.m <- summary(m)$coefficients["G", "Pr(>|t|)"]
#  return(c(Z.TSG.min = c(Z.f, Z.m)[which.min(c(p.f, p.m))],
#  p.TSG.min = min(p.f,p.m), rej.TSG.min = (min(p.f,p.m)<1-sqrt(1-alpha))))
  return(c(p.TSG.min = (min(p.f,p.m))))
}

# TSG.min(dat)

# 3. TSG.1: 
TSG.meta <- function(dat) {
  f <- lm(Y~G, data = dplyr::filter(dat, sex == 1))
  m <- lm(Y~G, data = dplyr::filter(dat, sex == 0))
  z.f <- summary(f)$coefficients["G", "t value"]
  z.m <- summary(m)$coefficients["G", "t value"]
  
  x.f <- summary(f)$coefficients["G", "Estimate"]
  x.m <- summary(m)$coefficients["G", "Estimate"]
  w.f = (1/(summary(f)$coefficients["G", "Std. Error"]))^2*1000/998
  w.m = (1/(summary(m)$coefficients["G", "Std. Error"]))^2*1000/998
  
  z.TSG.1 = sqrt(w.f/(w.f+w.m))*z.f+sqrt(w.m/(w.f+w.m))*z.m
  p.TSG.1 = 2*pnorm(abs(z.TSG.1), lower.tail = F)
  
  q.TSG.q = z.f^2+z.m^2
  p.TSG.q = pchisq(df = 2, q = q.TSG.q, lower.tail = FALSE)
  return(c( z.TSG.1 =  z.TSG.1, p.TSG.1 = p.TSG.1, 
            q.TSG.q =  q.TSG.q, p.TSG.q = p.TSG.q))
}

# TSG.meta(dat)  
 
# 5. TI.joint 

TI.2df <- function(dat){
  sig.f <- summary(lm(Y~G, data = dplyr::filter(dat, sex == 1)))$sigma
  sig.m <- summary(lm(Y~G, data = dplyr::filter(dat, sex == 0)))$sigma
  wt = 1/(ifelse(dat$sex==1, sig.f^2, sig.m^2))
  m1 = summary(lm(Y~G*sex, data = dat, weights = wt))
  V = m1$cov[c("G","G:sex"),c("G","G:sex")]
  beta = m1$coefficients[c("G","G:sex"),"Estimate"]
  wstat = t(beta)%*%solve(V)%*%beta
  p.TI.joint = pchisq(q = wstat, df = 2, lower.tail = FALSE)
  z.TI.joint = wstat
  return(c(z.TI.joint = z.TI.joint, p.TI.joint = p.TI.joint))
}

# TI.2df(dat) 

# 6. T.main 

T.main <- function(dat) {
  sig.f <- summary(lm(Y~G, data = dplyr::filter(dat, sex == 1)))$sigma
  sig.m <- summary(lm(Y~G, data = dplyr::filter(dat, sex == 0)))$sigma
  wt = 1/(ifelse(dat$sex==1, sig.f^2, sig.m^2))
  m1 = summary(lm(Y~G+sex, data = dat, weights = wt))
  p.T.main = m1$coefficients["G","Pr(>|t|)"]
  z.T.main = m1$coefficients["G","t value"]
  return(c(z.T.main = z.T.main, p.T.main = p.T.main))
}

#T.main(dat) 

param_generator <- function(scenario, maf = c(0.1, 0.25, 0.01), beta.f, beta.m,
                            k = c(0.5,1,1.5,2), n.f = 1000, r = c(0.5,1,2,5), err_dist = "norm",
                            model = "additive", response = "linear") {
  param = data.frame(scenario = scenario, 
                     expand.grid(maf = maf, n.f = n.f,k=k, r = r, 
                                 beta.f = beta.f,beta.m = beta.m,
                                 err_dist = err_dist, model = model,
                                 response = response))
  param$n.m = param$n.f*param$k
  dplyr::select(param,c("scenario", "maf", "n.f","n.m","k","r","beta.f","beta.m",
                        "err_dist", "model", "response"))
}

