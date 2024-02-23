source("simulation_library.R")

# set.seed(1000)
set.seed(20220428)

# read parameter
args = commandArgs(trailingOnly=TRUE)
i = 8
i = as.numeric(args[1])
print(paste("Parameter null setting i=",i))

param1 = param_generator(scenario = 1, n.f = 5000, k=1, r = c(0.5,1,2,5),
                         maf= c(0.05, 0.1, 0.25, 0.5), 
                         beta.f  = seq(0, 0.3, 0.03), beta.m = seq(0, 0.3, 0.03),
                         err_dist = c("norm"), 
                         model =  c("additive"),
                         response = c("linear")) %>%
  dplyr::filter(beta.m == beta.f)


param2 = param_generator(scenario = 2, n.f = 5000, k=1, r = c(0.5,1,2,5),
                         maf= c(0.05, 0.1, 0.25, 0.5), 
                         beta.f  = seq(0, 0.3, 0.03), beta.m = 0,
                         err_dist = c("norm"), 
                         model =  c("additive"),
                         response = c("linear")) 

param3 = param_generator(scenario = 3, n.f = 5000, k=1, r = c(0.5,1,2,5),
                         maf= c(0.05, 0.1, 0.25, 0.5), 
                         beta.f  = 0.15, beta.m = seq(0, 0.3, 0.06),
                         err_dist = c("norm"), 
                         model =  c("additive"),
                         response = c("linear"))

param4 = param_generator(scenario = 4, n.f = 5000, k=1, r = c(0.5,1,2,5),
                         maf= c(0.05, 0.1, 0.25, 0.5), 
                         beta.f  = 0.15, beta.m = -seq(0, 0.6, 0.06),
                         err_dist = c("norm"), 
                         model =  c("additive"),
                         response = c("linear"))

param_alt <- rbind(param1, param2, param3, param4)


# a = 0.05
a = 5e-8
# Changed from 50->100
# M = 1000
M = 10000
R = 1

mc.cores = min(dim(param_alt)[1], detectCores()-2,36)

print(paste("Number of cores:",mc.cores))

parameter = param_alt[i,]
print(paste("Param_alternative",i,"/",dim(param_alt)[1]))
print(parameter)
message(paste(Sys.time()," Start!", sep = ""))


power <- rep(0,6)
for (I in 1:R) {
  print(paste(Sys.time(),"Iteration I =",I,"/",R,"= R"))
  temp_power = data.frame(do.call("rbind",pbmclapply(1:M, mc.cores = mc.cores, function(j) {
#  if (j %% 10000 == 0) message(paste(Sys.time(),": Iteration ",j, "/", R, sep = ""))
  dat = generater(maf = parameter$maf, 
                  n.f = parameter$n.f,
                  n.m = parameter$n.m, 
                  r = parameter$r,
                  beta.f = parameter$beta.f,
                  beta.m = parameter$beta.m,
                  err_dist = parameter$err_dist,
                  model = parameter$model,
                  response = parameter$response)

  f <- lm(Y~G, data = dplyr::filter(dat, sex == 1))
  m <- lm(Y~G, data = dplyr::filter(dat, sex == 0))
  
  
  summary.f <- summary(f)
  summary.m <- summary(m)
  
  
  Z.f <- summary.f$coefficients["G", "t value"]
  p.f <- summary.f$coefficients["G", "Pr(>|t|)"]
  
  Z.m <- summary.m$coefficients["G", "t value"]
  p.m <- summary.m$coefficients["G", "Pr(>|t|)"]
  
  x.f <- summary.f$coefficients["G", "Estimate"]
  x.m <- summary.m$coefficients["G", "Estimate"]
  
  
  w.f = (1/summary.f$coefficients["G", "Std. Error"])^2
  w.m = (1/summary.m$coefficients["G", "Std. Error"])^2
  z.TSG.1 = sqrt(w.f/(w.f+w.m))*Z.f+sqrt(w.m/(w.f+w.m))*Z.m
  p.TSG.1 = 2*pnorm(abs(z.TSG.1), lower.tail = F)
  q.TSG.q = Z.f^2+Z.m^2
  p.TSG.q = pchisq(df = 2, q = q.TSG.q, lower.tail = FALSE)
  
  
  sig.f <- summary.f$sigma
  sig.m <- summary.m$sigma
  wt = 1/(ifelse(dat$sex==1, sig.f^2, sig.m^2))
  m1 = summary(lm(Y~G*sex, data = dat, weights = wt))
  V = m1$cov[c("G","G:sex"),c("G","G:sex")]
  beta = m1$coefficients[c("G","G:sex"),"Estimate"]
  z.TI.joint = t(beta)%*%solve(V)%*%beta
  p.TI.joint = pchisq(q = z.TI.joint, df = 2, lower.tail = FALSE)
  
  m0 = summary(lm(Y~G+sex, data = dat, weights = wt))
  p.T.main = m0$coefficients["G","Pr(>|t|)"]
  p = c(p.TSG.f = p.f,  p.TSG.min =  1 - (1 - min(p.f,p.m))^2, p.TSG.1 = p.TSG.1, 
    p.TSG.q = p.TSG.q, p.TI.joint = p.TI.joint,p.T.main = p.T.main)
  p<a
})))
power <-power + colSums(temp_power)
}

p_power <- data.frame(parameter, t(power/M))

saveRDS(p_power, file = paste("out/power/",parameter$err_dist, "/power_rep-", M, 
                           "_scenario-", parameter$scenario,
                           "_beta.f-", parameter$beta.f,
                           "_beta.m-", parameter$beta.m,
                           "_batch-",parameter$batch, 
                           "_k-",parameter$k,
                           "_r-",parameter$r,
                           "_err-", parameter$err_dist,
                           "_model-", parameter$model,
                           "_maf-", parameter$maf,
                           "_nf-", parameter$n.f,
                           "_response-", parameter$response,
                           "_20221120.rds", sep = ""))
