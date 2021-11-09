############
#
# This code was developed by Genie Leonova, Joe Breeden, and Tony Bellotti
#    for purposes os studying instabilities in Cox PH estimation via
#    simulated examples.
#
# Copyright (c) 2021 Leonova, Breeden, and Bellotti
#
############


rm(list = ls())
Sys.setlocale("LC_TIME", "C")
memory.limit(64000)

# library(matlib)
library(survival)
library(car)
library(data.table)
library(ggplot2)
library(scales)
library(grid)

setwd("/jbreeden/coxph")
source("ggplot_themes.R")

# Standard functions to generate risk factors
ln_baseline <- function(a) {0.0025*(1+2*dnorm(log(a*2), 1.8,1)/0.4)}

tfun60 <- function(t) {rnorm(t, 0,0.4) + 1*sin((1:t)/60*2*pi+pi)}

# Stationary AR(1) process for vintage
theta_ar1 <- 0.9
vfun_sar1 <- function(t) {
  v <- rep(0,t)
  for (ix in 2:t) {
    v[ix] <- theta_ar1 * v[ix-1] + rnorm(1,0,0.25)
  }
  v
}

pbin <- function(x){
  1 / (1 + exp(-x))
}

# Finds projection of vector x to vector y
Proj <- function(x, y) {
  if (length(x) != length(y)) {
    stop("x and y should be of the same length")
  }
  lambda <- sum(x*y) / sum(y^2)
  return(lambda*y)
}

# ------------------------- Full business cycle -------------------------

t_max <- 60
seed <- 2209
n <- 10000
scale <- 1.25
beta <- c(1, 1)
cp <- 0  #  1.05^(1/12)-1   Simulate without attrition
n_sims <- 100
alpha <- 0.01
afun <- ln_baseline
delta <- 0
vfun <- vfun_sar1
tfun <- tfun60
dr_target <- 0.02  # This is the target default rate
dr_tol <- 0.00001   # This is the tolerance for hitting the default rate

# Effect  #1 ------------------------------------------------------------------
# The approximate time of calculations  - about 18 hours.

survival_sim_new_effect1 <- function(n, seed, t_max, afun, cp, beta, delta, hscale.prev,
                                     estimation.method = c("all", "efron","breslow","exact"),
                                     dr_only = FALSE){
  set.seed(seed)
  
  # Simulate the start date of each observation (origination date)
  d <- ceiling(runif(n)*t_max)
  
  
  # Simulate TVC (Time-varying covariate):
  
  baseline.haz <- afun(1:t_max)
  
  Q <- matrix(nrow = n, ncol = t_max)
  for (i in 1:n) {
    q_i <- rnorm(t_max, mean=0, sd = 1)
    # Q[i,] <- q_i - Proj(q_i, log(baseline.haz))
    Q[i,] <- q_i - mean(q_i) - Proj(q_i - mean(q_i), log(baseline.haz) - mean(log(baseline.haz)))
    Q[i,] <- (sd(log(baseline.haz)) / sd(Q[i,])) * Q[i,]  ### It's an open question about how this should be set
  }

  # Threshold probablity for failure for each observation
  pf <- runif(n)
  
  # Set the thresholds for probability of attrition
  pa <- matrix(runif(t_max*n), ncol=t_max)
  
  sim1 <- function(hscale) {
    #print(hscale)
    
    X_delta <- t(apply(Q, MARGIN = 1, FUN = function(x){
      res <- (1 - delta) * x + delta * log(hscale*baseline.haz)
      return(res)
    }))
    
    # X_delta_std <- X_delta - mean(X_delta)
    # X_delta_std <- X_delta_std / sd(X_delta_std)
    #
    # test <- apply(X_delta, MARGIN = 1, FUN = function(x){
    #   crossprod(x, baseline.haz)
    # })
    # hist(test, plot = FALSE)
    
    
    # Keep track of time to event and record whether default or censorship;
    # Initialize as censored at last possible time point;
    def <- rep(FALSE,n)
    
    # The following implies that age=1 at the date of origination
    # a is the ending age, by default set to t_max minus origination date
    a <- t_max-d+1
    
    # Compute hazard at each age a_i;
    # get survival prob and simulate default
    chaz <- rep(0,n)
    for (a_i in 1:t_max) {
      
      # First determine cases that are censored
      # For sensored accounts ix, set end age to current age.
      ix <- which(pa[ , a_i] < cp & a > a_i)
      a[ix]<-a_i    ### This is switched off for now, so this should be empty set
      
      # Calculate cumulative hazard
      chaz <- chaz + hscale*baseline.haz[a_i] * exp( beta[1] * X_delta[, a_i])
      
      # Survival is related to cumulative hazard;
      # then determine defaults
      sp <- exp(-chaz)
      ix <- which(sp < pf & a > a_i)
      def[ix] <- TRUE
      a[ix] <- a_i
    }
    
    return(list(data.frame(a,def),X_delta))
  }
  
  err1 <- function(ln_hscale) {
    out <- sim1(exp(ln_hscale))  # We search on a log scale to improve the optimization
    
    # fraction of deafults
    dr <- sum(out[[1]]$def)/length(out[[1]]$def)
    
    dr.err <- abs(dr-dr_target)
    #print(paste(exp(ln_hscale), dr, dr.err))
    
    return(log(dr.err))
  }
  
#  if (hscale.prev == 1) {  # This is the default starting condition
    h.intvl <- log(c(0.0001, 100))
#  } else {
#    h.intvl <- log(c(hscale.prev*0.8, hscale.prev*1.2))
#  }
  h.opt <- optimize(err1, interval=h.intvl, tol=dr_tol)
  hscale <- exp(h.opt$minimum)
  out <- sim1(hscale)
  dr <- sum(out[[1]]$def)/length(out[[1]]$def)
  #print(dr)
  
  training.set <- data.frame(Loan.ID = 1:n, Age = out[[1]]$a, Default = out[[1]]$def)
  X_delta = out[[2]]
  
  # Now build Cox PH models from the simulated data
  cff <- list()
  if (estimation.method == "all") {
    cff[["efron"]] <- coxph(Surv(Age, Default) ~ tt(Loan.ID), data = training.set,
                                      tt = function(acct.id, t,...) {
                                        mapply(acct.id, t, FUN = function(i,j) X_delta[i, j] )
                                      }, 
                                      ties = "efron") 
    cff[["breslow"]] <- coxph(Surv(Age, Default) ~ tt(Loan.ID), data = training.set,
                                      tt = function(acct.id, t,...) {
                                        mapply(acct.id, t, FUN = function(i,j) X_delta[i, j] )
                                      }, 
                                      ties = "breslow") 
    cff[["exact"]] <- coxph(Surv(Age, Default) ~ tt(Loan.ID), data = training.set,
                                      tt = function(acct.id, t,...) {
                                        mapply(acct.id, t, FUN = function(i,j) X_delta[i, j] )
                                      }, 
                                      ties = "exact") 
    
  } else {
  cff[[estimation.method]] <- coxph(Surv(Age, Default) ~ tt(Loan.ID), data = training.set,
              tt = function(acct.id, t,...) {
                mapply(acct.id, t, FUN = function(i,j) X_delta[i, j] )
              }, 
              ties = estimation.method) 
  }
  
  # Estimate the Condition Number of a Matrix (kappa)
  
  tmp <- training.set[rep(1:nrow(training.set), times = training.set$Age),]
  tmp$Age.correct<- unlist(lapply(training.set$Age, FUN = function(x) {seq(1, x, by = 1)}))
  tmp$Default[tmp$Age > tmp$Age.correct] <- FALSE
  tmp$Age <- tmp$Age.correct
  tmp$Age.correct <- NULL
  tmp$hazard <- log(hscale*baseline.haz[tmp$Age])
  tmp$x <- mapply(tmp$Loan.ID, tmp$Age, FUN = function(i,j) X_delta[i, j] )
  #tmp$Age <- factor(tmp$Age)

  mm12 <- model.matrix(~ -1 + hazard + x, data = tmp)
  kappa <- kappa(mm12, exact = TRUE)
  #mm12 <- model.matrix(~ -1 + Age + x, data = tmp)
  #kappa <- kappa(mm12, exact = TRUE)
  
  coef <- list()
  ht <- list()
  pv <- list()
  for (cfn in names(cff)) {
    coef[[cfn]] <- matrix(c(cff[[cfn]]$coef[1],NA, NA),
                   ncol=3, byrow=TRUE)
    
    ht[[cfn]] <- paste0("tt(Loan.ID)=", beta[1])
    pv[[cfn]] <- matrix(c(linearHypothesis(cff[[cfn]], ht[[cfn]])$Pr[2],NA,NA),
                 ncol=3, byrow=TRUE)
  }
  
  data.frame(seed=seed, dr=dr, beta1 = beta[1], coef.efron=coef$efron[1,1], coef.breslow=coef$breslow[1,1], coef.exact=coef$exact[1,1], 
       pv.efron=pv$efron[1,1], pv.breslow=pv$breslow[1,1], pv.exact=pv$exact[1,1], `kappa` = kappa, hscale = exp(h.opt$minimum))
}

# start the clock
ptm <- proc.time()


delta_set <- c(seq(0, 0.5, by = 0.05), seq(0.52, 0.9, by = 0.02), seq(0.91, 0.99, by = 0.01))
full.output <- list()
res.efron <- data.table(Delta = delta_set, `Average absolute error` = as.numeric(NA), `Average error` = as.numeric(NA))
res.breslow <- data.table(Delta = delta_set, `Average absolute error` = as.numeric(NA), `Average error` = as.numeric(NA))
res.exact <- data.table(Delta = delta_set, `Average absolute error` = as.numeric(NA), `Average error` = as.numeric(NA))

hscale.prev <- 1
for (delta in delta_set) {
  print(paste("Delta =", delta))
  
  for (ix in 1:n_sims) {
    sim.res <- survival_sim_new_effect1(n = n, seed = seed + ix, t_max = t_max, hscale.prev, #scale = scale,
                                   afun = ln_baseline, cp = cp, beta = beta, delta = delta, dr_only=FALSE,
                                   estimation.method="all")
    hscale.prev <- sim.res$hscale
    if  (ix==1) R2 <- sim.res
    else R2 <- rbind(R2,sim.res)
  }
  full.output[[paste0("delta=", delta)]] <- R2
  
  res.efron[Delta == delta, c("Average absolute error", "Average error", "Kappa") := 
        list(mean(abs(R2$coef.efron - beta[1]))/beta[1], mean((R2$coef.efron - beta[1]))/beta[1],
             mean(R2$kappa))]
  res.breslow[Delta == delta, c("Average absolute error", "Average error", "Kappa") := 
              list(mean(abs(R2$coef.breslow - beta[1]))/beta[1], mean((R2$coef.breslow - beta[1]))/beta[1],
                   mean(R2$kappa))]
  res.exact[Delta == delta, c("Average absolute error", "Average error", "Kappa") := 
              list(mean(abs(R2$coef.exact - beta[1]))/beta[1], mean((R2$coef.exact - beta[1]))/beta[1],
                   mean(R2$kappa))]
}

# save(list = c("res.efron", "res.breslow", "res.exact", "full.output"), file = "new_sim_effect_1_all.RData")
save(list = c("res.efron", "res.breslow", "res.exact", "full.output"), file = "new_sim_effect_1_all_test1.RData")

# end the clock
print(proc.time() - ptm)



# Effect  #2 ------------------------------------------------------------------
# The approximate time of calculations  - about 13-14 hours.


survival_sim_new_effect2 <- function(n, seed, t_max, scale, afun, cp, beta, delta, hscale.prev,
                                     dr_only = FALSE){
  set.seed(seed)
  
  # Simulate the start date of each observation (origination date)
  d <- ceiling(runif(n)*t_max)
  
  # -------------------------------------
  # Simulate TVC (Time-varying covariate):
  
  baseline.haz <- afun(1:t_max)
  
  S <- matrix(nrow = n, ncol = t_max)
  Q1 <- matrix(nrow = n, ncol = t_max)
  Q2 <- matrix(nrow = n, ncol = t_max)
  for (i in 1:n) {
    s_i <- rnorm(t_max, mean=0, sd = 1) 
    # S[i,] <- s_i - Proj(s_i, log(baseline.haz))
    S[i,] <- s_i - mean(s_i) - Proj(s_i - mean(s_i), log(baseline.haz) - mean(log(baseline.haz)))
    
    q1_i <- rnorm(t_max, mean=0, sd = 1) 
    # Q1[i,] <- q1_i - Proj(q1_i, log(baseline.haz)) - Proj(q1_i, S[i, ])
    Q1[i,] <- q1_i - mean(q1_i) - Proj(q1_i - mean(q1_i), log(baseline.haz) - mean(log(baseline.haz))) - 
      Proj(q1_i - mean(q1_i), S[i, ])
  }
  
  sdS <- sd(S, na.rm=TRUE)
  sdQ1 <- sd(Q1, na.rm=TRUE)
  Q1 <- sdS / sdQ1 * Q1
    
  for (i in 1:n) {
    q2_i <- rnorm(t_max, mean=0, sd = 1) 
    # Q2[i,] <- q2_i - Proj(q2_i, log(baseline.haz)) - Proj(q2_i, S[i,]) - Proj(q2_i, Q1[i,])
    Q2[i,] <- q2_i - mean(q2_i) - Proj(q2_i - mean(q2_i), log(baseline.haz) - mean(log(baseline.haz))) - 
      Proj(q2_i - mean(q2_i), S[i,]) - Proj(q2_i - mean(q2_i), Q1[i,])
  }

  sdQ2 <- sd(Q2, na.rm=TRUE)
  Q2 <- sdS / sdQ2 * Q2
  
  X1_delta <- (1 - delta) * Q1 + delta * S
  X2_delta <- (1 - delta) * Q2 + delta * S
  
  # Threshold probablity for failure for each observation
  pf <- runif(n)
  
  # Set the thresholds for probability of attrition
  pa <- matrix(runif(t_max*n), ncol=t_max)
  
  sim2 <- function(hscale) {
    # Keep track of time to event and record whether default or censorship;
    # Initialize as censored at last possible time point;
    def <- rep(FALSE,n)
    
    # The following implies that age=1 at the date of origination
    a <- t_max-d+1
    
    # Compute hazard at each age a_i;
    # get survival prob and simulate default
    chaz <- rep(0,n)
    for (a_i in 1:t_max) {
      
      # First determine cases that are censored
      ix <- which(pa[ , a_i] < cp & a > a_i)
      a[ix]<-a_i
      
      # Calculate cumulative hazard
      chaz <- chaz + hscale*baseline.haz[a_i] * exp(beta[1] * X1_delta[, a_i] + beta[2] * X2_delta[, a_i])
      
      # Survival is related to cumulative hazard;
      # then determine defaults
      sp <- exp(-chaz)
      ix <- which(sp < pf & a > a_i)
      def[ix] <- TRUE
      a[ix] <- a_i    ### This is switched off for now, so this should be empty set
    }
    
    return(list(data.frame(a,def),X1_delta,X2_delta))
  }
  
  err2 <- function(ln_hscale) {
    out <- sim2(exp(ln_hscale))  # We search on a log scale to improve the optimization
    
    # fraction of deafults
    dr <- sum(out[[1]]$def)/length(out[[1]]$def)
    
    dr.err <- abs(dr-dr_target)
    #print(paste(exp(ln_hscale), dr, dr.err))
    
    return(log(dr.err))
  }
  
  if (hscale.prev == 1) {  # This is the default starting condition
    h.intvl <- log(c(0.01, 1))
  } else {
    h.intvl <- log(c(hscale.prev*0.8, hscale.prev*1.2))
  }
  h.opt <- optimize(err2, interval=h.intvl, tol=dr_tol)
  hscale <- exp(h.opt$minimum)
  out <- sim2(hscale = hscale)
  dr <- sum(out[[1]]$def)/length(out[[1]]$def)
  #print(dr)
  
  training.set <- data.frame(Loan.ID = 1:n, Loan.ID.duplicate = 1:n, Age = out[[1]]$a, Default = out[[1]]$def)
  X1_delta = out[[2]]
  X2_delta = out[[3]]
  
  # Now build Cox PH models from the simulated data
  
  c1 <- coxph(Surv(Age, Default) ~ tt(Loan.ID) + tt(Loan.ID.duplicate), data = training.set,
              tt = list(function(Loan.ID, t,...) {
                mapply(Loan.ID, t, FUN = function(i,j) X1_delta[i, j] )
              }, function(Loan.ID.duplicate, t,...) {
                mapply(Loan.ID.duplicate, t, FUN = function(i,j) X2_delta[i, j] )
              }) )
  
  
  # Estimate the Condition Number of a Matrix (kappa) -------------------------
  
  tmp <- training.set[rep(1:nrow(training.set), times = training.set$Age), ]
  tmp$Age.correct<- unlist(lapply(training.set$Age, FUN = function(x) {seq(1, x, by = 1)}))
  tmp$Default[tmp$Age > tmp$Age.correct] <- FALSE
  tmp$Age <- tmp$Age.correct
  tmp$Age.correct <- NULL
  tmp$x1 <- mapply(tmp$Loan.ID, tmp$Age, FUN = function(i,j) X1_delta[i, j] )
  tmp$x2 <- mapply(tmp$Loan.ID, tmp$Age, FUN = function(i,j) X2_delta[i, j] )
  tmp$hazard <- log(hscale*baseline.haz[tmp$Age])
  #tmp$Age <- factor(tmp$Age)
  
  mm12 <- model.matrix(~ -1 + hazard + x1 + x2, data = tmp)
  kappa <- kappa(mm12, exact = TRUE)
  #mm12 <- model.matrix(~ -1 + Age + x1 + x2, data = tmp)
  #kappa <- kappa(mm12, exact = TRUE)
  
  coef <- matrix(c(c1$coef[1], c1$coef[2]),
                 ncol=2, byrow=TRUE)
  
  # ht <- paste0("tt(Loan.ID)=", beta[1])
  # pv <- matrix(c(linearHypothesis(c1, ht)$Pr[2],NA,NA),
  #              ncol=3, byrow=TRUE)
  
  
  data.frame(seed, dr, beta1 = beta[1], beta2 = beta[2], coef, `kappa` = kappa, hscale = exp(h.opt$minimum))
}


# start the clock
ptm <- proc.time()


delta_set <- c(seq(0, 0.5, by = 0.05), seq(0.52, 0.9, by = 0.02), seq(0.91, 0.99, by = 0.01))
full.output <- list()
res <- data.table(Delta = delta_set)

hscale.prev <- 1
for (delta in delta_set) {
  print(paste("Delta =", delta))
  
  for (ix in 1:n_sims) {
    R1 <- survival_sim_new_effect2(n = n, seed = seed + ix, t_max = t_max, hscale.prev=hscale.prev, #scale = scale,
                                   afun = ln_baseline, cp = cp, beta = beta, delta = delta, dr_only=FALSE)
    hscale.prev <- R1$hscale
    if  (ix==1) R2 <- R1
    else R2 <- rbind(R2,R1)
  }
  full.output[[paste0("delta=", delta)]] <- R2
  
  res[Delta == delta, c("Average absolute error for X1", "Average error for X1",
                        "Average absolute error for X2", "Average error for X2",
                        "Kappa") :=
        list(mean(abs(R2$X1 - beta[1]))/beta[1], mean((R2$X1 - beta[1]))/beta[1], 
             mean(abs(R2$X2 - beta[2]))/beta[2], mean((R2$X2 - beta[2]))/beta[2],
             mean(R2$kappa))]
}

save(list = c("res", "full.output"), file = "new_sim_effect_2_test1.RData")

# end the clock
print(proc.time() - ptm)


# Effect  #2A -----------------------------------------------------------------
# The approximate time of calculations  - about 13-14 hours.


survival_sim_new_effect2A <- function(n, seed, t_max, scale, afun, cp, beta, delta, hscale.prev,
                                      dr_only = FALSE){
  set.seed(seed)
  
  # Simulate the start date of each observation (origination date)
  d <- ceiling(runif(n)*t_max)
  
  # Vintage effect over each year
  f_v <- vfun(t_max)
  v <- f_v[d]
  
  # -------------------------------------
  # Simulate TVC (Time-varying covariate):
  
  baseline.haz <- afun(1:t_max)
  
  Q1 <- rnorm(n, mean = 0, sd = sd(v))
  Q2 <- rnorm(n, mean = 0, sd = sd(v))

  X1_delta <- (1 - delta) * Q1 + delta * v
  X2_delta <- (1 - delta) * Q2 + delta * v
  
  # Threshold probablity for failure for each observation
  pf <- runif(n)
  
  # Set the thresholds for probability of attrition
  pa <- matrix(runif(t_max*n), ncol=t_max)
  
  sim2A <- function(hscale) {
    # Keep track of time to event and record whether default or censorship;
    # Initialize as censored at last possible time point;
    def <- rep(FALSE,n)
    
    # The following implies that age=1 at the date of origination
    a <- t_max-d+1
    
    # Compute hazard at each age a_i;
    # get survival prob and simulate default
    chaz <- rep(0,n)
    for (a_i in 1:t_max) {
      
      # First determine cases that are censored
      ix <- which(pa[ , a_i] < cp & a > a_i)
      a[ix]<-a_i
      
      # Calculate cumulative hazard
      chaz <- chaz + hscale*baseline.haz[a_i] * exp(beta[1] * X1_delta + beta[2] * X2_delta)
      
      # Survival is related to cumulative hazard;
      # then determine defaults
      sp <- exp(-chaz)
      ix <- which(sp < pf & a > a_i)
      def[ix] <- TRUE
      a[ix] <- a_i    ### This is switched off for now, so this should be empty set
    }
    
    return(list(data.frame(a,def),X1_delta,X2_delta))
  }
  
  err2 <- function(ln_hscale) {
    out <- sim2A(exp(ln_hscale))  # We search on a log scale to improve the optimization
    
    # fraction of deafults
    dr <- sum(out[[1]]$def)/length(out[[1]]$def)
    
    dr.err <- abs(dr-dr_target)
    #print(paste(exp(ln_hscale), dr, dr.err))
    
    return(log(dr.err))
  }
  
#  if (hscale.prev == 1) {  # This is the default starting condition
    h.intvl <- log(c(0.0001, 100))
#  } else {
#    h.intvl <- log(c(hscale.prev*0.8, hscale.prev*1.2))
#  }
  h.opt <- optimize(err2, interval=h.intvl, tol=dr_tol)
  hscale <- exp(h.opt$minimum)
  out <- sim2A(hscale = hscale)
  dr <- sum(out[[1]]$def)/length(out[[1]]$def)
  
  training.set <- data.frame(Age = out[[1]]$a, Default = out[[1]]$def, x1 = out[[2]], x2 = out[[3]])
  
  # Now build Cox PH models from the simulated data
  
  c1 <- coxph(Surv(Age, Default) ~ x1 + x2, data = training.set)
  
  
  # Estimate the Condition Number of a Matrix (kappa) -------------------------
  
  tmp <- training.set[rep(1:nrow(training.set), times = training.set$Age), ]
  tmp$Age.correct<- unlist(lapply(training.set$Age, FUN = function(x) {seq(1, x, by = 1)}))
  tmp$Default[tmp$Age > tmp$Age.correct] <- FALSE
  tmp$Age <- tmp$Age.correct
  tmp$Age.correct <- NULL
  tmp$hazard <- log(hscale*baseline.haz[tmp$Age])
  #tmp$Age <- factor(tmp$Age)
  mm12 <- model.matrix(~ -1 + hazard + x1 + x2, data = tmp)
  kappa <- kappa(mm12, exact = TRUE)
  
  coef <- matrix(c(c1$coef[1], c1$coef[2]),
                 ncol=2, byrow=TRUE)
  
  # ht <- paste0("tt(Loan.ID)=", beta[1])
  # pv <- matrix(c(linearHypothesis(c1, ht)$Pr[2],NA,NA),
  #              ncol=3, byrow=TRUE)
  
  
  data.frame(seed, dr, beta1 = beta[1], beta2 = beta[2], coef, `kappa` = kappa, hscale = exp(h.opt$minimum))
}


# start the clock
ptm <- proc.time()


delta_set <- c(seq(0, 0.5, by = 0.05), seq(0.52, 0.9, by = 0.02), seq(0.91, 0.99, by = 0.01))
full.output <- list()
res <- data.table(Delta = delta_set)

beta <- c(1, 1)

hscale.prev <- 1
for (delta in delta_set) {
  print(paste("Delta =", delta))
  
  for (ix in 1:n_sims) {
    R1 <- survival_sim_new_effect2A(n = n, seed = seed + ix, t_max = t_max, hscale.prev=hscale.prev, #scale = scale,
                                   afun = ln_baseline, cp = cp, beta = beta, delta = delta, dr_only=FALSE)
    hscale.prev <- R1$hscale
    if  (ix==1) R2 <- R1
    else R2 <- rbind(R2,R1)
  }
  full.output[[paste0("delta=", delta)]] <- R2
  
  res[Delta == delta, c("Average absolute error for X1", "Average error for X1",
                        "Average absolute error for X2", "Average error for X2",
                        "Kappa") :=
        list(mean(abs(R2$X1 - beta[1]))/beta[1], mean((R2$X1 - beta[1]))/beta[1], 
             mean(abs(R2$X2 - beta[2]))/beta[2], mean((R2$X2 - beta[2]))/beta[2],
             mean(R2$kappa))]
}

save(list = c("res", "full.output"), file = "new_sim_effect_2A.RData")

# end the clock
print(proc.time() - ptm)


# for (delta in seq(0.991, 0.999, by = 0.001)) {
#   print(paste("Delta =", delta))
#   
#   for (ix in 1:n_sims) {
#     R1 <- survival_sim_new_effect2A(n = n, seed = seed + ix, t_max = t_max, scale = scale,
#                                     afun = ln_baseline, vfun = vfun_sar1, cp = cp, beta = beta, delta = delta, dr_only=FALSE)
#     if  (ix==1) R2 <- R1
#     else R2 <- rbind(R2,R1)
#   }
#   full.output[[paste0("delta=", delta)]] <- R2
# }
# 
# 
# res0 <- data.table(Delta = seq(0.991, 0.999, by = 0.001))
# for (delta in seq(0.991, 0.999, by = 0.001)) {
#   
#   res0[Delta == delta,  c("Average absolute error for X1", "Average absolute error for X2") :=
#          list(mean(abs(full.output[[paste0("delta=", delta)]]$X1 - beta[1])), 
#               mean(abs(full.output[[paste0("delta=", delta)]]$X2 - beta[2])))]
#   
#   res0[Delta == delta, "Kappa" := mean(full.output[[paste0("delta=", delta)]]$kappa)]
# }
# 
# 
# res <- rbindlist(list(res, res0))
# save(list = c("res", "full.output"), file = "new_sim_effect_2A.RData")


# Effect  #2B -----------------------------------------------------------------
# The approximate time of calculations  - about 9 hours.

survival_sim_new_effect2B <- function(n, seed, t_max, scale, afun, tfun, cp, beta, delta, hscale.prev,
                                      dr_only = FALSE){
  set.seed(seed)
  
  # Simulate the start date of each observation (origination date)
  d <- ceiling(runif(n)*t_max)
  
  # -------------------------------------
  # Simulate TVC (Time-varying covariate):
  
  baseline.haz <- afun(1:t_max)
  
  # Calendar-time effect extending beyond data series for extrapolation
  f_t <- tfun(t_max)

  Q1_t <- rnorm(t_max, mean = 0, sd = 1)
  # Q1_t <- Q1_t - Proj(Q1_t, f_t)
  Q1_t <- Q1_t - mean(Q1_t) - Proj(Q1_t - mean(Q1_t), f_t - mean(f_t))
  Q1_t <- (sd(f_t) / sd(Q1_t)) * Q1_t

  Q2_t <- rnorm(t_max, mean = 0, sd = 1)
  # Q2_t <- Q2_t - Proj(Q2_t, f_t) - Proj(Q2_t, Q1_t)
  Q2_t <- Q2_t - mean(Q2_t) - Proj(Q2_t - mean(Q2_t), f_t - mean(f_t)) - 
    Proj(Q2_t - mean(Q2_t), Q1_t)
  Q2_t <- (sd(f_t) / sd(Q2_t)) * Q2_t
  
  X1_t_delta <- (1 - delta) * Q1_t + delta * f_t
  X2_t_delta <- (1 - delta) * Q2_t + delta * f_t
  
  # Threshold probablity for failure for each observation
  pf <- runif(n)
  
  # Set the thresholds for probability of attrition
  pa <- matrix(runif(t_max*n), ncol=t_max)
  
  sim2B <- function(hscale) {
    # Keep track of time to event and record whether default or censorship;
    # Initialize as censored at last possible time point;
    def <- rep(FALSE,n)
    
    # The following implies that age=1 at the date of origination.
    # This assignment resets the censoring calculations
    a <- t_max-d+1
    
    # Compute hazard at each age a_i;
    # get survival prob and simulate default
    chaz <- rep(0,n)
    for (a_i in 1:t_max) {

      # First determine cases that are censored
      ix <- which(pa[, a_i] < cp & a > a_i)
      a[ix]<-a_i
      
      # Calculate cumulative hazard
      chaz <- chaz + hscale*baseline.haz[a_i] * exp(beta[1] * X1_t_delta[d+a_i-1] + beta[2] * X2_t_delta[d+a_i-1])

      # Survival is related to cumulative hazard;
      # then determine defaults
      sp <- exp(-chaz)
      ix <- which(sp < pf & a > a_i)  # if the default probability is greater than the threshold and the account age is less than the max age (censored date)
      def[ix] <- TRUE
      a[ix] <- a_i    # Accounts that default are added to the censored list to prevent re-defaults
    }
    
    return(list(data.frame(a,def),X1_t_delta,X2_t_delta))
  }
  
  err2 <- function(ln_hscale) {
    out <- sim2B(exp(ln_hscale))  # We search on a log scale to improve the optimization
    
    # fraction of deafults
    dr <- sum(out[[1]]$def)/length(out[[1]]$def)
    
    dr.err <- abs(dr-dr_target)
    #print(paste(exp(ln_hscale), dr, dr.err))
    
    return(log(dr.err))
  }
  
  if (hscale.prev == 1) {  # This is the default starting condition
    h.intvl <- log(c(0.01, 1))
  } else {
    h.intvl <- log(c(hscale.prev*0.8, hscale.prev*1.2))
  }
  h.opt <- optimize(err2, interval=h.intvl, tol=dr_tol)
  hscale <- exp(h.opt$minimum)
  out <- sim2B(hscale = hscale)
  dr <- sum(out[[1]]$def)/length(out[[1]]$def)
  #print(paste(hscale,dr))
  
  training.set <- data.frame(Age = out[[1]]$a, Default = out[[1]]$def, Vintage.Date = d, 
                             Vintage.Date.dublicate = d)
  X1_t_delta = out[[2]]
  X2_t_delta = out[[3]]
  
  # Now build Cox PH models from the simulated data
  
  c1 <- coxph(Surv(Age, Default) ~ tt(Vintage.Date) + tt(Vintage.Date.dublicate), data = training.set,
              tt = list(function(Vintage.Date, t,...) X1_t_delta[Vintage.Date + t-1], 
                        function(Vintage.Date.dublicate, t,...) X2_t_delta[Vintage.Date.dublicate + t-1]) )
  

  # Estimate the Condition Number of a Matrix (kappa) -------------------------
  
  tmp <- training.set[rep(1:nrow(training.set), times = training.set$Age), ]
  tmp$Age.correct<- unlist(lapply(training.set$Age, FUN = function(x) {seq(1, x, by = 1)}))
  tmp$Default[tmp$Age > tmp$Age.correct] <- FALSE
  tmp$Age <- tmp$Age.correct
  tmp$Age.correct <- NULL
  
  tmp$x1 <- X1_t_delta[tmp$Age + tmp$Vintage.Date-1]
  tmp$x2 <- X2_t_delta[tmp$Age + tmp$Vintage.Date-1]
  tmp$hazard <- log(hscale*baseline.haz[tmp$Age])
  #tmp$Age <- factor(tmp$Age)
  mm12 <- model.matrix(~ -1 + hazard + x1 + x2, data = tmp)
  kappa <- kappa(mm12, exact = TRUE)
  
  coef <- matrix(c(c1$coef[1], c1$coef[2]),
                 ncol=2, byrow=TRUE)
  
  # ht <- paste0("tt(Loan.ID)=", beta[1])
  # pv <- matrix(c(linearHypothesis(c1, ht)$Pr[2],NA,NA),
  #              ncol=3, byrow=TRUE)
  
  
  data.frame(seed, dr, beta1 = beta[1], beta2 = beta[2], coef, `kappa` = kappa, hscale = exp(h.opt$minimum))
}

# start the clock
ptm <- proc.time()


delta_set <- c(seq(0, 0.5, by = 0.05), seq(0.52, 0.9, by = 0.02), seq(0.91, 0.99, by = 0.01))
full.output <- list()
res <- data.table(Delta = delta_set)

hscale.prev <- 1
for (delta in delta_set) {
  print(paste("Delta =", delta))
  
  for (ix in 1:n_sims) {
    R1 <- survival_sim_new_effect2B(n = n, seed = seed + ix, t_max = t_max, hscale.prev=hscale.prev, #scale = scale,
                                    afun = ln_baseline, tfun = tfun60, cp = cp, beta = beta, 
                                    delta = delta, dr_only=FALSE)
    hscale.prev <- R1$hscale
    if  (ix==1) R2 <- R1
    else R2 <- rbind(R2,R1)
  }
  full.output[[paste0("delta=", delta)]] <- R2
  
  res[Delta == delta, c("Average absolute error for X1", "Average error for X1",
                        "Average absolute error for X2", "Average error for X2",
                        "Kappa") :=
        list(mean(abs(R2$X1 - beta[1]))/beta[1], mean((R2$X1 - beta[1]))/beta[1], 
             mean(abs(R2$X2 - beta[2]))/beta[2], mean((R2$X2 - beta[2]))/beta[2],
             mean(R2$kappa))]
}

save(list = c("res", "full.output"), file = "new_sim_effect_2B_test1.RData")

# end the clock
print(proc.time() - ptm)


# Effect  #2C -----------------------------------------------------------------
# The approximate time of calculations  - about 115 hours which is about 4.8 days. 

d <- list()
baseline.haz <- list()
f_t <- list()
E1_t <- list()
E2_t <- list()
pf <- list()
pa <- list()

for (ix in 1:n_sims) {
  # Simulate the start date of each observation (origination date)
  d[[ix]] <- ceiling(runif(n)*t_max)
  
  # -------------------------------------
  # Simulate TVC (Time-varying covariate):
  
  baseline.haz[[ix]] <- afun(1:t_max)
  
  # Calendar-time effect extending beyond data series for extrapolation
  f_t[[ix]] <- tfun(t_max)
  
  E1_t[[ix]] <- rnorm(t_max, mean = 0, sd = 1)
  E1_t[[ix]] <- E1_t[[ix]] - Proj(E1_t[[ix]], f_t[[ix]]) - Proj(E1_t[[ix]], 1:t_max)
  E1_t[[ix]] <- (sd(f_t[[ix]]) / sd(E1_t[[ix]])) * E1_t[[ix]]
  
  E2_t[[ix]] <- rnorm(t_max, mean = 0, sd = 1)
  E2_t[[ix]] <- E2_t[[ix]] - Proj(E2_t[[ix]], f_t[[ix]]) - Proj(E2_t[[ix]], 1:t_max)
  E2_t[[ix]] <- (sd(f_t[[ix]]) / sd(E2_t[[ix]])) * E2_t[[ix]]
  
  # Threshold probablity for failure for each observation
  pf[[ix]] <- runif(n)
  
  # Set the thresholds for probability of attrition
  pa[[ix]] <- matrix(runif(t_max*n), ncol=t_max)
}


survival_sim_new_effect2C <- function(ix, n, seed, t_max, scale, afun, vfun, tfun, cp, beta, hscale.prev,
                                     delta_t, delta, dr_only = FALSE){
#  set.seed(seed)
  
  
  X1_t_delta <- (1 - delta) * E1_t[[ix]] + delta * f_t[[ix]] + delta_t * (1:t_max)
  X2_t_delta <- (1 - delta) * E2_t[[ix]] + delta * f_t[[ix]] + delta_t * (1:t_max)
  
  sim2C <- function(hscale) {
    # Keep track of time to event and record whether default or censorship;
    # Initialize as censored at last possible time point;
    def <- rep(FALSE,n)
    
    # The following implies that age=1 at the date of origination
    a <- t_max-d[[ix]]+1
    
    # Compute hazard at each age a_i;
    # get survival prob and simulate default
    chaz <- rep(0,n)
    for (a_i in 1:t_max) {
      
      # First determine cases that are censored
      ic <- which(pa[[ix]][, a_i] < cp & a > a_i)
      a[ic]<-a_i
      
      # Calculate cumulative hazard
      chaz <- chaz + hscale*baseline.haz[[ix]][a_i] * exp(beta[1] * X1_t_delta[d[[ix]]+a_i-1] + beta[2] * X2_t_delta[d[[ix]]])
      
      # Survival is related to cumulative hazard;
      # then determine defaults
      sp <- exp(-chaz)
      ic <- which(sp < pf[[ix]] & a > a_i)
      def[ic] <- TRUE
      a[ic] <- a_i
    }
    
    return(list(data.frame(a,def),X1_t_delta,X2_t_delta))
    
    sum(def)/length(def)
  }
  
  err2C <- function(ln_hscale) {
    out <- sim2C(exp(ln_hscale))  # We search on a log scale to improve the optimization
    
    # fraction of deafults
    dr <- sum(out[[1]]$def)/length(out[[1]]$def)
    
    dr.err <- abs(dr-dr_target)
    #print(paste(exp(ln_hscale), dr, dr.err))
    
    return(log(dr.err))
  }
  
  # if (hscale.prev == 1) {  # This is the default starting condition
  #   h.intvl <- log(c(0.01, 1))
  # } else {
  #   h.intvl <- log(c(hscale.prev*0.8, hscale.prev*1.2))
  # }
  
  h.rng <- c(0.0001, 1000)
  h.intvl <- log(h.rng) 
  dr <- NA
  
  for (counter in 1:4) { #No infinite loops
    h.opt <- optimize(err2C, interval=h.intvl, tol=dr_tol)
    hscale <- exp(h.opt$minimum)
    out <- sim2C(hscale = hscale)
    dr <- sum(out[[1]]$def)/length(out[[1]]$def)
    if (abs(dr-dr_target) > 0.02) {
      #print(paste(ix, toString(h.rng), dr, sep=": "))
      h.rng <- c(hscale/10000, hscale*10000)
      h.intvl <- log(h.rng) 
    } else {
      break
    }
  }
  
  if(abs(dr-dr_target) > 0.1) {
    print(paste(ix, toString(h.rng), dr, sep=": "))
    return(data.frame(seed, dr, beta1 = NA, beta2 = NA, X1=NA, X2=NA, `kappa` = NA, hscale = NA))
  }
  
  training.set <- data.frame(Age = out[[1]]$a, Default = out[[1]]$def, Vintage.Date = d[[ix]], 
                             Vintage.Date.dublicate = d[[ix]])
  X1_t_delta = out[[2]]
  X2_t_delta = out[[3]]
  
  # Now build Cox PH models from the simulated data
  
  c1 <- coxph(Surv(Age, Default) ~ tt(Vintage.Date) + tt(Vintage.Date.dublicate), data = training.set,
              tt = list(function(Vintage.Date, t,...) X1_t_delta[Vintage.Date + t-1], 
                        function(Vintage.Date.dublicate, t,...) X2_t_delta[Vintage.Date.dublicate + t-1]) )
  
  # Estimate the Condition Number of a Matrix (kappa) -------------------------
  
  tmp <- training.set[rep(1:nrow(training.set), times = training.set$Age), ]
  tmp$Age.correct<- unlist(lapply(training.set$Age, FUN = function(x) {seq(1, x, by = 1)}))
  tmp$Default[tmp$Age > tmp$Age.correct] <- FALSE
  tmp$Age <- tmp$Age.correct
  tmp$Age.correct <- NULL
  
  tmp$x1 <- X1_t_delta[tmp$Age + tmp$Vintage.Date - 1]
  tmp$x2 <- X2_t_delta[tmp$Age + tmp$Vintage.Date - 1]
  tmp$hazard <- log(hscale*baseline.haz[[ix]][tmp$Age])
  mm12 <- model.matrix(~ -1 + hazard + x1 + x2, data = tmp)
  kappa <- kappa(mm12, exact = TRUE)
  
  coef <- matrix(c(c1$coef[1], c1$coef[2]),
                 ncol=2, byrow=TRUE)
  
  # ht <- paste0("tt(Loan.ID)=", beta[1])
  # pv <- matrix(c(linearHypothesis(c1, ht)$Pr[2],NA,NA),
  #              ncol=3, byrow=TRUE)
  
  
  data.frame(seed, dr, beta1 = beta[1], beta2 = beta[2], coef, `kappa` = kappa, hscale = exp(h.opt$minimum))
}

delta_set0 <- c(seq(0, 0.99, by = 0.02))
delta_set <- seq(0.1, 0.3, by = 0.02)

# start the clock
ptm <- proc.time()

for (delta_t in delta_set) {
  
  full.output <- list()
  res <- as.data.table(expand.grid(Delta_t = delta_t, Delta = delta_set0))
  hscale.prev <- 1
  
  for (delta in delta_set0) {
    
    print(paste0("Delta_t=", delta_t, ", Delta=", delta))
    
    for (ix in 1:n_sims) {
      R1 <- survival_sim_new_effect2C(ix = ix, n = n, seed = seed + ix, t_max = t_max, hscale.prev=hscale.prev, #scale = scale,
                                     afun = ln_baseline, tfun = tfun60, vfun = vfun_sar1, 
                                     cp = cp, beta = beta, delta_t = delta_t, delta = delta, dr_only=FALSE)
      
      hscale.prev <- R1$hscale
      #print(paste(hscale.prev, R1$dr))
      if  (ix==1) R2 <- R1
      else R2 <- rbind(R2,R1)
    }
    full.output[[paste0("delta_t=", delta_t, ", delta=", delta)]] <- R2
    
    res[Delta_t == delta_t & Delta == delta, c("Average absolute error for X1", "Average error for X1",
                                                   "Average absolute error for X2", "Average error for X2",
                                                   "Kappa") :=
          list(mean(abs(R2$X1 - beta[1]))/beta[1], mean((R2$X1 - beta[1]))/beta[1], 
               mean(abs(R2$X2 - beta[2]))/beta[2], mean((R2$X2 - beta[2]))/beta[2],
               mean(R2$kappa))]
    
  }
  save(res, full.output, file = paste0("effect2C_delta_t_", delta_t, ".RData"))
}


# end the clock
print(proc.time() - ptm)


# Effect  #3A -----------------------------------------------------------------
# The approximate time of calculations  - about 115 hours which is about 4.8 days. 

# Dumping the initialization in the global space. Not good, but quick
# Doing this to make sure that each change in parameter settings tests the same set of functions

d <- list()
baseline.haz <- list()
f_t <- list()
f_v <- list()
pf <- list()
pa <- list()

set.seed(seed)  # Reset generator to the global seed

for (ix in 1:n_sims) {
  # Simulate the start date of each observation (origination date)
  d[[ix]] <- ceiling(runif(n)*t_max)
  
  # -------------------------------------
  # Simulate TVC (Time-varying covariate):
  
  baseline.haz[[ix]] <- 0.0025 # afun(1:t_max)
  
  # Calendar-time effect extending beyond data series for extrapolation
  f_t[[ix]] <- rev(vfun(t_max))  # We're using vfun instead of tfun to create a symmetric result
  f_t[[ix]] <- f_t[[ix]] - f_t[[ix]][t_max/2]
  
  # Vintage effect over each year
  f_v[[ix]] <- vfun(t_max)
  f_v[[ix]] <- f_v[[ix]] - f_v[[ix]][t_max/2]
  
  # Threshold probablity for failure for each observation
  pf[[ix]] <- runif(n)
  
  # Set the thresholds for probability of attrition
  pa[[ix]] <- matrix(runif(t_max*n), ncol=t_max)
}

# f_t.mat <- do.call(rbind, f_t)
# f_t.sd <- apply(f_t.mat, 2, sd)
# 
# f_v.mat <- do.call(rbind, f_v)
# f_v.sd <- apply(f_v.mat, 2, sd)
# 
# plot(f_v.sd, type='l')
# lines(f_t.sd, col=2)
# 
# plot(f_t[[1]], type='l', ylim=c(-2,2))
# for (ix in 2:n_sims) {
#   lines(rev(f_t[[ix]]))
# }


survival_sim_new_effect3A <- function(ix, n, seed, t_max, scale, afun, vfun, tfun, cp, beta, hscale.prev,
                                     delta_t, delta_v, dr_only = FALSE){
  
  #  set.seed(seed)  # Removed, because everything here is deterministic now.
  
  E_t <- f_t[[ix]] - Proj(f_t[[ix]], rep(1, t_max)) - Proj(f_t[[ix]], 1:t_max - Proj(1:t_max, rep(1, t_max)))
  E_t <- (sd(f_t[[ix]]) / sd(E_t)) * E_t
  
  G_v <- f_v[[ix]] - Proj(f_v[[ix]], rep(1, t_max)) - Proj(f_v[[ix]], 1:t_max - Proj(1:t_max, rep(1, t_max)))
  G_v <- (sd(f_v[[ix]]) / sd(G_v)) * G_v
  
  E_t_delta <- delta_t * (1:t_max) + E_t  
  G_v_delta <- delta_v * (1:t_max) + G_v
  
  # plot(E_t, type='l')
  # lines(G_v, col=2)
  # plot(E_t_delta, type='l')
  # lines(G_v_delta, col=2)
  
  sim3A <- function(hscale) {
    # Keep track of time to event and record whether default or censorship;
    # Initialize as censored at last possible time point;
    def <- rep(FALSE,n)
    
    # The following implies that age=1 at the date of origination
    a <- t_max-d[[ix]]+1
    
    # Compute hazard at each age a_i;
    # get survival prob and simulate default
    chaz <- rep(0,n)
    for (a_i in 1:t_max) {
      
      # First determine cases that are censored
      ic <- which(pa[[ix]][, a_i] < cp & a > a_i)
      a[ic]<-a_i
      
      # Calculate cumulative hazard
      chaz <- chaz + hscale*baseline.haz[[ix]][a_i] * exp(beta[1] * E_t_delta[d[[ix]]+a_i-1] + beta[2] * G_v_delta[d[[ix]]])
      
      # Survival is related to cumulative hazard;
      # then determine defaults
      sp <- exp(-chaz)
      ic <- which(sp < pf[[ix]] & a > a_i)
      def[ic] <- TRUE
      a[ic] <- a_i
    }
    
    return(list(data.frame(a,def),E_t_delta,G_v_delta))
  }
  
  err3A <- function(ln_hscale) {
    out <- sim3A(exp(ln_hscale))  # We search on a log scale to improve the optimization
    
    # fraction of deafults
    dr <- sum(out[[1]]$def)/length(out[[1]]$def)
    
    dr.err <- abs(dr-dr_target)
    #print(paste(exp(ln_hscale), dr, dr.err))
    
    return(log(dr.err))
  }
  
  if (hscale.prev == 1) {  # This is the default starting condition
    h.intvl <- log(c(0.0001, 100))
  } else {
    h.intvl <- log(c(hscale.prev/10, hscale.prev*10))
  }
  h.opt <- optimize(err3A, interval=h.intvl, tol=dr_tol)
  hscale <- exp(h.opt$minimum)
  out <- sim3A(hscale = hscale)
  dr <- sum(out[[1]]$def)/length(out[[1]]$def)
  #print(paste(hscale,dr))
  
  E_t_delta = out[[2]]
  G_v_delta = out[[3]]
  training.set <- data.frame(Age = out[[1]]$a, Default = out[[1]]$def, Vintage.Date = d[[ix]], 
                             x2 = G_v_delta[d[[ix]]])
  
  # Now build Cox PH models from the simulated data
  
  c1 <- coxph(Surv(Age, Default) ~ x2 + tt(Vintage.Date), data = training.set,
              tt = function(v, t,...) E_t_delta[v+t-1] )
  
  # Estimate the Condition Number of a Matrix (kappa) -------------------------
  
  tmp <- training.set[rep(1:nrow(training.set), times = training.set$Age), ]
  tmp$Age.correct<- unlist(lapply(training.set$Age, FUN = function(x) {seq(1, x, by = 1)}))
  tmp$Default[tmp$Age > tmp$Age.correct] <- FALSE
  tmp$Age <- tmp$Age.correct
  tmp$Age.correct <- NULL
  
  tmp$x1 <- E_t_delta[tmp$Age + tmp$Vintage.Date - 1]
  tmp$hazard <- log(hscale*baseline.haz[[ix]][tmp$Age])
  #tmp$Age <- factor(tmp$Age)
  mm12 <- model.matrix(~ -1 + hazard + x1 + x2, data = tmp)
  kappa <- kappa(mm12, exact = TRUE)
  
  coef <- matrix(c(c1$coef[1], c1$coef[2]),
                 ncol=2, byrow=TRUE)
  
  # ht <- paste0("tt(Loan.ID)=", beta[1])
  # pv <- matrix(c(linearHypothesis(c1, ht)$Pr[2],NA,NA),
  #              ncol=3, byrow=TRUE)
  
  
  data.frame(seed, dr, beta1 = beta[1], beta2 = beta[2], coef, `kappa` = kappa, hscale = exp(h.opt$minimum), E_t, G_v)
}

delta_set <- seq(0, 0.15, by = 0.01)

# start the clock
ptm <- proc.time()

#func_set <- list()
#finc <- 1

for (delta_t in delta_set) {
  
  full.output <- list()
  res <- as.data.table(expand.grid(Delta_t = delta_t, Delta_v = delta_set))
  hscale.prev <- 1
  
  for (delta_v in delta_set) {
    
    print(paste0("Delta_t=", delta_t, ", Delta_v=", delta_v))
    
    for (ix in 1:n_sims) {
      R1 <- survival_sim_new_effect3A(ix = ix, n = n, seed = seed + ix, t_max = t_max, hscale.prev=hscale.prev, #scale = scale,
                                     afun = ln_baseline, tfun = tfun60, vfun = vfun_sar1, 
                                     cp = cp, beta = beta, delta_t = delta_t, delta_v = delta_v, dr_only=FALSE)
      
      hscale.prev <- R1$hscale
      if  (ix==1) R2 <- R1
      else R2 <- rbind(R2,R1)
      
      #      func_set[[finc]] <- list(delta_t, delta_v, E_t=R1$E_t, G_V=R1$G_v)
      #      finc <- finc+1
    }
    full.output[[paste0("delta_t=", delta_t, ", delta_v=", delta_v)]] <- R2
    
    res[Delta_t == delta_t & Delta_v == delta_v, c("Average absolute error for X1", "Average error for X1",
                                                   "Average absolute error for X2", "Average error for X2",
                                                   "Kappa") :=
          list(mean(abs(R2$X1 - beta[1]))/beta[1], mean((R2$X1 - beta[1]))/beta[1], 
               mean(abs(R2$X2 - beta[2]))/beta[2], mean((R2$X2 - beta[2]))/beta[2],
               mean(R2$kappa))]
    
  }
  save(res, full.output, file = paste0("effect3A_delta_t_", delta_t, ".RData"))
}

#save(func_set, file="func_set.RData")

# end the clock
print(proc.time() - ptm)




# Effect  #3B -----------------------------------------------------------------
# The approximate time of calculations  - about 115 hours which is about 4.8 days. 

# Dumping the initialization in the global space. Not good, but quick
# Doing this to make sure that each change in parameter settings tests the same set of functions

d <- list()
baseline.haz <- list()
f_t <- list()
f_v <- list()
pf <- list()
pa <- list()

set.seed(seed)  # Reset generator to the global seed

for (ix in 1:n_sims) {
  # Simulate the start date of each observation (origination date)
  d[[ix]] <- ceiling(runif(n)*t_max)
  
  # -------------------------------------
  # Simulate TVC (Time-varying covariate):
  
  baseline.haz[[ix]] <- afun(1:t_max)
  
  # Calendar-time effect extending beyond data series for extrapolation
  f_t[[ix]] <- tfun(t_max)
  
  # Vintage effect over each year
  f_v[[ix]] <- vfun(t_max)

  # Threshold probablity for failure for each observation
  pf[[ix]] <- runif(n)
  
  # Set the thresholds for probability of attrition
  pa[[ix]] <- matrix(runif(t_max*n), ncol=t_max)
}


survival_sim_new_effect3B <- function(ix, n, seed, t_max, scale, afun, vfun, tfun, cp, beta, hscale.prev,
                                     delta_t, delta_v, dr_only = FALSE){
  
#  set.seed(seed)  # Removed, because everything here is deterministic now.
  
  # E_t <- f_t[[ix]] - Proj(f_t[[ix]], rep(1, t_max)) - Proj(f_t[[ix]], 1:t_max - Proj(1:t_max, rep(1, t_max)))
  E_t <- f_t[[ix]] - mean(f_t[[ix]]) - Proj(f_t[[ix]] - mean(f_t[[ix]]), 1:t_max - mean(1:t_max))
  E_t <- (sd(f_t[[ix]]) / sd(E_t)) * E_t
  
  # G_v <- f_v[[ix]] - Proj(f_v[[ix]], rep(1, t_max)) - Proj(f_v[[ix]], 1:t_max - Proj(1:t_max, rep(1, t_max)))
  G_v <- f_v[[ix]] - mean(f_v[[ix]]) - Proj(f_v[[ix]] - mean(f_v[[ix]]), 1:t_max - mean(1:t_max))
  G_v <- (sd(f_v[[ix]]) / sd(G_v)) * G_v
  
  E_t_delta <- delta_t * (1:t_max) + E_t  
  G_v_delta <- delta_v * (1:t_max) + G_v
  
  sim3B <- function(hscale) {
    # Keep track of time to event and record whether default or censorship;
    # Initialize as censored at last possible time point;
    def <- rep(FALSE,n)
    
    # The following implies that age=1 at the date of origination
    a <- t_max-d[[ix]]+1
    
    # Compute hazard at each age a_i;
    # get survival prob and simulate default
    chaz <- rep(0,n)
    for (a_i in 1:t_max) {
      
      # First determine cases that are censored
      ic <- which(pa[[ix]][, a_i] < cp & a > a_i)
      a[ic]<-a_i
      
      # Calculate cumulative hazard
      chaz <- chaz + hscale*baseline.haz[[ix]][a_i] * exp(beta[1] * E_t_delta[d[[ix]]+a_i-1] + beta[2] * G_v_delta[d[[ix]]])
      
      # Survival is related to cumulative hazard;
      # then determine defaults
      sp <- exp(-chaz)
      ic <- which(sp < pf[[ix]] & a > a_i)
      def[ic] <- TRUE
      a[ic] <- a_i
    }
    
    return(list(data.frame(a,def),E_t_delta,G_v_delta))
  }
  
  err3B <- function(ln_hscale) {
    out <- sim3B(exp(ln_hscale))  # We search on a log scale to improve the optimization
    
    # fraction of deafults
    dr <- sum(out[[1]]$def)/length(out[[1]]$def)
    
    dr.err <- abs(dr-dr_target)
    #print(paste(exp(ln_hscale), dr, dr.err))
    
    return(log(dr.err))
  }
  
  if (hscale.prev == 1) {  # This is the default starting condition
    h.intvl <- log(c(0.0001, 100))
  } else {
    h.intvl <- log(c(hscale.prev/10, hscale.prev*10))
  }
  h.opt <- optimize(err3B, interval=h.intvl, tol=dr_tol)
  hscale <- exp(h.opt$minimum)
  out <- sim3B(hscale = hscale)
  dr <- sum(out[[1]]$def)/length(out[[1]]$def)
  #print(paste(hscale,dr))
  
  E_t_delta = out[[2]]
  G_v_delta = out[[3]]
  training.set <- data.frame(Age = out[[1]]$a, Default = out[[1]]$def, Vintage.Date = d[[ix]], 
                             x2 = G_v_delta[d[[ix]]])
  
  # Now build Cox PH models from the simulated data
  
  c1 <- coxph(Surv(Age, Default) ~ x2 + tt(Vintage.Date), data = training.set,
              tt = function(v, t,...) E_t_delta[v+t-1] )
  
  # Estimate the Condition Number of a Matrix (kappa) -------------------------
  
  tmp <- training.set[rep(1:nrow(training.set), times = training.set$Age), ]
  tmp$Age.correct<- unlist(lapply(training.set$Age, FUN = function(x) {seq(1, x, by = 1)}))
  tmp$Default[tmp$Age > tmp$Age.correct] <- FALSE
  tmp$Age <- tmp$Age.correct
  tmp$Age.correct <- NULL
  
  tmp$x1 <- E_t_delta[tmp$Age + tmp$Vintage.Date - 1]
  tmp$hazard <- log(hscale*baseline.haz[[ix]][tmp$Age])
  #tmp$Age <- factor(tmp$Age)
  mm12 <- model.matrix(~ -1 + hazard + x1 + x2, data = tmp)
  kappa <- kappa(mm12, exact = TRUE)
  
  coef <- matrix(c(c1$coef[1], c1$coef[2]),
                 ncol=2, byrow=TRUE)
  
  # ht <- paste0("tt(Loan.ID)=", beta[1])
  # pv <- matrix(c(linearHypothesis(c1, ht)$Pr[2],NA,NA),
  #              ncol=3, byrow=TRUE)
  
  
  data.frame(seed, dr, beta1 = beta[1], beta2 = beta[2], coef, `kappa` = kappa, hscale = exp(h.opt$minimum), E_t, G_v)
}

delta_set <- seq(0, 0.15, by = 0.01)

# start the clock
ptm <- proc.time()

#func_set <- list()
#finc <- 1

for (delta_t in delta_set) {
  
  full.output <- list()
  res <- as.data.table(expand.grid(Delta_t = delta_t, Delta_v = delta_set))
  hscale.prev <- 1
  
  for (delta_v in delta_set) {
    
    print(paste0("Delta_t=", delta_t, ", Delta_v=", delta_v))
    
    for (ix in 1:n_sims) {
      R1 <- survival_sim_new_effect3B(ix = ix, n = n, seed = seed + ix, t_max = t_max, hscale.prev=hscale.prev, #scale = scale,
                                     afun = ln_baseline, tfun = tfun60, vfun = vfun_sar1, 
                                     cp = cp, beta = beta, delta_t = delta_t, delta_v = delta_v, dr_only=FALSE)
      
      hscale.prev <- R1$hscale
      if  (ix==1) R2 <- R1
      else R2 <- rbind(R2,R1)
      
#      func_set[[finc]] <- list(delta_t, delta_v, E_t=R1$E_t, G_V=R1$G_v)
#      finc <- finc+1
    }
    full.output[[paste0("delta_t=", delta_t, ", delta_v=", delta_v)]] <- R2
    
    res[Delta_t == delta_t & Delta_v == delta_v, c("Average absolute error for X1", "Average error for X1",
                                                   "Average absolute error for X2", "Average error for X2",
                                                   "Kappa") :=
          list(mean(abs(R2$X1 - beta[1]))/beta[1], mean((R2$X1 - beta[1]))/beta[1], 
               mean(abs(R2$X2 - beta[2]))/beta[2], mean((R2$X2 - beta[2]))/beta[2],
               mean(R2$kappa))]
    
  }
  save(res, full.output, file = paste0("effect3B_test1_100sim_delta_t_", delta_t, ".RData"))
}

#save(func_set, file="func_set.RData")

# end the clock
print(proc.time() - ptm)

E_t <- list()
G_v <- list()

for (ix in 1:n_sims) {
  E_t[[ix]] <- f_t[[ix]] - Proj(f_t[[ix]], rep(1, t_max)) - Proj(f_t[[ix]], 1:t_max - Proj(1:t_max, rep(1, t_max)))
  E_t[[ix]] <- (sd(f_t[[ix]]) / sd(E_t[[ix]])) * E_t[[ix]]
  
  G_v[[ix]] <- f_v[[ix]] - Proj(f_v[[ix]], rep(1, t_max)) - Proj(f_v[[ix]], 1:t_max - Proj(1:t_max, rep(1, t_max)))
  G_v[[ix]] <- (sd(f_v[[ix]]) / sd(G_v[[ix]])) * G_v[[ix]]
}

E.all <- do.call(rbind.data.frame, E_t)
colnames(E.all) <- paste0("t", 1:t_max)

G.all <- do.call(rbind.data.frame, G_v)
colnames(G.all) <- paste0("v", 1:t_max)

# Effect  #3A by Genie --------------------------------------------------------
# The approximate time of calculations  - about 115 hours which is about 4.8 days. 

#Test1: 
# E_t <- f_t[[ix]] - mean(f_t[[ix]]) - Proj(f_t[[ix]] - mean(f_t[[ix]]), 1:t_max - mean(1:t_max))
# G_v <- f_v[[ix]] - mean(f_v[[ix]]) - Proj(f_v[[ix]] - mean(f_v[[ix]]), 1:t_max - mean(1:t_max))

#Test2: Test1 + 
# true lifecycle --- baseline.haz[[ix]] <- afun(1:t_max) 

# Test3: Test1 + constant baseline + fixed vfun(t_max) (called only once)

# Test4: Test3 + comment 2 lines (has not run yet)
# f_t[[ix]] <- f_t[[ix]] - f_t[[ix]][t_max/2]
# f_v[[ix]] <- f_v[[ix]] - f_v[[ix]][t_max/2]

# Test4: switch E_t with G_v
# E_t_delta <- delta_t * (1:t_max) + G_v   
# G_v_delta <- delta_v * (1:t_max) + E_t

# Test5: Test3 + the following
# E_t_delta <- delta_t * (t_max:1) + E_t
# E_t <- f_t[[ix]] - mean(f_t[[ix]]) - Proj(f_t[[ix]] - mean(f_t[[ix]]), t_max:1 - mean(t_max:1))

# Test6: Test5 + true lifecycle  ( baseline.haz[[ix]] <- afun(1:t_max) )

# Dumping the initialization in the global space. Not good, but quick
# Doing this to make sure that each change in parameter settings tests the same set of functions

d <- list()
baseline.haz <- list()
f_t <- list()
f_v <- list()
pf <- list()
pa <- list()

set.seed(seed)  # Reset generator to the global seed

for (ix in 1:n_sims) {
  # Simulate the start date of each observation (origination date)
  d[[ix]] <- ceiling(runif(n)*t_max)
  
  # -------------------------------------
  # Simulate TVC (Time-varying covariate):
  
  baseline.haz[[ix]] <-  rep(0.0025, times = t_max) #  afun(1:t_max)
  
  credit.risk <- vfun(t_max)
  
  # Calendar-time effect extending beyond data series for extrapolation
  f_t[[ix]] <- rev(credit.risk)  # We're using vfun instead of tfun to create a symmetric result
  # f_t[[ix]] <- f_t[[ix]] - f_t[[ix]][t_max/2]
  
  # Vintage effect over each year
  f_v[[ix]] <- credit.risk
  # f_v[[ix]] <- f_v[[ix]] - f_v[[ix]][t_max/2]
  
  # Threshold probablity for failure for each observation
  pf[[ix]] <- runif(n)
  
  # Set the thresholds for probability of attrition
  pa[[ix]] <- matrix(runif(t_max*n), ncol=t_max)
}

# f_t.mat <- do.call(rbind, f_t)
# f_t.sd <- apply(f_t.mat, 2, sd)
# 
# f_v.mat <- do.call(rbind, f_v)
# f_v.sd <- apply(f_v.mat, 2, sd)
# 
# plot(f_v.sd, type='l')
# lines(f_t.sd, col=2)
# 
# plot(f_t[[1]], type='l', ylim=c(-2,2))
# for (ix in 2:n_sims) {
#   lines(rev(f_t[[ix]]))
# }


survival_sim_new_effect3A <- function(ix, n, seed, t_max, scale, afun, vfun, tfun, cp, beta, hscale.prev,
                                      delta_t, delta_v, dr_only = FALSE){
  
  #  set.seed(seed)  # Removed, because everything here is deterministic now.
  
  # E_t <- f_t[[ix]] - Proj(f_t[[ix]], rep(1, t_max)) - Proj(f_t[[ix]], 1:t_max - Proj(1:t_max, rep(1, t_max)))
  # E_t <- f_t[[ix]] - mean(f_t[[ix]]) - Proj(f_t[[ix]] - mean(f_t[[ix]]), 1:t_max - mean(1:t_max))
  E_t <- f_t[[ix]] - mean(f_t[[ix]]) - Proj(f_t[[ix]] - mean(f_t[[ix]]), t_max:1 - mean(t_max:1))
  E_t <- (sd(f_t[[ix]]) / sd(E_t)) * E_t
  
  # G_v <- f_v[[ix]] - Proj(f_v[[ix]], rep(1, t_max)) - Proj(f_v[[ix]], 1:t_max - Proj(1:t_max, rep(1, t_max)))
  G_v <- f_v[[ix]] - mean(f_v[[ix]]) - Proj(f_v[[ix]] - mean(f_v[[ix]]), 1:t_max - mean(1:t_max))
  G_v <- (sd(f_v[[ix]]) / sd(G_v)) * G_v
  
  # E_t_delta <- delta_t * (1:t_max) + E_t
  E_t_delta <- delta_t * (t_max:1) + E_t
  G_v_delta <- delta_v * (1:t_max) + G_v
  # E_t_delta <- delta_t * (1:t_max) + G_v   
  # G_v_delta <- delta_v * (1:t_max) + E_t
  
  # plot(E_t, type='l')
  # lines(G_v, col=2)
  # plot(E_t_delta, type='l')
  # lines(G_v_delta, col=2)
  
  sim3A <- function(hscale) {
    # Keep track of time to event and record whether default or censorship;
    # Initialize as censored at last possible time point;
    def <- rep(FALSE,n)
    
    # The following implies that age=1 at the date of origination
    a <- t_max-d[[ix]]+1
    
    # Compute hazard at each age a_i;
    # get survival prob and simulate default
    chaz <- rep(0,n)
    for (a_i in 1:t_max) {
      
      # First determine cases that are censored
      ic <- which(pa[[ix]][, a_i] < cp & a > a_i)
      a[ic]<-a_i
      
      # Calculate cumulative hazard
      chaz <- chaz + hscale*baseline.haz[[ix]][a_i] * exp(beta[1] * E_t_delta[d[[ix]]+a_i-1] + beta[2] * G_v_delta[d[[ix]]])
      
      # Survival is related to cumulative hazard;
      # then determine defaults
      sp <- exp(-chaz)
      ic <- which(sp < pf[[ix]] & a > a_i)
      def[ic] <- TRUE
      a[ic] <- a_i
    }
    
    return(list(data.frame(a,def),E_t_delta,G_v_delta))
  }
  
  err3A <- function(ln_hscale) {
    out <- sim3A(exp(ln_hscale))  # We search on a log scale to improve the optimization
    
    # fraction of deafults
    dr <- sum(out[[1]]$def)/length(out[[1]]$def)
    
    dr.err <- abs(dr-dr_target)
    #print(paste(exp(ln_hscale), dr, dr.err))
    
    return(log(dr.err))
  }
  
  if (hscale.prev == 1) {  # This is the default starting condition
    h.intvl <- log(c(0.0001, 100))
  } else {
    h.intvl <- log(c(hscale.prev/10, hscale.prev*10))
  }
  h.opt <- optimize(err3A, interval=h.intvl, tol=dr_tol)
  hscale <- exp(h.opt$minimum)
  out <- sim3A(hscale = hscale)
  dr <- sum(out[[1]]$def)/length(out[[1]]$def)
  #print(paste(hscale,dr))
  
  E_t_delta = out[[2]]
  G_v_delta = out[[3]]
  training.set <- data.frame(Age = out[[1]]$a, Default = out[[1]]$def, Vintage.Date = d[[ix]], 
                             x2 = G_v_delta[d[[ix]]])
  
  # Now build Cox PH models from the simulated data
  
  c1 <- coxph(Surv(Age, Default) ~ x2 + tt(Vintage.Date), data = training.set,
              tt = function(v, t,...) E_t_delta[v+t-1] )
  
  # Estimate the Condition Number of a Matrix (kappa) -------------------------
  
  tmp <- training.set[rep(1:nrow(training.set), times = training.set$Age), ]
  tmp$Age.correct<- unlist(lapply(training.set$Age, FUN = function(x) {seq(1, x, by = 1)}))
  tmp$Default[tmp$Age > tmp$Age.correct] <- FALSE
  tmp$Age <- tmp$Age.correct
  tmp$Age.correct <- NULL
  
  tmp$x1 <- E_t_delta[tmp$Age + tmp$Vintage.Date - 1]
  tmp$hazard <- log(hscale*baseline.haz[[ix]][tmp$Age])
  #tmp$Age <- factor(tmp$Age)
  mm12 <- model.matrix(~ -1 + hazard + x1 + x2, data = tmp)
  kappa <- kappa(mm12, exact = TRUE)
  
  coef <- matrix(c(c1$coef[1], c1$coef[2]),
                 ncol=2, byrow=TRUE)
  
  # ht <- paste0("tt(Loan.ID)=", beta[1])
  # pv <- matrix(c(linearHypothesis(c1, ht)$Pr[2],NA,NA),
  #              ncol=3, byrow=TRUE)
  
  
  data.frame(seed, dr, beta1 = beta[1], beta2 = beta[2], coef, `kappa` = kappa, hscale = exp(h.opt$minimum), E_t, G_v)
}

delta_set <- seq(0, 0.15, by = 0.01)

# start the clock
ptm <- proc.time()

#func_set <- list()
#finc <- 1

for (delta_t in delta_set) {
  
  full.output <- list()
  res <- as.data.table(expand.grid(Delta_t = delta_t, Delta_v = delta_set))
  hscale.prev <- 1
  
  for (delta_v in delta_set) {
    
    print(paste0("Delta_t=", delta_t, ", Delta_v=", delta_v))
    
    for (ix in 1:n_sims) {
      R1 <- survival_sim_new_effect3A(ix = ix, n = n, seed = seed + ix, t_max = t_max, hscale.prev=hscale.prev, #scale = scale,
                                      afun = ln_baseline, tfun = tfun60, vfun = vfun_sar1, 
                                      cp = cp, beta = beta, delta_t = delta_t, delta_v = delta_v, dr_only=FALSE)
      
      hscale.prev <- R1$hscale
      if  (ix==1) R2 <- R1
      else R2 <- rbind(R2,R1)
      
      #      func_set[[finc]] <- list(delta_t, delta_v, E_t=R1$E_t, G_V=R1$G_v)
      #      finc <- finc+1
    }
    full.output[[paste0("delta_t=", delta_t, ", delta_v=", delta_v)]] <- R2
    
    res[Delta_t == delta_t & Delta_v == delta_v, c("Average absolute error for X1", "Average error for X1",
                                                   "Average absolute error for X2", "Average error for X2",
                                                   "Kappa") :=
          list(mean(abs(R2$X1 - beta[1]))/beta[1], mean((R2$X1 - beta[1]))/beta[1], 
               mean(abs(R2$X2 - beta[2]))/beta[2], mean((R2$X2 - beta[2]))/beta[2],
               mean(R2$kappa))]
    
  }
  save(res, full.output, file = paste0("effect3A_test5_100sim_delta_t_", delta_t, ".RData"))
}

#save(func_set, file="func_set.RData")

# end the clock
print(proc.time() - ptm)

# ------- Plot from Genie -------------
E.all.ggdt <- as.data.table(E.all)
E.all.ggdt <- melt(E.all.ggdt, measure.vars = names(E.all.ggdt), variable.name = "date")
levels(E.all.ggdt$date) <- gsub("^t", "", levels(E.all.ggdt$date))
E.all.ggdt[, date := as.numeric(as.character(date))]

# 'mult' is the multiplier of the standard error
# this doesn't work
# ggplot(data=E.all.ggdt, aes(x=date, y=value)) +
#   stat_summary(fun.data ="mean_sdl", mult=1, geom = "smooth") + theme_bw()

G.all.ggdt <- as.data.table(G.all)
G.all.ggdt <- melt(G.all.ggdt, measure.vars = names(G.all.ggdt), variable.name = "date")
levels(G.all.ggdt$date) <- gsub("^v", "", levels(G.all.ggdt$date))
G.all.ggdt[, date := as.numeric(as.character(date))]

pdf("functions.pdf")
E.all.ggdt <- E.all.ggdt[, as.list(smean.sdl(value, mult = 1)), keyby = .(date)]
ggplot(data = E.all.ggdt, aes(x = date, y = Mean)) + 
  geom_line(size = 1, color = "blue") + 
  geom_ribbon(aes(ymax = Lower, ymin = Upper), alpha = 0.2)

G.all.ggdt <- G.all.ggdt[, as.list(smean.sdl(value, mult = 1)), keyby = .(date)]
ggplot(data = G.all.ggdt, aes(x = date, y = Mean)) + 
  geom_line(size = 1, color = "blue") + 
  geom_ribbon(aes(ymax = Lower, ymin = Upper), alpha = 0.2)
dev.off()


