#' Functions for testing model
#' @author Emma Louise Davis
#' @return Functions for use in temporal test sensitivity analysis

########
# Survival function, inputs alpha and gamma
# alpha <- exp(-intercept/scale))
# gamma <- 1/scale
# intercept and scale from survival::survreg, dist = "loglogistic"
# S(t) = Prob[T>t]

S_t <- function(t,alpha,gamma){
  out <- 1/(1+alpha*(t)^gamma)
  return(out)
}


########
# Calculating confidence intervals for binomial proportions

Binom_conf <- function(x,n,conf.level=0.95,lower=T){
  out <- c()
  if(length(n)==1) n <- rep(n,length(x)) 
  for(i in seq_along(x)){
    binomfit <- binom.test(x[i],n[i],alternative="two.sided",conf.level = conf.level)
    if(lower) out[i] <- binomfit$conf.int[1]
    if(!lower) out[i] <- binomfit$conf.int[2]
  }
  return(out)
}


###########################################################################
#
# Model Functions
# X = Time to first positive test (from infection)
# Z = X + Y = Time to last positive test (from infection)
# X < Z
#
###########################################################################

########
# Probability t = first test time
# P[X=t]

P_equals <- function(t,alpha,gamma){
  return(S_t(t-1,alpha,gamma) - S_t(t,alpha,gamma))
}

#######
# Probability t is between first and last test time
# P[X<t<X+Y]

P_between <- function(t,alpha_x,gamma_x,alpha_z,gamma_z){
  if(t %in% c(0,1)) return(0)
  if(t>1) {
    x <- 1:(t-1)
    seq_x <- P_equals(x,alpha_x,gamma_x)*S_t(t,alpha_z,gamma_z)/S_t(x,alpha_z,gamma_z)
    return(sum(seq_x))
  }
}

########
# Probability t = last test time
# P[Z=t]

P_last <- function(t,alpha_z,gamma_z){
  if(t %in% c(0,1)) return(0)
  if(t>1) {
    x <- 1:(t-1)
    seq_x <- P_equals(t,alpha_z,gamma_z)*P_equals(x,alpha_z,gamma_z)
    return(sum(seq_x))
  }
}

########
# Probability test result is positive at time t
# = P[X=t] + sens*P[X<t<Z] + P[Z=t]

P_pos <- function(t,alpha_x,gamma_x,alpha_z,gamma_z,sens){
  prob <- P_equals(t,alpha_z,gamma_z) + 
    sens*P_between(t,alpha_x,gamma_x,alpha_z,gamma_z) +
    P_last(t,alpha_z,gamma_z)
  return(prob)
}

###########################################################################
#
# Functions used in bootstrapping
#
###########################################################################

########
# Bootstrapping sampler

pickIDs <- function(N){
  return(sample(N,N,replace=T))
}

########
# Select data for bootstrapping sample

sampleData <- function(df,IDs){
  return(df[IDs,])
}

########
# Fit model

modelFit <- function(df){
  modelX <- survival::survreg(Surv(time_first, status_first) ~ 1,
                             data = df,
                             dist = "loglogistic"
  )
  modelZ <- survival::survreg(Surv(time_last, status_last) ~ 1,
                                   data = df,
                                   dist = "loglogistic"
  )
  return(list(modelX,modelZ))
}

########
# Get parameters from modelFit output

getModelParams <- function(model){
  intercept <- model$coeff[1]
  scale <- model$scale
  alpha <- as.numeric(exp(-intercept/scale))
  gamma <- 1/scale
  return(c(alpha,gamma))
}
