#' Case-resampling percentile bootstrap
#' @author Emma Louise Davis
#' @return Temporal sensitivity confidence intervals for all models

set.seed(20220412)
sims <- 10000
days <- 20
probs_nasal_LFT <- matrix(0,days+1,sims)
probs_throat_LFT <- matrix(0,days+1,sims)
probs_nasal_FFA <- matrix(0,days+1,sims)
probs_throat_FFA <- matrix(0,days+1,sims)
probs_nasal_PCR <- matrix(0,days+1,sims)
probs_throat_PCR <- matrix(0,days+1,sims)
probs_nasal_PCR_ext <- matrix(0,days+1,sims)
probs_throat_PCR_ext <- matrix(0,days+1,sims)
probs_nasal_PCR_weak_ext <- matrix(0,days+1,sims)
probs_throat_PCR_weak_ext <- matrix(0,days+1,sims)

for(i in 1:sims){  
  
  IDs <- pickIDs(N=18)
  
  sample_nasal_LFT <- sampleData(nasal_LFT,IDs)
  models <- modelFit(sample_nasal_LFT)
  paramsX <- getModelParams(models[[1]])
  paramsZ <- getModelParams(models[[2]])
  sensitivities <- binom.test(sum(num_pos_nasal_LFT[IDs]),sum(sample_nasal_LFT$duration),alternative="two.sided")
  sens_nasal_LFT_sim <- sensitivities$estimate
  for(t in 1:days){
    probs_nasal_LFT[t+1,i] <- P_pos(t,paramsX[1],paramsX[2],paramsZ[1],paramsZ[2],sens_nasal_LFT_sim)
  }
  
  sample_throat_LFT <- sampleData(throat_LFT,IDs)
  models <- modelFit(sample_throat_LFT)
  paramsX <- getModelParams(models[[1]])
  paramsZ <- getModelParams(models[[2]])
  sensitivities <- binom.test(sum(num_pos_throat_LFT[IDs]),sum(sample_throat_LFT$duration),alternative="two.sided")
  sens_throat_LFT_sim <- sensitivities$estimate
  for(t in 1:days){
    probs_throat_LFT[t+1,i] <- P_pos(t,paramsX[1],paramsX[2],paramsZ[1],paramsZ[2],sens_throat_LFT_sim)
  }
  
  sample_nasal_FFA <- sampleData(nasal_FFA,IDs)
  models <- modelFit(sample_nasal_FFA)
  paramsX <- getModelParams(models[[1]])
  paramsZ <- getModelParams(models[[2]])
  sensitivities <- binom.test(sum(num_pos_nasal_FFA[IDs]),sum(sample_nasal_FFA$duration),alternative="two.sided")
  sens_nasal_FFA_sim <- sensitivities$estimate
  for(t in 1:days){  
    probs_nasal_FFA[t+1,i] <- P_pos(t,paramsX[1],paramsX[2],paramsZ[1],paramsZ[2],sens_nasal_FFA_sim)
  }
  
  sample_throat_FFA <- sampleData(throat_FFA,IDs)
  models <- modelFit(sample_throat_FFA)
  paramsX <- getModelParams(models[[1]])
  paramsZ <- getModelParams(models[[2]])
  sensitivities <- binom.test(sum(num_pos_throat_FFA[IDs]),sum(sample_throat_FFA$duration),alternative="two.sided")
  sens_throat_FFA_sim <- sensitivities$estimate
  for(t in 1:days){
    probs_throat_FFA[t+1,i] <- P_pos(t,paramsX[1],paramsX[2],paramsZ[1],paramsZ[2],sens_throat_FFA_sim)
  }
  
  sample_nasal_PCR <- sampleData(nasal_PCR,IDs)
  models <- modelFit(sample_nasal_PCR)
  paramsX <- getModelParams(models[[1]])
  paramsZ <- getModelParams(models[[2]])
  sensitivities <- binom.test(sum(num_pos_nasal_PCR[IDs]),sum(sample_nasal_PCR$duration),alternative="two.sided")
  sens_nasal_PCR_sim <- sensitivities$estimate
  for(t in 1:days){
    probs_nasal_PCR[t+1,i] <- P_pos(t,paramsX[1],paramsX[2],paramsZ[1],paramsZ[2],sens_nasal_PCR_sim)
  }
  
  sample_throat_PCR <- sampleData(throat_PCR,IDs)
  models <- modelFit(sample_throat_PCR)
  paramsX <- getModelParams(models[[1]])
  paramsZ <- getModelParams(models[[2]])
  sensitivities <- binom.test(sum(num_pos_throat_PCR[IDs]),sum(sample_throat_PCR$duration),alternative="two.sided")
  sens_throat_PCR_sim <- sensitivities$estimate
  for(t in 1:days){
    probs_throat_PCR[t+1,i] <- P_pos(t,paramsX[1],paramsX[2],paramsZ[1],paramsZ[2],sens_throat_PCR_sim)
  }
  
  sample_nasal_PCR_ext <- sampleData(nasal_PCR_ext,IDs)
  models <- modelFit(sample_nasal_PCR_ext)
  paramsX <- getModelParams(models[[1]])
  paramsZ <- getModelParams(models[[2]])
  sensitivities <- binom.test(sum(num_pos_nasal_PCR[IDs]),sum(sample_nasal_PCR$duration),alternative="two.sided")
  sens_nasal_PCR_sim <- sensitivities$estimate
  for(t in 1:days){
    probs_nasal_PCR_ext[t+1,i] <- P_pos(t,paramsX[1],paramsX[2],paramsZ[1],paramsZ[2],sens_nasal_PCR_sim)
  }
  
  sample_throat_PCR_ext <- sampleData(throat_PCR_ext,IDs)
  models <- modelFit(sample_throat_PCR_ext)
  paramsX <- getModelParams(models[[1]])
  paramsZ <- getModelParams(models[[2]])
  sensitivities <- binom.test(sum(num_pos_throat_PCR[IDs]),sum(sample_throat_PCR$duration),alternative="two.sided")
  sens_throat_PCR_sim <- sensitivities$estimate
  for(t in 1:days){
    probs_throat_PCR_ext[t+1,i] <- P_pos(t,paramsX[1],paramsX[2],paramsZ[1],paramsZ[2],sens_throat_PCR_sim)
  }
  
  sample_nasal_PCR_weak_ext <- sampleData(nasal_PCR_weak_ext,IDs)
  sample_nasal_PCR_weak <- sampleData(nasal_PCR_weak,IDs)
  models <- modelFit(sample_nasal_PCR_weak_ext)
  paramsX <- getModelParams(models[[1]])
  paramsZ <- getModelParams(models[[2]])
  sensitivities <- binom.test(sum(num_pos_nasal_PCR_weak[IDs]),sum(sample_nasal_PCR_weak$duration),alternative="two.sided")
  sens_nasal_PCR_sim <- sensitivities$estimate
  for(t in 1:days){
    probs_nasal_PCR_weak_ext[t+1,i] <- P_pos(t,paramsX[1],paramsX[2],paramsZ[1],paramsZ[2],sens_nasal_PCR_sim)
  }
  
  sample_throat_PCR_weak_ext <- sampleData(throat_PCR_weak_ext,IDs)
  sample_throat_PCR_weak <- sampleData(throat_PCR_weak,IDs)
  models <- modelFit(sample_throat_PCR_weak_ext)
  paramsX <- getModelParams(models[[1]])
  paramsZ <- getModelParams(models[[2]])
  sensitivities <- binom.test(sum(num_pos_throat_PCR_weak[IDs]),sum(sample_throat_PCR_weak$duration),alternative="two.sided")
  sens_throat_PCR_sim <- sensitivities$estimate
  for(t in 1:days){
    probs_throat_PCR_weak_ext[t+1,i] <- P_pos(t,paramsX[1],paramsX[2],paramsZ[1],paramsZ[2],sens_throat_PCR_sim)
  }
  
}

lower_nasal_LFT <- c()
upper_nasal_LFT <- c()
lower_nasal_FFA <- c()
upper_nasal_FFA <- c()
lower_nasal_PCR <- c()
upper_nasal_PCR <- c()
lower_nasal_PCR_ext <- c()
upper_nasal_PCR_ext <- c()
lower_nasal_PCR_weak_ext <- c()
upper_nasal_PCR_weak_ext <- c()
lower_throat_LFT <- c()
upper_throat_LFT <- c()
lower_throat_FFA <- c()
upper_throat_FFA <- c()
lower_throat_PCR <- c()
upper_throat_PCR <- c()
lower_throat_PCR_ext <- c()
upper_throat_PCR_ext <- c()
lower_throat_PCR_weak_ext <- c()
upper_throat_PCR_weak_ext <- c()
for(t in 1:(days+1)){
  lower_nasal_LFT[t] <- quantile(probs_nasal_LFT[t,],0.025)
  upper_nasal_LFT[t] <- quantile(probs_nasal_LFT[t,],0.975)
  lower_nasal_FFA[t] <- quantile(probs_nasal_FFA[t,],0.025)
  upper_nasal_FFA[t] <- quantile(probs_nasal_FFA[t,],0.975)
  lower_nasal_PCR[t] <- quantile(probs_nasal_PCR[t,],0.025,na.rm=T)
  upper_nasal_PCR[t] <- quantile(probs_nasal_PCR[t,],0.975,na.rm=T)
  lower_nasal_PCR_ext[t] <- quantile(probs_nasal_PCR_ext[t,],0.025,na.rm=T)
  upper_nasal_PCR_ext[t] <- quantile(probs_nasal_PCR_ext[t,],0.975,na.rm=T)
  lower_nasal_PCR_weak_ext[t] <- quantile(probs_nasal_PCR_weak_ext[t,],0.025,na.rm=T)
  upper_nasal_PCR_weak_ext[t] <- quantile(probs_nasal_PCR_weak_ext[t,],0.975,na.rm=T)
  lower_throat_LFT[t] <- quantile(probs_throat_LFT[t,],0.025)
  upper_throat_LFT[t] <- quantile(probs_throat_LFT[t,],0.975)
  lower_throat_FFA[t] <- quantile(probs_throat_FFA[t,],0.025)
  upper_throat_FFA[t] <- quantile(probs_throat_FFA[t,],0.975)
  lower_throat_PCR[t] <- quantile(probs_throat_PCR[t,],0.025,na.rm=T)
  upper_throat_PCR[t] <- quantile(probs_throat_PCR[t,],0.975,na.rm=T)
  lower_throat_PCR_ext[t] <- quantile(probs_throat_PCR_ext[t,],0.025,na.rm=T)
  upper_throat_PCR_ext[t] <- quantile(probs_throat_PCR_ext[t,],0.975,na.rm=T)
  lower_throat_PCR_weak_ext[t] <- quantile(probs_throat_PCR_weak_ext[t,],0.025,na.rm=T)
  upper_throat_PCR_weak_ext[t] <- quantile(probs_throat_PCR_weak_ext[t,],0.975,na.rm=T)
}

# plot(0:days,upper_nasal_LFT,type="l")
# lines(0:days,lower_nasal_LFT)
# plot(0:days,upper_throat_LFT,type="l")
# lines(0:days,lower_throat_LFT)
# 
# plot(0:days,upper_nasal_FFA,type="l")
# lines(0:days,lower_nasal_FFA)
# plot(0:days,upper_throat_FFA,type="l")
# lines(0:days,lower_throat_FFA)
# 
# plot(0:days,upper_nasal_PCR,type="l")
# lines(0:days,lower_nasal_PCR)
# plot(0:days,upper_throat_PCR,type="l")
# lines(0:days,lower_nasal_PCR)

