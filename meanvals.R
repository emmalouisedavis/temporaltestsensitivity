#' Calculate temporal sensitivities
#' @author Emma Louise Davis
#' @return Mean temporal estimates of sensitivity for all models

probs_nasal_LFT <- 0
probs_throat_LFT <- 0
probs_nasal_FFA <- 0
probs_throat_FFA <- 0
probs_nasal_PCR <- 0
probs_throat_PCR <- 0
probs_nasal_PCR_ext <- 0
probs_throat_PCR_ext <- 0
probs_nasal_PCR_weak_ext <- 0
probs_throat_PCR_weak_ext <- 0

for(t in 1:20){
  probs_nasal_LFT[t+1] <- P_pos(t,alphaX_nasal_LFT,gammaX_nasal_LFT,
                                  alphaZ_nasal_LFT,gammaZ_nasal_LFT,sens_nasal_LFT)
  probs_throat_LFT[t+1] <- P_pos(t,alphaX_throat_LFT,gammaX_throat_LFT,
                                   alphaZ_throat_LFT,gammaZ_throat_LFT,sens_throat_LFT)
  probs_nasal_PCR[t+1] <- P_pos(t,alphaX_nasal_PCR,gammaX_nasal_PCR,
                                  alphaZ_nasal_PCR,gammaZ_nasal_PCR,sens_nasal_PCR)
  probs_throat_PCR[t+1] <- P_pos(t,alphaX_throat_PCR,gammaX_throat_PCR,
                                   alphaZ_throat_PCR,gammaZ_throat_PCR,sens_throat_PCR)
  probs_nasal_PCR_ext[t+1] <- P_pos(t,alphaX_nasal_PCR,gammaX_nasal_PCR,
                                      alphaZ_nasal_PCR_ext,gammaZ_nasal_PCR_ext,sens_nasal_PCR)
  probs_throat_PCR_ext[t+1] <- P_pos(t,alphaX_throat_PCR,gammaX_throat_PCR,
                                       alphaZ_throat_PCR_ext,gammaZ_throat_PCR_ext,sens_throat_PCR)
  probs_nasal_PCR_weak_ext[t+1] <- P_pos(t,alphaX_nasal_PCR_weak,gammaX_nasal_PCR_weak,
                                           alphaZ_nasal_PCR_weak_ext,gammaZ_nasal_PCR_weak_ext,sens_nasal_PCR_weak)
  probs_throat_PCR_weak_ext[t+1] <- P_pos(t,alphaX_throat_PCR_weak,gammaX_throat_PCR_weak,
                                            alphaZ_throat_PCR_weak_ext,gammaZ_throat_PCR_weak_ext,sens_throat_PCR_weak)
  probs_nasal_FFA[t+1] <- P_pos(t,alphaX_nasal_FFA,gammaX_nasal_FFA,
                                  alphaZ_nasal_FFA,gammaZ_nasal_FFA,sens_nasal_FFA)
  probs_throat_FFA[t+1] <- P_pos(t,alphaX_throat_FFA,gammaX_throat_FFA,
                                   alphaZ_throat_FFA,gammaZ_throat_FFA,sens_throat_FFA)
  
}