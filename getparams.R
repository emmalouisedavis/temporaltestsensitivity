#' Load data and get parameters
#' @author Emma Louise Davis
#' @return parameters for test sensitivity model

###########################################################################
#
# Original data: Human Challenge Study Supplementary Figure 1
# Days 1-14
# 18 participants
#
###########################################################################

testing <- as_tibble(read.csv('HumanChallengeDailyTests_noNAs.csv'))

# Index of patients (18 total) whom:
# treatment = TRUE (=1) or FALSE (=0)
treatment <- c(rep(1,6),rep(0,12))
# symptoms = TRUE (=1) or FALSE (=0)
symptoms <- c(rep(1,13),0,0,1,1,0)

# Lateral flow testing (positive=1, negative=0)
LFTtesting <- testing %>% filter(Assay == "LFA")
# FFA testing (strong positive=2, weak positive=1, negative=0)
FFAtesting <- testing %>% filter(Assay == "FFA")

# Version of FFA dataset with strong positive=1 and weak positive=negative=0
FFAtesting_strong <- FFAtesting %>% mutate(Day_1 = pmax(Day_1-1,0),
                                           Day_2 = pmax(Day_2-1,0),
                                           Day_3 = pmax(Day_3-1,0),
                                           Day_4 = pmax(Day_4-1,0),
                                           Day_5 = pmax(Day_5-1,0),
                                           Day_6 = pmax(Day_6-1,0),
                                           Day_7 = pmax(Day_7-1,0),
                                           Day_8 = pmax(Day_8-1,0),
                                           Day_9 = pmax(Day_9-1,0),
                                           Day_10 = pmax(Day_10-1,0),
                                           Day_11 = pmax(Day_11-1,0),
                                           Day_12 = pmax(Day_12-1,0),
                                           Day_13 = pmax(Day_13-1,0),
                                           Day_14 = pmax(Day_14-1,0))
# Version of FFA dataset with strong positive=weak positive=1, negative=0)
FFAtesting <- FFAtesting
FFAtesting$Day_1 <- FFAtesting$Day_1 - FFAtesting_strong$Day_1
FFAtesting$Day_2 <- FFAtesting$Day_2 - FFAtesting_strong$Day_2
FFAtesting$Day_3 <- FFAtesting$Day_3 - FFAtesting_strong$Day_3
FFAtesting$Day_4 <- FFAtesting$Day_4 - FFAtesting_strong$Day_4
FFAtesting$Day_5 <- FFAtesting$Day_5 - FFAtesting_strong$Day_5
FFAtesting$Day_6 <- FFAtesting$Day_6 - FFAtesting_strong$Day_6
FFAtesting$Day_7 <- FFAtesting$Day_7 - FFAtesting_strong$Day_7
FFAtesting$Day_8 <- FFAtesting$Day_8 - FFAtesting_strong$Day_8
FFAtesting$Day_9 <- FFAtesting$Day_9 - FFAtesting_strong$Day_9
FFAtesting$Day_10 <- FFAtesting$Day_10 - FFAtesting_strong$Day_10
FFAtesting$Day_11 <- FFAtesting$Day_11 - FFAtesting_strong$Day_11
FFAtesting$Day_12 <- FFAtesting$Day_12 - FFAtesting_strong$Day_12
FFAtesting$Day_13 <- FFAtesting$Day_13 - FFAtesting_strong$Day_13
FFAtesting$Day_14 <- FFAtesting$Day_14 - FFAtesting_strong$Day_14

# qPCR testing (strong positive=2, weak positive=1, negative=0)
PCRtesting <- testing %>% filter(Assay == "qPCR")

# Version of qPCR dataset with strong positive=1 and weak positive=negative=0
PCRtesting_strong <- PCRtesting %>% mutate(Day_1 = pmax(Day_1-1,0),
                                           Day_2 = pmax(Day_2-1,0),
                                           Day_3 = pmax(Day_3-1,0),
                                           Day_4 = pmax(Day_4-1,0),
                                           Day_5 = pmax(Day_5-1,0),
                                           Day_6 = pmax(Day_6-1,0),
                                           Day_7 = pmax(Day_7-1,0),
                                           Day_8 = pmax(Day_8-1,0),
                                           Day_9 = pmax(Day_9-1,0),
                                           Day_10 = pmax(Day_10-1,0),
                                           Day_11 = pmax(Day_11-1,0),
                                           Day_12 = pmax(Day_12-1,0),
                                           Day_13 = pmax(Day_13-1,0),
                                           Day_14 = pmax(Day_14-1,0))

# Version of qPCR dataset with strong positive=weak positive=1, negative=0)
PCRtesting_weak <- PCRtesting
PCRtesting_weak$Day_1 <- PCRtesting$Day_1 - PCRtesting_strong$Day_1
PCRtesting_weak$Day_2 <- PCRtesting$Day_2 - PCRtesting_strong$Day_2
PCRtesting_weak$Day_3 <- PCRtesting$Day_3 - PCRtesting_strong$Day_3
PCRtesting_weak$Day_4 <- PCRtesting$Day_4 - PCRtesting_strong$Day_4
PCRtesting_weak$Day_5 <- PCRtesting$Day_5 - PCRtesting_strong$Day_5
PCRtesting_weak$Day_6 <- PCRtesting$Day_6 - PCRtesting_strong$Day_6
PCRtesting_weak$Day_7 <- PCRtesting$Day_7 - PCRtesting_strong$Day_7
PCRtesting_weak$Day_8 <- PCRtesting$Day_8 - PCRtesting_strong$Day_8
PCRtesting_weak$Day_9 <- PCRtesting$Day_9 - PCRtesting_strong$Day_9
PCRtesting_weak$Day_10 <- PCRtesting$Day_10 - PCRtesting_strong$Day_10
PCRtesting_weak$Day_11 <- PCRtesting$Day_11 - PCRtesting_strong$Day_11
PCRtesting_weak$Day_12 <- PCRtesting$Day_12 - PCRtesting_strong$Day_12
PCRtesting_weak$Day_13 <- PCRtesting$Day_13 - PCRtesting_strong$Day_13
PCRtesting_weak$Day_14 <- PCRtesting$Day_14 - PCRtesting_strong$Day_14

###########################################################################
#
# Transcribed data from Human Challenge Study Supplementary Figure 1
# Values manually taken from data file: 'Data/HumanChallengeDailyTests_noNAs.csv'
#
###########################################################################

##################
# LFT: times of first positive for each patient (days since exposure)

time_firstpos_nasal_LFT <- c(4,3,8,3,3,4,6,6,4,2,3,5,3,3,4,4,5,4)
time_firstpos_throat_LFT <- c(4,3,4,3,3,3,4,4,4,3,4,5,4,4,7,8,2,3)
time_lastpos_nasal_LFT <-  c(12,12,13,12,8,8,11,12,10,11,10, 9,14,10,13,11, 9,10) # from infection
time_lastpos_throat_LFT <- c(12,12,13,12,7,8,10,11, 9, 8, 8,10,14, 9, 7, 8, 6, 9) # from infection
duration_nasal_LFT <- pmax(0,time_lastpos_nasal_LFT - time_firstpos_nasal_LFT - 1) # between first positive and last positive
duration_throat_LFT <- pmax(0,time_lastpos_throat_LFT - time_firstpos_throat_LFT - 1) # between first positive and last positive

# first positive test status (for survival analysis)
# If truncated (i.e. never positive) = 0
# If not truncated (i.e. at least one positive recorded) = 1
status_LFT_firstpos <- rep(1,18) # statuses are the same for throat and nasal for LFTs
# last positive test status (for survival analysis)
# If truncated (i.e. still positive at end of 14 days) = 0
# If not truncated (i.e. negative by the end of 14 days) = 1
status_LFT_lastpos <- status_LFT_firstpos
status_LFT_lastpos[13] <- 0 # Patient 13 still LFT positive (nasal and throat) at last sample


# FINAL LFT DATASETS FOR USE IN ANALYSIS:

nasal_LFT <- tibble(time_first = time_firstpos_nasal_LFT,time_last = time_lastpos_nasal_LFT,
                    duration = duration_nasal_LFT,status_first = status_LFT_firstpos,
                    status_last = status_LFT_lastpos,treatment,symptoms)
throat_LFT <- tibble(time_first = time_firstpos_throat_LFT,time_last = time_lastpos_throat_LFT,
                     duration = duration_throat_LFT,status_first = status_LFT_firstpos,
                     status_last = status_LFT_lastpos,treatment,symptoms)

##################
# Same process for FFA (considering all positives, not just strong positives)

time_firstpos_nasal_FFA <-  c(7,3,8,5,4,3,6,6,3,2,2,6,4,3,4,4,4,3)
time_firstpos_throat_FFA <- c(4,2,3,2,2,2,5,3,2,3,3,5,3,2,5,5,2,3)
time_lastpos_nasal_FFA <-  c(10,11,10,11,9,9,11,11,10,12,12,8,6,8,12,9,8,10) # from infection
time_lastpos_throat_FFA <- c(8, 4, 9, 10,5,9,9, 9, 10,8, 7, 6,5,4,11,5,6,8) # from infection
time_lastpos_any_FFA <- pmax(time_lastpos_nasal_FFA,time_lastpos_throat_FFA)
duration_nasal_FFA <- pmax(0,time_lastpos_nasal_FFA - time_firstpos_nasal_FFA - 1) # between first positive and last positive
duration_throat_FFA <- pmax(0,time_lastpos_throat_FFA - time_firstpos_throat_FFA - 1) # between first positive and last positive
status_FFA_firstpos <- rep(1,18) # statuses are the same for throat and nasal FFA
status_FFA_lastpos <- status_FFA_firstpos

# FINAL FFA DATASETS FOR USE IN ANALYSIS:

nasal_FFA <- tibble(time_first = time_firstpos_nasal_FFA,time_last = time_lastpos_nasal_FFA,
                    duration = duration_nasal_FFA,status_first = status_FFA_firstpos,
                    status_last = status_FFA_lastpos,treatment,symptoms)
throat_FFA <- tibble(time_first = time_firstpos_throat_FFA,time_last = time_lastpos_throat_FFA,
                     duration = duration_throat_FFA,status_first = status_FFA_firstpos,
                     status_last = status_FFA_lastpos,treatment,symptoms)

##################
# Time between first LFT positive and last FFA positive (i.e. duration of detectable + viable virus)

duration_nasal_LFTFFA <- pmax(0,time_lastpos_nasal_FFA - time_firstpos_nasal_LFT - 1) # between first LFT positive and last FFA positive
duration_throat_LFTFFA <- pmax(0,time_lastpos_throat_FFA - time_firstpos_throat_LFT - 1) # between first LFT positive and last FFA positive
duration_anyFFA_nasalLFT <- pmax(0,time_lastpos_any_FFA - time_firstpos_nasal_LFT - 1) 
duration_anyFFA_throatLFT <- pmax(0,time_lastpos_any_FFA - time_firstpos_throat_LFT - 1) 
duration_throat_LFTFFA[which(duration_throat_LFTFFA==0)] <- NA
duration_anyFFA_throatLFT[which(duration_anyFFA_throatLFT==0)] <- NA

nasal_LFTFFA <- tibble(time_first = time_firstpos_nasal_LFT,time_last = time_lastpos_nasal_FFA,
                       duration = duration_nasal_LFTFFA,status_first = status_LFT_firstpos,
                       status_last = status_FFA_lastpos,treatment,symptoms)
throat_LFTFFA <- tibble(time_first = time_firstpos_throat_LFT,time_last = time_lastpos_throat_FFA,
                        duration = duration_throat_LFTFFA,status_first = status_LFT_firstpos,
                        status_last = status_FFA_lastpos,treatment,symptoms)
anyFFA_nasalLFT <- tibble(time_first = time_firstpos_nasal_LFT,time_last = time_lastpos_any_FFA,
                          duration = duration_anyFFA_nasalLFT,status_first = status_LFT_firstpos,
                          status_last = status_FFA_lastpos,treatment,symptoms)
anyFFA_throatLFT <- tibble(time_first = time_firstpos_throat_LFT,time_last = time_lastpos_any_FFA,
                           duration = duration_anyFFA_throatLFT,status_first = status_LFT_firstpos,
                           status_last = status_FFA_lastpos,treatment,symptoms)


##################
# Same process for qPCR (considering all positives and only strong positives as two scenarios)

time_firstpos_nasal_PCR_weak <- c(4,3,4,2,4,2,5,3,2,2,2,5,1,2,3,4,3,3)
time_firstpos_nasal_PCR <- c(4,3,8,2,4,2,6,5,2,2,2,6,1,2,3,4,3,3)
time_firstpos_throat_PCR_weak <- c(2,2,2,2,2,2,3,2,2,2,2,3,2,2,3,2,1,2)
time_firstpos_throat_PCR <- c(3,2,3,2,2,2,3,3,2,2,2,3,3,2,3,2,2,2)

time_lastpos_nasal_PCR_weak <-  c(14,14,14,14,14,14,14,14,14,14,14,14,14,14,14,14,14,14) # from infection
time_lastpos_nasal_PCR <-       c(14,14,14,14,14,14,14,14,12,14,13,12,14,14,14,14,12,12) # from infection
time_lastpos_throat_PCR_weak <- c(14,14,14,14,14,11,13,14,14,13,14,14,14,14,14,14,14,14) # from infection
time_lastpos_throat_PCR <-      c(14,14,14,13,14,11,13,14,14,13,14,14,14,14,14,14,13,14) # from infection

duration_nasal_PCR <- pmax(0,time_lastpos_nasal_PCR - time_firstpos_nasal_PCR - 1) # between first positive and last positive
duration_throat_PCR <- pmax(0,time_lastpos_throat_PCR - time_firstpos_throat_PCR - 1) # between first positive and last positive
duration_nasal_PCR_weak <- pmax(0,time_lastpos_nasal_PCR_weak - time_firstpos_nasal_PCR_weak - 1) # between first positive and last positive
duration_throat_PCR_weak <- pmax(0,time_lastpos_throat_PCR_weak - time_firstpos_throat_PCR_weak - 1) # between first positive and last positive

status_nasal_PCR_firstpos <- rep(1,18) # statuses are the same for throat_PCR and nasal
status_nasal_PCR_lastpos <- rep(0,18)
status_nasal_PCR_lastpos[c(11,12,17,18)] <- 1 # considering only strong positives, weak positives treated as negative result (all patients still at least weakly positive on day 14)
status_throat_PCR_firstpos <- rep(1,18)
status_throat_PCR_lastpos <- rep(0,18)
status_throat_PCR_lastpos[c(4,6,7,10,17)] <- 1 # considering only strong positives, weak positives treated as negative result

status_nasal_PCR_weak_firstpos <- status_nasal_PCR_firstpos # statuses are the same for throat_PCR and nasal
status_nasal_PCR_weak_lastpos <- rep(0,18)
status_throat_PCR_weak_firstpos <- status_throat_PCR_firstpos # statuses are the same for throat_PCR and throat
status_throat_PCR_weak_lastpos <- rep(0,18)
status_throat_PCR_weak_lastpos[c(6,7,10)] <- 1

#####################################
# FINAL qPCR DATASETS FOR USE IN ANALYSIS (strong qPCR only, considering weak as negative)

nasal_PCR <- tibble(time_first = time_firstpos_nasal_PCR,time_last = time_lastpos_nasal_PCR,
                    duration = duration_nasal_PCR,status_first = status_nasal_PCR_firstpos,
                    status_last = status_nasal_PCR_lastpos,treatment,symptoms)
throat_PCR <- tibble(time_first = time_firstpos_throat_PCR,time_last = time_lastpos_throat_PCR,
                     duration = duration_throat_PCR,status_first = status_throat_PCR_firstpos,
                     status_last = status_throat_PCR_lastpos,treatment,symptoms)

nasal_PCR_weak <- tibble(time_first = time_firstpos_nasal_PCR_weak,time_last = time_lastpos_nasal_PCR_weak,
                    duration = duration_nasal_PCR_weak,status_first = status_nasal_PCR_weak_firstpos,
                    status_last = status_nasal_PCR_weak_lastpos,treatment,symptoms)
throat_PCR_weak <- tibble(time_first = time_firstpos_throat_PCR_weak,time_last = time_lastpos_throat_PCR_weak,
                     duration = duration_throat_PCR_weak,status_first = status_throat_PCR_weak_firstpos,
                     status_last = status_throat_PCR_weak_lastpos,treatment,symptoms)


###########################################################################
#
# Get model parameters for:
# X = Time to first positive test (from infection)
# Y = Time between first positive test and last positive test
# Z = X + Y = Time of last positive test (from infection)
#
###########################################################################

# X: LFT nasal 
nasal_LFT_modelX <- survival::survreg(Surv(time_first, status_first) ~ 1,
                                 data = nasal_LFT,
                                 dist = "loglogistic"
)

interceptX_nasal_LFT <- nasal_LFT_modelX$coeff[1]
scaleX_nasal_LFT <- nasal_LFT_modelX$scale

alphaX_nasal_LFT <- as.numeric(exp(-interceptX_nasal_LFT/scaleX_nasal_LFT))
gammaX_nasal_LFT <- 1/scaleX_nasal_LFT

# Z: LFT nasal 
nasal_LFT_modelZ <- survival::survreg(Surv(time_last, status_last) ~ 1,
                                  data = nasal_LFT,
                                  dist = "loglogistic"
)

interceptZ_nasal_LFT <- nasal_LFT_modelZ$coeff[1]
scaleZ_nasal_LFT <- nasal_LFT_modelZ$scale

alphaZ_nasal_LFT <- as.numeric(exp(-interceptZ_nasal_LFT/scaleZ_nasal_LFT))
gammaZ_nasal_LFT <- 1/scaleZ_nasal_LFT

# X: LFT throat
throat_LFT_modelX <- survival::survreg(Surv(time_first, status_first) ~ 1,
                                  data = throat_LFT,
                                  dist = "loglogistic"
)

interceptX_throat_LFT <- throat_LFT_modelX$coeff[1]
scaleX_throat_LFT <- throat_LFT_modelX$scale

alphaX_throat_LFT <- as.numeric(exp(-interceptX_throat_LFT/scaleX_throat_LFT))
gammaX_throat_LFT <- 1/scaleX_throat_LFT

# Z: LFT throat
throat_LFT_modelZ <- survival::survreg(Surv(time_last, status_last) ~ 1,
                                  data = throat_LFT,
                                  dist = "loglogistic"
)

interceptZ_throat_LFT <- throat_LFT_modelZ$coeff[1]
scaleZ_throat_LFT <- throat_LFT_modelZ$scale

alphaZ_throat_LFT <- as.numeric(exp(-interceptZ_throat_LFT/scaleZ_throat_LFT))
gammaZ_throat_LFT <- 1/scaleZ_throat_LFT

# X: FFA nasal (all positives, not just strong positives)
nasal_FFA_modelX <- survival::survreg(Surv(time_first, status_first) ~ 1,
                                      data = nasal_FFA,
                                      dist = "loglogistic"
)

interceptX_nasal_FFA <- nasal_FFA_modelX$coeff[1]
scaleX_nasal_FFA <- nasal_FFA_modelX$scale

alphaX_nasal_FFA <- as.numeric(exp(-interceptX_nasal_FFA/scaleX_nasal_FFA))
gammaX_nasal_FFA <- 1/scaleX_nasal_FFA

# Z: FFA nasal (all positives, not just strong positives)
nasal_FFA_modelZ <- survival::survreg(Surv(time_last, status_last) ~ 1,
                                      data = nasal_FFA,
                                      dist = "loglogistic"
)

interceptZ_nasal_FFA <- nasal_FFA_modelZ$coeff[1]
scaleZ_nasal_FFA <- nasal_FFA_modelZ$scale

alphaZ_nasal_FFA <- as.numeric(exp(-interceptZ_nasal_FFA/scaleZ_nasal_FFA))
gammaZ_nasal_FFA <- 1/scaleZ_nasal_FFA

# X: FFA throat (all positives, not just strong positives)
throat_FFA_modelX <- survival::survreg(Surv(time_first, status_first) ~ 1,
                                       data = throat_FFA,
                                       dist = "loglogistic"
)

interceptX_throat_FFA <- throat_FFA_modelX$coeff[1]
scaleX_throat_FFA <- throat_FFA_modelX$scale

alphaX_throat_FFA <- as.numeric(exp(-interceptX_throat_FFA/scaleX_throat_FFA))
gammaX_throat_FFA <- 1/scaleX_throat_FFA

# Z: FFA throat (all positives, not just strong positives)
throat_FFA_modelZ <- survival::survreg(Surv(time_last, status_last) ~ 1,
                                       data = throat_FFA,
                                       dist = "loglogistic"
)

interceptZ_throat_FFA <- throat_FFA_modelZ$coeff[1]
scaleZ_throat_FFA <- throat_FFA_modelZ$scale

alphaZ_throat_FFA <- as.numeric(exp(-interceptZ_throat_FFA/scaleZ_throat_FFA))
gammaZ_throat_FFA <- 1/scaleZ_throat_FFA

# X: qPCR nasal (just strong positives)
nasal_PCR_modelX <- survival::survreg(Surv(time_first, status_first) ~ 1,
                                      data = nasal_PCR,
                                      dist = "loglogistic"
)

interceptX_nasal_PCR <- nasal_PCR_modelX$coeff[1]
scaleX_nasal_PCR <- nasal_PCR_modelX$scale

alphaX_nasal_PCR <- as.numeric(exp(-interceptX_nasal_PCR/scaleX_nasal_PCR))
gammaX_nasal_PCR <- 1/scaleX_nasal_PCR

# Z: qPCR nasal (just strong positives)
nasal_PCR_modelZ <- survival::survreg(Surv(time_last, status_last) ~ 1,
                                      data = nasal_PCR,
                                      dist = "loglogistic"
)

interceptZ_nasal_PCR <- nasal_PCR_modelZ$coeff[1]
scaleZ_nasal_PCR <- nasal_PCR_modelZ$scale
 
alphaZ_nasal_PCR <- as.numeric(exp(-interceptZ_nasal_PCR/scaleZ_nasal_PCR))
gammaZ_nasal_PCR <- 1/scaleZ_nasal_PCR

# X: qPCR throat (just strong positives)
throat_PCR_modelX <- survival::survreg(Surv(time_first, status_first) ~ 1,
                                       data = throat_PCR,
                                       dist = "loglogistic"
)

interceptX_throat_PCR <- throat_PCR_modelX$coeff[1]
scaleX_throat_PCR <- throat_PCR_modelX$scale

alphaX_throat_PCR <- as.numeric(exp(-interceptX_throat_PCR/scaleX_throat_PCR))
gammaX_throat_PCR <- 1/scaleX_throat_PCR

# Z: qPCR throat (just strong positives)
throat_PCR_modelZ <- survival::survreg(Surv(time_last, status_last) ~ 1,
                                       data = throat_PCR,
                                       dist = "loglogistic"
)

interceptZ_throat_PCR <- throat_PCR_modelZ$coeff[1]
scaleZ_throat_PCR <- throat_PCR_modelZ$scale

alphaZ_throat_PCR <- as.numeric(exp(-interceptZ_throat_PCR/scaleZ_throat_PCR))
gammaZ_throat_PCR <- 1/scaleZ_throat_PCR

# X: qPCR nasal (including weak)
nasal_PCR_weak_modelX <- survival::survreg(Surv(time_first, status_first) ~ 1,
                                      data = nasal_PCR_weak,
                                      dist = "loglogistic"
)

interceptX_nasal_PCR_weak <- nasal_PCR_weak_modelX$coeff[1]
scaleX_nasal_PCR_weak <- nasal_PCR_weak_modelX$scale

alphaX_nasal_PCR_weak <- as.numeric(exp(-interceptX_nasal_PCR_weak/scaleX_nasal_PCR_weak))
gammaX_nasal_PCR_weak <- 1/scaleX_nasal_PCR_weak

# X: qPCR throat (including weak)
throat_PCR_weak_modelX <- survival::survreg(Surv(time_first, status_first) ~ 1,
                                       data = throat_PCR_weak,
                                       dist = "loglogistic"
)

interceptX_throat_PCR_weak <- throat_PCR_weak_modelX$coeff[1]
scaleX_throat_PCR_weak <- throat_PCR_weak_modelX$scale

alphaX_throat_PCR_weak <- as.numeric(exp(-interceptX_throat_PCR_weak/scaleX_throat_PCR_weak))
gammaX_throat_PCR_weak <- 1/scaleX_throat_PCR_weak

# Y: first LFT to last FFA (nasal)
nasal_LFTFFA_modelY <- survival::survreg(Surv(duration, status_last) ~ 1,
                                      data = nasal_LFTFFA,
                                      dist = "loglogistic"
)

interceptY_nasal_LFTFFA <- nasal_LFTFFA_modelY$coeff[1]
scaleY_nasal_LFTFFA <- nasal_LFTFFA_modelY$scale

alphaY_nasal_LFTFFA <- as.numeric(exp(-interceptY_nasal_LFTFFA/scaleY_nasal_LFTFFA))
gammaY_nasal_LFTFFA <- 1/scaleY_nasal_LFTFFA

# Y: first LFT to last FFA (throat)
throat_LFTFFA_modelY <- survival::survreg(Surv(duration, status_last) ~ 1,
                                         data = throat_LFTFFA,
                                         dist = "loglogistic"
)

interceptY_throat_LFTFFA <- throat_LFTFFA_modelY$coeff[1]
scaleY_throat_LFTFFA <- throat_LFTFFA_modelY$scale

alphaY_throat_LFTFFA <- as.numeric(exp(-interceptY_throat_LFTFFA/scaleY_throat_LFTFFA))
gammaY_throat_LFTFFA <- 1/scaleY_throat_LFTFFA

###########################################################################
#
# Test sensitivity during detectable period
#
###########################################################################

## LFT
num_neg_nasal_LFT <- c(0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0)
num_neg_throat_LFT <- c(0,0,2,0,0,0,0,0,0,2,2,0,3,0,0,0,2,0)
num_pos_nasal_LFT <- duration_nasal_LFT - num_neg_nasal_LFT
num_pos_throat_LFT <- duration_throat_LFT - num_neg_throat_LFT

sensitivities <- binom.test(sum(num_pos_nasal_LFT),sum(duration_nasal_LFT),alternative="two.sided")
sens_nasal_LFT <- sensitivities$estimate # 95% CI (0.9316141, 0.9976397)
sens_nasal_LFT_CI <- as.numeric(sensitivities$conf.int)
sensitivities <- binom.test(sum(num_pos_throat_LFT),sum(duration_throat_LFT),alternative="two.sided")
sens_throat_LFT <- sensitivities$estimate # 95% CI (0.7802336, 0.9335901)
sens_throat_LFT_CI <- as.numeric(sensitivities$conf.int)

## FFA (all +, including weak +)
num_neg_nasal_FFA <-  c(0,0,0,0,1,3,0,1,0,1,3,0,0,0,1,0,0,1)
num_neg_throat_FFA <- c(2,0,4,1,0,3,0,0,1,2,2,0,1,0,4,0,0,0)
num_pos_nasal_FFA <- duration_nasal_FFA - num_neg_nasal_FFA
num_pos_throat_FFA <- duration_throat_FFA - num_neg_throat_FFA

sensitivities <- binom.test(sum(num_pos_nasal_FFA),sum(duration_nasal_FFA),alternative="two.sided")
sens_nasal_FFA <- sensitivities$estimate # 95% CI (0.7726394, 0.9310879)
sens_nasal_FFA_CI <- as.numeric(sensitivities$conf.int)
sensitivities <- binom.test(sum(num_pos_throat_FFA),sum(duration_throat_FFA),alternative="two.sided")
sens_throat_FFA <- sensitivities$estimate # 95% CI (0.5331273, 0.7831306)
sens_throat_FFA_CI <- as.numeric(sensitivities$conf.int)

## PCR (split in strong + only || all, including weak +)
num_neg_nasal_PCR <- c(0,0,0,0,0,1,0,0,0,0,1,0,0,0,0,0,0,0)
num_neg_throat_PCR <- c(1,0,1,0,0,0,1,0,1,1,1,0,0,3,0,3,0,1)
num_pos_nasal_PCR <- duration_nasal_PCR - num_neg_nasal_PCR
num_pos_throat_PCR <- duration_throat_PCR - num_neg_throat_PCR

num_neg_nasal_PCR_weak <-  c(0,0,2,0,0,0, 0,0,0,0,0,0,0,0,0,0,0,0)
num_neg_throat_PCR_weak <- c(1,0,0,0,0,0, 1,0,0,0,0,0,0,2,0,1,0,0)
num_pos_nasal_PCR_weak <- duration_nasal_PCR_weak - num_neg_nasal_PCR_weak
num_pos_throat_PCR_weak <- duration_throat_PCR_weak - num_neg_throat_PCR_weak

sensitivities <- binom.test(sum(num_pos_nasal_PCR),sum(duration_nasal_PCR),alternative="two.sided")
sens_nasal_PCR <- sensitivities$estimate # 95% CI (0.9563815, 0.9985106)
sens_nasal_PCR_CI <- as.numeric(sensitivities$conf.int)
sensitivities <- binom.test(sum(num_pos_throat_PCR),sum(duration_throat_PCR),alternative="two.sided")
sens_throat_PCR <- sensitivities$estimate # 95% CI (0.8840546, 0.9624674)
sens_throat_PCR_CI <- as.numeric(sensitivities$conf.int)

sensitivities <- binom.test(sum(num_pos_nasal_PCR_weak),sum(duration_nasal_PCR_weak),alternative="two.sided")
sens_nasal_PCR_weak <- sensitivities$estimate # 95% CI (0.9604412 0.9986515)
sensitivities <- binom.test(sum(num_pos_throat_PCR_weak),sum(duration_throat_PCR_weak),alternative="two.sided")
sens_throat_PCR_weak <- sensitivities$estimate # 95% CI (0.9399696 0.9914466)

###########################################################################
#
# Extended Data Human Challenge Study Supplementary Figure 1 + Extended Figures 4 and 5
# Days 1-19
# 18 participants
#
###########################################################################

testing_ext <- as_tibble(read.csv('HumanChallengeDailyTests_Extended_19days.csv'))

PCRtesting_ext <- testing_ext %>% filter(Assay == "qPCR")

# Version of qPCR dataset with strong positive=1 and weak positive=negative=0
PCRtesting_strong_ext <- PCRtesting_ext %>% mutate(Day_1 = pmax(Day_1-1,0),
                                                   Day_2 = pmax(Day_2-1,0),
                                                   Day_3 = pmax(Day_3-1,0),
                                                   Day_4 = pmax(Day_4-1,0),
                                                   Day_5 = pmax(Day_5-1,0),
                                                   Day_6 = pmax(Day_6-1,0),
                                                   Day_7 = pmax(Day_7-1,0),
                                                   Day_8 = pmax(Day_8-1,0),
                                                   Day_9 = pmax(Day_9-1,0),
                                                   Day_10 = pmax(Day_10-1,0),
                                                   Day_11 = pmax(Day_11-1,0),
                                                   Day_12 = pmax(Day_12-1,0),
                                                   Day_13 = pmax(Day_13-1,0),
                                                   Day_14 = pmax(Day_14-1,0),
                                                   Day_15 = pmax(Day_15-1,0),
                                                   Day_16 = pmax(Day_16-1,0),
                                                   Day_17 = pmax(Day_17-1,0),
                                                   Day_18 = pmax(Day_18-1,0),
                                                   Day_19 = pmax(Day_19-1,0))
# Version of qPCR dataset with strong positive=weak positive=1, negative=0)
PCRtesting_weak_ext <- PCRtesting_ext
PCRtesting_weak_ext$Day_1 <- PCRtesting_ext$Day_1 - PCRtesting_strong_ext$Day_1
PCRtesting_weak_ext$Day_2 <- PCRtesting_ext$Day_2 - PCRtesting_strong_ext$Day_2
PCRtesting_weak_ext$Day_3 <- PCRtesting_ext$Day_3 - PCRtesting_strong_ext$Day_3
PCRtesting_weak_ext$Day_4 <- PCRtesting_ext$Day_4 - PCRtesting_strong_ext$Day_4
PCRtesting_weak_ext$Day_5 <- PCRtesting_ext$Day_5 - PCRtesting_strong_ext$Day_5
PCRtesting_weak_ext$Day_6 <- PCRtesting_ext$Day_6 - PCRtesting_strong_ext$Day_6
PCRtesting_weak_ext$Day_7 <- PCRtesting_ext$Day_7 - PCRtesting_strong_ext$Day_7
PCRtesting_weak_ext$Day_8 <- PCRtesting_ext$Day_8 - PCRtesting_strong_ext$Day_8
PCRtesting_weak_ext$Day_9 <- PCRtesting_ext$Day_9 - PCRtesting_strong_ext$Day_9
PCRtesting_weak_ext$Day_10 <- PCRtesting_ext$Day_10 - PCRtesting_strong_ext$Day_10
PCRtesting_weak_ext$Day_11 <- PCRtesting_ext$Day_11 - PCRtesting_strong_ext$Day_11
PCRtesting_weak_ext$Day_12 <- PCRtesting_ext$Day_12 - PCRtesting_strong_ext$Day_12
PCRtesting_weak_ext$Day_13 <- PCRtesting_ext$Day_13 - PCRtesting_strong_ext$Day_13
PCRtesting_weak_ext$Day_14 <- PCRtesting_ext$Day_14 - PCRtesting_strong_ext$Day_14
PCRtesting_weak_ext$Day_15 <- PCRtesting_ext$Day_15 - PCRtesting_strong_ext$Day_15
PCRtesting_weak_ext$Day_16 <- PCRtesting_ext$Day_16 - PCRtesting_strong_ext$Day_16
PCRtesting_weak_ext$Day_17 <- PCRtesting_ext$Day_17 - PCRtesting_strong_ext$Day_17
PCRtesting_weak_ext$Day_18 <- PCRtesting_ext$Day_18 - PCRtesting_strong_ext$Day_18
PCRtesting_weak_ext$Day_19 <- PCRtesting_ext$Day_19 - PCRtesting_strong_ext$Day_19

##################
# Same process for extended data qPCR (considering all positives and only strong positives as two scenarios)

time_lastpos_nasal_PCR_weak_ext <-  c(14,16,17,16,17,16,15,16,18,16,14,18,18,17,15,17,15,14) # from infection
time_lastpos_nasal_PCR_ext <-       c(14,15,17,16,17,16,14,14,17,15,13,17,17,16,15,17,15,12) # from infection
time_lastpos_throat_PCR_weak_ext <- c(15,16,19,14,18,15,15,16,19,13,16,17,18,17,15,17,14,15) # from infection
time_lastpos_throat_PCR_ext <-      c(14,14,18,13,17,15,13,15,17,13,16,17,17,17,15,17,13,14) # from infection

status_nasal_PCR_lastpos_ext  <-      c(1,1,1,0,1,1,1,1,1,1,1,1,1,1,1,0,0,1)
status_nasal_PCR_weak_lastpos_ext  <- c(1,0,1,0,1,1,0,0,1,0,1,0,0,0,1,0,0,1)
status_throat_PCR_lastpos_ext <-      c(1,1,1,1,1,1,1,1,1,1,0,1,1,0,1,0,1,1)
status_throat_PCR_weak_lastpos_ext <- c(0,0,0,1,0,1,0,0,0,1,0,1,0,0,1,0,1,0)

duration_nasal_PCR_ext <- pmax(0,time_lastpos_nasal_PCR_ext - time_firstpos_nasal_PCR - 1) # between first positive and last positive
duration_throat_PCR_ext <- pmax(0,time_lastpos_throat_PCR_ext - time_firstpos_throat_PCR - 1) # between first positive and last positive
duration_nasal_PCR_weak_ext <- pmax(0,time_lastpos_nasal_PCR_weak_ext - time_firstpos_nasal_PCR_weak - 1) # between first positive and last positive
duration_throat_PCR_weak_ext <- pmax(0,time_lastpos_throat_PCR_weak_ext - time_firstpos_throat_PCR_weak - 1) # between first positive and last positive

#####################################
# FINAL qPCR DATASETS FOR USE IN ANALYSIS (strong qPCR only, considering weak as negative)

nasal_PCR_ext <- tibble(time_first = time_firstpos_nasal_PCR,time_last = time_lastpos_nasal_PCR_ext,
                        status_first = status_nasal_PCR_firstpos,status_last = status_nasal_PCR_lastpos_ext,
                        treatment,symptoms)
throat_PCR_ext <- tibble(time_first = time_firstpos_throat_PCR,time_last = time_lastpos_throat_PCR_ext,
                         status_first = status_throat_PCR_firstpos,status_last = status_throat_PCR_lastpos_ext,
                         treatment,symptoms)

nasal_PCR_weak_ext <- tibble(time_first = time_firstpos_nasal_PCR_weak,time_last = time_lastpos_nasal_PCR_weak_ext,
                             status_first = status_nasal_PCR_weak_firstpos,status_last = status_nasal_PCR_weak_lastpos_ext,
                             treatment,symptoms)
throat_PCR_weak_ext <- tibble(time_first = time_firstpos_throat_PCR_weak,time_last = time_lastpos_throat_PCR_weak_ext,
                              status_first = status_throat_PCR_weak_firstpos,status_last = status_throat_PCR_weak_lastpos_ext,
                              treatment,symptoms)



# Z: qPCR nasal (just strong positives, extended data)
nasal_PCR_modelZ_ext <- survival::survreg(Surv(time_last, status_last) ~ 1,
                                          data = nasal_PCR_ext,
                                          dist = "loglogistic"
)

interceptZ_nasal_PCR_ext <- nasal_PCR_modelZ_ext$coeff[1]
scaleZ_nasal_PCR_ext <- nasal_PCR_modelZ_ext$scale

alphaZ_nasal_PCR_ext <- as.numeric(exp(-interceptZ_nasal_PCR_ext/scaleZ_nasal_PCR_ext))
gammaZ_nasal_PCR_ext <- 1/scaleZ_nasal_PCR_ext

# Z: qPCR throat (just strong positives, extended data)
throat_PCR_modelZ_ext <- survival::survreg(Surv(time_last, status_last) ~ 1,
                                           data = throat_PCR_ext,
                                           dist = "loglogistic"
)

interceptZ_throat_PCR_ext <- throat_PCR_modelZ_ext$coeff[1]
scaleZ_throat_PCR_ext <- throat_PCR_modelZ_ext$scale

alphaZ_throat_PCR_ext <- as.numeric(exp(-interceptZ_throat_PCR_ext/scaleZ_throat_PCR_ext))
gammaZ_throat_PCR_ext <- 1/scaleZ_throat_PCR_ext

# Z: qPCR nasal (including weak, extended data)
nasal_PCR_weak_modelZ_ext <- survival::survreg(Surv(time_last, status_last) ~ 1,
                                               data = nasal_PCR_weak_ext,
                                               dist = "loglogistic"
)

interceptZ_nasal_PCR_weak_ext <- nasal_PCR_weak_modelZ_ext$coeff[1]
scaleZ_nasal_PCR_weak_ext <- nasal_PCR_weak_modelZ_ext$scale

alphaZ_nasal_PCR_weak_ext <- as.numeric(exp(-interceptZ_nasal_PCR_weak_ext/scaleZ_nasal_PCR_weak_ext))
gammaZ_nasal_PCR_weak_ext <- 1/scaleZ_nasal_PCR_weak_ext

# Z: qPCR throat (including weak, extended data)
throat_PCR_weak_modelZ_ext <- survival::survreg(Surv(time_last, status_last) ~ 1,
                                                data = throat_PCR_weak_ext,
                                                dist = "loglogistic"
)

interceptZ_throat_PCR_weak_ext <- throat_PCR_weak_modelZ_ext$coeff[1]
scaleZ_throat_PCR_weak_ext <- throat_PCR_weak_modelZ_ext$scale

alphaZ_throat_PCR_weak_ext <- as.numeric(exp(-interceptZ_throat_PCR_weak_ext/scaleZ_throat_PCR_weak_ext))
gammaZ_throat_PCR_weak_ext <- 1/scaleZ_throat_PCR_weak_ext
