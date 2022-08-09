#' Construct Figures for Estimating LFT and qPCR test-sensitivity over time since infection from a human challenge study
#' @author Emma Louise Davis
#' @return ggplot objects

# Note: Run to generate figures after runscript.R
# Figures 1 and 2 returned as p1 and p2
# Figures S1-8 returned as pS1, pS1, ... , pS8

# Colorblind friendly palette:
cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

#####################################
# Figure 1: Time to First Positive
# Discrete Kaplan-Meier curves and continuous log-logistic survival curves
#####################################

fit_nasal_LFT <- survfit(Surv(time_firstpos_nasal_LFT,status_LFT_firstpos) ~ 1,conf.type="log-log") #conf.type = log vs log-log gives similar result
summary(fit_nasal_LFT)
pos.plot_nasal_LFT <- ggsurvplot(fit_nasal_LFT,data=nasal_LFT)

fit_throat_LFT <- survfit(Surv(time_firstpos_throat_LFT,status_LFT_firstpos) ~ 1,conf.type="log-log") #conf.type = log vs log-log gives similar result
summary(fit_throat_LFT)
pos.plot_throat_LFT <- ggsurvplot(fit_throat_LFT,data=throat_LFT)

fit_nasal_FFA <- survfit(Surv(time_firstpos_nasal_FFA,status_FFA_firstpos) ~ 1,conf.type="log-log") #conf.type = log vs log-log gives similar result
summary(fit_nasal_FFA)
pos.plot_nasal_FFA <- ggsurvplot(fit_nasal_FFA,data=nasal_FFA)

fit_throat_FFA <- survfit(Surv(time_firstpos_throat_FFA,status_FFA_firstpos) ~ 1,conf.type="log-log") #conf.type = log vs log-log gives similar result
summary(fit_throat_FFA)
pos.plot_throat_FFA <- ggsurvplot(fit_throat_FFA,data=throat_FFA)

fit_nasal_PCR <- survfit(Surv(time_firstpos_nasal_PCR,status_nasal_PCR_firstpos) ~ 1,conf.type="log-log") #conf.type = log vs log-log gives similar result
summary(fit_nasal_PCR)
pos.plot_nasal_PCR <- ggsurvplot(fit_nasal_PCR,data=nasal_PCR)

fit_throat_PCR <- survfit(Surv(time_firstpos_throat_PCR,status_nasal_PCR_firstpos) ~ 1,conf.type="log-log") #conf.type = log vs log-log gives similar result
summary(fit_throat_PCR)
pos.plot_throat_PCR <- ggsurvplot(fit_throat_PCR,data=throat_PCR)

fit_nasal_PCR_weak <- survfit(Surv(time_firstpos_nasal_PCR_weak,status_nasal_PCR_weak_firstpos) ~ 1,conf.type="log-log") #conf.type = log vs log-log gives similar result
summary(fit_nasal_PCR_weak)
pos.plot_nasal_PCR_weak <- ggsurvplot(fit_nasal_PCR_weak,data=nasal_PCR_weak)

fit_throat_PCR_weak <- survfit(Surv(time_firstpos_throat_PCR_weak,status_nasal_PCR_weak_firstpos) ~ 1,conf.type="log-log") #conf.type = log vs log-log gives similar result
summary(fit_throat_PCR_weak)
pos.plot_throat_PCR_weak <- ggsurvplot(fit_throat_PCR_weak,data=throat_PCR_weak)

df <- tibble(time=seq(0,8,0.1),y=S_t(seq(0,8,0.1),alphaX_nasal_LFT,gammaX_nasal_LFT))
p1C <- pos.plot_nasal_LFT$plot + geom_line(data=df,aes(x=time,y=y)) + 
  labs(title="Time to first positive LFT (nasal)",x="Days since infection",y="Proportion still negative",tag="C") + 
  theme_minimal(base_size=14) +
  theme(legend.position="none") +
  scale_colour_manual(values=cbPalette[c(2)]) +
  scale_fill_manual(values=cbPalette[c(2)])

df <- tibble(time=seq(0,8,0.1),y=S_t(seq(0,8,0.1),alphaX_throat_LFT,gammaX_throat_LFT))
p1D <- pos.plot_throat_LFT$plot + geom_line(data=df,aes(x=time,y=y)) + 
  labs(title="Time to first positive LFT (throat)",x="Days since infection",y="Proportion still negative",tag="D") + 
  theme_minimal(base_size=14) +
  theme(legend.position="none") +
  scale_colour_manual(values=cbPalette[c(3)]) +
  scale_fill_manual(values=cbPalette[c(3)]) 

df <- tibble(time=seq(0,8,0.1),y=S_t(seq(0,8,0.1),alphaX_nasal_FFA,gammaX_nasal_FFA))
p1E <- pos.plot_nasal_FFA$plot + geom_line(data=df,aes(x=time,y=y)) + 
  labs(title="Time to first positive FFA (nasal)",x="Days since infection",y="Proportion still negative",tag="E") + 
  theme_minimal(base_size=14) +
  theme(legend.position="none") +
  scale_colour_manual(values=cbPalette[c(8)]) +
  scale_fill_manual(values=cbPalette[c(8)])

df <- tibble(time=seq(0,8,0.1),y=S_t(seq(0,8,0.1),alphaX_throat_FFA,gammaX_throat_FFA))
p1F <- pos.plot_throat_FFA$plot + geom_line(data=df,aes(x=time,y=y)) + 
  labs(title="Time to first positive FFA (throat)",x="Days since infection",y="Proportion still negative",tag="F") + 
  theme_minimal(base_size=14) +
  theme(legend.position="none") +
  scale_colour_manual(values=cbPalette[c(7)]) +
  scale_fill_manual(values=cbPalette[c(7)]) +
  coord_cartesian(xlim=c(0,8), ylim=c(0,1)) +
  scale_x_continuous(breaks=c(0,2,4,6,8))

df <- tibble(time=seq(0,8,0.1),y=S_t(seq(0,8,0.1),alphaX_nasal_PCR_weak,gammaX_nasal_PCR_weak))
p1A <- pos.plot_nasal_PCR_weak$plot + geom_line(data=df,aes(x=time,y=y)) + 
  labs(title="Time to first positive qPCR (nasal)",x="Days since infection",y="Proportion still negative", tag="A") + 
  theme_minimal(base_size=14) +
  theme(legend.position="none") +
  scale_colour_manual(values=cbPalette[c(6)]) +
  scale_fill_manual(values=cbPalette[c(6)]) + 
  coord_cartesian(xlim=c(0,8), ylim=c(0,1)) +
  scale_x_continuous(breaks=c(0,2,4,6,8))

df <- tibble(time=seq(0,8,0.1),y=S_t(seq(0,8,0.1),alphaX_throat_PCR_weak,gammaX_throat_PCR_weak))
p1B <- pos.plot_throat_PCR_weak$plot + geom_line(data=df,aes(x=time,y=y)) + 
  labs(title="Time to first positive qPCR (throat)",x="Days since infection",y="Proportion still negative",tag="B") + 
  theme_minimal(base_size=14) +
  theme(legend.position="none") +
  scale_colour_manual(values=cbPalette[c(4)]) +
  scale_fill_manual(values=cbPalette[c(4)]) +
  coord_cartesian(xlim=c(0,8), ylim=c(0,1)) +
  scale_x_continuous(breaks=c(0,2,4,6,8))

p1 <- grid.arrange(p1A,p1B,p1C,p1D,p1E,p1F,nrow=3)

#####################################
# Figure 2: Temporal sensitivity model (LFT and FFA only)
# Error windows: 95% percentile case-resampling bootstrap
# Data + 95% binom CI errorbars
#####################################

LFTtesting <- LFTtesting %>% mutate(Day_0 = 0)
FFAtesting <- FFAtesting %>% mutate(Day_0 = 0)

modelLFT <- data.frame(days=0:20,nasal=probsB_nasal_LFT,throat=probsB_throat_LFT,
                       nasal_upper=upper_nasal_LFT,nasal_lower=lower_nasal_LFT,
                       throat_upper=upper_throat_LFT,throat_lower=lower_throat_LFT)
modelFFA <- data.frame(days=0:20,nasal=probsB_nasal_FFA,throat=probsB_throat_FFA,
                       nasal_upper=upper_nasal_FFA,nasal_lower=lower_nasal_FFA,
                       throat_upper=upper_throat_FFA,throat_lower=lower_throat_FFA)

p2A <- LFTtesting %>% group_by(Site) %>%
  summarise(across(c(Day_1,Day_2,Day_3,Day_4,Day_5,Day_6,Day_7,Day_8,Day_9,Day_10,Day_11,Day_12,Day_13,Day_14), ~ mean(.x, na.rm = TRUE))) %>%
  gather(var, value, -Site) %>% 
  spread(Site, value) %>%
  arrange(match(var, c("Day_1", "Day_2", "Day_3","Day_4","Day_5","Day_6","Day_7","Day_8","Day_9","Day_10","Day_11","Day_12","Day_13","Day_14"))) %>%
  mutate(var = 1:14) %>%
  mutate(Nose_l = Binom_conf(x=Nose*18,n=18,lower=T)) %>%
  mutate(Nose_u = Binom_conf(x=Nose*18,n=18,lower=F)) %>%
  mutate(Throat_l = Binom_conf(x=Throat*18,n=18,lower=T)) %>%
  mutate(Throat_u = Binom_conf(x=Throat*18,n=18,lower=F)) %>%
  ggplot() + geom_point(aes(x=var+0.05, y=Nose, colour="Nose",shape="Nose"),size=2) + geom_point(aes(x=var-0.05, y=Throat, col="Throat", shape="Throat"),size=2) +
  geom_errorbar(aes(x=var+0.05,ymin=Nose_l,ymax=Nose_u, colour="Nose"),width=0.5) +
  geom_errorbar(aes(x=var-0.05,ymin=Throat_l,ymax=Throat_u, colour="Throat"),width=0.5) +
  labs(x = "Days post infection", y = "Proportion LFT positive", color = "Swab site\n",fill="",tag="A") + 
  scale_colour_manual(name = "Swab site\n (data)",
                      labels = c("Nose","Throat"),
                      values = cbPalette[c(2,3)]) +   
  scale_shape_manual(name = "Swab site\n (data)",
                     labels = c("Nose","Throat"),
                     values = c(19, 17)) +
  geom_line(data=modelLFT, aes(x=days, y=nasal,colour="Nose",linetype="Nose")) +
  geom_ribbon(data=modelLFT, aes(x=days, ymin=nasal_lower, ymax=nasal_upper,fill="Nose"),alpha=0.25,show.legend = F) +
  geom_line(data=modelLFT, aes(x=days, y=throat,colour="Throat",linetype="Throat")) +
  geom_ribbon(data=modelLFT, aes(x=days, ymin=throat_lower, ymax=throat_upper,fill="Throat"),alpha=0.25,show.legend = F) +
  scale_linetype_manual(name = "Swab site\n (model)",
                        labels = c("Nose","Throat"),
                        values = c(1,2),guide = guide_legend(override.aes = list(colour=cbPalette[c(2,3)]) )) +
  scale_fill_manual(name = "Swab site\n (model)",
                    labels = c("Nose","Throat"),
                    values = cbPalette[c(2,3)]) +
  scale_x_continuous(expand=c(0,0), limits=c(0,15)) +
  scale_y_continuous(expand=c(0,0), limits=c(-2,2)) +
  coord_cartesian(xlim=c(0,15), ylim=c(-0.002,1.002)) +
  theme_minimal(base_size=18)

p2B <- FFAtesting %>% group_by(Site) %>%
  summarise(across(c(Day_1,Day_2,Day_3,Day_4,Day_5,Day_6,Day_7,Day_8,Day_9,Day_10,Day_11,Day_12,Day_13,Day_14), ~ mean(.x, na.rm = TRUE))) %>%
  gather(var, value, -Site) %>% 
  spread(Site, value) %>%
  arrange(match(var, c("Day_1", "Day_2", "Day_3","Day_4","Day_5","Day_6","Day_7","Day_8","Day_9","Day_10","Day_11","Day_12","Day_13","Day_14"))) %>%
  mutate(var = 1:14) %>%
  mutate(Nose_l = Binom_conf(x=Nose*18,n=18,lower=T)) %>%
  mutate(Nose_u = Binom_conf(x=Nose*18,n=18,lower=F)) %>%
  mutate(Throat_l = Binom_conf(x=Throat*18,n=18,lower=T)) %>%
  mutate(Throat_u = Binom_conf(x=Throat*18,n=18,lower=F)) %>%
  ggplot() + geom_point(aes(x=var+0.05, y=Nose, colour="Nose",shape="Nose"),size=2) + geom_point(aes(x=var-0.05, y=Throat, col="Throat", shape="Throat"),size=2) +
  geom_errorbar(aes(x=var+0.05,ymin=Nose_l,ymax=Nose_u, colour="Nose"),width=0.5) +
  geom_errorbar(aes(x=var-0.05,ymin=Throat_l,ymax=Throat_u, colour="Throat"),width=0.5) +
  labs(x = "Days post infection", y = "Proportion FFA positive", color = "Swab site\n",fill="",tag="B") + 
  scale_colour_manual(name = "Swab site\n (data)",
                      labels = c("Nose","Throat"),
                      values = cbPalette[c(8,7)]) +   
  scale_shape_manual(name = "Swab site\n (data)",
                     labels = c("Nose","Throat"),
                     values = c(19, 17)) +
  geom_line(data=modelFFA, aes(x=days, y=nasal,colour="Nose",linetype="Nose")) +
  geom_ribbon(data=modelFFA, aes(x=days, ymin=nasal_lower, ymax=nasal_upper,fill="Nose"),alpha=0.25,show.legend = F) +
  geom_line(data=modelFFA, aes(x=days, y=throat,colour="Throat",linetype="Throat")) +
  geom_ribbon(data=modelFFA, aes(x=days, ymin=throat_lower, ymax=throat_upper,fill="Throat"),alpha=0.25,show.legend = F) +
  scale_linetype_manual(name = "Swab site\n (model)",
                        labels = c("Nose","Throat"),
                        values = c(1,2),guide = guide_legend(override.aes = list(colour=cbPalette[c(8,7)]) )) +
  scale_fill_manual(name = "Swab site\n (model)",
                    labels = c("Nose","Throat"),
                    values = cbPalette[c(8,7)]) +
  scale_x_continuous(expand=c(0,0), limits=c(0,15)) +
  scale_y_continuous(expand=c(0,0), limits=c(-2,2)) +
  coord_cartesian(xlim=c(0,15), ylim=c(-0.002,1.002)) +
  theme_minimal(base_size=18)

p2 <- grid.arrange(p2A,p2B,ncol=2)

#####################################
# Figure S1: Log-failure odds plots 
# To test log-logistic survival time assumption
# Time to first positive
#####################################

# # Check log-logistic assumption (straight line -> good)

fit_nasal_LFT <- survfit(Surv(time_firstpos_nasal_LFT,status_LFT_firstpos) ~ 1,conf.type="log-log") #conf.type = log vs log-log gives similar result
fit_throat_LFT <- survfit(Surv(time_firstpos_throat_LFT,status_LFT_firstpos) ~ 1,conf.type="log-log") #conf.type = log vs log-log gives similar result
fit_nasal_FFA <- survfit(Surv(time_firstpos_nasal_FFA,status_FFA_firstpos) ~ 1,conf.type="log-log") #conf.type = log vs log-log gives similar result
fit_throat_FFA <- survfit(Surv(time_firstpos_throat_FFA,status_FFA_firstpos) ~ 1,conf.type="log-log") #conf.type = log vs log-log gives similar result
fit_nasal_PCR <- survfit(Surv(time_firstpos_nasal_PCR,status_nasal_PCR_firstpos) ~ 1,conf.type="log-log") #conf.type = log vs log-log gives similar result
fit_throat_PCR <- survfit(Surv(time_firstpos_throat_PCR,status_nasal_PCR_firstpos) ~ 1,conf.type="log-log") #conf.type = log vs log-log gives similar result

S_times <- fit_nasal_LFT$time
S_hat <- fit_nasal_LFT$surv
S_data <- tibble(time=S_times,surv=S_hat)
S1A <- ggplot(data=S_data,aes(x=log(time),y=log((1-surv)/surv))) + geom_point() +
  theme_minimal(base_size=14) +
  coord_cartesian(clip = 'off') +
  labs(title="Time to first positive LFT (nasal)",x="log(time)",y="log-failure odds",tag="A") 

S_times <- fit_throat_LFT$time
S_hat <- fit_throat_LFT$surv
S_data <- tibble(time=S_times,surv=S_hat)
S1B <- ggplot(data=S_data,aes(x=log(time),y=log((1-surv)/surv))) + geom_point() +
  theme_minimal(base_size=14) +
  coord_cartesian(clip = 'off') +
  labs(title="Time to first positive LFT (throat)",x="log(time)",y="log-failure odds",tag="B") 

S_times <- fit_nasal_FFA$time
S_hat <- fit_nasal_FFA$surv
S_data <- tibble(time=S_times,surv=S_hat)
S1C <- ggplot(data=S_data,aes(x=log(time),y=log((1-surv)/surv))) + geom_point() +
  theme_minimal(base_size=14) +
  coord_cartesian(clip = 'off') +
  labs(title="Time to first positive FFA (nasal)",x="log(time)",y="log-failure odds",tag="C") 

S_times <- fit_throat_FFA$time
S_hat <- fit_throat_FFA$surv
S_data <- tibble(time=S_times,surv=S_hat)
S1D <- ggplot(data=S_data,aes(x=log(time),y=log((1-surv)/surv))) + geom_point() +
  theme_minimal(base_size=14) +
  coord_cartesian(clip = 'off') +
  labs(title="Time to first positive FFA (throat)",x="log(time)",y="log-failure odds",tag="D") 

S_times <- fit_nasal_PCR$time
S_hat <- fit_nasal_PCR$surv
S_data <- tibble(time=S_times,surv=S_hat)
S1E <- ggplot(data=S_data,aes(x=log(time),y=log((1-surv)/surv))) + geom_point() +
  theme_minimal(base_size=14) +
  coord_cartesian(clip = 'off') +
  labs(title="Time to first positive qPCR (nasal)",x="log(time)",y="log-failure odds",tag="E") 

S_times <- fit_throat_PCR$time
S_hat <- fit_throat_PCR$surv
S_data <- tibble(time=S_times,surv=S_hat)
S1F <- ggplot(data=S_data,aes(x=log(time),y=log((1-surv)/surv))) + geom_point() +
  theme_minimal(base_size=14) +
  coord_cartesian(clip = 'off') +
  labs(title="Time to first positive qPCR (throat)",x="log(time)",y="log-failure odds",tag="F") 

pS1 <- grid.arrange(S1A,S1B,S1C,S1D,S1E,S1F,ncol=2)

#####################################
# Figure S2: Log-failure odds plots 
# To test log-logistic survival time assumption
# Time to last positive (LFT and FFA only)
#####################################

# # Check log-logistic assumption (straight line -> good)

fit_nasal_LFT <- survfit(Surv(time_lastpos_nasal_LFT,status_LFT_lastpos) ~ 1,conf.type="log-log") #conf.type = log vs log-log gives similar result
fit_throat_LFT <- survfit(Surv(time_lastpos_throat_LFT,status_LFT_lastpos) ~ 1,conf.type="log-log") #conf.type = log vs log-log gives similar result
fit_nasal_FFA <- survfit(Surv(time_lastpos_nasal_FFA,status_FFA_lastpos) ~ 1,conf.type="log-log") #conf.type = log vs log-log gives similar result
fit_throat_FFA <- survfit(Surv(time_lastpos_throat_FFA,status_FFA_lastpos) ~ 1,conf.type="log-log") #conf.type = log vs log-log gives similar result

S_times <- fit_nasal_LFT$time
S_hat <- fit_nasal_LFT$surv
S_data <- tibble(time=S_times,surv=S_hat)
S2A <- ggplot(data=S_data,aes(x=log(time),y=log((1-surv)/surv))) + geom_point() +
  theme_minimal(base_size=14) +
  coord_cartesian(clip = 'off') +
  labs(title="Time to last positive LFT (nasal)",x="log(time)",y="log-failure odds",tag="A") 

S_times <- fit_throat_LFT$time
S_hat <- fit_throat_LFT$surv
S_data <- tibble(time=S_times,surv=S_hat)
S2B <- ggplot(data=S_data,aes(x=log(time),y=log((1-surv)/surv))) + geom_point() +
  theme_minimal(base_size=14) +
  coord_cartesian(clip = 'off') +
  labs(title="Time to last positive LFT (throat)",x="log(time)",y="log-failure odds",tag="B") 

S_times <- fit_nasal_FFA$time
S_hat <- fit_nasal_FFA$surv
S_data <- tibble(time=S_times,surv=S_hat)
S2C <- ggplot(data=S_data,aes(x=log(time),y=log((1-surv)/surv))) + geom_point() +
  theme_minimal(base_size=14) +
  coord_cartesian(clip = 'off') +
  labs(title="Time to last positive FFA (nasal)",x="log(time)",y="log-failure odds",tag="C") 

S_times <- fit_throat_FFA$time
S_hat <- fit_throat_FFA$surv
S_data <- tibble(time=S_times,surv=S_hat)
S2D <- ggplot(data=S_data,aes(x=log(time),y=log((1-surv)/surv))) + geom_point() +
  theme_minimal(base_size=14) +
  coord_cartesian(clip = 'off') +
  labs(title="Time to last positive FFA (throat)",x="log(time)",y="log-failure odds",tag="D") 

pS2 <- grid.arrange(S2A,S2B,S2C,S2D,ncol=2)

#####################################
# Figure S3: Time to last positive (LFT and FFA only)
# Discrete Kaplan-Meier curves and continuous log-logistic survival curves
#####################################

fit_nasal_LFT <- survfit(Surv(time_lastpos_nasal_LFT,status_LFT_lastpos) ~ 1,conf.type="log-log") #conf.type = log vs log-log gives similar result
summary(fit_nasal_LFT)
neg.plot_nasal_LFT <- ggsurvplot(fit_nasal_LFT,data=nasal_LFT)

fit_throat_LFT <- survfit(Surv(time_lastpos_throat_LFT,status_LFT_lastpos) ~ 1,conf.type="log-log") #conf.type = log vs log-log gives similar result
summary(fit_throat_LFT)
neg.plot_throat_LFT <- ggsurvplot(fit_throat_LFT,data=throat_LFT)

fit_nasal_FFA <- survfit(Surv(time_lastpos_nasal_FFA,status_FFA_lastpos) ~ 1,conf.type="log-log") #conf.type = log vs log-log gives similar result
summary(fit_nasal_FFA)
neg.plot_nasal_FFA <- ggsurvplot(fit_nasal_FFA,data=nasal_FFA)

fit_throat_FFA <- survfit(Surv(time_lastpos_throat_FFA,status_FFA_lastpos) ~ 1,conf.type="log-log") #conf.type = log vs log-log gives similar result
summary(fit_throat_FFA)
neg.plot_throat_FFA <- ggsurvplot(fit_throat_FFA,data=throat_FFA)

df_neg <- tibble(time=seq(0,18,0.1),y=S_t(seq(0,18,0.1),alphaZ_nasal_LFT,gammaZ_nasal_LFT))
pS3A <- neg.plot_nasal_LFT$plot + geom_line(data=df_neg,aes(x=time,y=y)) + 
  labs(title="Time of last positive LFT (nasal)",x="Days since infection",y="Proportion still negative",tag="A") + 
  theme_minimal(base_size=14) +
  theme(legend.position="none") +
  scale_colour_manual(values=cbPalette[c(2)]) +
  scale_fill_manual(values=cbPalette[c(2)]) +
  coord_cartesian(xlim=c(0,18), ylim=c(0,1)) +
  scale_x_continuous(breaks=c(0,5,10,15))

df_neg <- tibble(time=seq(0,18,0.1),y=S_t(seq(0,18,0.1),alphaZ_throat_LFT,gammaZ_throat_LFT))
pS3B <- neg.plot_throat_LFT$plot + geom_line(data=df_neg,aes(x=time,y=y)) + 
  labs(title="Time of last positive LFT (throat)",x="Days since infection",y="Proportion still negative",tag="B") + 
  theme_minimal(base_size=14) +
  theme(legend.position="none") +
  scale_colour_manual(values=cbPalette[c(3)]) +
  scale_fill_manual(values=cbPalette[c(3)]) +
  coord_cartesian(xlim=c(0,18), ylim=c(0,1)) +
  scale_x_continuous(breaks=c(0,5,10,15))

df_neg <- tibble(time=seq(0,18,0.1),y=S_t(seq(0,18,0.1),alphaZ_nasal_FFA,gammaZ_nasal_FFA))
pS3C <- neg.plot_nasal_FFA$plot + geom_line(data=df_neg,aes(x=time,y=y)) + 
  labs(title="Time of last positive FFA (nasal)",x="Days since infection",y="Proportion still negative",tag="C") + 
  theme_minimal(base_size=14) +
  theme(legend.position="none") +
  scale_colour_manual(values=cbPalette[c(8)]) +
  scale_fill_manual(values=cbPalette[c(8)]) +
  coord_cartesian(xlim=c(0,18), ylim=c(0,1)) +
  scale_x_continuous(breaks=c(0,5,10,15))

df_neg <- tibble(time=seq(0,18,0.1),y=S_t(seq(0,18,0.1),alphaZ_throat_FFA,gammaZ_throat_FFA))
pS3D <- neg.plot_throat_FFA$plot + geom_line(data=df_neg,aes(x=time,y=y)) + 
  labs(title="Time of last positive FFA (throat)",x="Days since infection",y="Proportion still negative",tag="D") + 
  theme_minimal(base_size=14) +
  theme(legend.position="none") +
  scale_colour_manual(values=cbPalette[c(7)]) +
  scale_fill_manual(values=cbPalette[c(7)]) +
  coord_cartesian(xlim=c(0,18), ylim=c(0,1)) +
  scale_x_continuous(breaks=c(0,5,10,15))

pS3 <- grid.arrange(pS3A,pS3B,pS3C,pS3D,ncol=2)

#####################################
# Figure S4: Time from first LFT to last FFA
# Discrete Kaplan-Meier curves and continuous log-logistic survival curves
#####################################

fit_nasal_LFTFFA <- survfit(Surv(duration_nasal_LFTFFA,status_FFA_lastpos) ~ 1,conf.type="log-log") #conf.type = log vs log-log gives similar result
summary(fit_nasal_LFTFFA)
dur.plot_nasal_LFTFFA <- ggsurvplot(fit_nasal_LFTFFA,data=nasal_LFTFFA)

df <- tibble(time=seq(0,10,0.1),y=S_t(seq(0,10,0.1),alphaY_nasal_LFTFFA,gammaY_nasal_LFTFFA))
pS4A <- dur.plot_nasal_LFTFFA$plot + geom_line(data=df,aes(x=time,y=y)) + 
  labs(title="Time from first LFT positive to last FFA positive (nasal)",x="Days since first LFT positive",y="Proportion still FFA positive",tag="A") + 
  theme_minimal(base_size=14) +
  theme(legend.position="none") +
  scale_colour_manual(values=cbPalette[c(2)]) +
  scale_fill_manual(values=cbPalette[c(2)])  + 
  coord_cartesian(xlim=c(0,10)) +
  scale_x_continuous(breaks=0:10)


fit_throat_LFTFFA <- survfit(Surv(duration_throat_LFTFFA,status_FFA_lastpos) ~ 1,conf.type="log-log") #conf.type = log vs log-log gives similar result
summary(fit_throat_LFTFFA)
dur.plot_throat_LFTFFA <- ggsurvplot(fit_throat_LFTFFA,data=throat_LFTFFA)

df <- tibble(time=seq(0,10,0.1),y=S_t(seq(0,10,0.1),alphaY_throat_LFTFFA,gammaY_throat_LFTFFA))

pS4B <- dur.plot_throat_LFTFFA$plot + geom_line(data=df,aes(x=time,y=y)) + 
  labs(title="Time from first LFT positive to last FFA positive (throat)",x="Days since first LFT positive",y="Proportion still FFA positive",tag="B") + 
  theme_minimal(base_size=14) +
  theme(legend.position="none") +
  scale_colour_manual(values=cbPalette[c(3)]) +
  scale_fill_manual(values=cbPalette[c(3)]) + 
  coord_cartesian(xlim=c(0,10)) +
  scale_x_continuous(breaks=0:10)

pS4 <- grid.arrange(pS4A,pS4B,ncol=1)

#####################################
# Figure S5: Time to last strong qPCR positive
# Original data and extended data
# Discrete Kaplan-Meier curves and continuous log-logistic survival curves
#####################################

fit_nasal_PCR <- survfit(Surv(time_lastpos_nasal_PCR,status_nasal_PCR_lastpos) ~ 1,conf.type="log-log") #conf.type = log vs log-log gives similar result
summary(fit_nasal_PCR)
neg.plot_nasal_PCR <- ggsurvplot(fit_nasal_PCR,data=nasal_PCR)

fit_throat_PCR <- survfit(Surv(time_lastpos_throat_PCR,status_throat_PCR_lastpos) ~ 1,conf.type="log-log") #conf.type = log vs log-log gives similar result
summary(fit_throat_PCR)
neg.plot_throat_PCR <- ggsurvplot(fit_throat_PCR,data=throat_PCR)

fit_nasal_PCR_ext <- survfit(Surv(time_lastpos_nasal_PCR_ext,status_nasal_PCR_lastpos_ext) ~ 1,conf.type="log-log") #conf.type = log vs log-log gives similar result
summary(fit_nasal_PCR_ext)
neg.plot_nasal_PCR_ext <- ggsurvplot(fit_nasal_PCR_ext,data=nasal_PCR_ext)

fit_throat_PCR_ext <- survfit(Surv(time_lastpos_throat_PCR_ext,status_throat_PCR_lastpos_ext) ~ 1,conf.type="log-log") #conf.type = log vs log-log gives similar result
summary(fit_throat_PCR_ext)
neg.plot_throat_PCR_ext <- ggsurvplot(fit_throat_PCR_ext,data=throat_PCR_ext)

df_neg <- tibble(time=seq(0,28,0.1),y=S_t(seq(0,28,0.1),alphaZ_nasal_PCR,gammaZ_nasal_PCR))
pS5A <- neg.plot_nasal_PCR$plot + geom_line(data=df_neg,aes(x=time,y=y)) + 
  labs(title="Time of last strong positive qPCR (nasal)",x="Days since infection",y="Proportion still negative",tag="A") + 
  theme_minimal(base_size=14) +
  theme(legend.position="none") +
  scale_colour_manual(values=cbPalette[c(8)]) +
  scale_fill_manual(values=cbPalette[c(8)]) +
  coord_cartesian(xlim=c(0,28), ylim=c(0,1)) +
  scale_x_continuous(breaks=c(0,5,10,15,20,25))

df_neg <- tibble(time=seq(0,28,0.1),y=S_t(seq(0,28,0.1),alphaZ_throat_PCR,gammaZ_throat_PCR))
pS5B <- neg.plot_throat_PCR$plot + geom_line(data=df_neg,aes(x=time,y=y)) + 
  labs(title="Time of last strong positive qPCR (throat)",x="Days since infection",y="Proportion still negative",tag="B") + 
  theme_minimal(base_size=14) +
  theme(legend.position="none") +
  scale_colour_manual(values=cbPalette[c(7)]) +
  scale_fill_manual(values=cbPalette[c(7)]) +
  coord_cartesian(xlim=c(0,28), ylim=c(0,1)) +
  scale_x_continuous(breaks=c(0,5,10,15,20,25))

df_neg <- tibble(time=seq(0,28,0.1),y=S_t(seq(0,28,0.1),alphaZ_nasal_PCR_ext,gammaZ_nasal_PCR_ext))
pS5C <- neg.plot_nasal_PCR_ext$plot + geom_line(data=df_neg,aes(x=time,y=y)) + 
  labs(title="Time of last strong positive qPCR (nasal)\n extended data",x="Days since infection",y="Proportion still negative",tag="C") + 
  theme_minimal(base_size=14) +
  theme(legend.position="none") +
  scale_colour_manual(values=cbPalette[c(8)]) +
  scale_fill_manual(values=cbPalette[c(8)]) +
  coord_cartesian(xlim=c(0,28), ylim=c(0,1)) +
  scale_x_continuous(breaks=c(0,5,10,15,20,25))

df_neg <- tibble(time=seq(0,28,0.1),y=S_t(seq(0,28,0.1),alphaZ_throat_PCR_ext,gammaZ_throat_PCR_ext))
pS5D <- neg.plot_throat_PCR_ext$plot + geom_line(data=df_neg,aes(x=time,y=y)) + 
  labs(title="Time of last strong positive qPCR (throat)\n extended data",x="Days since infection",y="Proportion still negative",tag="D") + 
  theme_minimal(base_size=14) +
  theme(legend.position="none") +
  scale_colour_manual(values=cbPalette[c(7)]) +
  scale_fill_manual(values=cbPalette[c(7)]) +
  coord_cartesian(xlim=c(0,28), ylim=c(0,1)) +
  scale_x_continuous(breaks=c(0,5,10,15,20,25))

pS5 <- grid.arrange(pS5A,pS5B,pS5C,pS5D,ncol=2)

#####################################
# Figure S6: Time to last (any) qPCR positive
# Extended data only
# Discrete Kaplan-Meier curves and continuous log-logistic survival curves
#####################################

fit_nasal_PCR_weak_ext <- survfit(Surv(time_lastpos_nasal_PCR_weak_ext,status_nasal_PCR_weak_lastpos_ext) ~ 1,conf.type="log-log") #conf.type = log vs log-log gives similar result
summary(fit_nasal_PCR_weak_ext)
neg.plot_nasal_PCR_weak_ext <- ggsurvplot(fit_nasal_PCR_weak_ext,data=nasal_PCR_weak_ext)

fit_throat_PCR_weak_ext <- survfit(Surv(time_lastpos_throat_PCR_weak_ext,status_throat_PCR_weak_lastpos_ext) ~ 1,conf.type="log-log") #conf.type = log vs log-log gives similar result
summary(fit_throat_PCR_weak_ext)
neg.plot_throat_PCR_weak_ext <- ggsurvplot(fit_throat_PCR_weak_ext,data=throat_PCR_weak_ext)

df_neg <- tibble(time=seq(0,28,0.1),y=S_t(seq(0,28,0.1),alphaZ_nasal_PCR_weak_ext,gammaZ_nasal_PCR_weak_ext))
pS6A <- neg.plot_nasal_PCR_weak_ext$plot + geom_line(data=df_neg,aes(x=time,y=y)) + 
  labs(title="Time of last positive qPCR (nasal)\n extended data",x="Days since infection",y="Proportion still negative",tag="A") + 
  theme_minimal(base_size=14) +
  theme(legend.position="none") +
  scale_colour_manual(values=cbPalette[c(8)]) +
  scale_fill_manual(values=cbPalette[c(8)]) +
  coord_cartesian(xlim=c(0,28), ylim=c(0,1)) +
  scale_x_continuous(breaks=c(0,5,10,15,20,25))

df_neg <- tibble(time=seq(0,28,0.1),y=S_t(seq(0,28,0.1),alphaZ_throat_PCR_weak_ext,gammaZ_throat_PCR_weak_ext))

pS6B <- neg.plot_throat_PCR_weak_ext$plot + geom_line(data=df_neg,aes(x=time,y=y)) + 
  labs(title="Time of last positive qPCR (throat)\n extended data",x="Days since infection",y="Proportion still negative",tag="B") + 
  theme_minimal(base_size=14) +
  theme(legend.position="none") +
  scale_colour_manual(values=cbPalette[c(7)]) +
  scale_fill_manual(values=cbPalette[c(7)]) +
  coord_cartesian(xlim=c(0,28), ylim=c(0,1)) +
  scale_x_continuous(breaks=c(0,5,10,15,20,25))

pS6 <- grid.arrange(pS6A,pS6B,ncol=2)

#####################################
# Figure S7: Temporal sensitivity model (qPCR)
# Error windows: 95% percentile case-resampling bootstrap
# Data + 95% binom CI errorbars
# Strong, Strong (extended data) and Any (extended data)
#####################################

PCRtesting <- PCRtesting %>% mutate(Day_0 = 0)
PCRtesting_ext <- PCRtesting_strong_ext %>% mutate(Day_0 = 0)
PCRtesting_weak_ext <- PCRtesting_weak_ext %>% mutate(Day_0 = 0)

modelPCR <- data.frame(days=0:20,nasal=probsB_nasal_PCR,throat=probsB_throat_PCR,
                       nasal_upper=upper_nasal_PCR,nasal_lower=lower_nasal_PCR,
                       throat_upper=upper_throat_PCR,throat_lower=lower_throat_PCR)
modelPCR_ext <- data.frame(days=0:20,nasal=probsB_nasal_PCR_ext,throat=probsB_throat_PCR_ext,
                           nasal_upper=upper_nasal_PCR_ext,nasal_lower=lower_nasal_PCR_ext,
                           throat_upper=upper_throat_PCR_ext,throat_lower=lower_throat_PCR_ext)
modelPCR_weak_ext <- data.frame(days=0:20,nasal=probsB_nasal_PCR_weak_ext,throat=probsB_throat_PCR_weak_ext,
                                nasal_upper=upper_nasal_PCR_weak_ext,nasal_lower=lower_nasal_PCR_weak_ext,
                                throat_upper=upper_throat_PCR_weak_ext,throat_lower=lower_throat_PCR_weak_ext)

pS7A <- PCRtesting_strong %>% group_by(Site) %>%
  summarise(across(c(Day_1,Day_2,Day_3,Day_4,Day_5,Day_6,Day_7,Day_8,Day_9,Day_10,Day_11,Day_12,Day_13,Day_14), ~ mean(.x, na.rm = TRUE))) %>%
  gather(var, value, -Site) %>% 
  spread(Site, value) %>%
  arrange(match(var, c("Day_1", "Day_2", "Day_3","Day_4","Day_5","Day_6","Day_7","Day_8","Day_9","Day_10","Day_11","Day_12","Day_13","Day_14"))) %>%
  mutate(var = 1:14) %>%
  mutate(Nose_l = Binom_conf(x=Nose*18,n=18,lower=T)) %>%
  mutate(Nose_u = Binom_conf(x=Nose*18,n=18,lower=F)) %>%
  mutate(Throat_l = Binom_conf(x=Throat*18,n=18,lower=T)) %>%
  mutate(Throat_u = Binom_conf(x=Throat*18,n=18,lower=F)) %>%
  ggplot() + geom_point(aes(x=var+0.05, y=Nose, colour="Nose",shape="Nose"),size=2) + geom_point(aes(x=var-0.05, y=Throat, col="Throat", shape="Throat"),size=2) +
  geom_errorbar(aes(x=var+0.05,ymin=Nose_l,ymax=Nose_u, colour="Nose"),width=0.5) +
  geom_errorbar(aes(x=var-0.05,ymin=Throat_l,ymax=Throat_u, colour="Throat"),width=0.5) +
  labs(x = "Days post infection", y = "Proportion strong positive",title="qPCR positivity over time since exposure\n (strong positive only)", color = "Swab site\n",fill="",tag="A") + 
  scale_colour_manual(name = "Swab site\n (data)",
                      labels = c("Nose","Throat"),
                      values = cbPalette[c(6,4)]) +   
  scale_shape_manual(name = "Swab site\n (data)",
                     labels = c("Nose","Throat"),
                     values = c(19, 17)) +
  geom_line(data=modelPCR, aes(x=days, y=nasal,colour="Nose",linetype="Nose")) +
  geom_ribbon(data=modelPCR, aes(x=days, ymin=nasal_lower, ymax=nasal_upper,fill="Nose"),alpha=0.25,show.legend = F) +
  geom_line(data=modelPCR, aes(x=days, y=throat,colour="Throat",linetype="Throat")) +
  geom_ribbon(data=modelPCR, aes(x=days, ymin=throat_lower, ymax=throat_upper,fill="Throat"),alpha=0.25,show.legend = F) +
  scale_linetype_manual(name = "Swab site\n (model)",
                        labels = c("Nose","Throat"),
                        values = c(1,2),guide = guide_legend(override.aes = list(colour=cbPalette[c(6,4)]) )) +
  scale_fill_manual(name = "Swab site\n (model)",
                    labels = c("Nose","Throat"),
                    values = cbPalette[c(6,4)]) +
  scale_x_continuous(expand=c(0,0), limits=c(0,20)) +
  scale_y_continuous(expand=c(0,0), limits=c(-2,2)) +
  coord_cartesian(xlim=c(0,20), ylim=c(-0.002,1.002)) +
  theme_minimal(base_size=14)


pS7B <- PCRtesting_ext %>% group_by(Site) %>%
  summarise(across(c(Day_1,Day_2,Day_3,Day_4,Day_5,Day_6,Day_7,Day_8,Day_9,Day_10,
                     Day_11,Day_12,Day_13,Day_14,Day_15,Day_16,Day_17,Day_18,Day_19), ~ mean(.x, na.rm = TRUE))) %>%
  gather(var, value, -Site) %>% 
  spread(Site, value) %>%
  arrange(match(var, c("Day_1", "Day_2", "Day_3","Day_4","Day_5","Day_6","Day_7",
                       "Day_8","Day_9","Day_10","Day_11","Day_12","Day_13","Day_14",
                       "Day_15","Day_16","Day_17","Day_18","Day_19"))) %>%
  mutate(var = 1:19) %>%
  mutate(Nose_l = Binom_conf(x=Nose*c(rep(18,15),13,8,5,2),n=c(rep(18,15),13,8,5,2),lower=T)) %>%
  mutate(Nose_u = Binom_conf(x=Nose*c(rep(18,15),13,8,5,2),n=c(rep(18,15),13,8,5,2),lower=F)) %>%
  mutate(Throat_l = Binom_conf(x=Throat*c(rep(18,15),13,8,5,2),n=c(rep(18,15),13,8,5,2),lower=T)) %>%
  mutate(Throat_u = Binom_conf(x=Throat*c(rep(18,15),13,8,5,2),n=c(rep(18,15),13,8,5,2),lower=F)) %>%
  ggplot() + geom_point(aes(x=var+0.05, y=Nose, colour="Nose",shape="Nose"),size=2) + geom_point(aes(x=var-0.05, y=Throat, col="Throat", shape="Throat"),size=2) +
  geom_errorbar(aes(x=var+0.05,ymin=Nose_l,ymax=Nose_u, colour="Nose"),width=0.5) +
  geom_errorbar(aes(x=var-0.05,ymin=Throat_l,ymax=Throat_u, colour="Throat"),width=0.5) +
  labs(x = "Days post infection", y = "Proportion strong positive",title="qPCR positivity over time since exposure\n extended data (strong positive only)", color = "Swab site\n",fill="",tag="B") + 
  scale_colour_manual(name = "Swab site\n (data)",
                      labels = c("Nose","Throat"),
                      values = cbPalette[c(6,4)]) +   
  scale_shape_manual(name = "Swab site\n (data)",
                     labels = c("Nose","Throat"),
                     values = c(19, 17)) +
  geom_line(data=modelPCR_ext, aes(x=days, y=nasal,colour="Nose",linetype="Nose")) +
  geom_ribbon(data=modelPCR_ext, aes(x=days, ymin=nasal_lower, ymax=nasal_upper,fill="Nose"),alpha=0.25,show.legend = F) +
  geom_line(data=modelPCR_ext, aes(x=days, y=throat,colour="Throat",linetype="Throat")) +
  geom_ribbon(data=modelPCR_ext, aes(x=days, ymin=throat_lower, ymax=throat_upper,fill="Throat"),alpha=0.25,show.legend = F) +
  scale_linetype_manual(name = "Swab site\n (model)",
                        labels = c("Nose","Throat"),
                        values = c(1,2),guide = guide_legend(override.aes = list(colour=cbPalette[c(6,4)]) )) +
  scale_fill_manual(name = "Swab site\n (model)",
                    labels = c("Nose","Throat"),
                    values = cbPalette[c(6,4)]) +
  scale_x_continuous(expand=c(0,0), limits=c(0,20)) +
  scale_y_continuous(expand=c(0,0), limits=c(-2,2)) +
  coord_cartesian(xlim=c(0,20), ylim=c(-0.002,1.002)) +
  theme_minimal(base_size=14)


pS7C <- PCRtesting_weak_ext %>% group_by(Site) %>%
  summarise(across(c(Day_1,Day_2,Day_3,Day_4,Day_5,Day_6,Day_7,Day_8,Day_9,Day_10,
                     Day_11,Day_12,Day_13,Day_14,Day_15,Day_16,Day_17,Day_18,Day_19), ~ mean(.x, na.rm = TRUE))) %>%
  gather(var, value, -Site) %>% 
  spread(Site, value) %>%
  arrange(match(var, c("Day_1", "Day_2", "Day_3","Day_4","Day_5","Day_6","Day_7",
                       "Day_8","Day_9","Day_10","Day_11","Day_12","Day_13","Day_14",
                       "Day_15","Day_16","Day_17","Day_18","Day_19"))) %>%
  mutate(var = 1:19) %>%
  mutate(Nose_l = Binom_conf(x=Nose*c(rep(18,15),13,8,5,2),n=c(rep(18,15),13,8,5,2),lower=T)) %>%
  mutate(Nose_u = Binom_conf(x=Nose*c(rep(18,15),13,8,5,2),n=c(rep(18,15),13,8,5,2),lower=F)) %>%
  mutate(Throat_l = Binom_conf(x=Throat*c(rep(18,15),13,8,5,2),n=c(rep(18,15),13,8,5,2),lower=T)) %>%
  mutate(Throat_u = Binom_conf(x=Throat*c(rep(18,15),13,8,5,2),n=c(rep(18,15),13,8,5,2),lower=F)) %>%
  ggplot() + geom_point(aes(x=var+0.05, y=Nose, colour="Nose",shape="Nose"),size=2) + geom_point(aes(x=var-0.05, y=Throat, col="Throat", shape="Throat"),size=2) +
  geom_errorbar(aes(x=var+0.05,ymin=Nose_l,ymax=Nose_u, colour="Nose"),width=0.5) +
  geom_errorbar(aes(x=var-0.05,ymin=Throat_l,ymax=Throat_u, colour="Throat"),width=0.5) +
  labs(x = "Days post infection", y = "Proportion positive",title="qPCR positivity over time since exposure\n extended data", color = "Swab site\n",fill="",tag="C") + 
  scale_colour_manual(name = "Swab site\n (data)",
                      labels = c("Nose","Throat"),
                      values = cbPalette[c(6,4)]) +   
  scale_shape_manual(name = "Swab site\n (data)",
                     labels = c("Nose","Throat"),
                     values = c(19, 17)) +
  geom_line(data=modelPCR_weak_ext, aes(x=days, y=nasal,colour="Nose",linetype="Nose")) +
  geom_ribbon(data=modelPCR_weak_ext, aes(x=days, ymin=nasal_lower, ymax=nasal_upper,fill="Nose"),alpha=0.25,show.legend = F) +
  geom_line(data=modelPCR_weak_ext, aes(x=days, y=throat,colour="Throat",linetype="Throat")) +
  geom_ribbon(data=modelPCR_weak_ext, aes(x=days, ymin=throat_lower, ymax=throat_upper,fill="Throat"),alpha=0.25,show.legend = F) +
  scale_linetype_manual(name = "Swab site\n (model)",
                        labels = c("Nose","Throat"),
                        values = c(1,2),guide = guide_legend(override.aes = list(colour=cbPalette[c(6,4)]) )) +
  scale_fill_manual(name = "Swab site\n (model)",
                    labels = c("Nose","Throat"),
                    values = cbPalette[c(6,4)]) +
  scale_x_continuous(expand=c(0,0), limits=c(0,20)) +
  scale_y_continuous(expand=c(0,0), limits=c(-2,2)) +
  coord_cartesian(xlim=c(0,20), ylim=c(-0.002,1.002)) +
  theme_minimal(base_size=14)

pS7 <- grid.arrange(pS7A,pS7B,pS7C,ncol=1)

#####################################
# Figure S8: Temporal sensitivity model (nasal LFT)
# Error windows: 95% percentile case-resampling bootstrap
# Scaled for asymptomatic and swab administration route (self-swab, NHS, lab expert)
#####################################

# ratios from Dinnes et al (https://doi.org/10.1002/14651858.CD013705.pub2)
asym.ratio <- 58.1/72
selfswab <- 57.5
nhsswab <- 70
labswab <- 78.8
swab.ratio <- selfswab/labswab
nhs.ratio <- nhsswab/labswab

Day.no <- c(0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14)

pS8A <- LFTtesting %>% group_by(Site) %>%
  summarise(across(c(Day_0,Day_1,Day_2,Day_3,Day_4,Day_5,Day_6,Day_7,Day_8,Day_9,Day_10,Day_11,Day_12,Day_13,Day_14), ~ mean(.x, na.rm = TRUE))) %>%
  gather(var, value, -Site) %>% 
  spread(Site, value) %>%
  arrange(match(var, c("Day_0","Day_1", "Day_2", "Day_3","Day_4","Day_5","Day_6","Day_7","Day_8","Day_9","Day_10","Day_11","Day_12","Day_13","Day_14"))) %>%
  mutate(var = 0:14) %>%
  ggplot() + labs(title = "LFT positivity by symptoms\nnasal swab", x = "Days post infection", y = "Proportion LFT positive", color = "Swab site\n",fill="",tag="A") + 
  scale_colour_manual(name = "Symptoms",
                      labels = c("No","Yes"),
                      values = cbPalette[c(3,2)]) +   
  geom_line(data=modelLFT, aes(x=days, y=nasal*asym.ratio,colour="Asymptomatic",linetype="Asymptomatic")) +
  geom_ribbon(data=modelLFT, aes(x=days, ymin=nasal_lower*asym.ratio, ymax=nasal_upper*asym.ratio,fill="Asymptomatic"),alpha=0.25,show.legend = F) +
  geom_line(data=modelLFT, aes(x=days, y=nasal,colour="Symptomatic",linetype="Symptomatic")) +
  geom_ribbon(data=modelLFT, aes(x=days, ymin=nasal_lower, ymax=nasal_upper,fill="Nose"),alpha=0.25,show.legend = F) +
  scale_linetype_manual(name = "Symptoms",
                        labels = c("No","Yes"),
                        values = c(2,1),guide = guide_legend(override.aes = list(colour=cbPalette[c(3,2)]) )) +
  scale_fill_manual(name = "Symptoms",
                    labels = c("No","Yes"),
                    values = cbPalette[c(3,2)]) +
  scale_x_continuous(expand=c(0,0), limits=c(0,15)) +
  scale_y_continuous(expand=c(0,0), limits=c(-2,2)) +
  coord_cartesian(xlim=c(0,15), ylim=c(-0.002,1.002)) +
  theme_minimal(base_size=14)

pS8B <- LFTtesting %>% group_by(Site) %>%
  summarise(across(c(Day_0,Day_1,Day_2,Day_3,Day_4,Day_5,Day_6,Day_7,Day_8,Day_9,Day_10,Day_11,Day_12,Day_13,Day_14), ~ mean(.x, na.rm = TRUE))) %>%
  gather(var, value, -Site) %>% 
  spread(Site, value) %>%
  arrange(match(var, c("Day_0","Day_1", "Day_2", "Day_3","Day_4","Day_5","Day_6","Day_7","Day_8","Day_9","Day_10","Day_11","Day_12","Day_13","Day_14"))) %>%
  mutate(var = 0:14) %>%
  ggplot() + 
  labs(title="LFT positivity by swab type\nnasal swab",x = "Days post infection", y = "Proportion LFT positive", color = "Swab site\n",fill="",tag="B") + 
  geom_line(data=modelLFT, aes(x=days, y=nasal*swab.ratio,colour="Self-swab",linetype="Self-swab")) +
  geom_ribbon(data=modelLFT, aes(x=days, ymin=nasal_lower*swab.ratio, ymax=nasal_upper*swab.ratio,fill="Self-swab"),alpha=0.25,show.legend = F) +
  geom_line(data=modelLFT, aes(x=days, y=nasal,colour="1Lab expert",linetype="1Lab expert")) +
  geom_ribbon(data=modelLFT, aes(x=days, ymin=nasal_lower, ymax=nasal_upper,fill="1Lab expert"),alpha=0.25,show.legend = F) +
  geom_line(data=modelLFT, aes(x=days, y=nasal*nhs.ratio,colour="Health worker",linetype="Health worker")) +
  geom_ribbon(data=modelLFT, aes(x=days, ymin=nasal_lower*nhs.ratio, ymax=nasal_upper*nhs.ratio,fill="Health worker"),alpha=0.25,show.legend = F) +
  scale_linetype_manual(name = "Swab type",
                        labels = c("Lab expert","Health worker","Self-swab"),
                        values = c(1,2,3),guide = guide_legend(override.aes = list(colour=cbPalette[c(2,3,4)]))) +
  scale_fill_manual(name = "Swab type",
                    labels = c("Lab expert","Health worker","Self-swab"),
                    values = cbPalette[c(2,3,4)]) +
  scale_color_manual(name = "Swab type",
                     labels = c("Lab expert","Health worker","Self-swab"),
                     values = cbPalette[c(2,3,4)]) +
  scale_x_continuous(expand=c(0,0), limits=c(0,15)) +
  scale_y_continuous(expand=c(0,0), limits=c(-2,2)) +
  coord_cartesian(xlim=c(0,15), ylim=c(-0.002,1.002)) +
  theme_minimal(base_size=14)

pS8 <- grid.arrange(pS8A,pS8B,ncol=2)