##############################################################################
############################### Pre-processing ###############################
##############################################################################

library(medicaldata)
df <- medicaldata::opt
setwd("C:/Users/matte/Desktop/NonParam_OPT_Project")
source("DATA PREPROCESSING.R")

# NB: low birthweight if birthweight < 2500 g
#     pre-term birth if gestional age < 259 gg

##############################################################################
############################## SURVIVAL ANALYSIS #############################
##############################################################################

library(survival)
library(survminer)
library(dplyr) 
library(ggplot2)
library(knitr)
library(broom)

n <- dim(df)[1]
times <- df[,176]  # Delivery times; unit of measure = days
lost_to_FU <- which(is.na(df[,70]))
censoring <- rep(0,n)
censoring[lost_to_FU] <- 1
censoring <- factor(censoring, labels = (c('Event', 'Censor')))
groups <- df$Group  # or df$Diabetes to see the effect of other risk factor

threshold <- 240
times_threshold <- times[which(times <= threshold)]
censoring_threshold <- censoring[which(times <= threshold)]
groups_threshold <- df[which(times <= threshold),3]


##############################################################################
########################## KAPLAN-MEIER ESTIMATOR ############################
##############################################################################

KM <- survfit(Surv(times, censoring=="Event") ~ 1)
summary(KM)
kable(head(tidy(KM),20))
surv_median(KM)

plot(KM, conf.int = T, xlab='Time [days]', ylab = 'Survival Probability', col='red',
     main="Kaplan-Meier Curve for Lung Cancer Survival", xlim = range(times))

ggsurvplot(KM, data=as.data.frame(cbind(times,censoring)), risk.table = TRUE, risk.table.col = "strata", surv.median.line = "hv",
           ggtheme = theme_bw(), break.time.by=90,  title="Kaplan-Meier Curve Delivery Times",xlim = range(times))

ggsurvplot(KM, data=as.data.frame(cbind(times,censoring)), risk.table = TRUE, ggtheme = theme_bw(),
           break.time.by=90, fun='cumhaz', title="Cumulative Hazard Curve for Delivery Times",xlim = range(times))

# Conider the pre-term births only
KM_threshold <- survfit(Surv(times_threshold, censoring_threshold=="Event") ~ 1)
summary(KM_threshold)
kable(head(tidy(KM_threshold),20))
surv_median(KM_threshold)

plot(KM_threshold, conf.int = T, xlab='Time [days]', ylab = 'Survival Probability', col='red',
     main="Kaplan-Meier Curve for Lung Cancer Survival", xlim = range(times_threshold))

ggsurvplot(KM_threshold, data=as.data.frame(cbind(times_threshold,censoring_threshold)), risk.table = TRUE, risk.table.col = "strata", surv.median.line = "hv",
           ggtheme = theme_bw(), break.time.by=90,  title="Kaplan-Meier Curve Delivery Times",xlim = range(times_threshold))

ggsurvplot(KM_threshold, data=as.data.frame(cbind(times_threshold,censoring_threshold)), risk.table = TRUE, ggtheme = theme_bw(),
           break.time.by=90, fun='cumhaz', title="Cumulative Hazard Curve for Delivery Times",xlim = range(times_threshold))


##############################################################################
################### KAPLAN-MEIER ESTIMATOR FOR DIFFERENT GROUPS ##############
##############################################################################

KM_groups <- survfit(Surv(times, censoring=="Event") ~ groups)

ggsurvplot(KM_groups, data=as.data.frame(cbind(times,censoring)), conf.int = F, risk.table = FALSE, risk.table.col = "strata",xlim = range(times),  
           surv.median.line = "hv", ggtheme = theme_bw(), break.time.by=90,
           title="Kaplan-Meier Curves for Delivery times")

log_rank_test <- survdiff(Surv(times, censoring=="Event") ~ groups)
log_rank_test

hazard_ratio <- (log_rank_test$obs[1]/log_rank_test$exp[1])/(log_rank_test$obs[2]/log_rank_test$exp[2])
hazard_ratio

# Conider the pre-term births only
KM_groups_threshold <- survfit(Surv(times_threshold, censoring_threshold=="Event") ~ groups_threshold)

ggsurvplot(KM_groups_threshold, data=as.data.frame(cbind(times_threshold,censoring_threshold)), conf.int = F, risk.table = FALSE, risk.table.col = "strata",xlim = range(times),  
           surv.median.line = "hv", ggtheme = theme_bw(), break.time.by=90,
           title="Kaplan-Meier Curves for Delivery times")

log_rank_test_threshold <- survdiff(Surv(times_threshold, censoring_threshold=="Event") ~ groups_threshold)
log_rank_test_threshold

hazard_ratio_threshold <- (log_rank_test_threshold$obs[1]/log_rank_test_threshold$exp[1])/(log_rank_test_threshold$obs[2]/log_rank_test_threshold$exp[2])
hazard_ratio_threshold


##############################################################################
################################# COX MODEL ##################################
##############################################################################

cov1 <- df$Diabetes  # select covariates significant for the delivery time
cov2 <- df$Hypertension
cox <- coxph(Surv(times, censoring=="Event") ~ cov1 + cov2)
summary(cox)
x11()
ggforest(cox, data = as.data.frame(cbind(cov1,cov2,times,censoring)))

plot(survfit(cox, data=as.data.frame(cbind(cov1,cov2,times,censoring))), 
     col="darkorange2", lwd=2, lty=1,
     xlab='Time [days]', ylab='Survival Probability',
     main='Baseline estimated survival probability')
grid()

# Consider the pre-term births only
cov1_threshold <- df$Diabetes[which(times <= threshold)]  # select covariates significant for the delivery time (see permutational inference results)
cov2_threshold <- df$Hypertension[which(times <= threshold)]
cox_threshold <- coxph(Surv(times_threshold, censoring_threshold=="Event") ~ cov1_threshold + cov2_threshold)
summary(cox_threshold)
x11()
ggforest(cox_threshold, data = as.data.frame(cbind(cov1_threshold,cov2_threshold,times_threshold,censoring_threshold)))

plot(survfit(cox_threshold, data=as.data.frame(cbind(cov1_threshold,cov2_threshold,times_threshold,censoring_threshold))), 
     col="darkorange2", lwd=2, lty=1,
     xlab='Time [days]', ylab='Survival Probability',
     main='Baseline estimated survival probability')
grid()

# Goodness-of-fit
ggcoxdiagnostics(cox, type = "martingale")
ggcoxdiagnostics(cox, type = "deviance")
ggcoxdiagnostics(cox, type = "schoenfeld")
test.ph <- cox.zph(cox)
test.ph

ggcoxdiagnostics(cox_threshold, type = "martingale")
ggcoxdiagnostics(cox_threshold, type = "deviance")
ggcoxdiagnostics(cox_threshold, type = "schoenfeld")
test.ph_threshold <- cox.zph(cox_threshold)
test.ph_threshold