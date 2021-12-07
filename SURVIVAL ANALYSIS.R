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
groups <- df$Group # select categorical variable significant for delivery times (see permutational inference results)

# Apply thresholding to see only pre-term birth
threshold <- 230
times_threshold <- times[which(times <= threshold)]
censoring_threshold <- censoring[which(times <= threshold)]
groups_threshold <- groups[which(times <= threshold)]


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

### Conider the pre-term births only
KM_threshold <- survfit(Surv(times_threshold, censoring_threshold=="Event") ~ 1)
summary(KM_threshold)
kable(head(tidy(KM_threshold),20))
surv_median(KM_threshold)

plot(KM_threshold, conf.int = T, xlab='Time [days]', ylab = 'Survival Probability', col='red',
     main="Kaplan-Meier Curve for Delivery times", xlim = range(times_threshold))

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

hazard_ratio <- (log_rank_test$obs[2]/log_rank_test$exp[2])/(log_rank_test$obs[1]/log_rank_test$exp[1])
hazard_ratio
# NB: we have variables with levels NO/YES intead of YES/NO => we change the computation of the Hazard Ratio wrt what we have done in class

# Consider the pre-term births only
KM_groups_threshold <- survfit(Surv(times_threshold, censoring_threshold=="Event") ~ groups_threshold)

ggsurvplot(KM_groups_threshold, data=as.data.frame(cbind(times_threshold,censoring_threshold)), conf.int = F, risk.table = FALSE, risk.table.col = "strata",xlim = range(times),  
           surv.median.line = "hv", ggtheme = theme_bw(), break.time.by=90,
           title="Kaplan-Meier Curves for Delivery times")

log_rank_test_threshold <- survdiff(Surv(times_threshold, censoring_threshold=="Event") ~ groups_threshold)
log_rank_test_threshold

hazard_ratio_threshold <- (log_rank_test_threshold$obs[2]/log_rank_test_threshold$exp[2])/(log_rank_test_threshold$obs[1]/log_rank_test_threshold$exp[1])
hazard_ratio_threshold


##############################################################################
################################# COX MODEL ##################################
##############################################################################

cov1 <- df$BMI # select covariates significant for the delivery time
cox <- coxph(Surv(times, censoring=="Event") ~ cov1)
summary(cox)
x11()
ggforest(cox, data = as.data.frame(cbind(cov1,times,censoring)))

# Baseline survival curve
x11()
plot(survfit(cox, data=as.data.frame(cbind(cov1,cov2,times,censoring))), 
     col="darkorange2", lwd=2, lty=1,
     xlab='Time [days]', ylab='Survival Probability',
     main='Baseline estimated survival probability')
grid()

# Plot the survival curve for different values of a covariate
cov_values <- data.frame(cov1 = c(min(cov1[which(!is.na(cov1))]),mean(cov1[which(!is.na(cov1))]),max(cov1[which(!is.na(cov1))])))
KM_cov <- survfit(cox, newdata = cov_values)
x11()
plot(KM_cov, conf.int=F,
     col=c("dodgerblue2","navy","darkmagenta"), lwd=2, lty=1,
     xlab='Time [days]', ylab='Survival Probability',
     main='Adjusted Survival Probability Plot')
grid()
legend('bottomleft', c("min", "mean", "max"),
       lty=c(1,1,1), lwd=c(2,2,2), col=c("dodgerblue2","navy","darkmagenta"))

### Consider the pre-term births only
cov1_threshold <- df$BMI[which(times <= threshold)]  # select covariates significant for the delivery time (see permutational inference results)
cox_threshold <- coxph(Surv(times_threshold, censoring_threshold=="Event") ~ cov1_threshold)
summary(cox_threshold)
x11()
ggforest(cox_threshold, data = as.data.frame(cbind(cov1_threshold,times_threshold,censoring_threshold)))

# Plot baseline survival curve
plot(survfit(cox_threshold, data=as.data.frame(cbind(cov1_threshold,cov2_threshold,times_threshold,censoring_threshold))), 
     col="darkorange2", lwd=2, lty=1,
     xlab='Time [days]', ylab='Survival Probability',
     main='Baseline estimated survival probability')
grid()

# Plot the survival curve for different values of a covariate
cov_values_threshold <- data.frame(cov1_threshold = c(min(cov1_threshold[which(!is.na(cov1_threshold))]),mean(cov1_threshold[which(!is.na(cov1_threshold))]),max(cov1_threshold[which(!is.na(cov1_threshold))])))
KM_cov_threshold <- survfit(cox_threshold, newdata = cov_values_threshold)
x11()
plot(KM_cov_threshold, conf.int=F,
     col=c("dodgerblue2","navy","darkmagenta"), lwd=2, lty=1,
     xlab='Time [days]', ylab='Survival Probability',
     main='Adjusted Survival Probability Plot')
grid()
legend('bottomleft', c("min", "mean", "max"),
       lty=c(1,1,1), lwd=c(2,2,2), col=c("dodgerblue2","navy","darkmagenta"))


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


##############################################################################
############################### CONCLUSIONS ##################################
##############################################################################

# GROUP               THRESHOLD    PVALUE LOG-RANK TEST      HAZARD RATIO
# Treat/Control       <=230        0.03                      0.5 
# Hips Y/N            //           0.02                      0.83
# Pubblic Asstce Y/N  //           0.04                      1.16
# Hypertension Y/N    //           6e-06                     2.36
# Diabetes Y/N        //           2e-04                     2.06
# Use Tob Y/N         //           0.008                     1.32

# a) Looking at the pre-term births only, we have a significance difference for the survivals curve
#    in Control and Treatment groups. Moreover being in the Control group is a risk factor.
# b) Risk factors as hypertension, diabetes, tobacco, hisp give births previously. 
# c) Completed EDC, EDC necessary, previous pregnancies, alcohol, drug, Blank race do not give different survival curves

# COVARIATE       THRESHOLD      PVALUE COX MODEL        HAZARD RATIO
# AGE             //             0.934                   1.001
# BMI             //             0.0362                  1.011
# BL.GE           <=230          0.0465                  0.4356
# BL.PD.avg       <=230          0.0388                  0.4965
# V3.Calc.I       <=230          0.0594                  4.293 
# V3.Calc.I       <=235          0.0796                  3.22
# V3.Calc.I       <=260          0.101                   1.346

# a) BMI is a risk factor; Age is not a risk factor

# b) Gingival indexes and Bacteria5.perc in the COX model are protective factor wrt pre-term birth. 
#    This is strange and counter-intuitive, a possible explaination is given by the NA's and the point (ii) below.

# c) If we apply a threshold (i.e. looking at severe pre-term births only) we notice that the gingival factors
# at the baseline (BL) are still protective but they become risky when measured at visit 3 (V3).
# This is coherent with the treatment effect: 
# (i)  treated patients have lower indeces in V3 (reduction in the oral infection) and these indeces in V3 are risk factors 
#      ==> treatment helps to have less pre-term births (we have these risk factors smaller after treatment)
# (ii) Unfortunately we have a lot of NA's in V3 and V5, it's difficult to assess that they are risk factors in general.
#      We notice that ginigval indeces at the baseline BL are protective factors, indeed:
#      
#      higher indeces at BL ==> treatment ==> improvement of oral health ==> less pre-term births since it is protective.
#
#      Hence having a bad situation at BL (high indeces) results in a protective factor thanks to the treatment.


# Data analysis for gengival indeces
control.BL.GE <- df$BL.GE[which(df$Group == 'C')]
treat.BL.GE <- df$BL.GE[which(df$Group == 'T')]
control.BL.Calc.I <- df$BL.Calc.I[which(df$Group == 'C')]
treat.BL.Calc.I <- df$BL.Calc.I[which(df$Group == 'T')]
control.V3.GE <- df$V3.GE[which(df$Group == 'C')]
treat.V3.GE <- df$V3.GE[which(df$Group == 'T')]
control.V5.GE <- df$V5.GE[which(df$Group == 'C')]
treat.V5.GE <- df$V5.GE[which(df$Group == 'T')]
control.V5.PD.avg <- df$V5.PD.avg[which(df$Group == 'C')]
treat.V5.PD.avg <- df$V5.PD.avg[which(df$Group == 'T')]

summary(control.BL.GE)
summary(treat.BL.GE)
boxplot(control.BL.GE, treat.BL.GE)

summary(control.BL.Calc.I)
summary(treat.BL.Calc.I)
boxplot(control.BL.Calc.I, treat.BL.Calc.I)

summary(control.V3.GE)
summary(treat.V3.GE)
boxplot(control.V3.GE, treat.V3.GE)

summary(control.V5.GE)
summary(treat.V5.GE)
boxplot(control.V5.GE, treat.V5.GE)

summary(control.V5.PD.avg)
summary(treat.V5.PD.avg)
boxplot(control.V5.PD.avg, treat.V5.PD.avg)

# Treatment improve the oral health: indeed, at visit V3 and V5 (after some treatments), the indeces are lower!
# Moreover patients with high indeces at the BL are treated (treated group has higher indeces at the BL) --> we could assess it by a permutational test
