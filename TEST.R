##############################################################################
############################### Pre-processing ###############################
##############################################################################

library(medicaldata)
df <- medicaldata::opt
setwd("C:/Users/matte/Desktop/NonParam_OPT_Project")
source("DATA PREPROCESSING.R")

##############################################################################
################################ FUNCTIONS ###################################
##############################################################################

library(MASS)
library(rgl)
library(DepthProc)
library(hexbin)
library(packagefinder)
library(aplpack)
library(robustbase)
library(progress)

perm_t_test_mean = function(x,y,iter=1e3){
  T0=abs(mean(x)-mean(y))  
  T_stat=numeric(iter)
  x_pooled=c(x,y)
  n=length(x_pooled)
  n1=length(x)
  for(perm in 1:iter){
    # permutation:
    permutation <- sample(1:n)
    x_perm <- x_pooled[permutation]
    x1_perm <- x_perm[1:n1]
    x2_perm <- x_perm[(n1+1):n]
    # test statistic:
    T_stat[perm] <- abs(mean(x1_perm) - mean(x2_perm))
    
  }
  # p-value
  p_val <- sum(T_stat>=T0)/iter
  
  return(p_val)
}

perm_t_test_median = function(x,y,iter=1e3){
  T0=abs(median(x)-median(y))  
  T_stat=numeric(iter)
  x_pooled=c(x,y)
  n=length(x_pooled)
  n1=length(x)
  for(perm in 1:iter){
    # permutation:
    permutation <- sample(1:n)
    x_perm <- x_pooled[permutation]
    x1_perm <- x_perm[1:n1]
    x2_perm <- x_perm[(n1+1):n]
    # test statistic:
    T_stat[perm] <- abs(median(x1_perm) - median(x2_perm))
    
  }
  # p-value
  p_val <- sum(T_stat>=T0)/iter
  
  return(p_val)
}

perm_t_test_depth = function(x,y,iter=1e3){
  x_med=depthMedian(as.matrix(x),depth_params = list(method='Tukey'))
  y_med=depthMedian(as.matrix(y),depth_params = list(method='Tukey'))  
  T0=as.numeric((x_med-y_med) %*% (x_med-y_med))
  T_stat=numeric(iter)
  x_pooled=rbind(x,y)
  n=dim(x_pooled)[1]
  n1=dim(x)[1]
  for(perm in 1:iter){
    # permutation:
    permutation <- sample(1:n)
    x_perm <- x_pooled[permutation,]
    x1_perm <- x_perm[1:n1,]
    x2_perm <- x_perm[(n1+1):n,]
    # test statistic:
    x1_med=depthMedian(as.matrix(x1_perm),depth_params = list(method='Tukey'))
    x2_med=depthMedian(as.matrix(x2_perm),depth_params = list(method='Tukey'))  
    T_stat[perm]=as.numeric((x1_med-x2_med) %*% (x1_med-x2_med))
  }
  # p-value
  p_val <- sum(T_stat>=T0)/iter
  
  return(p_val)
}

perm_t_test_prop = function(x,y,iter=1e3){
  T_stat=numeric(iter)
  x_pooled=factor(c(x,y))
  n=length(x_pooled)
  n1=length(x)
  T0=abs(sum(x==levels(x)[2])/n1-sum(y==levels(y)[2])/(n-n1))  # Difference of percentages of the occurrences of a factor in the 2 groups 
  for(perm in 1:iter){
    # permutation:
    permutation <- sample(1:n)
    x_perm <- x_pooled[permutation]
    x1_perm <- x_perm[1:n1]
    x2_perm <- x_perm[(n1+1):n]
    # test statistic:
    T_stat[perm] <- abs(sum(x1_perm==levels(x1_perm)[1])/n1-sum(x2_perm==levels(x2_perm)[1])/(n-n1))
  }
  # p-value
  p_val <- sum(T_stat>=T0)/iter
  
  return(p_val)
}

perm_t_test_mean_multivariate = function(x,y,iter=1e3){
  T0=as.numeric((colMeans(x)-colMeans(y)) %*% (colMeans(x)-colMeans(y)))
  T_stat=numeric(iter)
  x_pooled=rbind(x,y)
  n=dim(x_pooled)[1]
  n1=dim(x)[1]
  for(perm in 1:iter){
    # permutation:
    permutation <- sample(1:n)
    x_perm <- x_pooled[permutation,]
    x1_perm <- x_perm[1:n1,]
    x2_perm <- x_perm[(n1+1):n,]
    # test statistic:
    T_stat[perm] <- as.numeric((colMeans(x1_perm)-colMeans(x2_perm)) %*% (colMeans(x1_perm)-colMeans(x2_perm)))  
  }
  # p-value
  p_val <- sum(T_stat>=T0)/iter
  
  return(p_val)
}

perm_anova = function(outcome,group,iter=1e3){
  fit <- aov(outcome ~ group)
  T0 <- summary(fit)[[1]][1,4]
  T_stat <- numeric(iter)
  n <- length(outcome)
  for(perm in 1:iter){
  # Permutation:
  permutation <- sample(1:n)
  outcome_perm <- outcome[permutation]
  fit_perm <- aov(outcome_perm ~ group)
  # Test statistic:
  T_stat[perm] <- summary(fit_perm)[[1]][1,4]
  }
  #pvalue
  p_val <- sum(T_stat>=T0)/iter
  return(p_val) 
}

perm_manova = function(outcome,group,iter=1e3){
  fit <- manova(as.matrix(outcome) ~ group)
  T0 <- -summary.manova(fit,test="Wilks")$stats[1,2]
  T_stat <- numeric(iter)
  n <- dim(outcome)[1]
  for(perm in 1:iter){
    # Permutation:
    permutation <- sample(1:n)
    outcome_perm <- outcome[permutation,]
    fit_perm <- manova(as.matrix(outcome_perm) ~ group)
    # Test statistic:
    T_stat[perm] <- -summary.manova(fit_perm,test="Wilks")$stats[1,2]
  }
  #pvalue
  p_val <- sum(T_stat>=T0)/iter
  return(p_val) 
}

##############################################################################
############################## T-test univariate #############################
##############################################################################

idx_Var = 34
idx_Group = 175
groups = levels(df[,idx_Group])

# discard NA
idx_NA_Var = which(is.na(df[,idx_Var]))
idx_NA_Group = which(is.na(df[,idx_Group]))
if (length(idx_NA_Var)==0 & length(idx_NA_Group)==0){
  Var_test = df[,idx_Var]
  Group_test = df[,idx_Group]
  x = Var_test[which(Group_test ==groups[1])]
  y = Var_test[which(Group_test ==groups[2])]
}
if(length(idx_NA_Var)!=0 & length(idx_NA_Group)==0) {
  Var_test = df[-idx_NA_Var,idx_Var]
  Group_test = df[-idx_NA_Var,idx_Group]
  x = Var_test[which(Group_test ==groups[1])]
  y = Var_test[which(Group_test ==groups[2])]
} 
if(length(idx_NA_Var)==0 & length(idx_NA_Group)!=0) {
  Var_test = df[-idx_NA_Group,idx_Var]
  Group_test = df[-idx_NA_Group,idx_Group]
  x = Var_test[which(Group_test ==groups[1])]
  y = Var_test[which(Group_test ==groups[2])]
} 
if(length(idx_NA_Var)!=0 & length(idx_NA_Group)!=0) {
  Var_test = df[-c(idx_NA_Var, idx_NA_Group),idx_Var]
  Group_test = df[-c(idx_NA_Var, idx_NA_Group),idx_Group]
  x = Var_test[which(Group_test ==groups[1])]
  y = Var_test[which(Group_test ==groups[2])]
}

# Set threshold for numerical variable if wanted
#threshold = 50
#x = x[which(x >= threshold)]
#y = y[which(y >= threshold)]

# Comparison with parametric case
#hist(x)
#hist(y)
#shapiro.test(x)
#shapiro.test(y)
#t.test(x,y)

# test
if(is.factor(df[,idx_Var])==TRUE){
  perm_t_test_prop(x,y) 
} else {
  x11()
  boxplot(x,main=groups[1])
  x11()
  boxplot(y,main=groups[2])
  
  perm_t_test_mean(x,y)
  perm_t_test_median(x,y)
  perm_anova(Var_test, Group_test)
}


##############################################################################
############################## T-test bivariate ##############################
##############################################################################

idx_Var = c(71,72)
idx_Group = 3
groups = levels(df[,idx_Group])

# discard NA
idx_NA_Var = c(which(is.na(df[,idx_Var[1]])),which(is.na(df[,idx_Var[2]])))
idx_NA_Var = sort(unique(idx_NA))
idx_NA_Group = which(is.na(df[,idx_Group]))
if (length(idx_NA_Var)==0 & length(idx_NA_Group)==0){
  Var_test = df[,idx_Var]
  Group_test = df[,idx_Group]
  x = Var_test[which(Group_test ==groups[1])]
  y = Var_test[which(Group_test ==groups[2])]
} 
if(length(idx_NA_Var)!=0 & length(idx_NA_Group)==0) {
  Var_test = df[-idx_NA_Var,idx_Var]
  Group_test = df[-idx_NA_Var,idx_Group]
  x = Var_test[which(Group_test ==groups[1])]
  y = Var_test[which(Group_test ==groups[2])]
} 
if(length(idx_NA_Var)==0 & length(idx_NA_Group)!=0) {
  Var_test = df[-idx_NA_Group,idx_Var]
  Group_test = df[-idx_NA_Group,idx_Group]
  x = Var_test[which(Group_test ==groups[1])]
  y = Var_test[which(Group_test ==groups[2])]
} 
if(length(idx_NA_Var)!=0 & length(idx_NA_Group)!=0) {
  Var_test = df[-c(idx_NA_Var, idx_NA_Group),idx_Var]
  Group_test = df[-c(idx_NA_Var, idx_NA_Group),idx_Group]
  x = Var_test[which(Group_test ==groups[1])]
  y = Var_test[which(Group_test ==groups[2])]
}


# Set threshold for both numerical variables if wanted
threshold1 = 240
threshold2 = 2000
x = x[which(x[,1] <= threshold1 & x[,2] <= threshold2),]
y = y[which(y[,1] <= threshold1 & y[,2] <= threshold2),]

# Set threshold for only 1 numerical variable if wanted
threshold = 250
idx_threshold = 1 # 1 or 2 in the bivariate case
x = x[which(x[,idx_threshold] <= threshold1),]
y = y[which(y[,idx_threshold] <= threshold1),]

# test
x11()
boxplot(x,main=groups[1])
x11()
boxplot(y,main=groups[2])
perm_t_test_mean_multivariate(x,y)
perm_t_test_depth(x,y,iter=250) # this could be slow
perm_manova(Var_test,Group_test)

##############################################################################
################################## RESULTS ###################################
##############################################################################

#    VAR                   THRESHOLD    GROUP                 TEST                PVALUE
# 1) N qualif. teeth (37); //           Control/Treat (3);    t_test_mean;        0.084/0.094;
# 2) Gestional age (71);   <240         Control/Treat (3);    t_test_mean/median; 0.071/0.061;  
# 3) BirthWeight (72);     <2000        Control/Treat (3);    t_test_mean/median; 0.128/0.07;  
# 4) Bacteria5 (155);      //           Control/Treat (3);    t_test_mean/median; 0/0;  
# 5) Bacteria5% (171);     //           Control/Treat (3);    t_test_mean/median; 0/0;  
# 6) GA&BW (71&72);        <240&<2000;  Control/Treat (3);    t_test_mean_multiv; 0.089;
# 7) Age (4);              //           Pre-term Y/N (70);    t_test_mean/median; 0.01/0.036;
# 8) Age (4);              //           Low-weight Y/N (175); t_test_mean/median; 0.152/0.54;
# 9) Black (5);            //           Pre-term Y/N (70);    t_test_prop;        0.085;
# 10) Black (5);           //           Low-weight Y/N (175); t_test_prop;        0.004;
# 11) Hisp (9);            //           Pre-term Y/N (70);    t_test_prop;        0.022;
# 12) Hisp (9);            //           Low-weight Y/N (175); t_test_prop;        0.138;
# 13) Pubblic Asstce (11); //           Pre-term Y/N (70);    t_test_prop;        0.092;
# 14) Pubblic Asstce (11); //           Low-weight Y/N (175); t_test_prop;        0.066;
# 15) Hypertension (12);   //           Pre-term Y/N (70);    t_test_prop;        0;
# 16) Hypertension (12);   //           Low-weight Y/N (175); t_test_prop;        0.003;
# 17) Diabetes (13);       //           Pre-term Y/N (70);    t_test_prop;        0.002;
# 18) Diabetes (13);       //           Low-weight Y/N (175); t_test_prop;        0.698;
# 19) BMI (15);            //           Pre-term Y/N (70);    t_test_mean;        0.005;
# 20) BMI (15);            //           Low-weight Y/N (175); t_test_mean;        0.49;
# 21) Use Tob (16);        //           Pre-term Y/N (70);    t_test_prop;        0.085;
# 22) Use Tob (16);        //           Low-weight Y/N (175); t_test_prop;        0.33;
# 23) Cigarettes/Day (17); //           Pre-term Y/N (70);    t_test_mean;        0.14;
# 24) Cigarettes/Day (17); //           Low-weight Y/N (175); t_test_mean;        0.008;
# 25) Use Alc (18);        //           Pre-term Y/N (70);    t_test_prop;        0.69;
# 26) Use Alc (18);        //           Low-weight Y/N (175); t_test_prop;        0.731;
# 27) Drinks/Day (19);     //           Pre-term Y/N (70);    t_test_mean;        0.56;
# 28) Drinks/Day (19);     //           Low-weight Y/N (175); t_test_mean;        0.633;
# 29) Drug Add (20);       //           Pre-term Y/N (70);    t_test_prop;        0.634;
# 30) Drug Add (20);       //           Low-weight Y/N (175); t_test_prop;        0.679;
# 31) Prev Preg (21);      //           Pre-term Y/N (70);    t_test_prop;        0.14;
# 32) Prev Preg (21);      //           Low-weight Y/N (175); t_test_prop;        0.422;
# 34) EDC necessary (33);  //           Pre-term Y/N (70);    t_test_prop;        0.741;
# 35) EDC necessary (33);  //           Low-weight Y/N (175); t_test_prop;        0.046;
# 33) Completed EDC (34);  //           Pre-term Y/N (70);    t_test_prop;        0.089;
# 34) Completed EDC (34);  //           Low-weight Y/N (175); t_test_prop;        0.633;
# 35) Bacteria5 (155);     //           Pre-term Y/N (70);    t_test_median;      0.067;
# 36) Bacteria5 (155);     //           Low-weight Y/N (175); t_test_median;      0.623;
# 37) Antibodies (173);    //           Pre-term Y/N (70);    t_test_mean/median; 0.036/0.015;
# 38) Antibodies (173);    //           Low-weight Y/N (175); t_test_mean/median; 0.0302/0.3;


##############################################################################
################################ CONCLUSIONS ################# ################
##############################################################################

# a) Our outcomes (Birthweight and gestional age) are differently distributed in the Control and Treatment groups,
# therefore there is an effect of the treatment on these risks (low weight and pre-term birth); this effect is positive
# since treated patients show higher birthweights and higher delivery times (see boxplots) ==> treatment reduces these risks.
# Notice thath, to see this difference, we have considered severe outcomes only!

# b) Some variables of dental health are naturally different in Control and Treat. grousp since treatment is required for severe patients 
# and it improves dental health reducing for example the amount of bacteria.

# c) We can see which features contribute to have pre-term birth and low birthweight from the pvalue above:
# e.g. diabetes and hypertension are higher in case of pre-term births, instead tobaccos/alcohol doesn't differ a lot for pre-term births but tobacco is significant for low birthweights;
# e.g. less bacteria and more antibodies prevent the pre-term birth risk.