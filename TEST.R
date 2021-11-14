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
  T0=abs(depthMedian(x,depth_params = list(method='Tukey'))-depthMedian(y,depth_params = list(method='Tukey')))  
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
    T_stat[perm] <- abs(depthMedian(x,depth_params = list(method='Tukey') - depthMedian(y,depth_params = list(method='Tukey')))
    
  }
  # p-value
  p_val <- sum(T_stat>=T0)/iter
  
  return(p_val)
}

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

##############################################################################
############################## T-test univariate #############################
##############################################################################

idx_Var = 71
idx_NA = which(is.na(df[,idx_Var]))
if (length(idx_NA)==0){
  Var_test = df[,idx_Var]
  Group_test = df[,3]
  x = Var_test[which(Group_test =='C')]
  y = Var_test[which(Group_test =='T')]
} else {
  Var_test = df[-idx_NA,idx_Var]
  Group_test = df[-idx_NA,3]
  x = Var_test[which(Group_test =='C')]
  y = Var_test[which(Group_test =='T')]
}

threshold = 2000 # 240 per delivery time, 2000 per birth weight
x = x[which(x <= threshold)]
y = y[which(y <= threshold)]

x11()
boxplot(x)
x11()
boxplot(y)

perm_t_test_mean(x,y)
perm_t_test_median(x,y)
perm_t_test_prop(x,y)

##############################################################################
############################## T-test bivariate ##############################
##############################################################################

idx_Var = c(71,72)
idx_NA = which(is.na(df[,idx_Var]))
if (length(idx_NA)==0){
  Var_test = df[,idx_Var]
  Group_test = df[,3]
  x = Var_test[which(Group_test =='C'),]
  y = Var_test[which(Group_test =='T'),]
} else {
  Var_test = df[-idx_NA,idx_Var]
  Group_test = df[-idx_NA,3]
  x = Var_test[which(Group_test =='C'),]
  y = Var_test[which(Group_test =='T'),]
}

# threshold = 2000 # 240 per delivery time, 2000 per birth weight
# x = x[which(x <= threshold)]
# y = y[which(y <= threshold)]

x11()
boxplot(x)
x11()
boxplot(y)

perm_t_test_mean(x,y)
perm_t_test_median(x,y)
perm_t_test_prop(x,y)