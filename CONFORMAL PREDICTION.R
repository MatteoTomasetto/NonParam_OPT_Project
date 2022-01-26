##############################################################################
############################### Pre-processing ###############################
##############################################################################

library(medicaldata)
df <- medicaldata::opt
source("DATA PREPROCESSING.R")

#LIBRARYES

library(dplyr) 
library(ggplot2)
library(knitr)
library(dbscan)
library(splines)
library(mgcv)
library(rgl)
#devtools::install_github(repo="ryantibs/conformal", subdir="conformalInference")
library(conformalInference)

head(df)
dim(df)
n <- dim(df)[1]
p <- dim(df)[2]
set.seed(06101998)


################################################################################
####                       FULL CONFORMAL PREDICTION                         ###
################################################################################

#1)Univariate Predictive Intervals variable for gestional age<240 in TREATMENT

# x = vector of univariate data
alpha <- 0.1 # set it properly 
#x.obs <- na.omit(df[,71])
#x.obs <- na.omit(df[which(df[,3]=='C'),71])
x.obs <- na.omit(df[which(df[,3]=='T' & df[,71]<240),71]) #GA nel gruppo treat<240
x.new.grid <- seq(min(x.obs) - 0.5*diff(range(x.obs)), max(x.obs) + 0.5*diff(range(x.obs)), length = 1000) 
p.value <- numeric(length(x.new.grid))

NC <- function(z.aug, i){
  # choose the wanted predictor:
  
  abs(z.aug[i] - mean(z.aug[-i]))        # sample mean as predictor
  # abs(z.aug[i] - median(z.aug[-i]))    # sample median as predictor
}

for(k in 1:length(x.new.grid)) {
  x.obs.aug <- c(x.obs, x.new.grid[k])
  scores <- numeric(length(x.obs.aug))
  for (i in 1:length(x.obs.aug)) {
    scores[i] <- NC(x.obs.aug, i)
  }
  p.value[k] <- sum(scores >= scores[length(x.obs.aug)])/(length(x.obs.aug))
}

# Plot the p-values
x11()
plot(x.new.grid, p.value, type='l', ylim=c(0,1))
abline(h=c(0,1))
abline(h=alpha, col='red', lty=2)
points(x.obs, numeric(length(x.obs)), pch=3)

# Compute the Prediction Interval
PI.grid <- x.new.grid[which(p.value >= alpha)]
PI <- c(min(PI.grid), max(PI.grid))
PI
plot(x.new.grid, p.value, type='l', ylim=c(0,1))
abline(h=c(0,1))
abline(h=alpha, col='red', lty=2)
points(x.obs, numeric(length(x.obs)), pch=3)
abline(v = PI, col='red')
points(PI.grid, numeric(length(PI.grid)), pch=16, col='red')
title('Severe cases of GA in Treatment group')


#2)Univariate Predictive Intervals variable gestional age<240 in Control
# x = vector of univariate data
alpha <- 0.1 # set it properly 
#x.obs <- na.omit(df[,71])
#x.obs <- na.omit(df[which(df[,3]=='C'),71])
x.obs <- na.omit(df[which(df[,3]=='C' & df[,71]<240),71]) #GA nel gruppo control<240
x.new.grid <- seq(min(x.obs) - 0.5*diff(range(x.obs)), max(x.obs) + 0.5*diff(range(x.obs)), length = 1000) 
p.value <- numeric(length(x.new.grid))

NC <- function(z.aug, i){
  # choose the wanted predictor:
  
  abs(z.aug[i] - mean(z.aug[-i]))        # sample mean as predictor
  # abs(z.aug[i] - median(z.aug[-i]))    # sample median as predictor
}

for(k in 1:length(x.new.grid)) {
  x.obs.aug <- c(x.obs, x.new.grid[k])
  scores <- numeric(length(x.obs.aug))
  for (i in 1:length(x.obs.aug)) {
    scores[i] <- NC(x.obs.aug, i)
  }
  p.value[k] <- sum(scores >= scores[length(x.obs.aug)])/(length(x.obs.aug))
}

# Plot the p-values
x11()
plot(x.new.grid, p.value, type='l', ylim=c(0,1))
abline(h=c(0,1))
abline(h=alpha, col='red', lty=2)
points(x.obs, numeric(length(x.obs)), pch=3)

# Compute the Prediction Interval
PI.grid <- x.new.grid[which(p.value >= alpha)]
PI <- c(min(PI.grid), max(PI.grid))
PI
plot(x.new.grid, p.value, type='l', ylim=c(0,1))
abline(h=c(0,1))
abline(h=alpha, col='red', lty=2)
points(x.obs, numeric(length(x.obs)), pch=3)
abline(v = PI, col='red')
points(PI.grid, numeric(length(PI.grid)), pch=16, col='red')
title('Severe cases of GA in Control group')

#3)Univariate Predictive Intervals variable birthweight<2000 in treatment
# x = vector of univariate data
alpha <- 0.1 # set it properly 
#x.obs <- na.omit(df[,71])
#x.obs <- na.omit(df[which(df[,3]=='C'),71])
x.obs <- na.omit(df[which(df[,3]=='T' & df[,72]<2000),72]) #weight nel gruppo treat<2000
x.new.grid <- seq(min(x.obs) - 0.5*diff(range(x.obs)), max(x.obs) + 0.5*diff(range(x.obs)), length = 1000) 
p.value <- numeric(length(x.new.grid))

NC <- function(z.aug, i){
  # choose the wanted predictor:
  
  abs(z.aug[i] - mean(z.aug[-i]))        # sample mean as predictor
  # abs(z.aug[i] - median(z.aug[-i]))    # sample median as predictor
}

for(k in 1:length(x.new.grid)) {
  x.obs.aug <- c(x.obs, x.new.grid[k])
  scores <- numeric(length(x.obs.aug))
  for (i in 1:length(x.obs.aug)) {
    scores[i] <- NC(x.obs.aug, i)
  }
  p.value[k] <- sum(scores >= scores[length(x.obs.aug)])/(length(x.obs.aug))
}

# Plot the p-values
x11()
plot(x.new.grid, p.value, type='l', ylim=c(0,1))
abline(h=c(0,1))
abline(h=alpha, col='red', lty=2)
points(x.obs, numeric(length(x.obs)), pch=3)

# Compute the Prediction Interval
PI.grid <- x.new.grid[which(p.value >= alpha)]
PI <- c(min(PI.grid), max(PI.grid))
PI
plot(x.new.grid, p.value, type='l', ylim=c(0,1))
abline(h=c(0,1))
abline(h=alpha, col='red', lty=2)
points(x.obs, numeric(length(x.obs)), pch=3)
abline(v = PI, col='red')
points(PI.grid, numeric(length(PI.grid)), pch=16, col='red')
title('Severe cases of Birthweight in Treatment group')

#4)Univariate Predictive Intervals variable birthweight<2000 in Control
# x = vector of univariate data
alpha <- 0.1 # set it properly 
#x.obs <- na.omit(df[,71])
#x.obs <- na.omit(df[which(df[,3]=='C'),71])
x.obs <- na.omit(df[which(df[,3]=='C' & df[,72]<2000),72]) #GA nel gruppo contr<2000
x.new.grid <- seq(min(x.obs) - 0.5*diff(range(x.obs)), max(x.obs) + 0.5*diff(range(x.obs)), length = 1000) 
p.value <- numeric(length(x.new.grid))

NC <- function(z.aug, i){
  # choose the wanted predictor:
  
  abs(z.aug[i] - mean(z.aug[-i]))        # sample mean as predictor
  # abs(z.aug[i] - median(z.aug[-i]))    # sample median as predictor
}

for(k in 1:length(x.new.grid)) {
  x.obs.aug <- c(x.obs, x.new.grid[k])
  scores <- numeric(length(x.obs.aug))
  for (i in 1:length(x.obs.aug)) {
    scores[i] <- NC(x.obs.aug, i)
  }
  p.value[k] <- sum(scores >= scores[length(x.obs.aug)])/(length(x.obs.aug))
}

# Plot the p-values
x11()
plot(x.new.grid, p.value, type='l', ylim=c(0,1))
abline(h=c(0,1))
abline(h=alpha, col='red', lty=2)
points(x.obs, numeric(length(x.obs)), pch=3)

# Compute the Prediction Interval
PI.grid <- x.new.grid[which(p.value >= alpha)]
PI <- c(min(PI.grid), max(PI.grid))
PI
plot(x.new.grid, p.value, type='l', ylim=c(0,1))
abline(h=c(0,1))
abline(h=alpha, col='red', lty=2)
points(x.obs, numeric(length(x.obs)), pch=3)
abline(v = PI, col='red')
points(PI.grid, numeric(length(PI.grid)), pch=16, col='red')
title('Severe cases of Birthweight in Control group')


# b) Multivariate Predictive Regions TREATMENT
# X = matrix of multivariate data (features are the columns)
alpha <- 0.1 # set it properly 
x.obs <- na.omit(df[which(df[,3]=='T' & df[,71]<240 & df[,72]<2000),c(71,72)])
x.obs <- scale(x.obs)
x1.new.grid <- seq(min(x.obs[,1]) - 0.25*diff(range(x.obs[,1])), max(x.obs[,1]) + 0.25*diff(range(x.obs[,1])), length = 20)
x2.new.grid <- seq(min(x.obs[,2]) - 0.25*diff(range(x.obs[,2])), max(x.obs[,2]) + 0.25*diff(range(x.obs[,2])), length = 20)
p.value <- matrix(nrow = length(x1.new.grid), ncol = length(x2.new.grid))

NC <- function(z.aug, i){
  # choose the wanted predictor:
  
  sum((z.aug[i,] - colMeans(z.aug[-i,]))^2) # Euclidean distance
  # as.numeric( as.matrix(z.aug[i,] - colMeans(z.aug[-i,])) %*% solve(cov(z.aug[-i,])) %*% as.matrix(t(z.aug[i,] - colMeans(z.aug[-i,]))) ) # Mahalanobis Distance
  # abs( z.aug[i,2] - sum(coefficients(lm(z.aug[-i,2]  ~ z.aug[-i,1]))*c(1, z.aug[i,1]))) # Regression if X has 2 columns
}

for(k in 1:length(x1.new.grid)) {
  for(h in 1:length(x2.new.grid)) {
    x.obs.aug <- rbind(x.obs, c(x1.new.grid[k],x2.new.grid[h]))
    scores <- numeric(dim(x.obs.aug)[1])
    for (i in 1:dim(x.obs.aug)[1]) {
      scores[i] <- NC(x.obs.aug, i)
    }
    p.value[k,h] <- sum(scores >= scores[dim(x.obs.aug)[1]])/(dim(x.obs.aug)[1])
    print(c(k,h))
  }
}

# Plot the p-values and the prediction region
x11()
image(x1.new.grid, x2.new.grid, p.value)
points(x.obs, pch=16)
contour(x1.new.grid, x2.new.grid, p.value, levels = alpha, add=T)
# Alternative plot
X=data.frame(x.obs)
x11()
x12_surface=expand.grid(x1.new.grid,x2.new.grid)
data_plot=cbind(p.value,x12_surface)
ggplot() + 
  scale_color_continuous()+
  geom_tile(data=data_plot, aes(Var1, Var2, fill= p.value)) +
  geom_point(data = data_plot, aes(Var1, Var2, fill= p.value)) +
  ylim(-5,5)+
  xlim(-5,5)
pred_set=x12_surface[p.value>alpha,]
poly_points=pred_set[chull(pred_set),]
ggplot() + 
  geom_tile(data=data_plot, aes(Var1, Var2, fill= p.value)) +
  geom_point(data=X, aes(Var1, Var2, fill= p.value)) + 
  geom_polygon(data=poly_points,aes(Var1,Var2),color='red',size=1,alpha=alpha)+
  ylim(-5,5)+
  xlim(-5,5)

# b) Multivariate Predictive Regions CONTROL
# X = matrix of multivariate data (features are the columns)
alpha <- 0.1 # set it properly 
x.obs <- na.omit(df[which(df[,3]=='C' & df[,71,]<240 & df[,72]<2000),c(71,72)])
x.obs <- scale(x.obs, center = false)
x1.new.grid <- seq(min(x.obs[,1]) - 0.25*diff(range(x.obs[,1])), max(x.obs[,1]) + 0.25*diff(range(x.obs[,1])), length = 20)
x2.new.grid <- seq(min(x.obs[,2]) - 0.25*diff(range(x.obs[,2])), max(x.obs[,2]) + 0.25*diff(range(x.obs[,2])), length = 20)
p.value <- matrix(nrow = length(x1.new.grid), ncol = length(x2.new.grid))

NC <- function(z.aug, i){
  # choose the wanted predictor:
  
  sum((z.aug[i,] - colMeans(z.aug[-i,]))^2) # Euclidean distance
  # as.numeric( as.matrix(z.aug[i,] - colMeans(z.aug[-i,])) %*% solve(cov(z.aug[-i,])) %*% as.matrix(t(z.aug[i,] - colMeans(z.aug[-i,]))) ) # Mahalanobis Distance
  # abs( z.aug[i,2] - sum(coefficients(lm(z.aug[-i,2]  ~ z.aug[-i,1]))*c(1, z.aug[i,1]))) # Regression if X has 2 columns
}

for(k in 1:length(x1.new.grid)) {
  for(h in 1:length(x2.new.grid)) {
    x.obs.aug <- rbind(x.obs, c(x1.new.grid[k],x2.new.grid[h]))
    scores <- numeric(dim(x.obs.aug)[1])
    for (i in 1:dim(x.obs.aug)[1]) {
      scores[i] <- NC(x.obs.aug, i)
    }
    p.value[k,h] <- sum(scores >= scores[dim(x.obs.aug)[1]])/(dim(x.obs.aug)[1])
    print(c(k,h))
  }
}

# Plot the p-values and the prediction region
x11()
image(x1.new.grid, x2.new.grid, p.value, xlim=c(0,300))
points(X, pch=16)
contour(x1.new.grid, x2.new.grid, p.value, levels = alpha, add=T)

