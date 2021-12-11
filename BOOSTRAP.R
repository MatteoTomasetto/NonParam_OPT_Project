##############################################################################
############################### Pre-processing ###############################
##############################################################################

library(medicaldata)
df <- medicaldata::opt
source("DATA PREPROCESSING.R")

# NB: low birthweight if birthweight < 2500 g
#     pre-term birth if gestional age < 259 gg

##############################################################################
################################### BOOSTRAP #################################
##############################################################################

### Bootstrap for Independent Samples ###
# Computing the bootstrap distribution of the difference of the two sample medians

#results for the difference of the medians in the two groups: significant differences: 103,101,67:64,62:54,52:48,44:43,35
#cio in cui differiscono sono parametri ortodontici e il numero di visite
group<-ifelse(df[,3] =='C',1,0)

set.seed(24021979)
B <- 10000

#possibili valori di i in corrispondenza delle variabili continue: 4,15,22,35-68,71:72,74:75,77,101:171
for(i in 71){
  control<-df[which(group== '1'),i]
  treat<-df[which(group== '0'),i]

  x1 <- na.omit(control)
  x2 <- na.omit(treat)

  # Plot data
  x11()
  par(mfrow=c(1,3))
  boxplot(x1, ylim=range(c(x1,x2)))
  title(main=names(df)[i], sub='control')
  boxplot(x2, ylim=range(c(x1,x2)))
  title(main=names(df)[i], sub='treatment')
  
  
  x1.obs <- x1
  x2.obs <- x2
  diff.Q2.obs <- quantile(x1, 0.50) - quantile(x2, 0.50)
  
  T.boot.diff.Q2 <- numeric(B)
  
  for(b in 1:B)
  {
    x1.b <- sample(x1.obs, replace = T)
    x2.b <- sample(x2.obs, replace = T)
    T.boot.diff.Q2[b] <- quantile(x1.b, 0.50) - quantile(x2.b, 0.50)
  }
  
  plot(ecdf(T.boot.diff.Q2), main='Sample Median control - Sample Median treatment')
  abline(v = diff.Q2.obs, lty=2)
  
  
  # RP intervals
  alpha <- 0.05
  
  right.quantile <- quantile(T.boot.diff.Q2, 1 - alpha/2)
  left.quantile  <- quantile(T.boot.diff.Q2, alpha/2)
  
  diff.Q2.obs
  right.quantile - diff.Q2.obs
  left.quantile  - diff.Q2.obs
  
  CI.RP <- c(diff.Q2.obs - (right.quantile - diff.Q2.obs), diff.Q2.obs - (left.quantile - diff.Q2.obs))
  CI.RP
  
  abline(v = CI.RP)
}



############################################################################################
#estimate of some quantities of varaibles 71 and 72 in the two groups (treatment and control)

bootstrap_quantile=function(x1, x2, col, perc, alpha){
  #x1,x2=the two groups
  #col=colonna della variabile di interesse
  #perc=percentuale del quantile da calcolare
  #alpha= level of the test
  x11()
  par(mfrow=c(1,2))
  boxplot(x1, ylim=range(c(x1,x2)))
  title(main=names(df)[col], sub='control')
  boxplot(x2, ylim=range(c(x1,x2)))
  title(main=names(df)[col], sub='treatment')

  x1.obs <- x1
  x2.obs <- x2

  
  med.obs1 <- quantile(x1,perc) 
  med.obs1
  med.obs2 <- quantile(x2, perc)
  med.obs2

  T.boot.med1 <- numeric(B)
  T.boot.med2 <- numeric(B)

  for(b in 1:B)
  {
    x1.b <- sample(x1.obs, replace = T)
    x2.b <- sample(x2.obs, replace = T)
    T.boot.med1[b] <- quantile(x1.b, perc) 
    T.boot.med2[b] <- quantile(x2.b, perc)
  }

  # RP intervals
  
  right.quantile1 <- quantile(T.boot.med1, 1 - alpha/2)
  left.quantile1  <- quantile(T.boot.med1, alpha/2)
  right.quantile2 <- quantile(T.boot.med2, 1 - alpha/2)
  left.quantile2  <- quantile(T.boot.med2, alpha/2)


  right.quantile1 - med.obs1
  left.quantile1  - med.obs1
  right.quantile2 - med.obs2
  left.quantile2  - med.obs2

  CI.RP1 <- c(med.obs1 - (right.quantile1 - med.obs1), med.obs1 - (left.quantile1 - med.obs1))
  CI.RP1
  CI.RP2 <- c(med.obs2 - (right.quantile2 - med.obs2), med.obs2 - (left.quantile2 - med.obs2))
  CI.RP2

  x11()
  par(mfrow=c(1,2))

  plot(ecdf(T.boot.med1), main='Sample Median control')
  abline(v = med.obs1, lty=2)
  abline(v = CI.RP1)

  plot(ecdf(T.boot.med2), main='Sample Median treatment')
  abline(v = med.obs2, lty=2)
  abline(v = CI.RP2)


}

bootstrap_mean=function(x1, x2, col, alpha){ 
  #x1,x2=the two groups
  #col=colonna della variabile di interesse
  #alpha=level of the test
  x11()
  par(mfrow=c(1,2))
  boxplot(x1, ylim=range(c(x1,x2)))
  title(main=names(df)[col], sub='control')
  boxplot(x2, ylim=range(c(x1,x2)))
  title(main=names(df)[col], sub='treatment')
  
  x1.obs <- x1
  x2.obs <- x2
  
  #mean
  med.obs1 <- mean(x1) 
  med.obs1
  med.obs2 <- mean(x2)
  med.obs2
  
  T.boot.med1 <- numeric(B)
  T.boot.med2 <- numeric(B)
  
  for(b in 1:B)
  {
    x1.b <- sample(x1.obs, replace = T)
    x2.b <- sample(x2.obs, replace = T)
    T.boot.med1[b] <- mean(x1.b) 
    T.boot.med2[b] <- mean(x2.b)
  }
  
  # RP intervals
  
  right.quantile1 <- quantile(T.boot.med1, 1 - alpha/2)
  left.quantile1  <- quantile(T.boot.med1, alpha/2)
  right.quantile2 <- quantile(T.boot.med2, 1 - alpha/2)
  left.quantile2  <- quantile(T.boot.med2, alpha/2)
  
  
  right.quantile1 - med.obs1
  left.quantile1  - med.obs1
  right.quantile2 - med.obs2
  left.quantile2  - med.obs2
  
  CI.RP1 <- c(med.obs1 - (right.quantile1 - med.obs1), med.obs1 - (left.quantile1 - med.obs1))
  CI.RP1
  CI.RP2 <- c(med.obs2 - (right.quantile2 - med.obs2), med.obs2 - (left.quantile2 - med.obs2))
  CI.RP2
  
  x11()
  par(mfrow=c(1,2))
  
  plot(ecdf(T.boot.med1), main='Sample Mean control')
  abline(v = med.obs1, lty=2)
  abline(v = CI.RP1)
  
  plot(ecdf(T.boot.med2), main='Sample Mean treatment')
  abline(v = med.obs2, lty=2)
  abline(v = CI.RP2)
  
  
}
# NB: low birthweight if birthweight < 2000 g
#     pre-term birth if gestional age < 240 gg

#two groups for gestional age
cd<-df[which(group== '1' & df[,71]<240),71]
td<-df[which(group== '0' & df[,71]<240),71]

#two groups for weights
cw<-df[which(group== '1' & df[,72]<2000),72]
tw<-df[which(group== '0' & df[,72]<2000),72]


bootstrap_quantile(cd,td,71,0.5,0.05)
bootstrap_quantile(cw,tw,72,0.5,0.05)
bootstrap_mean(cd,td,71,0.05)
bootstrap_mean(cw,tw,72,0.05)
