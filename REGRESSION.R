################################################################################
##########################  NON PARAMETRIC REGRESSION ################################
################################################################################

library(lme4)
library(ISLR2)
library(car)
library(mgcv)
library(rgl)
library(splines)
library(pbapply)
library(ISLR2)
library(car)
library(np)
library(splines)
library(fda)

# Function to build a logistic regression model with polynomials.
# INPUT:
#  y: response variable
#  x: covariate
#  deg: degree of the polynomial
#  grid_h: step size for the prediction grid
#  out: boolean, 1 if there are observations that should not be considered in regression
#  out_idx: indices of observations that should not be considered in regression
# OUTPUT: 
#  it returns the summary of the polynomial logistic regression till the 5th degree
logistic_1regression_poly = function(y, x, deg, grid_h = 0.5, out = 0, out_idx = NULL){
  
  if (out == 1){
    x <- x[-out_idx]
    y <- y[-out_idx]
  }
  
  x2 <- x[!is.na(y)]
  y2 <- y[!is.na(y)]
  y2 <- y2[!is.na(x2)]
  x2 <- x2[!is.na(x2)]
  y2 <- dummy(y2)
  
  
  m_list_logit = lapply(1:5,function(degree){glm(I(y2) ~ poly(x2,degree = degree),family = 'binomial')})
  
  age.grid = seq(range(x2)[1],range(x2)[2],by = grid_h)
  
  x11()
  preds = predict(m_list_logit[[deg]], list(x2 = age.grid), se=T)
  pfit = exp(preds$fit )/(1+ exp( preds$fit )) 
  se.bands.logit = cbind(preds$fit +2* preds$se.fit , preds$fit -2*
                           preds$se.fit)
  se.bands = exp(se.bands.logit)/(1+ exp(se.bands.logit))
  plot(x2, I(y2), xlim = range(age.grid), type = "n", ylim = c(0 ,.5))
  points (jitter(x2), I((y2)/4), cex = .5, pch = "|",
          col = " darkgrey ", main = 'Poly 4 Fit - Logistic')
  lines(age.grid, pfit, lwd = 2, col = " blue")
  matlines(age.grid, se.bands, lwd = 1, col = " blue",lty = 3)
  
  return(do.call(what = anova, c(list(test = "Chisq"), m_list_logit)))
  
}


# Function to build a logistic regression model with BSplines.
# INPUT:
#  y: response variable
#  x: covariate
#  knots: list of splines' knots 
#  deg: degree of the polynomial
#  grid_h: step size for the prediction grid
#  out: boolean, 1 if there are observations that should not be considered in regression
#  out_idx: indices of observations that should not be considered in regression
# OUTPUT: 
#  it returns the summary of the polynomial logistic regression till the 5th degree
logistic_1regression_bsplines = function(y, x, knots, deg, grid_h = 0.5, out = 0, out_idx = NULL){
  
  if (out == 1){
    x <- x[-out_idx]
    y <- y[-out_idx]
  }
  
  x2 <- x[!is.na(y)]
  y2 <- y[!is.na(y)]
  y2 <- y2[!is.na(x2)]
  x2 <- x2[!is.na(x2)]
  y2 <- dummy(y2)
  
  m_list_logit = lapply(1:5, function(degree){glm(I(y2) ~ bs(x2, knots = knots, degree = degree ),family='binomial')})
  
  age.grid = seq(range(x2)[1], range(x2)[2], by = grid_h)
  
  x11()
  preds = predict(m_list_logit[[deg]], list(x2 = age.grid), se=T)
  pfit = exp(preds$fit )/(1+ exp( preds$fit )) 
  se.bands.logit = cbind(preds$fit +2* preds$se.fit , preds$fit -2*
                           preds$se.fit)
  se.bands = exp(se.bands.logit)/(1+ exp(se.bands.logit))
  plot(x2, I(y2), xlim = range(age.grid), type = "n", ylim = c(0 ,.5))
  points (jitter(x2), I((y2)/4), cex = .5, pch = "|",
          col = " darkgrey ", main = 'Poly 4 Fit - Logistic')
  lines(age.grid, pfit, lwd = 2, col = " blue")
  matlines(age.grid, se.bands, lwd = 1, col = " blue",lty = 3)
  
  return(do.call(what = anova, c(list(test = "Chisq"), m_list_logit)))
  
}

# Function to build a logistic regression model with Smoothing Splines.
# INPUT:
#  y: response variable
#  x: covariate
#  grid_h: step size for the prediction grid
#  out: boolean, 1 if there are observations that should not be considered in regression
#  out_idx: indices of observations that should not be considered in regression
# OUTPUT: 
#  it returns the summary of the polynomial logistic regression till the 5th degree
logistic_1regression_Smoothingsplines = function(y, x, grid_h = 0.5, out = 0, out_idx = NULL){
  
  if (out == 1){
    x <- x[-out_idx]
    y <- y[-out_idx]
  }
  
  x2 <- x[!is.na(y)]
  y2 <- y[!is.na(y)]
  y2 <- y2[!is.na(x2)]
  x2 <- x2[!is.na(x2)]
  y2 <- dummy(y2)
  
  m_list_logit <- gam(I(y2) ~ s(x2, bs='cr'),method="GCV.Cp", family='binomial')
  
  age.grid = seq(range(x2)[1], range(x2)[2], by = grid_h)
  
  x11()
  preds = predict(m_list_logit, list(x2 = age.grid), se=T)
  pfit = exp(preds$fit )/(1+ exp( preds$fit )) 
  se.bands.logit = cbind(preds$fit +2* preds$se.fit , preds$fit -2*
                           preds$se.fit)
  se.bands = exp(se.bands.logit)/(1+ exp(se.bands.logit))
  plot(x2, I(y2), xlim = range(age.grid), type = "n", ylim = c(0 ,.5) , xlab = 'bl_GE', ylab = 'I(preterm)' )
  points (jitter(x2), I((y2)/4), cex = .5, pch = "|",
          col = " darkgrey ", main = 'Poly 4 Fit - Logistic')
  lines(age.grid, pfit, lwd = 2, col = " blue")
  matlines(age.grid, se.bands, lwd = 1, col = " blue",lty = 3)
  
  return(m_list_logit)
  
}

# Function to build a logistic regression model with Smoothing Splines.
# INPUT:
#  y: response variable
#  x: matrix of covariates (obtained by cbind(x1, x2, ...)). NB: its number of cols must be in [1,6]
#  grid_h: step size for the prediction grid
#  out: boolean, 1 if there are observations that should not be considered in regression
#  out_idx: indices of observations that should not be considered in regression
# OUTPUT: 
#  it returns the summary of the polynomial logistic regression till the 5th degree
logistic_1regression_Smoothingsplines_GAM = function(y, x, grid_h = 0.5, out = 0, out_idx = NULL){
  
  if (out == 1){
    x <- x[-out_idx,]
    y <- y[-out_idx]
  }
  
  data <- data.frame(y,x)
  p <- dim(x)[2]
  pp <- p+1
  dim(data)
  data <- na.omit(data)
  y2 <- data[,1]
  x2 <- data[,2:pp]
  
  y2 <- dummy(y2)
  
  if( dim(x2)[2] == 6){
  m_list_logit <- gam(I(y2) ~ s(x2[,1], bs = 'cr') + s(x2[,2], bs = 'cr') + s(x2[,3], bs = 'cr') + s(x2[,4], bs = 'cr') + s(x2[,5], bs = 'cr') + s(x2[,6], bs = 'cr'), family='binomial')}
  if( dim(x2)[2] == 5){
    m_list_logit <- gam(I(y2) ~ s(x2[,1], bs = 'cr') + s(x2[,2], bs = 'cr') + s(x2[,3], bs = 'cr') + s(x2[,4], bs = 'cr') + s(x2[,5], bs = 'cr'), family='binomial')}
  if( dim(x2)[2] == 4){
    m_list_logit <- gam(I(y2) ~ s(x2[,1], bs = 'cr') + s(x2[,2], bs = 'cr') + s(x2[,3], bs = 'cr') + s(x2[,4], bs = 'cr') , family='binomial')}
  if( dim(x2)[2] == 3){
    m_list_logit <- gam(I(y2) ~ s(x2[,1], bs = 'cr') + s(x2[,2], bs = 'cr') + s(x2[,3], bs = 'cr') , family='binomial')}
  if( dim(x2)[2] == 2){
    m_list_logit <- gam(I(y2) ~ s(x2[,1], bs = 'cr') + s(x2[,2], bs = 'cr') , family='binomial')}
  if( dim(x2)[2] == 1){
    m_list_logit <- gam(I(y2) ~ s(x2[,1], bs = 'cr') , family='binomial')}
  
  
  return(m_list_logit)
  
}


# we try to build a model for the dichotomous variable preterm 
# which indicates whether the pregnancy ended before 
# gestational age 37 weeks (259 days)
preterm <- df[,70]

#in order to do that, we used as covariates the periodontal variables from the first, the third and the fifth visits
bl_bacteria <- df[,136]
bl_ge <- df[,38]
bl_bleeding <- df[,39]
bl_pd <- df[,40]
bl_attachment <- df[,43]
bl_calculus <- df[,46]
bl_plaque <- df[,47]
v3_ge <- df[,48]
v3_bleeding <- df[,49]
v3_attachment <- df[,53]
v3_calculus <- df[,56]
v3_plaque <- df[,57]
v5_ge <- df[,58]
v5_bleeding <- df[,59]
v5_attachment <- df[,63]
v5_calculus <- df[,66]
v5_plaque <- df[,67]
v5_bacteria <- df[,146]
apgar <- df[,74]

# After many trials we find a good model which makes use of the Silness-Löe Gingival Index at baseline (first visit)
# We first applied a polynomial regression
logistic_1regression_poly(preterm,bl_ge,3,grid_h = 0.01)

# from the summary we see that also a 2 degree polynomial could be enough so we plot it
logistic_1regression_poly(preterm,bl_ge,2,grid_h = 0.01)

# we can see that there the probability of having a preterm pregnacy increases for 
# extreme values of the index
# that is because for extreme values we have very few data, so we pass to a bspline 
# regression in order to make a more precise prediction using knots were we have less data.
logistic_1regression_bsplines(preterm, bl_ge, knots = c(0.6,2.4), 3, grid_h = 0.001) 

# we notice that there could be a leverage effect on the rigth, due to the more extreme 
#point of bl_ge. So we try to refit the model without that points
x11()
boxplot(bl_ge)

range(bl_ge)
which(bl_ge >= 3)

logistic_1regression_bsplines(preterm, bl_ge, knots = c(0.6,2.4), 3, grid_h = 0.001, out = 1, out_idx = which(bl_ge >= 2.9)) 

# we notice that the probability of having a preterm birth increases for high value of bl_ge
# but we also notice that there is an (unexpected) oscillatory part in the middle that suggest that
# scores between (0.9, 1.1) have an higher probability wrt scores in (1.8, 2)
# that could be because at baseline, most of the patients to be treated were probably choosen
# among the ones with more severe periodontal conditions and so, after treatment they could 
# have improved their conditions and have prevent the preterm birth. 


# This is in some sense confirmed by the fact that this behaviour is less pronounced if we look at data from 
# the 5th visit. 
logistic_1regression_bsplines(preterm, v5_ge, knots = c(0.8,1.5,2.1), 3, grid_h = 0.001) 

# we also notice a strange behaviour of the curve of high value of v5_ge but we think is 
# due to the fact that there are very few data for high value of ge index in the last visit
which(v5_ge > 2.3)

# we try to use the smoothing splines 
mod <- logistic_1regression_Smoothingsplines(preterm, bl_ge, grid_h = 0.001)
anova(mod)

# we obtain a more smoothed function wrt bspline. Also, in this way we don't need to 
#take out the out_idx

# As a further step we tried to take into account the other covariates of the first 
#visit using a Generalized Additive Models (GAM)
x1 <- cbind(bl_ge, bl_attachment, bl_bacteria, bl_plaque, bl_calculus, bl_pd)

mod1 <- logistic_1regression_Smoothingsplines_GAM(preterm, x1, grid_h = 0.001)
anova(mod1)

# We see that is not a good model. Try to refit a model using only bl_ge and bl_pd,
# another model using bl_ge, bl_pd and bl_plaque
# and the last one using bl_ge and bl_plaque
x2 <- cbind(bl_ge, bl_pd )
x3 <- cbind(bl_ge, bl_bacteria)
x4 <- cbind(bl_ge, bl_plaque)
x5 <- cbind(bl_ge, bl_calculus)
x6 <- cbind(bl_ge, bl_bleeding)

mod2 <- logistic_1regression_Smoothingsplines_GAM(preterm, x2 , grid_h = 0.001)
anova(mod2)

mod3 <- logistic_1regression_Smoothingsplines_GAM(preterm, x3 , grid_h = 0.001, out = 1, out_idx = 94)
anova(mod3)

mod4 <- logistic_1regression_Smoothingsplines_GAM(preterm, x4 , grid_h = 0.001)
anova(mod4)

mod5 <- logistic_1regression_Smoothingsplines_GAM(preterm, x5 , grid_h = 0.001)
anova(mod5)

mod6 <- logistic_1regression_Smoothingsplines_GAM(preterm, x6 , grid_h = 0.001)
anova(mod6)


# now let's compare the model using the AIC index
AIC(mod, mod1, mod2, mod3, mod4, mod5, mod6)
