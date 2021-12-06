# USE THIS FILE AS TEMPLATE FOR THE PREPROCESSING STEP IN THE SCRIPTS

#install.packages("medicaldata")
library(medicaldata)
df <- medicaldata::opt
setwd("C:/Users/matte/Desktop/NonParam_OPT_Project")
#setwd("C:/Users/asia/Desktop/NonParam_OPT_Project")
#setwd("C:/Users/laura/Desktop/NonParam_OPT_Project")
source("DATA PREPROCESSING.R")

# NB: low birthweight if birthweight < 2500 g
#     pre-term birth if gestional age < 259 gg