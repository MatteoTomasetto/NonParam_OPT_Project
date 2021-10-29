#install.packages("medicaldata")
library(medicaldata)
df <- medicaldata::opt

################################################################################
##########################  DATA PRE-PROCESSING ################################
################################################################################

### MISSING VALUES CLEANING
blank = levels(df[,70])[1]

idx_NA = which(df[,70] == blank)
df[idx_NA,70] = NA
df[,70] = factor(df[,70])    # It removes " " from the levels
df[idx_NA,71] = NA           # NA if Gestional Age == time of the last visit 

for(i in 78:88){
idx_NA = which(df[,i] == blank)
df[idx_NA,i] = NA
df[,i] = factor(df[,i])
}

for(i in 89:100){
  df[,i] = factor(df[,i])
}

df[which(df[,103] == 100),103] = 0 # 0 if she does not miss any visits

for (i in 104:135){
  lev = levels(df[,i])
  df[,i] = as.numeric(lev[as.numeric(df[,i])])
}


### NEW VARIABLES
# Apgar_bin: Apgar_bin == 0 se neonato sano, Apgar_bin == 1 se neonato non sano
Apgar_bin = df[,74]
Apgar_bin[which(df[,74]>=7)] = 0  
Apgar_bin[which(df[,74]<7)] = 1 

# Antibodies: Antibodies = sum of the antibodies wrt different bacteria at visit 1
Antibodies = df[,104] + df[,105] +  df[,106] + df[,107] + df[,108] + df[,109] +  df[,110]
# Antibodies5: Antibodies5 = sum of the antibodies wrt different bacteria at visit 5
Antibodies5 = df[,120] + df[,121] +  df[,122] + df[,123] + df[,124] + df[,125] +  df[,126]

names(df)[145] = "Bacteria"
names(df)[155] = "Bacteria5"
names(df)[163] = "Bacteria%"
names(df)[171] = "Bacteria5%"

# idx_out <- [19, 34, 36]