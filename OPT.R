install.packages("medicaldata")
library(medicaldata)
df <- medicaldata::opt
pairs(df[,c(37,68)])
df[,139]
sum(is.na(df[,139]))
# idx_out <- [19, 34, 36]