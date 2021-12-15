library(ggplot2)
library(scales)
library(medicaldata)

### Boxplot for numerical variables
# Define x and y with values of a variable for 2 different groups
# In the following, variables and groups names have to be changed manually

g1 <- rep("No", length(x))
g2 <- rep("Yes", length(y))
df1 <- as.data.frame(cbind(c(x,y),c(g1,g2)))
names(df1)[1] <- "Antibodies"
names(df1)[2] <- "Low_weight"
attach(df1)
df1[,1] <- as.numeric(df1[,1])

bp <- ggplot(df1, aes(x=Low_weight, y=Antibodies, fill=Low_weight)) + 
     geom_boxplot()+
     geom_jitter(shape=16, position=position_jitter(0.2))
bp



### Pie Chart for categorical variables
# In the following, variables and groups names have to be changed manually

df1 <- data.frame(
  Hypertension = c("No", "Yes"),
  value = c(92/(92+11), 11/(11+92)) # percentages for the pie chart
)
attach(df1)
blank_theme <- theme_minimal()+
  theme(
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    panel.border = element_blank(),
    panel.grid=element_blank(),
    axis.ticks = element_blank(),
  )

pie <- ggplot(df1, aes(x="", y=value, fill=Hypertension))+
       geom_bar(width = 1, stat = "identity")+ 
       coord_polar("y", start=0)+
       scale_fill_brewer(palette="Dark2")+
       ggtitle("Preterm: YES")+
       blank_theme+
       theme(axis.text.x=element_blank())+
       geom_text(aes(y = value/3 + c(0, cumsum(value)[-length(value)]), 
                label = percent(value/100)), size=5)
pie
