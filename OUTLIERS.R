##############################################################################
############################### Pre-processing ###############################
##############################################################################

library(medicaldata)
df <- medicaldata::opt
source("DATA PREPROCESSING.R")

# NB: low birthweight if birthweight < 2500 g
#     pre-term birth if gestional age < 259 gg

##############################################################################
############################## OUTLIERS' ANALYSIS ############################
##############################################################################


#divido le donne in control and treatment
group<-ifelse(df[,3] =='C',1,0)
control<-df[which(group== '1'),71:72]
treat<-df[which(group== '0'),71:72]
prem<-df[which(df[,70]=='Yes'),]
no_prem<-df[which(df[,70]=='No '),]

x11()
bagplot_matrix<-aplpack::bagplot.pairs(na.omit(df[,71:72]))
#bagplot_matrix<-aplpack::bagplot.pairs(na.omit(df[,c(38,71)]))
#bagplot_matrix<-aplpack::bagplot.pairs(na.omit(df[,c(46,71)]))
#bagplot_matrix<-aplpack::bagplot.pairs(na.omit(df[,c(47,71)]))

#bagplots
x11()
bp<-bagplot(na.omit(df[,71:72]))
#bp<-bagplot(na.omit(df[,c(38,71)]))
#bp<-bagplot(na.omit(df[,c(46,71)]))
#bp<-bagplot(na.omit(df[,c(47,71)]))

#select the outliers (values)
bp$pxy.outlier #sono 40

#extract the indexes of the outliers
j=1
idx=numeric(dim(bp$pxy.outlier)[1])
for(i in 1:823){
  if(is.na(df[i,71])==FALSE && is.na(df[i,72])==FALSE &&
     df[i,71]==bp$pxy.outlier[j,1] && df[i,72] == bp$pxy.outlier[j,2]){
    idx[j]=i
    j=j+1
  }
}


#make some comparison with histograms
for(i in 1:37){
  x11()
  par(mfrow=c(1,2))
  if(!is.na(df[idx,i])){
    barplot(table(df[idx,i]))
    title(main=names(df)[i],sub='outliers')
  }
  if(!is.na(df[,i])){
    barplot(table(df[,i]))
    title(main=names(df)[i],sub='tutti')
  }
}



for(i in 68:103){
  x11()
  par(mfrow=c(1,2))
  if(!is.na(df[idx,i])){
    barplot(table(df[idx,i]))
    title(main=names(df)[i], sub='outliers')
  }
  if(!is.na(df[,i])){
    barplot(table(df[,i]))
    title(main=names(df)[i], sub='tutti')
  }
}


#outliers considerando colonne 71 e 72-> nati prematuri 
#102) le pazienti idonee a poche visite (1,2,3,4) sono principalmente negli outliers 
#101) la maggior parte ha partecipato a 0 o 1 visita dopo la baseline 
#76) la maggior parte ha sperimentato qualche evento avverso (circa il 35% del totale)
#70) la maggior parte ha terminato la gravidanza prima dei 259 giorni (parto prematuro)(circa il 34% del totale)
#69) la maggior parte delle elective abortion e dei non-live birth sono tra gli outliers
#21)la maggior parte ha avuto precendenti gravidanze
#20)la maggior parte non si droga
#18)la maggior parte non beve alcool
#16)la maggior parte non fuma
#13)la maggior parte non ha il diabete
#11)la maggior parte ha pubblica asssitenza
#8) nessuna asiatica
#5)la maggior parte sono nere
#3)no distinzione tra control & treatment tra gli outliers

############################################################################################################
##COLONNE 38-71
#la maggior parte degli outlier sono idonee a poche visite (col 102)
#tra gli outliers la maggior parte delle donne ha fatto poche visite (andamento opposto rix a tutti i dati)
#la maggior parte degli outliers ha avuto eventi avversi (andamento opposto rix a tutti i dati)
#la maggior parte degli outliers ha partorito prematuramente (circa il 40% del totale)
#tutte le elective abortion sono tra gli outliers e quasi tutte i nati morti


#COLONNE 46-71
#la maggior parte degli outliers sono idonei a poche visite 
#nessuno degli outlier fa 5 visite
#circa il 40% del totale ha avuto eventi avversi
#tutti hanno partorito prematuramente
#tutte le elective abortion sono tra gli outliers e quasi tutte i nati morti

#gli indici hanno piu o meno tutti gli stessi outliers
#gli indici alla quinta visita hanno meno outliers che a quella baseline