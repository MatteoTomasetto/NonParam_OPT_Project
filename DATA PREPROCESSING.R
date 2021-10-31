################################################################################
##########################  DATA PRE-PROCESSING ################################
################################################################################

### MISSING VALUES CLEANING

#col 9 (Hisp)-> vuoto è diventato NA perchè il paziente non ha voluto rispondere
df$Hisp[which(df$Hisp !='No ' & df$Hisp!='Yes')]=NA
df$Hisp<-as.character(df$Hisp)
df$Hisp<- as.factor(df$Hisp)

#col 14 (BL.Diab.Type): NA sostituiti con "no diab" 
df$BL.Diab.Type<-as.character(df$BL.Diab.Type)
df$BL.Diab.Type[is.na(df$BL.Diab.Type)] <- 'no diab'
df$BL.Diab.Type<- as.factor(df$BL.Diab.Type)

#col 16 (Use.Tob)-> vuoto è diventato NA perchè il paziente non ha voluto rispondere
df$Use.Tob[which(df$Use.Tob!='No ' & df$Use.Tob!='Yes')]=NA
df$Use.Tob<-as.character(df$Use.Tob)
df$Use.Tob<- as.factor(df$Use.Tob)


#col 17 (BL.Cig.Day): ho sostituito NA con 0 per chi non fuma e lasciato na per chi 
#fuma o non ha risposto
df$BL.Cig.Day[which(df$Use.Tob=='No ')]<- 0


#col 18 (Use.Alc)-> vuoto è diventato NA perchè il paziente non ha voluto rispondere
df$Use.Alc[which(df$Use.Alc!='No ' & df$Use.Alc!='Yes')]=NA
df$Use.Alc<-as.character(df$Use.Alc)
df$Use.Alc<- as.factor(df$Use.Alc)


#col 19 (BL.Drks.Day): ho sostituito NA con 0 per chi non beve e lasciato na per chi 
#beve o non ha risposto
df$BL.Drks.Day[which(df$Use.Alc=='No ')]<- 0
#NB: ci sono alcuni che usano alcol ma che bevono 0 drink al giorno-> cosa fare?


#col 20 (Drug.Add)-> vuoto è diventato NA perchè il paziente non ha voluto rispondere
df$Drug.Add[which(df$Drug.Add!='No ' & df$Drug.Add!='Yes')]=NA
df$Drug.Add<-as.character(df$Drug.Add)
df$Drug.Add<- as.factor(df$Drug.Add)

#col 22 (N.prev.preg): ho sostituito NA con 0 per chi non ha avuto gravidanze e 
#lasciato na per chi ne ha avute e non ha risposto
df$N.prev.preg[which(df$Prev.preg=='No ')]<- 0

#col 23 (Live.PTB): na sostituiti con no previous pregnancies
df$Live.PTB<-as.character(df$Live.PTB)
df$Live.PTB[which(df$Prev.preg == 'No ')]<-'NoPrevPreg'
df$Live.PTB<-as.factor(df$Live.PTB)

#col 24 (Any.stillbirth): na sostituiti con no previous pregnancies
df$Any.stillbirth<-as.character(df$Any.stillbirth)
df$Any.stillbirth[which(df$Prev.preg == 'No ')]<-'NoPrevPreg'
df$Any.stillbirth<-as.factor(df$Any.stillbirth)

#col 25 (Spont.ab): na sostituiti con no previous pregnancies e vuoti sostituiti con
#NA (vuoto= ha avuto precedenti gravidanze (col 21), ma qui non ha risposto)
df$Spont.ab<-as.character(df$Spont.ab)
df$Spont.ab[which(df$Prev.preg == 'No ')]<-'NoPrevPreg'
df$Spont.ab[which(df$Spont.ab != 'Yes' & df$Spont.ab != 'No ' & df$Prev.preg == 'Yes')]<-NA
df$Spont.ab<-as.factor(df$Spont.ab)

#col 26 (Induced.ab): na sostituiti con no previous pregnancies e vuoti sostituiti con
#NA (vuoto= ha avuto precedenti gravidanze (col 21), ma qui non ha risposto)
df$Induced.ab<-as.character(df$Induced.ab)
df$Induced.ab[which(df$Prev.preg == 'No ')]<-'NoPrevPreg'
df$Induced.ab[which(df$Induced.ab != 'Yes' & df$Induced.ab != 'No ' & df$Prev.preg == 'Yes')]<-NA
df$Induced.ab<-as.factor(df$Induced.ab)


#col 27 (Any.live.ptb.sb.sp.ab.in.ab): inserito un terzo livello per le donne che 
#non hanno mai avuto una gravidanza (NA sostituito con NoPrevPreg)
df$Any.live.ptb.sb.sp.ab.in.ab<-as.character(df$Any.live.ptb.sb.sp.ab.in.ab)
df$Any.live.ptb.sb.sp.ab.in.ab[which(df$Prev.preg == 'No ')] <- 'NoPrevPreg'
df$Any.live.ptb.sb.sp.ab.in.ab<-as.factor(df$Any.live.ptb.sb.sp.ab.in.ab)


#col 28 (N.living.kids): NA sostituiti con 0 per donne che non hanno figli, 
#NA lasciati NA per donne con figli che non hanno specificato quanti erano vivi 
#at baseline
df$N.living.kids[which(df$Prev.preg == 'No ')] <- 0


#col 29 (Tx.comp.): vuoto-> 'withdrawn' (ritirato) if group=T 
#NA->'Notherapy' if group=C =>a seconda del valore della varaibile 3 (Group). 
df$Tx.comp.<-as.character(df$Tx.comp.)
df$Tx.comp.[which(df$Group == 'T' & df$Tx.comp. != 'Yes' & df$Tx.comp. != 'No ' & df$Tx.comp. != 'Und')] <-'Withdrawn'
df$Tx.comp.[which(df$Group == 'C' )] <-'NoTherapy'
df$Tx.comp.<-as.factor(df$Tx.comp.)


#col 30 (Local.anes): na->notherapy; no-> distinzione tra chi non ha fatto l'anestesia locale
#e chi si è ritirato dal trattamento (informazione alla colonna 29)
df$Local.anes<-as.character(df$Local.anes)
df$Local.anes[which(df$Group == 'C')]<-'NoTherapy'
df$Local.anes[which(df$Tx.comp. == 'Withdrawn' & df$Local.anes == 'No ')]<-'Withdrawn'
df$Local.anes<-as.factor(df$Local.anes)


#col 31 (Topical.Anest): na->notherapy; no-> distinzione tra chi non ha fatto la topical anestesia 
#e chi si è ritirato dal trattamento (informazione alla colonna 29)
df$Topical.Anest<-as.character(df$Topical.Anest)
df$Topical.Anest[which(df$Group == 'C')]<-'NoTherapy'
df$Topical.Anest[which(df$Tx.comp. == 'Withdrawn' & df$Topical.Anest == 'No ')]<-'Withdrawn'
df$Topical.Anest<-as.factor(df$Topical.Anest)


#col 32 (Tx.time): NA-> distinzione tra chi non ha fatto la terapia->0
#e chi si è ritirato ->NA
sum(is.na(df$Tx.time))
df$Tx.time[which(df$Group == 'C')]<-0
df$Tx.time[which(df$Tx.comp. == 'Withdrawn' & df$Group == 'T')]<-NA






library(corrplot)

################## dalla 33 alla 68

#33
levels(df$EDC.necessary.)
sum(is.na(df$EDC.necessary.))
df$EDC.necessary.[which(df$EDC.necessary.== blank)]<-NA
df$EDC.necessary.=factor(df$EDC.necessary.)


#34 non serve

#35 gi� apposto

#36 non serve

#dalla 37 alla 68 apposto: i missing data sono gi� come NA

# guardo le correlazioni delle colonne delle visite (SEZIONE PERIODONOTAL SUMMARIES)
x11()
corrplot(cor(df[,37:68],use='na.or.complete'))
title('correlation of 37-68')
#vediamo che ci sono dei chiari pattern => possiamo trascurare qualche covariata
#zoomiamo il grafico
x11()
corrplot(cor(df[,37:47],use='na.or.complete'))
x11()
pairs(df[,37:47])
x11()
pairs(df[,48:58])
# vediamo che le terne di variabili _avg,_4,_5 sono spesso fortemente correlate.
#questo vale anche per le altre (dalla 38 alla 68) quindi potremmo decidere di tenere solo avg?
#lo vediamo anche dal pairs

#altre correlazioni tra qt� calcolato tra terza e quinta visita
correlazioni<-cor(df[,37:68],use='na.or.complete')

correlazioni['V3.Calc.I','V5.Calc.I'] #0.99
correlazioni['V3.CAL.avg','V5.CAL.avg'] #0.88
correlazioni['V3.PD.avg','V5.PD.avg'] #0.86










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

