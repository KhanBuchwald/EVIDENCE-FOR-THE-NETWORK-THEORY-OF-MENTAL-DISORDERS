## 1. Prepare SOPS

## Set working directory

setwd("")

## Read SOPS dataset

dat=read.table("sops01.txt", header=T)[2:4139,]

##  Order assessments from latest to earliest

dat2=dat[order(dat$interview_date, decreasing = T),]

##Select SOPS items required

dat3=dat2[,c(4,7,9,13,17,21,25,29:38)]

## Change empty string to NA

dat3[dat3==""]=NA

## Select observations with no missing data 

dat4=dat3[complete.cases(dat3),]

## Remove earlier observations

dat5=dat4[!duplicated(dat4$subjectkey),]

## Remove date column

dat6=dat5[,c(1,3:17)]

## Order according to participant ID

dat7=dat6[order(dat6$subjectkey),]

## Change variables to numeric

for(i in 2:ncol(dat7)){
  dat7[,i]=as.numeric(as.character(dat7[,i]))
}

## Select SOPS variables

dat7=dat7[,c(1:12)]

## Rewrite names of variables

names(dat7)[2:12]=c("P1","P2","P3","P4","P5","N1","N2","N3","N4","N5","N6")

## Change working directory

setwd("")

## Save dataset

write.csv(dat7, "NAPLS3_Data_Ready_For_Analysis.csv")

## Clear environment

rm(list = ls())

## 2. Prepare PANSS

## Change working directory

setwd("")

## Read data

dat=read.delim("panss01.txt", header=T)

##Select preintertvention data only

dat2=dat[dat$phase_ct=="Pre-Rand",]

## Order according to participant ID

dat2=dat2[order(dat2$subjectkey),]

##Select PANSS items

dat3=dat2[,c(4,239:268)]

## Change empty string to NA

dat3[dat3==""]=NA

## Select observations with no missing data 

dat4=dat3[complete.cases(dat3),]

## Select PANSS items

dat4=dat4[,c(1:15)]

## Rewrite names of variables

names(dat4)=c("subjectkey","P1","P2","P3","P4","P5","P6","P7","N1","N2","N3","N4","N5","N6","N7")

## Change working directory

setwd("")

## Save dataset

write.csv(dat4, "CATIE_Data_Ready_For_Analysis.csv")

## Clear environment

rm(list = ls())

## 3. Run structural equation models NAPLS3

library(lavaan)
library(bnlearn)
library(DiagrammeR)
library(rsvg)
library(DiagrammeRsvg)
library(semPlot)
library(lavaanPlot)

## Set working directory

setwd("")

## Read NAPLS3 preprared dataset

dat=read.csv("NAPLS3_Data_Ready_For_Analysis.csv", header = T)[,3:13]

##Convert all items to numeric

for(i in 1:11){
  dat[,i]=as.numeric(dat[,i])
} 

## Select and label 5 random groups

set.seed(999)
indices=cbind.data.frame(row=seq(1,795,1), rand=runif(795,0,10))
indices=indices[order(indices$rand),]
indices$group=c(rep(1,100), rep(2,100), rep(3,100),rep(4,100), rep(5,100), rep(0,295))

## Create five groups of train and test datasets

train.dat.1=dat[indices$row[indices$group==1],]
train.dat.2=dat[indices$row[indices$group==2],]
train.dat.3=dat[indices$row[indices$group==3],]
train.dat.4=dat[indices$row[indices$group==4],]
train.dat.5=dat[indices$row[indices$group==5],]

test.dat.1=dat[indices$row[indices$group!=1],]
test.dat.2=dat[indices$row[indices$group!=2],]
test.dat.3=dat[indices$row[indices$group!=3],]
test.dat.4=dat[indices$row[indices$group!=4],]
test.dat.5=dat[indices$row[indices$group!=5],]

## Run Hybrid Bayesian Network (On items)

Bn1=hc(train.dat.1)
Bn2=hc(train.dat.2)
Bn3=hc(train.dat.3)
Bn4=hc(train.dat.4)
Bn5=hc(train.dat.5)

## Obtain model specification

M1.1=modelstring(Bn1) 
M2.1=modelstring(Bn2)
M3.1=modelstring(Bn3)
M4.1=modelstring(Bn4)
M5.1=modelstring(Bn5)

## Change BN model string to Lavaan syntax function

changestring=function(M1){
  filtered_string <- gsub("\\[[^\\|]+\\]", "", M1)
  step1 <- gsub("]\\[", "]\n[", filtered_string)
  step2 <- gsub("\\[|\\]", "", step1)
  step3 <- gsub("\\|", "~", step2)
  final_format <- gsub(":", "+", step3)
  return(final_format)
}

## change model string

M1.2=changestring(M1.1)
M2.2=changestring(M2.1)
M3.2=changestring(M3.1)
M4.2=changestring(M4.1)
M5.2=changestring(M5.1)

## Write BN model in Lavaan syntax
## Run structural equation model on network

sys1.1=sem(model=M1.2, data=test.dat.1, estimator="MLR") 
sys1.2=sem(model=M2.2, data=test.dat.2, estimator="MLR")
sys1.3=sem(model=M3.2, data=test.dat.3, estimator="MLR")
sys1.4=sem(model=M4.2, data=test.dat.4, estimator="MLR")
sys1.5=sem(model=M5.2, data=test.dat.5, estimator="MLR")

## Write latent variable model in Lavaan syntax

Lavmod2.1='
Positive=~P1
Positive=~P2
Positive=~P3
Positive=~P4
Positive=~P5
Negative=~N1
Negative=~N2
Negative=~N3
Negative=~N4
Negative=~N5
Negative=~N6
'

Lavmod2.2='
Positive=~P1
Positive=~P2
Positive=~P3
Positive=~P4
Positive=~P5
Negative=~N1
Negative=~N2
Negative=~N3
Negative=~N4
Negative=~N5
Negative=~N6
Positive~~Negative
'

## Run structural equation model on latent variable model

sys2.1=cfa(model=Lavmod2.1, data=test.dat.1, estimator="MLR")
sys2.2=cfa(model=Lavmod2.1, data=test.dat.2, estimator="MLR")
sys2.3=cfa(model=Lavmod2.1, data=test.dat.3, estimator="MLR")
sys2.4=cfa(model=Lavmod2.1, data=test.dat.4, estimator="MLR")
sys2.5=cfa(model=Lavmod2.1, data=test.dat.5, estimator="MLR")
sys2.6=cfa(model=Lavmod2.2, data=test.dat.1, estimator="MLR") ## With covaraince term between latent variables
sys2.7=cfa(model=Lavmod2.2, data=test.dat.2, estimator="MLR") ## With covaraince term between latent variables
sys2.8=cfa(model=Lavmod2.2, data=test.dat.3, estimator="MLR") ## With covaraince term between latent variables
sys2.9=cfa(model=Lavmod2.2, data=test.dat.4, estimator="MLR") ## With covaraince term between latent variables
sys2.10=cfa(model=Lavmod2.2, data=test.dat.5, estimator="MLR") ## With covaraince term between latent variables

## Identify if covariance term between latent variables is significant

summary(sys2.6)$pe[12,]
summary(sys2.7)$pe[12,]
summary(sys2.8)$pe[12,]
summary(sys2.9)$pe[12,]
summary(sys2.10)$pe[12,]

## Write integrated model in Lavaan syntax

Lavmod3.1=paste(Lavmod2.2, M1.2, sep="\n")
Lavmod3.2=paste(Lavmod2.2, M2.2, sep="\n")
Lavmod3.3=paste(Lavmod2.2, M3.2, sep="\n")
Lavmod3.4=paste(Lavmod2.2, M4.2, sep="\n")
Lavmod3.5=paste(Lavmod2.2, M5.2, sep="\n")

## Run structural equation model on integrated model

sys3.1=cfa(model=Lavmod3.1, data=test.dat.1, estimator="MLR")
sys3.2=cfa(model=Lavmod3.2, data=test.dat.2, estimator="MLR")
sys3.3=cfa(model=Lavmod3.3, data=test.dat.3, estimator="MLR")
sys3.4=cfa(model=Lavmod3.4, data=test.dat.4, estimator="MLR")
sys3.5=cfa(model=Lavmod3.5, data=test.dat.5, estimator="MLR")

## Obtain likelihood, AIC, BIC, CFI, TLI, and RMSEA and write them into a csv file

Fit.Stat=cbind.data.frame(fit.sys1.1=fitmeasures(sys1.1), 
fit.sys1.2=fitmeasures(sys1.2),
fit.sys1.3=fitmeasures(sys1.3),
fit.sys1.4=fitmeasures(sys1.4),
fit.sys1.5=fitmeasures(sys1.5),
fit.sys1.6=fitmeasures(sys2.6),
fit.sys2.7=fitmeasures(sys2.7),
fit.sys2.8=fitmeasures(sys2.8),
fit.sys2.9=fitmeasures(sys2.9),
fit.sys2.10=fitmeasures(sys2.10),
fit.sys2.1=fitmeasures(sys3.1),
fit.sys3.2=fitmeasures(sys3.2),
fit.sys3.3=fitmeasures(sys3.3),
fit.sys3.4=fitmeasures(sys3.4),
fit.sys3.5=fitmeasures(sys3.5))

write.csv(Fit.Stat, "")

## Likelihood ratio test of latent variable model and integrated model (Sys2, Sys3)

anova(sys3.1,sys2.6,  test="LRT")
anova(sys3.2,sys2.7,  test="LRT")
anova(sys3.3,sys2.8,  test="LRT")
anova(sys3.4,sys2.9,  test="LRT")
anova(sys3.5,sys2.10,  test="LRT")

## Plot all intergrated SEMS

setwd("")

plot1=lavaanPlot(model=sys3.1, coefs=T, stars=c("regress"))
plot2=lavaanPlot(model=sys3.2, coefs=T, stars=c("regress"))
plot3=lavaanPlot(model=sys3.3, coefs=T, stars=c("regress"))
plot4=lavaanPlot(model=sys3.4, coefs=T, stars=c("regress"))
plot5=lavaanPlot(model=sys3.5, coefs=T, stars=c("regress"))

write(export_svg(plot1), file = "lavaanPlot1.svg")
write(export_svg(plot2), file = "lavaanPlot2.svg")
write(export_svg(plot3), file = "lavaanPlot3.svg")
write(export_svg(plot4), file = "lavaanPlot4.svg")
write(export_svg(plot5), file = "lavaanPlot5.svg")

rsvg_pdf("lavaanPlot1.svg", file = "lavaan_plot1.pdf", width = 600, height = 650)
rsvg_pdf("lavaanPlot2.svg", file = "lavaan_plot2.pdf", width = 600, height = 650)
rsvg_pdf("lavaanPlot3.svg", file = "lavaan_plot3.pdf", width = 600, height = 650)
rsvg_pdf("lavaanPlot4.svg", file = "lavaan_plot4.pdf", width = 600, height = 650)
rsvg_pdf("lavaanPlot5.svg", file = "lavaan_plot5.pdf", width = 600, height = 650)

## Clear environment

rm(list = ls())

## 4. Run structural equation models CATIE

## Set working directory

setwd("")

## Read CATIE preprared dataset

dat=read.csv("CATIE_Data_Ready_For_Analysis.csv", header = T)[,3:16]

## Convert all items to numeric

for(i in 1:14){
  dat[,i]=as.numeric(as.character(dat[,i]))
} 

## Select and label 5 random groups

set.seed(999)
indices=cbind.data.frame(row=seq(1,1446,1), rand=runif(1446,0,10))
indices=indices[order(indices$rand),]
indices$group=c(rep(1,180), rep(2,180), rep(3,180),rep(4,180), rep(5,180), rep(0,546))

train.dat.1=dat[indices$row[indices$group==1],]
train.dat.2=dat[indices$row[indices$group==2],]
train.dat.3=dat[indices$row[indices$group==3],]
train.dat.4=dat[indices$row[indices$group==4],]
train.dat.5=dat[indices$row[indices$group==5],]

test.dat.1=dat[indices$row[indices$group!=1],]
test.dat.2=dat[indices$row[indices$group!=2],]
test.dat.3=dat[indices$row[indices$group!=3],]
test.dat.4=dat[indices$row[indices$group!=4],]
test.dat.5=dat[indices$row[indices$group!=5],]

## Run Hybrid Bayesian Network (On items)

Bn1=hc(train.dat.1)
Bn2=hc(train.dat.2)
Bn3=hc(train.dat.3)
Bn4=hc(train.dat.4)
Bn5=hc(train.dat.5)

## Change BN model string to Lavaan syntax function

changestring=function(M1){
  filtered_string <- gsub("\\[[^\\|]+\\]", "", M1)
  step1 <- gsub("]\\[", "]\n[", filtered_string)
  step2 <- gsub("\\[|\\]", "", step1)
  step3 <- gsub("\\|", "~", step2)
  final_format <- gsub(":", "+", step3)
  return(final_format)
}

## Obtain model specification

M1.1=modelstring(Bn1) 
M2.1=modelstring(Bn2)
M3.1=modelstring(Bn3)
M4.1=modelstring(Bn4)
M5.1=modelstring(Bn5)

## change model string

M1.2=changestring(M1.1)
M2.2=changestring(M2.1)
M3.2=changestring(M3.1)
M4.2=changestring(M4.1)
M5.2=changestring(M5.1)

## Run structural equation model on network

sys1.1=sem(model=M1.2, data=test.dat.1, estimator="MLR") 
sys1.2=sem(model=M2.2, data=test.dat.2, estimator="MLR")
sys1.3=sem(model=M3.2, data=test.dat.3, estimator="MLR")
sys1.4=sem(model=M4.2, data=test.dat.4, estimator="MLR")
sys1.5=sem(model=M5.2, data=test.dat.5, estimator="MLR")

## Write latent variable model in Lavaan syntax

Lavmod2.1='
Positive=~P1
Positive=~P2
Positive=~P3
Positive=~P4
Positive=~P5
Positive=~P6
Positive=~P7
Negative=~N1
Negative=~N2
Negative=~N3
Negative=~N4
Negative=~N5
Negative=~N6
Negative=~N7
'

Lavmod2.2='
Positive=~P1
Positive=~P2
Positive=~P3
Positive=~P4
Positive=~P5
Positive=~P6
Positive=~P7
Negative=~N1
Negative=~N2
Negative=~N3
Negative=~N4
Negative=~N5
Negative=~N6
Negative=~N7
Positive~~Negative
'
sys2.1=cfa(model=Lavmod2.1, data=test.dat.1, estimator="MLR")
sys2.2=cfa(model=Lavmod2.1, data=test.dat.2, estimator="MLR")
sys2.3=cfa(model=Lavmod2.1, data=test.dat.3, estimator="MLR")
sys2.4=cfa(model=Lavmod2.1, data=test.dat.4, estimator="MLR")
sys2.5=cfa(model=Lavmod2.1, data=test.dat.5, estimator="MLR")
sys2.6=cfa(model=Lavmod2.2, data=test.dat.1, estimator="MLR") ## With covaraince term between latent variables
sys2.7=cfa(model=Lavmod2.2, data=test.dat.2, estimator="MLR") ## With covaraince term between latent variables
sys2.8=cfa(model=Lavmod2.2, data=test.dat.3, estimator="MLR") ## With covaraince term between latent variables
sys2.9=cfa(model=Lavmod2.2, data=test.dat.4, estimator="MLR") ## With covaraince term between latent variables
sys2.10=cfa(model=Lavmod2.2, data=test.dat.5, estimator="MLR") ## With covaraince term between latent variables

## Identify if covariance term between latent variables is significant

summary(sys2.6)$pe[15,]
summary(sys2.7)$pe[15,]
summary(sys2.8)$pe[15,]
summary(sys2.9)$pe[15,]
summary(sys2.10)$pe[15,]

## Write integrated model in Lavaan syntax

Lavmod3.1=paste(Lavmod2.2, M1.2, sep="\n")
Lavmod3.2=paste(Lavmod2.2, M2.2, sep="\n")
Lavmod3.3=paste(Lavmod2.2, M3.2, sep="\n")
Lavmod3.4=paste(Lavmod2.2, M4.2, sep="\n")
Lavmod3.5=paste(Lavmod2.2, M5.2, sep="\n")

## Run structural equation model on integrated model

sys3.1=cfa(model=Lavmod3.1, data=test.dat.1, estimator="MLR")
sys3.2=cfa(model=Lavmod3.2, data=test.dat.2, estimator="MLR")
sys3.3=cfa(model=Lavmod3.3, data=test.dat.3, estimator="MLR")
sys3.4=cfa(model=Lavmod3.4, data=test.dat.4, estimator="MLR")
sys3.5=cfa(model=Lavmod3.5, data=test.dat.5, estimator="MLR")

## Obtain likelihood, AIC, BIC, CFI, TLI, and RMSEA and write them into a csv file

Fit.Stat=cbind.data.frame(fit.sys1.1=fitmeasures(sys1.1), 
                          fit.sys1.2=fitmeasures(sys1.2),
                          fit.sys1.3=fitmeasures(sys1.3),
                          fit.sys1.4=fitmeasures(sys1.4),
                          fit.sys1.5=fitmeasures(sys1.5),
                          fit.sys1.6=fitmeasures(sys2.6),
                          fit.sys2.7=fitmeasures(sys2.7),
                          fit.sys2.8=fitmeasures(sys2.8),
                          fit.sys2.9=fitmeasures(sys2.9),
                          fit.sys2.10=fitmeasures(sys2.10),
                          fit.sys2.1=fitmeasures(sys3.1),
                          fit.sys3.2=fitmeasures(sys3.2),
                          fit.sys3.3=fitmeasures(sys3.3),
                          fit.sys3.4=fitmeasures(sys3.4),
                          fit.sys3.5=fitmeasures(sys3.5))

write.csv(Fit.Stat, "")

## Likelihood ratio test of latent variable model and integrated model (Sys2, Sys3)

anova(sys3.1,sys2.6,  test="LRT")
anova(sys3.2,sys2.7,  test="LRT")
anova(sys3.3,sys2.8,  test="LRT")
anova(sys3.4,sys2.9,  test="LRT")
anova(sys3.5,sys2.10,  test="LRT")

## Plot all intergrated SEMS

setwd("")

plot6=lavaanPlot(model=sys3.1, coefs=T, stars=c("regress"))
plot7=lavaanPlot(model=sys3.2, coefs=T, stars=c("regress"))
plot8=lavaanPlot(model=sys3.3, coefs=T, stars=c("regress"))
plot9=lavaanPlot(model=sys3.4, coefs=T, stars=c("regress"))
plot10=lavaanPlot(model=sys3.5, coefs=T, stars=c("regress"))

write(export_svg(plot6), file = "lavaanPlot6.svg")
write(export_svg(plot7), file = "lavaanPlot7.svg")
write(export_svg(plot8), file = "lavaanPlot8.svg")
write(export_svg(plot9), file = "lavaanPlot9.svg")
write(export_svg(plot10), file = "lavaanPlot10.svg")

rsvg_pdf("lavaanPlot6.svg", file = "lavaan_plot6.pdf", width = 600, height = 650)
rsvg_pdf("lavaanPlot7.svg", file = "lavaan_plot7.pdf", width = 600, height = 650)
rsvg_pdf("lavaanPlot8.svg", file = "lavaan_plot8.pdf", width = 600, height = 650)
rsvg_pdf("lavaanPlot9.svg", file = "lavaan_plot9.pdf", width = 600, height = 650)
rsvg_pdf("lavaanPlot10.svg", file = "lavaan_plot10.pdf", width = 600, height = 650)

## Clear environment

rm(list = ls())

## 5.  Descriptive Statistics NAPLS3 study

## Set working directory

setwd("")

## Read NAPLS3 preprared dataset

dat=read.csv("NAPLS3_Data_Ready_For_Analysis.csv", header = T)

## Change columns to numeric

for(i in 3:ncol(dat)){
  dat[,i]=as.numeric(dat[,i])
}

## Read demographics file

demo=read.delim("", header=T)

## Remove first row

demo2=demo[2:807,]

## Merge demographics with study participants

dat2=merge(dat, demo2, by="subjectkey", all.x=T)

## Frequency table of sex and ethnicity

table(dat2$sex)
table(dat2$ethnic_group)

## Proportions of each sex and ethnic group

round(prop.table(table(dat2$sex))*100,1)
round(prop.table(table(dat2$ethnic_group))*100,1)

## Calculate mean and sd for age

mean(as.numeric(dat2$interview_age)/12)
sd(as.numeric(dat2$interview_age)/12)

dat$Positive=rowSums(dat[,c(3:7)])
dat$Negative=rowSums(dat[,c(8:13)])
dat$Disorganised=rowSums(dat[,c(14:17)])

## Calculate summary statistics for each subscale

mean(dat$Positive)
mean(dat$Negative)
mean(dat$Disorganised)

## Calculate mean score on an item for each subscale

mean(dat$Positive)/5
mean(dat$Negative)/6
mean(dat$Disorganised)/4

## Calculate summary statistics for each subscale

median(dat$Positive)
median(dat$Negative)
median(dat$Disorganised)

sd(dat$Positive)
sd(dat$Negative)
sd(dat$Disorganised)

min(dat$Positive)
min(dat$Negative)
min(dat$Disorganised)

max(dat$Positive)
max(dat$Negative)
max(dat$Disorganised)

## Identify how many participants scored 3 or more on at least on item for each subscale

table(dat$P1>2|dat$P2>2|dat$P3>2|dat$P4>2|dat$P5>2)
table(dat$N1>2|dat$N2>2|dat$N3>2|dat$N4>2|dat$N5>2|dat$N6>2)
table(dat$D1>2|dat$D2>2|dat$D3>2|dat$D4>2)

## Identify how many participants scored 6 on at least on item for each subscale

table(dat$P1>5|dat$P2>5|dat$P3>5|dat$P4>5|dat$P5>5)
table(dat$N1>5|dat$N2>5|dat$N3>5|dat$N4>5|dat$N5>5|dat$N6>5)
table(dat$D1>5|dat$D2>5|dat$D3>5|dat$D4>5)

## Identify the proportion of participants who scored 6 on at least on item for each subscale

prop.table(table(dat$P1>5|dat$P2>5|dat$P3>5|dat$P4>5|dat$P5>5))
prop.table(table(dat$N1>5|dat$N2>5|dat$N3>5|dat$N4>5|dat$N5>5|dat$N6>5))
prop.table(table(dat$D1>5|dat$D2>5|dat$D3>5|dat$D4>5))

## Read dataset on diagnosis

diag=read.delim("", header=T)

## Removew first row

diag2=diag[2:nrow(diag),]

## Select study participants only

diag3=diag2[diag2$subjectkey%in%dat2$subjectkey,]

## Create diagnosis variable

diag3$X=rep(1, nrow(diag3))

## Change values of diagnosis variable to 0 if diagnostic code is 71.9 

diag3$X[diag3$scid_diag1=="71.9"]=0

## Change values of diagnosis variable to 0 if diagnostic code is an empty sting

diag3$X[diag3$scid_diag1==""]=0

## Order dataset according to date assessed and diagnosis status

diag3=diag3[order(diag3$interview_age, decreasing = T),]
diag3=diag3[order(diag3$X, decreasing = T),]

## Remove observations prior to and after a diagnosis was made and remove earlier observations in 
## those not diagnosed

diag4=diag3[!duplicated(diag3$subjectkey),]

## Frequency table of diagnosis

table(diag4$scid_diag1)

## Calculate proportion of people diagnosed

prop.table(table(diag4$scid_diag1))

## Clear environment

rm(list = ls())

## 6.  Descriptive Statistics CATIE study

## Change working directory

setwd("")

## Read CATIE preprared dataset

dat=read.csv("CATIE_Data_Ready_For_Analysis.csv", header=T)

## Calculate subscale scores

dat$positive=dat$P1+dat$P2+dat$P3+dat$P4+dat$P5+dat$P6+dat$P7
dat$negative=dat$N1+dat$N2+dat$N3+dat$N4+dat$N5+dat$N6+dat$N7
dat$general=dat$G1+dat$G2+dat$G3+dat$G4+dat$G5+dat$G6+dat$G7+dat$G8+dat$G9+dat$G10+dat$G11+dat$G12+dat$G13+dat$G14+
  dat$G15+dat$G16


## Calculate summary statistics for each subscale

mean(dat$positive)
median(dat$positive)
sd(dat$positive)
min(dat$positive)
max(dat$positive)

mean(dat$negative)
median(dat$negative)
sd(dat$negative)
min(dat$negative)
max(dat$negative)

mean(dat$general)
median(dat$general)
sd(dat$general)
min(dat$general)
max(dat$general)

## Create dataframe of subject ID for merge

dat2=as.data.frame(dat$subjectkey)

## Change column name of this dataset to "subjetID"

names(dat2)="subjectkey"

## Set working directory

setwd("")

## Read demographics dataset

demo=read.delim("demo01.txt", header=T)

## remove first row

demo2=demo[2:nrow(demo),]

## Select participants that are included in the prepared dataset

demo3=demo2[demo2$subjectkey%in%dat2$subjectkey,]

## Order by interview age

demo3=demo3[order(demo3$interview_age),]

## Selected earliest observation for each participant only

demo4=demo3[!duplicated(demo3$subjectkey),]

## Calculate mean age, SD of age, and number of missing values for age

mean(as.numeric(demo4$interview_age)/12, na.rm=T)
sd(as.numeric(demo4$interview_age)/12, na.rm=T)
table(is.na(demo4$interview_age))

## Frequency table for sex, proportions, and number of missing values of sex

table(demo4$sex)
prop.table(table(demo4$sex))
table(is.na(demo4$sex))

## Change empty string in Race to NA

demo4$race[demo4$race==""]=NA

## Change unknown in race to NA

demo4$race[demo4$race=="Unknown or not reported"]=NA

## Frequency table for race, proportions, and number of missing values of race 

table(demo4$race)
prop.table(table(demo4$race))
table(is.na(demo4$race))

## Read diagnosis dataset

SCID=read.delim("scid_ph01.txt", header = T)

## Remove first row

SCID2=SCID[2:nrow(SCID),]

## merge with prepared dataset

SCID3=merge(dat2,SCID2,by="subjectkey", all.x=T)

## Frequency table and proportions of people diagnosed in last 5 years

table(SCID3$scid01)
prop.table(table(SCID3$scid01))

## Clear environment

rm(list = ls())

## END