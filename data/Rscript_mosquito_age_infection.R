# This is an R script to replicate the analyses from the article entitled 
# "Mosquito aging modulates the development, virulence and transmission potential of pathogens" 
# by Bernard Somé, Edwige Guissou, Dari F. Da, Quentin Richard, Marc Choisy, Koudraogo Bienvenue Yameogo,
# Domombabele FdS Hien, Rakiswende S Yerbanga, Georges Anicet Ouedraogo, Kounbobr R Dabiré,
# Ramsès Djidjou-Demasse, Anna Cohuet, Thierry Lefèvre (2023)
# Code developed by Thierry Lefevre: thierry.lefevre@ird.fr

library(emmeans)
library(glmmTMB)
library(car)
library(coxme)
library(MASS)
####################################################################################
### Fig 2A: analysis of the effect of age on oocyst prevalence in An. coluzzii ####
###################################################################################

###### load data
ooc<-read.table("fig2_oocyst_data.txt",header=T,stringsAsFactors = T)

### metadata: dataframe "ooc" has 4 columns:
# 	mosquito_id : a unique code for each dissected mosquito (n= 605)
# 	age_class : a 3-level categorical variable corresponding to each of the three age group (12, 8 and 4-day-old)
# 	isolate : a 4-level categorical variable corresponding to each of the four gametocyte carriers (replicates) 
#   oocyst : number of developed oocyst in each mosquito midgut

ooc$infection <- ifelse(ooc$oocyst > 0 ,1,0) ## creation of the binomial "infection" variable (0: uninfected, 1: infected with at least one oocyst)

mod1<-glmmTMB(infection~age_class+(1|isolate),family=binomial,data=ooc)
Anova(mod1)
summary(emmeans(mod1,pairwise~age_class,type="response"),infer=TRUE)

#######################################################################################
### Fig 2B: analysis of the effect of age on oocyst intensity in An. coluzzii #########
#######################################################################################

ooc_inf<-subset(ooc,oocyst>0) ## select infected specimens only

mod2<-glmmTMB(oocyst~age_class+(1|isolate),data=ooc_inf,family=truncated_nbinom2(link = "log"))
Anova(mod2)
summary(emmeans(mod2,pairwise~age_class,type="response"),infer=TRUE)

#######################################################################################
### Fig 2C: analysis of the effect of age on sporozoite prevalence in An. coluzzii ####
#######################################################################################

###### load data
spz<-read.table("fig2_sporozoite_data.txt",header=T,stringsAsFactors = T)

### metadata: dataframe "spz" has 5 columns:
# 	mosquito_id : a unique code for each dissected mosquito (n= 729)
# 	age_class : a 3-level categorical variable corresponding to each of the three age group (12, 8 and 4-day-old)
# 	isolate : a 4-level categorical variable corresponding to each of the four gametocyte carriers (replicates) 
#   Ct : number of cycle during the qPCR (this is a proxy of the quantity of parasite DNA in salivary glands)
#   positive: a binary variable corresponding to the infection status of mosquitoes (1: presence of, 0: absence of sporozoite)

mod3<-glmmTMB(positive~age_class+(1|isolate),family=binomial,data=spz)
Anova(mod3)
summary(emmeans(mod3,pairwise~age_class,type="response"),infer=TRUE)

#######################################################################################
### Fig 2D: analysis of the effect of age on sporozoite intensity in An. coluzzii #####
#######################################################################################

spz_inf<-subset(spz,positive=="1") ## select infected specimens only

mod4<-glmmTMB(log(Ct)~age_class+(1|isolate),data=spz_inf)
Anova(mod4)
summary(emmeans(mod4,pairwise~age_class,type="response"),infer=TRUE)

#######################################################################################
### Fig 3: analysis of the effect of age and infection on mosquito survival      ######
#######################################################################################

###### load data
surva<-read.table("fig3_survival_data.txt",header=T,stringsAsFactors = T)

### metadata: dataframe "surva" has 8 columns:
# 	mosquito_id : a unique code for each tracked mosquito (n= 657)
#   cup: a code for each paper cup containing mosquitoes (n=48 cups of ~15 mosquitoes (range = 4-18))
# 	age_class : a 3-level categorical variable corresponding to each of the three age group (12, 8 and 4-day-old)
# 	isolate : a 4-level categorical variable corresponding to each of the four gametocyte carriers (replicates) 
#   exposure : a two-level categorical variable corresponding to whether mosquitoes received an infectious (yes) or uninfectious blood-meal (no)
#   infection_status: a 3-level categorical variable corresponding to whether mosquitoes were uninfected controls, became infected upon exposure or remained uninfected upon exposure
#   daysPI: the day (time post-infection) of mosquito death
#   censor: censoring indicator (1 indicates that the response is a time at death, 0 indicates that the individual was alive when last seen). because all mosquitoes were followed until death, this indicator is 1 for every mosquito

mod5<-coxme(Surv(daysPI,censor)~age_class*infection_status+(1|isolate)+(1|cup),data=surva)
Anova(mod5)
mod6<-coxme(Surv(daysPI,censor)~age_class+(1|isolate)+(1|cup),data=surva)
summary(emmeans(mod6, pairwise~age_class, type="response"), infer = TRUE) # multiple pairwise comparison

#####analysis focused on the old mosquitoes #########
old<-subset(surva,age_class=="12-day-old")
mod7<-coxme(Surv(daysPI,censor)~infection_status+(1|isolate)+(1|cup),data=old)
Anova(mod7)
library(rockchalk)
old$infection<-combineLevels(old$infection_status,levs = c("Uninfected control", "Exposed-uninfected"), newLabel = c("uninfected") )
mod8<-coxme(Surv(daysPI,censor)~infection+(1|isolate)+(1|cup),data=old)
summary(emmeans(mod8, pairwise~infection, type="response"), infer = TRUE) # multiple pairwise comparison

#######################################################################################
### Fig 4: analysis of the effect of age on the parasite's EIP                    #####
#######################################################################################

###### load data
eip<-read.table("fig4_EIP_data.txt",header=T,stringsAsFactors = T)

### metadata: dataframe "eip" has 7 columns:
# 	mosquito_id : a unique code for each dissected mosquito (n= 802)
# 	age_class : a 2-level categorical variable corresponding to each of the two age group (12 and 4-day-old)
# 	isolate : a 3-level categorical variable corresponding to each of the three gametocyte carriers (replicates) 
#   dpi :  day (time post-infection) at which mosquitoes were dissected
#   oocyst : number of developed oocyst in each mosquito midgut
#   broken_oocyst: number of ruptured oocysts 
#   spz: a binomial variable corresponding to the presence (1) or absence (0) of sporozoite in mosquito salivary glands

eip<-subset(eip,dpi>6) ## Mosquitoes dissected at 6 dpi were used to check the success of infection (fig.S3) not to estimate eip
eip<-subset(eip,dpi<15) #### There was no more mosquitoes after day 14 for isolate G (The results for each isolate separately until 17 dpibm are detailed in Figure S4)
eip$total_oocyst <- eip$oocyst+eip$broken_oocyst 
eip<-subset(eip,total_oocyst>0) ## Uninfected mosquitoes (from which no estimates of EIP could be derived), need to be excluded

#### Figure 4A  Proportion of infected mosquitoes with ruptured oocysts #########

eip$prev_broken <- ifelse(eip$broken_oocyst > 0 ,1,0) ## creation of the binomial variable "prev_broken" (1: presence of at least one ruptured oocyst, 0: absence of ruptured oocyst)

mod9<-glmmTMB(prev_broken~age_class*dpi+(1|isolate),family=binomial,data=eip)
Anova(mod9)

#### Figure 4B  Fraction of ruptured oocysts #########

mod10<-glmmTMB(cbind(eip$oocyst,eip$broken_oocyst)~age_class*dpi+(1|isolate),family=binomial,data=eip)
Anova(mod10)

#### Figure 4C  Fraction of ruptured oocysts #########
eip$spz <- as.factor(eip$spz)
mod11<-glmmTMB(spz~age_class*dpi+(1|isolate),family=binomial,data=eip)
Anova(mod11)

#######################################################################################
### Fig 5: analysis of the effect of age on vectorial capacity                    #####
#######################################################################################

###### load data
VC<-read.table("fig5_vectorial_capacity_data.txt",header=T,stringsAsFactors = T)

### metadata: dataframe "VC" has 3 columns:
# 	scenario : a 3-level categorical variable corresponding to each scenario
# 	age_class : a 2-level categorical variable corresponding to each of the two age group (12 and 4-day-old)
#   vectorial_capacity: the response variable simulated from the mathematical model

VC_S1<-subset(VC,scenario=="Scenario 1")
mod12<-glm.nb(vectorial_capacity~age_class,data=VC_S1)
Anova(mod12)
##### END -----
