###############################################################################################
#required R libraries are uploaded
library(deSolve)
library(dplyr)
library(FME)
library(reshape)
library(ggplot2)
library(gridExtra)
library(foreach)
library(doParallel)
library(DEoptim)
###############################################################################################

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
###########################################DATA################################################
#read the measured data
cal_data<-read.csv("./DB_concept/Hasan_Jolanta/hasan_jolanta.csv")

#Data contains aerobic and anaerobic incubations.
#Only aerobic incubation data are used. 
#Soils from Plesne catchment watershed (PL) and Certovo catchment watershed (CT)
#are evaluated separately

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
###########################################DB model############################################
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#Model parameters are estimated in python
source("Models/DB_cal_death_py.R")
#read the parameters
py_parsPL <- as.numeric(read.csv("./DB_concept/Hasan_Jolanta/PL_parameters.csv", header = F))
py_parsCT <- as.numeric(read.csv("./DB_concept/Hasan_Jolanta/CT_parameters.csv", header = F))

PL_cal_py<-DB_cal_death_py(dataset = cal_data[(cal_data$Soil=="PL" & cal_data$Status=="A"), ],
                           py_pars = py_parsPL)
CT_cal_py<-DB_cal_death_py(dataset = cal_data[(cal_data$Soil=="CT" & cal_data$Status=="A"), ],
                           py_pars = py_parsCT)

PL_cal_py$goodness$Gfit
CT_cal_py$goodness$Gfit

ggplot(PL_cal_py$goodness$Yhat, aes(time, obs))+geom_point(cex=6)+geom_line(data=PL_cal_py$simul, aes(time, value))+facet_wrap(~variable, scales="free")
ggplot(CT_cal_py$goodness$Yhat, aes(time, obs))+geom_point(cex=6)+geom_line(data=CT_cal_py$simul, aes(time, value))+facet_wrap(~variable, scales="free")

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
###########################################Monod model#########################################
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#Model parameters are estimated in python
source("Models/monod_py.R")
#read the parameters
PLmonod <- as.numeric(read.csv("./DB_concept/Hasan_Jolanta/PL_monod.csv", header = F))
CTmonod <- as.numeric(read.csv("./DB_concept/Hasan_Jolanta/CT_monod.csv", header = F))

PL_monod_out<-Monod_py(dataset = cal_data[(cal_data$Soil=="PL" & cal_data$Status=="A"), ],
                           py_pars = PLmonod)
CT_monod_out<-Monod_py(dataset = cal_data[(cal_data$Soil=="CT" & cal_data$Status=="A"), ],
                           py_pars = CTmonod)

PL_monod_out$goodness$Gfit
CT_monod_out$goodness$Gfit

ggplot(PL_monod_out$goodness$Yhat, aes(time, obs))+
  geom_point(cex=6)+geom_line(data=PL_monod_out$simul, aes(time, value))+
  geom_line(data=PL_cal_py$simul, aes(time, value), colour="red")+
  facet_wrap(~variable, scales="free")
ggplot(CT_monod_out$goodness$Yhat, aes(time, obs))+
  geom_point(cex=6)+geom_line(data=CT_monod_out$simul, aes(time, value))+
  geom_line(data=CT_cal_py$simul, aes(time, value), colour="red")+
  facet_wrap(~variable, scales="free")

#Calculate the significance of the difference
##F test
###Plesne
####Number of measurements
ntPL<-5+5+4+4+5+5+4+4

(sum(PL_monod_out$goodness$Gfit$SSres)-sum(PL_cal_py$goodness$Gfit$SSres))*(ntPL - 12)/
  sum(PL_cal_py$goodness$Gfit$SSres)/(4)
pf(q=(sum(PL_monod_out$goodness$Gfit$SSres)-sum(PL_cal_py$goodness$Gfit$SSres))*(ntPL - 12)/
     sum(PL_cal_py$goodness$Gfit$SSres)/(4), 
   df1=4, 
   df2=ntPL-12, 
   lower.tail=F)

###Certovo
####Number of measurements
ntCT<-5+5+4+4+5+5+4+4

(sum(CT_monod_out$goodness$Gfit$SSres)-sum(CT_cal_py$goodness$Gfit$SSres))*(ntCT - 11)/
  sum(CT_cal_py$goodness$Gfit$SSres)/(3)
pf(q=(sum(CT_monod_out$goodness$Gfit$SSres)-sum(CT_cal_py$goodness$Gfit$SSres))*(ntCT - 11)/
     sum(CT_cal_py$goodness$Gfit$SSres)/(3), 
   df1=3, 
   df2=ntCT-11, 
   lower.tail=F)

####Summary table
S1 <- data.frame(Fval = c((sum(PL_monod_out$goodness$Gfit$SSres)-sum(PL_cal_py$goodness$Gfit$SSres))*(ntPL - 12)/
                        sum(PL_cal_py$goodness$Gfit$SSres)/(4),
                      (sum(CT_monod_out$goodness$Gfit$SSres)-sum(CT_cal_py$goodness$Gfit$SSres))*(ntCT - 12)/
                        sum(CT_cal_py$goodness$Gfit$SSres)/(4)),
                pval = c(pf(q=(sum(PL_monod_out$goodness$Gfit$SSres)-sum(PL_cal_py$goodness$Gfit$SSres))*(ntPL - 12)/
                              sum(PL_cal_py$goodness$Gfit$SSres)/(4), 
                            df1=4, 
                            df2=ntPL-12, 
                            lower.tail=F),
                         pf(q=(sum(CT_monod_out$goodness$Gfit$SSres)-sum(CT_cal_py$goodness$Gfit$SSres))*(ntCT - 12)/
                              sum(CT_cal_py$goodness$Gfit$SSres)/(4), 
                            df1=4, 
                            df2=ntCT-12, 
                            lower.tail=F)))
print(S1)

#Calculate the significance of the difference
##Log likelihood ratio test
###Plesne
####for each pool separately
pchisq(q=2*(PL_cal_py$goodness$Gfit$ll-PL_monod_out$goodness$Gfit$ll), 
   df=4, 
   lower.tail=F)
####all pools together
pchisq(q=2*(sum(PL_cal_py$goodness$Gfit$ll)-sum(PL_monod_out$goodness$Gfit$ll)), 
       df=4, 
       lower.tail=F)

###Certovo
####for each pool separately
pchisq(q=2*(CT_cal_py$goodness$Gfit$ll-CT_monod_out$goodness$Gfit$ll), 
       df=4, 
       lower.tail=F)
####all pools together
pchisq(q=2*(sum(CT_cal_py$goodness$Gfit$ll)-sum(CT_monod_out$goodness$Gfit$ll)), 
       df=4, 
       lower.tail=F)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
###########################################MEND model##########################################
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#Model parameters are estimated in python
source("Models/mend_py.R")
#read the parameters
PLmend <- as.numeric(read.csv("./DB_concept/Hasan_Jolanta/PL_mend.csv", header = F))
CTmend <- as.numeric(read.csv("./DB_concept/Hasan_Jolanta/CT_mend.csv", header = F))

PL_mend_out<-Mend_py(dataset = cal_data[(cal_data$Soil=="PL" & cal_data$Status=="A"), ],
                       py_pars = PLmend)
CT_mend_out<-Mend_py(dataset = cal_data[(cal_data$Soil=="CT" & cal_data$Status=="A"), ],
                       py_pars = CTmend)

PL_mend_out$goodness$Gfit
CT_mend_out$goodness$Gfit

ggplot(PL_mend_out$goodness$Yhat, aes(time, obs))+
  geom_point(cex=6)+geom_line(data=PL_mend_out$simul, aes(time, value))+
  geom_line(data=PL_cal_py$simul, aes(time, value), colour="red")+
  facet_wrap(~variable, scales="free")
ggplot(CT_mend_out$goodness$Yhat, aes(time, obs))+
  geom_point(cex=6)+geom_line(data=CT_mend_out$simul, aes(time, value))+
  geom_line(data=CT_cal_py$simul, aes(time, value), colour="red")+
  facet_wrap(~variable, scales="free")

#Calculate the significance of the difference
##F test
###Plesne
####Number of measurements
ntPL<-5+5+4+4+5+5+4+4

(sum(PL_mend_out$goodness$Gfit$SSres)-sum(PL_cal_py$goodness$Gfit$SSres))*(ntPL - 12)/
  sum(PL_cal_py$goodness$Gfit$SSres)/(3)
pf(q=(sum(PL_mend_out$goodness$Gfit$SSres)-sum(PL_cal_py$goodness$Gfit$SSres))*(ntPL - 12)/
     sum(PL_cal_py$goodness$Gfit$SSres)/(3), 
   df1=3, 
   df2=ntPL-12, 
   lower.tail=F)

###Certovo
####Number of measurements
ntCT<-5+5+4+4+5+5+4+4

(sum(CT_mend_out$goodness$Gfit$SSres)-sum(CT_cal_py$goodness$Gfit$SSres))*(ntCT - 12)/
  sum(CT_cal_py$goodness$Gfit$SSres)/(3)
pf(q=(sum(CT_mend_out$goodness$Gfit$SSres)-sum(CT_cal_py$goodness$Gfit$SSres))*(ntCT - 12)/
     sum(CT_cal_py$goodness$Gfit$SSres)/(3), 
   df1=3, 
   df2=ntCT-12, 
   lower.tail=F)

####Summary table
S2 <- data.frame(Fval = c((sum(PL_mend_out$goodness$Gfit$SSres)-sum(PL_cal_py$goodness$Gfit$SSres))*(ntPL - 12)/
                            sum(PL_cal_py$goodness$Gfit$SSres)/(3),
                          (sum(CT_mend_out$goodness$Gfit$SSres)-sum(CT_cal_py$goodness$Gfit$SSres))*(ntCT - 12)/
                            sum(CT_cal_py$goodness$Gfit$SSres)/(3)),
                 pval = c(pf(q=(sum(PL_mend_out$goodness$Gfit$SSres)-sum(PL_cal_py$goodness$Gfit$SSres))*(ntPL - 12)/
                               sum(PL_cal_py$goodness$Gfit$SSres)/(3), 
                             df1=3, 
                             df2=ntPL-12, 
                             lower.tail=F),
                          pf(q=(sum(CT_mend_out$goodness$Gfit$SSres)-sum(CT_cal_py$goodness$Gfit$SSres))*(ntCT - 12)/
                               sum(CT_cal_py$goodness$Gfit$SSres)/(3), 
                             df1=3, 
                             df2=ntCT-12, 
                             lower.tail=F)))
print(S2)

#Calculate the significance of the difference
##Log likelihood ratio test
###Plesne
####for each pool separately
pchisq(q=2*(PL_cal_py$goodness$Gfit$ll-PL_mend_out$goodness$Gfit$ll), 
       df=3, 
       lower.tail=F)
####all pools together
pchisq(q=2*(sum(PL_cal_py$goodness$Gfit$ll)-sum(PL_mend_out$goodness$Gfit$ll)), 
       df=3, 
       lower.tail=F)

###Certovo
####for each pool separately
pchisq(q=2*(CT_cal_py$goodness$Gfit$ll-CT_mend_out$goodness$Gfit$ll), 
       df=3, 
       lower.tail=F)
####all pools together
pchisq(q=2*(sum(CT_cal_py$goodness$Gfit$ll)-sum(CT_mend_out$goodness$Gfit$ll)), 
       df=3, 
       lower.tail=F)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#############################################Figures ##########################################
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#Plesne
Pdata <- melt(cal_data[(cal_data$Soil=="PL" & cal_data$Status=="A"), 
                       c("time", "G12", "G13", "DOC12", "DOC13",
                         "Cmic12", "Cmic13", "CO212", "CO213")],
              id.vars = c("time"))

Pdata$facets = factor(Pdata$variable,
                      labels = c(
                        "Glucose^{12}",
                        "Glucose^{13}",
                        "DOC^{12}",
                        "DOC^{13}",
                        "MBC^{12}",
                        "MBC^{13}",
                        "CO[2]^{12}",
                        "CO[2]^{13}"
                      ))

DBsim<-PL_cal_py$simul[(PL_cal_py$simul$variable=="G_12C" |
                          PL_cal_py$simul$variable=="G_13C" |
                          PL_cal_py$simul$variable=="DOC_12C" |
                          PL_cal_py$simul$variable=="DOC_13C" |
                          PL_cal_py$simul$variable=="Cmic_12C" |
                          PL_cal_py$simul$variable=="Cmic_13C" |
                          PL_cal_py$simul$variable=="CO2_12C" |
                          PL_cal_py$simul$variable=="CO2_13C"), ]
DBsim$model<-"DB"
DBsim$facets = factor(DBsim$variable,
                      labels = c(
                        "Glucose^{12}",
                        "Glucose^{13}",
                        "DOC^{12}",
                        "DOC^{13}",
                        "CO[2]^{12}",
                        "CO[2]^{13}",
                        "MBC^{12}",
                        "MBC^{13}"
                      ))
Monodsim<-PL_monod_out$simul[(PL_monod_out$simul$variable=="G_12C" |
                                PL_monod_out$simul$variable=="G_13C" |
                                PL_monod_out$simul$variable=="DOC_12C" |
                                PL_monod_out$simul$variable=="DOC_13C" |
                                PL_monod_out$simul$variable=="Cmic_12C" |
                                PL_monod_out$simul$variable=="Cmic_13C" |
                                PL_monod_out$simul$variable=="CO2_12C" |
                                PL_monod_out$simul$variable=="CO2_13C"), ]
Monodsim$model<-"Monod"
Monodsim$facets = factor(Monodsim$variable,
                      labels = c(
                        "Glucose^{12}",
                        "Glucose^{13}",
                        "DOC^{12}",
                        "DOC^{13}",
                        "CO[2]^{12}",
                        "CO[2]^{13}",
                        "MBC^{12}",
                        "MBC^{13}"
                      ))
Mendsim<-PL_mend_out$simul[(PL_mend_out$simul$variable=="G_12C" |
                               PL_mend_out$simul$variable=="G_13C" |
                               PL_mend_out$simul$variable=="DOC_12C" |
                               PL_mend_out$simul$variable=="DOC_13C" |
                               PL_mend_out$simul$variable=="Cmic_12C" |
                               PL_mend_out$simul$variable=="Cmic_13C" |
                               PL_mend_out$simul$variable=="CO2_12C" |
                               PL_mend_out$simul$variable=="CO2_13C"), ]
Mendsim$model<-"MEND"
Mendsim$facets = factor(Mendsim$variable,
                         labels = c(
                           "Glucose^{12}",
                           "Glucose^{13}",
                           "DOC^{12}",
                           "DOC^{13}",
                           "CO[2]^{12}",
                           "CO[2]^{13}",
                           "MBC^{12}",
                           "MBC^{13}"
                         ))
Plsims<-rbind(DBsim, Monodsim, Mendsim)
ggplot(Pdata, aes(time, value))+geom_point(cex=6, pch=21, fill="grey")+
  geom_line(data = Plsims, aes(time, value, colour=model))+
  facet_wrap(~facets, scales="free_y", labeller = label_parsed, ncol=4)+theme_min
