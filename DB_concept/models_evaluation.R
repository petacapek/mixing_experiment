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

#Calculate the significance of the difference
##F test
###Plesne
####Number of measurements
ntPL<-5+5+4+4+5+5+4+4

(sum(PL_monod_out$goodness$Gfit$SSres)-sum(PL_cal_py$goodness$Gfit$SSres))*(ntPL - 12)/
  sum(PL_cal_py$goodness$Gfit$SSres)/(5)
pf(q=(sum(PL_monod_out$goodness$Gfit$SSres)-sum(PL_cal_py$goodness$Gfit$SSres))*(ntPL - 12)/
     sum(PL_cal_py$goodness$Gfit$SSres)/(5), 
   df1=5, 
   df2=ntPL-12, 
   lower.tail=F)

###Certovo
####Number of measurements
ntCT<-5+5+4+4+5+5+4+4

(sum(CT_monod_out$goodness$Gfit$SSres)-sum(CT_cal_py$goodness$Gfit$SSres))*(ntCT - 9)/
  sum(CT_cal_py$goodness$Gfit$SSres)/(3)
pf(q=(sum(CT_monod_out$goodness$Gfit$SSres)-sum(CT_cal_py$goodness$Gfit$SSres))*(ntCT - 9)/
     sum(CT_cal_py$goodness$Gfit$SSres)/(3), 
   df1=3, 
   df2=ntCT-9, 
   lower.tail=F)

####Summary table
S1 <- data.frame(Fval = c((sum(PL_monod_out$goodness$Gfit$SSres)-sum(PL_cal_py$goodness$Gfit$SSres))*(ntPL - 12)/
                        sum(PL_cal_py$goodness$Gfit$SSres)/(5),
                      (sum(CT_monod_out$goodness$Gfit$SSres)-sum(CT_cal_py$goodness$Gfit$SSres))*(ntCT - 9)/
                        sum(CT_cal_py$goodness$Gfit$SSres)/(3)),
                pval = c(pf(q=(sum(PL_monod_out$goodness$Gfit$SSres)-sum(PL_cal_py$goodness$Gfit$SSres))*(ntPL - 12)/
                              sum(PL_cal_py$goodness$Gfit$SSres)/(5), 
                            df1=5, 
                            df2=ntPL-12, 
                            lower.tail=F),
                         pf(q=(sum(CT_monod_out$goodness$Gfit$SSres)-sum(CT_cal_py$goodness$Gfit$SSres))*(ntCT - 11)/
                              sum(CT_cal_py$goodness$Gfit$SSres)/(4), 
                            df1=4, 
                            df2=ntCT-11, 
                            lower.tail=F)))
print(S1)

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

#Calculate the significance of the difference
##F test
###Plesne
####Number of measurements
ntPL<-5+5+4+4+5+5+4+4

(sum(PL_mend_out$goodness$Gfit$SSres)-sum(PL_cal_py$goodness$Gfit$SSres))*(ntPL - 12)/
  sum(PL_cal_py$goodness$Gfit$SSres)/(5)
pf(q=(sum(PL_mend_out$goodness$Gfit$SSres)-sum(PL_cal_py$goodness$Gfit$SSres))*(ntPL - 12)/
     sum(PL_cal_py$goodness$Gfit$SSres)/(5), 
   df1=5, 
   df2=ntPL-12, 
   lower.tail=F)

###Certovo
####Number of measurements
ntCT<-5+5+4+4+5+5+4+4

(sum(CT_mend_out$goodness$Gfit$SSres)-sum(CT_cal_py$goodness$Gfit$SSres))*(ntCT - 11)/
  sum(CT_cal_py$goodness$Gfit$SSres)/(4)
pf(q=(sum(CT_mend_out$goodness$Gfit$SSres)-sum(CT_cal_py$goodness$Gfit$SSres))*(ntCT - 11)/
     sum(CT_cal_py$goodness$Gfit$SSres)/(4), 
   df1=4, 
   df2=ntCT-11, 
   lower.tail=F)

####Summary table
S2 <- data.frame(Fval = c((sum(PL_mend_out$goodness$Gfit$SSres)-sum(PL_cal_py$goodness$Gfit$SSres))*(ntPL - 12)/
                            sum(PL_cal_py$goodness$Gfit$SSres)/(5),
                          (sum(CT_mend_out$goodness$Gfit$SSres)-sum(CT_cal_py$goodness$Gfit$SSres))*(ntCT - 11)/
                            sum(CT_cal_py$goodness$Gfit$SSres)/(5)),
                 pval = c(pf(q=(sum(PL_mend_out$goodness$Gfit$SSres)-sum(PL_cal_py$goodness$Gfit$SSres))*(ntPL - 12)/
                               sum(PL_cal_py$goodness$Gfit$SSres)/(5), 
                             df1=5, 
                             df2=ntPL-12, 
                             lower.tail=F),
                          pf(q=(sum(CT_mend_out$goodness$Gfit$SSres)-sum(CT_cal_py$goodness$Gfit$SSres))*(ntCT - 11)/
                               sum(CT_cal_py$goodness$Gfit$SSres)/(4), 
                             df1=4, 
                             df2=ntCT-11, 
                             lower.tail=F)))
print(S2)
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#############################################Figures ##########################################
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#Plesne - model simulations vs data
Pdata <- melt(cal_data[(cal_data$Soil=="PL" & cal_data$Status=="A"), 
                       c("time", "G12", "G13", "DOC12", "DOC13",
                         "Cmic12", "Cmic13", "CO212", "CO213")],
              id.vars = c("time"))

Pdata$facets = factor(Pdata$variable,
                      labels = c(
                        "{}^{12}~C~-~C[6]~H[12]~O[6]",
                        "{}^{13}~C~-~C[6]~H[12]~O[6]",
                        "{}^{12}~C~-~DOC",
                        "{}^{13}~C~-~DOC",
                        "{}^{12}~C~-~MBC",
                        "{}^{13}~C~-~MBC",
                        "{}^{12}~C~-~CO[2]",
                        "{}^{13}~C~-~CO[2]"
                      ))

DBsim<-PL_cal_py$simul[(PL_cal_py$simul$variable=="G_12C" |
                          PL_cal_py$simul$variable=="G_13C" |
                          PL_cal_py$simul$variable=="DOC_12C" |
                          PL_cal_py$simul$variable=="DOC_13C" |
                          PL_cal_py$simul$variable=="Cmic_12C" |
                          PL_cal_py$simul$variable=="Cmic_13C" |
                          PL_cal_py$simul$variable=="CO2_12C" |
                          PL_cal_py$simul$variable=="CO2_13C"), ]
DBsim$model<-"2-Pool MBC"
DBsim$facets = factor(DBsim$variable,
                      labels = c(
                        "{}^{12}~C~-~C[6]~H[12]~O[6]",
                        "{}^{13}~C~-~C[6]~H[12]~O[6]",
                        "{}^{12}~C~-~DOC",
                        "{}^{13}~C~-~DOC",
                        "{}^{12}~C~-~CO[2]",
                        "{}^{13}~C~-~CO[2]",
                        "{}^{12}~C~-~MBC",
                        "{}^{13}~C~-~MBC"
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
                        "{}^{12}~C~-~C[6]~H[12]~O[6]",
                        "{}^{13}~C~-~C[6]~H[12]~O[6]",
                        "{}^{12}~C~-~DOC",
                        "{}^{13}~C~-~DOC",
                        "{}^{12}~C~-~CO[2]",
                        "{}^{13}~C~-~CO[2]",
                        "{}^{12}~C~-~MBC",
                        "{}^{13}~C~-~MBC"
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
                           "{}^{12}~C~-~C[6]~H[12]~O[6]",
                           "{}^{13}~C~-~C[6]~H[12]~O[6]",
                           "{}^{12}~C~-~DOC",
                           "{}^{13}~C~-~DOC",
                           "{}^{12}~C~-~CO[2]",
                           "{}^{13}~C~-~CO[2]",
                           "{}^{12}~C~-~MBC",
                           "{}^{13}~C~-~MBC"
                         ))

Plsims<-rbind(DBsim, Monodsim, Mendsim)

texts<-data.frame(model = rep(c("MEND", "Monod", "2-Pool MBC"), each = 8),
                  variable = rep(c("G_12C", "G_13C", "DOC_12C", "DOC_13C",
                                   "CO2_12C", "CO2_13C", "Cmic_12C", "Cmic_13C"),
                                 times=3),
                  label = c(paste("R^2 ==", round(PL_mend_out$goodness$Gfit$R2[1], 2)),
                            paste("R^2 ==", round(PL_mend_out$goodness$Gfit$R2[2], 2)),
                            paste("R^2 ==", 0.00),
                            paste("R^2 ==", round(PL_mend_out$goodness$Gfit$R2[4], 2)),
                            paste("R^2 ==", round(PL_mend_out$goodness$Gfit$R2[5], 2)),
                            paste("R^2 ==", round(PL_mend_out$goodness$Gfit$R2[6], 2)),
                            paste("R^2 ==", round(PL_mend_out$goodness$Gfit$R2[7], 2)),
                            paste("R^2 ==", round(PL_mend_out$goodness$Gfit$R2[8], 2)),
                            paste("R^2 ==", round(PL_monod_out$goodness$Gfit$R2[1], 2)),
                            paste("R^2 ==", round(PL_monod_out$goodness$Gfit$R2[2], 2)),
                            paste("R^2 ==", 0.00),
                            paste("R^2 ==", round(PL_monod_out$goodness$Gfit$R2[4], 2)),
                            paste("R^2 ==", round(PL_monod_out$goodness$Gfit$R2[5], 2)),
                            paste("R^2 ==", round(PL_monod_out$goodness$Gfit$R2[6], 2)),
                            paste("R^2 ==", round(PL_monod_out$goodness$Gfit$R2[7], 2)),
                            paste("R^2 ==", round(PL_monod_out$goodness$Gfit$R2[8], 2)),
                            paste("R^2 ==", round(PL_cal_py$goodness$Gfit$R2[1], 2)),
                            paste("R^2 ==", round(PL_cal_py$goodness$Gfit$R2[2], 2)),
                            paste("R^2 ==", round(PL_cal_py$goodness$Gfit$R2[3], 2)),
                            paste("R^2 ==", round(PL_cal_py$goodness$Gfit$R2[4], 2)),
                            paste("R^2 ==", round(PL_cal_py$goodness$Gfit$R2[5], 2)),
                            paste("R^2 ==", round(PL_cal_py$goodness$Gfit$R2[6], 2)),
                            paste("R^2 ==", round(PL_cal_py$goodness$Gfit$R2[7], 2)),
                            paste("R^2 ==", round(PL_cal_py$goodness$Gfit$R2[8], 2))),
                  x=rep(c(55, 55, 14, 14, 55, 55, 55, 55,
                          55, 55, 14, 14, 55, 55, 55, 55,
                          55, 55, 19, 14, 49, 55, 56, 55)),
                  y=c(90, 2.5, 19.5, 0.225, 32, 0.65, 62, 0.85, 
                      75, 2.1, 18.5, 0.220, 22, 0.4, 59.5, 0.78,
                      60, 1.7, 17.5, 0.215, 12, 0.15, 57, 0.71))
texts$facets = factor(texts$variable,
                        labels = c(
                          "{}^{12}~C~-~MBC",
                          "{}^{13}~C~-~MBC",
                          "{}^{12}~C~-~CO[2]",
                          "{}^{13}~C~-~CO[2]",
                          "{}^{12}~C~-~DOC",
                          "{}^{13}~C~-~DOC",
                          "{}^{12}~C~-~C[6]~H[12]~O[6]",
                          "{}^{13}~C~-~C[6]~H[12]~O[6]"
                        ))

#5 per 12 inches
ggplot(Pdata, aes(time, value))+geom_point(cex=6, pch=21, fill="grey")+
  geom_line(data = Plsims, aes(time, value, colour=model, linetype=model), lwd=1.2)+
  facet_wrap(~facets, scales="free_y", labeller = label_parsed, ncol=4)+
  theme_min+xlab("Time (h)")+
  ylab(expression(paste("Pool size (", mu,"mol ", g(DW)^{-1}, ")")))+
  geom_text(data = texts, mapping = aes(x=x, y=y, label = label, colour=model), parse = T,
            fontface="bold", size = 5, show.legend = F)+
  scale_color_manual(values = c("black", "grey70", "grey30"))+
  scale_linetype_manual(values = c("solid", "F1", "dotdash"))+
  theme(legend.title = element_blank())
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#Certovo - model simulations vs data
Cdata <- melt(cal_data[(cal_data$Soil=="CT" & cal_data$Status=="A"), 
                       c("time", "G12", "G13", "DOC12", "DOC13",
                         "Cmic12", "Cmic13", "CO212", "CO213")],
              id.vars = c("time"))

Cdata$facets = factor(Cdata$variable,
                      labels = c(
                        "{}^{12}~C~-~C[6]~H[12]~O[6]",
                        "{}^{13}~C~-~C[6]~H[12]~O[6]",
                        "{}^{12}~C~-~DOC",
                        "{}^{13}~C~-~DOC",
                        "{}^{12}~C~-~MBC",
                        "{}^{13}~C~-~MBC",
                        "{}^{12}~C~-~CO[2]",
                        "{}^{13}~C~-~CO[2]"
                      ))

DBsim2<-CT_cal_py$simul[(CT_cal_py$simul$variable=="G_12C" |
                           CT_cal_py$simul$variable=="G_13C" |
                           CT_cal_py$simul$variable=="DOC_12C" |
                           CT_cal_py$simul$variable=="DOC_13C" |
                           CT_cal_py$simul$variable=="Cmic_12C" |
                           CT_cal_py$simul$variable=="Cmic_13C" |
                           CT_cal_py$simul$variable=="CO2_12C" |
                           CT_cal_py$simul$variable=="CO2_13C"), ]
DBsim2$model<-"2-Pool MBC"
DBsim2$facets = factor(DBsim2$variable,
                      labels = c(
                        "{}^{12}~C~-~C[6]~H[12]~O[6]",
                        "{}^{13}~C~-~C[6]~H[12]~O[6]",
                        "{}^{12}~C~-~DOC",
                        "{}^{13}~C~-~DOC",
                        "{}^{12}~C~-~CO[2]",
                        "{}^{13}~C~-~CO[2]",
                        "{}^{12}~C~-~MBC",
                        "{}^{13}~C~-~MBC"
                      ))
Monodsim2<-CT_monod_out$simul[(CT_monod_out$simul$variable=="G_12C" |
                                 CT_monod_out$simul$variable=="G_13C" |
                                 CT_monod_out$simul$variable=="DOC_12C" |
                                 CT_monod_out$simul$variable=="DOC_13C" |
                                 CT_monod_out$simul$variable=="Cmic_12C" |
                                 CT_monod_out$simul$variable=="Cmic_13C" |
                                 CT_monod_out$simul$variable=="CO2_12C" |
                                 CT_monod_out$simul$variable=="CO2_13C"), ]
Monodsim2$model<-"Monod"
Monodsim2$facets = factor(Monodsim2$variable,
                         labels = c(
                           "{}^{12}~C~-~C[6]~H[12]~O[6]",
                           "{}^{13}~C~-~C[6]~H[12]~O[6]",
                           "{}^{12}~C~-~DOC",
                           "{}^{13}~C~-~DOC",
                           "{}^{12}~C~-~CO[2]",
                           "{}^{13}~C~-~CO[2]",
                           "{}^{12}~C~-~MBC",
                           "{}^{13}~C~-~MBC"
                         ))
Mendsim2<-CT_mend_out$simul[(CT_mend_out$simul$variable=="G_12C" |
                               CT_mend_out$simul$variable=="G_13C" |
                               CT_mend_out$simul$variable=="DOC_12C" |
                               CT_mend_out$simul$variable=="DOC_13C" |
                               CT_mend_out$simul$variable=="Cmic_12C" |
                               CT_mend_out$simul$variable=="Cmic_13C" |
                               CT_mend_out$simul$variable=="CO2_12C" |
                               CT_mend_out$simul$variable=="CO2_13C"), ]
Mendsim2$model<-"MEND"
Mendsim2$facets = factor(Mendsim2$variable,
                        labels = c(
                          "{}^{12}~C~-~C[6]~H[12]~O[6]",
                          "{}^{13}~C~-~C[6]~H[12]~O[6]",
                          "{}^{12}~C~-~DOC",
                          "{}^{13}~C~-~DOC",
                          "{}^{12}~C~-~CO[2]",
                          "{}^{13}~C~-~CO[2]",
                          "{}^{12}~C~-~MBC",
                          "{}^{13}~C~-~MBC"
                        ))

CTsims<-rbind(DBsim2, Monodsim2, Mendsim2)

texts2<-data.frame(model = rep(c("MEND", "Monod", "2-Pool MBC"), each = 8),
                  variable = rep(c("G_12C", "G_13C", "DOC_12C", "DOC_13C",
                                   "CO2_12C", "CO2_13C", "Cmic_12C", "Cmic_13C"),
                                 times=3),
                  label = c(paste("R^2 ==", round(CT_mend_out$goodness$Gfit$R2[1], 2)),
                            paste("R^2 ==", round(CT_mend_out$goodness$Gfit$R2[2], 2)),
                            paste("R^2 ==", round(0.00, 2)),
                            paste("R^2 ==", round(0.00, 2)),
                            paste("R^2 ==", round(CT_mend_out$goodness$Gfit$R2[5], 2)),
                            paste("R^2 ==", round(CT_mend_out$goodness$Gfit$R2[6], 2)),
                            paste("R^2 ==", round(CT_mend_out$goodness$Gfit$R2[7], 2)),
                            paste("R^2 ==", round(CT_mend_out$goodness$Gfit$R2[8], 2)),
                            paste("R^2 ==", round(CT_monod_out$goodness$Gfit$R2[1], 2)),
                            paste("R^2 ==", round(CT_monod_out$goodness$Gfit$R2[2], 2)),
                            paste("R^2 ==", round(0.00, 2)),
                            paste("R^2 ==", round(0.00, 2)),
                            paste("R^2 ==", round(CT_monod_out$goodness$Gfit$R2[5], 2)),
                            paste("R^2 ==", round(CT_monod_out$goodness$Gfit$R2[6], 2)),
                            paste("R^2 ==", round(CT_monod_out$goodness$Gfit$R2[7], 2)),
                            paste("R^2 ==", round(CT_monod_out$goodness$Gfit$R2[8], 2)),
                            paste("R^2 ==", round(CT_cal_py$goodness$Gfit$R2[1], 2)),
                            paste("R^2 ==", round(CT_cal_py$goodness$Gfit$R2[2], 2)),
                            paste("R^2 ==", round(CT_cal_py$goodness$Gfit$R2[3], 2)),
                            paste("R^2 ==", round(0.00, 2)),
                            paste("R^2 ==", round(CT_cal_py$goodness$Gfit$R2[5], 2)),
                            paste("R^2 ==", round(CT_cal_py$goodness$Gfit$R2[6], 2)),
                            paste("R^2 ==", round(CT_cal_py$goodness$Gfit$R2[7], 2)),
                            paste("R^2 ==", round(CT_cal_py$goodness$Gfit$R2[8], 2))),
                  x=rep(c(55, 55, 53, 55, 55, 55, 55, 55,
                          55, 55, 53, 55, 55, 55, 55, 55,
                          55, 53, 58, 55, 55, 55, 55, 55)),
                  y=c(80, 2.2, 22.8, 0.25, 24, 0.55, 70, 0.91, 
                      65, 1.8, 21.8, 0.24, 16, 0.35, 65, 0.81,
                      50, 1.4, 20.8, 0.23, 8, 0.15, 60, 0.71))
texts2$facets = factor(texts2$variable,
                      labels = c(
                        "{}^{12}~C~-~MBC",
                        "{}^{13}~C~-~MBC",
                        "{}^{12}~C~-~CO[2]",
                        "{}^{13}~C~-~CO[2]",
                        "{}^{12}~C~-~DOC",
                        "{}^{13}~C~-~DOC",
                        "{}^{12}~C~-~C[6]~H[12]~O[6]",
                        "{}^{13}~C~-~C[6]~H[12]~O[6]"
                      ))

#5 per 12 inches
ggplot(Cdata, aes(time, value))+geom_point(cex=6, pch=21, fill="grey")+
  geom_line(data = CTsims, aes(time, value, colour=model, linetype=model), lwd=1.2)+
  facet_wrap(~facets, scales="free_y", labeller = label_parsed, ncol=4)+
  theme_min+xlab("Time (h)")+
  ylab(expression(paste("Pool size (", mu,"mol ", g(DW)^{-1}, ")")))+
  geom_text(data = texts2, mapping = aes(x=x, y=y, label = label, colour=model), parse = T,
            fontface="bold", size = 5, show.legend = F)+
  scale_color_manual(values = c("black", "grey70", "grey30"))+
  scale_linetype_manual(values = c("solid", "F1", "dotdash"))+
  theme(legend.title = element_blank())
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#Variability of the correction factor over time
kecvar<-data.frame(Soil = rep(c("Plesne", "Certovo"), each = 150),
                   Time = c(PL_cal_py$simul$time, CT_cal_py$simul$time),
                   kec = c(as.numeric(PL_cal_py$simul[(PL_cal_py$simul$variable=="Cmic_12C"), "value"])/
                             (as.numeric(PL_cal_py$simul[(PL_cal_py$simul$variable=="R_12C"), "value"])+
                             as.numeric(PL_cal_py$simul[(PL_cal_py$simul$variable=="S_12C"), "value"])),
                           as.numeric(CT_cal_py$simul[(CT_cal_py$simul$variable=="Cmic_12C"), "value"])/
                             (as.numeric(CT_cal_py$simul[(CT_cal_py$simul$variable=="R_12C"), "value"])+
                                as.numeric(CT_cal_py$simul[(CT_cal_py$simul$variable=="S_12C"), "value"]))))
#Variability of the CUE over time
plcue<-(PL_cal_py$simul[(PL_cal_py$simul$variable=="Cu_glucose"), "value"]*py_parsPL[["Ac_glucose"]]+
          PL_cal_py$simul[(PL_cal_py$simul$variable=="Cu_DOC"), "value"]*py_parsPL[["Ac_DOC"]]+
          PL_cal_py$simul[(PL_cal_py$simul$variable=="an"), "value"]*py_parsPL[["Yu"]])/
  (PL_cal_py$simul[(PL_cal_py$simul$variable=="Cu_glucose"), "value"]+
     PL_cal_py$simul[(PL_cal_py$simul$variable=="Cu_DOC"), "value"]+
     PL_cal_py$simul[(PL_cal_py$simul$variable=="an"), "value"])

ctcue<-(CT_cal_py$simul[(CT_cal_py$simul$variable=="Cu_glucose"), "value"]*py_parsCT[["Ac_glucose"]]+
          CT_cal_py$simul[(CT_cal_py$simul$variable=="Cu_DOC"), "value"]*py_parsCT[["Ac_DOC"]]+
          CT_cal_py$simul[(CT_cal_py$simul$variable=="an"), "value"]*py_parsCT[["Yu"]])/
  (CT_cal_py$simul[(CT_cal_py$simul$variable=="Cu_glucose"), "value"]+
     CT_cal_py$simul[(CT_cal_py$simul$variable=="Cu_DOC"), "value"]+
     CT_cal_py$simul[(CT_cal_py$simul$variable=="an"), "value"])


cuevar<-data.frame(Soil = rep(c("Plesne", "Certovo"), each = 150),
                   Time = c(PL_cal_py$simul$time, CT_cal_py$simul$time),
                   CUE = c(plcue, ctcue))
cuevar$kec<-kecvar$kec

grid.arrange(
  ggplot(subset(cuevar, Soil=="Plesne"), aes(x=Time))+
    geom_line(lwd=1.2, aes(y=CUE, color="CUE"), show.legend = F)+
    geom_line(lwd=1.2, aes(y=kec*4, color="kec"), show.legend = F)+
    theme_min+
    scale_y_continuous(sec.axis = sec_axis(~./4, 
                                           name = "Correction factor"))+
    facet_wrap(~Soil, scales="free")+
    ylab("Carbon use efficiency")+
    xlab("Time (h)")+
    theme(legend.title = element_blank(),
          axis.title.x = element_blank(),
          plot.margin = unit(c(0.1, 0.1, 0.1, 0.1), "in"))+
    scale_color_manual(values = c("black", "grey"),
                       name = " ", labels=c("CUE", "Correction factor")),
  ggplot(subset(cuevar, Soil=="Certovo"), aes(x=Time))+
    geom_line(lwd=1.2, aes(y=CUE, color="CUE"), show.legend=T)+
    geom_line(lwd=1.2, aes(y=kec*1.5, color="kec"), show.legend=T)+
    theme_min+
    scale_y_continuous(sec.axis = sec_axis(~./1.5, 
                                           name = "Correction factor"))+
    facet_wrap(~Soil, scales="free")+
    ylab("Carbon use efficiency")+
    xlab("Time (h)")+
    theme(legend.title = element_blank(),
          plot.margin = unit(c(0.1, 0.1, 0.1, 0.21), "in"),
          legend.position = c(0.4, 0.8))+
    scale_color_manual(values = c("black", "grey"),
                       name = " ", labels=c("CUE", "Correction factor")),
  ncol=1
)

kecvar2<-kecvar[(kecvar$Time==0.1 |
                   round(kecvar$Time, 0)==12 |
                   round(kecvar$Time, 0)==24 |
                   round(kecvar$Time, 0)==48 |
                   round(kecvar$Time, 0)==72), ]
kecvar2<-kecvar2[c(1, 3, 5, 7, 9, 10, 12, 14, 16, 18), ]

kecvar %>% group_by(Soil) %>% summarize(l = min(kec),
                                        m = max(kec))

