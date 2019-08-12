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
