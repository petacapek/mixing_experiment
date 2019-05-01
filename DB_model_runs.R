#This script is designed to use Dynamic Enenergy Budget theory numerical model
#to undesrtand patterns of carbon mass balance of the mixing experiment.
#The main idea is that the mass balance is biased by the change of conversion factor 
#between chloroform labile organic carbon and microbial biomass.

#Thus, the numerical model that accounts for changes in conversion factor is defined
#and calibrated against data from the experiment conducted with the same soils at the fine 
#temporal scale. The model parameters that explain the changes in conversion factor in time
#are used to fit the model on the data from mixing experiment.

#Model functions are defined in separate files and saved in the folder "Models". This script
#uploads these functions from source files. 

###############################################################################################
#required R labraries are uploaded
library(deSolve)
library(dplyr)
library(FME)
library(reshape)
library(ggplot2)
library(foreach)
library(doParallel)
library(DEoptim)
###############################################################################################

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#####################################Model calibration#########################################
########################################Jolanta data###########################################
#read the measured data
cal_data<-read.csv("./DB_concept/Hasan_Jolanta/hasan_jolanta.csv")

#Data contains aerobic and anaerobic incubations.
#Only aerobic incubation data are used. 
#The calibration is done separately for soils from Plesne catchment watershed (PL) 
#and Certovo catchment watershed (CT)

##################################Two different model formulations are tested and the best one is selected###################################################
#As source of soil organic carbon, two possibilities are tested - K2SO4 extractable DOC and water extractable DOC.
#Again, the best option is selected. The functions and results with K2SO4 extractable DOC are denoted
#by using letter "K" at the end of the function/output file name.  
#The functions and results with water extractable DOC are denoted
#by using letter "W" at the end of the function/output file name.  
###########################################K2SO4 extractable DOC##################################
#1. Death model - When Structures cannot be maintained, Structures are dying and releasing to DOC and Cres pools
source("Models/DB_cal_deathK.R")
PL_cal_deathK<-DB_cal_deathK(dataset = cal_data[(cal_data$Soil=="PL" & cal_data$Status=="A"), ])
CT_cal_deathK<-DB_cal_deathK(dataset = cal_data[(cal_data$Soil=="CT" & cal_data$Status=="A"), ])

PL_cal_deathK$pars
CT_cal_deathK$pars

PL_cal_deathK$goodness$Gfit
CT_cal_deathK$goodness$Gfit

ggplot(PL_cal_deathK$goodness$Yhat, aes(time, obs))+geom_point(cex=6)+geom_line(data=PL_cal_deathK$simul, aes(time, value))+facet_wrap(~variable, scales="free")
ggplot(CT_cal_deathK$goodness$Yhat, aes(time, obs))+geom_point(cex=6)+geom_line(data=CT_cal_deathK$simul, aes(time, value))+facet_wrap(~variable, scales="free")

#2. Respiration model - When Structures cannot be maintained, Structures are released as CO2 
source("Models/DB_cal_respK.R")
PL_cal_respK<-DB_cal_respK(dataset = cal_data[(cal_data$Soil=="PL" & cal_data$Status=="A"), ])
CT_cal_respK<-DB_cal_respK(dataset = cal_data[(cal_data$Soil=="CT" & cal_data$Status=="A"), ])

PL_cal_respK$pars
CT_cal_respK$pars

PL_cal_respK$goodness$Gfit
CT_cal_respK$goodness$Gfit

ggplot(PL_cal_respK$goodness$Yhat, aes(time, obs))+geom_point(cex=6)+geom_line(data=PL_cal_respK$simul, aes(time, value))+facet_wrap(~variable, scales="free")
ggplot(CT_cal_respK$goodness$Yhat, aes(time, obs))+geom_point(cex=6)+geom_line(data=CT_cal_respK$simul, aes(time, value))+facet_wrap(~variable, scales="free")

###########################################Water extractable DOC##################################
#1. Death model - When Structures cannot be maintained, Structures are dying and releasing to DOC and Cres pools
source("Models/DB_cal_deathW.R")
PL_cal_deathW<-DB_cal_deathW(dataset = cal_data[(cal_data$Soil=="PL" & cal_data$Status=="A"), ])
CT_cal_deathW<-DB_cal_deathW(dataset = cal_data[(cal_data$Soil=="CT" & cal_data$Status=="A"), ])

PL_cal_deathW$pars
CT_cal_deathW$pars

PL_cal_deathW$goodness$Gfit
CT_cal_deathW$goodness$Gfit

ggplot(PL_cal_deathW$goodness$Yhat, aes(time, obs))+geom_point(cex=6)+geom_line(data=PL_cal_deathW$simul, aes(time, value))+facet_wrap(~variable, scales="free")
ggplot(CT_cal_deathW$goodness$Yhat, aes(time, obs))+geom_point(cex=6)+geom_line(data=CT_cal_deathW$simul, aes(time, value))+facet_wrap(~variable, scales="free")

#2. Respiration model - When Structures cannot be maintained, Structures are released as CO2 
source("Models/DB_cal_respW.R")
PL_cal_respW<-DB_cal_respW(dataset = cal_data[(cal_data$Soil=="PL" & cal_data$Status=="A"), ])
CT_cal_respW<-DB_cal_respW(dataset = cal_data[(cal_data$Soil=="CT" & cal_data$Status=="A"), ])

PL_cal_respW$pars
CT_cal_respW$pars

PL_cal_respW$goodness$Gfit
CT_cal_respW$goodness$Gfit

ggplot(PL_cal_resp$goodness$Yhat, aes(time, obs))+geom_point(cex=6)+geom_line(data=PL_cal_resp$simul, aes(time, value))+facet_wrap(~variable, scales="free")
ggplot(CT_cal_resp$goodness$Yhat, aes(time, obs))+geom_point(cex=6)+geom_line(data=PL_cal_resp$simul, aes(time, value))+facet_wrap(~variable, scales="free")


##################################Possible isotope discrimination################################
#The simulations show that espetially 13C in microbial biomass is overestimated
#Therefore, the possible isotope fractionation during uptake is tested.
#New model structure with one additional parameter is defined and its correspondence with data is
#tested.
source("Models/DB_cal_deathK_D.R")
PL_cal_deathK_D<-DB_cal_deathK_D(dataset = cal_data[(cal_data$Soil=="PL" & cal_data$Status=="A"), ])
CT_cal_deathK_D<-DB_cal_deathK_D(dataset = cal_data[(cal_data$Soil=="CT" & cal_data$Status=="A"), ])

PL_cal_deathK_D$pars
CT_cal_deathK_D$pars

PL_cal_deathK_D$goodness$Gfit
CT_cal_deathK_D$goodness$Gfit

ggplot(PL_cal_deathK_D$goodness$Yhat, aes(time, obs))+geom_point(cex=6)+geom_line(data=PL_cal_deathK_D$simul, aes(time, value))+facet_wrap(~variable, scales="free")
ggplot(CT_cal_deathK_D$goodness$Yhat, aes(time, obs))+geom_point(cex=6)+geom_line(data=CT_cal_deathK_D$simul, aes(time, value))+facet_wrap(~variable, scales="free")

#Discrimination cannot explain the observed patterns. Therefore, it is not further considered

#################################################################################################
#################################################################################################
#######################################Model application#########################################
#################################################################################################
#################################################################################################
#Now, model is applied to understand the carbon mass balance in the mixing experiment
#Model parameters estimated in previous steps are used. K2SO4 extract is assumed as a soil carbon
#substrate and Structures that cannot be maintained are assumed to enter DOC and Cres pools.
#Mean estimates of fr and fs parameters are used without acknowledging the possible variability (these parameters are fixed)
#The rest of the parameters are estimated again but the lower and upper bounds of parameter
#estimates are used to define the parameter space.

#Measured data
mixing<-read.csv("data_mixing.csv")[c(1:160), ]

#Previously estimated parameters fr and fs
fractionsPL<-PL_cal_deathK$pars[c(10,11)]
fractionsCT<-CT_cal_deathK$pars[c(10,11)]

#Initial guess of parameters Vmax_glucose - Yu
initPL<-PL_cal_deathK$pars[c(1:9)]
initCT<-CT_cal_deathK$pars[c(1:9)]

#lower limits of parameter space
minpars<-pmin(summary(PL_cal_deathK$par_mcmc)["min", c(1:9)],
              summary(CT_cal_deathK$par_mcmc)["min", c(1:9)])
maxpars<-pmax(summary(PL_cal_deathK$par_mcmc)["max", c(1:9)],
              summary(CT_cal_deathK$par_mcmc)["max", c(1:9)])

#loading the testing function
# source("./Models/testf.R")
# 
# testf_out<-testf(dataset=mixing, fractionsPL = fractionsPL, fractionsCT = fractionsCT,
#                  initPL = initPL, initCT = initCT, minpar=minpars, maxpar = maxpars)

#loading the function
source("./Models/DB_mixing.R")

#defining number of cores
no_cors<-detectCores()
#creating cluster
cl<-makeCluster(no_cors)
#registering cluster
registerDoParallel(cl)

#function run
db_mixing_out<-DB_mixing(dataset=mixing, fractionsPL = fractionsPL, fractionsCT = fractionsCT,
                        initPL = initPL, initCT = initCT, minpar=minpars, maxpar = maxpars)

stopImplicitCluster()

#checking the goodness of correspondence
db_mixing_out$goodness1
db_mixing_out$goodness2
db_mixing_out$goodness3
db_mixing_out$goodness4

#The function cannot converge. Thus, litter horizons are run first separately as they seem to 
#converge without the problem
#loading the function
# source("./Models/DB_mixing_litter.R")
# 
# #defining number of cores
# no_cors<-detectCores()-1
# #creating cluster
# cl<-makeCluster(no_cors)
# #registering cluster
# registerDoParallel(cl)
# 
# #function run
# db_mixing_litter_out<-DB_mixing_litter(dataset=mixing, fractionsPL = fractionsPL, fractionsCT = fractionsCT,
#                          initPL = initPL, initCT = initCT, minpar=minpars, maxpar = maxpars)
# 
# stopImplicitCluster()
# 
# #checking the goodness of correspondence
# db_mixing_litter_out$goodness1
# db_mixing_litter_out$goodness2
# 
