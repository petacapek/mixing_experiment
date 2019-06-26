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
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~Python parameters~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
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

source("Models/DB_cal_death_pyRS.R")
#read the parameters
py_parsPLRS <- as.numeric(read.csv("./DB_concept/Hasan_Jolanta/PL_parametersRS.csv", header = F))
py_parsCTRS <- as.numeric(read.csv("./DB_concept/Hasan_Jolanta/CT_parametersRS.csv", header = F))

PL_cal_pyRS<-DB_cal_death_pyRS(dataset = cal_data[(cal_data$Soil=="PL" & cal_data$Status=="A"), ],
                           py_pars = py_parsPLRS)
CT_cal_pyRS<-DB_cal_death_pyRS(dataset = cal_data[(cal_data$Soil=="CT" & cal_data$Status=="A"), ],
                           py_pars = py_parsCTRS)

PL_cal_pyRS$goodness$Gfit
CT_cal_pyRS$goodness$Gfit

ggplot(PL_cal_pyRS$goodness$Yhat, aes(time, obs))+geom_point(cex=6)+geom_line(data=PL_cal_pyRS$simul, aes(time, value))+facet_wrap(~variable, scales="free")
ggplot(CT_cal_pyRS$goodness$Yhat, aes(time, obs))+geom_point(cex=6)+geom_line(data=CT_cal_pyRS$simul, aes(time, value))+facet_wrap(~variable, scales="free")

#Marstorp 1999
source("Models/DB_M_death_py.R")
Mdata<-read.csv("./DB_concept/Marstorp/marstorp1999.csv")
py_parsM <- as.numeric(read.csv("./DB_concept/Marstorp/M_parameters.csv", header = F))

M_py<-DB_M_death_py(dataset = Mdata, py_pars = py_parsM)

M_py$goodness$Gfit

ggplot(M_py$goodness$Yhat, aes(time, obs))+geom_point(cex=3)+geom_line(data=M_py$simul, aes(time, value))+facet_wrap(~variable, scales="free")
ggplot(Mdata, aes(Time, CLC))+geom_point(cex=6)
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#


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

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~Python parameters~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
source("Models/DB_cal_resp_py.R")
#read the parameters
py_parsPL_resp <- as.numeric(read.csv("./DB_concept/Hasan_Jolanta/PL_parameters_resp.csv", header = F))
py_parsCT_resp <- as.numeric(read.csv("./DB_concept/Hasan_Jolanta/CT_parameters_resp.csv", header = F))

PL_cal_py_resp<-DB_cal_resp_py(dataset = cal_data[(cal_data$Soil=="PL" & cal_data$Status=="A"), ],
                           py_pars = py_parsPL_resp)
CT_cal_py_resp<-DB_cal_resp_py(dataset = cal_data[(cal_data$Soil=="CT" & cal_data$Status=="A"), ],
                           py_pars = py_parsCT_resp)

PL_cal_py_resp$goodness$Gfit
CT_cal_py_resp$goodness$Gfit

ggplot(PL_cal_py_resp$goodness$Yhat, aes(time, obs))+geom_point(cex=6)+geom_line(data=PL_cal_py_resp$simul, aes(time, value))+facet_wrap(~variable, scales="free")
ggplot(CT_cal_py_resp$goodness$Yhat, aes(time, obs))+geom_point(cex=6)+geom_line(data=CT_cal_py_resp$simul, aes(time, value))+facet_wrap(~variable, scales="free")
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

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

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~Python parameters~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
source("Models/DB_W_death_py.R")
#read the parameters
py_parsPLW <- as.numeric(read.csv("./DB_concept/Hasan_Jolanta/PL_parametersW.csv", header = F))
py_parsCTW <- as.numeric(read.csv("./DB_concept/Hasan_Jolanta/CT_parametersW.csv", header = F))

PL_W_py<-DB_W_death_py(dataset = cal_data[(cal_data$Soil=="PL" & cal_data$Status=="A"), ],
                           py_pars = py_parsPLW)
CT_W_py<-DB_W_death_py(dataset = cal_data[(cal_data$Soil=="CT" & cal_data$Status=="A"), ],
                           py_pars = py_parsCTW)

PL_W_py$goodness$Gfit
CT_W_py$goodness$Gfit

ggplot(PL_W_py$goodness$Yhat, aes(time, obs))+geom_point(cex=6)+geom_line(data=PL_W_py$simul, aes(time, value))+facet_wrap(~variable, scales="free")
ggplot(CT_W_py$goodness$Yhat, aes(time, obs))+geom_point(cex=6)+geom_line(data=CT_W_py$simul, aes(time, value))+facet_wrap(~variable, scales="free")
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

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

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~Python parameters~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
source("Models/DB_W_resp_py.R")
#read the parameters
py_parsPL_respW <- as.numeric(read.csv("./DB_concept/Hasan_Jolanta/PL_parameters_respW.csv", header = F))
py_parsCT_respW <- as.numeric(read.csv("./DB_concept/Hasan_Jolanta/CT_parameters_respW.csv", header = F))

PL_W_py_resp<-DB_W_resp_py(dataset = cal_data[(cal_data$Soil=="PL" & cal_data$Status=="A"), ],
                               py_pars = py_parsPL_respW)
CT_W_py_resp<-DB_W_resp_py(dataset = cal_data[(cal_data$Soil=="CT" & cal_data$Status=="A"), ],
                               py_pars = py_parsCT_respW)

PL_W_py_resp$goodness$Gfit
CT_W_py_resp$goodness$Gfit

ggplot(PL_W_py_resp$goodness$Yhat, aes(time, obs))+geom_point(cex=6)+geom_line(data=PL_W_py_resp$simul, aes(time, value))+facet_wrap(~variable, scales="free")
ggplot(CT_W_py_resp$goodness$Yhat, aes(time, obs))+geom_point(cex=6)+geom_line(data=CT_W_py_resp$simul, aes(time, value))+facet_wrap(~variable, scales="free")
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#


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
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~Python parameters~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
source("Models/DB_cal_deathD_py.R")
#read the parameters
py_parsPLD <- as.numeric(read.csv("./DB_concept/Hasan_Jolanta/PL_parametersD.csv", header = F))
py_parsCTD <- as.numeric(read.csv("./DB_concept/Hasan_Jolanta/CT_parametersD.csv", header = F))

PL_calD_py<-DB_cal_deathD_py(dataset = cal_data[(cal_data$Soil=="PL" & cal_data$Status=="A"), ],
                           py_pars = py_parsPLD)
CT_calD_py<-DB_cal_deathD_py(dataset = cal_data[(cal_data$Soil=="CT" & cal_data$Status=="A"), ],
                           py_pars = py_parsCTD)

PL_calD_py$goodness$Gfit
CT_calD_py$goodness$Gfit

ggplot(PL_calD_py$goodness$Yhat, aes(time, obs))+geom_point(cex=6)+geom_line(data=PL_calD_py$simul, aes(time, value))+facet_wrap(~variable, scales="free")
ggplot(CT_calD_py$goodness$Yhat, aes(time, obs))+geom_point(cex=6)+geom_line(data=CT_calD_py$simul, aes(time, value))+facet_wrap(~variable, scales="free")
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#Respiratory losses
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~Python parameters~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
source("Models/DB_cal_respD_py.R")
#read the parameters
py_parsPLD2 <- as.numeric(read.csv("./DB_concept/Hasan_Jolanta/PL_parameters_respD2.csv", header = F))
py_parsCTD2 <- as.numeric(read.csv("./DB_concept/Hasan_Jolanta/CT_parameters_respD2.csv", header = F))

PL_calD2_py<-DB_cal_respD_py(dataset = cal_data[(cal_data$Soil=="PL" & cal_data$Status=="A"), ],
                             py_pars = py_parsPLD2)
CT_calD2_py<-DB_cal_respD_py(dataset = cal_data[(cal_data$Soil=="CT" & cal_data$Status=="A"), ],
                             py_pars = py_parsCTD2)

PL_calD2_py$goodness$Gfit
CT_calD2_py$goodness$Gfit

ggplot(PL_calD2_py$goodness$Yhat, aes(time, obs))+geom_point(cex=6)+geom_line(data=PL_calD2_py$simul, aes(time, value))+facet_wrap(~variable, scales="free")
ggplot(CT_calD2_py$goodness$Yhat, aes(time, obs))+geom_point(cex=6)+geom_line(data=CT_calD2_py$simul, aes(time, value))+facet_wrap(~variable, scales="free")
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
##################################Possible C sorption################################
cal_data %>% filter(Status=="A") %>% group_by(Soil, time) %>% summarize(DOC = DOC12 + DOC13,
                                                                        WOC = WOC12 + WOC13)
mixing %>% group_by(Plesne, cas, znaceni, horizont) %>% summarize(DOC = mean(DOC2, na.rm=T),
                                                                 WOC = mean(DOC, na.rm = T))

source("Models/DB_cal_death_sorp_py.R")
#read the parameters
py_parsPL_sorp <- as.numeric(read.csv("./DB_concept/Hasan_Jolanta/PL_parameters_sorp.csv", header = F))
py_parsCT_sorp <- as.numeric(read.csv("./DB_concept/Hasan_Jolanta/CT_parameters_sorp.csv", header = F))

PL_cal_sorp_py<-DB_cal_death_sorp_py(dataset = cal_data[(cal_data$Soil=="PL" & cal_data$Status=="A"), ],
                             py_pars = py_parsPL_sorp)
CT_cal_sorp_py<-DB_cal_death_sorp_py(dataset = cal_data[(cal_data$Soil=="CT" & cal_data$Status=="A"), ],
                             py_pars = py_parsCT_sorp)

PL_cal_sorp_py$goodness$Gfit
CT_cal_sorp_py$goodness$Gfit

ggplot(PL_cal_sorp_py$goodness$Yhat, aes(time, obs))+geom_point(cex=6)+geom_line(data=PL_cal_sorp_py$simul, aes(time, value))+facet_wrap(~variable, scales="free")
ggplot(CT_cal_sorp_py$goodness$Yhat, aes(time, obs))+geom_point(cex=6)+geom_line(data=CT_cal_sorp_py$simul, aes(time, value))+facet_wrap(~variable, scales="free")



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

#Plotting the correspondence
ggplot(db_mixing_out$Yhat_treatments[db_mixing_out$Yhat_treatments$horizon=="Litter", ], 
       aes(value, obs))+geom_point(aes(colour=horizon))+
  facet_wrap(Treatment~variable, scales="free")+geom_abline(intercept=0, slope = 1)

#Plotting the parameters
##in long format
Pars<-melt(db_mixing_out$pars_all, id.vars = c("Plesne", "Certovo", "horizon", "id"))
##Add uncertainty
Pars$sd<-melt(cbind(db_mixing_out$pars_all[, c("Plesne", "Certovo", "horizon", "id")],
                    as.data.frame(rbind(summary(db_mixing_out$pars_raw[[1]]$mcmc_out)["sd", ],
                                        summary(db_mixing_out$pars_raw[[2]]$mcmc_out)["sd", ],
                                        summary(db_mixing_out$pars_raw[[3]]$mcmc_out)["sd", ],
                                        summary(db_mixing_out$pars_raw[[4]]$mcmc_out)["sd", ],
                                        summary(db_mixing_out$pars_raw[[5]]$mcmc_out)["sd", ],
                                        summary(db_mixing_out$pars_raw[[6]]$mcmc_out)["sd", ],
                                        summary(db_mixing_out$pars_raw[[7]]$mcmc_out)["sd", ],
                                        summary(db_mixing_out$pars_raw[[8]]$mcmc_out)["sd", ],
                                        summary(db_mixing_out$pars_raw[[9]]$mcmc_out)["sd", ],
                                        summary(db_mixing_out$pars_raw[[10]]$mcmc_out)["sd", ]))), 
      id.vars = c("Plesne", "Certovo", "horizon", "id"))[, "value"]

ggplot(Pars, aes(Plesne, value))+geom_point(cex=6, aes(colour=horizon))+
  geom_errorbar(aes(ymin=value-sd, ymax=value+sd, color=horizon))+
  facet_wrap(horizon~variable, scales="free")

#f parameter can be fixed across all treatments, but separate estimates are applied 
#for litter and organic soil 
mean(db_mixing_out$pars_all[db_mixing_out$pars_all$horizon=="Litter", "f"])
mean(db_mixing_out$pars_all[db_mixing_out$pars_all$horizon=="Organic soil", "f"])

########################################Fixing f parameter################################################
#loading the function
source("./Models/DB_mixing_f.R")

#defining number of cores
no_cors<-detectCores()
#creating cluster
cl<-makeCluster(no_cors)
#registering cluster
registerDoParallel(cl)

#function run
db_mixing_out_f<-DB_mixing_f(dataset=mixing, fractionsPL = fractionsPL, fractionsCT = fractionsCT,
                             initPL = initPL[-8], initCT = initCT[-8], minpar=minpars[-8], maxpar = maxpars[-8])

stopImplicitCluster()

#checking the goodness of correspondence
db_mixing_out_f$goodness1
db_mixing_out_f$goodness2
db_mixing_out_f$goodness3
db_mixing_out_f$goodness4

#Plotting the correspondence
ggplot(db_mixing_out_f$Yhat_treatments[db_mixing_out_f$Yhat_treatments$horizon=="Organic soil", ], 
       aes(value, obs))+geom_point(aes(colour=horizon))+
  facet_wrap(Treatment~variable, scales="free")+geom_abline(intercept=0, slope = 1)

#Plotting the parameters
##in long format
Pars_f<-melt(db_mixing_out_f$pars_all, id.vars = c("Plesne", "Certovo", "horizon", "id"))
##Add uncertainty
Pars_f$sd<-melt(cbind(db_mixing_out_f$pars_all[, c("Plesne", "Certovo", "horizon", "id")],
                    as.data.frame(rbind(summary(db_mixing_out_f$pars_raw[[1]]$mcmc_out)["sd", ],
                                        summary(db_mixing_out_f$pars_raw[[2]]$mcmc_out)["sd", ],
                                        summary(db_mixing_out_f$pars_raw[[3]]$mcmc_out)["sd", ],
                                        summary(db_mixing_out_f$pars_raw[[4]]$mcmc_out)["sd", ],
                                        summary(db_mixing_out_f$pars_raw[[5]]$mcmc_out)["sd", ],
                                        summary(db_mixing_out_f$pars_raw[[6]]$mcmc_out)["sd", ],
                                        summary(db_mixing_out_f$pars_raw[[7]]$mcmc_out)["sd", ],
                                        summary(db_mixing_out_f$pars_raw[[8]]$mcmc_out)["sd", ],
                                        summary(db_mixing_out_f$pars_raw[[9]]$mcmc_out)["sd", ],
                                        summary(db_mixing_out_f$pars_raw[[10]]$mcmc_out)["sd", ]))), 
              id.vars = c("Plesne", "Certovo", "horizon", "id"))[, "value"]

ggplot(Pars_f, aes(Plesne, value))+geom_point(cex=6, aes(colour=horizon))+
  geom_errorbar(aes(ymin=value-sd, ymax=value+sd, color=horizon))+
  facet_wrap(horizon~variable, scales="free")

#Plotting the simulations
Simul_f<-db_mixing_out_f$simul

Simul_f %>% filter(horizon=="Organic soil" & variable=="Cmic_12C") %>% 
  group_by(time, variable, Plesne, Certovo, Treatment) %>% 
  summarize(value = mean(value)) %>%
  ggplot(aes(time, value))+geom_line(aes(colour = Treatment))+
  facet_wrap(variable~Plesne, scales="free")

######################################Organic soil only#######################################
#Use model parameters from calibration phase - same soils
#Difference is theoretically given by the initial amount of Reserves
names(py_parsPLRS)<-c("Ac_glucose", "Vmaxg", "Kmg", 
                    "Ac_DOC", "Vmax", "Km",
                    "mr", "f", "Yu", "fs", "fr", "RSinit")
names(py_parsCTRS)<-c("Ac_glucose", "Vmaxg", "Kmg", 
                    "Ac_DOC", "Vmax", "Km",
                    "mr", "f", "Yu", "fs", "fr", "RSinit")

#loading the function
source("./Models/DB_mixing_organic.R")

#defining number of cores
no_cors<-detectCores()-2
#creating cluster
cl<-makeCluster(no_cors)
#registering cluster
registerDoParallel(cl)

#function run
db_mixing_organic<-DB_mixing_organic(dataset=mixing, initPL = py_parsPLRS, initCT = py_parsCTRS)

stopImplicitCluster()

#checking the goodness of correspondence
db_mixing_organic$goodness1
db_mixing_organic$goodness2
db_mixing_organic$goodness3
db_mixing_organic$goodness4