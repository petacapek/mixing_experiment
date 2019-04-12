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
source("Models/DB_cal.R")
PL_cal<-DB_cal(dataset = cal_data[(cal_data$Soil=="PL" & cal_data$Status=="A"), ])
CT_cal<-DB_cal(dataset = cal_data[(cal_data$Soil=="CT" & cal_data$Status=="A"), ])
