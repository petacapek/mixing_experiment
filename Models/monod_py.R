Monod_py<-function(dataset, py_pars){
  
  #model is defined here
  monod_model<-function(time, state, pars){
    with(as.list(c(state, pars)),{
      #Abbreviations:
      #States:
        #MBC_12C - 12C in microbial biomass
        #MBC_13C - 13C in microbial biomass
        #Cmic_12C - 12C chloroform labile organic carbon from microbial biomass 
        #Cmic_13C - 13C chloroform labile organic carbon from microbial biomass
        #G_12C - 12C glucose
        #G_13C - 13C glucose
        #DOC_12C - 12C dissolved organic carbon
        #DOC_13C - 13C dissolved organic carbon
        #Cres_12C - 12C in residual pool - unsoluble part of death cells
        #Cres_13C - 13C in residual pool - unsoluble part of death cells
        #CO2_12C - cumulative 12C-CO2 in headspace
        #CO2_13C - cumulative 13C-CO2 in headspace
      #Fluxes:
        #Cu_glucose - uptake rate of the glucose
        #Cu_DOC - uptake rate of the DOC
        #r_12C - 12C-CO2 production rate
        #r_13C - 13C-CO2 production rate
      #Scaling factors:
        #Gatm - 13C atm% of glucose
        #DOCatm - 13C atm% of dissolved organic carbon
      #Model parameters:
        #Vmax_glucose - maximum velocity constant for glucose uptake rate
        #Vmax_DOC - maximum velocity constant for dissolved organic carbon uptake rate
        #Km_glucose - affinity constant of glucose uptake
        #Km_DOC - affinity constant of dissolved organic carbon uptake
        #CUE - carbon use efficiency
        #kb - microbial biomass decay rate
        #kec - chloroform labile part of microbial biomass carbon
      
      #Equations:
      #Uptake rate - we assume that glucose is preferred substrate for uptake/growth 
      #Organic carbon uptake rate - glucose
      Cu_glucose=Vmaxg*(MBC_12C+MBC_13C)*(G_12C+G_13C)/(G_12C+G_13C+Kmg*(1+(DOC_12C+DOC_13C)/Km))#
      #Organic carbon uptake rate - DOC
      Cu_DOC=Vmax*(MBC_12C+MBC_13C)*(DOC_12C+DOC_13C)/(DOC_12C+DOC_13C+Km*(1+(G_12C+G_13C)/Kmg))#
      
      #Scaling factors:
      Gatm=G_13C/(G_12C+G_13C)
      DOCatm=DOC_13C/(DOC_12C+DOC_13C)
      
      #Chloroform labile organic carbon is part of MBC
      Cmic_12C=kec*MBC_12C
      Cmic_13C=kec*MBC_13C
      
      #States
      dMBC_12C<-CUE*Cu_glucose*(1-Gatm)+CUE*Cu_DOC*(1-DOCatm)-kb*MBC_12C
      dMBC_13C<-CUE*Cu_glucose*Gatm+CUE*Cu_DOC*DOCatm-kb*MBC_13C
      dG_12C<--Cu_glucose*(1-Gatm)
      dG_13C<--Cu_glucose*Gatm
      dDOC_12C<--Cu_DOC*(1-DOCatm)+kb*kec*MBC_12C
      dDOC_13C<--Cu_DOC*DOCatm+kb*kec*MBC_13C
      dCres_12C<-kb*(1-kec)*MBC_12C
      dCres_13C<-kb*(1-kec)*MBC_13C
      dCO2_12C<-(1-CUE)*Cu_glucose*(1-Gatm)+(1-CUE)*Cu_DOC*(1-DOCatm)
      dCO2_13C<-(1-CUE)*Cu_glucose*Gatm+(1-CUE)*Cu_DOC*DOCatm
      
      return(list(c(dMBC_12C, dMBC_13C, 
                    dG_12C, dG_13C, 
                    dDOC_12C, dDOC_13C, 
                    dCres_12C, dCres_13C,
                    dCO2_12C, dCO2_13C), Cmic_12C=Cmic_12C, Cmic_13C=Cmic_13C))
      
    })
  }
  #################################################################################################
  #Goodness of model simulation is calculated here for each variable separately
  #To do so, "good" function is defined. This function pretty much similar to "cost" function
  good<-function(x){
    
    par<-x
    names(par)<-c("Vmaxg", "Kmg", 
                  "Vmax", "Km",
                  "CUE", "kb", "kec")
    
    #Extracting initial concentration of state variables from data
    MBC_12Cinit=as.numeric(dataset[1, "Cmic12init"])/par[["kec"]]
    MBC_13Cinit=as.numeric(dataset[1, "Cmic13init"])/par[["kec"]]
    G_12Cinit=as.numeric(dataset[1, "G12init"])
    G_13Cinit=as.numeric(dataset[1, "G13init"])
    DOC_12Cinit=as.numeric(dataset[1, "DOC12init"])
    DOC_13Cinit=as.numeric(dataset[1, "DOC13init"])
    #Cres and CO2 pools are initialy 0
    
    #time of the sampling
    t_sampling=as.numeric(dataset$time)
    
    #model simulation is run here
    yhat_all<-as.data.frame(ode(y=c(MBC_12C=MBC_12Cinit, MBC_13C=MBC_13Cinit,
                                    G_12C=G_12Cinit, G_13C=G_13Cinit,
                                    DOC_12C=DOC_12Cinit, DOC_13C=DOC_13Cinit,
                                    Cres_12C=0, Cres_13C=0, CO2_12C=0, CO2_13C=0), 
                                parms=par, monod_model, times=t_sampling))
    
    #variables that were measured in the experiment are extracted
    yhat<-select(yhat_all, c("time", "G_12C", "G_13C", "DOC_12C", "DOC_13C", "CO2_12C", "CO2_13C", "Cmic_12C", "Cmic_13C"))
    
    #convert the simulated dataset into long format data frame
    Yhat<-melt(yhat, id.vars = "time")
    
    #add the measured data
    Yhat$obs<-c(as.numeric(dataset[,"G12"]), as.numeric(dataset[,"G13"]),
                as.numeric(dataset[,"DOC12"]), as.numeric(dataset[,"DOC13"]),
                as.numeric(dataset[,"CO212"]), as.numeric(dataset[,"CO213"]),
                as.numeric(dataset[,"Cmic12"]), as.numeric(dataset[,"Cmic13"]))
    
    #several goodness of fit metrics are calulated for each variable
    #number of data points for each variable
    n=length(t_sampling)
    
    Gfit<-Yhat %>% group_by(variable) %>% summarise(SSres=sum(((obs-value)^2), na.rm = T), 
                                                    SStot=sum(((obs-mean(obs, na.rm = T))^2), na.rm = T),
                                                    ll=-sum(((obs-value)^2), na.rm = T)/2/(sd(obs, na.rm = T)^2))
    Gfit$R2<-with(Gfit, 1-SSres/SStot)
    Gfit$N<-length(x)
    Gfit$AIC<-with(Gfit, 2*N-2*ll)
    
    rsq_out<-list(Yhat=Yhat, Gfit=Gfit)
    
    return(rsq_out)
  }
  
  #The goodness of correspondence between model simulation and measured data is calculated 
  #using "good" function here
  fit<-good(py_pars)
  
  #For nice plots, model simulation is run with calibrated model parameters
  #at fine temporal scale and the simulation is stored in "simul" data frame
  #First, initial concentration of state variables are defined
  py_pars<-py_pars
  names(py_pars)<-c("Vmaxg", "Kmg", 
                    "Vmax", "Km",
                    "CUE", "kb", "kec")
  
  MBC_12Cinit=as.numeric(dataset[1, "Cmic12init"])/py_pars[["kec"]]
  MBC_13Cinit=as.numeric(dataset[1, "Cmic13init"])/py_pars[["kec"]]
  G_12Cinit=as.numeric(dataset[1, "G12init"])
  G_13Cinit=as.numeric(dataset[1, "G13init"])
  DOC_12Cinit=as.numeric(dataset[1, "DOC12init"])
  DOC_13Cinit=as.numeric(dataset[1, "DOC13init"])
  #Cres and CO2 pools are initialy 0
  
  #time of the sampling
  t_simul=seq(from=min(dataset$time), to=max(dataset$time), length.out = 150)
  
  #model simulation is run here
  simul<-as.data.frame(ode(y=c(MBC_12C=MBC_12Cinit, MBC_13C=MBC_13Cinit,
                               G_12C=G_12Cinit, G_13C=G_13Cinit,
                               DOC_12C=DOC_12Cinit, DOC_13C=DOC_13Cinit,
                               Cres_12C=0, Cres_13C=0, CO2_12C=0, CO2_13C=0),
                           parms=py_pars, monod_model, times=t_simul))
  Simul<-melt(simul, id.vars = "time")
  
  f_out<-list(goodness=fit,
              simul=Simul)
  
  return(f_out)
}