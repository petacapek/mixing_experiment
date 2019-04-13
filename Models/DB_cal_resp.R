DB_cal_resp<-function(dataset){
  
  #model is defined here
  db_model<-function(time, state, pars){
    with(as.list(c(state, pars)),{
      #Abbreviations:
      #States:
        #R_12C - 12C in Reserves
        #R_13C - 13C in Reserves
        #S_12C - 12C in Structures 
        #S_13C - 13C in Structures 
        #Cmic_12C - 12C chloroform labile organic carbon from microbial biomass 
        #Cmic_13C - 13C chloroform labile organic carbon from microbial biomass
        #G_12C - 12C glucose
        #G_13C - 13C glucose
        #DOC_12C - 12C dissolved organic carbon
        #DOC_13C - 13C dissolved organic carbon
        #CO2_12C - cumulative 12C-CO2 in headspace
        #CO2_13C - cumulative 13C-CO2 in headspace
      #Fluxes:
        #Cu_glucose - uptake rate of the glucose
        #Cu_DOC - uptake rate of the DOC
        #m - maintnance rate of Structures 
        #an - mobilization rate of Reserves available for growth
        #r_12C - 12C-CO2 production rate
        #r_13C - 13C-CO2 production rate
      #Scaling factors:
        #Ratm - 13C atm% of Reserves pool
        #Satm - 13C atm% of Structures pool
        #Gatm - 13C atm% of glucose
        #DOCatm - 13C atm% of dissolved organic carbon
      #Model parameters:
        #Vmax_glucose - maximum velocity constant for glucose uptake rate
        #Vmax_DOC - maximum velocity constant for dissolved organic carbon uptake rate
        #Km_glucose - affinity constant of glucose uptake
        #Km_DOC - affinity constant of dissolved organic carbon uptake
        #Ac_glucose - assimilation efficiency of glucose
        #Ac_DOC - assimilation efficiency of dissolved organic carbon
        #mr - maintnance rate constant
        #f - reserves mobilization rate constant
        #Yu - growth yield 
        #fr - chloroform labile part of Reserves
        #fs - chloroform labile part of Structures
      
      #Equations:
      #Uptake rate - we assume that glucose is preferred substrate for uptake/growth 
      #Organic carbon uptake rate - glucose
      Cu_glucose=Vmax_glucose*(S_12C+S_13C)*(G_12C+G_13C)/(Km_glucose+(G_12C+G_13C)*(1+(DOC_12C+DOC_13C)/Km_DOC))
      #Organic carbon uptake rate - DOC
      Cu_DOC=Vmax_DOC*(S_12C+S_13C)*(DOC_12C+DOC_13C)/(Km_DOC+(DOC_12C+DOC_13C)*(1+(G_12C+G_13C)/Km_glucose))
      
      #maintnance of Structures
      m=mr*(S_12C+S_13C)
      
      #mobilization rate of Reserves available for growth
      an=f*(R_12C+R_13C)-m
      
      #Scaling factors:
      Ratm=R_13C/(R_12C+R_13C)
      Satm=S_13C/(S_12C+S_13C)
      Gatm=G_13C/(G_12C+G_13C)
      DOCatm=DOC_13C/(DOC_12C+DOC_13C)
      
      #Respiration rate
      #Respiration rate is composed of three processes - organic carbon assimilation associated respiration,
      #growth associated respiration, and maintnance respiration.
      #Growth associated respiration occur only when growth is realized.
      #Maintnance respiration is defined by the concentration of Structures. When Reserves do not contain enough carbon to cover 
      #Structures maintnance, required difference is cover by respiration from Structures.
      r_12C<-(1-Ac_glucose)*Cu_glucose*(1-Gatm)+(1-Ac_DOC)*Cu_DOC*(1-DOCatm)+pmax((1-Yu)*an*(1-Ratm), 0)+ifelse(an>0, m*(1-Ratm), f*R_12C-an*(1-Satm))
      r_13C<-(1-Ac_glucose)*Cu_glucose*Gatm+(1-Ac_DOC)*Cu_DOC*DOCatm+pmax((1-Yu)*an*Ratm, 0)+ifelse(an>0, m*Ratm, f*R_13C-an*Satm)
      
      #Chloroform labile organic carbon is part of Reserves and part of Structures
      Cmic_12C=fr*R_12C+fs*S_12C
      Cmic_13C=fr*R_13C+fs*S_13C
      
      #States
      #If Reserves contain enough C, growth (i.e. Structures production) can be realized.
      #In opposite case, Structures cannot increase.
      #If C in Reserves is insufficient to cover maintnance of Structures (i.e. an is negative),
      #respective amount of Structures are lost to DOC and Cres pool.
      #The partitioning of C lost from Structures between DOC and Cres pool is controlled by the fs parameter.
      dR_12C<-Ac_glucose*Cu_glucose*(1-Gatm)+Ac_DOC*Cu_DOC*(1-DOCatm)-f*R_12C
      dR_13C<-Ac_glucose*Cu_glucose*Gatm+Ac_DOC*Cu_DOC*DOCatm-f*R_13C
      dS_12C<-pmax(an*Yu*(1-Ratm), 0)+pmin(0, an*(1-Satm))
      dS_13C<-pmax(an*Yu*Ratm, 0)+pmin(0, an*Satm)
      dG_12C<--Cu_glucose*(1-Gatm)
      dG_13C<--Cu_glucose*Gatm
      dDOC_12C<--Cu_DOC*(1-DOCatm)
      dDOC_13C<--Cu_DOC*DOCatm
      dCO2_12C<-r_12C
      dCO2_13C<-r_13C
      
      return(list(c(dR_12C, dR_13C, 
                    dS_12C, dS_13C, 
                    dG_12C, dG_13C, 
                    dDOC_12C, dDOC_13C, 
                    dCO2_12C, dCO2_13C), Cmic_12C=Cmic_12C, Cmic_13C=Cmic_13C))
      
    })
  }
  #define names of parameters
  parnames<-c("Vmax_glucose", "Vmax_DOC", 
              "Km_glucose", "Km_DOC", 
              "Ac_glucose", "Ac_DOC",
              "mr", "f", "Yu", "fr", "fs", "Rinit", "Ratm_init")
  
  #Minimization ("cost") function is defined here
  cost<-function(x){
    
    par<-x[1:length(parnames)]
    names(par)<-parnames
    
    #Extracting initial concentration of state variables from data
    #The initial abundance of Reserves and Strctures in microbial biomass is not know.
    #Therefore initial concentration of Reserves ("Rinit") is estimated together with all model parameters.
    #It is also not known if the Reserves and Structures have the same isotopic signal (i.e. Ratm/Satm).
    #Therefore, one more parameter - Ratm_init, has to be defined and estimated.
    R_12Cinit=par[["Rinit"]]*(1-par[["Ratm_init"]])
    R_13Cinit=par[["Rinit"]]*par[["Ratm_init"]]
    S_12Cinit=(as.numeric(dataset[1, "Cmic12init"])-par[["fr"]]*R_12Cinit)/par[["fs"]]
    S_13Cinit=(as.numeric(dataset[1, "Cmic13init"])-par[["fr"]]*R_13Cinit)/par[["fs"]]
    G_12Cinit=as.numeric(dataset[1, "G12init"])
    G_13Cinit=as.numeric(dataset[1, "G13init"])
    DOC_12Cinit=as.numeric(dataset[1, "DOC12init"])
    DOC_13Cinit=as.numeric(dataset[1, "DOC13init"])
    #Cres and CO2 pools are initialy 0
    
    #time of the sampling
    t_sampling=as.numeric(dataset$time)
    
    #model simulation is run here
    yhat_all<-as.data.frame(ode(y=c(R_12C=R_12Cinit, R_13C=R_13Cinit,
                                    S_12C=S_12Cinit, S_13C=S_13Cinit,
                                    G_12C=G_12Cinit, G_13C=G_13Cinit,
                                    DOC_12C=DOC_12Cinit, DOC_13C=DOC_13Cinit,
                                    CO2_12C=0, CO2_13C=0), 
                                parms=par, db_model, times=t_sampling))
    
    #variables that were measured in the experiment are extracted
    yhat<-select(yhat_all, c("time", "G_12C", "G_13C", "DOC_12C", "DOC_13C", "CO2_12C", "CO2_13C", "Cmic_12C", "Cmic_13C"))
    
    #convert the simulated dataset into long format data frame
    Yhat<-melt(yhat, id.vars = "time")
    
    #add the measured data
    Yhat$obs<-c(as.numeric(dataset[,"G12"]), as.numeric(dataset[,"G13"]),
                as.numeric(dataset[,"DOC12"]), as.numeric(dataset[,"DOC13"]),
                as.numeric(dataset[,"CO212"]), as.numeric(dataset[,"CO213"]),
                as.numeric(dataset[,"Cmic12"]), as.numeric(dataset[,"Cmic13"]))
    
    #add the weighting factor
    #I want to have the weighting factor to be proportional to mean of the particular variables 
    yweights<-Yhat %>% group_by(variable) %>% summarise(weights=mean(value, na.rm = T))
    
    #match yweights with the Yhat data frame
    Yhat<-merge(Yhat, yweights, by.x=c("variable"), by.y=c("variable"))
    
    #now, the root mean square error is calculated
    RMSE<-as.numeric(Yhat %>% summarise(RMSE=sum((((value-obs)/mean(weights))^2), na.rm = T)))
    
    return(RMSE)
  }
  
  #Model parameters are calibrated here
  #The calibration is done in two consecutive steps
  #First, approximate model parameters are estimated using MCMC algorithm
  #approximate parameter estimation is done by MCMC method
  Rinit_guess=as.numeric(dataset[1, "Cmic12init"])
  Ratm_init_guess=as.numeric(dataset[1, "Cmic13init"])/(as.numeric(dataset[1, "Cmic12init"])+as.numeric(dataset[1, "Cmic13init"]))
  
  par_mcmc<-modMCMC(f=cost, 
                    p=c(Vmax_glucose=0.1, Vmax_DOC=0.05, 
                        Km_glucose=10, Km_DOC=10, 
                        Ac_glucose=0.9, Ac_DOC=0.8,
                        mr=1e-4, f=1e-3, Yu=0.6, fr=0.8, fs=0.3, 
                        Rinit=0.1*Rinit_guess, Ratm_init=Ratm_init_guess), 
                    lower=c(Vmax_glucose=1e-4, Vmax_DOC=1e-4, 
                            Km_glucose=1e-4, Km_DOC=1e-4, 
                            Ac_glucose=0, Ac_DOC=0,
                            mr=1e-8, f=1e-8, Yu=0, fr=0, fs=0, 
                            Rinit=1e-3*Rinit_guess, Ratm_init=0.5*Ratm_init_guess),
                    upper=c(Vmax_glucose=1e2, Vmax_DOC=1e2, 
                            Km_glucose=1e3, Km_DOC=1e3, 
                            Ac_glucose=1, Ac_DOC=1,
                            mr=1e2, f=1e2, Yu=1, fr=1, fs=1, 
                            Rinit=0.95*Rinit_guess, Ratm_init=1.2*Ratm_init_guess), niter=10000)
  
  #lower and upper limits for parameters estimates are extracted
  pl<-summary(par_mcmc)["min",]
  pu<-summary(par_mcmc)["max",]
  
  #these limits are used to find global optimum by DEoptim algorithm
  opt_par<-DEoptim(fn=cost, lower=pl, upper=pu, 
                   control = c(itermax = 10000, steptol = 50, reltol = 1e-8, 
                               trace=FALSE, strategy=3, NP=250))
  
  #################################################################################################
  #Goodness of model simulation is calculated here for each variable separately
  #To do so, "good" function is defined. This function pretty much similar to "cost" function
  good<-function(x){
    
    par<-x[1:length(parnames)]
    names(par)<-parnames
    
    #Extracting initial concentration of state variables from data
    #The initial abundance of Reserves and Strctures in microbial biomass is not know.
    #Therefore initial concentration of Reserves ("Rinit") is estimated together with all model parameters.
    #It is also not known if the Reserves and Structures have the same isotopic signal (i.e. Ratm/Satm).
    #Therefore, one more parameter - Ratm_init, has to be defined and estimated.
    R_12Cinit=par[["Rinit"]]*(1-par[["Ratm_init"]])
    R_13Cinit=par[["Rinit"]]*par[["Ratm_init"]]
    S_12Cinit=(as.numeric(dataset[1, "Cmic12init"])-par[["fr"]]*R_12Cinit)/par[["fs"]]
    S_13Cinit=(as.numeric(dataset[1, "Cmic13init"])-par[["fr"]]*R_13Cinit)/par[["fs"]]
    G_12Cinit=as.numeric(dataset[1, "G12init"])
    G_13Cinit=as.numeric(dataset[1, "G13init"])
    DOC_12Cinit=as.numeric(dataset[1, "DOC12init"])
    DOC_13Cinit=as.numeric(dataset[1, "DOC13init"])
    #Cres and CO2 pools are initialy 0
    
    #time of the sampling
    t_sampling=as.numeric(dataset$time)
    
    #model simulation is run here
    yhat_all<-as.data.frame(ode(y=c(R_12C=R_12Cinit, R_13C=R_13Cinit,
                                    S_12C=S_12Cinit, S_13C=S_13Cinit,
                                    G_12C=G_12Cinit, G_13C=G_13Cinit,
                                    DOC_12C=DOC_12Cinit, DOC_13C=DOC_13Cinit,
                                    CO2_12C=0, CO2_13C=0), 
                                parms=par, db_model, times=t_sampling))
    
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
                                                    ll=-n*log(2*mean(obs, na.rm = T))/2-n*log(sd(obs, na.rm = T)^2)/2-sum(((obs-value)^2), na.rm = T)/2/(sd(obs, na.rm = T)^2))
    Gfit$R2<-with(Gfit, 1-SSres/SStot)
    Gfit$N<-length(x)
    Gfit$AIC<-with(Gfit, 2*N-2*ll)
    
    rsq_out<-list(Yhat=Yhat, Gfit=Gfit)
    
    return(rsq_out)
  }
  
  #The goodness of correspondence between model simulation and measured data is calculated 
  #using "good" function here
  fit<-good(opt_par$optim$bestmem)
  
  #For nice plots, model simulation is run with calibrated model parameters
  #at fine temporal scale and the simulation is stored in "simul" data frame
  #First, initial concentration of state variables are defined
  R_12Cinit=opt_par$optim$bestmem[["Rinit"]]*(1-opt_par$optim$bestmem[["Ratm_init"]])
  R_13Cinit=opt_par$optim$bestmem[["Rinit"]]*opt_par$optim$bestmem[["Ratm_init"]]
  S_12Cinit=(as.numeric(dataset[1, "Cmic12init"])-opt_par$optim$bestmem[["fr"]]*R_12Cinit)/opt_par$optim$bestmem[["fs"]]
  S_13Cinit=(as.numeric(dataset[1, "Cmic13init"])-opt_par$optim$bestmem[["fr"]]*R_13Cinit)/opt_par$optim$bestmem[["fs"]]
  G_12Cinit=as.numeric(dataset[1, "G12init"])
  G_13Cinit=as.numeric(dataset[1, "G13init"])
  DOC_12Cinit=as.numeric(dataset[1, "DOC12init"])
  DOC_13Cinit=as.numeric(dataset[1, "DOC13init"])
  #Cres and CO2 pools are initialy 0
  
  #time of the sampling
  t_simul=seq(from=min(dataset$time), to=max(dataset$time), length.out = 150)
  
  #model simulation is run here
  simul<-as.data.frame(ode(y=c(R_12C=R_12Cinit, R_13C=R_13Cinit,
                               S_12C=S_12Cinit, S_13C=S_13Cinit,
                               G_12C=G_12Cinit, G_13C=G_13Cinit,
                               DOC_12C=DOC_12Cinit, DOC_13C=DOC_13Cinit,
                               CO2_12C=0, CO2_13C=0),
                           parms=opt_par$optim$bestmem, db_model, times=t_simul))
  Simul<-melt(simul, id.vars = "time")
  
  #All important calulations are stored in the "f_out" list and returned
  #1. best model parameter
  #2. goodness of correspondence
  #3. MCMC output
  #4. simulation
  
  f_out<-list(pars=opt_par$optim$bestmem,
              goodness=fit,
              par_mcmc=par_mcmc,
              simul=Simul)
  
  return(f_out)
}