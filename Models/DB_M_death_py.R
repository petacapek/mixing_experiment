DB_M_death_py<-function(dataset, py_pars){
  
  #model is defined here
  db_model<-function(time, state, pars){
    with(as.list(c(state, pars)),{
      
      #Equations:
      #Uptake rate - we assume that glucose is preferred substrate for uptake/growth 
      #Organic carbon uptake rate - glucose
      Cu_glucose=Vmaxg*S*G/(G+Kmg*(1+DOC/Km))#
      #Organic carbon uptake rate - DOC
      Cu_DOC=Vmax*S*DOC/(DOC+Km*(1+G/Kmg))#
      
      #maintnance of Structures
      m=mr*S
      
      #mobilization rate of Reserves available for growth
      an=f*R-m
      
      #Respiration rate
      #Respiration rate is composed of three processes - organic carbon assimilation associated respiration,
      #growth associated respiration, and maintnance respiration.
      #Growth associated respiration occur only when growth is realized.
      #Maintnance respiration is defined by the concentration of Structures, but at the same time, mobilization rate of Reserves.
      #When an is negative, all C mobilized from pool of Reserves are respired, which doesn't need to correspond with 
      #the maintnance requirements. In that case, maintnance respiration is lower than it should be and Structers are dying.
      r<-(1-Ac_glucose)*Cu_glucose+(1-Ac_DOC)*Cu_DOC+pmax((1-Yu)*an, 0)+ifelse(an>0, m, f*R)
      
      #Chloroform labile organic carbon is part of Reserves and part of Structures
      #DNA is part of structures
      Cmic=fr*R+fs*S
      DNA=fd*S
      
      #States
      #If Reserves contain enough C, growth (i.e. Structures production) can be realized.
      #In opposite case, Structures cannot increase.
      #If C in Reserves is insufficient to cover maintnance of Structures (i.e. an is negative),
      #respective amount of Structures are lost to DOC and Cres pool.
      #The partitioning of C lost from Structures between DOC and Cres pool is controlled by the fs parameter.
      dR<-Ac_glucose*Cu_glucose+Ac_DOC*Cu_DOC-f*R
      dS<-pmax(an*Yu, 0)+pmin(0, an/mr)
      dG<--Cu_glucose
      dDOC<--Cu_DOC-pmin(0, an/mr*fs)
      dCres<--pmin(0, an/mr*(1-fs))
      dCO2<-r
      
      return(list(c(dR, dS, dG, dDOC, dCres, dCO2), 
                  Cmic=Cmic, DNA=DNA, r=r))
      
    })
  }
  #################################################################################################
  #Goodness of model simulation is calculated here for each variable separately
  #To do so, "good" function is defined. This function pretty much similar to "cost" function
  good<-function(x){
    
    par<-x
    names(par)<-c("Ac_glucose", "Vmaxg", "Kmg", 
                  "Ac_DOC", "Vmax", "Km",
                  "mr", "f", "Yu", "fs", "fr", "fd", "DOCinit")
    
    #Extracting initial concentration of state variables from data
    #The initial abundance of Reserves and Strctures in microbial biomass is not know.
    #Therefore initial concentration of Reserves ("Rinit") is estimated together with all model parameters.
    #It is also not known if the Reserves and Structures have the same isotopic signal (i.e. Ratm/Satm).
    #Therefore, one more parameter - Ratm_init, has to be defined and estimated.
    Sinit=(as.numeric(dataset[1, "DNAinit"]))/par[["fd"]]
    Ginit=as.numeric(dataset[1, "Glinit"])
    DOCinit=par[["DOCinit"]]
    CLCinit=as.numeric(dataset[1, "CLCinit"])
    Rinit=(CLCinit-Sinit*par[["fs"]])/par[["fr"]]
    #Cres and CO2 pools are initialy 0
    
    #time of the sampling
    t_sampling=as.numeric(dataset$Time)
    
    #model simulation is run here
    yhat_all<-as.data.frame(ode(y=c(R=Rinit, S=Sinit, G=Ginit, DOC=DOCinit, Cres=0, CO2=0), 
                                parms=par[1:12], db_model, times=t_sampling))
    
    #variables that were measured in the experiment are extracted
    yhat<-select(yhat_all, c("time", "G", "r", "Cmic", "DNA"))
    
    #convert the simulated dataset into long format data frame
    Yhat<-melt(yhat, id.vars = "time")
    
    #add the measured data
    Yhat$obs<-c(as.numeric(dataset[,"Gl"]), as.numeric(dataset[,"CO2"]),
                as.numeric(dataset[,"CLC"]), as.numeric(dataset[,"DNA"]))
    
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
  names(py_pars)<-c("Ac_glucose", "Vmaxg", "Kmg", 
                    "Ac_DOC", "Vmax", "Km",
                    "mr", "f", "Yu", "fs", "fr", "fd", "DOCinit")
  
  Sinit=(as.numeric(dataset[1, "DNAinit"]))/py_pars[["fd"]]
  Ginit=as.numeric(dataset[1, "Glinit"])
  DOCinit=py_pars[["DOCinit"]]
  CLCinit=as.numeric(dataset[1, "CLCinit"])
  Rinit=(CLCinit-Sinit*py_pars[["fs"]])/py_pars[["fr"]]
  #Cres and CO2 pools are initialy 0
  
  #time of the sampling
  t_simul=seq(from=min(dataset$Time), to=max(dataset$Time), length.out = 150)
  
  #model simulation is run here
  simul<-as.data.frame(ode(y=c(R=Rinit, S=Sinit, G=Ginit, DOC=DOCinit, Cres=0, CO2=0), 
                           parms=py_pars[1:12], db_model, times=t_simul))
    
  Simul<-melt(simul, id.vars = "time")
  
  f_out<-list(goodness=fit,
              simul=Simul)
  
  return(f_out)
}