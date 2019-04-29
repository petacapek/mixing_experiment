DB_mixing_litter<-function(dataset, fractionsPL, fractionsCT, initPL, initCT, minpar, maxpar){
  #Some information that is needed for model simulations are pasted into the dataframe
  #1. Identification number of treatments
  dataset$id<-rep(c(rep(seq(1,5), times=4),#O horizons - unlabelled
                rep(seq(6,10), times=4),#A horizons - unlabelled
                rep(seq(11,15), times=4),#O horizons - labelled
                rep(seq(16,20), times=4)), times=2)#A horizons - labelled
  #2. Estimated fractions of Reserves and Structures in chloroform labile organic carbon
    #that are provided in "fractionsPL" and "fractionsCT" vectors entering the function
  dataset$fr<-dataset$Plesne*fractionsPL[["fr"]]+dataset$Certovo*fractionsCT[["fr"]]
  dataset$fs<-dataset$Plesne*fractionsPL[["fs"]]+dataset$Certovo*fractionsCT[["fs"]]
  #3. Rename "horizont"
  dataset$horizon<-ifelse(dataset$horizont=="O", "Litter", "Organic soil")

  #Two models are defined. The only difference between these models is presence/absence
  #of the glucose pool that is available for microbial consumption.
  #db_modelun is used to simulate tretaments with no glucose addition and
  #db_modell is used to simulate tretaments with glucose addition
  #1. Without glucose addition
  db_modelun<-function(time, state, pars){
    with(as.list(c(state, pars)),{
      #Abbreviations:
      #States:
      #R_12C - 12C in Reserves
      #R_13C - 13C in Reserves
      #S_12C - 12C in Structures
      #S_13C - 13C in Structures
      #Cmic_12C - 12C chloroform labile organic carbon from microbial biomass
      #Cmic_13C - 13C chloroform labile organic carbon from microbial biomass
      #DOC_12C - 12C dissolved organic carbon
      #DOC_13C - 13C dissolved organic carbon
      #Cres_12C - 12C in residual pool - unsoluble part of death cells
      #Cres_13C - 13C in residual pool - unsoluble part of death cells
      #CO2_12C - cumulative 12C-CO2 in headspace
      #CO2_13C - cumulative 13C-CO2 in headspace
      #Fluxes:
      #Cu_DOC - uptake rate of the DOC
      #m - maintnance rate of Structures
      #an - mobilization rate of Reserves available for growth
      #r_12C - 12C-CO2 production rate
      #r_13C - 13C-CO2 production rate
      #Scaling factors:
      #Ratm - 13C atm% of Reserves pool
      #Satm - 13C atm% of Structures pool
      #DOCatm - 13C atm% of dissolved organic carbon
      #Model parameters:
      #Vmax_DOC - maximum velocity constant for dissolved organic carbon uptake rate
      #Km_DOC - affinity constant of dissolved organic carbon uptake
      #Ac_DOC - assimilation efficiency of dissolved organic carbon
      #mr - maintnance rate constant
      #f - reserves mobilization rate constant
      #Yu - growth yield
      #fr - chloroform labile part of Reserves
      #fs - chloroform labile part of Structures

      #Equations:
      #Uptake rate
      #Organic carbon uptake rate - DOC
      Cu_DOC=Vmax_DOC*(S_12C+S_13C)*(DOC_12C+DOC_13C)/(Km_DOC+(DOC_12C+DOC_13C))#

      #maintnance of Structures
      m=mr*(S_12C+S_13C)

      #mobilization rate of Reserves available for growth
      an=f*(R_12C+R_13C)-m

      #Scaling factors:
      Ratm=R_13C/(R_12C+R_13C)
      Satm=S_13C/(S_12C+S_13C)
      DOCatm=DOC_13C/(DOC_12C+DOC_13C)

      #Respiration rate
      #Respiration rate is composed of three processes - organic carbon assimilation associated respiration,
      #growth associated respiration, and maintnance respiration.
      #Growth associated respiration occur only when growth is realized.
      #Maintnance respiration is defined by the concentration of Structures, but at the same time, mobilization rate of Reserves.
      #When an is negative, all C mobilized from pool of Reserves are respired, which doesn't need to correspond with
      #the maintnance requirements. In that case, maintnance respiration is lower than it should be and Structers are dying.
      r_12C<-(1-Ac_DOC)*Cu_DOC*(1-DOCatm)+pmax((1-Yu)*an*(1-Ratm), 0)+ifelse(an>0, m*(1-Ratm), f*R_12C)
      r_13C<-(1-Ac_DOC)*Cu_DOC*DOCatm+pmax((1-Yu)*an*Ratm, 0)+ifelse(an>0, m*Ratm, f*R_13C)

      #Chloroform labile organic carbon is part of Reserves and part of Structures
      Cmic_12C=fr*R_12C+fs*S_12C
      Cmic_13C=fr*R_13C+fs*S_13C

      #States
      #If Reserves contain enough C, growth (i.e. Structures production) can be realized.
      #In opposite case, Structures cannot increase.
      #If C in Reserves is insufficient to cover maintnance of Structures (i.e. an is negative),
      #respective amount of Structures are lost to DOC and Cres pool.
      #The partitioning of C lost from Structures between DOC and Cres pool is controlled by the fs parameter.
      dR_12C<-Ac_DOC*Cu_DOC*(1-DOCatm)-f*R_12C
      dR_13C<-Ac_DOC*Cu_DOC*DOCatm-f*R_13C
      dS_12C<-pmax(an*Yu*(1-Ratm), 0)+pmin(0, an/mr*(1-Satm))
      dS_13C<-pmax(an*Yu*Ratm, 0)+pmin(0, an/mr*Satm)
      dDOC_12C<--Cu_DOC*(1-DOCatm)-pmin(0, an/mr*(1-Satm)*fs)
      dDOC_13C<--Cu_DOC*DOCatm-pmin(0, an/mr*Satm*fs)
      dCres_12C<--pmin(0, an/mr*(1-Satm)*(1-fs))
      dCres_13C<--pmin(0, an/mr*Satm*(1-fs))
      dCO2_12C<-r_12C
      dCO2_13C<-r_13C

      return(list(c(dR_12C, dR_13C,
                    dS_12C, dS_13C,
                    dDOC_12C, dDOC_13C,
                    dCres_12C, dCres_13C,
                    dCO2_12C, dCO2_13C), Cmic_12C=Cmic_12C, Cmic_13C=Cmic_13C))

    })
  }

  #2. With glucose addition
  db_modell<-function(time, state, pars){
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
      #Cres_12C - 12C in residual pool - unsoluble part of death cells
      #Cres_13C - 13C in residual pool - unsoluble part of death cells
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
      Cu_glucose=Vmax_glucose*(S_12C+S_13C)*(G_12C+G_13C)/((G_12C+G_13C)+Km_glucose*(1+(DOC_12C+DOC_13C)/Km_DOC))#
      #Organic carbon uptake rate - DOC
      Cu_DOC=Vmax_DOC*(S_12C+S_13C)*(DOC_12C+DOC_13C)/((DOC_12C+DOC_13C)+Km_DOC*(1+(G_12C+G_13C)/Km_glucose))#

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
      #Maintnance respiration is defined by the concentration of Structures, but at the same time, mobilization rate of Reserves.
      #When an is negative, all C mobilized from pool of Reserves are respired, which doesn't need to correspond with
      #the maintnance requirements. In that case, maintnance respiration is lower than it should be and Structers are dying.
      r_12C<-(1-Ac_glucose)*Cu_glucose*(1-Gatm)+(1-Ac_DOC)*Cu_DOC*(1-DOCatm)+pmax((1-Yu)*an*(1-Ratm), 0)+ifelse(an>0, m*(1-Ratm), f*R_12C)
      r_13C<-(1-Ac_glucose)*Cu_glucose*Gatm+(1-Ac_DOC)*Cu_DOC*DOCatm+pmax((1-Yu)*an*Ratm, 0)+ifelse(an>0, m*Ratm, f*R_13C)

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
      dS_12C<-pmax(an*Yu*(1-Ratm), 0)+pmin(0, an/mr*(1-Satm))
      dS_13C<-pmax(an*Yu*Ratm, 0)+pmin(0, an/mr*Satm)
      dG_12C<--Cu_glucose*(1-Gatm)
      dG_13C<--Cu_glucose*Gatm
      dDOC_12C<--Cu_DOC*(1-DOCatm)-pmin(0, an/mr*(1-Satm)*fs)
      dDOC_13C<--Cu_DOC*DOCatm-pmin(0, an/mr*Satm*fs)
      dCres_12C<--pmin(0, an/mr*(1-Satm)*(1-fs))
      dCres_13C<--pmin(0, an/mr*Satm*(1-fs))
      dCO2_12C<-r_12C
      dCO2_13C<-r_13C

      return(list(c(dR_12C, dR_13C,
                    dS_12C, dS_13C,
                    dG_12C, dG_13C,
                    dDOC_12C, dDOC_13C,
                    dCres_12C, dCres_13C,
                    dCO2_12C, dCO2_13C), Cmic_12C=Cmic_12C, Cmic_13C=Cmic_13C))

    })
  }
  #define names of parameters
  #parameters fr and fs are not estimated, but copied to the parameter vector from dataset
  parnames<-c("Vmax_glucose", "Vmax_DOC",
              "Km_glucose", "Km_DOC",
              "Ac_glucose", "Ac_DOC",
              "mr", "f", "Yu", "Rinit", "Ratm_init")#fr and fs will be supplied additionally

  #All treatments are analysed sequentially. Parallel for loop computing function is defined here.
  epar<-function(unlabelled, labelled){

    #First, minimization ("cost") function is defined
    cost<-function(x){

      par<-x[1:length(parnames)]
      names(par)<-parnames

      #separate vector of parameters into "model parameters" and "initial states parameters"
      #model parameters are further separated into parameters for "db_modelun" (mparun) and "db_modell" (mparl)
      #model parameters
      mparl<-c(par[1:9], "fr"=as.numeric(unlabelled[1, "fr"]), "fs"=as.numeric(unlabelled[1, "fs"]))
      mparun<-mparl[-c(1, 3, 5)]
      #initial states parameters
      ipar<-par[10:11]

      #Extracting initial concentration of state variables from data
      #The initial abundance of Reserves and Strctures in microbial biomass is not known.
      #Therefore initial concentration of Reserves ("Rinit") is estimated together with all model parameters.
      #It is also not known if the Reserves and Structures have the same isotopic signal (i.e. Ratm/Satm).
      #Therefore, one more parameter - Ratm_init, has to be defined and estimated.
      #Every experimental treatment has 4 replicates, which differ in intial states a little bit.
      #However, the model parameters should be same for all replicates.
      #Therefore, 4 initial states are defined and 4 model simulations are run.
      #These are merged together and the differences between simulations and measurements are minimized
      #First set of initial parameters. The fast initial glucose sorption is accounted for.
      R_12Cinit1=ipar[["Rinit"]]*(1-ipar[["Ratm_init"]])
      R_13Cinit1=ipar[["Rinit"]]*ipar[["Ratm_init"]]
      S_12Cinit1=(as.numeric(unlabelled[1, "Cmic"])*(1-as.numeric(unlabelled[1, "Cmicatm"]))-mparl[["fr"]]*R_12Cinit1)/mparl[["fs"]]
      S_13Cinit1=(as.numeric(unlabelled[1, "Cmic"])*as.numeric(unlabelled[1, "Cmicatm"])-mparl[["fr"]]*R_13Cinit1)/mparl[["fs"]]
      G_12Cinit1=as.numeric(labelled[1, "Glcorrected"])*(1-0.06296353)
      G_13Cinit1=as.numeric(labelled[1, "Glcorrected"])*0.06296353
      DOC_12Cinit1=as.numeric(unlabelled[1, "DOC2"])*(1-as.numeric(unlabelled[1, "DOCatm"]))
      DOC_13Cinit1=as.numeric(unlabelled[1, "DOC2"])*as.numeric(unlabelled[1, "DOCatm"])
      #Second set of initial parameters. Parameter "Rinit" is scaled by a ratio between the
      #microbial biomass in replicate 1 and 2.
      R_12Cinit2=ipar[["Rinit"]]*(1-ipar[["Ratm_init"]])*(as.numeric(unlabelled[2, "Cmic"])/as.numeric(unlabelled[1, "Cmic"]))
      R_13Cinit2=ipar[["Rinit"]]*ipar[["Ratm_init"]]*(as.numeric(unlabelled[2, "Cmic"])/as.numeric(unlabelled[1, "Cmic"]))
      S_12Cinit2=(as.numeric(unlabelled[2, "Cmic"])*(1-as.numeric(unlabelled[2, "Cmicatm"]))-mparl[["fr"]]*R_12Cinit2)/mparl[["fs"]]
      S_13Cinit2=(as.numeric(unlabelled[2, "Cmic"])*as.numeric(unlabelled[2, "Cmicatm"])-mparl[["fr"]]*R_13Cinit2)/mparl[["fs"]]
      G_12Cinit2=as.numeric(labelled[2, "Glcorrected"])*(1-0.06296353)
      G_13Cinit2=as.numeric(labelled[2, "Glcorrected"])*0.06296353
      DOC_12Cinit2=as.numeric(unlabelled[2, "DOC2"])*(1-as.numeric(unlabelled[2, "DOCatm"]))
      DOC_13Cinit2=as.numeric(unlabelled[2, "DOC2"])*as.numeric(unlabelled[2, "DOCatm"])
      #Third set of initial parameters. Parameter "Rinit" is scaled by a ratio between the
      #microbial biomass in replicate 1 and 3.
      R_12Cinit3=ipar[["Rinit"]]*(1-ipar[["Ratm_init"]])*(as.numeric(unlabelled[3, "Cmic"])/as.numeric(unlabelled[1, "Cmic"]))
      R_13Cinit3=ipar[["Rinit"]]*ipar[["Ratm_init"]]*(as.numeric(unlabelled[3, "Cmic"])/as.numeric(unlabelled[1, "Cmic"]))
      S_12Cinit3=(as.numeric(unlabelled[3, "Cmic"])*(1-as.numeric(unlabelled[3, "Cmicatm"]))-mparl[["fr"]]*R_12Cinit3)/mparl[["fs"]]
      S_13Cinit3=(as.numeric(unlabelled[3, "Cmic"])*as.numeric(unlabelled[3, "Cmicatm"])-mparl[["fr"]]*R_13Cinit3)/mparl[["fs"]]
      G_12Cinit3=as.numeric(labelled[3, "Glcorrected"])*(1-0.06296353)
      G_13Cinit3=as.numeric(labelled[3, "Glcorrected"])*0.06296353
      DOC_12Cinit3=as.numeric(unlabelled[3, "DOC2"])*(1-as.numeric(unlabelled[3, "DOCatm"]))
      DOC_13Cinit3=as.numeric(unlabelled[3, "DOC2"])*as.numeric(unlabelled[3, "DOCatm"])
      #Fourth set of initial parameters. Parameter "Rinit" is scaled by a ratio between the
      #microbial biomass in replicate 1 and 4.
      R_12Cinit4=ipar[["Rinit"]]*(1-ipar[["Ratm_init"]])*(as.numeric(unlabelled[4, "Cmic"])/as.numeric(unlabelled[1, "Cmic"]))
      R_13Cinit4=ipar[["Rinit"]]*ipar[["Ratm_init"]]*(as.numeric(unlabelled[4, "Cmic"])/as.numeric(unlabelled[1, "Cmic"]))
      S_12Cinit4=(as.numeric(unlabelled[4, "Cmic"])*(1-as.numeric(unlabelled[4, "Cmicatm"]))-mparl[["fr"]]*R_12Cinit4)/mparl[["fs"]]
      S_13Cinit4=(as.numeric(unlabelled[4, "Cmic"])*as.numeric(unlabelled[4, "Cmicatm"])-mparl[["fr"]]*R_13Cinit4)/mparl[["fs"]]
      G_12Cinit4=as.numeric(labelled[4, "Glcorrected"])*(1-0.06296353)
      G_13Cinit4=as.numeric(labelled[4, "Glcorrected"])*0.06296353
      DOC_12Cinit4=as.numeric(unlabelled[4, "DOC2"])*(1-as.numeric(unlabelled[4, "DOCatm"]))
      DOC_13Cinit4=as.numeric(unlabelled[4, "DOC2"])*as.numeric(unlabelled[4, "DOCatm"])
      #Cres and CO2 pools are initialy 0

      #8 model simulation is run. 4 simulations for unlabelled treatments ("db_modelun") and 4 for labelled tratments ("db_modell").
      #Thus, it is expected that the only difference between labelled and unlabelled tretaments is
      #the presence/absence of glucose as a source of organic carbon.
      #All model simulations are merged together. This has one big advantage - it increases the number
      #of data the model is calibrated against.
      #1. unlabelled treatments
      yhat_1un<-as.data.frame(ode(y=c(R_12C=R_12Cinit1, R_13C=R_13Cinit1,
                                      S_12C=S_12Cinit1, S_13C=S_13Cinit1,
                                      DOC_12C=DOC_12Cinit1, DOC_13C=DOC_13Cinit1,
                                      Cres_12C=0, Cres_13C=0, CO2_12C=0, CO2_13C=0),
                                  parms=mparun, db_modelun, times=seq(0,50)))
      yhat_2un<-as.data.frame(ode(y=c(R_12C=R_12Cinit2, R_13C=R_13Cinit2,
                                      S_12C=S_12Cinit2, S_13C=S_13Cinit2,
                                      DOC_12C=DOC_12Cinit2, DOC_13C=DOC_13Cinit2,
                                      Cres_12C=0, Cres_13C=0, CO2_12C=0, CO2_13C=0),
                                  parms=mparun, db_modelun, times=seq(0,50)))
      yhat_3un<-as.data.frame(ode(y=c(R_12C=R_12Cinit3, R_13C=R_13Cinit3,
                                      S_12C=S_12Cinit3, S_13C=S_13Cinit3,
                                      DOC_12C=DOC_12Cinit3, DOC_13C=DOC_13Cinit3,
                                      Cres_12C=0, Cres_13C=0, CO2_12C=0, CO2_13C=0),
                                  parms=mparun, db_modelun, times=seq(0,50)))
      yhat_4un<-as.data.frame(ode(y=c(R_12C=R_12Cinit4, R_13C=R_13Cinit4,
                                      S_12C=S_12Cinit4, S_13C=S_13Cinit4,
                                      DOC_12C=DOC_12Cinit4, DOC_13C=DOC_13Cinit4,
                                      Cres_12C=0, Cres_13C=0, CO2_12C=0, CO2_13C=0),
                                  parms=mparun, db_modelun, times=seq(0,50)))
      #2. labelled treatments
      yhat_1l<-as.data.frame(ode(y=c(R_12C=R_12Cinit1, R_13C=R_13Cinit1,
                                     S_12C=S_12Cinit1, S_13C=S_13Cinit1,
                                     G_12C=G_12Cinit1, G_13C=G_13Cinit1,
                                     DOC_12C=DOC_12Cinit1, DOC_13C=DOC_13Cinit1,
                                     Cres_12C=0, Cres_13C=0, CO2_12C=0, CO2_13C=0),
                                 parms=mparl, db_modell, times=seq(0,50)))
      yhat_2l<-as.data.frame(ode(y=c(R_12C=R_12Cinit2, R_13C=R_13Cinit2,
                                     S_12C=S_12Cinit2, S_13C=S_13Cinit2,
                                     G_12C=G_12Cinit2, G_13C=G_13Cinit2,
                                     DOC_12C=DOC_12Cinit2, DOC_13C=DOC_13Cinit2,
                                     Cres_12C=0, Cres_13C=0, CO2_12C=0, CO2_13C=0),
                                 parms=mparl, db_modell, times=seq(0,50)))
      yhat_3l<-as.data.frame(ode(y=c(R_12C=R_12Cinit3, R_13C=R_13Cinit3,
                                     S_12C=S_12Cinit3, S_13C=S_13Cinit3,
                                     G_12C=G_12Cinit3, G_13C=G_13Cinit3,
                                     DOC_12C=DOC_12Cinit3, DOC_13C=DOC_13Cinit3,
                                     Cres_12C=0, Cres_13C=0, CO2_12C=0, CO2_13C=0),
                                 parms=mparl, db_modell, times=seq(0,50)))
      yhat_4l<-as.data.frame(ode(y=c(R_12C=R_12Cinit4, R_13C=R_13Cinit4,
                                     S_12C=S_12Cinit4, S_13C=S_13Cinit4,
                                     G_12C=G_12Cinit4, G_13C=G_13Cinit4,
                                     DOC_12C=DOC_12Cinit4, DOC_13C=DOC_13Cinit4,
                                     Cres_12C=0, Cres_13C=0, CO2_12C=0, CO2_13C=0),
                                 parms=mparl, db_modell, times=seq(0,50)))
      #select only times 0 and 48 at which the measurements were done
      yhat_1un<-yhat_1un %>% filter(time==0 | time==48)
      yhat_2un<-yhat_2un %>% filter(time==0 | time==48)
      yhat_3un<-yhat_3un %>% filter(time==0 | time==48)
      yhat_4un<-yhat_4un %>% filter(time==0 | time==48)
      yhat_1l<-yhat_1l %>% filter(time==0 | time==48)
      yhat_2l<-yhat_2l %>% filter(time==0 | time==48)
      yhat_3l<-yhat_3l %>% filter(time==0 | time==48)
      yhat_4l<-yhat_4l %>% filter(time==0 | time==48)

      #combine all unlabelled simulations together
      yhat_allun<-rbind(yhat_1un, yhat_2un, yhat_3un, yhat_4un)
      #order them by time
      yhat_allun<-yhat_allun[order(yhat_allun$time), ]

      #variables that were measured in the experiment are extracted
      yhatun<-select(yhat_allun, c("time", "DOC_12C", "DOC_13C", "CO2_12C", "CO2_13C", "Cmic_12C", "Cmic_13C"))

      #convert the simulated dataset into long format data frame
      Yhatun<-melt(yhatun, id.vars = "time")

      #extract the relevant measurements from the dataset
      obsun<-data.frame(time = yhatun$time,
                        DOC_12C = c(as.numeric(unlabelled[c(1:4), "DOC2"])*(1-as.numeric(unlabelled[c(1:4), "DOCatm"])),
                                    as.numeric(unlabelled[c(5:8), "DOC2"])*(1-as.numeric(unlabelled[c(5:8), "DOCatm"]))),
                        DOC_13C = c(as.numeric(unlabelled[c(1:4), "DOC2"])*as.numeric(unlabelled[c(1:4), "DOCatm"]),
                                    as.numeric(unlabelled[c(5:8), "DOC2"])*as.numeric(unlabelled[c(5:8), "DOCatm"])),
                        CO2_12C = c(rep(0, times=4),
                                    as.numeric(unlabelled[c(5:8), "CCO2c"])*(1-as.numeric(unlabelled[c(5:8), "CO2atm"]))),
                        CO2_13C = c(rep(0, times=4),
                                    as.numeric(unlabelled[c(5:8), "CCO2c"])*as.numeric(unlabelled[c(5:8), "CO2atm"])),
                        Cmic_12C = c(as.numeric(unlabelled[c(1:4), "Cmic"])*(1-as.numeric(unlabelled[c(1:4), "Cmicatm"])),
                                     as.numeric(unlabelled[c(5:8), "Cmic"])*(1-as.numeric(unlabelled[c(5:8), "Cmicatm"]))),
                        Cmic_13C = c(as.numeric(unlabelled[c(1:4), "Cmic"])*as.numeric(unlabelled[c(1:4), "Cmicatm"]),
                                     as.numeric(unlabelled[c(5:8), "Cmic"])*as.numeric(unlabelled[c(5:8), "Cmicatm"])))

      #convert dataset with measurements to long format and copy the "value" column to Yhat data frame
      Yhatun$obs<-melt(obsun, id.vars = "time")[, "value"]

      #combine all labelled simulations together
      yhat_alll<-rbind(yhat_1l, yhat_2l, yhat_3l, yhat_4l)
      #order them by time
      yhat_alll<-yhat_alll[order(yhat_alll$time), ]

      #variables that were measured in the experiment are extracted
      yhatl<-select(yhat_alll, c("time", "G_12C", "G_13C", "DOC_12C", "DOC_13C", "CO2_12C", "CO2_13C", "Cmic_12C", "Cmic_13C"))

      #convert the simulated dataset into long format data frame
      Yhatl<-melt(yhatl, id.vars = "time")

      #extract the relevant measurements from the dataset
      obsl<-data.frame(time = yhatl$time,
                       G_12C = c(as.numeric(labelled[c(1:4), "Glcorrected"])*(1-0.06296353),
                                 as.numeric(labelled[c(5:8), "DOCg"])*(1-0.06296353)),
                       G_13C = c(as.numeric(labelled[c(1:4), "Glcorrected"])*0.06296353,
                                 as.numeric(labelled[c(5:8), "DOCg"])*0.06296353),
                       DOC_12C = c(as.numeric(labelled[c(1:4), "DOC2"])*(1-as.numeric(labelled[c(1:4), "DOCatm"])),
                                   (as.numeric(labelled[c(5:8), "DOC2"])-as.numeric(labelled[c(5:8), "DOCg"]))*(1-as.numeric(unlabelled[c(5:8), "DOCatm"]))),
                       DOC_13C = c(as.numeric(labelled[c(1:4), "DOC2"])*as.numeric(labelled[c(1:4), "DOCatm"]),
                                   (as.numeric(labelled[c(5:8), "DOC2"])-as.numeric(labelled[c(5:8), "DOCg"]))*as.numeric(unlabelled[c(5:8), "DOCatm"])),
                       CO2_12C = c(rep(0, times=4),
                                   as.numeric(labelled[c(5:8), "CCO2c"])*(1-as.numeric(labelled[c(5:8), "CO2atm"]))),
                       CO2_13C = c(rep(0, times=4),
                                   as.numeric(labelled[c(5:8), "CCO2c"])*as.numeric(labelled[c(5:8), "CO2atm"])),
                       Cmic_12C = c(as.numeric(labelled[c(1:4), "Cmic"])*(1-as.numeric(labelled[c(1:4), "Cmicatm"])),
                                    as.numeric(labelled[c(5:8), "Cmic"])*(1-as.numeric(labelled[c(5:8), "Cmicatm"]))),
                       Cmic_13C = c(as.numeric(labelled[c(1:4), "Cmic"])*as.numeric(labelled[c(1:4), "Cmicatm"]),
                                    as.numeric(labelled[c(5:8), "Cmic"])*as.numeric(labelled[c(5:8), "Cmicatm"])))

      #convert dataset with measurements to long format and copy the "value" column to Yhat data frame
      Yhatl$obs<-melt(obsl, id.vars = "time")[, "value"]

      #merge all simulations together
      Yhat<-rbind(Yhatun, Yhatl)

      #add the weighting factor
      #I want to have the weighting factor to be proportional to mean of the particular variables
      yweights<-Yhat %>% group_by(variable) %>% summarise(weights=mean(value, na.rm = T))

      #match yweights with the Yhat data frame
      Yhat<-merge(Yhat, yweights, by.x=c("variable"), by.y=c("variable"))

      #now, the root mean square error is calculated
      RMSE<-as.numeric(Yhat %>% summarise(RMSE=sum((((value-obs)/mean(weights))^2), na.rm = T)))

      return(RMSE)
    }

    #Second, model parameters are estimated
    #The calibration is done in two consecutive steps
    #First, approximate model parameters are estimated using MCMC algorithm
    #approximate parameter estimation is done by MCMC method
    #Initial guess and upper/lower limits of parameters estimates are taken from
    #previous calibration iteration and supplied in function call as "minpar", "maxpar",
    #"initPL" and "initCT"
    Rinit_guess=as.numeric(unlabelled[1, "Cmic"])*(1-as.numeric(unlabelled[1, "Cmicatm"]))
    Ratm_init_guess=as.numeric(unlabelled[1, "Cmicatm"])

    mpar_guess<-c(as.numeric(unlabelled[1, "Plesne"])*initPL+as.numeric(unlabelled[1, "Certovo"])*initCT,
                  0.1*Rinit_guess, Ratm_init_guess)
    llimit<-c(as.numeric(minpar), 1e-3*Rinit_guess, 0.5*Ratm_init_guess)
    ulimit<-c(as.numeric(maxpar), 0.95*Rinit_guess, 1.2*Ratm_init_guess)

    names(mpar_guess)<-parnames
    names(llimit)<-parnames
    names(ulimit)<-parnames

    par_mcmc<-modMCMC(f=cost,
                      p=mpar_guess,
                      lower=llimit,
                      upper=ulimit, niter=10000)

    #lower and upper limits for parameters estimates are extracted
    pl<-summary(par_mcmc)["min",]
    pu<-summary(par_mcmc)["max",]

    #these limits are used to find global optimum by DEoptim algorithm
    opt_par<-DEoptim(fn=cost, lower=pl, upper=pu,
                     control = c(itermax = 10000, steptol = 50, reltol = 1e-8,
                                 trace=FALSE, strategy=3, NP=250))
    out1<-list(mcmc_out=par_mcmc,
                   pars=opt_par)

    return(out1)
  }

  #Now, the parameters are estimated for each treatment
  epar_out<-foreach(i=seq(1:5), .combine=list, .multicombine = TRUE,
                    .packages=c("FME", "dplyr", "DEoptim", "reshape")) %dopar% {

                      epar(unlabelled = dataset[dataset$id==i, ],
                           labelled = dataset[dataset$id==i+10, ])

                    }

  #parameters are further extracted and merged into one data frame to allow
  #subsequent calculation of goodness of simulations.
  pars_all<-as.data.frame(rbind(epar_out[[1]]$pars$optim$bestmem,
                                epar_out[[2]]$pars$optim$bestmem,
                                epar_out[[3]]$pars$optim$bestmem,
                                epar_out[[4]]$pars$optim$bestmem,
                                epar_out[[5]]$pars$optim$bestmem))
                                #epar_out[[6]]$pars$optim$bestmem,
                                #epar_out[[7]]$pars$optim$bestmem,
                                #epar_out[[8]]$pars$optim$bestmem,
                                #epar_out[[9]]$pars$optim$bestmem,
                                #epar_out[[10]]$pars$optim$bestmem))
  #add labels
  # pars_all$Plesne<-rep(c(1, 0.75, 0.5, 0.25, 0), times=2)
  # pars_all$Certovo<-rep(c(0, 0.25, 0.5, 0.75, 1), times=2)
  # pars_all$horizon<-c(rep("Litter", times=5), rep("Organic soil", times=5))
  # pars_all$id<-seq(1,10)

  pars_all$Plesne<-rep(c(1, 0.75, 0.5, 0.25, 0), times=1)
  pars_all$Certovo<-rep(c(0, 0.25, 0.5, 0.75, 1), times=1)
  pars_all$horizon<-rep("Litter", times=5)
  pars_all$id<-seq(1,5)

  #################################################################################################
  #Goodness of model simulations is calculated here for each treatment and each variable separately
  #To do so, "good" function is defined. This function pretty much similar to "cost" function
  good<-function(x, unlabelled, labelled){

    par<-x[1:length(parnames)]
    names(par)<-parnames

    #separate vector of parameters into "model parameters" and "initial states parameters"
    #model parameters are further separated into parameters for "db_modelun" (mparun) and "db_modell" (mparl)
    #model parameters
    mparl<-c(par[1:9], "fr"=as.numeric(unlabelled[1, "fr"]), "fs"=as.numeric(unlabelled[1, "fs"]))
    mparun<-mparl[-c(1, 3, 5)]
    #initial states parameters
    ipar<-par[10:11]

    #Extracting initial concentration of state variables from data
    #The initial abundance of Reserves and Strctures in microbial biomass is not known.
    #Therefore initial concentration of Reserves ("Rinit") is estimated together with all model parameters.
    #It is also not known if the Reserves and Structures have the same isotopic signal (i.e. Ratm/Satm).
    #Therefore, one more parameter - Ratm_init, has to be defined and estimated.
    #Every experimental treatment has 4 replicates, which differ in intial states a little bit.
    #However, the model parameters should be same for all replicates.
    #Therefore, 4 initial states are defined and 4 model simulations are run.
    #These are merged together and the differences between simulations and measurements are minimized
    #First set of initial parameters. The fast initial glucose sorption is accounted for.
    R_12Cinit1=ipar[["Rinit"]]*(1-ipar[["Ratm_init"]])
    R_13Cinit1=ipar[["Rinit"]]*ipar[["Ratm_init"]]
    S_12Cinit1=(as.numeric(unlabelled[1, "Cmic"])*(1-as.numeric(unlabelled[1, "Cmicatm"]))-mparl[["fr"]]*R_12Cinit1)/mparl[["fs"]]
    S_13Cinit1=(as.numeric(unlabelled[1, "Cmic"])*as.numeric(unlabelled[1, "Cmicatm"])-mparl[["fr"]]*R_13Cinit1)/mparl[["fs"]]
    G_12Cinit1=as.numeric(labelled[1, "Glcorrected"])*(1-0.06296353)
    G_13Cinit1=as.numeric(labelled[1, "Glcorrected"])*0.06296353
    DOC_12Cinit1=as.numeric(unlabelled[1, "DOC2"])*(1-as.numeric(unlabelled[1, "DOCatm"]))
    DOC_13Cinit1=as.numeric(unlabelled[1, "DOC2"])*as.numeric(unlabelled[1, "DOCatm"])
    #Second set of initial parameters. Parameter "Rinit" is scaled by a ratio between the
    #microbial biomass in replicate 1 and 2.
    R_12Cinit2=ipar[["Rinit"]]*(1-ipar[["Ratm_init"]])*(as.numeric(unlabelled[2, "Cmic"])/as.numeric(unlabelled[1, "Cmic"]))
    R_13Cinit2=ipar[["Rinit"]]*ipar[["Ratm_init"]]*(as.numeric(unlabelled[2, "Cmic"])/as.numeric(unlabelled[1, "Cmic"]))
    S_12Cinit2=(as.numeric(unlabelled[2, "Cmic"])*(1-as.numeric(unlabelled[2, "Cmicatm"]))-mparl[["fr"]]*R_12Cinit2)/mparl[["fs"]]
    S_13Cinit2=(as.numeric(unlabelled[2, "Cmic"])*as.numeric(unlabelled[2, "Cmicatm"])-mparl[["fr"]]*R_13Cinit2)/mparl[["fs"]]
    G_12Cinit2=as.numeric(labelled[2, "Glcorrected"])*(1-0.06296353)
    G_13Cinit2=as.numeric(labelled[2, "Glcorrected"])*0.06296353
    DOC_12Cinit2=as.numeric(unlabelled[2, "DOC2"])*(1-as.numeric(unlabelled[2, "DOCatm"]))
    DOC_13Cinit2=as.numeric(unlabelled[2, "DOC2"])*as.numeric(unlabelled[2, "DOCatm"])
    #Third set of initial parameters. Parameter "Rinit" is scaled by a ratio between the
    #microbial biomass in replicate 1 and 3.
    R_12Cinit3=ipar[["Rinit"]]*(1-ipar[["Ratm_init"]])*(as.numeric(unlabelled[3, "Cmic"])/as.numeric(unlabelled[1, "Cmic"]))
    R_13Cinit3=ipar[["Rinit"]]*ipar[["Ratm_init"]]*(as.numeric(unlabelled[3, "Cmic"])/as.numeric(unlabelled[1, "Cmic"]))
    S_12Cinit3=(as.numeric(unlabelled[3, "Cmic"])*(1-as.numeric(unlabelled[3, "Cmicatm"]))-mparl[["fr"]]*R_12Cinit3)/mparl[["fs"]]
    S_13Cinit3=(as.numeric(unlabelled[3, "Cmic"])*as.numeric(unlabelled[3, "Cmicatm"])-mparl[["fr"]]*R_13Cinit3)/mparl[["fs"]]
    G_12Cinit3=as.numeric(labelled[3, "Glcorrected"])*(1-0.06296353)
    G_13Cinit3=as.numeric(labelled[3, "Glcorrected"])*0.06296353
    DOC_12Cinit3=as.numeric(unlabelled[3, "DOC2"])*(1-as.numeric(unlabelled[3, "DOCatm"]))
    DOC_13Cinit3=as.numeric(unlabelled[3, "DOC2"])*as.numeric(unlabelled[3, "DOCatm"])
    #Fourth set of initial parameters. Parameter "Rinit" is scaled by a ratio between the
    #microbial biomass in replicate 1 and 4.
    R_12Cinit4=ipar[["Rinit"]]*(1-ipar[["Ratm_init"]])*(as.numeric(unlabelled[4, "Cmic"])/as.numeric(unlabelled[1, "Cmic"]))
    R_13Cinit4=ipar[["Rinit"]]*ipar[["Ratm_init"]]*(as.numeric(unlabelled[4, "Cmic"])/as.numeric(unlabelled[1, "Cmic"]))
    S_12Cinit4=(as.numeric(unlabelled[4, "Cmic"])*(1-as.numeric(unlabelled[4, "Cmicatm"]))-mparl[["fr"]]*R_12Cinit4)/mparl[["fs"]]
    S_13Cinit4=(as.numeric(unlabelled[4, "Cmic"])*as.numeric(unlabelled[4, "Cmicatm"])-mparl[["fr"]]*R_13Cinit4)/mparl[["fs"]]
    G_12Cinit4=as.numeric(labelled[4, "Glcorrected"])*(1-0.06296353)
    G_13Cinit4=as.numeric(labelled[4, "Glcorrected"])*0.06296353
    DOC_12Cinit4=as.numeric(unlabelled[4, "DOC2"])*(1-as.numeric(unlabelled[4, "DOCatm"]))
    DOC_13Cinit4=as.numeric(unlabelled[4, "DOC2"])*as.numeric(unlabelled[4, "DOCatm"])
    #Cres and CO2 pools are initialy 0

    #8 model simulation is run. 4 simulations for unlabelled treatments ("db_modelun") and 4 for labelled tratments ("db_modell").
    #Thus, it is expected that the only difference between labelled and unlabelled tretaments is
    #the presence/absence of glucose as a source of organic carbon.
    #All model simulations are merged together.
    #1. unlabelled treatments
    yhat_1un<-as.data.frame(ode(y=c(R_12C=R_12Cinit1, R_13C=R_13Cinit1,
                                    S_12C=S_12Cinit1, S_13C=S_13Cinit1,
                                    DOC_12C=DOC_12Cinit1, DOC_13C=DOC_13Cinit1,
                                    Cres_12C=0, Cres_13C=0, CO2_12C=0, CO2_13C=0),
                                parms=mparun, db_modelun, times=seq(0,50)))
    yhat_2un<-as.data.frame(ode(y=c(R_12C=R_12Cinit2, R_13C=R_13Cinit2,
                                    S_12C=S_12Cinit2, S_13C=S_13Cinit2,
                                    DOC_12C=DOC_12Cinit2, DOC_13C=DOC_13Cinit2,
                                    Cres_12C=0, Cres_13C=0, CO2_12C=0, CO2_13C=0),
                                parms=mparun, db_modelun, times=seq(0,50)))
    yhat_3un<-as.data.frame(ode(y=c(R_12C=R_12Cinit3, R_13C=R_13Cinit3,
                                    S_12C=S_12Cinit3, S_13C=S_13Cinit3,
                                    DOC_12C=DOC_12Cinit3, DOC_13C=DOC_13Cinit3,
                                    Cres_12C=0, Cres_13C=0, CO2_12C=0, CO2_13C=0),
                                parms=mparun, db_modelun, times=seq(0,50)))
    yhat_4un<-as.data.frame(ode(y=c(R_12C=R_12Cinit4, R_13C=R_13Cinit4,
                                    S_12C=S_12Cinit4, S_13C=S_13Cinit4,
                                    DOC_12C=DOC_12Cinit4, DOC_13C=DOC_13Cinit4,
                                    Cres_12C=0, Cres_13C=0, CO2_12C=0, CO2_13C=0),
                                parms=mparun, db_modelun, times=seq(0,50)))
    #2. labelled treatments
    yhat_1l<-as.data.frame(ode(y=c(R_12C=R_12Cinit1, R_13C=R_13Cinit1,
                                   S_12C=S_12Cinit1, S_13C=S_13Cinit1,
                                   G_12C=G_12Cinit1, G_13C=G_13Cinit1,
                                   DOC_12C=DOC_12Cinit1, DOC_13C=DOC_13Cinit1,
                                   Cres_12C=0, Cres_13C=0, CO2_12C=0, CO2_13C=0),
                               parms=mparl, db_modell, times=seq(0,50)))
    yhat_2l<-as.data.frame(ode(y=c(R_12C=R_12Cinit2, R_13C=R_13Cinit2,
                                   S_12C=S_12Cinit2, S_13C=S_13Cinit2,
                                   G_12C=G_12Cinit2, G_13C=G_13Cinit2,
                                   DOC_12C=DOC_12Cinit2, DOC_13C=DOC_13Cinit2,
                                   Cres_12C=0, Cres_13C=0, CO2_12C=0, CO2_13C=0),
                               parms=mparl, db_modell, times=seq(0,50)))
    yhat_3l<-as.data.frame(ode(y=c(R_12C=R_12Cinit3, R_13C=R_13Cinit3,
                                   S_12C=S_12Cinit3, S_13C=S_13Cinit3,
                                   G_12C=G_12Cinit3, G_13C=G_13Cinit3,
                                   DOC_12C=DOC_12Cinit3, DOC_13C=DOC_13Cinit3,
                                   Cres_12C=0, Cres_13C=0, CO2_12C=0, CO2_13C=0),
                               parms=mparl, db_modell, times=seq(0,50)))
    yhat_4l<-as.data.frame(ode(y=c(R_12C=R_12Cinit4, R_13C=R_13Cinit4,
                                   S_12C=S_12Cinit4, S_13C=S_13Cinit4,
                                   G_12C=G_12Cinit4, G_13C=G_13Cinit4,
                                   DOC_12C=DOC_12Cinit4, DOC_13C=DOC_13Cinit4,
                                   Cres_12C=0, Cres_13C=0, CO2_12C=0, CO2_13C=0),
                               parms=mparl, db_modell, times=seq(0,50)))
    #select only times 0 and 48 at which the measurements were done
    yhat_1un<-yhat_1un %>% filter(time==0 | time==48)
    yhat_2un<-yhat_2un %>% filter(time==0 | time==48)
    yhat_3un<-yhat_3un %>% filter(time==0 | time==48)
    yhat_4un<-yhat_4un %>% filter(time==0 | time==48)
    yhat_1l<-yhat_1l %>% filter(time==0 | time==48)
    yhat_2l<-yhat_2l %>% filter(time==0 | time==48)
    yhat_3l<-yhat_3l %>% filter(time==0 | time==48)
    yhat_4l<-yhat_4l %>% filter(time==0 | time==48)

    #combine all unlabelled simulations together
    yhat_allun<-rbind(yhat_1un, yhat_2un, yhat_3un, yhat_4un)
    #order them by time
    yhat_allun<-yhat_allun[order(yhat_allun$time), ]

    #variables that were measured in the experiment are extracted
    yhatun<-select(yhat_allun, c("time", "DOC_12C", "DOC_13C", "CO2_12C", "CO2_13C", "Cmic_12C", "Cmic_13C"))

    #convert the simulated dataset into long format data frame
    Yhatun<-melt(yhatun, id.vars = "time")

    #extract the relevant measurements from the dataset
    obsun<-data.frame(time = yhatun$time,
                      DOC_12C = c(as.numeric(unlabelled[c(1:4), "DOC2"])*(1-as.numeric(unlabelled[c(1:4), "DOCatm"])),
                                  as.numeric(unlabelled[c(5:8), "DOC2"])*(1-as.numeric(unlabelled[c(5:8), "DOCatm"]))),
                      DOC_13C = c(as.numeric(unlabelled[c(1:4), "DOC2"])*as.numeric(unlabelled[c(1:4), "DOCatm"]),
                                  as.numeric(unlabelled[c(5:8), "DOC2"])*as.numeric(unlabelled[c(5:8), "DOCatm"])),
                      CO2_12C = c(rep(0, times=4),
                                  as.numeric(unlabelled[c(5:8), "CCO2c"])*(1-as.numeric(unlabelled[c(5:8), "CO2atm"]))),
                      CO2_13C = c(rep(0, times=4),
                                  as.numeric(unlabelled[c(5:8), "CCO2c"])*as.numeric(unlabelled[c(5:8), "CO2atm"])),
                      Cmic_12C = c(as.numeric(unlabelled[c(1:4), "Cmic"])*(1-as.numeric(unlabelled[c(1:4), "Cmicatm"])),
                                   as.numeric(unlabelled[c(5:8), "Cmic"])*(1-as.numeric(unlabelled[c(5:8), "Cmicatm"]))),
                      Cmic_13C = c(as.numeric(unlabelled[c(1:4), "Cmic"])*as.numeric(unlabelled[c(1:4), "Cmicatm"]),
                                   as.numeric(unlabelled[c(5:8), "Cmic"])*as.numeric(unlabelled[c(5:8), "Cmicatm"])))

    #convert dataset with measurements to long format and copy the "value" column to Yhat data frame
    Yhatun$obs<-melt(obsun, id.vars = "time")[, "value"]

    #combine all labelled simulations together
    yhat_alll<-rbind(yhat_1l, yhat_2l, yhat_3l, yhat_4l)
    #order them by time
    yhat_alll<-yhat_alll[order(yhat_alll$time), ]

    #variables that were measured in the experiment are extracted
    yhatl<-select(yhat_alll, c("time", "G_12C", "G_13C", "DOC_12C", "DOC_13C", "CO2_12C", "CO2_13C", "Cmic_12C", "Cmic_13C"))

    #convert the simulated dataset into long format data frame
    Yhatl<-melt(yhatl, id.vars = "time")

    #extract the relevant measurements from the dataset
    obsl<-data.frame(time = yhatl$time,
                     G_12C = c(as.numeric(labelled[c(1:4), "Glcorrected"])*(1-0.06296353),
                               as.numeric(labelled[c(5:8), "DOCg"])*(1-0.06296353)),
                     G_13C = c(as.numeric(labelled[c(1:4), "Glcorrected"])*0.06296353,
                               as.numeric(labelled[c(5:8), "DOCg"])*0.06296353),
                     DOC_12C = c(as.numeric(labelled[c(1:4), "DOC2"])*(1-as.numeric(labelled[c(1:4), "DOCatm"])),
                                 (as.numeric(labelled[c(5:8), "DOC2"])-as.numeric(labelled[c(5:8), "DOCg"]))*(1-as.numeric(unlabelled[c(5:8), "DOCatm"]))),
                     DOC_13C = c(as.numeric(labelled[c(1:4), "DOC2"])*as.numeric(labelled[c(1:4), "DOCatm"]),
                                 (as.numeric(labelled[c(5:8), "DOC2"])-as.numeric(labelled[c(5:8), "DOCg"]))*as.numeric(unlabelled[c(5:8), "DOCatm"])),
                     CO2_12C = c(rep(0, times=4),
                                 as.numeric(labelled[c(5:8), "CCO2c"])*(1-as.numeric(labelled[c(5:8), "CO2atm"]))),
                     CO2_13C = c(rep(0, times=4),
                                 as.numeric(labelled[c(5:8), "CCO2c"])*as.numeric(labelled[c(5:8), "CO2atm"])),
                     Cmic_12C = c(as.numeric(labelled[c(1:4), "Cmic"])*(1-as.numeric(labelled[c(1:4), "Cmicatm"])),
                                  as.numeric(labelled[c(5:8), "Cmic"])*(1-as.numeric(labelled[c(5:8), "Cmicatm"]))),
                     Cmic_13C = c(as.numeric(labelled[c(1:4), "Cmic"])*as.numeric(labelled[c(1:4), "Cmicatm"]),
                                  as.numeric(labelled[c(5:8), "Cmic"])*as.numeric(labelled[c(5:8), "Cmicatm"])))

    #convert dataset with measurements to long format and copy the "value" column to Yhat data frame
    Yhatl$obs<-melt(obsl, id.vars = "time")[, "value"]

    #merge all simulations together
    Yhat<-rbind(Yhatun, Yhatl)

    #add labels
    Yhat$Plesne<-c(rep(as.numeric(unlabelled$Plesne[1]), times = nrow(Yhatun)),
                   rep(as.numeric(unlabelled$Plesne[1]), times = nrow(Yhatl)))
    Yhat$Certovo<-c(rep(as.numeric(unlabelled$Certovo[1]), times = nrow(Yhatun)),
                   rep(as.numeric(unlabelled$Certovo[1]), times = nrow(Yhatl)))
    Yhat$horizon<-c(rep(unlabelled$horizon[1], times = nrow(Yhatun)),
                    rep(unlabelled$horizon[1], times = nrow(Yhatl)))
    Yhat$Treatment<-c(rep("Unlabelled", times = nrow(Yhatun)),
                      rep("Labelled", times = nrow(Yhatl)))

    #several goodness of fit metrics are calulated for each variable
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
  #using "good" function for each treatment in parallel
  good_out<-foreach(i=seq(1:5), .combine=list, .multicombine = TRUE,
                    .packages=c("FME", "dplyr", "DEoptim", "reshape")) %dopar% {

                      good(x=as.numeric(pars_all[pars_all$id==i, c(1:11)]),
                           unlabelled = dataset[dataset$id==i, ],
                           labelled = dataset[dataset$id==i+10, ])

                    }

  #all measurements and simulations are combined together to calculate the goodness
  #of correspondence across all treatments, but for each variable separately
  Yhat_treatments<-as.data.frame(rbind(good_out[[1]]$Yhat, good_out[[2]]$Yhat,
                                       good_out[[3]]$Yhat, good_out[[4]]$Yhat,
                                       good_out[[5]]$Yhat))
                                       # good_out[[6]]$Yhat,
                                       # good_out[[7]]$Yhat, good_out[[8]]$Yhat,
                                       # good_out[[9]]$Yhat, good_out[[10]]$Yhat))
                                       #
  #goodness of fit metrics are calulated in several steps
  #1. First, goodness of correspondence is calculated across all treatments
  Gfit_treatments<-Yhat_treatments %>% group_by(variable) %>% summarise(SSres=sum(((obs-value)^2), na.rm = T),
                                                                        SStot=sum(((obs-mean(obs, na.rm = T))^2), na.rm = T),
                                                                        ll=-sum(((obs-value)^2), na.rm = T)/2/(sd(obs, na.rm = T)^2))
  Gfit_treatments$R2<-with(Gfit_treatments, 1-SSres/SStot)
  Gfit_treatments$N<-rep(110, times=nrow(Gfit_treatments))
  Gfit_treatments$AIC<-with(Gfit_treatments, 2*N-2*ll)
  #2. Second, goodness of correspondence is calculated for unlabelled and labelled treatments separately
  Gfit_labels<-Yhat_treatments %>% group_by(variable, Treatment) %>% summarise(SSres=sum(((obs-value)^2), na.rm = T),
                                                                               SStot=sum(((obs-mean(obs, na.rm = T))^2), na.rm = T),
                                                                               ll=-sum(((obs-value)^2), na.rm = T)/2/(sd(obs, na.rm = T)^2))
  Gfit_labels$R2<-with(Gfit_labels, 1-SSres/SStot)
  Gfit_labels$N<-rep(110, times=nrow(Gfit_labels))
  Gfit_labels$AIC<-with(Gfit_labels, 2*N-2*ll)
  # #3. Third, goodness of correspondence is calculated for each horizon separately
  # Gfit_horizon<-Yhat_treatments %>% group_by(variable, horizon) %>% summarise(SSres=sum(((obs-value)^2), na.rm = T),
  #                                                                             SStot=sum(((obs-mean(obs, na.rm = T))^2), na.rm = T),
  #                                                                             ll=-sum(((obs-value)^2), na.rm = T)/2/(sd(obs, na.rm = T)^2))
  # Gfit_horizon$R2<-with(Gfit_horizon, 1-SSres/SStot)
  # Gfit_horizon$N<-rep(110, times=nrow(Gfit_horizon))
  # Gfit_horizon$AIC<-with(Gfit_horizon, 2*N-2*ll)
  # #4. Fourth, goodness of correspondence is calculated for each horizon and unlabelled/labelled treatment separately
  # Gfit_each<-Yhat_treatments %>% group_by(variable, horizon, Treatment) %>% summarise(SSres=sum(((obs-value)^2), na.rm = T),
  #                                                                                     SStot=sum(((obs-mean(obs, na.rm = T))^2), na.rm = T),
  #                                                                                     ll=-sum(((obs-value)^2), na.rm = T)/2/(sd(obs, na.rm = T)^2))
  # Gfit_each$R2<-with(Gfit_each, 1-SSres/SStot)
  # Gfit_each$N<-rep(110, times=nrow(Gfit_each))
  # Gfit_each$AIC<-with(Gfit_each, 2*N-2*ll)

  #For nice plots, model simulation is run with calibrated model parameters
  #at fine temporal scale. The simulation is done in parallel again using function "simul".
  simul<-function(x, unlabelled, labelled){

    par<-x[1:length(parnames)]
    names(par)<-parnames

    #separate vector of parameters into "model parameters" and "initial states parameters"
    #model parameters are further separated into parameters for "db_modelun" (mparun) and "db_modell" (mparl)
    #model parameters
    mparl<-c(par[1:9], "fr"=as.numeric(unlabelled[1, "fr"]), "fs"=as.numeric(unlabelled[1, "fs"]))
    mparun<-mparl[-c(1, 3, 5)]
    #initial states parameters
    ipar<-par[10:11]

    #Extracting initial concentration of state variables from data
    #The initial abundance of Reserves and Strctures in microbial biomass is not known.
    #Therefore initial concentration of Reserves ("Rinit") is estimated together with all model parameters.
    #It is also not known if the Reserves and Structures have the same isotopic signal (i.e. Ratm/Satm).
    #Therefore, one more parameter - Ratm_init, has to be defined and estimated.
    #Every experimental treatment has 4 replicates, which differ in intial states a little bit.
    #However, the model parameters should be same for all replicates.
    #Therefore, 4 initial states are defined and 4 model simulations are run.
    #These are merged together and the differences between simulations and measurements are minimized
    #First set of initial parameters. The fast initial glucose sorption is accounted for.
    R_12Cinit1=ipar[["Rinit"]]*(1-ipar[["Ratm_init"]])
    R_13Cinit1=ipar[["Rinit"]]*ipar[["Ratm_init"]]
    S_12Cinit1=(as.numeric(unlabelled[1, "Cmic"])*(1-as.numeric(unlabelled[1, "Cmicatm"]))-mparl[["fr"]]*R_12Cinit1)/mparl[["fs"]]
    S_13Cinit1=(as.numeric(unlabelled[1, "Cmic"])*as.numeric(unlabelled[1, "Cmicatm"])-mparl[["fr"]]*R_13Cinit1)/mparl[["fs"]]
    G_12Cinit1=as.numeric(labelled[1, "Glcorrected"])*(1-0.06296353)
    G_13Cinit1=as.numeric(labelled[1, "Glcorrected"])*0.06296353
    DOC_12Cinit1=as.numeric(unlabelled[1, "DOC2"])*(1-as.numeric(unlabelled[1, "DOCatm"]))
    DOC_13Cinit1=as.numeric(unlabelled[1, "DOC2"])*as.numeric(unlabelled[1, "DOCatm"])
    #Second set of initial parameters. Parameter "Rinit" is scaled by a ratio between the
    #microbial biomass in replicate 1 and 2.
    R_12Cinit2=ipar[["Rinit"]]*(1-ipar[["Ratm_init"]])*(as.numeric(unlabelled[2, "Cmic"])/as.numeric(unlabelled[1, "Cmic"]))
    R_13Cinit2=ipar[["Rinit"]]*ipar[["Ratm_init"]]*(as.numeric(unlabelled[2, "Cmic"])/as.numeric(unlabelled[1, "Cmic"]))
    S_12Cinit2=(as.numeric(unlabelled[2, "Cmic"])*(1-as.numeric(unlabelled[2, "Cmicatm"]))-mparl[["fr"]]*R_12Cinit2)/mparl[["fs"]]
    S_13Cinit2=(as.numeric(unlabelled[2, "Cmic"])*as.numeric(unlabelled[2, "Cmicatm"])-mparl[["fr"]]*R_13Cinit2)/mparl[["fs"]]
    G_12Cinit2=as.numeric(labelled[2, "Glcorrected"])*(1-0.06296353)
    G_13Cinit2=as.numeric(labelled[2, "Glcorrected"])*0.06296353
    DOC_12Cinit2=as.numeric(unlabelled[2, "DOC2"])*(1-as.numeric(unlabelled[2, "DOCatm"]))
    DOC_13Cinit2=as.numeric(unlabelled[2, "DOC2"])*as.numeric(unlabelled[2, "DOCatm"])
    #Third set of initial parameters. Parameter "Rinit" is scaled by a ratio between the
    #microbial biomass in replicate 1 and 3.
    R_12Cinit3=ipar[["Rinit"]]*(1-ipar[["Ratm_init"]])*(as.numeric(unlabelled[3, "Cmic"])/as.numeric(unlabelled[1, "Cmic"]))
    R_13Cinit3=ipar[["Rinit"]]*ipar[["Ratm_init"]]*(as.numeric(unlabelled[3, "Cmic"])/as.numeric(unlabelled[1, "Cmic"]))
    S_12Cinit3=(as.numeric(unlabelled[3, "Cmic"])*(1-as.numeric(unlabelled[3, "Cmicatm"]))-mparl[["fr"]]*R_12Cinit3)/mparl[["fs"]]
    S_13Cinit3=(as.numeric(unlabelled[3, "Cmic"])*as.numeric(unlabelled[3, "Cmicatm"])-mparl[["fr"]]*R_13Cinit3)/mparl[["fs"]]
    G_12Cinit3=as.numeric(labelled[3, "Glcorrected"])*(1-0.06296353)
    G_13Cinit3=as.numeric(labelled[3, "Glcorrected"])*0.06296353
    DOC_12Cinit3=as.numeric(unlabelled[3, "DOC2"])*(1-as.numeric(unlabelled[3, "DOCatm"]))
    DOC_13Cinit3=as.numeric(unlabelled[3, "DOC2"])*as.numeric(unlabelled[3, "DOCatm"])
    #Fourth set of initial parameters. Parameter "Rinit" is scaled by a ratio between the
    #microbial biomass in replicate 1 and 4.
    R_12Cinit4=ipar[["Rinit"]]*(1-ipar[["Ratm_init"]])*(as.numeric(unlabelled[4, "Cmic"])/as.numeric(unlabelled[1, "Cmic"]))
    R_13Cinit4=ipar[["Rinit"]]*ipar[["Ratm_init"]]*(as.numeric(unlabelled[4, "Cmic"])/as.numeric(unlabelled[1, "Cmic"]))
    S_12Cinit4=(as.numeric(unlabelled[4, "Cmic"])*(1-as.numeric(unlabelled[4, "Cmicatm"]))-mparl[["fr"]]*R_12Cinit4)/mparl[["fs"]]
    S_13Cinit4=(as.numeric(unlabelled[4, "Cmic"])*as.numeric(unlabelled[4, "Cmicatm"])-mparl[["fr"]]*R_13Cinit4)/mparl[["fs"]]
    G_12Cinit4=as.numeric(labelled[4, "Glcorrected"])*(1-0.06296353)
    G_13Cinit4=as.numeric(labelled[4, "Glcorrected"])*0.06296353
    DOC_12Cinit4=as.numeric(unlabelled[4, "DOC2"])*(1-as.numeric(unlabelled[4, "DOCatm"]))
    DOC_13Cinit4=as.numeric(unlabelled[4, "DOC2"])*as.numeric(unlabelled[4, "DOCatm"])
    #Cres and CO2 pools are initialy 0

    #8 model simulation is run. 4 simulations for unlabelled treatments ("db_modelun") and 4 for labelled tratments ("db_modell").
    #Thus, it is expected that the only difference between labelled and unlabelled tretaments is
    #the presence/absence of glucose as a source of organic carbon.
    #All model simulations are merged together.
    #1. unlabelled treatments
    yhat_1un<-as.data.frame(ode(y=c(R_12C=R_12Cinit1, R_13C=R_13Cinit1,
                                    S_12C=S_12Cinit1, S_13C=S_13Cinit1,
                                    DOC_12C=DOC_12Cinit1, DOC_13C=DOC_13Cinit1,
                                    Cres_12C=0, Cres_13C=0, CO2_12C=0, CO2_13C=0),
                                parms=mparun, db_modelun, times=seq(0,50)))
    yhat_2un<-as.data.frame(ode(y=c(R_12C=R_12Cinit2, R_13C=R_13Cinit2,
                                    S_12C=S_12Cinit2, S_13C=S_13Cinit2,
                                    DOC_12C=DOC_12Cinit2, DOC_13C=DOC_13Cinit2,
                                    Cres_12C=0, Cres_13C=0, CO2_12C=0, CO2_13C=0),
                                parms=mparun, db_modelun, times=seq(0,50)))
    yhat_3un<-as.data.frame(ode(y=c(R_12C=R_12Cinit3, R_13C=R_13Cinit3,
                                    S_12C=S_12Cinit3, S_13C=S_13Cinit3,
                                    DOC_12C=DOC_12Cinit3, DOC_13C=DOC_13Cinit3,
                                    Cres_12C=0, Cres_13C=0, CO2_12C=0, CO2_13C=0),
                                parms=mparun, db_modelun, times=seq(0,50)))
    yhat_4un<-as.data.frame(ode(y=c(R_12C=R_12Cinit4, R_13C=R_13Cinit4,
                                    S_12C=S_12Cinit4, S_13C=S_13Cinit4,
                                    DOC_12C=DOC_12Cinit4, DOC_13C=DOC_13Cinit4,
                                    Cres_12C=0, Cres_13C=0, CO2_12C=0, CO2_13C=0),
                                parms=mparun, db_modelun, times=seq(0,50)))
    #2. labelled treatments
    yhat_1l<-as.data.frame(ode(y=c(R_12C=R_12Cinit1, R_13C=R_13Cinit1,
                                   S_12C=S_12Cinit1, S_13C=S_13Cinit1,
                                   G_12C=G_12Cinit1, G_13C=G_13Cinit1,
                                   DOC_12C=DOC_12Cinit1, DOC_13C=DOC_13Cinit1,
                                   Cres_12C=0, Cres_13C=0, CO2_12C=0, CO2_13C=0),
                               parms=mparl, db_modell, times=seq(0,50)))
    yhat_2l<-as.data.frame(ode(y=c(R_12C=R_12Cinit2, R_13C=R_13Cinit2,
                                   S_12C=S_12Cinit2, S_13C=S_13Cinit2,
                                   G_12C=G_12Cinit2, G_13C=G_13Cinit2,
                                   DOC_12C=DOC_12Cinit2, DOC_13C=DOC_13Cinit2,
                                   Cres_12C=0, Cres_13C=0, CO2_12C=0, CO2_13C=0),
                               parms=mparl, db_modell, times=seq(0,50)))
    yhat_3l<-as.data.frame(ode(y=c(R_12C=R_12Cinit3, R_13C=R_13Cinit3,
                                   S_12C=S_12Cinit3, S_13C=S_13Cinit3,
                                   G_12C=G_12Cinit3, G_13C=G_13Cinit3,
                                   DOC_12C=DOC_12Cinit3, DOC_13C=DOC_13Cinit3,
                                   Cres_12C=0, Cres_13C=0, CO2_12C=0, CO2_13C=0),
                               parms=mparl, db_modell, times=seq(0,50)))
    yhat_4l<-as.data.frame(ode(y=c(R_12C=R_12Cinit4, R_13C=R_13Cinit4,
                                   S_12C=S_12Cinit4, S_13C=S_13Cinit4,
                                   G_12C=G_12Cinit4, G_13C=G_13Cinit4,
                                   DOC_12C=DOC_12Cinit4, DOC_13C=DOC_13Cinit4,
                                   Cres_12C=0, Cres_13C=0, CO2_12C=0, CO2_13C=0),
                               parms=mparl, db_modell, times=seq(0,50)))

    #combine all unlabelled simulations together
    yhat_allun<-rbind(yhat_1un, yhat_2un, yhat_3un, yhat_4un)
    #order them by time
    yhat_allun<-yhat_allun[order(yhat_allun$time), ]

    #convert the simulated dataset into long format data frame
    Yhatun<-melt(yhat_allun, id.vars = "time")

    #combine all labelled simulations together
    yhat_alll<-rbind(yhat_1l, yhat_2l, yhat_3l, yhat_4l)
    #order them by time
    yhat_alll<-yhat_alll[order(yhat_alll$time), ]

    #convert the simulated dataset into long format data frame
    Yhatl<-melt(yhat_alll, id.vars = "time")

    #merge all simulations together
    Yhat<-rbind(Yhatun, Yhatl)

    #add labels
    Yhat$Plesne<-c(rep(as.numeric(unlabelled$Plesne[1]), times = nrow(Yhatun)),
                   rep(as.numeric(unlabelled$Plesne[1]), times = nrow(Yhatl)))
    Yhat$Certovo<-c(rep(as.numeric(unlabelled$Certovo[1]), times = nrow(Yhatun)),
                    rep(as.numeric(unlabelled$Certovo[1]), times = nrow(Yhatl)))
    Yhat$horizon<-c(rep(unlabelled$horizon[1], times = nrow(Yhatun)),
                    rep(unlabelled$horizon[1], times = nrow(Yhatl)))
    Yhat$Treatment<-c(rep("Unlabelled", times = nrow(Yhatun)),
                      rep("Labelled", times = nrow(Yhatl)))
    return(Yhat)
  }

  #First simulation is run
  simul_out<-simul(x=as.numeric(pars_all[pars_all$id==1, c(1:11)]),
                   unlabelled = dataset[dataset$id==1, ],
                   labelled = dataset[dataset$id==11, ])

  #All other simulations are run and appended to existing "simul_out"
  for(i in 2:5){

    simul_out<-rbind(simul_out, simul(x=as.numeric(pars_all[pars_all$id==i, c(1:11)]),
                                      unlabelled = dataset[dataset$id==i, ],
                                      labelled = dataset[dataset$id==i+10, ]))
  }


  #All important calulations are stored in the "f_out" list and returned
  #1. best model parameter
  #2. goodness of correspondence across all treatments
  #3. goodness of correspondence for unlabelled and labelled treatments separately
  #4. goodness of correspondence for each horizon separately
  #5. goodness of correspondence for each horizon and unlabelled/labelled treatment separately
  #6. simulated and measured data
  #7. simulations
  #8. raw data with complete information about model parameters estimations
  #9. detail information about the goodness of correspondence between simulations
     #and measurements for each treatment and each variable

  f_out<-list(pars_all=pars_all,
              goodness1=Gfit_treatments,
              goodness2=Gfit_labels,
              #goodness3=Gfit_horizon,
              #goodness4=Gfit_each,
              Yhat_treatments=Yhat_treatments,
              simul=simul_out,
              pars_raw=epar_out,
              goodness_raw=good_out)

  return(f_out)
}
