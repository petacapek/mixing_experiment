DEB_picina12<-function(parameters, data){
  
  #define the model
  deriv<-function(time, state, pars){
    
    with(as.list(c(state, pars)),{
      
      R<-matrix(state[1:5], ncol = 1)
      S<-matrix(state[6:10], ncol = 1)
      DOC<-matrix(state[11:15], ncol = 1)
      
      #C uptake
      Cu=Vmax*S*DOC/(Km+DOC)
      
      #udrzovani
      m=S*m0
      
      #mobilizace rezerv
      an=f*R-m
      
      #respirace
      r=m+pmax(an*(1-Yu), 0)+(1-Ac)*Cu
      
      #microbialni uhlik je neco z rezerv a neco ze struktur
      Cmic=fr*R+fs*S
      
      dR<-Ac*Cu-f*R
      dS<-pmax(an*Yu,0)+pmin(0, an)
      dDOC<--Cu
      
      
      return(list(c(dR, dS, dDOC), Cmic=Cmic, r=r))
      
    })
  }
  
  
  #define names of parameters
  parnames<-c("Ac", "Vmax", "Km", "m0", "f", "Yu", "fr", "fs", "Rinit")
  
    #defining goodness of fit function 
    rsq_ode<-function(x, odeset){
      
      par<-x[1:length(parnames)]
      
      names(par)<-parnames
      
      Cmicinit=93.6
      DOCinit=315
      
      #initial R
      Rinit<-par[["Rinit"]]
      
      #initial S is then
      Sinit<-(Cmicinit-Rinit*par[["fr"]])/par[["fs"]]
     
      
      #first, pars dependent output from ode is matched with measured values
      yhat_all<-as.data.frame(ode(y=c(R1=Rinit, R2=Rinit, R3=Rinit, R4=Rinit, R5=Rinit,
                                      S1=Sinit, S2=Sinit, S3=Sinit, S4=Sinit, S5=Sinit,
                                      DOC1=315, DOC2=315, DOC3=315, DOC4=315, DOC5=315), 
                                  parms=par[1:8], 
                                  deriv, 
                                  times=c(0, 8, 16, 24, 48, 72)))
      #select time and the measured variables 
      yhat<-select(yhat_all, c("time", "Cmic1", "Cmic2", "Cmic3","Cmic4", "Cmic5",
                               "DOC1", "DOC2", "DOC3", "DOC4", "DOC5",
                               "r1", "r2", "r3", "r4", "r5"))
      
      #reformat to long format data frame
      Yhat<-melt(yhat, id.vars = "time")
      
      
      #add the measured data to a data frame
      Yhat$obs<-c(as.numeric(odeset[c(1,6,11,16, 21, 26), c("Cmic")]),
                  as.numeric(odeset[c(2,7,12,17, 22, 27), c("Cmic")]),
                  as.numeric(odeset[c(3,8,13,18, 23, 28), c("Cmic")]),
                  as.numeric(odeset[c(4,9,14,19, 24, 29), c("Cmic")]),
                  as.numeric(odeset[c(5,10,15,20, 25, 30), c("Cmic")]),
                  as.numeric(odeset[c(1,6,11,16, 21, 26), c("DOC")]),
                  as.numeric(odeset[c(2,7,12,17, 22, 27), c("DOC")]),
                  as.numeric(odeset[c(3,8,13,18, 23, 28), c("DOC")]),
                  as.numeric(odeset[c(4,9,14,19, 24, 29), c("DOC")]),
                  as.numeric(odeset[c(5,10,15,20, 25, 30), c("DOC")]),
                  as.numeric(odeset[c(1,6,11,16, 21, 26), c("r")]),
                  as.numeric(odeset[c(2,7,12,17, 22, 27), c("r")]),
                  as.numeric(odeset[c(3,8,13,18, 23, 28), c("r")]),
                  as.numeric(odeset[c(4,9,14,19, 24, 29), c("r")]),
                  as.numeric(odeset[c(5,10,15,20, 25, 30), c("r")]))
      
      Yhat$variable2<-substr(Yhat$variable, 1,1)
      
      
      #rsquared calculation for each variable
      Gfit<-Yhat %>% group_by(variable2) %>% summarise(SSres=sum(((obs-value)^2), na.rm = T), 
                                                      SStot=sum(((obs-mean(obs, na.rm = T))^2), na.rm = T),
                                                      ll=-sum(((obs-value)^2), na.rm = T)/2/(sd(obs, na.rm = T)^2))
      Gfit$R2<-with(Gfit, 1-SSres/SStot)
      Gfit$N<-length(x)
      Gfit$AIC<-with(Gfit, 2*N-2*ll)
      Gfit$Soil<-odeset$Soil[1]
      
      rsq_out<-list(Yhat=Yhat, Gfit=Gfit)
      
      return(rsq_out)
      
    }
    
    
    fit<-rsq_ode(x=parameters, odeset = data)
    
  return(fit)
}