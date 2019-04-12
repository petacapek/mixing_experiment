AWB_picina12<-function(parameters, data){
  
  #define the model
  deriv<-function(time, state, pars){
    
    with(as.list(c(state, pars)),{
      
      Cb<-matrix(state[1:5], ncol = 1)
      DOC<-matrix(state[6:10], ncol = 1)
      
      #C uptake
      Cu=Vmax*Cb*DOC/(Km+DOC)
      
      dCb<-CUE*Cu-k*Cb
      dDOC<--Cu+a*k*Cb
      
      return(list(c(dCb, dDOC), r=(1-CUE)*Cu, Cmic=Cb*kec))
      
    })
  }
  
  
  #define names of parameters
  parnames<-c("Vmax", "Km", "CUE", "k", "a", "kec")
  
    #defining goodness of fit function 
    rsq_ode<-function(x, odeset){
      
      par<-x[1:length(parnames)]
      
      names(par)<-parnames
      
      Cmicinit=93.6/par[["kec"]]
      DOCinit=315
      
      #first, pars dependent output from ode is matched with measured values
      yhat_all<-as.data.frame(ode(y=c(Cb1=Cmicinit, Cb2=Cmicinit, Cb3=Cmicinit, Cb4=Cmicinit, Cb5=Cmicinit,
                                      DOC1=315, DOC2=315, DOC3=315, DOC4=315, DOC5=315), 
                                  parms=par, 
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