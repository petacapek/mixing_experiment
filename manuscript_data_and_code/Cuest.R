Cuest<-function(data){
  #Define microbial model
  cfun<-function(time, state, pars){
    
    with(as.list(c(state, pars)),{
      
      dCmic<--a*k2*Cmic+CUE*k1*Cmic*C
      dC<-a*k2*Cmic-CUE*k1*Cmic*C
      dCO2<-(1-CUE)*k1*Cmic*C
      
      return(list(c(dCmic, dC, dCO2)))
      
    })
  }
  #Define estim function
  estim<-function(d){
    cost<-function(x){
      par<-x
      names(par)<-c("k1", "k2", "a", "CUE")
      yhat <- as.data.frame(ode(y=c(Cmic=d$Cmic[1], C=d$DOC[1], CO2=0), parms=par, 
                                cfun, times=c(0, 48)))
      Yhat<-melt(yhat, id.vars="time")
      Yhat$obs<-c(d$Cmic[1], d$Cmic[2], d$DOC[1], d$DOC[2], 0, d$CCO2[2])
      
      yweights<-Yhat %>% group_by(variable) %>% summarise(weights=mean(value, na.rm = T))
      
      #match with the Yhat data frame
      Yhat<-merge(Yhat, yweights, by.x=c("variable"), by.y=c("variable"))
      
      #now, the root mean square error is calculated
      NRMSE<-as.numeric(Yhat %>% group_by(variable) %>% summarise(NRMSE=sum((((value-obs)/mean(weights))^2), na.rm = T)) %>%
                          summarise(NRMSE=sum(NRMSE)))
      
      return(NRMSE)
    }
    good<-function(x){
      par<-x
      names(par)<-c("k1", "k2", "a", "CUE")
      yhat <- as.data.frame(ode(y=c(Cmic=d$Cmic[1], C=d$DOC[1], CO2=0), parms=par, 
                                cfun, times=c(0, 48)))
      Yhat<-melt(yhat, id.vars="time")
      Yhat$obs<-c(d$Cmic[1], d$Cmic[2], d$DOC[1], d$DOC[2], 0, d$CCO2[2])
      
      Yhat$Plesne<-rep(d$Plesne[1], times=nrow(Yhat))
      Yhat$Horizon<-rep(d$Horizon[1], times=nrow(Yhat))
      
      return(Yhat)
    }
    
    pt<-modMCMC(f=cost, p=c(k1=0.1, k2=0.1, a=0.5, CUE=0.5),
                lower=c(k1=1e-5, k2=1e-5, a=0, CUE=0),
                upper=c(k1=10, k2=10, a=1, CUE=0.8), niter=10000)
    #lower and upper limits for parameters are extracted
    pl<-summary(pt)["min",]
    pu<-summary(pt)["max",]
    
    #these limits are used to find global optimum by DEoptim
    opt_par<-DEoptim(fn=cost, lower=pl, upper=pu,
                     control = c(itermax = 10000, steptol = 50, reltol = 1e-8,
                                 trace=FALSE, strategy=3, NP=250))
    
    #goodness of fit
    fit<-good(opt_par$optim$bestmem)
    
    #best parameters
    p<-opt_par$optim$bestmem
    names(p)<-c("k1", "k2", "a", "CUE")
    
    #return list with opt_par and par_prof
    estim_out<-list(pars=p, fit=fit)
    
    return(estim_out)
  }
  
  res<-foreach(i=1:40, .combine=list, .multicombine = TRUE,
               .packages=c("FME", "dplyr", "DEoptim", "reshape")) %dopar% {
                 
                 estim(d=data[c(i, i+80),])
                 
               }
  
  return(res)
}