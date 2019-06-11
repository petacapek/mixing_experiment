enzymes_calc<-function(filepath, Sconc){
  
  require(openxlsx)
  require(reshape)
  require(minpack.lm)
  
  #reading the file data 
  #there are three measurments in time, each is loaded separately
  
  #################################################################################
  #first the script has to read sample identifier
  IDs<-as.vector(read.xlsx(xlsxFile=filepath, 
                 sheet=2, 
                 rows = c(3:6),colNames = F, 
                 cols = c(16))[,1])
  
  #times t1 and t2
  t<-as.vector(read.xlsx(xlsxFile=filepath, 
                           sheet=2, 
                           rows = c(24,34),colNames = F, 
                           cols = c(19))[,1])
  #and dry weights
  DW<-as.vector(read.xlsx(xlsxFile=filepath, 
                          sheet=1, 
                          rows = c(3:6),colNames = F, 
                          cols = c(4))[,1])
  
  #time 0
  dat_t0<-as.data.frame(t(read.xlsx(xlsxFile=filepath, sheet=2, 
                                    rows = c(14:18),colNames = F, 
                    cols = c(2:13))))
  
  ##specify the enzymes, time of measurement, sample identification and substrate concentration
  ##we need long format data
  colnames(dat_t0)<-rep(c("glu", "cel", "phos", "leu", "chit"), times=ncol(dat_t0)/5)
  dat_t0$time<-rep(0.5, times=nrow(dat_t0))
  #what sample it is
  dat_t0$Sample<-rep(IDs, each=3)
  #and what is the DW of the sample
  dat_t0$DW<-rep(DW, each=3)
  #long format
  Dat_t0<-melt(dat_t0, id.vars = c("time", "Sample", "DW"))
  colnames(Dat_t0)<-c("time", "Sample", "DW", "E", "measure")
  
  #concentration of added substrate - this depends on the sample
  #therefore for loop has to be used
  #First, empty column (conc) is created in the data frame
  Dat_t0$conc<-vector("numeric", length = nrow(Dat_t0))
  elist<-c("glu", "cel", "phos", "leu", "chit")
  
  #than, for loop is running
  for(i in IDs){
    
    for(n in elist){
      
      Dat_t0[(Dat_t0$Sample==i & Dat_t0$E==n), "conc"]<-Sconc[(Sconc$Sample==i & Sconc$E==n), "conc"]
      
      
    }
    
    
  }
    
  #define the substrate
  Dat_t0$Substrate<-ifelse(Dat_t0$E=="chit", "amc", "muf")
  
  #time 1
  dat_t1<-as.data.frame(t(read.xlsx(xlsxFile=filepath, sheet=2, colNames=F, rows = c(24:28),
                      cols = c(2:13))))
  
  ##specify the enzymes (colnames), time of measurement and substrate concentration
  colnames(dat_t1)<-c("glu", "cel", "phos", "leu", "chit")
  dat_t1$time<-rep(t[1], times=nrow(dat_t1))
  #what sample it is
  dat_t1$Sample<-dat_t0$Sample
  #and what is the DW of the sample
  dat_t1$DW<-rep(DW, each=3)
  #long format
  Dat_t1<-melt(dat_t1, id.vars = c("time", "Sample", "DW"))
  colnames(Dat_t1)<-c("time", "Sample", "DW", "E", "measure")
  Dat_t1$conc<-Dat_t0$conc
  Dat_t1$Substrate<-Dat_t0$Substrate
  
  #time 2
  dat_t2<-as.data.frame(t(read.xlsx(xlsxFile=filepath, sheet=2, colNames = F, rows = c(34:38),
                      cols = c(2:13))))
  
  ##specify the enzymes (colnames), time of measurement and substrate concentration
  colnames(dat_t2)<-c("glu", "cel", "phos", "leu", "chit")
  dat_t2$time<-rep(t[2], times=nrow(dat_t2))
  #what sample it is
  dat_t2$Sample<-dat_t0$Sample
  #and what is the DW of the sample
  dat_t2$DW<-rep(DW, each=3)
  #long format
  Dat_t2<-melt(dat_t2, id.vars = c("time", "Sample", "DW"))
  colnames(Dat_t2)<-c("time", "Sample", "DW", "E", "measure")
  Dat_t2$conc<-Dat_t0$conc
  Dat_t2$Substrate<-Dat_t0$Substrate
  
  #################################################################################
  
  #after reading the data itself, calibration data are loaded
  #every sample (three time points together) has own calibration data
  #there are two calibration data types
  #first are calibration data of MUF product (glu. cell, phos and leu) (muf is released when functional group is detached by enzyme action - muf:functional group = 1:1)
    #this data will be labeled as muf_cal0, muf_cal1 and muf_cal2
  #second are calibration data of AMC product (chit)
    #this data will be labeled as amc_cal0, amc_cal1 and amc_cal2
  
  ##################################################################################
  muf_cal0<-as.data.frame(t(read.xlsx(xlsxFile=filepath, sheet=2, rows = c(19,20),
                                     cols = c(2:13), colNames = F)))
  
  muf_cal0$Sample<-rep(IDs, each=3)
  muf_cal0<-melt(muf_cal0, id.vars = 'Sample')
  muf_cal0<-muf_cal0[order(muf_cal0$Sample),]
  muf_cal0$conc<-rep(c(rep(0, 3), 1, 5, 10), times=length(IDs))
  
  
  ###
  muf_cal1<-as.data.frame(t(read.xlsx(xlsxFile=filepath, sheet=2, rows = c(29,30),
                                      cols = c(2:13), colNames = F)))
  
  muf_cal1$Sample<-rep(IDs, each=3)
  muf_cal1<-melt(muf_cal1, id.vars = 'Sample')
  muf_cal1<-muf_cal1[order(muf_cal1$Sample),]
  muf_cal1$conc<-rep(c(rep(0, 3), 1, 5, 10), times=length(IDs))
  
  ###
  muf_cal2<-as.data.frame(t(read.xlsx(xlsxFile=filepath, sheet=2, rows = c(39,40),
                                      cols = c(2:13), colNames = F)))
  
  muf_cal2$Sample<-rep(IDs, each=3)
  muf_cal2<-melt(muf_cal2, id.vars = 'Sample')
  muf_cal2<-muf_cal2[order(muf_cal2$Sample),]
  muf_cal2$conc<-rep(c(rep(0, 3), 1, 5, 10), times=length(IDs))
  
  ####################
  amc_cal0<-as.data.frame(t(read.xlsx(xlsxFile=filepath, sheet=2, rows = c(19,21),
                                      cols = c(2:13), colNames = F)))
  
  amc_cal0$Sample<-rep(IDs, each=3)
  amc_cal0<-melt(amc_cal0, id.vars = 'Sample')
  amc_cal0<-amc_cal0[order(amc_cal0$Sample),]
  amc_cal0$conc<-rep(c(rep(0, 3), 1, 5, 10), times=length(IDs))
  
  
  ###
  amc_cal1<-as.data.frame(t(read.xlsx(xlsxFile=filepath, sheet=2, rows = c(29, 31),
                                      cols = c(2:13), colNames = F)))
  
  amc_cal1$Sample<-rep(IDs, each=3)
  amc_cal1<-melt(amc_cal1, id.vars = 'Sample')
  amc_cal1<-amc_cal1[order(amc_cal1$Sample),]
  amc_cal1$conc<-rep(c(rep(0, 3), 1, 5, 10), times=length(IDs))
  
  
  ###
  amc_cal2<-as.data.frame(t(read.xlsx(xlsxFile=filepath, sheet=2, rows = c(39,41),
                                      cols = c(2:13), colNames = F)))
  
  amc_cal2$Sample<-rep(IDs, each=3)
  amc_cal2<-melt(amc_cal2, id.vars = 'Sample')
  amc_cal2<-amc_cal2[order(amc_cal2$Sample),]
  amc_cal2$conc<-rep(c(rep(0, 3), 1, 5, 10), times=length(IDs))
  #######################################################################################
  #Based on the calibration data, linear regression models are calculated for both products (muf and amc)
  #for each sample
  #we need to know the model coeficients only
  #they are abeled as muf_coef0 - muf_coef2 and amc_coef0 - amc_coef2
  #these models will be used to calculate the product concentration later on
  ########################################################################################
  
  #muf substrate
  #Time 0
  #data drame with all calibration data is created
  muf_coef0<-data.frame(Sample=rep(IDs), Intercept=vector("numeric", length = length(IDs)),
                       Slope=vector("numeric", length = length(IDs)))
  
  #model coefficients are estimated for each sample 
  for(i in IDs){
    
    #Intercept
    muf_coef0[muf_coef0$Sample==i, "Intercept"]<-coef(lm(conc~value, data = muf_cal0[muf_cal0$Sample==i,]))[1]
    #Slope
    muf_coef0[muf_coef0$Sample==i, "Slope"]<-coef(lm(conc~value, data = muf_cal0[muf_cal0$Sample==i,]))[2]
  
  }
  
  #Time 1
  #data drame with all calibration data is created
  muf_coef1<-data.frame(Sample=rep(IDs), Intercept=vector("numeric", length = length(IDs)),
                        Slope=vector("numeric", length = length(IDs)))
  
  #model coefficients are estimated for each sample 
  for(i in IDs){
    
    #Intercept
    muf_coef1[muf_coef1$Sample==i, "Intercept"]<-coef(lm(conc~value, data = muf_cal1[muf_cal1$Sample==i,]))[1]
    #Slope
    muf_coef1[muf_coef1$Sample==i, "Slope"]<-coef(lm(conc~value, data = muf_cal1[muf_cal1$Sample==i,]))[2]
    
  }
  
  #Time 2
  #data drame with all calibration data is created
  muf_coef2<-data.frame(Sample=rep(IDs), Intercept=vector("numeric", length = length(IDs)),
                        Slope=vector("numeric", length = length(IDs)))
  
  #model coefficients are estimated for each sample 
  for(i in IDs){
    
    #Intercept
    muf_coef2[muf_coef2$Sample==i, "Intercept"]<-coef(lm(conc~value, data = muf_cal2[muf_cal2$Sample==i,]))[1]
    #Slope
    muf_coef2[muf_coef2$Sample==i, "Slope"]<-coef(lm(conc~value, data = muf_cal2[muf_cal2$Sample==i,]))[2]
    
  }
    
  #the same is done with calibration data for amc substrate
  #Time 0
  #data drame with all calibration data is created
  amc_coef0<-data.frame(Sample=rep(IDs), Intercept=vector("numeric", length = length(IDs)),
                        Slope=vector("numeric", length = length(IDs)))
  
  #model coefficients are estimated for each sample 
  for(i in IDs){
    
    #Intercept
    amc_coef0[amc_coef0$Sample==i, "Intercept"]<-coef(lm(conc~value, data = amc_cal0[amc_cal0$Sample==i,]))[1]
    #Slope
    amc_coef0[amc_coef0$Sample==i, "Slope"]<-coef(lm(conc~value, data = amc_cal0[amc_cal0$Sample==i,]))[2]
    
  }
  
  #Time 1
  #data drame with all calibration data is created
  amc_coef1<-data.frame(Sample=rep(IDs), Intercept=vector("numeric", length = length(IDs)),
                        Slope=vector("numeric", length = length(IDs)))
  
  #model coefficients are estimated for each sample 
  for(i in IDs){
    
    #Intercept
    amc_coef1[amc_coef1$Sample==i, "Intercept"]<-coef(lm(conc~value, data = amc_cal1[amc_cal1$Sample==i,]))[1]
    #Slope
    amc_coef1[amc_coef1$Sample==i, "Slope"]<-coef(lm(conc~value, data = amc_cal1[amc_cal1$Sample==i,]))[2]
    
  }
  
  #Time 2
  #data drame with all calibration data is created
  amc_coef2<-data.frame(Sample=rep(IDs), Intercept=vector("numeric", length = length(IDs)),
                        Slope=vector("numeric", length = length(IDs)))
  
  #model coefficients are estimated for each sample 
  for(i in IDs){
    
    #Intercept
    amc_coef2[amc_coef2$Sample==i, "Intercept"]<-coef(lm(conc~value, data = amc_cal2[amc_cal2$Sample==i,]))[1]
    #Slope
    amc_coef2[amc_coef2$Sample==i, "Slope"]<-coef(lm(conc~value, data = amc_cal2[amc_cal2$Sample==i,]))[2]
    
  }
  
  
  
  #######################################################################################
  
  #############################################################################################
  #now, product concentration is calculated using stored linear regression models
  #notice that the linearity is required
  #the final concentration is umol of product per gram of dry soil ([(250/200] - dilution when suspension is mixed with substrate solution),
  #[(200/1e6) - from conc in umol/l to amount umol in 200 ul], [(50e3/200) - scale to 50 ml of total suspension],
  #[(0.5*DW) - to gram of dry soil]
  #substrate concentration is also recalculated per gram of dry soil
  #all negative values are again replaced by 0
  #calculation is done separately for each Sample
  #################################################################################################
  #time point 0
  #new data frame with new column is created
  final_t0<-Dat_t0
  final_t0$product<-vector("numeric", length = nrow(final_t0))
  
  #calculation is done for each sample
  #for muf substrate
  for(i in IDs){
    
    final_t0[(final_t0$Sample==i & final_t0$E!="chit" ), "product"]<-(final_t0[(final_t0$Sample==i & final_t0$E!="chit" ), "measure"]*
      as.numeric(muf_coef0[muf_coef0$Sample==i, "Slope"])+as.numeric(muf_coef0[muf_coef0$Sample==i, "Intercept"]))*
      (250/200)*(200/1e6)*(50e3/200)/(0.5*final_t0[(final_t0$Sample==i & final_t0$E!="chit" ), "DW"])
  }
  #and for amc
  for(i in IDs){
    
    final_t0[(final_t0$Sample==i & final_t0$E=="chit" ), "product"]<-(final_t0[(final_t0$Sample==i & final_t0$E=="chit" ), "measure"]*
                                                                        as.numeric(amc_coef0[amc_coef0$Sample==i, "Slope"])+as.numeric(amc_coef0[amc_coef0$Sample==i, "Intercept"]))*
      (250/200)*(200/1e6)*(50e3/200)/(0.5*final_t0[(final_t0$Sample==i & final_t0$E=="chit" ), "DW"])
  }
  
  #negative values are replaced by 0
  final_t0[final_t0$product<0, "product"]<-c(0)
  
  #further, concentration of added substrate needs to be recalculated per g of dry weight as well
  final_t0$conc<-with(final_t0, conc*(200/1e6)*(50e3/200)/(0.5*DW))
  
  
  
  #time point 1
  #new data frame with new column is created
  final_t1<-Dat_t1
  final_t1$product<-vector("numeric", length = nrow(final_t1))
  
  #calculation is done for each sample
  #for muf substrate
  for(i in IDs){
    
    final_t1[(final_t1$Sample==i & final_t1$E!="chit" ), "product"]<-(final_t1[(final_t1$Sample==i & final_t1$E!="chit" ), "measure"]*
                                                                        as.numeric(muf_coef1[muf_coef1$Sample==i, "Slope"])+as.numeric(muf_coef1[muf_coef1$Sample==i, "Intercept"]))*
      (250/200)*(200/1e6)*(50e3/200)/(0.5*final_t1[(final_t1$Sample==i & final_t1$E!="chit" ), "DW"])
  }
  #and for amc
  for(i in IDs){
    
    final_t1[(final_t1$Sample==i & final_t1$E=="chit" ), "product"]<-(final_t1[(final_t1$Sample==i & final_t1$E=="chit" ), "measure"]*
                                                                        as.numeric(amc_coef1[amc_coef1$Sample==i, "Slope"])+as.numeric(amc_coef1[amc_coef1$Sample==i, "Intercept"]))*
      (250/200)*(200/1e6)*(50e3/200)/(0.5*final_t1[(final_t1$Sample==i & final_t1$E=="chit" ), "DW"])
  }
  
  #negative values are replaced by 0
  final_t1[final_t1$product<0, "product"]<-c(0)
  
  #further, concentration of added substrate needs to be recalculated per g of dry weight as well
  final_t1$conc<-with(final_t1, conc*(200/1e6)*(50e3/200)/(0.5*DW))
  
  
  #time point 2
  #new data frame with new column is created
  final_t2<-Dat_t2
  final_t2$product<-vector("numeric", length = nrow(final_t2))
  
  #calculation is done for each sample
  #for muf substrate
  for(i in IDs){
    
    final_t2[(final_t2$Sample==i & final_t2$E!="chit" ), "product"]<-(final_t2[(final_t2$Sample==i & final_t2$E!="chit" ), "measure"]*
                                                                        as.numeric(muf_coef2[muf_coef2$Sample==i, "Slope"])+as.numeric(muf_coef2[muf_coef2$Sample==i, "Intercept"]))*
      (250/200)*(200/1e6)*(50e3/200)/(0.5*final_t2[(final_t2$Sample==i & final_t2$E!="chit" ), "DW"])
  }
  #and for amc
  for(i in IDs){
    
    final_t2[(final_t2$Sample==i & final_t2$E=="chit" ), "product"]<-(final_t2[(final_t2$Sample==i & final_t2$E=="chit" ), "measure"]*
                                                                        as.numeric(amc_coef2[amc_coef2$Sample==i, "Slope"])+as.numeric(amc_coef2[amc_coef2$Sample==i, "Intercept"]))*
      (250/200)*(200/1e6)*(50e3/200)/(0.5*final_t2[(final_t2$Sample==i & final_t2$E=="chit" ), "DW"])
  }
  
  #negative values are replaced by 0
  final_t2[final_t2$product<0, "product"]<-c(0)
  
  #further, concentration of added substrate needs to be recalculated per g of dry weight as well
  final_t2$conc<-with(final_t2, conc*(200/1e6)*(50e3/200)/(0.5*DW))
  
  ######################################################################################################
  #the last calculation is the enzymatic activity
  #the calculated values have following units - Vmax - µmol/h/µmol(Substrate), Km - µmols
  
  #creat new data frame with all variables I need to know for the calculation
  activity<-rbind(final_t0[,c("Sample", "E", "time", "product", "conc")],
                  final_t1[,c("Sample", "E", "time", "product", "conc")],
                  final_t2[,c("Sample", "E", "time", "product", "conc")])
  #data transformation for linear regression
  activity$y<-with(activity, product/time)
  activity$x<-with(activity, log(1-product/conc)/time)
  
  
  #enzyme activity is calculated for each sample and each enzyme
  #summary table is created first
  res<-data.frame(Sample=rep(IDs, each=5),
                  E=rep(unique(activity$E, times=length(IDs))))
  res$Vmax<-vector("numeric", length = nrow(res))
  res$Km<-vector("numeric", length = nrow(res))              
  
  #calculation
  for(i in IDs){
    
    for(n in unique(activity$E)){
      
      res[(res$Sample==i & res$E==n), "Vmax"]<-lmodel2(y~x, data=activity[(activity$Sample==i & activity$E==n), ], 
                                                       nperm = 999)$regression.results[2,2]
      res[(res$Sample==i & res$E==n), "Km"]<-(-lmodel2(y~x, data=activity[(activity$Sample==i & activity$E==n), ],
                                                     nperm = 999)$regression.results[2,3])
    }
  }
  
  #transform to short format
  res_finala<-reshape(res[, c("Sample", "E", "Vmax")], idvar = "Sample", timevar = "E", direction = "wide")
  colnames(res_finala)<-c("Sample", "gluVmax", "celVmax", "phosVmax", "leuVmax", "chitVmax")
  res_finalb<-reshape(res[, c("Sample", "E", "Km")], idvar = "Sample", timevar = "E", direction = "wide")
  colnames(res_finalb)<-c("Sample", "gluKm", "celKm", "phosKm", "leuKm", "chitKm")
  
  res_final<-cbind(res_finala, res_finalb[,-1])
  return(res_final)
}

