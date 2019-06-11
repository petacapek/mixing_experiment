############################################Biomass data##########################################
##Cmic
mich2 %>% group_by(Plesne, znaceni, horizont, Legend) %>% summarize(y.sd=sd(Cmic), y=mean(Cmic)) %>%
  ggplot(aes(factor(Plesne), y))+geom_point(pch=21, cex=6, aes(fill=Legend))+
  facet_grid(znaceni~horizont, scales="free")+theme_min+theme(legend.position = c(0.2, 0.88))+
  geom_errorbar(aes(ymin=y-y.sd, ymax=y+y.sd), width=0.1, lwd=0.5)+
  ylab(expression(paste("Microbial biomass carbon ( ", mu, "mol ",g^{-1}, ")" )))+
  xlab("Plesne : Certovo mixing ratio")+scale_fill_manual(values = c("white", "grey"))

##Nmic
mich2 %>% group_by(Plesne, znaceni, horizont, Legend) %>% summarize(y.sd=sd(Nmic), y=mean(Nmic)) %>%
ggplot(aes(factor(Plesne), y))+geom_point(pch=21, cex=6, aes(fill=Legend))+
  geom_errorbar(aes(ymin=y-y.sd, ymax=y+y.sd), width=0.1, lwd=0.5)+
  facet_grid(znaceni~horizont, scales="free")+theme_min+theme(legend.position = c(0.2, 0.88))+
  ylab(expression(paste("Microbial biomass nitrogen ( ", mu, "mol ",g^{-1}, ")" )))+
  xlab("Plesne : Certovo mixing ratio")+scale_fill_manual(values = c("white", "grey"))

##Pmic
mich2 %>% group_by(Plesne, znaceni, horizont, Legend) %>% summarize(y.sd=sd(Pmic), y=mean(Pmic)) %>%
  ggplot(aes(factor(Plesne), y))+geom_point(pch=21, cex=6, aes(fill=Legend))+
  geom_errorbar(aes(ymin=y-y.sd, ymax=y+y.sd), width=0.1, lwd=0.5)+
  facet_grid(znaceni~horizont, scales="free")+theme_min+theme(legend.position = c(0.17, 0.88))+
  ylab(expression(paste("Microbial biomass phosphorus ( ", mu, "mol ",g^{-1}, ")" )))+
  xlab("Plesne : Certovo mixing ratio")+scale_fill_manual(values = c("white", "grey"))

##C/N mic
mich2 %>% group_by(Plesne, znaceni, horizont, Legend) %>% summarize(y.sd=sd(Cmic/Nmic), y=mean(Cmic/Nmic)) %>%
ggplot(aes(factor(Plesne), y))+geom_point(pch=21, cex=6, aes(fill=Legend))+
  geom_errorbar(aes(ymin=y-y.sd, ymax=y+y.sd), width=0.1, lwd=0.5)+
  facet_grid(znaceni~horizont, scales="free")+theme_min+theme(legend.position = c(0.33, 0.9),
                                                              legend.key.size = unit(0.3, "in"))+
  ylab(expression(paste("Microbial C/N (mol/mol)" )))+
  xlab("Plesne : Certovo mixing ratio")+scale_fill_manual(values = c("white", "grey"))

##C/P mic
mich2 %>% group_by(Plesne, znaceni, horizont, Legend) %>% summarize(y.sd=sd(Cmic/Pmic), y=mean(Cmic/Pmic)) %>%
  ggplot(aes(factor(Plesne), y))+geom_point(pch=21, cex=6, aes(fill=Legend))+
  geom_errorbar(aes(ymin=y-y.sd, ymax=y+y.sd), width=0.1, lwd=0.5)+
  facet_grid(znaceni~horizont)+theme_min+theme(legend.position = c(0.17, 0.88),
                                                              legend.key.size = unit(0.3, "in"))+
  ylab(expression(paste("Microbial C/P (mol/mol)" )))+
  xlab("Plesne : Certovo mixing ratio")+scale_fill_manual(values = c("white", "grey"))

##N/P mic
mich2 %>% group_by(Plesne, znaceni, horizont, Legend) %>% summarize(y.sd=sd(Nmic/Pmic), y=mean(Nmic/Pmic)) %>%
  ggplot(aes(factor(Plesne), y))+geom_point(pch=21, cex=6, aes(fill=Legend))+
  geom_errorbar(aes(ymin=y-y.sd, ymax=y+y.sd), width=0.1, lwd=0.5)+
  facet_grid(znaceni~horizont, scales="free")+theme_min+theme(legend.position = c(0.17, 0.88),
                                                              legend.key.size = unit(0.3, "in"),
                                                              legend.title = element_blank())+
  ylab(expression(paste("Microbial N/P (mol/mol)" )))+
  xlab("Plesne : Certovo mixing ratio")+scale_fill_manual(values = c("white", "grey"))

####zmeny
mb.zmeny<-mich2[1:80, c("Plesne", "Certovo", "horizont", "znaceni")]
mb.zmeny$dCmic<-mich2[81:160, "Cmic"]-mich2[1:80, "Cmic"]
mb.zmeny$dNmic<-mich2[81:160, "Nmic"]-mich2[1:80, "Nmic"]
mb.zmeny$dPmic<-mich2[81:160, "Pmic"]-mich2[1:80, "Pmic"]
mb.zmeny$O2<-mich2[81:160, "O2"]
mb.zmeny$CO2<-mich2[81:160, "CCO2c"]
mb.zmeny$DOC_all<-mich2[1:80, "DOC"]+mich2[1:80, "Gl"]

##Cmic
mb.zmeny %>% group_by(Plesne, znaceni, horizont) %>% summarize(y.sd=sd(dCmic), y=mean(dCmic)) %>%
ggplot(aes(factor(Plesne), y))+geom_point(cex=6, pch=21, aes(fill=znaceni))+
  facet_grid(.~horizont, scales="free")+theme_min+theme(legend.position = c(0.2, 0.2))+
  geom_errorbar(aes(ymin=y-y.sd, ymax=y+y.sd), width=0.1, lwd=0.5)+
  ylab(expression(paste(Delta, "Microbial biomass carbon ( ", mu, "mol ",g^{-1}, ")" )))+
  xlab("Plesne : Certovo mixing ratio")+scale_fill_manual(values = c("white", "grey"))+
  geom_hline(yintercept = 0, lwd=1.5)+theme(legend.position = c(0.15, 0.2),
                                            legend.key.size = unit(0.3, "in"),
                                            legend.title = element_blank())

mb.zmeny %>% group_by(Plesne, znaceni, horizont) %>% summarize(y.sd=sd(dCmic), y=mean(dCmic),
                                                               x.sd=sd(O2), x=mean(O2)) %>%
  ggplot(aes(x, y))+geom_point(cex=6, pch=21, aes(fill=znaceni))+
  facet_wrap(znaceni~horizont, scales="free")+theme_min+theme(legend.position = c(0.2, 0.2))+
  geom_errorbar(aes(ymin=y-y.sd, ymax=y+y.sd), width=0.1, lwd=0.5)+
  geom_errorbarh(aes(xmin=x-x.sd, xmax=x+x.sd), width=0.1, lwd=0.5)+
  ylab(expression(paste(Delta, "Microbial biomass carbon ( ", mu, "mol ",g^{-1}, ")" )))+
  xlab("Oxygen concentration (%)")+scale_fill_manual(values = c("white", "grey"))+
  geom_hline(yintercept = 0, lwd=1.5)+theme(legend.position = c(0.15, 0.2),
                                            legend.key.size = unit(0.3, "in"),
                                            legend.title = element_blank())

ggplot(mb.zmeny, aes(O2, dCmic*DOC_all))+geom_point(cex=6, pch=21, aes(fill=znaceni))+
  facet_wrap(~horizont, scales="free")+theme_min

##Nmic
mb.zmeny %>% group_by(Plesne, znaceni, horizont) %>% summarize(y.sd=sd(dNmic), y=mean(dNmic)) %>%
  ggplot(aes(factor(Plesne), y))+geom_point(cex=6, pch=21, aes(fill=znaceni))+
  facet_grid(.~horizont, scales="free")+theme_min+theme(legend.position = c(0.2, 0.2))+
  geom_errorbar(aes(ymin=y-y.sd, ymax=y+y.sd), width=0.1, lwd=0.5)+
  ylab(expression(paste(Delta, "Microbial biomass nitrogen ( ", mu, "mol ",g^{-1}, ")" )))+
  xlab("Plesne : Certovo mixing ratio")+scale_fill_manual(values = c("white", "grey"))+
  geom_hline(yintercept = 0, lwd=1.5)+theme(legend.position = c(0.15, 0.2),
                                            legend.key.size = unit(0.3, "in"),
                                            legend.title = element_blank())
##Pmic
mb.zmeny %>% group_by(Plesne, znaceni, horizont) %>% summarize(y.sd=sd(dPmic), y=mean(dPmic)) %>%
  ggplot(aes(factor(Plesne), y))+geom_point(cex=6, pch=21, aes(fill=znaceni))+
  facet_grid(.~horizont, scales="free")+theme_min+theme(legend.position = c(0.2, 0.2))+
  geom_errorbar(aes(ymin=y-y.sd, ymax=y+y.sd), width=0.1, lwd=0.5)+
  ylab(expression(paste(Delta, "Microbial biomass phosphorus ( ", mu, "mol ",g^{-1}, ")" )))+
  xlab("Plesne : Certovo mixing ratio")+scale_fill_manual(values = c("white", "grey"))+
  geom_hline(yintercept = 0, lwd=1.5)+theme(legend.position = c(0.15, 0.2),
                                            legend.key.size = unit(0.3, "in"),
                                            legend.title = element_blank())
##C/N mic
mb.zmeny[-59, ] %>% group_by(Plesne, znaceni, horizont) %>% summarize(y.sd=sd(dCmic/dNmic), y=mean(dCmic/dNmic)) %>%
  ggplot(aes(factor(Plesne), y))+geom_point(cex=6, pch=21, aes(fill=znaceni))+
  facet_grid(.~horizont, scales="free")+theme_min+theme(legend.position = c(0.2, 0.2))+
  geom_errorbar(aes(ymin=y-y.sd, ymax=y+y.sd), width=0.1, lwd=0.5)+
  ylab(expression(paste(Delta, "C/N biomass ( ", mu, "mol ",g^{-1}, ")" )))+
  xlab("Plesne : Certovo mixing ratio")+scale_fill_manual(values = c("white", "grey"))+
  geom_hline(yintercept = 0, lwd=1.5)+theme(legend.position = c(0.7, 0.2),
                                            legend.key.size = unit(0.3, "in"),
                                            legend.title = element_blank())

##C/P mic
mb.zmeny[-c(73, 31, 60, 78), ] %>% group_by(Plesne, znaceni, horizont) %>% summarize(y.sd=sd(dCmic/dPmic), y=mean(dCmic/dPmic)) %>%
  ggplot(aes(factor(Plesne), y))+geom_point(cex=6, pch=21, aes(fill=znaceni))+
  facet_grid(.~horizont, scales="free")+theme_min+theme(legend.position = c(0.2, 0.2))+
  geom_errorbar(aes(ymin=y-y.sd, ymax=y+y.sd), width=0.1, lwd=0.5)+
  ylab(expression(paste(Delta, "C/P biomass ( ", mu, "mol ",g^{-1}, ")" )))+
  xlab("Plesne : Certovo mixing ratio")+scale_fill_manual(values = c("white", "grey"))+
  geom_hline(yintercept = 0, lwd=1.5)+theme(legend.position = c(0.2, 0.8),
                                            legend.key.size = unit(0.3, "in"),
                                            legend.title = element_blank())

##N/P mic
mb.zmeny[-c(73, 60, 78), ] %>% group_by(Plesne, znaceni, horizont) %>% summarize(y.sd=sd(dNmic/dPmic), y=mean(dNmic/dPmic)) %>%
  ggplot(aes(factor(Plesne), y))+geom_point(cex=6, pch=21, aes(fill=znaceni))+
  facet_grid(.~horizont, scales="free")+theme_min+theme(legend.position = c(0.2, 0.2))+
  geom_errorbar(aes(ymin=y-y.sd, ymax=y+y.sd), width=0.1, lwd=0.5)+
  ylab(expression(paste(Delta, "N/P biomass ( ", mu, "mol ",g^{-1}, ")" )))+
  xlab("Plesne : Certovo mixing ratio")+scale_fill_manual(values = c("white", "grey"))+
  geom_hline(yintercept = 0, lwd=1.5)+theme(legend.position = c(0.2, 0.8),
                                            legend.key.size = unit(0.3, "in"),
                                            legend.title = element_blank())

#################################################Molecular data##################################################
#Bacteria
b<-t(read.xlsx(xlsxFile = c("./Chemical_and_microbial_data/michaci_otu_tables.xlsx"), sheet = 1,
             colNames = F)[-1, -c(1, 42:49)])
rownames(b)<-as.numeric(read.xlsx(xlsxFile = c("./Chemical_and_microbial_data/michaci_otu_tables.xlsx"), sheet = 1,
                                  colNames = F)[1,-c(1, 42:49)])
#Fungi
f<-t(read.xlsx(xlsxFile = c("./Chemical_and_microbial_data/michaci_otu_tables.xlsx"), sheet = 2,
             colNames = F)[-1, -c(1, 42:49)])
rownames(f)<-as.numeric(read.xlsx(xlsxFile = c("./Chemical_and_microbial_data/michaci_otu_tables.xlsx"), sheet = 2,
                                  colNames = F)[1,-c(1, 42:49)])
colnames(f)<-t(read.xlsx(xlsxFile = c("./Chemical_and_microbial_data/michaci_otu_tables.xlsx"), sheet = 2,
                                  colNames = F)[-1,1])
rownames(b)==rownames(f)

##labels
lb<-read.xlsx(xlsxFile = c("./Chemical_and_microbial_data/michaci_otu_tables.xlsx"), sheet = 3)
lb$SampleID<-as.factor(lb$SampleID)
lb$SampleID<-factor(lb$SampleID, levels = rownames(b))
lb<-lb[order(lb$SampleID), ]

##Calculating distance matrix
library(vegan)
###Bacteria
bnorm<-decostand(b, method=c("normalize"))
bdist<-vegdist(bnorm, method = "jaccard")

bad1<-adonis(bdist~horizon, lb)
bad1

bad2<-adonis(bdist~horizon+PL/horizon, lb)
bad2

bmatrix<-as.matrix(bdist)

#Intercept 1 - Plesne - O - 0
#PL0:CTO = 0.75 : 0.25
PLO1<-bmatrix[intersect(which(lb$PL==0),which(lb$horizon=="O")),
               intersect(which(lb$PL==1),which(lb$horizon=="O"))]
PLO1[PLO1==0]<-NA

#PL0:CTO = 0.75 : 0.25
PLO75<-bmatrix[intersect(which(lb$PL==0),which(lb$horizon=="O")),
              intersect(which(lb$PL==0.75),which(lb$horizon=="O"))]
PLO75[PLO75==0]<-NA
#PL0:CTO = 0.5 : 0.5
PLO5<-bmatrix[intersect(which(lb$PL==0),which(lb$horizon=="O")),
               intersect(which(lb$PL==0.5),which(lb$horizon=="O"))]
PLO5[PLO5==0]<-NA
#PL0:CTO = 0.25 : 0.75
PLO25<-bmatrix[intersect(which(lb$PL==0),which(lb$horizon=="O")),
              intersect(which(lb$PL==0.25),which(lb$horizon=="O"))]
PLO25[PLO25==0]<-NA
#PL0:CTO = 0 : 1
PLO0<-bmatrix[intersect(which(lb$PL==0),which(lb$horizon=="O")),
               intersect(which(lb$PL==0),which(lb$horizon=="O"))]
PLO0[PLO0==0]<-NA
#Intercept 2 Plesne - A - 0
#PLA:CTA = 1 : 0
PLA1<-bmatrix[intersect(which(lb$PL==0),which(lb$horizon=="A")),
              intersect(which(lb$PL==1),which(lb$horizon=="A"))]
PLA1[PLA1==0]<-NA
#PLA:CTA = 0.75 : 0.25
PLA75<-bmatrix[intersect(which(lb$PL==0),which(lb$horizon=="A")),
              intersect(which(lb$PL==0.75),which(lb$horizon=="A"))]
PLO75[PLO75==0]<-NA
#PLA:CTA = 0.5 : 0.5
PLA5<-bmatrix[intersect(which(lb$PL==0),which(lb$horizon=="A")),
               intersect(which(lb$PL==0.5),which(lb$horizon=="A"))]
PLA5[PLA5==0]<-NA
#PLA:CTA = 0.75 : 0.25
PLA25<-bmatrix[intersect(which(lb$PL==0),which(lb$horizon=="A")),
              intersect(which(lb$PL==0.25),which(lb$horizon=="A"))]
PLA25[PLA25==0]<-NA
#PLA:CTA = 0 : 1
PLA0<-bmatrix[intersect(which(lb$PL==0),which(lb$horizon=="A")),
               intersect(which(lb$PL==0),which(lb$horizon=="A"))]
PLA0[PLA0==0]<-NA

bdiff<-data.frame(Distance = c(mean(log(PLO0/PLO1), na.rm = T), 
                               mean(log(PLO0/PLO75), na.rm = T),  
                               mean(log(PLO0/PLO5), na.rm = T),  
                               mean(log(PLO0/PLO25), na.rm = T),
                               mean(log(PLO0/PLO0), na.rm = T),  
                               mean(log(PLA0/PLA1), na.rm = T), 
                               mean(log(PLA0/PLA75), na.rm = T), 
                               mean(log(PLA0/PLA5), na.rm = T),
                               mean(log(PLA0/PLA25), na.rm = T), 
                               mean(log(PLA0/PLA0), na.rm = T)),
                  Distance.sd = c(sd(log(PLO0/PLO1), na.rm = T), 
                                  sd(log(PLO0/PLO75), na.rm = T),  
                                  sd(log(PLO0/PLO5), na.rm = T),  
                                  sd(log(PLO0/PLO25), na.rm = T),
                                  sd(log(PLO0), na.rm = T),  
                                  sd(log(PLA0/PLA1), na.rm = T), 
                                  sd(log(PLA0/PLA75), na.rm = T), 
                                  sd(log(PLA0/PLA5), na.rm = T),
                                  sd(log(PLA0/PLA25), na.rm = T), 
                                  sd(log(PLA0), na.rm = T)),
                  Plesne = rep(rev(seq(0,1, by = 0.25)), times = 2),
                  Horizon = c(rep("O", 5), rep("A", 5)))

ggplot(bdiff, aes(factor(Plesne), Distance))+geom_point(cex=6, pch=21)+
  facet_grid(.~Horizon, scales = "free")+theme_min+
  geom_errorbar(aes(ymin=Distance-Distance.sd, ymax=Distance+Distance.sd), width=0.1, lwd=0.5)+
  ylab(expression(paste("16S Distance")))+
  xlab("Plesne : Certovo mixing ratio")+theme(legend.position = c(0.2, 0.8),
                                            legend.key.size = unit(0.3, "in"),
                                            legend.title = element_blank())

###Fungi
fnorm<-decostand(f, method=c("normalize"))
fdist<-vegdist(fnorm, method = "jaccard")

fad1<-adonis(fdist~horizon, lb)
fad1

fad2<-adonis(fdist~horizon+PL/horizon, lb)
fad2

fmatrix<-as.matrix(fdist)

#Intercept 1 - Plesne - O - 0
#PL0:CTO = 0.75 : 0.25
fPLO1<-fmatrix[intersect(which(lb$PL==0),which(lb$horizon=="O")),
              intersect(which(lb$PL==1),which(lb$horizon=="O"))]
fPLO1[fPLO1==0]<-NA

#PL0:CTO = 0.75 : 0.25
fPLO75<-fmatrix[intersect(which(lb$PL==0),which(lb$horizon=="O")),
               intersect(which(lb$PL==0.75),which(lb$horizon=="O"))]
fPLO75[fPLO75==0]<-NA
#PL0:CTO = 0.5 : 0.5
fPLO5<-fmatrix[intersect(which(lb$PL==0),which(lb$horizon=="O")),
              intersect(which(lb$PL==0.5),which(lb$horizon=="O"))]
fPLO5[fPLO5==0]<-NA
#PL0:CTO = 0.25 : 0.75
fPLO25<-fmatrix[intersect(which(lb$PL==0),which(lb$horizon=="O")),
               intersect(which(lb$PL==0.25),which(lb$horizon=="O"))]
fPLO25[fPLO25==0]<-NA
#PL0:CTO = 0 : 1
fPLO0<-fmatrix[intersect(which(lb$PL==0),which(lb$horizon=="O")),
              intersect(which(lb$PL==0),which(lb$horizon=="O"))]
fPLO0[fPLO0==0]<-NA
#Intercept 2 Plesne - A - 0
#PLA:CTA = 1 : 0
fPLA1<-fmatrix[intersect(which(lb$PL==0),which(lb$horizon=="A")),
              intersect(which(lb$PL==1),which(lb$horizon=="A"))]
fPLA1[fPLA1==0]<-NA
#PLA:CTA = 0.75 : 0.25
fPLA75<-fmatrix[intersect(which(lb$PL==0),which(lb$horizon=="A")),
               intersect(which(lb$PL==0.75),which(lb$horizon=="A"))]
fPLO75[fPLO75==0]<-NA
#PLA:CTA = 0.5 : 0.5
fPLA5<-fmatrix[intersect(which(lb$PL==0),which(lb$horizon=="A")),
              intersect(which(lb$PL==0.5),which(lb$horizon=="A"))]
fPLA5[fPLA5==0]<-NA
#PLA:CTA = 0.75 : 0.25
fPLA25<-fmatrix[intersect(which(lb$PL==0),which(lb$horizon=="A")),
               intersect(which(lb$PL==0.25),which(lb$horizon=="A"))]
fPLA25[fPLA25==0]<-NA
#PLA:CTA = 0 : 1
fPLA0<-fmatrix[intersect(which(lb$PL==0),which(lb$horizon=="A")),
              intersect(which(lb$PL==0),which(lb$horizon=="A"))]
fPLA0[fPLA0==0]<-NA

fdiff<-data.frame(Distance = c(mean(log(fPLO0/fPLO1), na.rm = T), 
                               mean(log(fPLO0/fPLO75), na.rm = T),  
                               mean(log(fPLO0/fPLO5), na.rm = T),  
                               mean(log(fPLO0/fPLO25), na.rm = T),
                               mean(log(fPLO0/fPLO0), na.rm = T),  
                               mean(log(fPLA0/fPLA1), na.rm = T), 
                               mean(log(fPLA0/fPLA75), na.rm = T), 
                               mean(log(fPLA0/fPLA5), na.rm = T),
                               mean(log(fPLA0/fPLA25), na.rm = T), 
                               mean(log(fPLA0/fPLA0), na.rm = T)),
                  Distance.sd = c(sd(log(fPLO0/fPLO1), na.rm = T), 
                                  sd(log(fPLO0/fPLO75), na.rm = T),  
                                  sd(log(fPLO0/fPLO5), na.rm = T),  
                                  sd(log(fPLO0/fPLO25), na.rm = T),
                                  sd(log(fPLO0), na.rm = T),  
                                  sd(log(fPLA0/fPLA1), na.rm = T), 
                                  sd(log(fPLA0/fPLA75), na.rm = T), 
                                  sd(log(fPLA0/fPLA5), na.rm = T),
                                  sd(log(fPLA0/fPLA25), na.rm = T), 
                                  sd(log(fPLA0), na.rm = T)),
                  Plesne = rep(rev(seq(0,1, by = 0.25)), times = 2),
                  Horizon = c(rep("O", 5), rep("A", 5)))

ggplot(fdiff, aes(factor(Plesne), Distance))+geom_point(cex=6, pch=21)+
  facet_grid(.~Horizon, scales = "free")+theme_min+
  geom_errorbar(aes(ymin=Distance-Distance.sd, ymax=Distance+Distance.sd), width=0.1, lwd=0.5)+
  ylab(expression(paste("ITS Distance")))+
  xlab("Plesne : Certovo mixing ratio")+theme(legend.position = c(0.2, 0.8),
                                              legend.key.size = unit(0.3, "in"),
                                              legend.title = element_blank())
