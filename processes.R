##########################################procesy#################################################
mich2<-mixing
mich2$Legend<-ifelse(mich2$cas==0, "Before incubation", "After incubation")
mich2$znaceni<-ifelse(mich2$znaceni=="JO", "Labelled", "Unlabelled")

##CO2
mich2[mich2$cas==1, ] %>% group_by(Plesne, znaceni, horizont) %>% summarize(y.sd=sd(CCO2c), y=mean(CCO2c)) %>%
ggplot(aes(factor(Plesne), y))+geom_point(pch=21, cex=6)+
  geom_errorbar(aes(ymin=y-y.sd, ymax=y+y.sd), width=0.1, lwd=0.5)+
  facet_grid(znaceni~horizont, scales="free")+theme_min+theme(legend.position = c(0.2, 0.88),
                                                              legend.key.size = unit(0.3, "in"),
                                                              legend.title = element_blank())+
  ylab(expression(paste("Cumulative ",CO[2]," ( ", mu, "mol ",g^{-1}, ")" )))+
  xlab("Plesne : Certovo mixing ratio")+scale_fill_manual(values = c("white", "grey"))

# ##obrat
# model_res$pars %>% group_by(Plesne, znaceni, horizont) %>% summarize(y.sd=sd(Vmax), y=mean(Vmax)) %>%
#   ggplot(aes(factor(Plesne), y))+geom_point(pch=21, cex=6)+
#   geom_errorbar(aes(ymin=y-y.sd, ymax=y+y.sd), width=0.1, lwd=0.5)+
#   facet_grid(znaceni~horizont, scales="free")+theme_min+theme(legend.position = c(0.2, 0.88),
#                                                               legend.key.size = unit(0.3, "in"),
#                                                               legend.title = element_blank())+
#  ylab(expression(paste(C[MB], " turnover rate ( ",day^{-1}, ")" )))+
#  xlab("Plesne : Certovo mixing ratio")+scale_fill_manual(values = c("white", "grey"))

##CO2 glukoza
mich2[(mich2$cas==1 & mich2$znaceni=="Labelled"), ] %>% group_by(Plesne, znaceni, horizont) %>% summarize(y.sd=sd(C.CO2.g), y=mean(C.CO2.g)) %>%
  ggplot(aes(factor(Plesne), y))+geom_point(pch=21, cex=6)+
  geom_errorbar(aes(ymin=y-y.sd, ymax=y+y.sd), width=0.1, lwd=0.5)+
  facet_grid(znaceni~horizont, scales="free")+theme_min+theme(legend.position = c(0.2, 0.88),
                                                              legend.key.size = unit(0.3, "in"),
                                                              legend.title = element_blank())+
  ylab(expression(paste("Cumulative glucose  ",CO[2]," ( ", mu, "mol ",g^{-1}, ")" )))+
  xlab("Plesne : Certovo mixing ratio")+scale_fill_manual(values = c("white", "grey"))

##Cmic glukoza
mich2[(mich2$cas==1 & mich2$znaceni=="Labelled"), ] %>% group_by(Plesne, znaceni, horizont) %>% summarize(y.sd=sd(Cmicg), y=mean(Cmicg)) %>%
  ggplot(aes(factor(Plesne), y))+geom_point(pch=21, cex=6)+
  geom_errorbar(aes(ymin=y-y.sd, ymax=y+y.sd), width=0.1, lwd=0.5)+
  facet_wrap(znaceni~horizont, scales="free")+theme_min+theme(legend.position = c(0.2, 0.88),
                                                              legend.key.size = unit(0.3, "in"),
                                                              legend.title = element_blank())+
  ylab(expression(paste(C[MB], " glucose ( ", mu, "mol ",g^{-1}, ")" )))+
  xlab("Plesne : Certovo mixing ratio")+scale_fill_manual(values = c("white", "grey"))

#Cmic glukoza vzhledem ke ztrate biomasy
XX<-mich2[(mich2$cas==1 & mich2$znaceni=="Labelled"), ]
XX$dCmb<-as.numeric(mb.zmeny[mb.zmeny$znaceni=="Labelled", "dCmic"])
XX$Cmic0<-mich2[(mich2$cas==0 & mich2$znaceni=="Labelled"), "Cmic"]
XX$dNH4<-as.numeric(ziv.zmeny[ziv.zmeny$znaceni=="Labelled", "dNH4"])
XX$dNO3<-as.numeric(ziv.zmeny[ziv.zmeny$znaceni=="Labelled", "dNO3"])
XX$dPO4<-as.numeric(ziv.zmeny[ziv.zmeny$znaceni=="Labelled", "dPO4"])
XX$dDON<-as.numeric(ziv.zmeny[ziv.zmeny$znaceni=="Labelled", "dDON"])
XX$dDOC<-as.numeric(ziv.zmeny[ziv.zmeny$znaceni=="Labelled", "dDOC"])

ggplot(XX, aes(Cmicg, CCO2c))+geom_point(cex=6)+theme_min+
  facet_wrap(.~Plesne, scales="free")

summary(lm(log(Cmicg)~Cmic0, XX, subset = Plesne == 1))
summary(lm(log(Cmicg)~Cmic0, XX, subset = Plesne == 0.75))
summary(lm(log(Cmicg)~Cmic0, XX, subset = Plesne == 0.5))
summary(lm(log(Cmicg)~Cmic0, XX, subset = Plesne == 0.25))
summary(lm(log(Cmicg)~Cmic0, XX, subset = Plesne == 0))

ggplot(XX, aes(CCO2c, -dNH4))+geom_point(cex=6)+theme_min+
  facet_wrap(.~Plesne, scales="free")

summary(lm(log(-dNH4)~Cmic0, XX, subset = Plesne == 1))
summary(lm(log(-dNH4)~Cmic0, XX, subset = Plesne == 0.75))
summary(lm(log(-dNH4)~Cmic0, XX, subset = Plesne == 0.5))
summary(lm(log(-dNH4)~Cmic0, XX, subset = Plesne == 0.25))
summary(lm(log(-dNH4)~Cmic0, XX, subset = Plesne == 0))

ggplot(XX, aes(CCO2c, -dNO3))+geom_point(cex=6)+theme_min+
  facet_wrap(.~Plesne, scales="free")
ggplot(XX, aes(CCO2c, -dPO4))+geom_point(cex=6)+theme_min+
  facet_wrap(.~Plesne, scales="free")
ggplot(XX, aes(CCO2c, -dDON))+geom_point(cex=6)+theme_min+
  facet_wrap(.~Plesne, scales="free")
ggplot(XX, aes(-dNO3, -dDOC))+geom_point(cex=6)+theme_min#+
  facet_wrap(.~Plesne, scales="free")
  
#Cmic glukoza na jednotku puvodni Cmic
mich0<-mich2[mich2$cas==0, ]
mich0$CmicgCmic<-mich2[mich2$cas==1, "Cmicg"]/mich0$Cmic

mich0[mich0$znaceni=="Labelled", ] %>% group_by(Plesne, horizont) %>% 
  summarize(y.sd=sd(CmicgCmic), y=mean(CmicgCmic)) %>%
  ggplot(aes(factor(Plesne), y))+geom_point(pch=21, cex=6)+
  geom_errorbar(aes(ymin=y-y.sd, ymax=y+y.sd), width=0.1, lwd=0.5)+
  facet_wrap(~horizont, scales="free")+theme_min+theme(legend.position = c(0.2, 0.88),
                                                              legend.key.size = unit(0.3, "in"),
                                                              legend.title = element_blank())+
  ylab(expression(paste("Relative amount of glucose derived microbial biomass")))+
  xlab("Plesne : Certovo mixing ratio")

##Priming 
mich2$Priming<-numeric(length = nrow(mich2))
mich2[(mich2$cas==1 & mich2$znaceni=="Labelled"),"Priming" ]<-(mich2[(mich2$cas==1 & mich2$znaceni=="Labelled"),"CCO2.s" ]-mich2[(mich2$cas==1 & mich2$znaceni=="Unlabelled"),"CCO2.s" ])/mich2[(mich2$cas==1 & mich2$znaceni=="Unlabelled"),"CCO2.s" ]

mich2[(mich2$cas==1 & mich2$znaceni=="Labelled"), ] %>% group_by(Plesne, znaceni, horizont) %>% summarize(y.sd=sd(Priming), y=mean(Priming)) %>%
  ggplot(aes(factor(Plesne), y))+geom_point(pch=21, cex=6)+
  geom_errorbar(aes(ymin=y-y.sd, ymax=y+y.sd), width=0.1, lwd=0.5)+
  facet_wrap(znaceni~horizont, scales="free")+theme_min+theme(legend.position = c(0.2, 0.88),
                                                              legend.key.size = unit(0.3, "in"),
                                                              legend.title = element_blank())+
  ylab(expression(paste("Relative priming effect")))+
  xlab("Plesne : Certovo mixing ratio")+scale_fill_manual(values = c("white", "grey"))

##CUE 
mich2[(mich2$cas==1 & mich2$znaceni=="Labelled"), ] %>% group_by(Plesne, znaceni, horizont) %>% summarize(y.sd=sd(Cmicg/(CCO2.g+Cmicg)), 
                                                                                                          y=mean(Cmicg/(CCO2.g+Cmicg))) %>%
  ggplot(aes(factor(Plesne), y))+geom_point(pch=21, cex=6)+
  geom_errorbar(aes(ymin=y-y.sd, ymax=y+y.sd), width=0.1, lwd=0.5)+
  facet_wrap(znaceni~horizont, scales="free")+theme_min+theme(legend.position = c(0.2, 0.88),
                                                              legend.key.size = unit(0.3, "in"),
                                                              legend.title = element_blank())+
  ylab(expression(paste(CUE)))+
  xlab("Plesne : Certovo mixing ratio")+scale_fill_manual(values = c("white", "grey"))


mich2[(mich2$cas==1 & mich2$znaceni=="Labelled"), ] %>% group_by(Plesne, znaceni, horizont) %>% summarize(y.sd=sd(CUE), y=mean(CUE),
                                                                                                          x.sd=sd(pH), x=mean(pH)) %>%
  ggplot(aes(x,y))+geom_point(cex=6, pch=21)+theme_min+geom_errorbar(width=0.1, aes(ymin=y-y.sd,
                                                                                    ymax=y+y.sd))+
  geom_errorbarh(width=0.1, aes(xmin=x-x.sd, xmax=x+x.sd))+ylab("CUE")+xlab("pH")

##C residualni 
mich2[(mich2$cas==1 & mich2$znaceni=="Labelled"), ] %>% group_by(Plesne, znaceni, horizont) %>% summarize(y.sd=sd(100-Gl.found), y=mean(100-Gl.found)) %>%
  ggplot(aes(factor(Plesne), y))+geom_point(pch=21, cex=6)+
  geom_errorbar(aes(ymin=y-y.sd, ymax=y+y.sd), width=0.1, lwd=0.5)+
  facet_wrap(~horizont)+theme_min+theme(legend.position = c(0.2, 0.88),
                                                              legend.key.size = unit(0.3, "in"),
                                                              legend.title = element_blank())+
  ylab(expression(paste("Residual C (%)")))+
  xlab("Plesne : Certovo mixing ratio")+scale_fill_manual(values = c("white", "grey"))

ziv.zmeny$Cres<-(100-mich2[(mich2$cas==1), "Gl.found"])*mich2[(mich2$cas==1), "Gl"]/100

mich2[(mich2$cas==1 & mich2$znaceni=="Labelled"), ] %>% group_by(Plesne, horizont) %>% 
  summarize(y.sd=sd((100-Gl.found)), y=mean((100-Gl.found)), x.sd=sd(bg), x=mean(bg)) %>%
  ggplot(aes(x, y))+geom_point(pch=21, cex=6, aes(fill=horizont))+#facet_wrap(~horizont, scales="free")+
  geom_errorbar(aes(ymin=y-y.sd, ymax=y+y.sd), width=0.1, lwd=0.5)+
  geom_errorbarh(aes(xmin=x-x.sd, xmax=x+x.sd), width=0.1, lwd=0.5)+
  facet_wrap(~horizont)+
  theme_min+theme(legend.position = c(0.2, 0.88),
                                        legend.key.size = unit(0.3, "in"),
                                        legend.title = element_blank())+
  ylab(expression(paste("Residual C (%)")))+
  xlab("Plesne : Certovo mixing ratio")+scale_fill_manual(values = c("white", "grey"))

#Nalezena glukoza
mich2[(mich2$cas==1 & mich2$znaceni=="Labelled"), ] %>% group_by(Plesne, horizont) %>% 
  summarize(y.sd=sd((100-Glfound)), y=mean((100-Glfound))) %>%
  ggplot(aes(factor(Plesne), y))+geom_point(pch=21, cex=6, aes(fill=horizont), show.legend = F)+
  facet_wrap(~horizont, scales="free")+
  geom_errorbar(aes(ymin=y-y.sd, ymax=y+y.sd), width=0.1, lwd=0.5)+
  theme_min+theme(legend.position = c(0.2, 0.88),
                  legend.key.size = unit(0.3, "in"),
                  legend.title = element_blank())+
  ylab(expression(paste("Residual C (%)")))+
  xlab("Plesne : Certovo mixing ratio")+scale_fill_manual(values = c("white", "grey"))

##13C
##CO2
mich2 %>% group_by(Plesne, znaceni, horizont, Legend) %>% summarize(y.sd=sd(CO2atm), y=mean(CO2atm)) %>%
  ggplot(aes(factor(Plesne), y))+geom_point(pch=21, cex=6)+
  geom_errorbar(aes(ymin=y-y.sd, ymax=y+y.sd), width=0.1, lwd=0.5)+
  facet_grid(znaceni~horizont, scales="free")+theme_min+theme(legend.position = c(0.2, 0.88),
                                                              legend.key.size = unit(0.3, "in"),
                                                              legend.title = element_blank())+
  ylab(expression(paste(CO[2]," (at%)" )))+
  xlab("Plesne : Certovo mixing ratio")+scale_fill_manual(values = c("white", "grey"))

##DOC
mich2 %>% group_by(Plesne, znaceni, horizont, Legend) %>% summarize(y.sd=sd(DOCatm), y=mean(DOCatm)) %>%
  ggplot(aes(factor(Plesne), y))+geom_point(pch=21, cex=6, aes(fill=Legend))+
  geom_errorbar(aes(ymin=y-y.sd, ymax=y+y.sd), width=0.1, lwd=0.5)+
  facet_wrap(znaceni~horizont, scales="free")+theme_min+theme(legend.position = c(0.2, 0.88),
                                                              legend.key.size = unit(0.3, "in"),
                                                              legend.title = element_blank())+
  ylab(expression(paste("DOC (at%)" )))+
  xlab("Plesne : Certovo mixing ratio")+scale_fill_manual(values = c("white", "grey"))

##Cmic
mich2 %>% group_by(Plesne, znaceni, horizont, Legend) %>% summarize(y.sd=sd(Cmicatm), y=mean(Cmicatm)) %>%
  ggplot(aes(factor(Plesne), y))+geom_point(pch=21, cex=6, aes(fill=Legend))+
  geom_errorbar(aes(ymin=y-y.sd, ymax=y+y.sd), width=0.1, lwd=0.5)+
  facet_wrap(znaceni~horizont, scales="free")+theme_min+theme(legend.position = c(0.2, 0.88),
                                                              legend.key.size = unit(0.3, "in"),
                                                              legend.title = element_blank())+
  ylab(expression(paste("Cmic (at%)" )))+
  xlab("Plesne : Certovo mixing ratio")+scale_fill_manual(values = c("white", "grey"))

michizo<-mich2[, c("Plesne", "Certovo", "horizont", "Legend", "znaceni", "CO2atm",
                   "DOCatm", "Cmicatm")]

Michizo<-melt(michizo, id.vars = c("Plesne", "Certovo", "horizont", "Legend", "znaceni"))
Michizo$delta<-with(Michizo, (value/(1-value)/0.011237-1)*1000)

Michizo[-173, ] %>% filter(znaceni!="Labelled") %>%
  group_by(Plesne, znaceni, horizont, variable) %>% summarize(y.sd=sd(delta, na.rm = T), 
                                                              y=mean(delta, na.rm = T)) %>%
  ggplot(aes(factor(Plesne), y))+geom_point(cex=6, pch=21, aes(fill=variable))+
  geom_errorbar(aes(ymin=y-y.sd, ymax=y+y.sd), width=0.1, lwd=0.5)+
  facet_wrap(znaceni~horizont, scales="free")+theme_min+theme(legend.position = "top",
                                                              legend.key.size = unit(0.3, "in"),
                                                              legend.title = element_blank())+
  ylab(expression(paste(delta^{13},C)))+
  xlab("Plesne : Certovo mixing ratio")

Michizo[(Michizo$Plesne==0.5 & Michizo$znaceni=="Unlabelled" &
           Michizo$horizont=="O" & Michizo$Legend=="Before incubation" &
           Michizo$variable=="DOCatm"), ]

#Zmena signalu biomasy a DOC
ziv.zmeny$dDOCatm<-mich2[81:160, "DOCatm"]-mich2[1:80, "DOCatm"]
ziv.zmeny$dCmicCO2atm<-mich2[81:160, "Cmicatm"]+mich2[81:160, "CO2atm"]
ziv.zmeny$dCmic<-mb.zmeny$dCmic
ziv.zmeny$CO2c<-mich2[81:160, "CCO2c"]
ziv.zmeny$CO2atm<-mich2[81:160, "CO2atm"]
ziv.zmeny$Cmicatm<-mich2[1:80, "Cmicatm"]
ziv.zmeny$DOCatm<-mich2[1:80, "DOCatm"]


ziv.zmeny %>% filter(znaceni=="Unlabelled") %>%
  group_by(Plesne, znaceni, horizont) %>% summarize(y.sd=sd(dDOCatm), y=mean(dDOCatm)) %>%
  ggplot(aes(factor(Plesne), y))+geom_point(cex=6, pch=21, aes(fill=znaceni))+
  facet_grid(.~horizont, scales="free")+theme_min+theme(legend.position = c(0.2, 0.2))+
  geom_errorbar(aes(ymin=y-y.sd, ymax=y+y.sd), width=0.1, lwd=0.5)+
  theme(legend.position = c(0.15, 0.2),legend.key.size = unit(0.3, "in"),
        legend.title = element_blank())+
  ylab(expression(paste(Delta~DOC, " (atm%)" )))+
  xlab("Plesne : Certovo mixing ratio")+scale_fill_manual(values = c("white", "grey"))+
  geom_hline(yintercept = 0, lwd=1.5)+theme(legend.title = element_blank())

ggplot(ziv.zmeny,  aes(dDOC2*DOCatm, CO2atm))+geom_point(cex=6, aes(colour=factor(Plesne)))+
  facet_wrap(znaceni~horizont, scales="free")+
  stat_smooth(method=lm)

####Ztraty
#####C
mich1<-mich2[mich2$cas==1, ]
ztraty<-mich2[mich2$cas==0, c("Plesne", "Certovo", "horizont", "znaceni", "cas", "Legend")]
ztraty$CCO2<-mich1[, "CCO2"]
ztraty$Cmic<-mich0[, "Cmic"]-mich1[, "Cmic"]
ztraty$DOC<-mich0[, "DOC2"]-mich1[, "DOC2"]

Ztraty<-melt(ztraty, id.vars = c("Plesne", "Certovo", "horizont", "znaceni", "cas", "Legend"))

Ztraty %>%  group_by(Plesne, znaceni, horizont, variable) %>% summarize(y.sd=sd(value), y=mean(value)) %>%
  ggplot(aes(factor(Plesne), y))+geom_point(cex=6, pch=21, aes(fill=variable))+
  facet_wrap(znaceni~horizont, scales="free")+theme_min+
  geom_errorbar(aes(ymin=y-y.sd, ymax=y+y.sd), width=0.1, lwd=0.5)+
  theme(legend.position = c(0.1, 0.4),legend.key.size = unit(0.3, "in"),
        legend.title = element_blank())+
  ylab(expression(paste("C losses (", mu, "mol ", g^{-1}, ")")))+
  xlab("Plesne : Certovo mixing ratio")+
  theme(legend.title = element_blank())

######N
nz<-mich2[mich2$cas==0, c("Plesne", "Certovo", "horizont", "znaceni", "cas", "Legend")]
nz$Nmic<-mich0[, "Nmic"]-mich1[, "Nmic"]
nz$DON<-mich0[, "DON"]-mich1[, "DON"]
nz$NH4<-mich0[, "NH42"]-mich1[, "NH42"]
#nz$NO3<-mich0[, "NO3"]-mich1[, "NO3"]

#Mineralizace dusiku
nz$leuVmax<-mich2e[mich2e$cas==0, "leuVmax"]
nz$leuKm<-mich2e[mich2e$cas==0, "leuKm"]
nz$chitVmax<-mich2e[mich2e$cas==0, "chitVmax"]
nz$chitKm<-mich2e[mich2e$cas==0, "chitKm"]
nz$DON0<-mich0[, "DON"]

nz$DONmin<-with(nz, chitVmax*DON0/(chitKm+DON0)*48+leuVmax*DON0/(leuKm+DON0)*48)

ggplot(nz, aes(DONmin, DON))+geom_point(cex=6, aes(colour=horizont))+
  facet_wrap(znaceni~horizont, scales="free")

Nz<-melt(nz, id.vars = c("Plesne", "Certovo", "horizont", "znaceni", "cas", "Legend"))

Nz %>%  group_by(Plesne, znaceni, horizont, variable) %>% summarize(y.sd=sd(value, na.rm = T), 
                                                                    y=mean(value, na.rm = T)) %>%
  ggplot(aes(factor(Plesne), y))+geom_point(cex=6, pch=21, aes(fill=variable))+
  facet_wrap(znaceni~horizont, scales="free")+theme_min+
  geom_errorbar(aes(ymin=y-y.sd, ymax=y+y.sd), width=0.1, lwd=0.5)+
  theme(legend.position = c(0.1, 0.4),legend.key.size = unit(0.3, "in"),
        legend.title = element_blank())+
  ylab(expression(paste("N losses (", mu, "mol ", g^{-1}, ")")))+
  xlab("Plesne : Certovo mixing ratio")+
  theme(legend.title = element_blank())

######P
pz<-mich2[mich2$cas==0, c("Plesne", "Certovo", "horizont", "znaceni", "cas", "Legend")]
pz$Pmic<-mich0[, "Pmic"]-mich1[, "Pmic"]
pz$PO4<-mich0[, "PO4"]-mich1[, "PO4"]
pz$PO42<-mich0[, "PO42"]-mich1[, "PO42"]
#nz$NO3<-mich0[, "NO3"]-mich1[, "NO3"]

Pz<-melt(pz, id.vars = c("Plesne", "Certovo", "horizont", "znaceni", "cas", "Legend"))

Pz %>%  group_by(Plesne, znaceni, horizont, variable) %>% summarize(y.sd=sd(value, na.rm = T), 
                                                                    y=mean(value, na.rm = T)) %>%
  ggplot(aes(factor(Plesne), y))+geom_point(cex=6, pch=21, aes(fill=variable))+
  facet_wrap(znaceni~horizont, scales="free")+theme_min+
  geom_errorbar(aes(ymin=y-y.sd, ymax=y+y.sd), width=0.1, lwd=0.5)+
  theme(legend.position = c(0.1, 0.4),legend.key.size = unit(0.3, "in"),
        legend.title = element_blank())+
  ylab(expression(paste("P losses (", mu, "mol ", g^{-1}, ")")))+
  xlab("Plesne : Certovo mixing ratio")+
  theme(legend.title = element_blank())


#Bilance fosforu
##Potreba fosforu
mich1$Sorption<-ifelse(mich1$Sorption>1, 1, mich1$Sorption)
pz$Pdem<-mich1[,"CCO2c"]*mich1[, "CUE"]/(1-mich1[, "CUE"])/(mich1[, "Cmic"]/mich1[, "Pmic"])
pz$Sorption<-mich1$Sorption
summary(pz$Pdem)

ggplot(pz[pz$znaceni=="Unlabelled", ], aes(1-Sorption, PO42))+geom_point(cex=6, aes(colour=factor(Plesne)))+
  stat_smooth(method = lm, se = F, formula = y~poly(x,2))+facet_wrap(~horizont, scales="free")

summary(lm(PO42~Pdem, pz[pz$znaceni=="Unlabelled", ]))

#Mineralizace fosforu
pz$phosVmax<-mich2e[mich2e$cas==0, "phosVmax"]
pz$phosKm<-mich2e[mich2e$cas==0, "phosKm"]

pz$Penzrel<-with(pz, Pmic*phosVmax/(phosKm+Pmic)*48)

ggplot(pz[pz$znaceni=="Unlabelled", ], aes(Pdem, PO42+Penzrel))+geom_point(cex=6, aes(colour=horizont))+
  facet_wrap(~horizont, scales = "free")

#Rychlost obratu
Obrat<-mich0[, c("Plesne", "Certovo", "horizont", "znaceni", "Legend")]
Obrat$dCmic<-mich0$Cmic-mich1$Cmic
Obrat$rust<-mich1$CCO2c*mich1$CUE/(1-mich1$CUE)
Obrat$obrat<-(Obrat$dCmic+Obrat$rust)/mich0$Cmic/48

Obrat %>% filter(znaceni=="Unlabelled") %>% group_by(Plesne, horizont) %>% 
  summarize(y.sd=sd(obrat, na.rm = T),
            y=mean(obrat, na.rm = T)) %>%
  ggplot(aes(factor(Plesne), y))+geom_point(cex=6, pch=21)+
  facet_wrap(~horizont, scales="free")+theme_min+
  geom_errorbar(aes(ymin=y-y.sd, ymax=y+y.sd), width=0.1, lwd=0.5)+
  theme(legend.position = c(0.1, 0.4),legend.key.size = unit(0.3, "in"),
        legend.title = element_blank())+
  ylab(expression(paste("Microbial biomass turnover rate")))+
  xlab("Plesne : Certovo mixing ratio")+
  theme(legend.title = element_blank())

#Correlation with microbial community composition
bdiff$Legend<-"Bacteria"
fdiff$Legend<-"Fungi"

Community<-rbind(bdiff, fdiff)

Community$obrat<-c(as.data.frame(Obrat %>% filter(znaceni=="Unlabelled") %>%
                                 group_by(Plesne, horizont) %>%
                                 summarize(mean(obrat)))[,3],
                   as.data.frame(Obrat %>% filter(znaceni=="Unlabelled") %>%
                                   group_by(Plesne, horizont) %>%
                                   summarize(mean(obrat)))[,3])
Community$obratsd<-c(as.data.frame(Obrat %>% filter(znaceni=="Unlabelled") %>%
                                   group_by(Plesne, horizont) %>%
                                   summarize(sd(obrat)))[,3],
                   as.data.frame(Obrat %>% filter(znaceni=="Unlabelled") %>%
                                   group_by(Plesne, horizont) %>%
                                   summarize(sd(obrat)))[,3])

ggplot(Community, aes(Distance, obrat))+geom_point(cex=6, aes(colour=Legend))+
  #geom_errorbar(aes(ymin = obrat-obratsd, ymax=obrat+obratsd))+
  #geom_errorbarh(aes(xmin=Distance-Distance.sd, xmax=Distance+Distance.sd))+
  facet_wrap(Horizon~Legend, scales="free")

Community$CUE<-c(as.data.frame(mich1 %>% filter(znaceni=="Labelled") %>%
                                 group_by(Plesne, horizont) %>%
                                 summarize(mean(CUE)))[,3],
                 as.data.frame(mich1 %>% filter(znaceni=="Labelled") %>%
                                 group_by(Plesne, horizont) %>%
                                 summarize(mean(CUE)))[,3])

Community$CUEsd<-c(as.data.frame(mich1 %>% filter(znaceni=="Labelled") %>%
                                 group_by(Plesne, horizont) %>%
                                 summarize(sd(CUE)))[,3],
                 as.data.frame(mich1 %>% filter(znaceni=="Labelled") %>%
                                 group_by(Plesne, horizont) %>%
                                 summarize(sd(CUE)))[,3])

ggplot(Community, aes(Distance, CUE))+geom_point(cex=6, aes(colour=Legend))+
  geom_errorbar(aes(ymin = CUE-CUEsd, ymax=CUE+CUEsd))+
  geom_errorbarh(aes(xmin=Distance-Distance.sd, xmax=Distance+Distance.sd))+
  facet_wrap(Legend~Horizon, scales = "free")


#Sorpce podel Kani
CTsorp<-read.csv("clipboard", header = F)
PLsorp<-read.csv("clipboard", header = F)

colnames(CTsorp)<-c("Padd", "Psorp")
colnames(PLsorp)<-c("Padd", "Psorp")

CTsorp$Padd<-CTsorp$Padd*0.03*1000
PLsorp$Padd<-PLsorp$Padd*0.03*1000

CTiso<-nls(Psorp~Xmax*Padd/(Kx+Padd), CTsorp, start = list(Xmax=20, Kx=10))
summary(CTiso)
PLiso<-nls(Psorp~Xmax*Padd/(Kx+Padd), PLsorp, start = list(Xmax=20, Kx=10))
summary(PLiso)


pzA<-pz[pz$horizont=="A", ]
pzA$PO4add<-mich0[mich0$horizont=="A", "PO42"]

pzA$Psorp<-with(pzA, coef(CTiso)[1]*PO4add/(coef(CTiso)[2]+PO4add)*Certovo+
                  coef(PLiso)[1]*PO4add/(coef(PLiso)[2]+PO4add)*Plesne)

ggplot(pzA, aes(Psorp, PO42))+geom_point(cex=6, aes(colour=znaceni))
