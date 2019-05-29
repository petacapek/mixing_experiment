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

##Priming 
mich2$Priming<-numeric(length = nrow(mich2))
mich2[(mich2$cas==1 & mich2$znaceni=="Labelled"),"Priming" ]<-mich2[(mich2$cas==1 & mich2$znaceni=="Labelled"),"C.CO2.s" ]-mich2[(mich2$cas==1 & mich2$znaceni=="Unlabelled"),"C.CO2.s" ]

mich2[(mich2$cas==1 & mich2$znaceni=="Labelled"), ] %>% group_by(Plesne, znaceni, horizont) %>% summarize(y.sd=sd(Priming), y=mean(Priming)) %>%
  ggplot(aes(factor(Plesne), y))+geom_point(pch=21, cex=6)+
  geom_errorbar(aes(ymin=y-y.sd, ymax=y+y.sd), width=0.1, lwd=0.5)+
  facet_wrap(znaceni~horizont, scales="free")+theme_min+theme(legend.position = c(0.2, 0.88),
                                                              legend.key.size = unit(0.3, "in"),
                                                              legend.title = element_blank())+
  ylab(expression(paste(CO[2]," soil surplus ( ", mu, "mol ",g^{-1}, ")")))+
  xlab("Plesne : Certovo mixing ratio")+scale_fill_manual(values = c("white", "grey"))

##CUE 
mich2[(mich2$cas==1 & mich2$znaceni=="Labelled"), ] %>% group_by(Plesne, znaceni, horizont) %>% summarize(y.sd=sd(CUE), y=mean(CUE)) %>%
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
  group_by(Plesne, znaceni, horizont, Legend, variable) %>% summarize(y.sd=sd(delta), y=mean(delta)) %>%
  ggplot(aes(factor(Plesne), y))+geom_point(cex=6, aes(colour=Legend, shape=variable))+
  geom_errorbar(aes(ymin=y-y.sd, ymax=y+y.sd), width=0.1, lwd=0.5)+
  facet_wrap(znaceni~horizont, scales="free")+theme_min+theme(legend.position = "top",
                                                              legend.key.size = unit(0.3, "in"),
                                                              legend.title = element_blank())+
  ylab(expression(paste("Cmic (at%)" )))+
  xlab("Plesne : Certovo mixing ratio")+scale_fill_manual(values = c("white", "grey"))

Michizo[(Michizo$Plesne==0.5 & Michizo$znaceni=="Unlabelled" &
           Michizo$horizont=="O" & Michizo$Legend=="Before incubation" &
           Michizo$variable=="DOCatm"), ]

#Zmena signalu biomasy a DOC
ziv.zmeny$dDOCatm<-mich2[81:160, "DOCatm"]-mich2[1:80, "DOCatm"]

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

ggplot(ziv.zmeny[ziv.zmeny$znaceni=="Unlabelled", ],
       aes(dDOC2, dDOCatm))+geom_point(cex=6)+
  facet_wrap(.~horizont, scales="free")
