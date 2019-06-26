##########################################procesy#################################################
mich2<-mixing
mich2$Legend<-ifelse(mich2$cas==0, "Before incubation", "After incubation")
mich2$znaceni<-ifelse(mich2$znaceni=="JO", "Labelled", "Unlabelled")
mich2$DOP<-read.xlsx(xlsxFile = c("./Data_overview/prehled.xlsx"), 1)[, 39]


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


zmeny_vse %>% filter(znaceni=="Unlabelled") %>% group_by(Plesne, Horizon) %>% 
  summarize(y.sd=sd(CUE2), y=mean(CUE2)) %>%
  ggplot(aes(factor(Plesne),y))+geom_point(cex=6, pch=21, aes(fill=Horizon))+
  theme_min+geom_errorbar(width=0.1, aes(ymin=y-y.sd,ymax=y+y.sd))+
  scale_fill_manual(values = c("white", "grey"))+
  theme(legend.title = element_blank(),
        legend.position = c(0.8, 0.15))+
  ylim(0.4, 0.9)+
  ylab("CUE")+xlab("Plesne : Certovo mixing ratio")+
  annotate("text", label="Horizon***", 0.5, 0.5, hjust="left", size=6, fontface="italic")+
  annotate("text", label="Plesne : Certovo**", 0.5, 0.45, hjust="left", size=6, fontface="italic")+
  ggtitle("A")

anova(glm(CUE2~horizont+Plesne/horizont, zmeny_vse, subset = znaceni=="Unlabelled",
          family=Gamma), test="F")

zmeny_vse %>% filter(znaceni=="Unlabelled") %>% group_by(Plesne, Horizon) %>% 
  summarize(y.sd=sd(obrat2, na.rm = T), y=mean(obrat2, na.rm = T)) %>%
  ggplot(aes(factor(Plesne),y))+geom_point(cex=6, pch=21, aes(fill=Horizon), show.legend = F)+
  theme_min+geom_errorbar(width=0.1, aes(ymin=y-y.sd,ymax=y+y.sd))+
  scale_fill_manual(values = c("white", "grey"))+
  theme(legend.title = element_blank(),
        legend.position = c(0.9, 0.15))+
  ylim(0, 0.09)+
  ylab(expression(paste("turnover rate (", h^{-1}, ")")))+
  xlab("Plesne : Certovo mixing ratio")+
  annotate("text", label="Horizon n.s.", 0.5, 0.005, hjust="left", size=6, fontface="italic")+
  annotate("text", label="Plesne : Certovo n.s", 0.5, 0, hjust="left", size=6, fontface="italic")+
  ggtitle("B")

anova(glm(obrat2~horizont+Plesne/horizont, zmeny_vse, subset = znaceni=="Unlabelled",
          family=Gamma), test="F")

zmeny_vse %>% filter(znaceni=="Unlabelled") %>% group_by(Plesne, Horizon) %>% 
  summarize(y.sd=sd(CO2un, na.rm = T), y=mean(CO2un, na.rm = T)) %>%
  ggplot(aes(factor(Plesne),y))+geom_point(cex=6, pch=21, aes(fill=Horizon), show.legend = T)+
  theme_min+geom_errorbar(width=0.1, aes(ymin=y-y.sd,ymax=y+y.sd))+
  scale_fill_manual(values = c("white", "grey"))+
  theme(legend.title = element_blank(),
        legend.position = c(0.8, 0.15))+
  ylim(0, 25)+
  ylab(expression(paste("Cumulative ",CO[2]," loss  (", mu, "mol ", g^{-1}, ")")))+
  xlab("Plesne : Certovo mixing ratio")+
  annotate("text", label="Horizon***", 0.5, 24, hjust="left", size=6, fontface="italic")+
  annotate("text", label="Plesne : Certovo**", 0.5, 22, hjust="left", size=6, fontface="italic")

anova(glm(CO2un~horizont+Plesne/horizont, zmeny_vse, subset = znaceni=="Unlabelled",
          family=Gamma), test="F")

mich2[(mich2$cas==1 & mich2$znaceni=="Labelled"), ] %>% group_by(Plesne, znaceni, horizont) %>% summarize(y.sd=sd(Cmicg/(CCO2.g+Cmicg)), 
                                                                                                          y=mean(Cmicg/(CCO2.g+Cmicg))) %>%
  ggplot(aes(factor(Plesne), y))+geom_point(pch=21, cex=6)+
  geom_errorbar(aes(ymin=y-y.sd, ymax=y+y.sd), width=0.1, lwd=0.5)+
  facet_wrap(znaceni~horizont, scales="free")+theme_min+theme(legend.position = c(0.2, 0.88),
                                                              legend.key.size = unit(0.3, "in"),
                                                              legend.title = element_blank())+
  ylab(expression(paste(CUE)))+
  xlab("Plesne : Certovo mixing ratio")+scale_fill_manual(values = c("white", "grey"))


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

Michizo$Horizon<-Michizo$horizont
Michizo$Horizon<-ifelse(Michizo$horizont=="O", "Litter", "Organic topsoil")

Michizo %>% filter(znaceni!="Labelled" & Legend=="After incubation") %>%
  group_by(Plesne, Horizon, variable) %>% 
  summarize(y.sd=sd(delta, na.rm = T), y=mean(delta, na.rm = T)) %>%
  ggplot(aes(factor(Plesne), y))+geom_point(cex=6, pch=21, aes(fill=variable))+
  geom_errorbar(aes(ymin=y-y.sd, ymax=y+y.sd), width=0.1, lwd=0.5)+
  facet_grid(.~Horizon)+theme_min+
  theme(legend.position = c(0.1, 0.85),
        legend.title = element_blank(),
        legend.text.align = 0)+
  ylab(expression(paste(delta^{13},C)))+
  xlab("Plesne : Certovo mixing ratio")+
  scale_fill_manual(values = c("black", "grey", "white"),
                    name = '',
                    labels = expression(C-CO[2], K[2]~SO[4]-EC, MBC))

anova(glm(value~variable+Plesne/variable, 
          data=subset(Michizo, znaceni=="Unlabelled" & Legend=="After incubation" & horizont=="A" & variable!="CO2atm"),
          family = Gamma), test="F")

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
ztraty$CCO2<-zmeny_vse$CO2un
ztraty$Cmic<-(mich0[, "Cmic"]-mich1[, "Cmic"])/0.38
ztraty$DOC<-mich0[, "DOC"]-mich1[, "DOC"]

Ztraty<-melt(ztraty, id.vars = c("Plesne", "Certovo", "horizont", "znaceni", "cas", "Legend"))
Ztraty$Horizon<-ifelse(Ztraty$horizont=="O", "Litter", "Organic topsoil")

Ztraty %>% filter(znaceni=="Unlabelled") %>% group_by(Plesne, Horizon, variable) %>% 
  summarize(y.sd=sd(value), y=mean(value)) %>%
  ggplot(aes(factor(Plesne), y))+geom_point(cex=6, pch=21, aes(fill=variable))+
  facet_wrap(~Horizon, scales="free")+theme_min+
  geom_errorbar(aes(ymin=y-y.sd, ymax=y+y.sd), width=0.1, lwd=0.5)+
  theme(legend.position = c(0.65, 0.8),legend.key.size = unit(0.3, "in"),
        legend.title = element_blank(),
        legend.text.align = 0)+
  ylab(expression(paste("C losses (", mu, "mol ", g^{-1}, ")")))+
  xlab("Plesne : Certovo mixing ratio")+
  theme(legend.title = element_blank())+
  geom_hline(yintercept = 0)+
  scale_fill_manual(values = c("black", "grey", "white"),
                    name = '',
                    labels = expression(C-CO[2], MBC, WEC))

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
nz$Ndem<-mich1[,"CCO2c"]*mich1[, "CUE"]/(1-mich1[, "CUE"])/(mich1[, "Cmic"]/mich1[, "Nmic"])


nz$DONmin<-with(nz, chitVmax*DON0/(chitKm+DON0)*48+leuVmax*DON0/(leuKm+DON0)*48)

ggplot(nz, aes(Ndem, DON+NH4))+geom_point(cex=6, aes(colour=horizont))+
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
pz$Pdem<-mich1[,"CCO2c"]*mich1[, "CUE"]/(1-mich1[, "CUE"])/(mich0[, "Cmic"]/mich0[, "Pmic"])
pz$Sorption<-mich1$Sorption
pz$pH<-mich0$pH
pz$Psorbed<-mich0$PO42-(mich0$PO42*exp(-mich0$Sorption/0.75*48))
pz$Psorbed2<--mich0$PO42-mich0$PO4*exp(-mich0$Sorption/0.75*48)


ggplot(pz[pz$znaceni=="Unlabelled", ], aes(Pdem, PO42))+geom_point(cex=6, aes(colour=horizont))+
  facet_wrap(~horizont, scales="free")+stat_smooth(method=lm)


summary(lm(PO42~Psorbed, pz[pz$znaceni=="Unlabelled", ], subset = horizont=="O"))
summary(lm(PO42~Psorbed, pz[pz$znaceni=="Unlabelled", ], subset = horizont=="A"))

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
Obrat$rust<-with(zmeny_vse, CO2un*CUE2/(1-CUE2))
Obrat$obrat<-(Obrat$dCmic+Obrat$rust)/mich0$Cmic/48
Obrat$obrat2<-(-log(Obrat$dCmic/0.38+Obrat$rust)+log(mich0$Cmic/0.38))/48
Obrat$obrat3<-(-log(Obrat$dCmic/0.38+Obrat$rust)+log(mich0$Cmic/0.38))/48
Obrat$Priming<-c(mich2[(mich2$cas==1 & mich2$znaceni=="Labelled"),"Priming" ], 
                 mich2[(mich2$cas==1 & mich2$znaceni=="Labelled"),"Priming" ])
Obrat$Sorption<-mich0$Sorption

ggplot(Obrat, aes(Sorption, rust))+geom_point(cex=6, aes(color=horizont))+
  facet_wrap(horizont~znaceni, scales="free")

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

Community$Priming<-c(as.data.frame(mich2 %>% filter(znaceni=="Labelled") %>%
                                 group_by(Plesne, horizont) %>%
                                 summarize(mean(Priming, na.rm = T)))[,3],
                     as.data.frame(mich2 %>% filter(znaceni=="Labelled") %>%
                                     group_by(Plesne, horizont) %>%
                                     summarize(mean(Priming, na.rm = T)))[,3])

Community$Primingsd<-c(as.data.frame(mich2 %>% filter(znaceni=="Labelled") %>%
                                   group_by(Plesne, horizont) %>%
                                   summarize(sd(Priming, na.rm = T)))[,3],
                   as.data.frame(mich2 %>% filter(znaceni=="Labelled") %>%
                                   group_by(Plesne, horizont) %>%
                                   summarize(sd(Priming, na.rm = T)))[,3])

ggplot(Community, aes(Distance, CUE))+geom_point(cex=6, aes(colour=Legend))+
  geom_errorbar(aes(ymin = CUE-CUEsd, ymax=CUE+CUEsd))+
  geom_errorbarh(aes(xmin=Distance-Distance.sd, xmax=Distance+Distance.sd))+
  facet_wrap(Legend~Horizon, scales = "free")

ggplot(Community, aes(Distance, Priming))+geom_point(cex=6, aes(colour=Horizon))+
  geom_errorbar(aes(ymin = Priming-Primingsd, ymax=Priming+Primingsd))+
  geom_errorbarh(aes(xmin=Distance-Distance.sd, xmax=Distance+Distance.sd))+
  facet_wrap(~Legend, scales = "free")

ggplot(mich2, aes(1-Sorption, Glfound))+geom_point(cex=6)

##########################################Mineralizace dusiku##############################################
#Stechiometrickej pristup
zmeny_vse$Dekompozice<-with(zmeny_vse, CO2un+CUE2*CO2un/(1-CUE2))
zmeny_vse$dNH4_pred<-with(zmeny_vse, Dekompozice*(DON/DOC-(Nmic/0.54)*CUE2/(Cmic/0.38)))

zmeny_vse %>% filter(znaceni=="Unlabelled") %>% group_by(Plesne, horizont) %>% 
  summarize(y.sd=sd(dNH4), y=mean(dNH4)) %>%
  ggplot(aes(factor(Plesne), y))+geom_point(cex=6, pch=21, fill="grey")+
  facet_wrap(~horizont, scales="free")+theme_min+theme(legend.position = c(0.2, 0.2))+
  geom_errorbar(aes(ymin=y-y.sd, ymax=y+y.sd), width=0.1, lwd=0.5)+
  theme(legend.position = c(0.15, 0.2),legend.key.size = unit(0.3, "in"),
        legend.title = element_blank())+
  ylab(expression(paste(Delta~NH[4], " ( ", mu, "mol ",g^{-1}, ")" )))+
  xlab("Plesne : Certovo mixing ratio")+scale_fill_manual(values = c("white", "grey"))+
  geom_hline(yintercept = 0, lwd=1.5)+theme(legend.title = element_blank())

zmeny_vse %>% filter(znaceni=="Unlabelled") %>% group_by(Plesne, horizont) %>% 
  summarize(y.sd=sd(dNH4_pred), y=mean(dNH4_pred)) %>%
  ggplot(aes(factor(Plesne), y))+geom_point(cex=6, pch=21, fill="grey")+
  facet_wrap(~horizont, scales="free")+theme_min+theme(legend.position = c(0.2, 0.2))+
  geom_errorbar(aes(ymin=y-y.sd, ymax=y+y.sd), width=0.1, lwd=0.5)+
  theme(legend.position = c(0.15, 0.2),legend.key.size = unit(0.3, "in"),
        legend.title = element_blank())+
  ylab(expression(paste(Delta~NH[4], " ( ", mu, "mol ",g^{-1}, ")" )))+
  xlab("Plesne : Certovo mixing ratio")+scale_fill_manual(values = c("white", "grey"))+
  geom_hline(yintercept = 0, lwd=1.5)+theme(legend.title = element_blank())

#Urcim z izotopu proporci biomasy v CO2
zmeny_vse$CO2atm<-michizo[(michizo$znaceni=="Unlabelled" & michizo$Legend=="After incubation"), "CO2atm"]
zmeny_vse$deltaCO2<-with(zmeny_vse, (CO2atm/(1-CO2atm)/0.011237-1)*1000)
zmeny_vse$Cmicatm<-michizo[(michizo$znaceni=="Unlabelled" & michizo$Legend=="After incubation"), "Cmicatm"]
zmeny_vse$deltaCmic<-with(zmeny_vse, (Cmicatm/(1-Cmicatm)/0.011237-1)*1000)
zmeny_vse$deltaBase<-mean(zmeny_vse[(zmeny_vse$horizont=="O" & zmeny_vse$Plesne==0), "deltaCO2"], na.rm=T)
zmeny_vse$Cmic.prop<-with(zmeny_vse, (deltaCO2-deltaBase)/(deltaCmic-deltaBase))

zmeny_vse %>% filter(znaceni=="Unlabelled") %>% group_by(Plesne, horizont) %>% 
  summarize(y.sd=sd(Cmic.prop), y=mean(Cmic.prop)) %>%
  ggplot(aes(factor(Plesne), y))+geom_point(cex=6, pch=21, fill="grey")+
  facet_wrap(~horizont, scales="free")+theme_min+theme(legend.position = c(0.2, 0.2))+
  geom_errorbar(aes(ymin=y-y.sd, ymax=y+y.sd), width=0.1, lwd=0.5)+
  theme(legend.position = c(0.15, 0.2),legend.key.size = unit(0.3, "in"),
        legend.title = element_blank())+
  xlab("Plesne : Certovo mixing ratio")+scale_fill_manual(values = c("white", "grey"))+
  geom_hline(yintercept = 0, lwd=1.5)+theme(legend.title = element_blank())


#Predpokladam 0 u Acek
zmeny_vse[(zmeny_vse$znaceni=="Unlabelled" & zmeny_vse$horizont=="A"), "Cmic.prop"]<-0
#A nic nesmi byt vetsi nez 0
zmeny_vse[(zmeny_vse$Cmic.prop>1), "Cmic.prop"]<-1
zmeny_vse[(zmeny_vse$Cmic.prop<0), "Cmic.prop"]<-0

zmeny_vse %>% filter(znaceni=="Unlabelled") %>% group_by(Plesne, horizont) %>% 
  summarize(y.sd=sd(Cmic.prop), y=mean(Cmic.prop)) %>%
  ggplot(aes(factor(Plesne), y))+geom_point(cex=6, pch=21, fill="grey")+
  facet_wrap(~horizont, scales="free")+theme_min+theme(legend.position = c(0.2, 0.2))+
  geom_errorbar(aes(ymin=y-y.sd, ymax=y+y.sd), width=0.1, lwd=0.5)+
  theme(legend.position = c(0.15, 0.2),legend.key.size = unit(0.3, "in"),
        legend.title = element_blank())+
  xlab("Plesne : Certovo mixing ratio")+scale_fill_manual(values = c("white", "grey"))+
  geom_hline(yintercept = 0, lwd=1.5)+theme(legend.title = element_blank())

#Opravim odhad mineralizace N
zmeny_vse$dNH4_pred2<-with(zmeny_vse, Dekompozice*(((1-Cmic.prop)*(DON/DOC)+Cmic.prop*(Nmic*0.38/Cmic/0.54))-Nmic*CUE2*0.38/Cmic/0.54))
zmeny_vse$dNH4_pred3<-with(zmeny_vse, Dekompozice*(((1-Cmic.prop)*(DONb/DOCb)+Cmic.prop*(Nmic*0.38/Cmic/0.54))-Nmic*CUE2*0.38/Cmic/0.54))

grid.arrange(
  zmeny_vse %>% filter(znaceni=="Unlabelled") %>% group_by(Plesne, Horizon) %>% 
  summarize(y.sd=sd(dNH4_pred), y=mean(dNH4_pred),
            y2.sd=sd(dNH4+dNO3), y2=mean(dNH4+dNO3)) %>%
  ggplot(aes(Plesne, y))+geom_point(cex=6, pch=21, fill="dodgerblue3")+
  facet_grid(.~Horizon, scales="free")+theme_min+theme(legend.position = c(0.2, 0.2))+
  geom_errorbar(aes(ymin=y-y.sd, ymax=y+y.sd), width=0.1, lwd=0.5)+
  theme(legend.position = c(0.15, 0.2),legend.key.size = unit(0.3, "in"),
        legend.title = element_blank(),
        axis.title.x = element_blank())+
  ylab(expression(paste(Delta, " Mineral nitrogen (", mu, "mol ",g^{-1}, ")" )))+
  xlab("Plesne : Certovo mixing ratio")+scale_fill_manual(values = c("grey", "white"))+
  geom_hline(yintercept = 0, lwd=1)+theme(legend.title = element_blank())+
  geom_point(cex=6, pch=21, aes(Plesne, y2, fill=Horizon), show.legend = F)+
  geom_errorbar(aes(ymin=y2-y2.sd, ymax=y2+y2.sd), width=0.1, lwd=0.5)+ylim(-6.5, 1.5)+
  ggtitle("A)"),
  zmeny_vse %>% filter(znaceni=="Unlabelled") %>% group_by(Plesne, Horizon) %>% 
  summarize(y.sd=sd(dNH4_pred2), y=mean(dNH4_pred2),
            y2.sd=sd(dNH4+dNO3), y2=mean(dNH4+dNO3)) %>%
  ggplot(aes(Plesne, y))+geom_point(cex=6, pch=21, fill="dodgerblue3")+
  facet_grid(.~Horizon, scales="free")+theme_min+theme(legend.position = c(0.2, 0.2))+
  geom_errorbar(aes(ymin=y-y.sd, ymax=y+y.sd), width=0.1, lwd=0.5)+
  theme(legend.position = c(0.15, 0.2),legend.key.size = unit(0.3, "in"),
        legend.title = element_blank())+
  ylab(expression(paste(Delta, " Mineral nitrogen (", mu, "mol ",g^{-1}, ")" )))+
  xlab("Plesne : Certovo mixing ratio")+scale_fill_manual(values = c("grey", "white"))+
  geom_hline(yintercept = 0, lwd=1)+theme(legend.title = element_blank())+
  geom_point(cex=6, pch=21, aes(Plesne, y2, fill=Horizon), show.legend = F)+
  geom_errorbar(aes(ymin=y2-y2.sd, ymax=y2+y2.sd), width=0.1, lwd=0.5)+ylim(-6.5, 1.5)+
  ggtitle("B)"), nrow=2)

zmeny_vse %>% filter(znaceni=="Unlabelled") %>% group_by(Plesne, horizont) %>% 
  summarize(y.sd=sd(dNH4_pred2), y=mean(dNH4_pred2),
            y2.sd=sd(dNH4+dNO3), y2=mean(dNH4+dNO3)) %>%
  ggplot(aes(factor(Plesne), y))+geom_point(cex=6, pch=21, fill="grey")+
  facet_wrap(~horizont, scales="free")+theme_min+theme(legend.position = c(0.2, 0.2))+
  geom_errorbar(aes(ymin=y-y.sd, ymax=y+y.sd), width=0.1, lwd=0.5)+
  theme(legend.position = c(0.15, 0.2),legend.key.size = unit(0.3, "in"),
        legend.title = element_blank())+
  ylab(expression(paste(Delta~NH[4], " ( ", mu, "mol ",g^{-1}, ")" )))+
  xlab("Plesne : Certovo mixing ratio")+scale_fill_manual(values = c("white", "grey"))+
  geom_hline(yintercept = 0, lwd=1.5)+theme(legend.title = element_blank())+
  geom_point(cex=6, pch=21, aes(factor(Plesne), y2))+
  geom_errorbar(aes(ymin=y2-y2.sd, ymax=y2+y2.sd), width=0.1, lwd=0.5)


zmeny_vse %>% filter(znaceni=="Unlabelled" & horizont=="A") %>% group_by(Plesne) %>% 
  summarize(y=mean(dNH4+dNO3))


ggplot(zmeny_vse[zmeny_vse$znaceni=="Unlabelled", ], aes(dNH4+dNO3, dNH4_pred2))+geom_point(cex=6, pch=21, aes(fill=horizont))+
  theme_min+geom_abline(intercept = 0, slope = 1)+facet_wrap(~horizont, scales="free")
ggplot(zmeny_vse[zmeny_vse$znaceni=="Unlabelled", ], aes(dNH4, dNH4_pred))+geom_point(cex=6, pch=21, aes(fill=horizont))+
  theme_min+geom_abline(intercept = 0, slope = 1)

summary(lm(dNH4~dNH4_pred, zmeny_vse, subset = znaceni=="Unlabelled"))
summary(lm(dNH4~dNH4_pred2, zmeny_vse, subset = znaceni=="Unlabelled"))


#Cim je danej nesouhlas
zmeny_vse$dNH4_pred2_res<-with(zmeny_vse, (dNH4_pred2-dNH4)^2)

zmeny_vse2<-zmeny_vse
zmeny_vse2$dNH4_pred0.1<-with(zmeny_vse2, (CO2un+CO2un*(CUE2-0.1)/(1-CUE2-0.1))*(((1-Cmic.prop)*(DON/DOC)+Cmic.prop*(Nmic*0.38/Cmic/0.54))-Nmic*(CUE2-0.1)*0.38/Cmic/0.54))

with(subset(zmeny_vse2, Horizon=="Organic topsoil" & znaceni=="Unlabelled"), dNH4_pred2-dNH4_pred0.1)
##########################################Mineralizace fosforu##############################################
#Stechiometrickej pristup
zmeny_vse$dPO4_pred<-with(zmeny_vse, Dekompozice*(DOP/DOC-(Pmic/0.4)*CUE2/(Cmic/0.38)))

zmeny_vse %>% filter(znaceni=="Unlabelled") %>% group_by(Plesne, horizont) %>% 
  summarize(y.sd=sd(dPO42), y=mean(dPO42)) %>%
  ggplot(aes(factor(Plesne), y))+geom_point(cex=6, pch=21, fill="grey")+
  facet_wrap(~horizont, scales="free")+theme_min+theme(legend.position = c(0.2, 0.2))+
  geom_errorbar(aes(ymin=y-y.sd, ymax=y+y.sd), width=0.1, lwd=0.5)+
  theme(legend.position = c(0.15, 0.2),legend.key.size = unit(0.3, "in"),
        legend.title = element_blank())+
  ylab(expression(paste(Delta~PO[4], " ( ", mu, "mol ",g^{-1}, ")" )))+
  xlab("Plesne : Certovo mixing ratio")+scale_fill_manual(values = c("white", "grey"))+
  geom_hline(yintercept = 0, lwd=1.5)+theme(legend.title = element_blank())

zmeny_vse %>% filter(znaceni=="Unlabelled") %>% group_by(Plesne, horizont) %>% 
  summarize(y.sd=sd(dPO4_pred), y=mean(dPO4_pred)) %>%
  ggplot(aes(factor(Plesne), y))+geom_point(cex=6, pch=21, fill="grey")+
  facet_wrap(~horizont, scales="free")+theme_min+theme(legend.position = c(0.2, 0.2))+
  geom_errorbar(aes(ymin=y-y.sd, ymax=y+y.sd), width=0.1, lwd=0.5)+
  theme(legend.position = c(0.15, 0.2),legend.key.size = unit(0.3, "in"),
        legend.title = element_blank())+
  ylab(expression(paste(Delta~PO[4], " ( ", mu, "mol ",g^{-1}, ")" )))+
  xlab("Plesne : Certovo mixing ratio")+scale_fill_manual(values = c("white", "grey"))+
  geom_hline(yintercept = 0, lwd=1.5)+theme(legend.title = element_blank())

zmeny_vse %>% filter(znaceni=="Unlabelled") %>% group_by(Plesne, Horizon) %>% 
  summarize(y.sd=sd(dPO4_pred), y=mean(dPO4_pred),
            y2.sd=sd(dPO4), y2=mean(dPO4)) %>%
  ggplot(aes(Plesne, y))+geom_point(cex=6, pch=21, fill="dodgerblue3")+
  facet_grid(.~Horizon, scales="free")+theme_min+theme(legend.position = c(0.2, 0.2))+
  geom_errorbar(aes(ymin=y-y.sd, ymax=y+y.sd), width=0.1, lwd=0.5)+
  theme(legend.position = c(0.15, 0.2),legend.key.size = unit(0.3, "in"),
        legend.title = element_blank(),
        axis.title.x = element_blank())+
  ylab(expression(paste(Delta, " Mineral phosphorus (", mu, "mol ",g^{-1}, ")" )))+
  xlab("Plesne : Certovo mixing ratio")+scale_fill_manual(values = c("grey", "white"))+
  geom_hline(yintercept = 0, lwd=1)+theme(legend.title = element_blank())+
  geom_point(cex=6, pch=21, aes(Plesne, y2, fill=Horizon), show.legend = F)+
  geom_errorbar(aes(ymin=y2-y2.sd, ymax=y2+y2.sd), width=0.1, lwd=0.5)+
  ggtitle("A)")


#Opravim odhad mineralizace P
zmeny_vse$dPO4_pred2<-with(zmeny_vse, Dekompozice*(((1-Cmic.prop)*(DOP/DOC)+Cmic.prop*(Pmic*0.38/Cmic/0.4))-Pmic*CUE2*0.38/Cmic/0.4))

grid.arrange(
  zmeny_vse %>% filter(znaceni=="Unlabelled") %>% group_by(Plesne, Horizon) %>% 
    summarize(y.sd=sd(dPO4_pred), y=mean(dPO4_pred),
              y2.sd=sd(dPO4), y2=mean(dPO4)) %>%
    ggplot(aes(Plesne, y))+geom_point(cex=6, pch=21, fill="dodgerblue3")+
    facet_grid(.~Horizon, scales="free")+theme_min+theme(legend.position = c(0.2, 0.2))+
    geom_errorbar(aes(ymin=y-y.sd, ymax=y+y.sd), width=0.1, lwd=0.5)+
    theme(legend.position = c(0.15, 0.2),legend.key.size = unit(0.3, "in"),
          legend.title = element_blank(),
          axis.title.x = element_blank())+
    ylab(expression(paste(Delta, " Mineral phosphorus (", mu, "mol ",g^{-1}, ")" )))+
    xlab("Plesne : Certovo mixing ratio")+scale_fill_manual(values = c("grey", "white"))+
    geom_hline(yintercept = 0, lwd=1)+theme(legend.title = element_blank())+
    geom_point(cex=6, pch=21, aes(Plesne, y2, fill=Horizon), show.legend = F)+
    geom_errorbar(aes(ymin=y2-y2.sd, ymax=y2+y2.sd), width=0.1, lwd=0.5)+
    ggtitle("A)"),
  zmeny_vse %>% filter(znaceni=="Unlabelled") %>% group_by(Plesne, Horizon) %>% 
  summarize(y.sd=sd(dPO4_pred2), y=mean(dPO4_pred2),
            y2.sd=sd(dPO4), y2=mean(dPO4)) %>%
  ggplot(aes(Plesne, y))+geom_point(cex=6, pch=21, fill="dodgerblue3")+
  facet_grid(.~Horizon, scales="free")+theme_min+theme(legend.position = c(0.2, 0.2))+
  geom_errorbar(aes(ymin=y-y.sd, ymax=y+y.sd), width=0.1, lwd=0.5)+
  theme(legend.position = c(0.15, 0.2),legend.key.size = unit(0.3, "in"),
        legend.title = element_blank(),
        axis.title.x = element_blank())+
  ylab(expression(paste(Delta, " Mineral phosphorus (", mu, "mol ",g^{-1}, ")" )))+
  xlab("Plesne : Certovo mixing ratio")+scale_fill_manual(values = c("grey", "white"))+
  geom_hline(yintercept = 0, lwd=1)+theme(legend.title = element_blank())+
  geom_point(cex=6, pch=21, aes(Plesne, y2, fill=Horizon), show.legend = F)+
  geom_errorbar(aes(ymin=y2-y2.sd, ymax=y2+y2.sd), width=0.1, lwd=0.5)+
  ggtitle("B)"),
  zmeny_vse %>% filter(znaceni=="Unlabelled") %>% group_by(Plesne, Horizon) %>% 
    summarize(y.sd=sd(-PO4sorbed), y=mean(-PO4sorbed),
              y2.sd=sd(dPO4), y2=mean(dPO4)) %>%
    ggplot(aes(Plesne, y))+geom_point(cex=6, pch=21, fill="dodgerblue3")+
    facet_grid(.~Horizon, scales="free")+theme_min+theme(legend.position = c(0.2, 0.2))+
    geom_errorbar(aes(ymin=y-y.sd, ymax=y+y.sd), width=0.1, lwd=0.5)+
    theme(legend.position = c(0.15, 0.2),legend.key.size = unit(0.3, "in"))+
    ylab(expression(paste(Delta, " Mineral phosphorus (", mu, "mol ",g^{-1}, ")" )))+
    xlab("Plesne : Certovo mixing ratio")+scale_fill_manual(values = c("grey", "white"))+
    geom_hline(yintercept = 0, lwd=1)+theme(legend.title = element_blank())+
    geom_point(cex=6, pch=21, aes(Plesne, y2, fill=Horizon), show.legend = F)+
    geom_errorbar(aes(ymin=y2-y2.sd, ymax=y2+y2.sd), width=0.1, lwd=0.5)+
    ggtitle("C)"), ncol=1)


#Zkusim bez prispevku polyphosphatu
FEdata<-read.delim("clipboard")
ggplot(FEdata, aes(CP, u0))+geom_point(cex=6, pch=21, aes(fill=Position))+
  theme_min

summary(lm(CP~u0, FEdata, subset = Position=="Top"))

#Opravim odhad mineralizace P
zmeny_vse$dPO4_pred3<-with(zmeny_vse, Dekompozice*((DOP/DOC)-1/75.76))

zmeny_vse %>% filter(znaceni=="Unlabelled") %>% group_by(Plesne, Horizon) %>% 
  summarize(y.sd=sd(dPO4_pred3), y=mean(dPO4_pred3),
            y2.sd=sd(dPO4), y2=mean(dPO4)) %>%
  ggplot(aes(Plesne, y))+geom_point(cex=6, pch=21, fill="dodgerblue3")+
  facet_grid(.~Horizon, scales="free")+theme_min+theme(legend.position = c(0.2, 0.2))+
  geom_errorbar(aes(ymin=y-y.sd, ymax=y+y.sd), width=0.1, lwd=0.5)+
  theme(legend.position = c(0.15, 0.2),legend.key.size = unit(0.3, "in"),
        legend.title = element_blank(),
        axis.title.x = element_blank())+
  ylab(expression(paste(Delta, " Mineral phosphorus (", mu, "mol ",g^{-1}, ")" )))+
  xlab("Plesne : Certovo mixing ratio")+scale_fill_manual(values = c("grey", "white"))+
  geom_hline(yintercept = 0, lwd=1)+theme(legend.title = element_blank())+
  geom_point(cex=6, pch=21, aes(Plesne, y2, fill=Horizon), show.legend = F)+
  geom_errorbar(aes(ymin=y2-y2.sd, ymax=y2+y2.sd), width=0.1, lwd=0.5)+
  ggtitle("C)")

#Opravim odhad mineralizace P
zmeny_vse$CPcr<-mich2[c(1:80), "CPcr"]
zmeny_vse$dPO4_pred4<-with(zmeny_vse, Dekompozice*((DOP/DOC)-1/CPcr))

zmeny_vse %>% filter(znaceni=="Unlabelled") %>% group_by(Plesne, Horizon) %>% 
  summarize(y.sd=sd(dPO4_pred4), y=mean(dPO4_pred4),
            y2.sd=sd(dPO4), y2=mean(dPO4)) %>%
  ggplot(aes(Plesne, y))+geom_point(cex=6, pch=21, fill="dodgerblue3")+
  facet_grid(.~Horizon, scales="free")+theme_min+theme(legend.position = c(0.2, 0.2))+
  geom_errorbar(aes(ymin=y-y.sd, ymax=y+y.sd), width=0.1, lwd=0.5)+
  theme(legend.position = c(0.15, 0.2),legend.key.size = unit(0.3, "in"),
        legend.title = element_blank(),
        axis.title.x = element_blank())+
  ylab(expression(paste(Delta, " Mineral phosphorus (", mu, "mol ",g^{-1}, ")" )))+
  xlab("Plesne : Certovo mixing ratio")+scale_fill_manual(values = c("grey", "white"))+
  geom_hline(yintercept = 0, lwd=1)+theme(legend.title = element_blank())+
  geom_point(cex=6, pch=21, aes(Plesne, y2, fill=Horizon), show.legend = F)+
  geom_errorbar(aes(ymin=y2-y2.sd, ymax=y2+y2.sd), width=0.1, lwd=0.5)+
  ggtitle("D)")

#Mineralizace organickyho fosforu
library(lamW)
zmeny_vse$dDOP_pred<-with(zmeny_vse, phosKm*lambertW0(DOP/phosKm*exp(DOP/phosKm-phosVmax*48))-DOP)

zmeny_vse %>% filter(znaceni=="Unlabelled") %>% group_by(Plesne, horizont) %>% 
  summarize(y.sd=sd(dDOP_pred), y=mean(dDOP_pred),
            y2.sd=sd(dDOP), y2=mean(dDOP)) %>%
  ggplot(aes(factor(Plesne), y))+geom_point(cex=6, pch=21, fill="grey")+
  facet_wrap(~horizont, scales="free")+theme_min+theme(legend.position = c(0.2, 0.2))+
  geom_errorbar(aes(ymin=y-y.sd, ymax=y+y.sd), width=0.1, lwd=0.5)+
  theme(legend.position = c(0.15, 0.2),legend.key.size = unit(0.3, "in"),
        legend.title = element_blank())+
  ylab(expression(paste(Delta~PO[4], " ( ", mu, "mol ",g^{-1}, ")" )))+
  xlab("Plesne : Certovo mixing ratio")+scale_fill_manual(values = c("white", "grey"))+
  geom_hline(yintercept = 0, lwd=1.5)+theme(legend.title = element_blank())+
  geom_point(cex=6, pch=21, aes(factor(Plesne), y2))+
  geom_errorbar(aes(ymin=y2-y2.sd, ymax=y2+y2.sd), width=0.1, lwd=0.5)

zmeny_vse$dPO4_pred6<-with(zmeny_vse, (DOP-phosKm*lambertW0(DOP/phosKm*exp(DOP/phosKm-phosVmax*48)))-PO4sorbed)

zmeny_vse %>% filter(znaceni=="Unlabelled") %>% group_by(Plesne, horizont) %>% 
  summarize(y.sd=sd(dPO4_pred6), y=mean(dPO4_pred6),
            y2.sd=sd(dPO4), y2=mean(dPO4)) %>%
  ggplot(aes(factor(Plesne), y))+geom_point(cex=6, pch=21, fill="grey")+
  facet_wrap(~horizont, scales="free")+theme_min+theme(legend.position = c(0.2, 0.2))+
  geom_errorbar(aes(ymin=y-y.sd, ymax=y+y.sd), width=0.1, lwd=0.5)+
  theme(legend.position = c(0.15, 0.2),legend.key.size = unit(0.3, "in"),
        legend.title = element_blank())+
  ylab(expression(paste(Delta~PO[4], " ( ", mu, "mol ",g^{-1}, ")" )))+
  xlab("Plesne : Certovo mixing ratio")+scale_fill_manual(values = c("white", "grey"))+
  geom_hline(yintercept = 0, lwd=1.5)+theme(legend.title = element_blank())+
  geom_point(cex=6, pch=21, aes(factor(Plesne), y2))+
  geom_errorbar(aes(ymin=y2-y2.sd, ymax=y2+y2.sd), width=0.1, lwd=0.5)


#mineralizace organickyho fosforu sofistikovanejc
#definuju model
Pmodel<-function(time, state, pars){
  with(as.list(c(state, pars)),{
    dPON<--Vmax*PON/(Km+PON)+Pmic*k
    dPO4<-Vmax*PON/(Km+PON)-PO4*Sorption
    dPmic<--Pmic*k
    return(list(c(dPON, dPO4, dPmic)))
  })
}


zmeny_vse$dPO4_pred5<-vector("numeric", length = 80)

for(i in 1:80){
  zmeny_vse$dPO4_pred5[i]<-
    (as.data.frame(ode(y=c(PON=zmeny_vse$DOP[i], PO4=zmeny_vse$PO4[i], Pmic = zmeny_vse$Pmic[i]/0.4), 
                    parms=c(Vmax = zmeny_vse$phosVmax[i],
                            Km = zmeny_vse$phosVmax[i],
                            Sorption = zmeny_vse$Sorption[i]/0.75,
                            k = -zmeny_vse$dPmic[i]/0.4/48), 
                    Pmodel, times=seq(0,48)))[49, 3])-zmeny_vse$PO4[i]
}

zmeny_vse %>% filter(znaceni=="Unlabelled" & dPO4_pred5>-0.15) %>% group_by(Plesne, horizont) %>% 
  summarize(y.sd=sd(dPO4_pred5), y=mean(dPO4_pred5),
            y2.sd=sd(dPO4), y2=mean(dPO4)) %>%
  ggplot(aes(factor(Plesne), y))+geom_point(cex=6, pch=21, fill="grey")+
  facet_wrap(~horizont, scales="free")+theme_min+theme(legend.position = c(0.2, 0.2))+
  geom_errorbar(aes(ymin=y-y.sd, ymax=y+y.sd), width=0.1, lwd=0.5)+
  theme(legend.position = c(0.15, 0.2),legend.key.size = unit(0.3, "in"),
        legend.title = element_blank())+
  ylab(expression(paste(Delta~PO[4], " ( ", mu, "mol ",g^{-1}, ")" )))+
  xlab("Plesne : Certovo mixing ratio")+scale_fill_manual(values = c("white", "grey"))+
  geom_hline(yintercept = 0, lwd=1.5)+theme(legend.title = element_blank())+
  geom_point(cex=6, pch=21, aes(factor(Plesne), y2))+
  geom_errorbar(aes(ymin=y2-y2.sd, ymax=y2+y2.sd), width=0.1, lwd=0.5)

ggplot(zmeny_vse[zmeny_vse$znaceni=="Unlabelled", ], aes(dPO4_pred5, dPO4))+
  geom_point(cex=6, pch=21, aes(fill=horizont))+theme_min+geom_abline(intercept = 0, slope = 1)

zmeny_vse %>% filter(znaceni=="Unlabelled") %>% group_by(Plesne, horizont) %>% 
  summarize(y.sd=sd(-PO4sorbed), y=mean(-PO4sorbed),
            y2.sd=sd(dPO4), y2=mean(dPO4)) %>%
  ggplot(aes(factor(Plesne), y))+geom_point(cex=6, pch=21, fill="grey")+
  facet_wrap(~horizont, scales="free")+theme_min+theme(legend.position = c(0.2, 0.2))+
  geom_errorbar(aes(ymin=y-y.sd, ymax=y+y.sd), width=0.1, lwd=0.5)+
  theme(legend.position = c(0.15, 0.2),legend.key.size = unit(0.3, "in"),
        legend.title = element_blank())+
  ylab(expression(paste(Delta~PO[4], " ( ", mu, "mol ",g^{-1}, ")" )))+
  xlab("Plesne : Certovo mixing ratio")+scale_fill_manual(values = c("white", "grey"))+
  geom_hline(yintercept = 0, lwd=1.5)+theme(legend.title = element_blank())+
  geom_point(cex=6, pch=21, aes(factor(Plesne), y2))+
  geom_errorbar(aes(ymin=y2-y2.sd, ymax=y2+y2.sd), width=0.1, lwd=0.5)

zmeny_vse$dPO4_pred6<-with(zmeny_vse, -dPmic*phosVmax/(-dPmic+phosKm)*48-PO4sorbed)

zmeny_vse %>% filter(znaceni=="Unlabelled") %>% group_by(Plesne, horizont) %>% 
  summarize(y.sd=sd(-PO4sorbed), y=mean(-PO4sorbed),
            y2.sd=sd(dPO4), y2=mean(dPO4)) %>%
  ggplot(aes(factor(Plesne), y))+geom_point(cex=6, pch=21, fill="grey")+
  facet_wrap(~horizont, scales="free")+theme_min+theme(legend.position = c(0.2, 0.2))+
  geom_errorbar(aes(ymin=y-y.sd, ymax=y+y.sd), width=0.1, lwd=0.5)+
  theme(legend.position = c(0.15, 0.2),legend.key.size = unit(0.3, "in"),
        legend.title = element_blank())+
  ylab(expression(paste(Delta~PO[4], " ( ", mu, "mol ",g^{-1}, ")" )))+
  xlab("Plesne : Certovo mixing ratio")+scale_fill_manual(values = c("white", "grey"))+
  geom_hline(yintercept = 0, lwd=1.5)+theme(legend.title = element_blank())+
  geom_point(cex=6, pch=21, aes(factor(Plesne), y2))+
  geom_errorbar(aes(ymin=y2-y2.sd, ymax=y2+y2.sd), width=0.1, lwd=0.5)

############################################Statistiky#################################################
#Cmic init
anova(glm(Cmic~horizont*Plesne/horizont, data = mich2[(mich2$znaceni=="Unlabelled" &
                                                  mich2$cas==0), ], family = Gamma),
      test="F")

anova(glm(Cmic~horizont+Plesne/horizont+cas/Plesne/horizont, data = mich2[(mich2$znaceni=="Unlabelled" ), ], family = Gamma),
      test="F")

#Nmic init
anova(glm(Nmic~horizont*Plesne/horizont, data = mich2[(mich2$znaceni=="Unlabelled" &
                                                mich2$cas==0), ], family = Gamma),
      test="F")

anova(glm(Nmic~horizont+Plesne/horizont+cas/Plesne/horizont, data = mich2[(mich2$znaceni=="Unlabelled"), ], family = Gamma),
      test="F")

#Pmic init
anova(glm(Pmic~horizont*Plesne/horizont, data = mich2[(mich2$znaceni=="Unlabelled" &
                                                         mich2$cas==0), ], family = Gamma),
      test="F")

anova(glm(Pmic~horizont+Plesne/horizont+cas/Plesne/horizont, data = mich2[(mich2$znaceni=="Unlabelled" ), ], family = Gamma),
      test="F")

anova(glm(PO4~horizont+Plesne/horizont+cas/Plesne/horizont, data = mich2[(mich2$znaceni=="Unlabelled" ), ], family = Gamma),
      test="F")

#DOC init
anova(glm(DOC~horizont*Plesne/horizont, data = mich2[(mich2$znaceni=="Unlabelled" &
                                                mich2$cas==0), ], family = Gamma),
      test="F")
#NH4 init
anova(glm(NH4~horizont*Plesne/horizont, data = mich2[(mich2$znaceni=="Unlabelled" &
                                               mich2$cas==0), ], family = Gamma),
      test="F")

#DON init
anova(glm(DON~horizont*Plesne/horizont, data = mich2[(mich2$znaceni=="Unlabelled" &
                                               mich2$cas==0), ], family = Gamma),
      test="F")

#DOP init
anova(glm(DOP~horizont*Plesne/horizont, data = mich2[(mich2$znaceni=="Unlabelled" &
                                                        mich2$cas==0), ], family = Gamma),
      test="F")

#NO3 init
anova(glm(NO3~horizont*Plesne/horizont, data = mich2[(mich2$znaceni=="Unlabelled" &
                                               mich2$cas==0), ], family = Gamma),
      test="F")

#CNmic init
mich2$CNmic<-with(mich2, Cmic*0.54/Nmic/0.38)
anova(glm(CNmic~horizont*Plesne/horizont, data = mich2[(mich2$znaceni=="Unlabelled" &
                                               mich2$cas==0), ], family = Gamma),
      test="F")

#CPmic init
mich2$CPmic<-with(mich2, Cmic*0.4/Pmic/0.38)
anova(glm(CPmic~horizont*Plesne/horizont, data = mich2[(mich2$znaceni=="Unlabelled" &
                                                          mich2$cas==0), ], family = Gamma),
      test="F")

#CPmic init
anova(glm(DOC/DOP~horizont*Plesne/horizont, data = mich2[(mich2$znaceni=="Unlabelled" &
                                                              mich2$cas==0), ], family = Gamma),
      test="F")

#DOC/N init
anova(glm(DOC/DON~horizont*Plesne/horizont, data = mich2[(mich2$znaceni=="Unlabelled" &
                                                 mich2$cas==0), ], family = Gamma),
      test="F")

write.csv(Obrat, (".obrat.csv"))

anova(glm((NH4+NO3)~horizont+Plesne/horizont+cas/Plesne/horizont, data = mich2[(mich2$znaceni=="Unlabelled" ), ], family = Gamma),
      test="F")

anova(glm((NH4+NO3)~Plesne+cas/Plesne, 
          data = mich2[(mich2$znaceni=="Unlabelled"), ], family = Gamma,
          subset = horizont=="O"),
      test="F")

anova(glm((NH4+NO3)~Plesne+cas/Plesne, 
          data = mich2[(mich2$znaceni=="Unlabelled"), ], family = Gamma,
          subset = horizont=="A"),
      test="F")

#13C
anova(glm(CO2atm~Plesne, data = mich2[(mich2$znaceni=="Unlabelled" & mich2$horizont=="O" & 
                                   mich2$cas==1), ], family = Gamma),
      test="F")

Aizo<-mich2[(mich2$znaceni=="Unlabelled" & mich2$horizont=="A" &
               mich2$cas==1), c("Plesne", "CO2atm", "DOCatm", "Cmicatm")]
Aizo<-melt(Aizo, id.vars = "Plesne")

anova(glm(value~Plesne+variable/Plesne, data = Aizo, family = Gamma),
     test="F")

summary(zmeny_vse[zmeny_vse$znaceni=="Unlabelled", "dNH4"]+
          zmeny_vse[zmeny_vse$znaceni=="Unlabelled", "dNO3"])

ggplot(zmeny_vse[(zmeny_vse$horizont=="O" & zmeny_vse$znaceni=="Unlabelled"), ],
       aes(dCmic, Cmicatm))+geom_point()
