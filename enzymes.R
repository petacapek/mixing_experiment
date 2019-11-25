##########################################enzymy#################################################
##beta-glucosidase
mich2 %>% group_by(Plesne, znaceni, horizont, Legend) %>% summarize(y.sd=sd(bg), y=mean(bg)) %>%
  ggplot(aes(factor(Plesne), y))+geom_point(pch=21, cex=6, aes(fill=Legend))+
  facet_grid(znaceni~horizont, scales="free")+theme_min+theme(legend.position = c(0.17, 0.88),
                                                              legend.key.size = unit(0.3, "in"),
                                                              legend.title = element_blank())+
  geom_errorbar(aes(ymin=y-y.sd, ymax=y+y.sd), width=0.1, lwd=0.5)+
  ylab(expression(paste(beta,"-glucosidase ( ",mu, "mol ",g^{-1}~h^{-1}, ")" )))+
  xlab("Plesne : Certovo mixing ratio")+scale_fill_manual(values = c("white", "grey"))

##cellobiosidase
mich2 %>% group_by(Plesne, znaceni, horizont, Legend) %>% summarize(y.sd=sd(cell), y=mean(cell)) %>%
  ggplot(aes(factor(Plesne), y))+geom_point(pch=21, cex=6, aes(fill=Legend))+
  facet_grid(znaceni~horizont, scales="free")+theme_min+theme(legend.position = c(0.17, 0.9),
                                                              legend.key.size = unit(0.3, "in"),
                                                              legend.title = element_blank())+
  geom_errorbar(aes(ymin=y-y.sd, ymax=y+y.sd), width=0.1, lwd=0.5)+
  ylab(expression(paste("cellobiosidase ( ",mu, "mol ",g^{-1}~h^{-1}, ")" )))+
  xlab("Plesne : Certovo mixing ratio")+scale_fill_manual(values = c("white", "grey"))

##alanin-aminopeptidase
mich2 %>% group_by(Plesne, znaceni, horizont, Legend) %>% summarize(y.sd=sd(ala), y=mean(ala)) %>%
  ggplot(aes(factor(Plesne), y))+geom_point(pch=21, cex=6, aes(fill=Legend))+
  facet_grid(znaceni~horizont, scales="free")+theme_min+theme(legend.position = c(0.2, 0.92),
                                                              legend.key.size = unit(0.3, "in"),
                                                              legend.title = element_blank())+
  geom_errorbar(aes(ymin=y-y.sd, ymax=y+y.sd), width=0.1, lwd=0.5)+
  ylab(expression(paste("alanin aminopeptidase ( ",mu, "mol ",g^{-1}~h^{-1}, ")" )))+
  xlab("Plesne : Certovo mixing ratio")+scale_fill_manual(values = c("white", "grey"))

##chitinase
mich2 %>% group_by(Plesne, znaceni, horizont, Legend) %>% summarize(y.sd=sd(chit), y=mean(chit)) %>%
  ggplot(aes(factor(Plesne), y))+geom_point(pch=21, cex=6, aes(fill=Legend))+
  facet_grid(znaceni~horizont, scales="free")+theme_min+theme(legend.position = c(0.8, 0.1),
                                                              legend.key.size = unit(0.3, "in"),
                                                              legend.title = element_blank())+
  geom_errorbar(aes(ymin=y-y.sd, ymax=y+y.sd), width=0.1, lwd=0.5)+
  ylab(expression(paste("chitinase ( ", mu, "mol ",g^{-1}~h^{-1}, ")" )))+
  xlab("Plesne : Certovo mixing ratio")+scale_fill_manual(values = c("white", "grey"))

##phosphatase
mich2 %>% group_by(Plesne, znaceni, horizont, Legend) %>% summarize(y.sd=sd(phos), y=mean(phos)) %>%
  ggplot(aes(factor(Plesne), y))+geom_point(pch=21, cex=6, aes(fill=Legend))+
  facet_grid(znaceni~horizont, scales="free")+theme_min+theme(legend.position = c(0.2, 0.9),
                                                              legend.key.size = unit(0.3, "in"),
                                                              legend.title = element_blank())+
  geom_errorbar(aes(ymin=y-y.sd, ymax=y+y.sd), width=0.1, lwd=0.5)+
  ylab(expression(paste("phosphatase ( ", mu, "mol ",g^{-1}~h^{-1}, ")" )))+
  xlab("Plesne : Certovo mixing ratio")+scale_fill_manual(values = c("white", "grey"))

##phenoloxidase
mich2 %>% group_by(Plesne, znaceni, horizont, Legend) %>% summarize(y.sd=sd(PHO), y=mean(PHO)) %>%
  ggplot(aes(factor(Plesne), y))+geom_point(pch=21, cex=6, aes(fill=Legend))+
  facet_grid(znaceni~horizont, scales="free")+theme_min+theme(legend.position = c(0.7, 0.9),
                                                              legend.key.size = unit(0.3, "in"),
                                                              legend.title = element_blank())+
  geom_errorbar(aes(ymin=y-y.sd, ymax=y+y.sd), width=0.1, lwd=0.5)+
  ylab(expression(paste("phenoloxidase ( ", mu, "mol ",g^{-1}~h^{-1}, ")" )))+
  xlab("Plesne : Certovo mixing ratio")+scale_fill_manual(values = c("white", "grey"))

##peroxidase
mich2 %>% group_by(Plesne, znaceni, horizont, Legend) %>% summarize(y.sd=sd(PER), y=mean(PER)) %>%
  ggplot(aes(factor(Plesne), y))+geom_point(pch=21, cex=6, aes(fill=Legend))+
  facet_grid(znaceni~horizont, scales="free")+theme_min+theme(legend.position = c(0.7, 0.9),
                                                              legend.key.size = unit(0.3, "in"),
                                                              legend.title = element_blank())+
  geom_errorbar(aes(ymin=y-y.sd, ymax=y+y.sd), width=0.1, lwd=0.5)+
  ylab(expression(paste("peroxidase ( ", mu, "mol ",g^{-1}~h^{-1}, ")" )))+
  xlab("Plesne : Certovo mixing ratio")+scale_fill_manual(values = c("white", "grey"))

##################################pomery####################################################
##Ec/En
mich2 %>% group_by(Plesne, znaceni, horizont, Legend) %>% summarize(y.sd=sd((bg+cell)/(ala+chit)), y=mean((BG+CELL)/(ALA+CHIT))) %>%
  ggplot(aes(factor(Plesne), y))+geom_point(pch=21, cex=6, aes(fill=Legend))+
  facet_grid(znaceni~horizont, scales="free")+theme_min+theme(legend.position = c(0.66, 0.95),
                                                              legend.key.size = unit(0.3, "in"),
                                                              legend.title = element_blank())+
  geom_errorbar(aes(ymin=y-y.sd, ymax=y+y.sd), width=0.1, lwd=0.5)+
  ylab(expression(paste(E[C]/E[N])))+
  xlab("Plesne : Certovo mixing ratio")+scale_fill_manual(values = c("white", "grey"))

##Ec/Ep
mich2 %>% group_by(Plesne, znaceni, horizont, Legend) %>% summarize(y.sd=sd((bg+cell)/(phos)), y=mean((BG+CELL)/(PH))) %>%
  ggplot(aes(factor(Plesne), y))+geom_point(pch=21, cex=6, aes(fill=Legend))+
  facet_grid(znaceni~horizont, scales="free")+theme_min+theme(legend.position = c(0.7, 0.93),
                                                              legend.key.size = unit(0.3, "in"),
                                                              legend.title = element_blank())+
  geom_errorbar(aes(ymin=y-y.sd, ymax=y+y.sd), width=0.1, lwd=0.5)+
  ylab(expression(paste(E[C]/E[P])))+
  xlab("Plesne : Certovo mixing ratio")+scale_fill_manual(values = c("white", "grey"))

##En/Ep
mich2 %>% group_by(Plesne, znaceni, horizont, Legend) %>% summarize(y.sd=sd((ala+chit)/(phos)), y=mean((ALA+CHIT)/(PH))) %>%
  ggplot(aes(factor(Plesne), y))+geom_point(pch=21, cex=6, aes(fill=Legend))+
  facet_grid(znaceni~horizont, scales="free")+theme_min+theme(legend.position = c(0.67, 0.9),
                                                              legend.key.size = unit(0.3, "in"),
                                                              legend.title = element_blank())+
  geom_errorbar(aes(ymin=y-y.sd, ymax=y+y.sd), width=0.1, lwd=0.5)+
  ylab(expression(paste(E[N]/E[P])))+
  xlab("Plesne : Certovo mixing ratio")+scale_fill_manual(values = c("white", "grey"))

##############################zmeny##########################################################
####zmeny
e.zmeny<-mich2[1:80, c("Plesne", "Certovo", "horizont", "znaceni")]
e.zmeny$dBG<-mich2[81:160, "bg"]-mich2[1:80, "bg"]
e.zmeny$dCELL<-mich2[81:160, "cell"]-mich2[1:80, "cell"]
e.zmeny$dALA<-mich2[81:160, "ala"]-mich2[1:80, "ala"]
e.zmeny$dCHIT<-mich2[81:160, "chit"]-mich2[1:80, "chit"]
e.zmeny$dPH<-mich2[81:160, "phos"]-mich2[1:80, "phos"]
e.zmeny$dPHO<-mich2[81:160, "PHO"]-mich2[1:80, "PHO"]
e.zmeny$dPER<-mich2[81:160, "PER"]-mich2[1:80, "PER"]

##BG
e.zmeny %>% group_by(Plesne, znaceni, horizont) %>% summarize(y.sd=sd(dBG), y=mean(dBG)) %>%
  ggplot(aes(factor(Plesne), y))+geom_point(cex=6, pch=21, aes(fill=znaceni))+
  facet_grid(.~horizont, scales="free")+theme_min+theme(legend.position = c(0.2, 0.2))+
  geom_errorbar(aes(ymin=y-y.sd, ymax=y+y.sd), width=0.1, lwd=0.5)+
  theme(legend.position = c(0.15, 0.2),legend.key.size = unit(0.3, "in"),
       legend.title = element_blank())+
  ylab(expression(paste(Delta~beta, "-glucosidase ( ", mu, "mol ",g^{-1}~h^{-1}, ")" )))+
  xlab("Plesne : Certovo mixing ratio")+scale_fill_manual(values = c("white", "grey"))+
  geom_hline(yintercept = 0, lwd=1.5)

##CELL
e.zmeny %>% group_by(Plesne, znaceni, horizont) %>% summarize(y.sd=sd(dCELL), y=mean(dCELL)) %>%
  ggplot(aes(factor(Plesne), y))+geom_point(cex=6, pch=21, aes(fill=znaceni))+
  facet_grid(.~horizont, scales="free")+theme_min+theme(legend.position = c(0.2, 0.2))+
  geom_errorbar(aes(ymin=y-y.sd, ymax=y+y.sd), width=0.1, lwd=0.5)+
  theme(legend.position = c(0.15, 0.2),legend.key.size = unit(0.3, "in"),
        legend.title = element_blank())+
  ylab(expression(paste(Delta, " cellobiosidase ( ", mu, "mol ",g^{-1}~h^{-1}, ")" )))+
  xlab("Plesne : Certovo mixing ratio")+scale_fill_manual(values = c("white", "grey"))+
  geom_hline(yintercept = 0, lwd=1.5)

##ALA
e.zmeny %>% group_by(Plesne, znaceni, horizont) %>% summarize(y.sd=sd(dALA), y=mean(dALA)) %>%
  ggplot(aes(factor(Plesne), y))+geom_point(cex=6, pch=21, aes(fill=znaceni))+
  facet_grid(.~horizont, scales="free")+theme_min+theme(legend.position = c(0.2, 0.2))+
  geom_errorbar(aes(ymin=y-y.sd, ymax=y+y.sd), width=0.1, lwd=0.5)+
  theme(legend.position = c(0.15, 0.2),legend.key.size = unit(0.3, "in"),
        legend.title = element_blank())+
  ylab(expression(paste(Delta, " alaninaminopeptidase ( ", mu, "mol ",g^{-1}~h^{-1}, ")" )))+
  xlab("Plesne : Certovo mixing ratio")+scale_fill_manual(values = c("white", "grey"))+
  geom_hline(yintercept = 0, lwd=1.5)+theme(legend.title = element_blank())

##CHIT
e.zmeny %>% group_by(Plesne, znaceni, horizont) %>% summarize(y.sd=sd(dCHIT), y=mean(dCHIT)) %>%
  ggplot(aes(factor(Plesne), y))+geom_point(cex=6, pch=21, aes(fill=znaceni))+
  facet_grid(.~horizont, scales="free")+theme_min+theme(legend.position = c(0.2, 0.2))+
  geom_errorbar(aes(ymin=y-y.sd, ymax=y+y.sd), width=0.1, lwd=0.5)+
  theme(legend.position = c(0.15, 0.2),legend.key.size = unit(0.3, "in"),
        legend.title = element_blank())+
  ylab(expression(paste(Delta, " leucinaminopeptidase ( ", mu, "mol ",g^{-1}~h^{-1}, ")" )))+
  xlab("Plesne : Certovo mixing ratio")+scale_fill_manual(values = c("white", "grey"))+
  geom_hline(yintercept = 0, lwd=1.5)+theme(legend.title = element_blank())

##PH
e.zmeny %>% group_by(Plesne, znaceni, horizont) %>% summarize(y.sd=sd(dPH), y=mean(dPH)) %>%
  ggplot(aes(factor(Plesne), y))+geom_point(cex=6, pch=21, aes(fill=znaceni))+
  facet_grid(.~horizont, scales="free")+theme_min+theme(legend.position = c(0.2, 0.2))+
  geom_errorbar(aes(ymin=y-y.sd, ymax=y+y.sd), width=0.1, lwd=0.5)+
  theme(legend.position = c(0.15, 0.2),legend.key.size = unit(0.3, "in"),
        legend.title = element_blank())+
  ylab(expression(paste(Delta, " phosphatase ( ", mu, "mol ",g^{-1}~h^{-1}, ")" )))+
  xlab("Plesne : Certovo mixing ratio")+scale_fill_manual(values = c("white", "grey"))+
  geom_hline(yintercept = 0, lwd=1.5)+theme(legend.title = element_blank())

######################################pomery s biomasou###################################
##Ec/Cmic
mich2 %>% group_by(Plesne, znaceni, horizont, Legend) %>% summarize(y.sd=sd((bg+cell)/Cmic), y=mean((bg+cell)/Cmic)) %>%
  ggplot(aes(factor(Plesne), y))+geom_point(pch=21, cex=6, aes(fill=Legend))+
  facet_grid(znaceni~horizont, scales="free")+theme_min+theme(legend.position = c(0.7, 0.9),
                                                              legend.key.size = unit(0.3, "in"),
                                                              legend.title = element_blank())+
  geom_errorbar(aes(ymin=y-y.sd, ymax=y+y.sd), width=0.1, lwd=0.5)+
  ylab(expression(paste(E[C]/C[MB])))+
  xlab("Plesne : Certovo mixing ratio")+scale_fill_manual(values = c("white", "grey"))

##En/Nmic
mich2 %>% group_by(Plesne, znaceni, horizont, Legend) %>% summarize(y.sd=sd((ala+chit)/Nmic), y=mean((ala+chit)/Nmic)) %>%
  ggplot(aes(factor(Plesne), y))+geom_point(pch=21, cex=6, aes(fill=Legend))+
  facet_grid(znaceni~horizont, scales="free")+theme_min+theme(legend.position = c(0.7, 0.9),
                                                              legend.key.size = unit(0.3, "in"),
                                                              legend.title = element_blank())+
  geom_errorbar(aes(ymin=y-y.sd, ymax=y+y.sd), width=0.1, lwd=0.5)+
  ylab(expression(paste(E[N]/N[MB])))+
  xlab("Plesne : Certovo mixing ratio")+scale_fill_manual(values = c("white", "grey"))

##Ep/Pmic
mich2 %>% group_by(Plesne, znaceni, horizont, Legend) %>% summarize(y.sd=sd((phos)/Pmic), y=mean((phos)/Pmic)) %>%
  ggplot(aes(factor(Plesne), y))+geom_point(pch=21, cex=6, aes(fill=Legend))+
  facet_grid(znaceni~horizont, scales="free")+theme_min+theme(legend.position = c(0.25, 0.9),
                                                              legend.key.size = unit(0.3, "in"),
                                                              legend.title = element_blank())+
  geom_errorbar(aes(ymin=y-y.sd, ymax=y+y.sd), width=0.1, lwd=0.5)+
  ylab(expression(paste(E[P]/P[MB])))+
  xlab("Plesne : Certovo mixing ratio")+scale_fill_manual(values = c("white", "grey"))

#####################################Duggleby 1995################################################
source("enzyme_calculation_function.R")
library(lmodel2)

Sconc<-read.xlsx(xlsxFile = c("./Enzymes_data/MUF/Sconcentrations.xlsx"))

#Files 
fenz<-list.files(path = "./Enzymes_data/MUF/")[-c(46:48)]
fenz<-paste0("./Enzymes_data/MUF/", fenz)

#First file
Enz<-enzymes_calc(filepath = fenz[1], Sconc = Sconc)

#The rest of the files
for(i in 2:length(fenz)){
  Enz<-rbind(Enz, enzymes_calc(filepath = fenz[i], Sconc = Sconc))
}

Enz2<-Enz

Enz2[Enz2$gluVmax<0, "gluVmax"]<-0
Enz2[Enz2$celVmax<0, "celVmax"]<-0
Enz2[Enz2$phosVmax<0, "phosVmax"]<-0
Enz2[Enz2$leuVmax<0, "leuVmax"]<-0
Enz2[Enz2$chitVmax<0, "chitVmax"]<-0

Enz$Sample<-as.factor(Enz$Sample)
Enz$Sample<-factor(Enz$Sample, levels = unique(Sconc$Sample))
Enz<-Enz[order(Enz$Sample), ]

Enz_complete<-rbind(Enz2[1:40, ], Enz2[1:40, ], Enz2[41:120, ])

mich2e<-cbind(mich2, Enz_complete[,-1])

###Plot Vmaxes
##beta-glucosidase
mich2e %>% group_by(Plesne, znaceni, horizont, Legend) %>% summarize(y.sd=sd(gluVmax, na.rm = T), 
                                                                     y=mean(gluVmax, na.rm = T)) %>%
  ggplot(aes(factor(Plesne), y))+geom_point(pch=21, cex=6, aes(fill=Legend))+
  facet_grid(znaceni~horizont, scales="free")+theme_min+theme(legend.position = "top",
                                                              legend.key.size = unit(0.3, "in"),
                                                              legend.title = element_blank())+
  geom_errorbar(aes(ymin=y-y.sd, ymax=y+y.sd), width=0.1, lwd=0.5)+
  ylab(expression(paste(V[MAX],~beta,"-glucosidase ( ",mu, "mol ",g^{-1}~h^{-1}, ")" )))+
  xlab("Plesne : Certovo mixing ratio")+scale_fill_manual(values = c("white", "grey"))

##cellobiosidase
mich2e %>% group_by(Plesne, znaceni, horizont, Legend) %>% summarize(y.sd=sd(celVmax, na.rm = T), 
                                                                     y=mean(celVmax, na.rm = T)) %>%
  ggplot(aes(factor(Plesne), y))+geom_point(pch=21, cex=6, aes(fill=Legend))+
  facet_grid(znaceni~horizont, scales="free")+theme_min+theme(legend.position = c(0.17, 0.9),
                                                              legend.key.size = unit(0.3, "in"),
                                                              legend.title = element_blank())+
  geom_errorbar(aes(ymin=y-y.sd, ymax=y+y.sd), width=0.1, lwd=0.5)+
  ylab(expression(paste(V[MAX]," cellobiosidase ( ",mu, "mol ",g^{-1}~h^{-1}, ")" )))+
  xlab("Plesne : Certovo mixing ratio")+scale_fill_manual(values = c("white", "grey"))

##leucin-aminopeptidase
mich2e %>% group_by(Plesne, znaceni, horizont, Legend) %>% summarize(y.sd=sd(leuVmax), y=mean(leuVmax)) %>%
  ggplot(aes(factor(Plesne), y))+geom_point(pch=21, cex=6, aes(fill=Legend))+
  facet_grid(znaceni~horizont, scales="free")+theme_min+theme(legend.position = c(0.2, 0.92),
                                                              legend.key.size = unit(0.3, "in"),
                                                              legend.title = element_blank())+
  geom_errorbar(aes(ymin=y-y.sd, ymax=y+y.sd), width=0.1, lwd=0.5)+
  ylab(expression(paste(V[MAX]," leucin aminopeptidase ( ",mu, "mol ",g^{-1}~h^{-1}, ")" )))+
  xlab("Plesne : Certovo mixing ratio")+scale_fill_manual(values = c("white", "grey"))

##chitinase
mich2e %>% group_by(Plesne, znaceni, horizont, Legend) %>% summarize(y.sd=sd(chitVmax), y=mean(chitVmax)) %>%
  ggplot(aes(factor(Plesne), y))+geom_point(pch=21, cex=6, aes(fill=Legend))+
  facet_grid(znaceni~horizont, scales="free")+theme_min+theme(legend.position = "top",
                                                              legend.key.size = unit(0.3, "in"),
                                                              legend.title = element_blank())+
  geom_errorbar(aes(ymin=y-y.sd, ymax=y+y.sd), width=0.1, lwd=0.5)+
  ylab(expression(paste(V[MAX]," chitinase ( ", mu, "mol ",g^{-1}~h^{-1}, ")" )))+
  xlab("Plesne : Certovo mixing ratio")+scale_fill_manual(values = c("white", "grey"))

##phosphatase
mich2e %>% group_by(Plesne, znaceni, horizont, Legend) %>% summarize(y.sd=sd(phosVmax), y=mean(phosVmax)) %>%
  ggplot(aes(factor(Plesne), y))+geom_point(pch=21, cex=6, aes(fill=Legend))+
  facet_grid(znaceni~horizont, scales="free")+theme_min+theme(legend.position = "top",
                                                              legend.key.size = unit(0.3, "in"),
                                                              legend.title = element_blank())+
  geom_errorbar(aes(ymin=y-y.sd, ymax=y+y.sd), width=0.1, lwd=0.5)+
  ylab(expression(paste(V[MAX]," phosphatase ( ", mu, "mol ",g^{-1}~h^{-1}, ")" )))+
  xlab("Plesne : Certovo mixing ratio")+scale_fill_manual(values = c("white", "grey"))+
  geom_hline(yintercept = 0)

##phosphatase Km
mich2e %>% group_by(Plesne, znaceni, horizont, Legend) %>% summarize(y.sd=sd(phosKm), y=mean(phosKm)) %>%
  ggplot(aes(factor(Plesne), y))+geom_point(pch=21, cex=6, aes(fill=Legend))+
  facet_grid(znaceni~horizont, scales="free")+theme_min+theme(legend.position = c(0.2, 0.9),
                                                              legend.key.size = unit(0.3, "in"),
                                                              legend.title = element_blank())+
  geom_errorbar(aes(ymin=y-y.sd, ymax=y+y.sd), width=0.1, lwd=0.5)+
  ylab(expression(paste(K[M]," phosphatase ( ", mu, "mol ",g^{-1}~h^{-1}, ")" )))+
  xlab("Plesne : Certovo mixing ratio")+scale_fill_manual(values = c("white", "grey"))



PKMs<-as.data.frame(mich2e %>% group_by(Plesne, znaceni, horizont, Legend) %>% summarize(y=mean(phosVmax)))
PKMs$PO4<-as.data.frame(mich2e %>% group_by(Plesne, znaceni, horizont, Legend) %>% summarize(mean(Cmic)))[,5]

ggplot(PKMs[PKMs$Legend=="Before incubation", ], aes(PO4, y))+geom_point(cex=6, aes(colour=Legend))+
  theme_min+facet_wrap(~horizont, scales = "free")+stat_smooth(method = lm)


mich2e %>% filter(znaceni=="Unlabelled" & Legend=="Before incubation") %>% 
  group_by(Plesne, horizont) %>% summarize(y.sd=sd(phosVmax), 
                                           y=mean(phosVmax),
                                           x.sd = sd(NH42/PO42), x= mean(NH42/PO42)) %>%
  ggplot(aes(x, y))+geom_point(pch=21, cex=6, fill="grey")+
  theme_min+theme(legend.position = c(0.2, 0.9),
                                                              legend.key.size = unit(0.3, "in"),
                                                              legend.title = element_blank())+
  geom_errorbar(aes(ymin=y-y.sd, ymax=y+y.sd), width=0.1, lwd=0.5)+
  geom_errorbarh(aes(xmin=x-x.sd, xmax=x+x.sd), width=0.1, lwd=0.5)+
  ylab(expression(paste(K[M]," phosphatase ( ", mu, "mol ",g^{-1}~h^{-1}, ")" )))+
  xlab("Plesne : Certovo mixing ratio")+scale_fill_manual(values = c("white", "grey"))+
  stat_smooth(method = lm, se=F, colour="black")

mich2e %>% filter(znaceni=="Unlabelled" & Legend=="Before incubation") %>% 
  group_by(Plesne, horizont) %>% summarize(y.sd=sd(phosVmax), 
                                           y=mean(phosVmax),
                                           x.sd = sd(Nmic/Pmic/NH4*PO4), x= mean(Nmic/Pmic/NH4*PO4)) %>%
  ggplot(aes(x, y))+geom_point(pch=21, cex=6, fill="grey")+
  theme_min+theme(legend.position = c(0.2, 0.9),
                  legend.key.size = unit(0.3, "in"),
                  legend.title = element_blank())+
  geom_errorbar(aes(ymin=y-y.sd, ymax=y+y.sd), width=0.1, lwd=0.5)+
  geom_errorbarh(aes(xmin=x-x.sd, xmax=x+x.sd), width=0.1, lwd=0.5)+
  ylab(expression(paste(K[M]," phosphatase ( ", mu, "mol ",g^{-1}~h^{-1}, ")" )))+
  xlab("Plesne : Certovo mixing ratio")+scale_fill_manual(values = c("white", "grey"))+
  stat_smooth(method = lm, se=F, colour="black")

ggplot(mich2e[(mich2e$znaceni=="Unlabelled" & mich2e$Legend=="Before incubation"), ],
       aes(Nmic/Pmic/NH4*PO4, phosVmax))+geom_point(cex=6, aes(colour=horizont))+theme_min


#Potreba 